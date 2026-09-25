#include "WgsValidator.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonObject>
#include <QProcess>

#include <algorithm>
#include <random>
#include <vector>

#include "Downloader.h"
#include "EnsemblLookup.h"
#include "TextTable.h"
#include "ToolLocator.h"
#include "WgsPipeline.h"

namespace WgsValidator {

namespace {

const char *const kChrom = "22";
const char *const kFastaUrl =
    "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz";
const int kReadLen = 150;
const int kDepth = 30;

struct Target {
    qint64 pos;       // 1-based
    QString rsid;     // empty for a ClinVar-only position
    char ref;
};

struct Truth {
    QString rsid;
    char ref, alt;
    QString genotype;   // sorted pair for het
    QString zygosity;   // hom_ref / het / hom_alt
    char hap[2];
};

QByteArray loadChromosome(const QString &path)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        return {};
    QByteArray seq;
    seq.reserve(52 * 1024 * 1024);
    while (!f.atEnd()) {
        QByteArray line = f.readLine();
        if (line.startsWith('>'))
            continue;
        line = line.trimmed().toUpper();
        seq += line;
    }
    return seq;
}

std::vector<Target> pickTargets(const QByteArray &seq, const QString &dataDir, int limit,
                                const Reporter &reporter)
{
    std::vector<std::pair<qint64, QString>> raw;
    const QJsonObject lookup = EnsemblLookup::loadLookup(dataDir + "/rsid_positions_grch37.json");
    for (auto it = lookup.begin(); it != lookup.end(); ++it) {
        const QJsonObject info = it.value().toObject();
        if (info.value("chrom").toString() == QLatin1String(kChrom))
            raw.push_back({info.value("pos").toString().toLongLong(), it.key()});
    }
    QFile cv(dataDir + "/clinvar_alleles.tsv");
    if (cv.open(QIODevice::ReadOnly)) {
        TsvTable t(&cv);
        const int cChrom = t.column("chrom"), cPos = t.column("pos"), cRef = t.column("ref"), cAlt = t.column("alt");
        int seen = 0;
        while (seen < limit && t.next()) {
            if (t.field(cChrom) != QLatin1String(kChrom))
                continue;
            if (t.field(cRef).size() != 1 || t.field(cAlt).size() != 1)
                continue;
            raw.push_back({t.field(cPos).toLongLong(), QString()});
            seen++;
        }
    }
    std::sort(raw.begin(), raw.end());
    raw.erase(std::unique(raw.begin(), raw.end()), raw.end());

    std::vector<Target> usable;
    for (const auto &[pos, rsid] : raw) {
        const qint64 start = pos - 1 - kReadLen, end = pos - 1 + kReadLen;
        if (start < 0 || end >= seq.size())
            continue;
        if (QByteArray::fromRawData(seq.constData() + start, int(end - start)).contains('N'))
            continue;
        usable.push_back({pos, rsid, seq[int(pos - 1)]});
    }
    int named = 0;
    for (const Target &t : usable)
        if (!t.rsid.isEmpty())
            named++;
    reporter.log(QStringLiteral("%1 usable loci on chr%2 (%3 with an rsID, %4 positional)")
                     .arg(usable.size()).arg(QLatin1String(kChrom)).arg(named).arg(usable.size() - named));
    return usable;
}

// Reads from two haplotypes. Every variant inside a read is applied, not
// just the one it was drawn for: chr22 targets sit as close as 30bp, and a
// read that only carried its own variant would drag a neighbouring hom-alt
// site down to half alt.
QHash<qint64, Truth> simulate(const QByteArray &seq, const std::vector<Target> &targets,
                              const QString &fastqPath, unsigned seed)
{
    std::mt19937_64 rng(seed);
    static const char bases[] = "ACGT";
    QHash<qint64, Truth> truth;
    std::vector<qint64> positions;

    for (size_t i = 0; i < targets.size(); ++i) {
        const Target &t = targets[i];
        char alt;
        do {
            alt = bases[std::uniform_int_distribution<int>(0, 3)(rng)];
        } while (alt == t.ref);
        Truth tr;
        tr.rsid = t.rsid;
        tr.ref = t.ref;
        tr.alt = alt;
        switch (i % 3) {
        case 0: tr.zygosity = "hom_ref"; tr.hap[0] = t.ref; tr.hap[1] = t.ref; tr.genotype = QString(2, QChar(t.ref)); break;
        case 1: tr.zygosity = "het"; tr.hap[0] = t.ref; tr.hap[1] = alt;
                tr.genotype = QString(QChar(qMin(t.ref, alt))) + QChar(qMax(t.ref, alt)); break;
        default: tr.zygosity = "hom_alt"; tr.hap[0] = alt; tr.hap[1] = alt; tr.genotype = QString(2, QChar(alt)); break;
        }
        truth[t.pos] = tr;
        positions.push_back(t.pos);
    }
    std::sort(positions.begin(), positions.end());

    std::vector<QByteArray> reads;
    for (const Target &t : targets) {
        for (int copy = 0; copy < kDepth; ++copy) {
            const int offset = std::uniform_int_distribution<int>(20, kReadLen - 21)(rng);
            const qint64 start = t.pos - 1 - offset;
            if (start < 0 || start + kReadLen > seq.size())
                continue;
            const int hap = copy % 2;
            QByteArray read = seq.mid(int(start), kReadLen);
            auto lo = std::lower_bound(positions.begin(), positions.end(), start + 1);
            for (auto it = lo; it != positions.end() && *it <= start + kReadLen; ++it)
                read[int(*it - 1 - start)] = truth[*it].hap[hap];
            reads.push_back("@sim:" + QByteArray::number(t.pos) + ':' + QByteArray::number(copy) +
                            ":h" + QByteArray::number(hap) + '\n' + read + "\n+\n" +
                            QByteArray(read.size(), 'I') + '\n');
        }
    }
    std::shuffle(reads.begin(), reads.end(), rng);
    QFile f(fastqPath);
    f.open(QIODevice::WriteOnly | QIODevice::Truncate);
    for (const QByteArray &r : reads)
        f.write(r);
    return truth;
}

bool runTool(const QString &tool, const QStringList &args, const Reporter &reporter, QString *error,
             const QString &stdoutFile = QString())
{
    const QString exe = ToolLocator::find(tool);
    if (exe.isEmpty()) {
        if (error)
            *error = QStringLiteral("%1 not found").arg(tool);
        return false;
    }
    reporter.log(QStringLiteral("$ %1 %2").arg(tool, args.join(' ')));
    QProcess p;
    if (!stdoutFile.isEmpty())
        p.setStandardOutputFile(stdoutFile, QIODevice::Truncate);
    p.start(exe, args);
    p.waitForFinished(-1);
    if (p.exitStatus() != QProcess::NormalExit || p.exitCode() != 0) {
        if (error)
            *error = QStringLiteral("%1 failed: %2").arg(tool, QString::fromUtf8(p.readAllStandardError()).trimmed());
        return false;
    }
    return true;
}

} // namespace

QString referenceDir(const QString &root) { return root + "/reference"; }

bool isSetUp(const QString &root)
{
    return QFileInfo(referenceDir(root) + "/chr22.mmi").exists() &&
           QFileInfo(referenceDir(root) + "/chr22.fa").exists();
}

bool setup(const QString &root, const Reporter &reporter, QString *error)
{
    const QString dir = referenceDir(root);
    QDir().mkpath(dir);
    const QString gz = dir + "/chr22.fa.gz", fa = dir + "/chr22.fa", mmi = dir + "/chr22.mmi";
    if (!QFileInfo(fa).exists()) {
        if (!QFileInfo(gz).exists()) {
            reporter.log(QStringLiteral("downloading %1").arg(QLatin1String(kFastaUrl)));
            if (!Downloader::download(QUrl(QLatin1String(kFastaUrl)), gz, reporter, error))
                return false;
        }
        reporter.log(QStringLiteral("decompressing -> %1").arg(fa));
        QProcess gunzip;
        gunzip.setStandardOutputFile(fa, QIODevice::Truncate);
        gunzip.start("gzip", {"-dc", gz});
        gunzip.waitForFinished(-1);
        if (gunzip.exitCode() != 0) {
            if (error)
                *error = QStringLiteral("gzip failed");
            return false;
        }
    }
    if (!QFileInfo(fa + ".fai").exists() && !runTool("samtools", {"faidx", fa}, reporter, error))
        return false;
    if (!QFileInfo(mmi).exists() && !runTool("minimap2", {"-x", "sr", "-d", mmi, fa}, reporter, error))
        return false;
    reporter.log(QStringLiteral("ready: %1").arg(fa));
    return true;
}

bool run(const Options &options, const Reporter &reporter, QStringList *problems, QString *error)
{
    if (!isSetUp(options.root)) {
        if (error)
            *error = QStringLiteral("chr22 sandbox not set up (run setup first)");
        return false;
    }
    const QString work = options.root + "/wgs_work_validate";
    QDir().mkpath(work);

    reporter.log(QStringLiteral("loading %1/chr22.fa").arg(referenceDir(options.root)));
    const QByteArray seq = loadChromosome(referenceDir(options.root) + "/chr22.fa");
    if (seq.isEmpty()) {
        if (error)
            *error = QStringLiteral("could not read chr22.fa");
        return false;
    }
    const std::vector<Target> targets = pickTargets(seq, options.root + "/data", options.clinvarTargets, reporter);
    const QString fastq = work + "/simulated.fastq";
    const QHash<qint64, Truth> truth = simulate(seq, targets, fastq, options.seed);
    int het = 0, homAlt = 0, homRef = 0;
    for (const Truth &t : truth) {
        if (t.zygosity == "het") het++;
        else if (t.zygosity == "hom_alt") homAlt++;
        else homRef++;
    }
    reporter.log(QStringLiteral("simulated %1 MB of reads: %2 het, %3 hom_alt, %4 hom_ref")
                     .arg(QFileInfo(fastq).size() / 1e6, 0, 'f', 1).arg(het).arg(homAlt).arg(homRef));

    WgsConfig cfg = WgsConfig::defaults(options.root);
    cfg.refGenome = referenceDir(options.root) + "/chr22.fa";
    cfg.refMmi = referenceDir(options.root) + "/chr22.mmi";
    cfg.workDir = work;
    cfg.skipQc = true;
    cfg.skipAnalysis = true;
    cfg.lookupRsids = false;
    cfg.backupExisting = false;
    cfg.keepIntermediates = true;
    cfg.outputGenome = work + "/genome_from_wgs.txt";
    if (options.threads > 0)
        cfg.threads = options.threads;
    // The BED is rebuilt each time so a changed ClinVar or lookup is reflected.
    QFile::remove(work + "/target_regions.bed");

    WgsPipeline pipeline(cfg, reporter);
    QString bam, vcf;
    if (!pipeline.alignAndSort(fastq, &bam) || !pipeline.callVariants(bam, &vcf)) {
        if (error)
            *error = QStringLiteral("pipeline step failed");
        return false;
    }
    if (pipeline.convert(vcf, cfg.outputGenome) < 0) {
        if (error)
            *error = QStringLiteral("conversion failed");
        return false;
    }

    // Read back what the pipeline produced.
    QHash<qint64, QPair<QString, QString>> called;   // pos -> (rsid, genotype)
    QFile g(cfg.outputGenome);
    if (g.open(QIODevice::ReadOnly)) {
        while (!g.atEnd()) {
            const QByteArray line = g.readLine().trimmed();
            if (line.startsWith('#'))
                continue;
            const QList<QByteArray> parts = line.split('\t');
            if (parts.size() >= 4 && parts[1] == kChrom)
                called[parts[2].toLongLong()] = {QString::fromUtf8(parts[0]), QString::fromUtf8(parts[3])};
        }
    }

    QStringList diffs;
    QList<qint64> keys = truth.keys();
    std::sort(keys.begin(), keys.end());
    for (qint64 pos : keys) {
        const Truth &want = truth[pos];
        const bool have = called.contains(pos);
        if (want.zygosity == "hom_ref") {
            // `bcftools call -v` emits variant sites only.
            if (have)
                diffs << QStringLiteral("%1:%2 is hom-ref but was called %3").arg(kChrom).arg(pos).arg(called[pos].second);
            continue;
        }
        if (!have) {
            diffs << QStringLiteral("%1:%2 %3 %4>%5 was not called at all")
                        .arg(kChrom).arg(pos).arg(want.zygosity).arg(want.ref).arg(want.alt);
            continue;
        }
        QString got = called[pos].second;
        std::sort(got.begin(), got.end());
        if (got != want.genotype)
            diffs << QStringLiteral("%1:%2 expected %3 got %4 (%5)").arg(kChrom).arg(pos)
                        .arg(want.genotype, called[pos].second, want.zygosity);
        else if (!want.rsid.isEmpty() && called[pos].first != want.rsid)
            diffs << QStringLiteral("%1:%2 expected rsID %3 got %4").arg(kChrom).arg(pos).arg(want.rsid, called[pos].first);
    }
    if (problems)
        *problems = diffs;
    reporter.stage(diffs.isEmpty() ? QStringLiteral("Validation passed") : QStringLiteral("Validation failed"), 100);
    if (diffs.isEmpty())
        reporter.log(QStringLiteral("ok: %1/%1 loci round-tripped FASTQ -> BAM -> VCF -> 23andMe with the expected genotype")
                         .arg(truth.size()));
    else
        reporter.log(QStringLiteral("FAILED: %1 of %2 loci differ").arg(diffs.size()).arg(truth.size()));
    return diffs.isEmpty();
}

} // namespace WgsValidator
