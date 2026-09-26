#include "WgsPipeline.h"
#include "TextTable.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QProcessEnvironment>
#include <QThread>

#include "EnsemblLookup.h"
#include "TargetRegions.h"
#include "ToolLocator.h"
#include "VcfConverter.h"

namespace {

QString envOr(const char *name, const QString &fallback)
{
    const QString v = QProcessEnvironment::systemEnvironment().value(QLatin1String(name));
    return v.isEmpty() ? fallback : v;
}

} // namespace

WgsConfig WgsConfig::defaults(const QString &root)
{
    WgsConfig c;
    c.root = root;
    c.dataDir = root + "/data";
    c.reportsDir = root + "/reports";
    const QString refDir = envOr("GH_REF_DIR", root + "/reference");
    c.workDir = envOr("GH_WGS_WORK", root + "/wgs_work");
    c.refGenome = envOr("GH_REF_GENOME", refDir + "/human_g1k_v37.fasta");
    c.refMmi = envOr("GH_REF_MMI", refDir + "/human_g1k_v37.mmi");
    c.rsidLookup = envOr("GH_RSID_LOOKUP", c.dataDir + "/rsid_positions_grch37.json");
    c.outputGenome = c.dataDir + "/genome.txt";
    c.analysisBinary = ToolLocator::analysisBinary(root);
    c.threads = QThread::idealThreadCount();
    return c;
}

WgsPipeline::WgsPipeline(const WgsConfig &config, const Reporter &reporter)
    : m_cfg(config), m_rep(reporter)
{
    if (m_cfg.threads <= 0)
        m_cfg.threads = QThread::idealThreadCount();
}

QStringList WgsPipeline::prerequisites() const
{
    QStringList problems;
    const QStringList absent = ToolLocator::missing(ToolLocator::requiredWgsTools());
    if (!absent.isEmpty())
        problems << QStringLiteral("%1 not installed. %2")
                        .arg(absent.join(", "), ToolLocator::installHint());
    if (!QFileInfo(m_cfg.refGenome).exists())
        problems << QStringLiteral("reference genome not found: %1 (run: make setup)").arg(m_cfg.refGenome);
    if (!m_cfg.skipAnalysis && !QFileInfo(m_cfg.analysisBinary).isExecutable())
        problems << QStringLiteral("analysis binary not found: %1 (run: make -C c)").arg(m_cfg.analysisBinary);
    return problems;
}

QString WgsPipeline::fail(const QString &what)
{
    m_error = what;
    m_rep.log(QStringLiteral("ERROR: %1").arg(what));
    return what;
}

bool WgsPipeline::waitFor(QProcess &p)
{
    while (!p.waitForFinished(200)) {
        if (m_rep.cancelled()) {
            p.kill();
            p.waitForFinished(2000);
            fail(QStringLiteral("cancelled"));
            return false;
        }
        // Drain stderr as it arrives so long steps show progress.
        const QByteArray err = p.readAllStandardError();
        for (const QByteArray &line : err.split('\n'))
            if (!line.trimmed().isEmpty())
                m_rep.log(QString::fromUtf8(line.trimmed()));
    }
    const QByteArray err = p.readAllStandardError();
    for (const QByteArray &line : err.split('\n'))
        if (!line.trimmed().isEmpty())
            m_rep.log(QString::fromUtf8(line.trimmed()));
    return p.exitStatus() == QProcess::NormalExit && p.exitCode() == 0;
}

bool WgsPipeline::runProcess(const QString &program, const QStringList &args,
                             QString *capturedStdout, bool logStdout)
{
    const QString exe = ToolLocator::find(program).isEmpty() ? program : ToolLocator::find(program);
    m_rep.log(QStringLiteral("$ %1 %2").arg(program, args.join(' ')));
    QProcess p;
    p.setProcessEnvironment(ToolLocator::environmentFor(ToolLocator::wgsTools()));
    p.start(exe, args);
    if (!p.waitForStarted(10000)) {
        fail(QStringLiteral("cannot start %1: %2").arg(program, p.errorString()));
        return false;
    }
    const bool ok = waitFor(p);
    const QString out = QString::fromUtf8(p.readAllStandardOutput());
    if (capturedStdout)
        *capturedStdout = out;
    else if (logStdout)
        for (const QString &line : out.split('\n'))
            if (!line.trimmed().isEmpty())
                m_rep.log(QStringLiteral("  %1").arg(line.trimmed()));
    if (!ok && m_error.isEmpty())
        fail(QStringLiteral("%1 exited with code %2").arg(program).arg(p.exitCode()));
    return ok;
}

bool WgsPipeline::runPipe(const QString &program1, const QStringList &args1,
                          const QString &program2, const QStringList &args2)
{
    m_rep.log(QStringLiteral("$ %1 %2 | %3 %4").arg(program1, args1.join(' '), program2, args2.join(' ')));
    QProcess p1, p2;
    const QProcessEnvironment env = ToolLocator::environmentFor(ToolLocator::wgsTools());
    p1.setProcessEnvironment(env);
    p2.setProcessEnvironment(env);
    p1.setStandardOutputProcess(&p2);
    p1.start(ToolLocator::find(program1), args1);
    p2.start(ToolLocator::find(program2), args2);
    if (!p1.waitForStarted(10000) || !p2.waitForStarted(10000)) {
        fail(QStringLiteral("cannot start %1 | %2").arg(program1, program2));
        return false;
    }
    const bool ok1 = waitFor(p1);
    const bool ok2 = waitFor(p2);
    if (!ok1 || !ok2) {
        if (m_error.isEmpty())
            fail(QStringLiteral("%1 | %2 failed (exit %3, %4)")
                     .arg(program1, program2).arg(p1.exitCode()).arg(p2.exitCode()));
        return false;
    }
    return true;
}

bool WgsPipeline::qc(const QString &fastq, QString *trimmed)
{
    const QString out = m_cfg.workDir + "/trimmed.fastq.gz";
    const QString json = m_cfg.workDir + "/fastp.json";
    m_rep.log(QStringLiteral("Step 2/6: Quality control (fastp)"));
    if (!runProcess("fastp", {"-i", fastq, "-o", out, "-j", json, "-h", m_cfg.workDir + "/fastp.html",
                              "-w", QString::number(qMin(m_cfg.threads, 16)), "--length_required", "50"}))
        return false;
    QFile jf(json);
    if (jf.open(QIODevice::ReadOnly)) {
        const QJsonObject summary = QJsonDocument::fromJson(jf.readAll()).object().value("summary").toObject();
        const double before = summary.value("before_filtering").toObject().value("total_reads").toDouble();
        const double after = summary.value("after_filtering").toObject().value("total_reads").toDouble();
        const double q30 = summary.value("after_filtering").toObject().value("q30_rate").toDouble();
        if (before > 0)
            m_rep.log(QStringLiteral("  Reads: %L1 -> %L2 (%3% passed, Q30: %4%)")
                          .arg(qint64(before)).arg(qint64(after))
                          .arg(after / before * 100, 0, 'f', 1).arg(q30 * 100, 0, 'f', 1));
    }
    *trimmed = out;
    return true;
}

bool WgsPipeline::alignAndSort(const QString &fastq, QString *bam)
{
    const QString out = m_cfg.workDir + "/aligned.bam";
    const bool haveIndex = QFileInfo(m_cfg.refMmi).exists();
    const QString ref = haveIndex ? m_cfg.refMmi : m_cfg.refGenome;
    if (haveIndex)
        m_rep.log(QStringLiteral("  Using pre-built index: %1").arg(QFileInfo(m_cfg.refMmi).fileName()));
    const int alignThreads = qMax(1, m_cfg.threads - 4);
    const int sortThreads = qMin(4, qMax(1, m_cfg.threads / 3));

    m_rep.log(QStringLiteral("Step 3/6: Alignment + sorting (minimap2 -> samtools sort)"));
    m_rep.stage(QStringLiteral("Aligning"), 10);
    if (!runPipe("minimap2", {"-a", "-x", "sr", "-t", QString::number(alignThreads), "--secondary=no", ref, fastq},
                 "samtools", {"sort", "-@", QString::number(sortThreads), "-m", "4G", "-o", out}))
        return false;
    m_rep.log(QStringLiteral("  Indexing BAM"));
    if (!runProcess("samtools", {"index", "-@", QString::number(m_cfg.threads), out}))
        return false;
    QString flagstat;
    if (!runProcess("samtools", {"flagstat", out}, &flagstat))
        return false;
    for (const QString &line : flagstat.split('\n'))
        if (line.contains("in total") || line.contains("mapped (") || line.contains("primary"))
            m_rep.log(QStringLiteral("  %1").arg(line.trimmed()));
    *bam = out;
    return true;
}

bool WgsPipeline::buildTargets(QString *bed)
{
    const QString path = m_cfg.workDir + "/target_regions.bed";
    if (QFileInfo(path).exists()) {
        *bed = path;
        return true;
    }
    QString err;
    const qint64 n = TargetRegions::build(m_cfg.dataDir + "/clinvar_alleles.tsv", m_cfg.rsidLookup,
                                          path, m_rep, &err);
    if (n < 0) {
        fail(err);
        return false;
    }
    *bed = n > 0 ? path : QString();
    return true;
}

bool WgsPipeline::callVariants(const QString &bam, QString *vcf)
{
    const QString out = m_cfg.workDir + "/variants.vcf.gz";
    QStringList mpileup{"mpileup", "-f", m_cfg.refGenome, "-q", "20", "-Q", "20",
                        "--threads", QString::number(qMax(1, m_cfg.threads / 2)), "-a", "FORMAT/DP"};
    if (!m_cfg.full) {
        QString bed;
        if (!buildTargets(&bed))
            return false;
        if (!bed.isEmpty()) {
            mpileup << "-T" << bed;
            m_rep.log(QStringLiteral("  Using targeted calling (ClinVar + SNP positions only)"));
        }
    }
    mpileup << bam;
    m_rep.log(QStringLiteral("Step 4/6: Variant calling (bcftools mpileup + call)"));
    m_rep.stage(QStringLiteral("Calling variants"), 55);
    if (!runPipe("bcftools", mpileup,
                 "bcftools", {"call", "-m", "-v", "--ploidy", "GRCh37",
                              "--threads", QString::number(qMax(1, m_cfg.threads / 4)), "-Oz", "-o", out}))
        return false;
    if (!runProcess("bcftools", {"index", out}))
        return false;
    QString stats;
    if (!runProcess("bcftools", {"stats", out}, &stats))
        return false;
    for (const QString &line : stats.split('\n')) {
        if (line.startsWith("SN") && line.contains("number of records"))
            m_rep.log(QStringLiteral("  Variants called: %1").arg(line.section('\t', -1).trimmed()));
        if (line.startsWith("SN") && line.contains("number of SNPs"))
            m_rep.log(QStringLiteral("  SNPs: %1").arg(line.section('\t', -1).trimmed()));
    }
    *vcf = out;
    return true;
}

qint64 WgsPipeline::convert(const QString &vcf, const QString &genomeOut)
{
    m_rep.stage(QStringLiteral("Converting"), 80);
    QString err;
    m_rep.log(QStringLiteral("Step 5/6: Converting the VCF to the array layout"));
    const qint64 n = VcfConverter::convert(vcf, genomeOut, m_cfg.rsidLookup, m_rep, &err);
    if (n < 0)
        fail(err);
    return n;
}

bool WgsPipeline::analyse(const QString &genome)
{
    m_rep.log(QStringLiteral("Step 6/6: Running genetic health analysis"));
    m_rep.stage(QStringLiteral("Analysing"), 90);
    QDir().mkpath(m_cfg.reportsDir);
    QStringList args{"--data", m_cfg.dataDir, "--html", m_cfg.reportsDir + "/GENETIC_HEALTH_REPORT.html", "--quiet"};
    if (!m_cfg.subjectName.isEmpty())
        args << "--name" << m_cfg.subjectName;
    args << genome;
    if (!runProcess(m_cfg.analysisBinary, args))
        return false;
    QString json;
    if (!runProcess(m_cfg.analysisBinary, {"--data", m_cfg.dataDir, "--json", "--quiet", genome}, &json))
        return false;
    QFile jf(m_cfg.reportsDir + "/comprehensive_results.json");
    if (jf.open(QIODevice::WriteOnly | QIODevice::Truncate))
        jf.write(json.toUtf8());
    m_rep.log(QStringLiteral("  Report: %1/GENETIC_HEALTH_REPORT.html").arg(m_cfg.reportsDir));
    return true;
}

bool WgsPipeline::run(WgsOutcome *outcome, QString *error)
{
    m_error.clear();
    const QStringList problems = prerequisites();
    if (!problems.isEmpty()) {
        if (error)
            *error = problems.join("; ");
        for (const QString &p : problems)
            m_rep.log(QStringLiteral("ERROR: %1").arg(p));
        return false;
    }
    QDir().mkpath(m_cfg.workDir);

    if (m_cfg.backupExisting && QFileInfo(m_cfg.outputGenome).exists()) {
        const QString backup = m_cfg.outputGenome.section('.', 0, -2) + ".txt.bak";
        QFile::remove(backup);
        QFile::copy(m_cfg.outputGenome, backup);
        m_rep.log(QStringLiteral("Backed up existing genome -> %1").arg(QFileInfo(backup).fileName()));
    }

    if (m_cfg.lookupRsids) {
        QString ids;
        if (!runProcess(m_cfg.analysisBinary, {"--list-rsids"}, &ids))
            return false;
        const QStringList rsids = outputLines(ids);
        m_rep.log(QStringLiteral("Step 1/6: rsID positions (Ensembl GRCh37, cached)"));
        m_rep.stage(QStringLiteral("Looking up rsID positions"), 5);
        if (!EnsemblLookup::update(m_cfg.rsidLookup, rsids, m_rep, nullptr)) {
            if (error)
                *error = QStringLiteral("rsID lookup failed");
            return false;
        }
    }

    QString input = m_cfg.fastq;
    if (m_cfg.skipQc || ToolLocator::find("fastp").isEmpty()) {
        if (ToolLocator::find("fastp").isEmpty())
            m_rep.log(QStringLiteral("Skipping QC (fastp not installed)"));
    } else if (!qc(m_cfg.fastq, &input)) {
        if (error) *error = m_error;
        return false;
    }

    WgsOutcome out;
    if (!alignAndSort(input, &out.bam) || !callVariants(out.bam, &out.vcf)) {
        if (error) *error = m_error;
        return false;
    }
    out.variants = convert(out.vcf, m_cfg.outputGenome);
    if (out.variants < 0) {
        if (error) *error = m_error;
        return false;
    }
    out.genome = m_cfg.outputGenome;

    if (!m_cfg.skipAnalysis && !analyse(out.genome)) {
        if (error) *error = m_error;
        return false;
    }

    if (!m_cfg.keepIntermediates) {
        m_rep.log(QStringLiteral("Cleaning up intermediate files..."));
        QDir work(m_cfg.workDir);
        for (const QString &f : work.entryList({"trimmed.*", "fastp.*"}, QDir::Files))
            work.remove(f);
    }
    m_rep.stage(QStringLiteral("Done"), 100);
    m_rep.log(QStringLiteral("Pipeline complete: %L1 SNPs extracted -> %2").arg(out.variants).arg(out.genome));
    if (outcome)
        *outcome = out;
    return true;
}
