#include "TestGenome.h"
#include "TextTable.h"

#include <QCryptographicHash>
#include <QDateTime>
#include <QDir>
#include <QEventLoop>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMap>
#include <QSet>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QProcess>
#include <QRegularExpression>
#include <QTemporaryDir>
#include <QThread>
#include <QTimer>

#include <algorithm>

#include "ChainLift.h"
#include "Downloader.h"
#include "EnsemblLookup.h"
#include "RemoteFasta.h"
#include "ToolLocator.h"

namespace TestGenome {

namespace {

const char *const kBase = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/";

struct Site {
    qint64 pos;
    QString rsid;
    QString ref;              // from the VCF, or Ensembl for hom-ref sites
    QStringList alts;
    QString gt;               // "0|1", "1", ...
    qint64 offset = 0;        // pos - record POS, for a record spanning the site
    bool inVcf = false;
};

int chromOrder(const QString &c)
{
    bool ok;
    const int n = c.toInt(&ok);
    if (ok) return n;
    if (c == "X") return 23;
    if (c == "Y") return 24;
    if (c == "MT") return 25;
    return 26;
}

const char *const kSeqGRCh37 = "https://grch37.rest.ensembl.org/sequence/region/human";
const char *const kSeqGRCh38 = "https://rest.ensembl.org/sequence/region/human";

// Reference sequence for a batch of regions ("chrom:start..end:1"), via an
// Ensembl sequence endpoint (POST, up to 50 regions), keyed by the query.
QMap<QString, QString> referenceBases(QNetworkAccessManager &nam, const char *endpoint,
                                      const QStringList &regions, QString *error)
{
    QMap<QString, QString> out;
    for (int i = 0; i < regions.size(); i += 50) {
        QJsonArray arr;
        for (const QString &r : regions.mid(i, 50))
            arr.append(r);
        QByteArray body;
        QString err;
        for (int attempt = 0; attempt < 5 && body.isEmpty(); attempt++) {
            if (attempt)
                QThread::msleep(1500 << qMin(attempt - 1, 3));
            QNetworkRequest req{QUrl(QLatin1String(endpoint))};
            req.setHeader(QNetworkRequest::ContentTypeHeader, "application/json");
            req.setRawHeader("Accept", "application/json");
            QNetworkReply *reply = nam.post(req, QJsonDocument(QJsonObject{{"regions", arr}}).toJson(QJsonDocument::Compact));
            QEventLoop loop;
            QTimer timeout;
            timeout.setSingleShot(true);
            QObject::connect(&timeout, &QTimer::timeout, reply, &QNetworkReply::abort);
            QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
            timeout.start(90000);
            loop.exec();
            if (reply->error() == QNetworkReply::NoError)
                body = reply->readAll();
            else
                err = reply->errorString();
            reply->deleteLater();
        }
        if (body.isEmpty()) {
            if (error)
                *error = QStringLiteral("Ensembl sequence request failed: %1").arg(err);
            return out;
        }
        for (const QJsonValue &v : QJsonDocument::fromJson(body).array()) {
            const QJsonObject o = v.toObject();
            out[o.value("query").toString()] = o.value("seq").toString().toUpper();
        }
    }
    return out;
}

QString reverseComplement(const QString &seq)
{
    QString out;
    for (int i = seq.size() - 1; i >= 0; i--) {
        switch (seq.at(i).toLatin1()) {
        case 'A': out += 'T'; break;
        case 'C': out += 'G'; break;
        case 'G': out += 'C'; break;
        case 'T': out += 'A'; break;
        default:  out += 'N'; break;
        }
    }
    return out;
}

// Every rsID the analysis reads, from the C binary.
bool listRsids(const Options &options, QStringList *rsids, QString *error)
{
    QProcess list;
    list.start(options.analysisBinary, {"--list-rsids"});
    list.waitForFinished(-1);
    if (list.exitCode() != 0) {
        if (error)
            *error = QStringLiteral("could not run %1 --list-rsids").arg(options.analysisBinary);
        return false;
    }
    *rsids = outputLines(QString::fromUtf8(list.readAllStandardOutput()));
    return true;
}

// The rsID -> GRCh37 position lookup for `rsids`, refreshed from Ensembl
// for ids it lacks. A packaged build starts from the copy compiled in, so
// the Ensembl host being unreachable only matters for ids added since; it
// is fatal only when nothing at all is positioned, since a sample built
// from no positions would be an empty file.
bool positionLookup(const QString &lookupPath, const QStringList &rsids, const Reporter &reporter,
                    QJsonObject *lookup, QString *error)
{
    QString lookupError;
    const bool updated = EnsemblLookup::update(lookupPath, rsids, reporter, nullptr, &lookupError);
    *lookup = EnsemblLookup::loadLookup(lookupPath);
    int positioned = 0;
    for (const QString &r : rsids)
        if (lookup->contains(r))
            positioned++;
    if (positioned == 0) {
        if (error)
            *error = QStringLiteral("no GRCh37 position is known for any of the %1 rsIDs: %2")
                         .arg(rsids.size()).arg(updated ? QStringLiteral("the lookup is empty") : lookupError);
        return false;
    }
    if (!updated)
        reporter.log(QStringLiteral("  Warning: %1; continuing with the %2 positions on file")
                         .arg(lookupError).arg(positioned));
    return true;
}

// Reference bases on `options.sample` for the ids that have a GRCh37
// position, from the assembly sequence itself: `samtools faidx` on Ensembl's
// indexed FASTA, with GRCh37 positions lifted to GRCh38 through Ensembl's
// chain file. A site on the reverse strand in GRCh38 is complemented so the
// file stays in GRCh37 orientation, which is what the tables use. Ids the
// chain cannot place are returned in `unresolved`.
bool referenceFromAssembly(const Options &options, const QString &samtools, const QJsonObject &lookup,
                           const QStringList &rsids, const Reporter &reporter,
                           QMap<QString, QString> *refBase, QStringList *unresolved, QString *error)
{
    ChainLift chain;
    if (options.sample == "GRCh38") {
        const QString chainPath = options.dataDir + "/ensembl/GRCh37_to_GRCh38.chain.gz";
        if (!QFile::exists(chainPath)) {
            QDir().mkpath(QFileInfo(chainPath).path());
            reporter.log(QStringLiteral("  downloading GRCh37_to_GRCh38.chain.gz from Ensembl"));
            if (!Downloader::download(QUrl(QLatin1String(RemoteFasta::kChainGRCh37ToGRCh38)), chainPath, reporter, error))
                return false;
        }
        if (!chain.load(chainPath, error))
            return false;
    }

    QStringList regions;
    QMap<QString, QString> regionOf;
    QStringList flipped;
    for (const QString &r : rsids) {
        const QJsonObject info = lookup.value(r).toObject();
        if (info.isEmpty())
            continue;
        const QString chrom = info.value("chrom").toString();
        const qint64 pos = info.value("pos").toString().toLongLong();
        if (options.sample == "GRCh38") {
            ChainLift::Lifted l;
            if (!chain.lift(chrom, pos, &l)) {
                unresolved->append(r);
                continue;
            }
            regionOf[r] = RemoteFasta::region(l.chrom, l.pos);
            if (l.reverse)
                flipped << r;
        } else {
            regionOf[r] = RemoteFasta::region(chrom, pos);
        }
        regions << regionOf[r];
    }

    reporter.stage(QStringLiteral("Fetching %1 bases from the assembly (%2 regions)")
                       .arg(options.sample).arg(regions.size()), 10);
    const char *url = options.sample == "GRCh38" ? RemoteFasta::kGRCh38 : RemoteFasta::kGRCh37;
    QString err;
    const QMap<QString, QString> seq = RemoteFasta::bases(
        samtools, url, regions, reporter, &err, 4,
        options.dataDir + "/test_genome_cache/fasta_" + options.sample + ".json");
    if (!err.isEmpty()) {
        if (error) *error = err;
        return false;
    }
    for (auto it = regionOf.begin(); it != regionOf.end(); ++it) {
        QString b = seq.value(it.value());
        if (b.size() != 1) {
            unresolved->append(it.key());
            continue;
        }
        if (flipped.contains(it.key()))
            b = reverseComplement(b);
        (*refBase)[it.key()] = b;
    }
    if (!flipped.isEmpty())
        reporter.log(QStringLiteral("  %1 site(s) on the opposite strand in %2, written in GRCh37 orientation: %3")
                         .arg(flipped.size()).arg(options.sample, flipped.join(", ")));
    return true;
}

// The same through the Ensembl REST API: the first allele of a variation
// record's mapping is the reference base on that build's forward strand.
// Used when samtools is not installed, and for the few sites the chain file
// does not cover (MSMB rs10993994 lies in a region GRCh38 inverted). The
// GRCh38 host goes through spells of failing, hence the retries and the
// per-id fallback; ids that still fail are left as no-calls and named.
bool referenceFromRest(const Options &options, const QJsonObject &lookup, const QStringList &rsids,
                       QNetworkAccessManager &nam, const Reporter &reporter,
                       QMap<QString, QString> *refBase, QString *error)
{
    const char *endpoint = options.sample == "GRCh38" ? EnsemblLookup::kEndpointGRCh38
                                                      : EnsemblLookup::kEndpoint;
    QMap<QString, QString> buildSite;   // rsid -> "chrom:pos" on the build
    const int batchSize = 100;
    for (int i = 0; i < rsids.size(); i += batchSize) {
        if (reporter.cancelled()) {
            if (error) *error = "cancelled";
            return false;
        }
        const QStringList batch = rsids.mid(i, batchSize);
        reporter.stage(QStringLiteral("Fetching %1 reference alleles from Ensembl (%2-%3 of %4)")
                           .arg(options.sample).arg(i + 1).arg(qMin(i + batchSize, rsids.size())).arg(rsids.size()),
                       5 + 85 * i / qMax(1, rsids.size()));
        QJsonObject results;
        QString err;
        bool ok = false;
        for (int attempt = 0; attempt < 6 && !ok; attempt++) {
            if (attempt) {
                reporter.log(QStringLiteral("  retrying batch %1 (%2)").arg(i / batchSize + 1).arg(err));
                QThread::msleep(2000 << qMin(attempt - 1, 3));
            }
            ok = EnsemblLookup::postVariationBatch(nam, endpoint, batch, &results, &err);
        }
        if (!ok) {
            reporter.log(QStringLiteral("  batch POST keeps failing (%1); fetching the %2 ids one by one")
                             .arg(err).arg(batch.size()));
            results = QJsonObject();
            QStringList failed;
            for (const QString &id : batch) {
                if (reporter.cancelled()) {
                    if (error) *error = "cancelled";
                    return false;
                }
                QJsonObject one;
                bool got = false;
                for (int attempt = 0; attempt < 4 && !got; attempt++) {
                    if (attempt)
                        QThread::msleep(1000 << attempt);
                    got = EnsemblLookup::getVariation(nam, endpoint, id, &one, &err);
                }
                if (got)
                    results[id] = one;
                else
                    failed << id;
                QThread::msleep(70);   // Ensembl's 15 requests/second
            }
            if (!failed.isEmpty())
                reporter.log(QStringLiteral("  Warning: no %1 record could be fetched for %2; written as no-calls")
                                 .arg(options.sample, failed.join(", ")));
        }
        const QSet<QString> requested(batch.begin(), batch.end());
        for (auto it = results.begin(); it != results.end(); ++it) {
            const QJsonObject info = it.value().toObject();
            const QJsonObject mapping = EnsemblLookup::primaryMapping(info);
            if (mapping.isEmpty())
                continue;
            const QString ref = mapping.value("allele_string").toString().section('/', 0, 0);
            QSet<QString> names{it.key(), info.value("name").toString()};
            for (const QJsonValue &syn : info.value("synonyms").toArray())
                names.insert(syn.toString());
            for (const QString &name : names)
                if (requested.contains(name)) {
                    (*refBase)[name] = ref;
                    buildSite[name] = mapping.value("chrom").toString() + ":"
                        + QString::number(qint64(mapping.value("start").toDouble()));
                }
        }
        QThread::msleep(200);
    }

    // Strand: compare an 11-base flank on both builds; where the GRCh38
    // flank is the reverse complement of the GRCh37 one, complement the base.
    if (options.sample != "GRCh37" && !buildSite.isEmpty()) {
        QStringList regions37, regions38;
        QMap<QString, QString> key37, key38;
        for (auto it = buildSite.begin(); it != buildSite.end(); ++it) {
            const QJsonObject l = lookup.value(it.key()).toObject();
            if (l.isEmpty())
                continue;
            const qint64 p37 = l.value("pos").toString().toLongLong();
            const qint64 p38 = it.value().section(':', 1).toLongLong();
            key37[it.key()] = QStringLiteral("%1:%2..%3:1").arg(l.value("chrom").toString()).arg(p37 - 5).arg(p37 + 5);
            key38[it.key()] = QStringLiteral("%1:%2..%3:1").arg(it.value().section(':', 0, 0)).arg(p38 - 5).arg(p38 + 5);
            regions37 << key37[it.key()];
            regions38 << key38[it.key()];
        }
        QString e37, e38;
        const QMap<QString, QString> flank37 = referenceBases(nam, kSeqGRCh37, regions37, &e37);
        const QMap<QString, QString> flank38 = referenceBases(nam, kSeqGRCh38, regions38, &e38);
        if (!e37.isEmpty() || !e38.isEmpty()) {
            if (error)
                *error = QStringLiteral("strand check failed: %1 %2").arg(e37, e38);
            return false;
        }
        QStringList flipped;
        for (auto it = key37.begin(); it != key37.end(); ++it) {
            const QString a = flank37.value(it.value()), b = flank38.value(key38[it.key()]);
            if (a.size() != 11 || b.size() != 11 || a == b)
                continue;
            if (b == reverseComplement(a)) {
                (*refBase)[it.key()] = reverseComplement(refBase->value(it.key()));
                flipped << it.key();
            }
        }
        if (!flipped.isEmpty())
            reporter.log(QStringLiteral("  %1 site(s) on the opposite strand in %2, written in GRCh37 orientation: %3")
                             .arg(flipped.size()).arg(options.sample, flipped.join(", ")));
    }
    return true;
}

// Apply one chromosome's bcftools rows to its sites. A site is covered by
// the record starting at it, or failing that by an earlier record whose REF
// spans it. Returns how many sites got a record.
int applyRows(const QByteArray &output, QList<Site> &sites)
{
    QMap<qint64, QStringList> rows;
    for (const QByteArray &line : outputLines(output)) {
        const QList<QByteArray> f = line.split('\t');
        if (f.size() >= 5)
            rows[f[1].toLongLong()] = QStringList{QString::fromUtf8(f[2]), QString::fromUtf8(f[3]),
                                                  QString::fromUtf8(f[4]).trimmed()};
    }
    int hits = 0;
    for (Site &s : sites) {
        auto it = rows.find(s.pos);
        if (it == rows.end()) {
            auto prev = rows.lowerBound(s.pos);
            while (prev != rows.begin()) {
                --prev;
                if (prev.key() + prev.value()[0].size() > s.pos) {
                    it = prev;
                    break;
                }
                if (s.pos - prev.key() > 64)
                    break;
            }
        }
        if (it == rows.end())
            continue;
        const QStringList &r = it.value();
        s.ref = r[0];
        s.alts = r[1].split(',');
        s.gt = r[2];
        s.offset = s.pos - it.key();
        s.inVcf = true;
        hits++;
    }
    return hits;
}

QList<qint64> positionsOf(const QList<Site> &sites)
{
    QList<qint64> out;
    for (const Site &s : sites)
        out << s.pos;
    return out;
}

// The reference assembly as a sample. Positions stay GRCh37 so the
// annotation data keeps matching by coordinate; the bases are the build's.
bool buildReference(const Options &options, const Reporter &reporter, Stats *stats, QString *error)
{
    const QString output = options.output.isEmpty()
        ? options.dataDir + "/genome_reference_" + options.sample + ".txt" : options.output;

    reporter.stage(QStringLiteral("Collecting positions"), 2);
    QStringList rsids;
    if (!listRsids(options, &rsids, error))
        return false;
    const QString lookupPath = options.dataDir + "/rsid_positions_grch37.json";
    QJsonObject lookup;
    if (!positionLookup(lookupPath, rsids, reporter, &lookup, error))
        return false;

    QMap<QString, QList<Site>> byChrom;
    QStringList placed;
    Stats st;
    for (const QString &r : rsids) {
        const QJsonObject info = lookup.value(r).toObject();
        if (info.isEmpty())
            continue;
        Site s;
        s.pos = info.value("pos").toString().toLongLong();
        s.rsid = r;
        byChrom[info.value("chrom").toString()] << s;
        placed << r;
        st.positions++;
    }

    QNetworkAccessManager nam;
    QMap<QString, QString> refBase;   // rsid -> the build's reference base, GRCh37 orientation
    const QString samtools = ToolLocator::find("samtools");
    if (!samtools.isEmpty()) {
        QStringList leftover;
        if (!referenceFromAssembly(options, samtools, lookup, placed, reporter, &refBase, &leftover, error))
            return false;
        if (!leftover.isEmpty()) {
            reporter.log(QStringLiteral("  %1 site(s) the chain file does not place; asking the Ensembl REST API: %2")
                             .arg(leftover.size()).arg(leftover.join(", ")));
            QString err;
            if (!referenceFromRest(options, lookup, leftover, nam, reporter, &refBase, &err))
                reporter.log(QStringLiteral("  Warning: %1; left as no-calls").arg(err));
        }
    } else {
        reporter.log(QStringLiteral("  samtools not found (%1); using the Ensembl REST API instead")
                         .arg(ToolLocator::installHint()));
        if (!referenceFromRest(options, lookup, placed, nam, reporter, &refBase, error))
            return false;
    }
    if (refBase.isEmpty()) {
        if (error)
            *error = QStringLiteral("no reference base could be fetched for any of the %1 positions; "
                                    "nothing written (is the network reachable?)").arg(st.positions);
        return false;
    }

    QFile out(output);
    if (!out.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (error)
            *error = QStringLiteral("cannot write %1").arg(output);
        return false;
    }
    out.write(QStringLiteral(
        "# The %1 reference assembly as a sample: its own base at each of the %2 positions\n"
        "# this analysis reads, homozygous. Positions are GRCh37 so the annotation data\n"
        "# matches by coordinate; the bases are %1's (Ensembl%3). Built by gh-data test-genome.\n"
        "# rsid\tchromosome\tposition\tgenotype\n")
        .arg(options.sample).arg(st.positions)
        .arg(samtools.isEmpty() ? " REST" : options.sample == "GRCh38" ? " assembly FASTA, chain-lifted" : " assembly FASTA")
        .toUtf8());
    QStringList chroms = byChrom.keys();
    std::sort(chroms.begin(), chroms.end(), [](const QString &a, const QString &b) { return chromOrder(a) < chromOrder(b); });
    for (const QString &chrom : chroms) {
        auto &sites = byChrom[chrom];
        std::sort(sites.begin(), sites.end(), [](const Site &a, const Site &b) { return a.pos < b.pos; });
        const bool haploid = chrom == "MT" || chrom == "Y";
        for (const Site &s : sites) {
            const QString base = refBase.value(s.rsid);
            QString genotype = "--";
            if (base.size() == 1 && QString("ACGT").contains(base)) {
                genotype = haploid ? base : base + base;
                st.homRef++;
            } else {
                st.noCall++;   // no record on this build, or an indel/repeat site
            }
            out.write((s.rsid + '\t' + chrom + '\t' + QString::number(s.pos) + '\t' + genotype + '\n').toUtf8());
        }
    }
    reporter.log(QStringLiteral("wrote %1: %2 reference genotypes, %3 no-calls")
                     .arg(output).arg(st.homRef).arg(st.noCall));
    reporter.stage(QStringLiteral("Done"), 100);
    if (stats)
        *stats = st;
    return true;
}

} // namespace

bool isReferenceSample(const QString &sample)
{
    return sample == "GRCh37" || sample == "GRCh38";
}

QByteArray chromCacheHeader(const QList<qint64> &positions)
{
    QByteArray text;
    for (qint64 p : positions)
        text += QByteArray::number(p) + ',';
    const QByteArray digest = QCryptographicHash::hash(text, QCryptographicHash::Sha1).toHex();
    return "#positions " + QByteArray::number(positions.size()) + " " + digest + "\n";
}

bool chromCacheMatches(const QByteArray &text, const QList<qint64> &positions, QByteArray *body)
{
    const QByteArray header = chromCacheHeader(positions);
    if (!text.startsWith(header))
        return false;
    if (body)
        *body = text.mid(header.size());
    return true;
}

QString genotypeFromCall(const QString &ref, const QStringList &alts, const QString &gt, qint64 offset)
{
    QStringList alleles{ref};
    alleles += alts;
    if (offset < 0 || offset >= ref.size())
        return QString();
    QString g;
    for (const QString &i : gt.split(QRegularExpression("[|/]"))) {
        bool ok;
        const int k = i.toInt(&ok);
        if (!ok || k < 0 || k >= alleles.size())
            return QString();
        const QString &allele = alleles[k];
        // The reference allele is the reference sequence whatever the record
        // describes; an alternate allele only names a base at this offset
        // when it is base-for-base the same length as REF.
        if (k != 0 && allele.size() != ref.size())
            return QString();
        g += allele.at(int(offset));
    }
    return g;
}

QString vcfUrl(const QString &chrom)
{
    bool ok;
    const int n = chrom.toInt(&ok);
    if (ok && n >= 1 && n <= 22)
        return QLatin1String(kBase) + QStringLiteral("ALL.chr%1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz").arg(n);
    if (chrom == "X")
        return QLatin1String(kBase) + "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz";
    if (chrom == "Y")
        return QLatin1String(kBase) + "ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz";
    if (chrom == "MT")
        return QLatin1String(kBase) + "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz";
    return QString();
}

bool build(const Options &options, const Reporter &reporter, Stats *stats, QString *error)
{
    if (isReferenceSample(options.sample))
        return buildReference(options, reporter, stats, error);

    const QString bcftools = ToolLocator::find("bcftools");
    if (bcftools.isEmpty()) {
        if (error)
            *error = QStringLiteral("bcftools is not installed. %1").arg(ToolLocator::installHint());
        return false;
    }
    const QString output = options.output.isEmpty()
        ? options.dataDir + "/genome_1000g_" + options.sample + ".txt" : options.output;

    // 1. Every rsID the analysis reads, with a GRCh37 position.
    reporter.stage(QStringLiteral("Collecting positions"), 2);
    QStringList rsids;
    if (!listRsids(options, &rsids, error))
        return false;
    const QString lookupPath = options.dataDir + "/rsid_positions_grch37.json";
    QJsonObject lookup;
    if (!positionLookup(lookupPath, rsids, reporter, &lookup, error))
        return false;

    QMap<QString, QList<Site>> byChrom;   // chrom -> sites
    for (const QString &r : rsids) {
        const QJsonObject info = lookup.value(r).toObject();
        if (info.isEmpty())
            continue;
        Site s;
        s.pos = info.value("pos").toString().toLongLong();
        s.rsid = r;
        byChrom[info.value("chrom").toString()] << s;
    }
    Stats st;
    for (auto &sites : byChrom) {
        std::sort(sites.begin(), sites.end(), [](const Site &a, const Site &b) { return a.pos < b.pos; });
        st.positions += sites.size();
    }
    reporter.log(QStringLiteral("%1 positions across %2 chromosomes for sample %3")
                     .arg(st.positions).arg(byChrom.size()).arg(options.sample));

    // 2. The sample's calls at those positions: one indexed query per
    // chromosome, several in flight at once because each spends most of its
    // time waiting on the remote index rather than transferring data.
    QTemporaryDir tmp;
    struct Query {
        QString chrom;
        QProcess *proc;
    };
    // bcftools caches each remote file's .tbi index in its working
    // directory; keep those out of the repository root.
    const QString indexDir = options.dataDir + "/1000g_index";
    QDir().mkpath(indexDir);
    QStringList chromsToQuery;
    for (auto it = byChrom.begin(); it != byChrom.end(); ++it) {
        if (vcfUrl(it.key()).isEmpty())
            reporter.log(QStringLiteral("  chr%1: no 1000 Genomes file, %2 site(s) left as no-calls")
                             .arg(it.key()).arg(it->size()));
        else
            chromsToQuery << it.key();
    }
    QSet<QString> absent;   // chromosomes the sample is not in at all

    // Chromosomes fetched by an earlier, interrupted run.
    const QString cacheDir = options.dataDir + "/test_genome_cache/" + options.sample;
    QDir().mkpath(cacheDir);
    QStringList pending;
    for (const QString &chrom : chromsToQuery) {
        QFile cf(cacheDir + "/chr" + chrom + ".tsv");
        QByteArray body;
        if (cf.open(QIODevice::ReadOnly) && chromCacheMatches(cf.readAll(), positionsOf(byChrom[chrom]), &body)) {
            if (body.startsWith("#absent"))
                absent.insert(chrom);
            else
                applyRows(body, byChrom[chrom]);
            continue;
        }
        pending << chrom;
    }
    if (pending.size() < chromsToQuery.size())
        reporter.log(QStringLiteral("  %1 of %2 chromosomes already fetched (resuming)")
                         .arg(chromsToQuery.size() - pending.size()).arg(chromsToQuery.size()));

    const int parallel = 6;
    int done = chromsToQuery.size() - pending.size();
    for (int batch = 0; batch < pending.size(); batch += parallel) {
        if (reporter.cancelled()) {
            if (error) *error = "cancelled";
            return false;
        }
        QList<Query> running;
        for (const QString &chrom : pending.mid(batch, parallel)) {
            const QString regions = tmp.filePath("regions_" + chrom + ".tsv");
            QFile rf(regions);
            rf.open(QIODevice::WriteOnly);
            for (const Site &s : byChrom[chrom])
                rf.write((chrom + '\t' + QString::number(s.pos) + '\n').toUtf8());
            rf.close();
            auto *q = new QProcess;
            q->setWorkingDirectory(indexDir);
            q->start(bcftools, {"query", "-R", regions, "-s", options.sample,
                                "-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", vcfUrl(chrom)});
            running.append({chrom, q});
        }
        reporter.stage(QStringLiteral("Querying 1000 Genomes (chromosome %1 of %2)")
                           .arg(done + 1).arg(chromsToQuery.size()),
                       5 + 80 * done / chromsToQuery.size());
        // Each query is cached the moment it completes, so a stop loses at
        // most the queries in flight.
        QList<Query> waiting = running;
        while (!waiting.isEmpty()) {
            if (reporter.cancelled()) {
                for (const Query &qr : running) {
                    if (qr.proc->state() != QProcess::NotRunning) {
                        qr.proc->kill();
                        qr.proc->waitForFinished(2000);
                    }
                    delete qr.proc;
                }
                if (error) *error = "cancelled";
                return false;
            }
            QList<Query> still;
            for (const Query &qr : waiting) {
                if (!qr.proc->waitForFinished(200)) {
                    still.append(qr);
                    continue;
                }
                const QString stderrText = QString::fromUtf8(qr.proc->readAllStandardError()).trimmed();
                QFile cf(cacheDir + "/chr" + qr.chrom + ".tsv");
                if (qr.proc->exitStatus() != QProcess::NormalExit || qr.proc->exitCode() != 0) {
                    // A sample absent from a chromosome's file -- a female and
                    // chrY -- is no data, not a failure.
                    if (stderrText.contains("not found in the header")) {
                        reporter.log(QStringLiteral("  chr%1: %2 is not in this file (no data; %3 site(s) left as no-calls)")
                                         .arg(qr.chrom, options.sample).arg(byChrom[qr.chrom].size()));
                        absent.insert(qr.chrom);
                        if (cf.open(QIODevice::WriteOnly | QIODevice::Truncate))
                            cf.write(chromCacheHeader(positionsOf(byChrom[qr.chrom])) + "#absent\n");
                        done++;
                        continue;
                    }
                    if (error)
                        *error = QStringLiteral("bcftools query failed for chr%1: %2").arg(qr.chrom, stderrText);
                    for (const Query &other : running) {
                        if (other.proc->state() != QProcess::NotRunning)
                            other.proc->kill();
                        delete other.proc;
                    }
                    return false;
                }
                const QByteArray output = qr.proc->readAllStandardOutput();
                const int hits = applyRows(output, byChrom[qr.chrom]);
                if (cf.open(QIODevice::WriteOnly | QIODevice::Truncate))
                    cf.write(chromCacheHeader(positionsOf(byChrom[qr.chrom])) + output);
                reporter.log(QStringLiteral("  chr%1: %2 of %3 sites are variant in the cohort")
                                 .arg(qr.chrom).arg(hits).arg(byChrom[qr.chrom].size()));
                done++;
                reporter.stage(QStringLiteral("Querying 1000 Genomes (%1 of %2 chromosomes)")
                                   .arg(done).arg(chromsToQuery.size()),
                               5 + 80 * done / chromsToQuery.size());
            }
            waiting = still;
        }
        for (const Query &qr : running)
            delete qr.proc;
    }

    // 3. Reference bases for the rest.
    reporter.stage(QStringLiteral("Fetching reference bases"), 88);
    QStringList regions;
    for (auto it = byChrom.begin(); it != byChrom.end(); ++it)
        for (const Site &s : *it)
            if (!s.inVcf && !vcfUrl(it.key()).isEmpty() && !absent.contains(it.key()))
                regions << QStringLiteral("%1:%2..%2:1").arg(it.key()).arg(s.pos);
    QNetworkAccessManager nam;
    QString seqError;
    QMap<QString, QString> bases;   // "chrom:pos..pos:1" -> base
    const QString samtools = ToolLocator::find("samtools");
    if (!samtools.isEmpty()) {
        // The assembly FASTA over HTTPS: no API, no rate limit.
        QStringList faidxRegions;
        for (const QString &r : regions)
            faidxRegions << RemoteFasta::region(r.section(':', 0, 0), r.section(':', 1).section("..", 0, 0).toLongLong());
        const QMap<QString, QString> got = RemoteFasta::bases(
            samtools, RemoteFasta::kGRCh37, faidxRegions, reporter, &seqError, 4,
            options.dataDir + "/test_genome_cache/fasta_GRCh37.json");
        for (int i = 0; i < regions.size(); i++)
            if (got.contains(faidxRegions[i]))
                bases[regions[i]] = got[faidxRegions[i]];
    } else {
        bases = referenceBases(nam, kSeqGRCh37, regions, &seqError);
    }
    if (seqError == "cancelled") {
        if (error) *error = seqError;
        return false;
    }
    if (!seqError.isEmpty())
        reporter.log(QStringLiteral("  Warning: %1").arg(seqError));

    // 4. Write the file.
    QFile out(output);
    if (!out.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (error)
            *error = QStringLiteral("cannot write %1").arg(output);
        return false;
    }
    out.write(QStringLiteral(
        "# 1000 Genomes Project phase 3, sample %1, GRCh37.\n"
        "# Genotypes extracted at the %2 positions this analysis reads; sites with no\n"
        "# variant in the 2,504-sample cohort are homozygous reference (base from Ensembl).\n"
        "# Open-access, consented data. Built by gh-data test-genome on %3.\n"
        "# rsid\tchromosome\tposition\tgenotype\n")
        .arg(options.sample).arg(st.positions).arg(QDateTime::currentDateTime().toString(Qt::ISODate)).toUtf8());

    QStringList chroms = byChrom.keys();
    std::sort(chroms.begin(), chroms.end(), [](const QString &a, const QString &b) { return chromOrder(a) < chromOrder(b); });
    for (const QString &chrom : chroms) {
        const bool haploid = chrom == "MT" || chrom == "Y";
        for (const Site &s : byChrom[chrom]) {
            QString genotype = "--";
            if (s.inVcf) {
                const QString g = genotypeFromCall(s.ref, s.alts, s.gt, s.offset);
                if (!g.isEmpty()) {
                    genotype = g;
                    st.fromVcf++;
                } else {
                    st.noCall++;
                }
            } else {
                const QString base = bases.value(QStringLiteral("%1:%2..%2:1").arg(chrom).arg(s.pos));
                if (base.size() == 1 && QString("ACGT").contains(base)) {
                    genotype = haploid ? base : base + base;
                    st.homRef++;
                } else {
                    st.noCall++;
                }
            }
            out.write((s.rsid + '\t' + chrom + '\t' + QString::number(s.pos) + '\t' + genotype + '\n').toUtf8());
        }
    }
    out.close();
    reporter.stage(QStringLiteral("Done"), 100);
    reporter.log(QStringLiteral("wrote %1: %2 genotypes from the 1000 Genomes calls, %3 homozygous reference, %4 no-calls")
                     .arg(output).arg(st.fromVcf).arg(st.homRef).arg(st.noCall));
    if (stats)
        *stats = st;
    return true;
}

} // namespace TestGenome
