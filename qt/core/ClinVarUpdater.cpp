#include "ClinVarUpdater.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonObject>
#include <QLocale>

#include "DataVersions.h"
#include "Downloader.h"
#include "GzipStream.h"
#include "TextTable.h"

namespace ClinVarUpdater {

const char *const kUrl =
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz";

int goldStars(const QString &reviewStatus)
{
    // Same table, same order: the first key contained in the status wins.
    static const struct { const char *key; int stars; } table[] = {
        {"practice guideline", 4},
        {"reviewed by expert panel", 4},
        {"criteria provided, multiple submitters, no conflicts", 3},
        {"criteria provided, multiple submitters, conflicting interpretations", 2},
        {"criteria provided, conflicting interpretations", 2},
        {"criteria provided, single submitter", 1},
        {"no assertion for the individual variant", 0},
        {"no assertion criteria provided", 0},
        {"no assertion provided", 0},
    };
    const QString low = reviewStatus.trimmed().toLower();
    for (const auto &e : table)
        if (low.contains(QLatin1String(e.key)))
            return e.stars;
    return 0;
}

bool process(const QString &gzPath, const QString &tsvPath, Stats *stats,
             const Reporter &reporter, QString *error)
{
    GzipStream in(gzPath);
    if (!in.open(QIODevice::ReadOnly)) {
        if (error)
            *error = in.errorText();
        return false;
    }
    QFile out(tsvPath);
    if (!out.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (error)
            *error = QStringLiteral("cannot write %1: %2").arg(tsvPath, out.errorString());
        return false;
    }

    TsvTable table(&in);
    if (!table.ok()) {
        if (error)
            *error = QStringLiteral("empty or unreadable archive");
        return false;
    }

    static const QStringList outputFields = {
        "chrom", "pos", "ref", "alt", "clinical_significance",
        "review_status", "gold_stars", "all_traits", "symbol",
        "inheritance_modes", "hgvs_p", "hgvs_c",
        "molecular_consequence", "xrefs",
    };
    out.write(tsvRow(outputFields));

    const int cAssembly = table.column("Assembly");
    const int cChrom = table.column("Chromosome");
    const int cStart = table.column("Start");
    const int cRef = table.column("ReferenceAlleleVCF");
    const int cAlt = table.column("AlternateAlleleVCF");
    const int cReview = table.column("ReviewStatus");
    const int cSig = table.column("ClinicalSignificance");
    const int cPheno = table.column("PhenotypeList");
    const int cGene = table.column("GeneSymbol");
    const int cOrigin = table.column("OriginSimple");
    // variant_summary has no ProteinChange or HGVS(c) column; the Python
    // writes "" for both when the columns are absent, and they always are.
    const int cHgvsP = table.column("ProteinChange");
    const int cHgvsC = table.column("HGVS(c)");
    const int cCons = table.column("MolecularConsequence");
    const int cXrefs = table.column("OtherIDs");

    Stats local;
    const qint64 inputSize = QFileInfo(gzPath).size();
    while (table.next()) {
        local.rowsProcessed++;
        if ((local.rowsProcessed & 0xFFFF) == 0) {
            if (reporter.cancelled()) {
                if (error)
                    *error = QStringLiteral("cancelled");
                return false;
            }
            reporter.stage(QStringLiteral("Processing ClinVar (%L1 rows)").arg(local.rowsProcessed),
                           inputSize > 0 ? int(qMin<qint64>(99, table.bytesRead() * 30 / (inputSize * 10))) : 0);
        }
        if (table.field(cAssembly) != QLatin1String("GRCh37"))
            continue;
        const QString chrom = table.field(cChrom), pos = table.field(cStart);
        if (chrom.isEmpty() || pos.isEmpty())
            continue;
        const QString ref = table.field(cRef), alt = table.field(cAlt);
        if (ref.isEmpty() || alt.isEmpty())
            continue;
        const QString review = table.field(cReview);
        out.write(tsvRow({
            chrom, pos, ref, alt, table.field(cSig), review,
            QString::number(goldStars(review)), table.field(cPheno), table.field(cGene),
            // "inheritance_modes" is filled from OriginSimple, as the Python
            // does: variant_summary carries no inheritance column at all.
            table.field(cOrigin),
            cHgvsP >= 0 ? table.field(cHgvsP) : QString(),
            cHgvsC >= 0 ? table.field(cHgvsC) : QString(),
            table.field(cCons), table.field(cXrefs),
        }));
        local.written++;
    }
    out.close();
    if (!in.errorText().isEmpty()) {
        if (error)
            *error = in.errorText();
        return false;
    }
    if (stats)
        *stats = local;
    return true;
}

bool update(const QString &dataDir, const Reporter &reporter, QString *error)
{
    QDir().mkpath(dataDir);
    const QString gz = dataDir + "/variant_summary.txt.gz";
    const QString tsv = dataDir + "/clinvar_alleles.tsv";

    reporter.log(QStringLiteral(">>> Downloading ClinVar from %1").arg(QLatin1String(kUrl)));
    if (!Downloader::download(QUrl(QLatin1String(kUrl)), gz, reporter, error))
        return false;

    reporter.log(QStringLiteral(">>> Processing ClinVar data (filtering to GRCh37 SNPs)..."));
    Stats stats;
    if (!process(gz, tsv, &stats, reporter, error)) {
        QFile::remove(gz);
        return false;
    }
    QFile::remove(gz);

    QJsonObject versions = DataVersions::load(dataDir);
    QJsonObject cv;
    cv["updated"] = DataVersions::nowIso();
    cv["source"] = QLatin1String(kUrl);
    cv["total_variants_processed"] = double(stats.rowsProcessed);
    cv["grch37_variants_written"] = double(stats.written);
    cv["file"] = QStringLiteral("clinvar_alleles.tsv");
    versions["clinvar"] = cv;
    DataVersions::save(dataDir, versions);

    reporter.log(QStringLiteral("    Total ClinVar rows processed: %L1").arg(stats.rowsProcessed));
    reporter.log(QStringLiteral("    GRCh37 variants written: %L1").arg(stats.written));
    reporter.log(QStringLiteral("    Output: %1").arg(tsv));
    reporter.stage(QStringLiteral("ClinVar updated"), 100);
    return true;
}

} // namespace ClinVarUpdater
