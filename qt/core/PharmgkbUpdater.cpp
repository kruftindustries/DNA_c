#include "PharmgkbUpdater.h"
#include "ZipArchive.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonObject>
#include <QLocale>
#include <QTemporaryDir>

#include "DataVersions.h"
#include "Downloader.h"
#include "TextTable.h"

namespace PharmgkbUpdater {

const char *const kUrl = "https://api.clinpgx.org/v1/download/file/data/summaryAnnotations.zip";

namespace {

const char *const kIdColumns[] = {"Clinical Annotation ID", "Summary Annotation ID"};

struct Member {
    const char *suffix;   // matched against the member's base name
    const char *outName;
};
const Member kMembers[] = {
    {"annotations.tsv", "clinical_annotations.tsv"},
    {"ann_alleles.tsv", "clinical_ann_alleles.tsv"},
};

} // namespace

QString resolveMember(const QStringList &names, const QString &suffix)
{
    for (const QString &n : names) {
        const QString base = n.section('/', -1);
        if (base.endsWith(suffix))
            return n;
    }
    return QString();
}

bool update(const QString &dataDir, const Reporter &reporter, QString *error)
{
    QDir().mkpath(dataDir);
    reporter.log(QStringLiteral(">>> Downloading ClinPGx annotations from %1").arg(QLatin1String(kUrl)));

    QTemporaryDir tmp;
    const QString archive = tmp.filePath("summaryAnnotations.zip");
    if (!Downloader::download(QUrl(QLatin1String(kUrl)), archive, reporter, error))
        return false;

    reporter.log(QStringLiteral(">>> Extracting annotation tables"));
    const QStringList names = ZipArchive::names(archive, error);
    if (names.isEmpty())
        return false;

    QString created = resolveMember(names, "CREATED.txt");
    if (created.isEmpty())
        for (const QString &n : names)
            if (n.startsWith("CREATED_")) { created = n; break; }
    QString release = QStringLiteral("unknown");
    if (!created.isEmpty())
        release = QFileInfo(created).completeBaseName().replace("CREATED_", "");

    QJsonArray written;
    for (const Member &m : kMembers) {
        const QString member = resolveMember(names, QLatin1String(m.suffix));
        if (member.isEmpty()) {
            if (error)
                *error = QStringLiteral("no member ending in '%1' in %2; archive contents: %3")
                             .arg(QLatin1String(m.suffix), QLatin1String(kUrl), names.join(", "));
            return false;
        }
        const QString outPath = dataDir + '/' + QLatin1String(m.outName);
        if (!ZipArchive::extract(archive, member, outPath, error))
            return false;
        written.append(QLatin1String(m.outName));
        reporter.log(QStringLiteral("    %1 -> %2").arg(member, outPath));
    }

    QJsonObject versions = DataVersions::load(dataDir);
    QJsonObject pg;
    pg["updated"] = DataVersions::nowIso();
    pg["source"] = QLatin1String(kUrl);
    pg["release"] = release;
    pg["files"] = written;
    versions["pharmgkb"] = pg;
    DataVersions::save(dataDir, versions);
    reporter.log(QStringLiteral("    ClinPGx release: %1").arg(release));

    QStringList problems;
    const bool ok = validate(dataDir, reporter, &problems);
    if (!ok && error)
        *error = problems.join("; ");
    return ok;
}

bool validate(const QString &dataDir, const Reporter &reporter, QStringList *problems)
{
    struct Check { const char *file; QStringList expected; };
    const Check checks[] = {
        {"clinical_annotations.tsv", {"Variant/Haplotypes", "Gene", "Drug(s)"}},
        {"clinical_ann_alleles.tsv", {"Genotype/Allele", "Annotation Text"}},
    };
    QStringList errors;
    for (const Check &c : checks) {
        const QString path = dataDir + '/' + QLatin1String(c.file);
        QFile f(path);
        if (!f.open(QIODevice::ReadOnly)) {
            errors << QStringLiteral("Missing: %1").arg(path);
            reporter.log(QStringLiteral("    ERROR: %1 not found").arg(path));
            continue;
        }
        TsvTable table(&f);
        QStringList missing;
        for (const QString &h : c.expected)
            if (!table.hasColumn(h))
                missing << h;
        bool hasId = false;
        for (const char *id : kIdColumns)
            hasId = hasId || table.hasColumn(QLatin1String(id));
        if (!hasId)
            missing << QStringLiteral("Clinical Annotation ID or Summary Annotation ID");
        if (!missing.isEmpty()) {
            errors << QStringLiteral("Missing columns in %1: %2").arg(QLatin1String(c.file), missing.join(", "));
            reporter.log(QStringLiteral("    WARNING: Missing columns: %1").arg(missing.join(", ")));
            continue;
        }
        qint64 rows = 0;
        while (table.next())
            rows++;
        reporter.log(QStringLiteral("    %1: OK (%L2 rows)").arg(QLatin1String(c.file)).arg(rows));
    }

    if (errors.isEmpty()) {
        QJsonObject versions = DataVersions::load(dataDir);
        QJsonObject pg = versions.value("pharmgkb").toObject();
        pg["validated"] = DataVersions::nowIso();
        if (!pg.contains("files"))
            pg["files"] = QJsonArray{"clinical_annotations.tsv", "clinical_ann_alleles.tsv"};
        versions["pharmgkb"] = pg;
        DataVersions::save(dataDir, versions);
        reporter.log(QStringLiteral("    ClinPGx validation passed. Metadata updated."));
    } else {
        reporter.log(QStringLiteral("    ClinPGx validation failed: %1 error(s)").arg(errors.size()));
    }
    if (problems)
        *problems = errors;
    return errors.isEmpty();
}

} // namespace PharmgkbUpdater
