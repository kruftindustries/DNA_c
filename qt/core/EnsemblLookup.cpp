#include "EnsemblLookup.h"

#include <QEventLoop>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QThread>
#include <QTimer>

#include "DataVersions.h"

namespace EnsemblLookup {

const char *const kEndpoint = "https://grch37.rest.ensembl.org/variation/homo_sapiens";
const char *const kEndpointGRCh38 = "https://rest.ensembl.org/variation/homo_sapiens";

QJsonObject loadLookup(const QString &path)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        return {};
    return QJsonDocument::fromJson(f.readAll()).object();
}

bool saveLookup(const QString &path, const QJsonObject &lookup)
{
    // Same two-space layout as json.dump(indent=2); the C loader reads the
    // shape, not the whitespace, but keeping it stable keeps diffs small.
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate))
        return false;
    f.write(DataVersions::indentedJson(lookup));
    return true;
}

QJsonObject primaryMapping(const QJsonObject &info)
{
    static const QSet<QString> valid = [] {
        QSet<QString> s;
        for (int i = 1; i <= 22; ++i)
            s.insert(QString::number(i));
        s.insert("X"); s.insert("Y"); s.insert("MT");
        return s;
    }();
    for (const QJsonValue &m : info.value("mappings").toArray()) {
        QJsonObject mo = m.toObject();
        QString chrom = mo.value("seq_region_name").toString();
        chrom.replace("chr", "");
        if (valid.contains(chrom)) {
            mo["chrom"] = chrom;
            return mo;
        }
    }
    return QJsonObject();
}

int mergeBatch(const QJsonObject &results, const QSet<QString> &requested, QJsonObject &lookup)
{
    int gained = 0;
    for (auto it = results.begin(); it != results.end(); ++it) {
        const QJsonObject info = it.value().toObject();
        const QJsonObject mapping = primaryMapping(info);
        if (mapping.isEmpty())
            continue;
        QJsonObject position;
        position["chrom"] = mapping.value("chrom").toString();
        position["pos"] = QString::number(qint64(mapping.value("start").toDouble()));
        QSet<QString> names{it.key(), info.value("name").toString()};
        for (const QJsonValue &s : info.value("synonyms").toArray())
            names.insert(s.toString());
        for (const QString &name : names) {
            if (requested.contains(name)) {
                lookup[name] = position;
                gained++;
            }
        }
    }
    return gained;
}

bool postVariationBatch(QNetworkAccessManager &nam, const char *endpoint, const QStringList &ids,
                        QJsonObject *results, QString *error)
{
    QNetworkRequest req{QUrl(QLatin1String(endpoint))};
    req.setHeader(QNetworkRequest::ContentTypeHeader, "application/json");
    req.setRawHeader("Accept", "application/json");
    QJsonArray arr;
    for (const QString &id : ids)
        arr.append(id);
    QNetworkReply *reply = nam.post(req, QJsonDocument(QJsonObject{{"ids", arr}}).toJson(QJsonDocument::Compact));

    QEventLoop loop;
    QTimer timeout;
    timeout.setSingleShot(true);
    QObject::connect(&timeout, &QTimer::timeout, reply, &QNetworkReply::abort);
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    timeout.start(60000);
    loop.exec();

    bool ok = reply->error() == QNetworkReply::NoError;
    if (ok) {
        QJsonParseError perr;
        const QJsonDocument doc = QJsonDocument::fromJson(reply->readAll(), &perr);
        ok = perr.error == QJsonParseError::NoError && doc.isObject();
        if (ok)
            *results = doc.object();
        else if (error)
            *error = perr.errorString();
    } else if (error) {
        *error = reply->errorString();
    }
    reply->deleteLater();
    return ok;
}

bool getVariation(QNetworkAccessManager &nam, const char *endpoint, const QString &id,
                  QJsonObject *result, QString *error)
{
    QNetworkRequest req{QUrl(QLatin1String(endpoint) + "/" + id + "?content-type=application/json")};
    req.setRawHeader("Accept", "application/json");
    QNetworkReply *reply = nam.get(req);
    QEventLoop loop;
    QTimer timeout;
    timeout.setSingleShot(true);
    QObject::connect(&timeout, &QTimer::timeout, reply, &QNetworkReply::abort);
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    timeout.start(60000);
    loop.exec();
    bool ok = reply->error() == QNetworkReply::NoError;
    if (ok) {
        QJsonParseError perr;
        const QJsonDocument doc = QJsonDocument::fromJson(reply->readAll(), &perr);
        ok = perr.error == QJsonParseError::NoError && doc.isObject();
        if (ok)
            *result = doc.object();
        else if (error)
            *error = perr.errorString();
    } else if (error) {
        *error = reply->errorString();
    }
    reply->deleteLater();
    return ok;
}

bool update(const QString &lookupPath, const QStringList &rsids, const Reporter &reporter,
            QStringList *unresolved)
{
    reporter.log(QStringLiteral("Building the rsID position lookup"));
    QJsonObject lookup = loadLookup(lookupPath);
    QStringList missing;
    for (const QString &r : rsids)
        if (!lookup.contains(r))
            missing << r;

    if (missing.isEmpty()) {
        reporter.log(QStringLiteral("  Lookup already complete (%1 rsIDs)").arg(lookup.size()));
    } else {
        reporter.log(QStringLiteral("  %1 rsIDs cached, %2 to query from Ensembl...")
                         .arg(lookup.size()).arg(missing.size()));
        QNetworkAccessManager nam;
        const int batchSize = 200;
        for (int i = 0; i < missing.size(); i += batchSize) {
            if (reporter.cancelled())
                return false;
            const QStringList batch = missing.mid(i, batchSize);
            QJsonObject results;
            QString err;
            bool ok = false;
            for (int attempt = 0; attempt < 3 && !ok; ++attempt) {
                ok = postVariationBatch(nam, kEndpoint, batch, &results, &err);
                if (!ok && attempt < 2)
                    QThread::msleep(1000u << (attempt + 1));
            }
            if (!ok) {
                reporter.log(QStringLiteral("  Warning: API batch %1 failed: %2")
                                 .arg(i / batchSize + 1).arg(err));
                continue;
            }
            QSet<QString> requested(batch.begin(), batch.end());
            mergeBatch(results, requested, lookup);
            QThread::msleep(200);
        }
    }

    QStringList left;
    for (const QString &r : rsids)
        if (!lookup.contains(r))
            left << r;
    if (!left.isEmpty())
        reporter.log(QStringLiteral("  WARNING: %1 rsID(s) have no GRCh37 position and will not "
                                    "be called: %2").arg(left.size()).arg(left.join(", ")));
    if (unresolved)
        *unresolved = left;

    if (!saveLookup(lookupPath, lookup))
        return false;
    reporter.log(QStringLiteral("  Saved %1 rsID positions (%2 requested)").arg(lookup.size()).arg(rsids.size()));
    return true;
}

} // namespace EnsemblLookup
