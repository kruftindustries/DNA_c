#include "RemoteFasta.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QList>
#include <QSaveFile>
#include <QPair>
#include <QProcess>
#include <QThread>

namespace RemoteFasta {

const char *const kGRCh37 =
    "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/14/ensembl/2013_09/genome/unmasked.fa.bgz";
const char *const kGRCh38 =
    "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome/unmasked.fa.bgz";
const char *const kChainGRCh37ToGRCh38 =
    "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome/assembly_mapping/GRCh37_to_GRCh38.chain.gz";

QString region(const QString &chrom, qint64 pos, int flank)
{
    return QStringLiteral("%1:%2-%3").arg(chrom).arg(pos - flank).arg(pos + flank);
}

QMap<QString, QString> parseFasta(const QByteArray &text)
{
    QMap<QString, QString> out;
    QString name;
    for (const QByteArray &raw : text.split('\n')) {
        const QByteArray line = raw.trimmed();
        if (line.isEmpty())
            continue;
        if (line.startsWith('>')) {
            name = QString::fromLatin1(line.mid(1));
            out[name] = QString();
        } else if (!name.isEmpty()) {
            out[name] += QString::fromLatin1(line).toUpper();
        }
    }
    return out;
}

namespace {

QMap<QString, QString> loadCache(const QString &path)
{
    QMap<QString, QString> out;
    QFile f(path);
    if (path.isEmpty() || !f.open(QIODevice::ReadOnly))
        return out;
    const QJsonObject o = QJsonDocument::fromJson(f.readAll()).object();
    for (auto it = o.begin(); it != o.end(); ++it)
        out[it.key()] = it.value().toString();
    return out;
}

void saveCache(const QString &path, const QMap<QString, QString> &cache)
{
    if (path.isEmpty())
        return;
    QDir().mkpath(QFileInfo(path).path());
    QJsonObject o;
    for (auto it = cache.begin(); it != cache.end(); ++it)
        o[it.key()] = it.value();
    QSaveFile f(path);
    if (f.open(QIODevice::WriteOnly)) {
        f.write(QJsonDocument(o).toJson(QJsonDocument::Compact));
        f.commit();
    }
}

// Wait for every process, polling so a cancellation kills them promptly.
// Returns false when cancelled.
bool waitAll(const QList<QProcess *> &procs, const Reporter &reporter)
{
    for (QProcess *p : procs) {
        while (!p->waitForFinished(500)) {
            if (reporter.cancelled()) {
                for (QProcess *q : procs)
                    if (q->state() != QProcess::NotRunning) {
                        q->kill();
                        q->waitForFinished(2000);
                    }
                return false;
            }
        }
    }
    return true;
}

} // namespace

QMap<QString, QString> bases(const QString &samtools, const QString &url, const QStringList &regions,
                             const Reporter &reporter, QString *error, int parallel, const QString &cachePath)
{
    QMap<QString, QString> cache = loadCache(cachePath);
    QMap<QString, QString> out;
    QStringList missing;
    for (const QString &r : regions) {
        if (cache.contains(r))
            out[r] = cache.value(r);
        else
            missing << r;
    }
    if (!cachePath.isEmpty() && missing.size() < regions.size())
        reporter.log(QStringLiteral("  %1 of %2 regions already fetched (resuming)")
                         .arg(regions.size() - missing.size()).arg(regions.size()));

    const int chunk = 60;
    QList<QStringList> chunks;
    for (int i = 0; i < missing.size(); i += chunk)
        chunks.append(missing.mid(i, chunk));

    for (int batch = 0; batch < chunks.size(); batch += parallel) {
        if (reporter.cancelled()) {
            if (error) *error = "cancelled";
            return out;
        }
        // A chunk that fails (the server refusing one of the index
        // downloads is the usual way) is retried on its own before it
        // counts as an error.
        QList<QStringList> pending = chunks.mid(batch, parallel);
        for (int attempt = 0; attempt < 4 && !pending.isEmpty(); attempt++) {
            if (attempt)
                QThread::msleep(2000 << (attempt - 1));
            QList<QPair<QStringList, QProcess *>> running;
            QList<QProcess *> procs;
            for (const QStringList &c : pending) {
                auto *p = new QProcess;
                p->start(samtools, QStringList{"faidx", url} + c);
                running.append({c, p});
                procs.append(p);
            }
            const bool finished = waitAll(procs, reporter);
            pending.clear();
            QString lastError;
            for (auto &r : running) {
                QProcess *p = r.second;
                if (finished && p->exitStatus() == QProcess::NormalExit && p->exitCode() == 0) {
                    const QMap<QString, QString> got = parseFasta(p->readAllStandardOutput());
                    for (auto it = got.begin(); it != got.end(); ++it) {
                        out[it.key()] = it.value();
                        cache[it.key()] = it.value();
                    }
                } else {
                    lastError = QString::fromUtf8(p->readAllStandardError()).trimmed();
                    pending.append(r.first);
                }
                delete p;
            }
            saveCache(cachePath, cache);
            if (!finished) {
                if (error) *error = "cancelled";
                return out;
            }
            if (!pending.isEmpty()) {
                reporter.log(QStringLiteral("  %1 chunk(s) failed (%2); retrying").arg(pending.size()).arg(lastError));
                if (attempt == 3 && error)
                    *error = QStringLiteral("samtools faidx failed: %1").arg(lastError);
            }
        }
        reporter.log(QStringLiteral("  %1 of %2 regions fetched").arg(qMin((batch + parallel) * chunk, missing.size())).arg(missing.size()));
    }
    return out;
}

} // namespace RemoteFasta
