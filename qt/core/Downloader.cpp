#include "Downloader.h"

#include <QEventLoop>
#include <QFile>
#include <QFileInfo>
#include <QLocale>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QTimer>

namespace Downloader {

const char *const kUserAgent =
    "genetic-health-pipeline/1.0 (+https://github.com/unbalancedparentheses/DNA)";

bool download(const QUrl &url, const QString &dest, const Reporter &reporter, QString *error)
{
    QFile out(dest);
    if (!out.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (error)
            *error = QStringLiteral("cannot write %1: %2").arg(dest, out.errorString());
        return false;
    }

    QNetworkAccessManager nam;
    QNetworkRequest req(url);
    req.setHeader(QNetworkRequest::UserAgentHeader, QString::fromLatin1(kUserAgent));
    req.setAttribute(QNetworkRequest::RedirectPolicyAttribute,
                     QNetworkRequest::NoLessSafeRedirectPolicy);
    QNetworkReply *reply = nam.get(req);

    QEventLoop loop;
    bool cancelled = false;
    qint64 received = 0, total = -1;
    QObject::connect(reply, &QNetworkReply::readyRead, &loop, [&] {
        out.write(reply->readAll());
    });
    QObject::connect(reply, &QNetworkReply::downloadProgress, &loop, [&](qint64 r, qint64 t) {
        received = r;
        total = t;
        const QString size = t > 0
            ? QStringLiteral("%1 of %2").arg(QLocale().formattedDataSize(r), QLocale().formattedDataSize(t))
            : QLocale().formattedDataSize(r);
        reporter.stage(QStringLiteral("Downloading %1").arg(QFileInfo(dest).fileName()),
                       t > 0 ? int(r * 100 / t) : 0);
        Q_UNUSED(size);
    });
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    QTimer poll;
    QObject::connect(&poll, &QTimer::timeout, &loop, [&] {
        if (reporter.cancelled()) {
            cancelled = true;
            reply->abort();
        }
    });
    poll.start(250);
    loop.exec();

    out.write(reply->readAll());
    out.close();
    const bool ok = !cancelled && reply->error() == QNetworkReply::NoError;
    if (!ok && error)
        *error = cancelled ? QStringLiteral("cancelled") : reply->errorString();
    reply->deleteLater();
    if (ok)
        reporter.log(QStringLiteral("    Downloaded %1 (%2)")
                         .arg(dest, QLocale().formattedDataSize(received)));
    else
        QFile::remove(dest);
    return ok;
}

} // namespace Downloader
