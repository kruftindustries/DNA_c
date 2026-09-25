#include "DataVersions.h"

#include <QDateTime>
#include <QFile>
#include <QJsonDocument>

namespace DataVersions {

QJsonObject load(const QString &dataDir)
{
    QFile f(dataDir + "/data_versions.json");
    if (!f.open(QIODevice::ReadOnly))
        return {};
    return QJsonDocument::fromJson(f.readAll()).object();
}

bool save(const QString &dataDir, const QJsonObject &versions)
{
    QFile f(dataDir + "/data_versions.json");
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate))
        return false;
    f.write(indentedJson(versions));
    return true;
}

QByteArray indentedJson(const QJsonObject &object)
{
    // Python writes indent=2 in insertion order; QJsonDocument writes four
    // spaces with sorted keys. Key order is not something any reader
    // depends on; the indentation is kept at two so the files stay the shape
    // they were.
    QByteArray text = QJsonDocument(object).toJson(QJsonDocument::Indented);
    if (text.endsWith('\n'))
        text.chop(1);
    QList<QByteArray> lines;
    for (const QByteArray &line : text.split('\n')) {
        int spaces = 0;
        while (spaces < line.size() && line[spaces] == ' ')
            spaces++;
        lines << QByteArray(spaces / 2, ' ') + line.mid(spaces);
    }
    return lines.join('\n');
}

QString nowIso()
{
    return QDateTime::currentDateTime().toString(QStringLiteral("yyyy-MM-ddTHH:mm:ss.zzz")) + "000";
}

} // namespace DataVersions
