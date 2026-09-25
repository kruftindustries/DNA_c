// data/data_versions.json: what release of each annotation source is
// installed. Same file the Python updater writes and the GUI reads.
#pragma once

#include <QJsonObject>
#include <QString>

namespace DataVersions {

QJsonObject load(const QString &dataDir);
bool save(const QString &dataDir, const QJsonObject &versions);

// datetime.now().isoformat() -- local time with six fractional digits.
QString nowIso();

// QJsonDocument's four-space output re-indented to json.dump(indent=2)'s two,
// without the trailing newline json.dump does not write.
QByteArray indentedJson(const QJsonObject &object);

} // namespace DataVersions
