// Read members of a zip archive in-process (vendored miniz), so fetching the
// ClinPGx bundle needs no unzip binary.
#pragma once

#include <QString>
#include <QStringList>

namespace ZipArchive {

// Every member name, or an empty list with *error set.
QStringList names(const QString &archive, QString *error);

// Extract one member (exact name) to a file.
bool extract(const QString &archive, const QString &member, const QString &dest, QString *error);

} // namespace ZipArchive
