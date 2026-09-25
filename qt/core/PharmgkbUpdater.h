// ClinPGx (formerly PharmGKB): download summaryAnnotations.zip and install
// the two annotation tables under the names the analysis reads. Mirrors
// update_data.update_pharmgkb / validate_pharmgkb. The archive is read with
// the unzip binary; Qt has no zip support of its own.
#pragma once

#include <QString>
#include <QStringList>

#include "Reporter.h"

namespace PharmgkbUpdater {

extern const char *const kUrl;

// Archive member whose base name ends with `suffix`, or empty.
QString resolveMember(const QStringList &names, const QString &suffix);

bool update(const QString &dataDir, const Reporter &reporter, QString *error);

// Column checks on the installed files; records a validated timestamp on
// success. `problems` receives one line per failure.
bool validate(const QString &dataDir, const Reporter &reporter, QStringList *problems);

} // namespace PharmgkbUpdater
