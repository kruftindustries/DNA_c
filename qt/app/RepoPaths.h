// Where the analysis lives: the C binary and the data directory. Resolved
// once from the application's location and overridable through settings, so
// a packaged build and a build tree both work.
#pragma once

#include <QString>

namespace RepoPaths {

// The checkout root: the nearest ancestor of the executable (or of the
// working directory) containing c/src/main.c. Empty if none.
QString root();

QString dataDir();          // <root>/data unless overridden
QString analysisBinary();   // the engine beside the app, at the root or in c/bin, unless overridden

void setDataDir(const QString &path);
void setAnalysisBinary(const QString &path);

} // namespace RepoPaths
