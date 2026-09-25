// Where the analysis lives: the C binary and the data directory. Resolved
// once from the application's location and overridable through settings, so
// a packaged build and a build tree both work.
#pragma once

#include <QString>

namespace RepoPaths {

// The working root: the checkout around the executable (or the working
// directory), else the per-user application data directory -- see
// core/Locations.h. Overridable in Settings.
QString root();

QString dataDir();          // <root>/data unless overridden
QString analysisBinary();   // the engine beside the app, at the root or in c/bin, unless overridden

void setDataDir(const QString &path);
void setAnalysisBinary(const QString &path);

} // namespace RepoPaths
