// Finds the external tools the WGS pipeline shells out to (minimap2,
// samtools, bcftools, fastp), and the analysis engine.
//
// Search order: the `tools` folder beside this executable (the Windows
// release ships the three required tools there), then PATH, then the
// package-manager prefixes a Finder- or Explorer-launched app does not see
// on PATH (/opt/homebrew/bin, /usr/local/bin, /opt/local/bin, ~/.local/bin).
#pragma once

#include <QProcessEnvironment>
#include <QString>
#include <QStringList>

namespace ToolLocator {

// Absolute path, or empty when not found anywhere.
QString find(const QString &tool);

// The bundled tools folder beside this executable (may not exist).
QString bundledDir();

// Directories to prepend to PATH so every tool that was found is on it.
QStringList extraPathDirs(const QStringList &tools);

// The process environment with those directories prepended.
QProcessEnvironment environmentFor(const QStringList &tools);

// The tools the FASTQ pipeline needs, required ones first.
QStringList wgsTools();
QStringList requiredWgsTools();   // without fastp, which is optional

// Which of `tools` cannot be found at all.
QStringList missing(const QStringList &tools);

// One sentence telling the user how to get the missing tools on this
// platform: the package-manager command on Linux and macOS; on Windows,
// where they ship in the release zip, how to restore them.
QString installHint();

// The analysis engine for a checkout root: beside this executable, at the
// root (where the top-level build puts it), in c/bin (the c/Makefile build),
// then on PATH. Empty if none.
QString analysisBinary(const QString &root);

} // namespace ToolLocator
