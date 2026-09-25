// Finds the external tools the WGS pipeline shells out to (minimap2,
// samtools, bcftools, fastp) on PATH, and the analysis engine.
#pragma once

#include <QProcessEnvironment>
#include <QString>
#include <QStringList>

namespace ToolLocator {

// Absolute path, or empty when not found anywhere.
QString find(const QString &tool);

// Directories to prepend to PATH so every tool that was found is on it.
QStringList extraPathDirs(const QStringList &tools);

// The process environment with those directories prepended.
QProcessEnvironment environmentFor(const QStringList &tools);

// The tools the FASTQ pipeline needs, required ones first.
QStringList wgsTools();
QStringList requiredWgsTools();   // without fastp, which is optional

// Which of `tools` cannot be found at all.
QStringList missing(const QStringList &tools);

// The one-line install command for everything the pipeline needs, for
// this platform's package manager.
QString installHint();

// The analysis engine for a checkout root: beside this executable, at the
// root (where the top-level build puts it), in c/bin (the c/Makefile build),
// then on PATH. Empty if none.
QString analysisBinary(const QString &root);

} // namespace ToolLocator
