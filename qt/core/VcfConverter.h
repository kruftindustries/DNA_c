// VCF -> the four-column genotype file the analysis reads. Mirrors
// wgs_pipeline.step_convert: single-nucleotide variants only, the rsID from
// the VCF's ID column, else from the position lookup, else a positional
// identifier the loader can resolve later.
#pragma once

#include <QString>

#include "Reporter.h"

namespace VcfConverter {

// Returns the number of genotypes written, or -1 with *error.
qint64 convert(const QString &vcfPath, const QString &genomeOut, const QString &lookupJson,
               const Reporter &reporter, QString *error);

} // namespace VcfConverter
