// The BED of positions worth calling: every ClinVar row plus every rsID with
// a known position. Calling only there is what makes a WGS run take minutes
// rather than hours. Mirrors wgs_pipeline._build_target_regions.
#pragma once

#include <QString>

#include "Reporter.h"

namespace TargetRegions {

// Writes `bedPath`; returns the number of positions, or -1 with *error.
qint64 build(const QString &clinvarTsv, const QString &lookupJson, const QString &bedPath,
             const Reporter &reporter, QString *error);

} // namespace TargetRegions
