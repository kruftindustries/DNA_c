// ClinVar: download variant_summary.txt.gz and reduce it to the GRCh37
// clinvar_alleles.tsv the analysis reads. Mirrors update_data.update_clinvar,
// including the columns it writes and how it maps review status to stars;
// tests/test_qt_parity.py checks the two processors produce identical bytes.
#pragma once

#include <QString>

#include "Reporter.h"

namespace ClinVarUpdater {

extern const char *const kUrl;

struct Stats {
    qint64 rowsProcessed = 0;
    qint64 written = 0;
};

int goldStars(const QString &reviewStatus);

// The filtering step alone: gz in, tsv out. False with *error on failure.
bool process(const QString &gzPath, const QString &tsvPath, Stats *stats,
             const Reporter &reporter, QString *error);

// Download, process, delete the archive, record the release in
// data_versions.json.
bool update(const QString &dataDir, const Reporter &reporter, QString *error);

} // namespace ClinVarUpdater
