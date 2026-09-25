// The GRCh37 reference the WGS pipeline aligns against: human_g1k_v37.fasta
// with its samtools index and minimap2 index, ~10 GB in reference/. Fetched
// and indexed from inside the app (or `gh-data reference`), so no Makefile
// or shell step is needed; each step is skipped when its output exists, so
// an interrupted setup resumes at the step it stopped in.
#pragma once

#include <QString>

#include "Reporter.h"

namespace ReferenceGenome {

extern const char *const kUrl;   // 1000 Genomes' human_g1k_v37.fasta.gz (~900 MB)

struct Status {
    QString dir, fasta, fai, mmi;
    qint64 fastaBytes = -1, mmiBytes = -1;   // -1 when absent
    bool hasFasta = false, hasFai = false, hasMmi = false;
    bool ready() const { return hasFasta && hasFai && hasMmi; }
};

Status status(const QString &root);

// Download, decompress, index. Needs samtools and minimap2 on PATH.
bool setup(const QString &root, const Reporter &reporter, QString *error);

} // namespace ReferenceGenome
