// End-to-end check of the WGS toolchain against a truth set, natively.
//
// The port of scripts/validate_wgs.py: a single chromosome (chr22, ~10 MB)
// instead of the 10 GB reference, reads simulated from two haplotypes whose
// genotypes were chosen here, then the real pipeline steps -- align, call,
// convert -- and a comparison of what came out with what went in.
#pragma once

#include <QString>
#include <QStringList>

#include "Reporter.h"

namespace WgsValidator {

struct Options {
    QString root;
    int clinvarTargets = 150;
    int threads = 0;
    unsigned seed = 1;
};

QString referenceDir(const QString &root);
bool isSetUp(const QString &root);

// Fetch and index chr22. Idempotent.
bool setup(const QString &root, const Reporter &reporter, QString *error);

// Simulate, run, compare. `problems` lists each locus that disagreed.
bool run(const Options &options, const Reporter &reporter, QStringList *problems, QString *error);

} // namespace WgsValidator
