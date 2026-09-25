// Bases from a bgzipped, indexed FASTA over HTTPS, through `samtools faidx`,
// which fetches the .fai/.gzi indexes and then only the blocks it needs.
// Ensembl publishes every assembly this way (unmasked.fa.bgz beside its
// .fai and .gzi), so a reference base at any position is a range request
// rather than an API call. Each region costs a round trip, so the regions
// are split across several samtools processes.
#pragma once

#include <QByteArray>
#include <QMap>
#include <QString>
#include <QStringList>

#include "Reporter.h"

namespace RemoteFasta {

extern const char *const kGRCh37;   // Ensembl GRCh37.p13 unmasked.fa.bgz
extern const char *const kGRCh38;   // Ensembl GRCh38.p14 unmasked.fa.bgz
extern const char *const kChainGRCh37ToGRCh38;

// "chrom:start-end" (1-based, inclusive) -> upper-case sequence. Regions
// samtools could not fetch are absent from the result; *error is set only
// when a process fails outright, or to "cancelled". With `cachePath`, a
// JSON map of regions already fetched is read first and extended after
// every batch, so an interrupted run resumes where it stopped.
QMap<QString, QString> bases(const QString &samtools, const QString &url, const QStringList &regions,
                             const Reporter &reporter, QString *error, int parallel = 4,
                             const QString &cachePath = QString());

// Parse `samtools faidx` output: ">region" lines followed by sequence lines.
QMap<QString, QString> parseFasta(const QByteArray &text);

QString region(const QString &chrom, qint64 pos, int flank = 0);

} // namespace RemoteFasta
