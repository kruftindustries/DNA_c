// A real, public test genome: one 1000 Genomes Project sample, extracted at
// the positions this analysis reads and written in 23andMe layout.
//
// The phase 3 VCFs are indexed, so bcftools pulls just the needed rows over
// HTTPS -- a few megabytes of range requests rather than the 50 GB of files.
// Sites the cohort has no variant at are homozygous reference; the reference
// base for those comes from Ensembl, so the file carries a call at every
// position the way a genotyping array would. NA12878 (HG001, CEU) is the
// default: the most-studied human genome, openly consented.
//
// The sample may also be a reference assembly, "GRCh37" or "GRCh38": the
// build's own base at every position, homozygous. That is a null test of the
// curated tables -- each variant the report then claims is a statement about
// the reference genome, and must be a known reference-is-the-variant
// haplotype (CYP3A5*3 on GRCh37) or a table error -- and the GRCh37 mtDNA is
// rCRS, so its haplogroup is H by definition. The two builds differ at the
// sites GRCh38 corrected from a rare to the major allele.
#pragma once

#include <QString>

#include "Reporter.h"

namespace TestGenome {

struct Options {
    QString dataDir;          // rsid_positions_grch37.json lives here
    QString analysisBinary;   // for --list-rsids
    QString sample = QStringLiteral("NA12878");   // or "GRCh37" / "GRCh38"
    QString output;           // default: <dataDir>/genome_1000g_<sample>.txt
                              //      or <dataDir>/genome_reference_<build>.txt
};

struct Stats {
    int positions = 0;        // rsIDs with a GRCh37 position
    int fromVcf = 0;          // genotype taken from the 1000 Genomes call
    int homRef = 0;           // no cohort variant: reference base from Ensembl
    int noCall = 0;           // indel or unresolved: written as "--"
};

bool build(const Options &options, const Reporter &reporter, Stats *stats, QString *error);

// True for the reference-assembly pseudo-samples, which need no bcftools.
bool isReferenceSample(const QString &sample);

// A build can be stopped and resumed: every chromosome's 1000 Genomes rows
// and every fetched reference base are kept under <dataDir>/test_genome_cache/
// as soon as they arrive, and the next run reads them instead of asking
// again. A chromosome cache carries the positions it was queried for, so a
// changed rsID list invalidates it.
QByteArray chromCacheHeader(const QList<qint64> &positions);
bool chromCacheMatches(const QByteArray &text, const QList<qint64> &positions, QByteArray *body);

// The 1000 Genomes file for a chromosome label, or empty when there is none.
QString vcfUrl(const QString &chrom);

// The genome-file genotype for one position from a VCF record's REF, ALT
// list and sample GT, where `offset` is the position's distance from the
// record's POS. Empty when the call cannot be written as bases: a missing
// GT, or a non-reference allele of a different length than REF (an indel).
// A homozygous-reference call at any record is the reference sequence, and
// a multi-base record whose alleles are all the same length (1000 Genomes
// writes m.2706 as AA>GA,GG,GC) yields the base at the offset.
QString genotypeFromCall(const QString &ref, const QStringList &alts, const QString &gt, qint64 offset);

} // namespace TestGenome
