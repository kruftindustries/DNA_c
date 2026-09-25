// Detects what kind of file the user handed us, from its contents.
//
// The extension is not trusted: vendors ship 23andMe and AncestryDNA exports
// both as .txt, and the two layouts differ in a way that used to be read
// silently wrong (AncestryDNA's split allele columns parsed as haploid
// calls). The same per-row rule the loaders use decides here.
#pragma once

#include <QString>

enum class GenomeFormat {
    Unknown,
    TwentyThreeAndMe,   // rsid chrom pos genotype
    AncestryDNA,        // rsid chrom pos allele1 allele2
    Vcf,
    Fastq,
    GzipCompressed,     // gzip magic; a FASTQ if the name says so
};

struct GenomeSniff {
    GenomeFormat format = GenomeFormat::Unknown;
    QString description;   // one line for the UI
    int dataRows = 0;      // genotype rows seen in the sampled lines
    int noCalls = 0;
    int sampledLines = 0;
};

QString formatName(GenomeFormat format);

// True for the layouts the analysis binary reads directly.
bool isArrayFormat(GenomeFormat format);

// Read up to maxLines lines and classify. Never throws; an unreadable file
// comes back Unknown with the reason in `description`.
GenomeSniff sniffGenome(const QString &path, int maxLines = 4000);

// The classifier itself, exposed for tests: one tab-split row.
enum class RowKind { Ok23andMe, OkAncestry, NoCall, Invalid, Short };
RowKind classifyRow(const QStringList &parts);
