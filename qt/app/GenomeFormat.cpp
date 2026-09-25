#include "GenomeFormat.h"

#include <QFile>
#include <QFileInfo>
#include <QTextStream>

namespace {

bool isBases(const QString &s)
{
    if (s.isEmpty() || s.size() > 2)
        return false;
    for (QChar c : s)
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            return false;
    return true;
}

} // namespace

QString formatName(GenomeFormat format)
{
    switch (format) {
    case GenomeFormat::TwentyThreeAndMe: return QStringLiteral("23andMe raw data");
    case GenomeFormat::AncestryDNA:      return QStringLiteral("AncestryDNA raw data");
    case GenomeFormat::Vcf:              return QStringLiteral("VCF");
    case GenomeFormat::Fastq:            return QStringLiteral("FASTQ (sequencing reads)");
    case GenomeFormat::GzipCompressed:   return QStringLiteral("gzip-compressed");
    case GenomeFormat::Unknown:          break;
    }
    return QStringLiteral("unrecognised");
}

bool isArrayFormat(GenomeFormat format)
{
    return format == GenomeFormat::TwentyThreeAndMe || format == GenomeFormat::AncestryDNA;
}

RowKind classifyRow(const QStringList &parts)
{
    // Mirrors loading.py:parse_genome_row and gh_genome_parse_row.
    if (parts.size() >= 5 && parts[3].size() == 1 && parts[4].size() == 1) {
        if (parts[3] == "0" || parts[4] == "0")
            return RowKind::NoCall;
        return isBases(parts[3] + parts[4]) ? RowKind::OkAncestry : RowKind::Invalid;
    }
    if (parts.size() >= 4) {
        if (parts[3] == "--")
            return RowKind::NoCall;
        return isBases(parts[3]) ? RowKind::Ok23andMe : RowKind::Invalid;
    }
    return RowKind::Short;
}

GenomeSniff sniffGenome(const QString &path, int maxLines)
{
    GenomeSniff out;
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) {
        out.description = QStringLiteral("cannot open: %1").arg(f.errorString());
        return out;
    }

    const QByteArray head = f.peek(4);
    if (head.size() >= 2 && quint8(head[0]) == 0x1f && quint8(head[1]) == 0x8b) {
        const QString name = QFileInfo(path).fileName().toLower();
        if (name.endsWith(".fastq.gz") || name.endsWith(".fq.gz")) {
            out.format = GenomeFormat::Fastq;
            out.description = QStringLiteral("gzip-compressed FASTQ");
        } else {
            out.format = GenomeFormat::GzipCompressed;
            out.description = QStringLiteral("gzip-compressed file; decompress it first");
        }
        return out;
    }

    QTextStream in(&f);
    int rows23 = 0, rowsAnc = 0, vcfHeader = 0;
    bool firstNonEmpty = true;
    for (int i = 0; i < maxLines && !in.atEnd(); ++i) {
        const QString line = in.readLine();
        out.sampledLines++;
        if (line.isEmpty())
            continue;
        if (firstNonEmpty) {
            firstNonEmpty = false;
            if (line.startsWith("##fileformat=VCF")) {
                out.format = GenomeFormat::Vcf;
                out.description = QStringLiteral("VCF; convert with the WGS pipeline's "
                                                 "conversion step before analysis");
                return out;
            }
            if (line.startsWith('@') && !line.startsWith("@rs")) {
                // A FASTQ record is @name / bases / + / qualities.
                const QString bases = in.readLine();
                const QString plus = in.readLine();
                if (plus.startsWith('+') && !bases.isEmpty()) {
                    out.format = GenomeFormat::Fastq;
                    out.description = QStringLiteral("FASTQ sequencing reads; run the WGS pipeline");
                    return out;
                }
            }
        }
        if (line.startsWith('#')) {
            if (line.startsWith("#CHROM"))
                vcfHeader++;
            continue;
        }
        switch (classifyRow(line.split(QLatin1Char('\t')))) {
        case RowKind::Ok23andMe: rows23++; out.dataRows++; break;
        case RowKind::OkAncestry: rowsAnc++; out.dataRows++; break;
        case RowKind::NoCall: out.noCalls++; break;
        default: break;
        }
    }

    if (vcfHeader) {
        out.format = GenomeFormat::Vcf;
        out.description = QStringLiteral("VCF without a fileformat line");
        return out;
    }
    if (rowsAnc && rowsAnc >= rows23) {
        out.format = GenomeFormat::AncestryDNA;
        out.description = QStringLiteral("AncestryDNA layout (allele1/allele2 columns): "
                                         "%1 genotypes, %2 no-calls in the first %3 lines")
                              .arg(out.dataRows).arg(out.noCalls).arg(out.sampledLines);
    } else if (rows23) {
        out.format = GenomeFormat::TwentyThreeAndMe;
        out.description = QStringLiteral("23andMe layout: %1 genotypes, %2 no-calls "
                                         "in the first %3 lines")
                              .arg(out.dataRows).arg(out.noCalls).arg(out.sampledLines);
    } else {
        out.description = QStringLiteral("no genotype rows found in the first %1 lines")
                              .arg(out.sampledLines);
    }
    return out;
}
