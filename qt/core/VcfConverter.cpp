#include "VcfConverter.h"

#include <QFile>
#include <QHash>
#include <QJsonObject>

#include "EnsemblLookup.h"
#include "GzipStream.h"
#include "TextTable.h"

namespace VcfConverter {

namespace {

bool isBases(const QByteArray &s)
{
    if (s.isEmpty())
        return false;
    for (char c : s)
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            return false;
    return true;
}

} // namespace

qint64 convert(const QString &vcfPath, const QString &genomeOut, const QString &lookupJson,
               const Reporter &reporter, QString *error)
{

    QHash<QByteArray, QByteArray> posToRsid;
    const QJsonObject lookup = EnsemblLookup::loadLookup(lookupJson);
    for (auto it = lookup.begin(); it != lookup.end(); ++it) {
        const QJsonObject info = it.value().toObject();
        posToRsid[(info.value("chrom").toString() + ':' + info.value("pos").toString()).toUtf8()] =
            it.key().toUtf8();
    }

    QIODevice *in;
    GzipStream gz(vcfPath);
    QFile plain(vcfPath);
    if (vcfPath.endsWith(".gz")) {
        if (!gz.open(QIODevice::ReadOnly)) {
            if (error)
                *error = gz.errorText();
            return -1;
        }
        in = &gz;
    } else {
        if (!plain.open(QIODevice::ReadOnly)) {
            if (error)
                *error = plain.errorString();
            return -1;
        }
        in = &plain;
    }

    QFile out(genomeOut);
    if (!out.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (error)
            *error = QStringLiteral("cannot write %1").arg(genomeOut);
        return -1;
    }
    out.write("# rsid\tchromosome\tposition\tgenotype\n");

    LineReader lines(in);
    QByteArray line;
    qint64 count = 0;
    while (lines.readLine(line)) {
        if (line.startsWith('#'))
            continue;
        // line.split('\t', 10): at most eleven fields, the last holding the rest.
        QList<QByteArray> fields;
        int start = 0;
        for (int splits = 0; splits < 10; ++splits) {
            const int tab = line.indexOf('\t', start);
            if (tab < 0)
                break;
            fields << line.mid(start, tab - start);
            start = tab + 1;
        }
        fields << line.mid(start);
        if (fields.size() < 10)
            continue;

        const QByteArray chrom = fields[0], pos = fields[1], vid = fields[2];
        const QByteArray ref = fields[3], alt = fields[4];
        const QList<QByteArray> alts = alt.split(',');
        bool snv = ref.size() == 1;
        for (const QByteArray &a : alts)
            if (a.size() != 1)
                snv = false;
        // A biallelic indel is written the way an array would call it: D
        // for the shorter allele, I for the longer. Several alternate
        // alleles, or alleles of equal length (an MNV), have no such
        // spelling and are dropped as before.
        const bool indel = !snv && alts.size() == 1 && ref.size() != alts[0].size() &&
                           isBases(ref) && isBases(alts[0]);
        if (!snv && !indel)
            continue;

        const QList<QByteArray> fmt = fields[8].split(':');
        const QList<QByteArray> sample = fields[9].split(':');
        const int gtIdx = fmt.indexOf("GT");
        if (gtIdx < 0 || gtIdx >= sample.size())
            continue;
        const QByteArray gt = sample[gtIdx];
        const char sep = gt.contains('/') ? '/' : '|';
        const QList<QByteArray> parts = gt.split(sep);
        if (parts.size() != 2 || parts[0].contains('.') || parts[1].contains('.'))
            continue;

        QList<QByteArray> alleles{ref};
        alleles += alts;
        if (indel) {
            const bool refShorter = ref.size() < alts[0].size();
            alleles = {refShorter ? "D" : "I", refShorter ? "I" : "D"};
        }
        bool ok1 = false, ok2 = false;
        const int i1 = parts[0].toInt(&ok1), i2 = parts[1].toInt(&ok2);
        if (!ok1 || !ok2 || i1 < 0 || i2 < 0 || i1 >= alleles.size() || i2 >= alleles.size())
            continue;
        const QByteArray genotype = alleles[i1] + alleles[i2];

        QByteArray rsid;
        const QByteArray posKey = chrom + ':' + pos;
        // dbSNP places an indel on the first changed base, one past the
        // VCF anchor the caller reports it at.
        const QByteArray nextKey = chrom + ':' + QByteArray::number(pos.toLongLong() + 1);
        if (vid != "." && vid.startsWith("rs"))
            rsid = vid;
        else if (posToRsid.contains(posKey))
            rsid = posToRsid.value(posKey);
        else if (indel && posToRsid.contains(nextKey))
            rsid = posToRsid.value(nextKey);
        else
            rsid = "chr" + chrom + '_' + pos;

        out.write(rsid + '\t' + chrom + '\t' + pos + '\t' + genotype + '\n');
        count++;
    }
    reporter.log(QStringLiteral("  Converted %L1 variants (SNVs and D/I indel calls) to 23andMe format -> %2")
                     .arg(count).arg(genomeOut));
    return count;
}

} // namespace VcfConverter
