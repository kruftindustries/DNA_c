#include <QtTest>
#include <QProcess>

#include "Reporter.h"
#include "VcfConverter.h"

class VcfConvertTest : public QObject {
    Q_OBJECT
    QTemporaryDir dir;

    QString write(const QString &name, const QByteArray &content)
    {
        const QString p = dir.filePath(name);
        QFile f(p);
        f.open(QIODevice::WriteOnly);
        f.write(content);
        return p;
    }

    static QByteArray vcf()
    {
        return "##fileformat=VCFv4.2\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
               "22\t100\trs1\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:20\n"      // het, ID given
               "22\t200\t.\tC\tT\t50\tPASS\t.\tGT\t1/1\n"              // hom alt, from lookup
               "22\t300\t.\tG\tA\t50\tPASS\t.\tGT\t1|1\n"              // phased, positional
               "22\t400\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1\n"             // deletion: D/I call, rsID from pos+1
               "22\t450\t.\tA\tATT\t50\tPASS\t.\tGT\t1/1\n"            // insertion: II
               "22\t460\t.\tAT\tA,ATT\t50\tPASS\t.\tGT\t1/2\n"         // multi-allelic indel: dropped
               "22\t470\t.\tAT\tGC\t50\tPASS\t.\tGT\t0/1\n"            // MNV: dropped
               "22\t500\t.\tA\tG,T\t50\tPASS\t.\tGT\t1/2\n"            // multiallelic SNV
               "22\t600\t.\tA\tG\t50\tPASS\t.\tGT\t./.\n"              // no call: dropped
               "22\t700\t.\tA\tG\t50\tPASS\t.\tDP:GT\t9:0/0\n";        // GT not first
    }

private slots:
    void plainVcf()
    {
        const QString in = write("v.vcf", vcf());
        // rs401 sits on the first deleted base, one past the VCF anchor at 400.
        const QString lookup = write("l.json", "{\"rs200\": {\"chrom\": \"22\", \"pos\": \"200\"},"
                                                " \"rs401\": {\"chrom\": \"22\", \"pos\": \"401\"}}");
        const QString out = dir.filePath("g.txt");
        QString err;
        QCOMPARE(VcfConverter::convert(in, out, lookup, Reporter(), &err), qint64(7));
        QFile f(out);
        f.open(QIODevice::ReadOnly);
        const QStringList lines = QString::fromUtf8(f.readAll()).split('\n', Qt::SkipEmptyParts);
        QCOMPARE(lines[0], QString("# rsid\tchromosome\tposition\tgenotype"));
        QCOMPARE(lines[1], QString("rs1\t22\t100\tAG"));
        QCOMPARE(lines[2], QString("rs200\t22\t200\tTT"));
        QCOMPARE(lines[3], QString("chr22_300\t22\t300\tAA"));
        QCOMPARE(lines[4], QString("rs401\t22\t400\tID"));
        QCOMPARE(lines[5], QString("chr22_450\t22\t450\tII"));
        QCOMPARE(lines[6], QString("chr22_500\t22\t500\tGT"));
        QCOMPARE(lines[7], QString("chr22_700\t22\t700\tAA"));
    }

    void gzippedVcf()
    {
        const QString in = write("v2.vcf", vcf());
        QProcess gz;
        gz.start("gzip", {"-f", in});
        QVERIFY(gz.waitForFinished());
        const QString out = dir.filePath("g2.txt");
        QString err;
        QCOMPARE(VcfConverter::convert(in + ".gz", out, QString(), Reporter(), &err), qint64(7));
    }
};

QTEST_GUILESS_MAIN(VcfConvertTest)
#include "tst_vcfconvert.moc"
