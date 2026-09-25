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
               "22\t400\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1\n"             // indel: dropped
               "22\t500\t.\tA\tG,T\t50\tPASS\t.\tGT\t1/2\n"            // multiallelic SNV
               "22\t600\t.\tA\tG\t50\tPASS\t.\tGT\t./.\n"              // no call: dropped
               "22\t700\t.\tA\tG\t50\tPASS\t.\tDP:GT\t9:0/0\n";        // GT not first
    }

private slots:
    void plainVcf()
    {
        const QString in = write("v.vcf", vcf());
        const QString lookup = write("l.json", "{\"rs200\": {\"chrom\": \"22\", \"pos\": \"200\"}}");
        const QString out = dir.filePath("g.txt");
        QString err;
        QCOMPARE(VcfConverter::convert(in, out, lookup, Reporter(), &err), qint64(5));
        QFile f(out);
        f.open(QIODevice::ReadOnly);
        const QStringList lines = QString::fromUtf8(f.readAll()).split('\n', Qt::SkipEmptyParts);
        QCOMPARE(lines[0], QString("# rsid\tchromosome\tposition\tgenotype"));
        QCOMPARE(lines[1], QString("rs1\t22\t100\tAG"));
        QCOMPARE(lines[2], QString("rs200\t22\t200\tTT"));
        QCOMPARE(lines[3], QString("chr22_300\t22\t300\tAA"));
        QCOMPARE(lines[4], QString("chr22_500\t22\t500\tGT"));
        QCOMPARE(lines[5], QString("chr22_700\t22\t700\tAA"));
    }

    void gzippedVcf()
    {
        const QString in = write("v2.vcf", vcf());
        QProcess gz;
        gz.start("gzip", {"-f", in});
        QVERIFY(gz.waitForFinished());
        const QString out = dir.filePath("g2.txt");
        QString err;
        QCOMPARE(VcfConverter::convert(in + ".gz", out, QString(), Reporter(), &err), qint64(5));
    }
};

QTEST_GUILESS_MAIN(VcfConvertTest)
#include "tst_vcfconvert.moc"
