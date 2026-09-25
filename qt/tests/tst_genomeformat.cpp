#include <QtTest>

#include "GenomeFormat.h"

class GenomeFormatTest : public QObject {
    Q_OBJECT

    QTemporaryDir dir;

    QString write(const QString &name, const QByteArray &content)
    {
        const QString path = dir.filePath(name);
        QFile f(path);
        f.open(QIODevice::WriteOnly);
        f.write(content);
        return path;
    }

private slots:
    void twentyThreeAndMe()
    {
        const GenomeSniff s = sniffGenome(write("g.txt",
            "# rsid\tchromosome\tposition\tgenotype\n"
            "rsid\tchromosome\tposition\tgenotype\n"
            "rs1\t1\t100\tAG\nrs2\t1\t200\t--\nrs3\tX\t300\tT\n"));
        QCOMPARE(s.format, GenomeFormat::TwentyThreeAndMe);
        QCOMPARE(s.dataRows, 2);
        QCOMPARE(s.noCalls, 1);
        QVERIFY(isArrayFormat(s.format));
    }

    void ancestryDna()
    {
        const GenomeSniff s = sniffGenome(write("a.txt",
            "#AncestryDNA raw data download\n"
            "rsid\tchromosome\tposition\tallele1\tallele2\n"
            "rs1\t1\t100\tA\tG\nrs2\t23\t200\tC\tC\nrs3\t1\t300\t0\t0\nrs4\t1\t400\tA\t0\n"));
        QCOMPARE(s.format, GenomeFormat::AncestryDNA);
        QCOMPARE(s.dataRows, 2);
        QCOMPARE(s.noCalls, 2);
        QVERIFY(isArrayFormat(s.format));
    }

    void vcf()
    {
        const GenomeSniff s = sniffGenome(write("v.vcf",
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n1\t100\trs1\tA\tG\n"));
        QCOMPARE(s.format, GenomeFormat::Vcf);
        QVERIFY(!isArrayFormat(s.format));
    }

    void fastq()
    {
        const GenomeSniff s = sniffGenome(write("r.fastq",
            "@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nACGT\n+\nIIII\n"));
        QCOMPARE(s.format, GenomeFormat::Fastq);
    }

    void gzipByName()
    {
        const QByteArray magic("\x1f\x8b\x08\x00garbage", 11);
        QCOMPARE(sniffGenome(write("r.fastq.gz", magic)).format, GenomeFormat::Fastq);
        QCOMPARE(sniffGenome(write("x.txt.gz", magic)).format, GenomeFormat::GzipCompressed);
    }

    void headerOnlyFileNamed()
    {
        const GenomeSniff s = sniffGenome(write("empty.txt",
            "# The GRCh38 reference assembly as a sample\n# rsid\tchromosome\tposition\tgenotype\n"));
        QCOMPARE(s.format, GenomeFormat::Unknown);
        QVERIFY2(s.description.contains("only 2 comment line"), qPrintable(s.description));
        QCOMPARE(sniffGenome(write("zero.txt", "")).description, QString("empty file"));
    }

    void missingFile()
    {
        const GenomeSniff s = sniffGenome(dir.filePath("nope.txt"));
        QCOMPARE(s.format, GenomeFormat::Unknown);
        QVERIFY(s.description.startsWith("cannot open"));
    }

    void rowClassifier()
    {
        QCOMPARE(classifyRow({"rs1", "1", "100", "AG"}), RowKind::Ok23andMe);
        QCOMPARE(classifyRow({"rs1", "1", "100", "--"}), RowKind::NoCall);
        QCOMPARE(classifyRow({"rs1", "1", "100", "XY"}), RowKind::Invalid);
        QCOMPARE(classifyRow({"rs1", "1", "100"}), RowKind::Short);
        QCOMPARE(classifyRow({"rs1", "1", "100", "A", "G"}), RowKind::OkAncestry);
        QCOMPARE(classifyRow({"rs1", "1", "100", "0", "0"}), RowKind::NoCall);
        QCOMPARE(classifyRow({"rs1", "1", "100", "A", "0"}), RowKind::NoCall);
        QCOMPARE(classifyRow({"rs1", "1", "100", "N", "N"}), RowKind::Invalid);
        // The vendor header is not a genotype and must not look like one.
        QCOMPARE(classifyRow({"rsid", "chromosome", "position", "allele1", "allele2"}),
                 RowKind::Invalid);
    }
};

QTEST_GUILESS_MAIN(GenomeFormatTest)
#include "tst_genomeformat.moc"
