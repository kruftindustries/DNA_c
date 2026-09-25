#include <QtTest>

#include "TestGenome.h"

using TestGenome::genotypeFromCall;

class TestGenomeTest : public QObject {
    Q_OBJECT

private slots:
    void snv()
    {
        QCOMPARE(genotypeFromCall("G", {"A"}, "0|1", 0), QString("GA"));
        QCOMPARE(genotypeFromCall("G", {"A"}, "1|1", 0), QString("AA"));
        QCOMPARE(genotypeFromCall("G", {"A", "T"}, "2/0", 0), QString("TG"));
        // Haploid calls (MT, Y) are a single base.
        QCOMPARE(genotypeFromCall("G", {"A"}, "0", 0), QString("G"));
    }

    void multiBaseRecord()
    {
        // 1000 Genomes writes m.2706 as an MNP: AA > GA,GG,GC. The reference
        // call is the reference base at the offset -- the H marker.
        QCOMPARE(genotypeFromCall("AA", {"GA", "GG", "GC"}, "0", 0), QString("A"));
        QCOMPARE(genotypeFromCall("AA", {"GA", "GG", "GC"}, "1", 0), QString("G"));
        QCOMPARE(genotypeFromCall("AA", {"GA", "GG", "GC"}, "3", 1), QString("C"));
        // m.12372 as GA > AA,AG: the site is the first base.
        QCOMPARE(genotypeFromCall("GA", {"AA", "AG"}, "0", 0), QString("G"));
        QCOMPARE(genotypeFromCall("GA", {"AA", "AG"}, "2", 1), QString("G"));
    }

    void indels()
    {
        // Homozygous reference at an indel record is still the reference base.
        QCOMPARE(genotypeFromCall("C", {"CT"}, "0|0", 0), QString("CC"));
        QCOMPARE(genotypeFromCall("CTT", {"C"}, "0|0", 2), QString("TT"));
        // A carried indel allele has no base at this position.
        QCOMPARE(genotypeFromCall("C", {"CT"}, "0|1", 0), QString());
        QCOMPARE(genotypeFromCall("CTT", {"C"}, "1|1", 1), QString());
        // Copy-number and symbolic alleles cannot be written either.
        QCOMPARE(genotypeFromCall("T", {"<CN0>", "C"}, "0|1", 0), QString());
        QCOMPARE(genotypeFromCall("T", {"<CN0>", "C"}, "0|2", 0), QString("TC"));
    }

    void chromosomeCache()
    {
        const QList<qint64> positions{100, 2000, 30000};
        const QByteArray header = TestGenome::chromCacheHeader(positions);
        QVERIFY(header.startsWith("#positions 3 "));
        QByteArray body;
        QVERIFY(TestGenome::chromCacheMatches(header + "1\t100\tA\tG\t0|1\n", positions, &body));
        QCOMPARE(body, QByteArray("1\t100\tA\tG\t0|1\n"));
        // A different position list (an rsID added) invalidates the cache.
        QVERIFY(!TestGenome::chromCacheMatches(header + "x", {100, 2000, 30000, 40000}, nullptr));
        QVERIFY(!TestGenome::chromCacheMatches("garbage", positions, nullptr));
    }

    void unusable()
    {
        QCOMPARE(genotypeFromCall("G", {"A"}, ".", 0), QString());
        QCOMPARE(genotypeFromCall("G", {"A"}, "./.", 0), QString());
        QCOMPARE(genotypeFromCall("G", {"A"}, "0|5", 0), QString());
        QCOMPARE(genotypeFromCall("G", {"A"}, "0|1", 1), QString());
    }
};

QTEST_APPLESS_MAIN(TestGenomeTest)
#include "tst_testgenome.moc"
