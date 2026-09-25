#include <QtTest>

#include "ChainLift.h"
#include "RemoteFasta.h"

class ChainLiftTest : public QObject {
    Q_OBJECT

    // Two chains: a forward one with a gap, and a reverse-strand one. The
    // "chr"-named duplicate and the patch chain must be ignored.
    static QByteArray chainText()
    {
        return
            "chain 100 1 1000 + 100 160 1 2000 + 500 570 1\n"
            "20 10 20\n"
            "30\n"
            "\n"
            "chain 90 10 1000 + 300 320 10 5000 - 100 120 2\n"
            "20\n"
            "\n"
            "chain 80 chr1 1000 + 100 160 chr1 2000 + 500 570 3\n"
            "20 10 20\n"
            "30\n"
            "\n"
            "chain 70 HG7_PATCH 1000 + 0 10 1 2000 + 0 10 4\n"
            "10\n";
    }

private slots:
    void forwardWithGap()
    {
        ChainLift lift;
        QBuffer buf;
        buf.setData(chainText());
        buf.open(QIODevice::ReadOnly);
        QString err;
        QVERIFY2(lift.parse(buf, &err), qPrintable(err));
        QCOMPARE(lift.blockCount(), 3);   // the bare-named chains only

        ChainLift::Lifted l;
        QVERIFY(lift.lift("1", 101, &l));          // first base of block 1: t=100 -> q=500
        QCOMPARE(l.chrom, QString("1"));
        QCOMPARE(l.pos, qint64(501));
        QVERIFY(!l.reverse);
        QVERIFY(lift.lift("1", 120, &l));          // last base of block 1
        QCOMPARE(l.pos, qint64(520));
        QVERIFY(!lift.lift("1", 125, &l));         // in the 10-base gap
        QVERIFY(lift.lift("1", 131, &l));          // block 2: t=130 -> q=540
        QCOMPARE(l.pos, qint64(541));
        QVERIFY(!lift.lift("1", 161, &l));         // past the chain
        QVERIFY(!lift.lift("2", 101, &l));         // no chain
    }

    void reverseStrand()
    {
        ChainLift lift;
        QBuffer buf;
        buf.setData(chainText());
        buf.open(QIODevice::ReadOnly);
        QVERIFY(lift.parse(buf, nullptr));
        ChainLift::Lifted l;
        // t=300..320 maps to q=100..120 on the minus strand of a 5000-base
        // target: t=300 -> reverse coordinate 100 -> forward 5000-1-100.
        QVERIFY(lift.lift("10", 301, &l));
        QCOMPARE(l.chrom, QString("10"));
        QVERIFY(l.reverse);
        QCOMPARE(l.pos, qint64(5000 - 100));
        QVERIFY(lift.lift("10", 320, &l));
        QCOMPARE(l.pos, qint64(5000 - 119));
    }

    void fastaParsing()
    {
        const auto got = RemoteFasta::parseFasta(">7:99672916-99672916\nT\n>1:10-21\nacgtac\ngtacgt\n");
        QCOMPARE(got.value("7:99672916-99672916"), QString("T"));
        QCOMPARE(got.value("1:10-21"), QString("ACGTACGTACGT"));
        QCOMPARE(RemoteFasta::region("MT", 2706), QString("MT:2706-2706"));
        QCOMPARE(RemoteFasta::region("10", 100, 5), QString("10:95-105"));
    }
};

QTEST_APPLESS_MAIN(ChainLiftTest)
#include "tst_chainlift.moc"
