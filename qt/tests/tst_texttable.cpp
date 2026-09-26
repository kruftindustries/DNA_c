#include <QtTest>
#include <QBuffer>

#include "TextTable.h"

// The reader and writer against Python's csv module behaviour on the cases
// that a naive split gets wrong.
class TextTableTest : public QObject {
    Q_OBJECT

    static QList<QStringList> readAll(const QByteArray &data)
    {
        QBuffer buf;
        buf.setData(data);
        buf.open(QIODevice::ReadOnly);
        TsvReader r(&buf);
        QList<QStringList> rows;
        QStringList f;
        while (r.readRecord(f))
            rows << f;
        return rows;
    }

private slots:
    void outputLinesStripCarriageReturns()
    {
        const QStringList ids = outputLines(QString("rs1\r\nrs2\r\n\r\nrs3\n"));
        QCOMPARE(ids, (QStringList{"rs1", "rs2", "rs3"}));
        QCOMPARE(outputLines(QString("rs1\nrs2")), (QStringList{"rs1", "rs2"}));
        const QList<QByteArray> rows = outputLines(QByteArray("22\t100\tA\tG\t0|1\r\n"));
        QCOMPARE(rows.size(), 1);
        QVERIFY(!rows[0].endsWith('\r'));
    }

    void plainRows()
    {
        const auto rows = readAll("a\tb\tc\r\n1\t2\t3\n");
        QCOMPARE(rows.size(), 2);
        QCOMPARE(rows[1], QStringList({"1", "2", "3"}));
    }

    void quotedFieldWithTabAndNewline()
    {
        const auto rows = readAll("x\t\"has\ttab\r\nand line\"\ty\nnext\t1\t2\n");
        QCOMPARE(rows.size(), 2);
        QCOMPARE(rows[0][1], QString("has\ttab\nand line"));
        QCOMPARE(rows[0][2], QString("y"));
        QCOMPARE(rows[1][0], QString("next"));
    }

    void doubledQuoteInsideQuoted()
    {
        const auto rows = readAll("\"say \"\"hi\"\"\"\tz\n");
        QCOMPARE(rows[0][0], QString("say \"hi\""));
    }

    void quoteInMiddleIsLiteral()
    {
        // Only a quote at the very start of a field opens quoting.
        const auto rows = readAll("5' UTR\tb\n");
        QCOMPARE(rows[0][0], QString("5' UTR"));
        const auto rows2 = readAll("ab\"cd\tb\n");
        QCOMPARE(rows2[0][0], QString("ab\"cd"));
    }

    void textAfterClosingQuoteIsKept()
    {
        const auto rows = readAll("\"q\"tail\tb\n");
        QCOMPARE(rows[0][0], QString("qtail"));
    }

    void blankLinesSkippedAndNoTrailingNewline()
    {
        const auto rows = readAll("a\tb\n\n\nc\td");
        QCOMPARE(rows.size(), 2);
        QCOMPARE(rows[1], QStringList({"c", "d"}));
    }

    void emptyFieldsPreserved()
    {
        const auto rows = readAll("a\t\t\n");
        QCOMPARE(rows[0], QStringList({"a", "", ""}));
    }

    void table()
    {
        QBuffer buf;
        buf.setData(QByteArray("chrom\tpos\tname\n1\t100\tx\n2\t200\ty\n"));
        buf.open(QIODevice::ReadOnly);
        TsvTable t(&buf);
        QVERIFY(t.ok());
        QCOMPARE(t.column("pos"), 1);
        QVERIFY(t.next());
        QCOMPARE(t.field("name"), QString("x"));
        QCOMPARE(t.field(9), QString());
        QVERIFY(t.next());
        QVERIFY(!t.next());
    }

    void writerQuotesMinimally()
    {
        QCOMPARE(tsvRow({"a", "b"}), QByteArray("a\tb\r\n"));
        QCOMPARE(tsvRow({"has\ttab", "plain"}), QByteArray("\"has\ttab\"\tplain\r\n"));
        QCOMPARE(tsvRow({"say \"hi\""}), QByteArray("\"say \"\"hi\"\"\"\r\n"));
        QCOMPARE(tsvRow({"line\nbreak"}), QByteArray("\"line\nbreak\"\r\n"));
        QCOMPARE(tsvRow({"", ""}), QByteArray("\t\r\n"));
        QCOMPARE(tsvRow({"5' UTR"}), QByteArray("5' UTR\r\n"));
    }

    void roundTrip()
    {
        const QStringList fields{"a\tb", "quote\"d", "multi\nline", "plain"};
        const auto rows = readAll(tsvRow(fields));
        QCOMPARE(rows.size(), 1);
        QCOMPARE(rows[0], fields);
    }

    void lineReader()
    {
        QBuffer buf;
        buf.setData(QByteArray("one\r\ntwo\nthree"));
        buf.open(QIODevice::ReadOnly);
        LineReader lr(&buf);
        QByteArray l;
        QVERIFY(lr.readLine(l)); QCOMPARE(l, QByteArray("one"));
        QVERIFY(lr.readLine(l)); QCOMPARE(l, QByteArray("two"));
        QVERIFY(lr.readLine(l)); QCOMPARE(l, QByteArray("three"));
        QVERIFY(!lr.readLine(l));
    }
};

QTEST_GUILESS_MAIN(TextTableTest)
#include "tst_texttable.moc"
