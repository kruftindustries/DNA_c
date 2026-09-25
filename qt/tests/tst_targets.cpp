#include <QtTest>

#include "Reporter.h"
#include "TargetRegions.h"

class TargetsTest : public QObject {
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

private slots:
    void sortedDedupedBed()
    {
        // A quoted field spanning lines must not yield a bogus row.
        const QString cv = write("cv.tsv",
            "chrom\tpos\tref\tall_traits\r\n"
            "X\t50\tA\tplain\r\n"
            "2\t10\tA\t\"multi\r\nline\ttrait\"\r\n"
            "1\t300\tA\tp\r\n"
            "1\t300\tA\tduplicate\r\n"
            "MT\t7\tA\tp\r\n"
            "1\tnotanumber\tA\tp\r\n"
            "\t5\tA\tp\r\n");
        const QString lookup = write("l.json",
            "{\"rs1\": {\"chrom\": \"1\", \"pos\": \"20\"}, \"rs2\": {\"chrom\": \"22\", \"pos\": \"9\"}}");
        const QString bed = dir.filePath("t.bed");
        QString err;
        QCOMPARE(TargetRegions::build(cv, lookup, bed, Reporter(), &err), qint64(6));
        QFile f(bed);
        f.open(QIODevice::ReadOnly);
        const QString text = QString::fromUtf8(f.readAll());
        // (chrom.zfill(2), pos): 01 < 02 < 0X < 22 < MT
        QCOMPARE(text, QString("1\t19\t20\n1\t299\t300\n2\t9\t10\nX\t49\t50\n22\t8\t9\nMT\t6\t7\n"));
    }
};

QTEST_GUILESS_MAIN(TargetsTest)
#include "tst_targets.moc"
