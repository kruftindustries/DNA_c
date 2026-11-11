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
            "chrom\tpos\tref\talt\tall_traits\r\n"
            "X\t50\tA\tG\tplain\r\n"
            "7\t117199646\tCTT\t-\tF508del, covers the anchor base before it too\r\n"
            "2\t10\tA\tG\t\"multi\r\nline\ttrait\"\r\n"
            "1\t300\tA\tG\tp\r\n"
            "1\t300\tA\tG\tduplicate\r\n"
            "MT\t7\tA\tG\tp\r\n"
            "1\tnotanumber\tA\tG\tp\r\n"
            "\t5\tA\tG\tp\r\n");
        const QString lookup = write("l.json",
            "{\"rs1\": {\"chrom\": \"1\", \"pos\": \"20\"}, \"rs2\": {\"chrom\": \"22\", \"pos\": \"9\"}}");
        const QString bed = dir.filePath("t.bed");
        QString err;
        QCOMPARE(TargetRegions::build(cv, lookup, bed, Reporter(), &err), qint64(8));
        QFile f(bed);
        f.open(QIODevice::ReadOnly);
        const QString text = QString::fromUtf8(f.readAll());
        // (chrom.zfill(2), pos): 01 < 02 < 0X < 22 < MT
        QCOMPARE(text, QString("1\t19\t20\n1\t299\t300\n2\t9\t10\n7\t117199644\t117199645\n7\t117199645\t117199646\nX\t49\t50\n22\t8\t9\nMT\t6\t7\n"));
    }
};

QTEST_GUILESS_MAIN(TargetsTest)
#include "tst_targets.moc"
