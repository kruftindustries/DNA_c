#include <QtTest>
#include <QJsonArray>
#include <QJsonObject>

#include "FindingsModel.h"

class FindingsTest : public QObject {
    Q_OBJECT

    static QJsonArray sample()
    {
        auto f = [](const char *rsid, const char *gene, const char *cat, int mag, const char *desc) {
            return QJsonObject{{"rsid", rsid}, {"gene", gene}, {"category", cat},
                               {"genotype", "AG"}, {"status", "reduced"},
                               {"description", desc}, {"magnitude", mag}, {"note", ""}};
        };
        return QJsonArray{f("rs1801133", "MTHFR", "Methylation", 3, "Reduced MTHFR"),
                          f("rs762551", "CYP1A2", "Drug Metabolism", 1, "Fast caffeine"),
                          f("rs4244285", "CYP2C19", "Drug Metabolism", 4, "Poor CYP2C19")};
    }

private slots:
    void modelShape()
    {
        FindingsModel m;
        m.setFindings(sample());
        QCOMPARE(m.rowCount(), 3);
        QCOMPARE(m.columnCount(), int(FindingsModel::ColumnCount));
        QCOMPARE(m.data(m.index(0, FindingsModel::Gene)).toString(), QString("MTHFR"));
        QCOMPARE(m.data(m.index(2, FindingsModel::Magnitude)).toInt(), 4);
        QCOMPARE(m.categories(), QStringList({"Drug Metabolism", "Methylation"}));
    }

    void filters()
    {
        FindingsModel m;
        m.setFindings(sample());
        FindingsFilter p;
        p.setSourceModel(&m);
        QCOMPARE(p.rowCount(), 3);

        p.setMinimumMagnitude(3);
        QCOMPARE(p.rowCount(), 2);

        p.setCategory("Drug Metabolism");
        QCOMPARE(p.rowCount(), 1);
        QCOMPARE(p.data(p.index(0, FindingsModel::Gene)).toString(), QString("CYP2C19"));

        p.setCategory(QString());
        p.setMinimumMagnitude(0);
        p.setSearchText("caffeine");
        QCOMPARE(p.rowCount(), 1);
        QCOMPARE(p.data(p.index(0, FindingsModel::Rsid)).toString(), QString("rs762551"));
    }

    void sortsMagnitudeNumerically()
    {
        FindingsModel m;
        m.setFindings(sample());
        FindingsFilter p;
        p.setSourceModel(&m);
        p.sort(FindingsModel::Magnitude, Qt::DescendingOrder);
        QCOMPARE(p.data(p.index(0, FindingsModel::Magnitude)).toInt(), 4);
        QCOMPARE(p.data(p.index(2, FindingsModel::Magnitude)).toInt(), 1);
    }

    void tsvExport()
    {
        FindingsModel m;
        m.setFindings(sample());
        const QString tsv = FindingsModel::toTsv(m.findings());
        const QStringList lines = tsv.split('\n', Qt::SkipEmptyParts);
        QCOMPARE(lines.size(), 4);
        QVERIFY(lines[0].startsWith("rsid\tgene\t"));
        QVERIFY(lines[1].contains("\tMTHFR\t"));
    }
};

QTEST_GUILESS_MAIN(FindingsTest)
#include "tst_findings.moc"
