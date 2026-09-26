#include <QtTest>
#include <QJsonArray>
#include <QJsonDocument>

#include "EnsemblLookup.h"

class EnsemblTest : public QObject {
    Q_OBJECT
    QTemporaryDir dir;

private slots:
    // A data directory without the lookup starts from the copy compiled
    // in: that is what a packaged build has, and it is what lets the test
    // genomes and reference samples be built when the Ensembl GRCh37 host
    // is down (which used to leave an empty sample file behind).
    void missingLookupSeededFromResource()
    {
        const QString path = dir.filePath("fresh/rsid_positions_grch37.json");
        QVERIFY(!QFile::exists(path));
        const QJsonObject lookup = EnsemblLookup::loadLookup(path);
        QVERIFY(QFile::exists(path));
        QVERIFY(lookup.size() > 500);
        QCOMPARE(lookup["rs429358"].toObject()["chrom"].toString(), QString("19"));
        QCOMPARE(lookup["rs429358"].toObject()["pos"].toString(), QString("45411941"));
        // and the copy is writable, so an update can extend it
        QVERIFY(QFileInfo(path).isWritable());
        QVERIFY(EnsemblLookup::saveLookup(path, lookup));
    }

    // A lookup an earlier build saved empty (every Ensembl batch failed)
    // is seeded like a missing one; a lookup with entries is left alone.
    void emptyLookupSeededOwnEntriesKept()
    {
        const QString path = dir.filePath("partial/rsid_positions_grch37.json");
        QDir().mkpath(dir.filePath("partial"));
        QJsonObject mine;
        mine.insert("rs429358", QJsonObject{{"chrom", "19"}, {"pos", "1"}});
        QVERIFY(EnsemblLookup::saveLookup(path, mine));
        QCOMPARE(EnsemblLookup::loadLookup(path), mine);

        const QString empty = dir.filePath("empty/rsid_positions_grch37.json");
        QDir().mkpath(dir.filePath("empty"));
        QVERIFY(EnsemblLookup::saveLookup(empty, QJsonObject()));
        QVERIFY(EnsemblLookup::loadLookup(empty).size() > 500);
    }

    void mergedRsidStoredUnderRequestedId()
    {
        const QJsonObject results = QJsonDocument::fromJson(R"({
            "rs193302994": {"name": "rs193302994", "synonyms": ["rs3088309", "rs377204770"],
                            "mappings": [{"seq_region_name": "MT", "start": 15452}]}
        })").object();
        QJsonObject lookup;
        QCOMPARE(EnsemblLookup::mergeBatch(results, {"rs3088309"}, lookup), 1);
        QCOMPARE(lookup.keys(), QStringList{"rs3088309"});
        QCOMPARE(lookup["rs3088309"].toObject()["pos"].toString(), QString("15452"));
        QVERIFY(!lookup.contains("rs193302994"));
    }

    void patchScaffoldSkipped()
    {
        const QJsonObject results = QJsonDocument::fromJson(R"({
            "rs56392308": {"name": "rs56392308", "synonyms": ["rs8176750"],
                           "mappings": [{"seq_region_name": "HG79_PATCH", "start": 136131204},
                                        {"seq_region_name": "9", "start": 136131057}]}
        })").object();
        QJsonObject lookup;
        EnsemblLookup::mergeBatch(results, {"rs8176750"}, lookup);
        QCOMPARE(lookup["rs8176750"].toObject()["chrom"].toString(), QString("9"));
        QCOMPARE(lookup["rs8176750"].toObject()["pos"].toString(), QString("136131057"));
    }

    void noUsableMappingLeftOut()
    {
        const QJsonObject results = QJsonDocument::fromJson(R"({
            "rs999": {"name": "rs999", "synonyms": [],
                      "mappings": [{"seq_region_name": "HG79_PATCH", "start": 1}]}
        })").object();
        QJsonObject lookup;
        QCOMPARE(EnsemblLookup::mergeBatch(results, {"rs999"}, lookup), 0);
        QVERIFY(lookup.isEmpty());
    }

    void saveMatchesPythonLayout()
    {
        QJsonObject lookup;
        lookup["rs1"] = QJsonObject{{"chrom", "1"}, {"pos", "100"}};
        const QString path = dir.filePath("l.json");
        QVERIFY(EnsemblLookup::saveLookup(path, lookup));
        QFile f(path);
        f.open(QIODevice::ReadOnly);
        QCOMPARE(QString::fromUtf8(f.readAll()),
                 QString("{\n  \"rs1\": {\n    \"chrom\": \"1\",\n    \"pos\": \"100\"\n  }\n}"));
        QCOMPARE(EnsemblLookup::loadLookup(path), lookup);
    }
};

QTEST_GUILESS_MAIN(EnsemblTest)
#include "tst_ensembl.moc"
