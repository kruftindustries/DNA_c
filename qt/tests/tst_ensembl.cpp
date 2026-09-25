#include <QtTest>
#include <QJsonArray>
#include <QJsonDocument>

#include "EnsemblLookup.h"

class EnsemblTest : public QObject {
    Q_OBJECT
    QTemporaryDir dir;

private slots:
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
