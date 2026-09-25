#include <QtTest>
#include <QJsonArray>
#include <QSignalSpy>

#include "AnalysisRunner.h"
#include "RepoPaths.h"

// Runs the actual C binary end to end: an AncestryDNA-layout genome in, JSON
// results and an HTML report out. Skipped when the binary is not built.
class AnalysisRunnerTest : public QObject {
    Q_OBJECT

private slots:
    void initTestCase()
    {
        QCoreApplication::setOrganizationName("genetic-health-test");
        QCoreApplication::setApplicationName("genetic-health-qt-test");
        if (RepoPaths::analysisBinary().isEmpty())
            QSKIP("analysis binary not built (make -C c)");
        if (!QFileInfo(RepoPaths::dataDir() + "/clinical_annotations.tsv").exists())
            QSKIP("data directory not populated");
    }

    void ancestryGenomeEndToEnd()
    {
        QTemporaryDir out;
        const QString genome = out.filePath("g.txt");
        QFile f(genome);
        QVERIFY(f.open(QIODevice::WriteOnly));
        f.write("#AncestryDNA raw data download\n"
                "rsid\tchromosome\tposition\tallele1\tallele2\n"
                "rs429358\t19\t45411941\tC\tT\n"
                "rs7412\t19\t45412079\tC\tC\n"
                "rs4680\t22\t19951271\tA\tG\n");
        f.close();

        AnalysisRunner runner;
        QSignalSpy done(&runner, &AnalysisRunner::finished);
        AnalysisRunner::Request req;
        req.genomePath = genome;
        req.dataDir = RepoPaths::dataDir();
        req.subjectName = "Runner Test";
        req.outputDir = out.filePath("reports");
        runner.start(req);

        QVERIFY(done.wait(120000));
        QVERIFY2(done.first().at(0).toBool(), qPrintable(done.first().at(1).toString()));

        const QJsonObject results = runner.results();
        QCOMPARE(results.value("apoe").toObject().value("apoe_type").toString(), QString("e3/e4"));
        QCOMPARE(results.value("summary").toObject().value("total_snps").toInt(), 3);
        QVERIFY(!results.value("findings").toArray().isEmpty());

        QFile html(runner.reportPath());
        QVERIFY(html.open(QIODevice::ReadOnly));
        const QString text = QString::fromUtf8(html.readAll());
        QVERIFY(text.contains("Genetic Health Report — Runner Test"));
        QVERIFY(text.contains("APOE e3/e4"));
    }
};

QTEST_GUILESS_MAIN(AnalysisRunnerTest)
#include "tst_analysisrunner.moc"
