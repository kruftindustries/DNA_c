#include <QtTest>
#include <QProcess>

#include "ClinVarUpdater.h"
#include "Reporter.h"

class ClinVarProcessTest : public QObject {
    Q_OBJECT
    QTemporaryDir dir;

    QString writeGz(const QString &name, const QByteArray &content)
    {
        const QString p = dir.filePath(name);
        QFile f(p);
        f.open(QIODevice::WriteOnly);
        f.write(content);
        f.close();
        QProcess gz;
        gz.start("gzip", {"-f", p});
        gz.waitForFinished();
        return p + ".gz";
    }

private slots:
    void goldStars()
    {
        QCOMPARE(ClinVarUpdater::goldStars("practice guideline"), 4);
        QCOMPARE(ClinVarUpdater::goldStars("  Reviewed by Expert Panel "), 4);
        QCOMPARE(ClinVarUpdater::goldStars("criteria provided, multiple submitters, no conflicts"), 3);
        QCOMPARE(ClinVarUpdater::goldStars("criteria provided, conflicting interpretations"), 2);
        QCOMPARE(ClinVarUpdater::goldStars("criteria provided, single submitter"), 1);
        QCOMPARE(ClinVarUpdater::goldStars("no assertion criteria provided"), 0);
        QCOMPARE(ClinVarUpdater::goldStars("something new"), 0);
    }

    void filtersAndColumns()
    {
        const QByteArray header =
            "#AlleleID\tGeneSymbol\tClinicalSignificance\tPhenotypeList\tOriginSimple\tAssembly\t"
            "Chromosome\tStart\tReviewStatus\tOtherIDs\tReferenceAlleleVCF\tAlternateAlleleVCF\tMolecularConsequence\n";
        const QByteArray rows =
            "1\tBRCA1\tPathogenic\tBreast cancer\tgermline\tGRCh37\t17\t41276045\tpractice guideline\tOMIM:1\tA\tG\tmissense\n"
            "2\tBRCA1\tPathogenic\tBreast cancer\tgermline\tGRCh38\t17\t43124000\tpractice guideline\tOMIM:1\tA\tG\tmissense\n"   // wrong assembly
            "3\tX\tBenign\tp\tgermline\tGRCh37\t\t100\tno assertion provided\t\tA\tG\t\n"                                       // no chrom
            "4\tX\tBenign\tp\tgermline\tGRCh37\t1\t100\tno assertion provided\t\t\tG\t\n"                                       // no ref
            "5\tHFE\t\"Pathogenic; risk factor\"\t\"Trait with\ttab\"\tgermline\tGRCh37\t6\t26093141\tcriteria provided, single submitter\t\tG\tA\t\n";
        const QString gz = writeGz("vs.txt", header + rows);
        const QString out = dir.filePath("out.tsv");
        ClinVarUpdater::Stats stats;
        QString err;
        QVERIFY2(ClinVarUpdater::process(gz, out, &stats, Reporter(), &err), qPrintable(err));
        QCOMPARE(stats.rowsProcessed, qint64(5));
        QCOMPARE(stats.written, qint64(2));

        QFile f(out);
        f.open(QIODevice::ReadOnly);
        const QByteArray text = f.readAll();
        const QList<QByteArray> lines = text.split('\n');
        QCOMPARE(lines[0], QByteArray("chrom\tpos\tref\talt\tclinical_significance\treview_status\tgold_stars\t"
                                      "all_traits\tsymbol\tinheritance_modes\thgvs_p\thgvs_c\tmolecular_consequence\txrefs\r"));
        QCOMPARE(lines[1], QByteArray("17\t41276045\tA\tG\tPathogenic\tpractice guideline\t4\tBreast cancer\tBRCA1\tgermline\t\t\tmissense\tOMIM:1\r"));
        // The quoted input field is unquoted on read and re-quoted only where needed on write.
        QCOMPARE(lines[2], QByteArray("6\t26093141\tG\tA\tPathogenic; risk factor\tcriteria provided, single submitter\t1\t"
                                      "\"Trait with\ttab\"\tHFE\tgermline\t\t\t\t\r"));
    }
};

QTEST_GUILESS_MAIN(ClinVarProcessTest)
#include "tst_clinvarprocess.moc"
