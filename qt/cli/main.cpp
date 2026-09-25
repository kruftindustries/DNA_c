// gh-data: the data updater and WGS pipeline as a console tool.
//
// The same core the GUI uses, driven from the command line, so the two
// Python modules it replaces (update_data.py, wgs_pipeline.py) have a
// scriptable successor and the parity tests have something to run.
#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QJsonObject>
#include <QLocale>
#include <QThread>

#include <cstdio>

#include "ClinVarUpdater.h"
#include "DataVersions.h"
#include "EnsemblLookup.h"
#include "PharmgkbUpdater.h"
#include "Reporter.h"
#include "TargetRegions.h"
#include "TestGenome.h"
#include "ReferenceGenome.h"
#include "ToolLocator.h"
#include "VcfConverter.h"
#include "WgsPipeline.h"
#include "WgsValidator.h"

namespace {

int usage()
{
    std::fputs(
        "gh-data - annotation data updates and the WGS pipeline\n"
        "\n"
        "Usage:\n"
        "  gh-data clinvar   [--data DIR]          download + process ClinVar\n"
        "  gh-data pharmgkb  [--data DIR]          download ClinPGx annotations\n"
        "  gh-data validate-pharmgkb [--data DIR]  check the installed ClinPGx files\n"
        "  gh-data status    [--data DIR]          show installed releases\n"
        "  gh-data clinvar-process IN.gz OUT.tsv   the ClinVar filter alone (for tests)\n"
        "  gh-data wgs FASTQ [--name N] [--threads T] [--skip-qc] [--full]\n"
        "                    [--skip-analysis] [--keep-intermediates] [--output PATH]\n"
        "  gh-data reference [--status]           fetch + index the GRCh37 reference (~12 GB)\n"
        "  gh-data wgs-validate [--setup] [--clinvar-targets N] [--threads T]\n"
        "  gh-data build-targets CLINVAR.tsv LOOKUP.json OUT.bed\n"
        "  gh-data convert-vcf VCF OUT.txt [LOOKUP.json]\n"
        "  gh-data lookup-rsids LOOKUP.json rsID...\n"
        "  gh-data test-genome [--sample NA12878|GRCh37|GRCh38] [--output PATH]\n"
        "                    a public 1000 Genomes genome at the analysed positions, or a\n"
        "                    reference assembly's own bases as a null-test sample\n"
        "  gh-data tools                           where the WGS tools were found\n"
        "\n"
        "The repository root is found from the executable's location; --root DIR overrides it.\n",
        stderr);
    return 2;
}

QString findRoot()
{
    QString dir = QCoreApplication::applicationDirPath();
    for (int i = 0; i < 8; ++i) {
        if (QFileInfo(dir + "/c/src/main.c").exists())
            return QDir(dir).absolutePath();
        QDir d(dir);
        if (!d.cdUp())
            break;
        dir = d.absolutePath();
    }
    return QDir::currentPath();
}

QString takeOption(QStringList &args, const QString &name, const QString &fallback = QString())
{
    const int i = args.indexOf(name);
    if (i < 0 || i + 1 >= args.size())
        return fallback;
    const QString v = args[i + 1];
    args.removeAt(i + 1);
    args.removeAt(i);
    return v;
}

bool takeFlag(QStringList &args, const QString &name)
{
    return args.removeAll(name) > 0;
}

int status(const QString &dataDir)
{
    const QJsonObject v = DataVersions::load(dataDir);
    std::printf("Data directory: %s\n", qPrintable(dataDir));
    const QFileInfo cv(dataDir + "/clinvar_alleles.tsv");
    std::printf("\nClinVar (clinvar_alleles.tsv): %s\n", cv.exists()
        ? qPrintable(QStringLiteral("EXISTS (%1)").arg(QLocale().formattedDataSize(cv.size()))) : "MISSING");
    if (v.contains("clinvar")) {
        const QJsonObject c = v.value("clinvar").toObject();
        std::printf("    Last updated: %s\n    GRCh37 variants: %lld\n", qPrintable(c.value("updated").toString()),
                    qlonglong(c.value("grch37_variants_written").toDouble()));
    }
    std::printf("\nClinPGx:\n");
    for (const char *f : {"clinical_annotations.tsv", "clinical_ann_alleles.tsv"}) {
        const QFileInfo fi(dataDir + '/' + f);
        std::printf("    %s: %s\n", f, fi.exists()
            ? qPrintable(QStringLiteral("EXISTS (%1)").arg(QLocale().formattedDataSize(fi.size()))) : "MISSING");
    }
    if (v.contains("pharmgkb")) {
        const QJsonObject p = v.value("pharmgkb").toObject();
        std::printf("    Release: %s\n    Last updated: %s\n    Last validated: %s\n",
                    qPrintable(p.value("release").toString("unknown")),
                    qPrintable(p.value("updated").toString("unknown")),
                    qPrintable(p.value("validated").toString("unknown")));
    }
    return 0;
}

} // namespace

int main(int argc, char **argv)
{
    QCoreApplication app(argc, argv);
    QStringList args = app.arguments().mid(1);
    if (args.isEmpty())
        return usage();

    const QString root = takeOption(args, "--root", findRoot());
    const QString dataDir = takeOption(args, "--data", root + "/data");
    const Reporter rep = stdoutReporter();
    const QString cmd = args.takeFirst();
    QString error;

    if (cmd == "clinvar") {
        if (ClinVarUpdater::update(dataDir, rep, &error))
            return 0;
        std::fprintf(stderr, "error: %s\n", qPrintable(error));
        return 1;
    }
    if (cmd == "pharmgkb") {
        if (PharmgkbUpdater::update(dataDir, rep, &error))
            return 0;
        std::fprintf(stderr, "error: %s\n", qPrintable(error));
        return 1;
    }
    if (cmd == "validate-pharmgkb") {
        QStringList problems;
        return PharmgkbUpdater::validate(dataDir, rep, &problems) ? 0 : 1;
    }
    if (cmd == "status")
        return status(dataDir);
    if (cmd == "clinvar-process") {
        if (args.size() != 2)
            return usage();
        ClinVarUpdater::Stats stats;
        if (!ClinVarUpdater::process(args[0], args[1], &stats, rep, &error)) {
            std::fprintf(stderr, "error: %s\n", qPrintable(error));
            return 1;
        }
        std::printf("processed %lld rows, wrote %lld\n", qlonglong(stats.rowsProcessed), qlonglong(stats.written));
        return 0;
    }
    if (cmd == "wgs") {
        WgsConfig cfg = WgsConfig::defaults(root);
        cfg.dataDir = dataDir;
        cfg.subjectName = takeOption(args, "--name");
        cfg.threads = takeOption(args, "--threads", "0").toInt();
        cfg.outputGenome = takeOption(args, "--output", cfg.outputGenome);
        cfg.skipQc = takeFlag(args, "--skip-qc");
        cfg.full = takeFlag(args, "--full");
        cfg.skipAnalysis = takeFlag(args, "--skip-analysis");
        cfg.keepIntermediates = takeFlag(args, "--keep-intermediates");
        if (args.size() != 1)
            return usage();
        cfg.fastq = QFileInfo(args[0]).absoluteFilePath();
        if (!QFileInfo(cfg.fastq).isFile()) {
            std::fprintf(stderr, "FASTQ file not found: %s\n", qPrintable(cfg.fastq));
            return 1;
        }
        std::printf("WGS Pipeline: FASTQ -> Genetic Health Reports\nInput:   %s (%.1f GB)\nThreads: %d\nMode:    %s\n",
                    qPrintable(QFileInfo(cfg.fastq).fileName()), QFileInfo(cfg.fastq).size() / 1e9,
                    cfg.threads > 0 ? cfg.threads : QThread::idealThreadCount(),
                    cfg.full ? "Full genome-wide variant calling" : "Targeted (ClinVar + SNP db positions only)");
        WgsPipeline pipeline(cfg, rep);
        WgsOutcome out;
        if (pipeline.run(&out, &error))
            return 0;
        std::fprintf(stderr, "error: %s\n", qPrintable(error));
        return 1;
    }
    if (cmd == "reference") {
        const bool statusOnly = takeFlag(args, "--status");
        if (!args.isEmpty())
            return usage();
        const ReferenceGenome::Status st = ReferenceGenome::status(root);
        std::printf("Reference directory: %s\n  %s: %s\n  %s: %s\n  %s: %s\n",
                    qPrintable(st.dir),
                    qPrintable(QFileInfo(st.fasta).fileName()), st.hasFasta ? "present" : "missing",
                    qPrintable(QFileInfo(st.fai).fileName()), st.hasFai ? "present" : "missing",
                    qPrintable(QFileInfo(st.mmi).fileName()), st.hasMmi ? "present" : "missing");
        if (statusOnly || st.ready())
            return st.ready() ? 0 : 1;
        if (!ReferenceGenome::setup(root, rep, &error)) {
            std::fprintf(stderr, "error: %s\n", qPrintable(error));
            return 1;
        }
        return 0;
    }
    if (cmd == "wgs-validate") {
        const bool doSetup = takeFlag(args, "--setup");
        WgsValidator::Options opt;
        opt.root = root;
        opt.clinvarTargets = takeOption(args, "--clinvar-targets", "150").toInt();
        opt.threads = takeOption(args, "--threads", "0").toInt();
        if (doSetup || !WgsValidator::isSetUp(root)) {
            if (!WgsValidator::setup(root, rep, &error)) {
                std::fprintf(stderr, "error: %s\n", qPrintable(error));
                return 1;
            }
            if (doSetup)
                return 0;
        }
        QStringList problems;
        if (WgsValidator::run(opt, rep, &problems, &error))
            return 0;
        for (const QString &p : problems.mid(0, 20))
            std::printf("  %s\n", qPrintable(p));
        if (!error.isEmpty())
            std::fprintf(stderr, "error: %s\n", qPrintable(error));
        return 1;
    }
    if (cmd == "build-targets") {
        if (args.size() != 3)
            return usage();
        return TargetRegions::build(args[0], args[1], args[2], rep, &error) >= 0 ? 0 : 1;
    }
    if (cmd == "convert-vcf") {
        if (args.size() < 2 || args.size() > 3)
            return usage();
        return VcfConverter::convert(args[0], args[1], args.value(2), rep, &error) >= 0 ? 0 : 1;
    }
    if (cmd == "lookup-rsids") {
        if (args.size() < 2)
            return usage();
        const QString path = args.takeFirst();
        QStringList unresolved;
        return EnsemblLookup::update(path, args, rep, &unresolved) ? 0 : 1;
    }
    if (cmd == "tools") {
        for (const QString &t : ToolLocator::wgsTools()) {
            const QString p = ToolLocator::find(t);
            std::printf("%-10s %s\n", qPrintable(t), p.isEmpty() ? "NOT FOUND" : qPrintable(p));
        }
        const QStringList absent = ToolLocator::missing(ToolLocator::requiredWgsTools());
        if (!absent.isEmpty())
            std::printf("\nRequirements not installed (%s). Install them with:\n  %s\n",
                        qPrintable(absent.join(", ")), qPrintable(ToolLocator::installHint()));
        return absent.isEmpty() ? 0 : 1;
    }
    if (cmd == "test-genome") {
        TestGenome::Options opt;
        opt.dataDir = dataDir;
        opt.analysisBinary = ToolLocator::analysisBinary(root);
        opt.sample = takeOption(args, "--sample", opt.sample);
        opt.output = takeOption(args, "--output");
        if (!args.isEmpty())
            return usage();
        TestGenome::Stats stats;
        if (TestGenome::build(opt, rep, &stats, &error))
            return 0;
        std::fprintf(stderr, "error: %s\n", qPrintable(error));
        return 1;
    }
    return usage();
}
