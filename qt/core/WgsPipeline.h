// FASTQ -> genotype file -> report, natively.
//
// The steps are the Python pipeline's: optional fastp QC, minimap2 into
// samtools sort, samtools index, an rsID position lookup, targeted bcftools
// calling over the ClinVar and rsID positions, VCF conversion, then the C
// analysis binary. Each external step is a QProcess; the two pipes use
// QProcess's own chaining rather than a shell.
#pragma once

#include <QProcess>
#include <QString>
#include <QStringList>

#include "Reporter.h"

struct WgsConfig {
    QString fastq;
    QString subjectName;
    int threads = 0;              // 0 = all cores
    bool skipQc = false;
    bool full = false;            // genome-wide calling instead of targeted
    bool keepIntermediates = false;
    bool skipAnalysis = false;
    bool lookupRsids = true;      // query Ensembl for positions not cached
    bool backupExisting = true;   // copy an existing output to .bak first

    // Locations. defaults() fills these from the repository root and the
    // same environment overrides the Python honours (GH_REF_DIR, GH_WGS_WORK,
    // GH_REF_GENOME, GH_REF_MMI, GH_RSID_LOOKUP).
    QString root, dataDir, reportsDir, workDir;
    QString refGenome, refMmi, rsidLookup, outputGenome, analysisBinary;

    static WgsConfig defaults(const QString &root);
};

struct WgsOutcome {
    QString bam, vcf, genome;
    qint64 variants = 0;
};

class WgsPipeline {
public:
    WgsPipeline(const WgsConfig &config, const Reporter &reporter);

    // Everything the Python's check_prereqs checked, as a list of problems.
    QStringList prerequisites() const;

    bool run(WgsOutcome *outcome, QString *error);

    // Individual steps, public so the validator can run a subset.
    bool alignAndSort(const QString &fastq, QString *bam);
    bool callVariants(const QString &bam, QString *vcf);
    qint64 convert(const QString &vcf, const QString &genomeOut);

private:
    bool qc(const QString &fastq, QString *trimmed);
    bool buildTargets(QString *bed);
    bool analyse(const QString &genome);
    bool runProcess(const QString &program, const QStringList &args,
                    QString *capturedStdout = nullptr, bool logStdout = true);
    bool runPipe(const QString &program1, const QStringList &args1,
                 const QString &program2, const QStringList &args2);
    bool waitFor(QProcess &p);
    QString fail(const QString &what);

    WgsConfig m_cfg;
    Reporter m_rep;
    QString m_error;
};
