#include "WgsPage.h"

#include <QCheckBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontDatabase>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMessageBox>
#include <QLocale>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QSpinBox>
#include <QTableWidget>
#include <QThread>
#include <QVBoxLayout>

#include "ReferenceGenome.h"
#include "RepoPaths.h"
#include "ToolLocator.h"
#include "WgsPipeline.h"
#include "WgsValidator.h"

WgsPage::WgsPage(QWidget *parent) : QWidget(parent)
{
    auto *layout = new QVBoxLayout(this);

    m_tools = new QTableWidget(0, 2);
    m_tools->setHorizontalHeaderLabels({"Tool", "Found at"});
    m_tools->horizontalHeader()->setStretchLastSection(true);
    m_tools->verticalHeader()->setVisible(false);
    m_tools->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_tools->setSelectionMode(QAbstractItemView::NoSelection);
    m_tools->setMaximumHeight(150);
    layout->addWidget(m_tools);

    m_requirements = new QLabel;
    m_requirements->setWordWrap(true);
    m_requirements->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(m_requirements);

    m_reference = new QLabel;
    m_reference->setWordWrap(true);
    layout->addWidget(m_reference);

    auto *form = new QFormLayout;
    auto *fastqRow = new QHBoxLayout;
    m_fastq = new QLineEdit;
    auto *browse = new QPushButton("Browse…");
    fastqRow->addWidget(m_fastq, 1);
    fastqRow->addWidget(browse);
    form->addRow("FASTQ", fastqRow);
    m_name = new QLineEdit;
    form->addRow("Subject name", m_name);
    m_threads = new QSpinBox;
    m_threads->setRange(1, 256);
    m_threads->setValue(QThread::idealThreadCount());
    form->addRow("Threads", m_threads);
    auto *flags = new QHBoxLayout;
    m_skipQc = new QCheckBox("Skip fastp QC");
    m_full = new QCheckBox("Genome-wide calling (slow)");
    m_keep = new QCheckBox("Keep intermediates");
    flags->addWidget(m_skipQc);
    flags->addWidget(m_full);
    flags->addWidget(m_keep);
    flags->addStretch();
    form->addRow(QString(), flags);
    layout->addLayout(form);

    auto *bar = new QHBoxLayout;
    m_validate = new QPushButton("Validate toolchain (chr22)");
    m_setupRef = new QPushButton("Set up GRCh37 reference…");
    m_setupRef->setToolTip("Download human_g1k_v37.fasta.gz from the 1000 Genomes site (~900 MB), "
                           "decompress it and build the samtools and minimap2 indexes (~12 GB in "
                           "reference/, ~10-20 minutes). Resumes if interrupted.");
    m_run = new QPushButton("Run WGS pipeline");
    auto *refreshButton = new QPushButton("Refresh");
    m_cancel = new QPushButton("Cancel");
    m_cancel->setEnabled(false);
    bar->addWidget(m_validate);
    bar->addWidget(m_setupRef);
    bar->addWidget(m_run);
    bar->addWidget(refreshButton);
    bar->addStretch();
    bar->addWidget(m_cancel);
    layout->addLayout(bar);

    auto *prog = new QHBoxLayout;
    m_stage = new QLabel("Idle");
    m_progress = new QProgressBar;
    m_progress->setRange(0, 100);
    prog->addWidget(m_stage);
    prog->addWidget(m_progress, 1);
    layout->addLayout(prog);

    m_log = new QPlainTextEdit;
    m_log->setReadOnly(true);
    m_log->setMaximumBlockCount(50000);
    m_log->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    layout->addWidget(m_log, 1);

    connect(browse, &QPushButton::clicked, this, [this] {
        const QString p = QFileDialog::getOpenFileName(this, "Choose a FASTQ", m_fastq->text(),
                                                       "FASTQ (*.fastq *.fq *.fastq.gz *.fq.gz);;All files (*)");
        if (!p.isEmpty())
            m_fastq->setText(p);
    });
    connect(m_validate, &QPushButton::clicked, this, &WgsPage::validateToolchain);
    connect(m_run, &QPushButton::clicked, this, &WgsPage::runPipeline);
    connect(m_setupRef, &QPushButton::clicked, this, &WgsPage::setupReference);
    connect(refreshButton, &QPushButton::clicked, this, &WgsPage::refresh);
    connect(m_cancel, &QPushButton::clicked, &m_job, &JobRunner::cancel);
    connect(&m_job, &JobRunner::logLine, m_log, &QPlainTextEdit::appendPlainText);
    connect(&m_job, &JobRunner::stage, this, [this](const QString &s, int pct) {
        m_stage->setText(s);
        m_progress->setValue(pct);
    });
    connect(&m_job, &JobRunner::finished, this, [this](bool ok, const QString &error) {
        m_log->appendPlainText(ok ? "(done)" : QStringLiteral("(failed: %1)").arg(error));
        m_stage->setText(ok ? "Done" : "Failed");
        if (ok)
            m_progress->setValue(100);
        setBusy(false);
        refresh();
        if (ok && m_lastWasRun && !m_full->isChecked()) {
            const QString reports = RepoPaths::root() + "/reports";
            emit reportReady(reports + "/GENETIC_HEALTH_REPORT.html", reports + "/comprehensive_results.json");
        }
    });

    refresh();
}

void WgsPage::setFastqPath(const QString &path) { m_fastq->setText(path); }

void WgsPage::refresh()
{
    m_tools->setRowCount(0);
    bool allRequired = true;
    const QStringList tools = ToolLocator::wgsTools();
    for (const QString &t : tools) {
        const QString path = ToolLocator::find(t);
        const int row = m_tools->rowCount();
        m_tools->insertRow(row);
        m_tools->setItem(row, 0, new QTableWidgetItem(t));
        m_tools->setItem(row, 1, new QTableWidgetItem(
            path.isEmpty() ? (t == "fastp" ? "not found (QC step will be skipped)" : "NOT FOUND") : path));
        if (path.isEmpty() && t != "fastp")
            allRequired = false;
    }
    m_tools->resizeColumnsToContents();

    if (!allRequired)
        m_requirements->setText(QStringLiteral(
            "<span style=\"color:palette(highlight)\"><b>Requirements not installed.</b></span> "
            "Install them with:<br><code>%1</code>").arg(ToolLocator::installHint()));
    else
        m_requirements->clear();
    m_requirements->setVisible(!m_requirements->text().isEmpty());

    const QString root = RepoPaths::root();
    const ReferenceGenome::Status ref = ReferenceGenome::status(root);
    const bool chr22 = WgsValidator::isSetUp(root);
    QString refText;
    if (ref.ready())
        refText = QStringLiteral("GRCh37 present (%1, index %2)")
                      .arg(QLocale().formattedDataSize(ref.fastaBytes), QLocale().formattedDataSize(ref.mmiBytes));
    else if (ref.hasFasta)
        refText = QStringLiteral("GRCh37 downloaded, indexes missing — press <i>Set up GRCh37 reference</i> to finish");
    else
        refText = QStringLiteral("not set up — <i>Set up GRCh37 reference</i> fetches and indexes it (~12 GB)");
    m_reference->setText(QStringLiteral("Reference genome: %1. chr22 validation sandbox: %2.")
        .arg(refText, chr22 ? "ready" : "not set up (the validation button fetches it, ~10 MB)"));
    const bool busy = m_job.isRunning();
    m_run->setEnabled(allRequired && ref.ready() && !busy);
    m_run->setToolTip(!allRequired ? "minimap2, samtools and bcftools are all required"
                      : !ref.ready() ? QStringLiteral("the GRCh37 reference is needed first: %1").arg(ref.fasta)
                                     : QString());
    m_setupRef->setEnabled(allRequired && !ref.ready() && !busy);
    m_validate->setEnabled(allRequired && !busy);
}

void WgsPage::validateToolchain()
{
    const QString root = RepoPaths::root();
    if (root.isEmpty())
        return;
    m_lastWasRun = false;
    const int threads = m_threads->value();
    runJob("Validating toolchain", [root, threads](const Reporter &r, QString *e) {
        if (!WgsValidator::isSetUp(root) && !WgsValidator::setup(root, r, e))
            return false;
        WgsValidator::Options opt;
        opt.root = root;
        opt.threads = threads;
        QStringList problems;
        const bool ok = WgsValidator::run(opt, r, &problems, e);
        for (const QString &p : problems.mid(0, 20))
            r.log("  " + p);
        return ok;
    });
}

void WgsPage::runPipeline()
{
    startPipeline(false);
}

void WgsPage::setupReference()
{
    const QString root = RepoPaths::root();
    if (root.isEmpty())
        return;
    m_lastWasRun = false;
    runJob("Setting up the GRCh37 reference", [root](const Reporter &r, QString *e) {
        return ReferenceGenome::setup(root, r, e);
    });
}

void WgsPage::startForFile(const QString &fastq, const QString &name)
{
    m_fastq->setText(fastq);
    if (!name.isEmpty())
        m_name->setText(name);
    refresh();
    if (m_job.isRunning())
        return;
    if (!m_requirements->text().isEmpty()) {
        m_log->appendPlainText("The WGS tools are not installed; see above for the install command.");
        return;
    }
    if (ReferenceGenome::status(RepoPaths::root()).ready()) {
        startPipeline(false);
        return;
    }
    const auto answer = QMessageBox::question(
        this, "Reference genome needed",
        "Sequencing reads are aligned against the GRCh37 reference, which is not set up yet.\n\n"
        "Download it now (~900 MB, about 12 GB on disk once indexed, 10-20 minutes) and then run "
        "the pipeline on this file?",
        QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes);
    if (answer == QMessageBox::Yes)
        startPipeline(true);
}

void WgsPage::startPipeline(bool setupReferenceFirst)
{
    const QString fastq = m_fastq->text().trimmed();
    if (fastq.isEmpty() || !QFileInfo(fastq).isFile()) {
        m_log->appendPlainText("choose a FASTQ file first");
        return;
    }
    WgsConfig cfg = WgsConfig::defaults(RepoPaths::root());
    cfg.dataDir = RepoPaths::dataDir();
    cfg.analysisBinary = RepoPaths::analysisBinary();
    cfg.fastq = fastq;
    cfg.subjectName = m_name->text().trimmed();
    cfg.threads = m_threads->value();
    cfg.skipQc = m_skipQc->isChecked();
    cfg.full = m_full->isChecked();
    cfg.keepIntermediates = m_keep->isChecked();
    m_lastWasRun = true;
    runJob(setupReferenceFirst ? "Setting up the reference, then running the WGS pipeline"
                               : "Running WGS pipeline",
           [cfg, setupReferenceFirst](const Reporter &r, QString *e) {
        if (setupReferenceFirst && !ReferenceGenome::status(cfg.root).ready()
                && !ReferenceGenome::setup(cfg.root, r, e))
            return false;
        WgsPipeline pipeline(cfg, r);
        WgsOutcome out;
        return pipeline.run(&out, e);
    });
}

void WgsPage::runJob(const QString &title, const JobRunner::Job &job)
{
    if (m_job.isRunning())
        return;
    m_log->clear();
    m_stage->setText(title);
    m_progress->setValue(0);
    setBusy(true);
    m_job.start(job);
}

void WgsPage::setBusy(bool busy)
{
    m_validate->setEnabled(!busy);
    m_run->setEnabled(!busy);
    m_setupRef->setEnabled(!busy);
    m_cancel->setEnabled(busy);
    if (!busy)
        refresh();   // the reference may have just appeared
}
