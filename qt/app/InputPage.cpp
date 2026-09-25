#include "InputPage.h"

#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSettings>
#include <QVBoxLayout>

#include "GenomeFormat.h"
#include "RepoPaths.h"
#include "TestGenome.h"
#include "ToolLocator.h"

InputPage::InputPage(QWidget *parent) : QWidget(parent)
{
    auto *layout = new QVBoxLayout(this);
    auto *form = new QFormLayout;

    auto *genomeRow = new QHBoxLayout;
    m_genome = new QLineEdit;
    m_genome->setPlaceholderText("23andMe or AncestryDNA raw data, or a WGS-converted file");
    auto *browse = new QPushButton("Browse…");
    genomeRow->addWidget(m_genome, 1);
    genomeRow->addWidget(browse);
    form->addRow("Genome file", genomeRow);

    m_format = new QLabel;
    m_format->setWordWrap(true);
    form->addRow("Detected", m_format);

    m_name = new QLineEdit;
    m_name->setPlaceholderText("optional; appears in the report title");
    form->addRow("Subject name", m_name);

    m_dataDir = new QLineEdit(RepoPaths::dataDir());
    form->addRow("Data directory", m_dataDir);

    m_outDir = new QLineEdit(RepoPaths::root().isEmpty() ? QString() : RepoPaths::root() + "/reports");
    form->addRow("Output directory", m_outDir);

    layout->addLayout(form);

    auto *buttons = new QHBoxLayout;
    m_run = new QPushButton("Run analysis");
    m_run->setEnabled(false);
    buttons->addWidget(m_run);
    buttons->addStretch();
    layout->addLayout(buttons);

    auto *testRow = new QHBoxLayout;
    m_testGenome = new QPushButton("Get a public test genome…");
    m_testSample = new QComboBox;
    m_testSample->addItem("1000 Genomes NA12878", "NA12878");
    m_testSample->addItem("GRCh37 reference (null test)", "GRCh37");
    m_testSample->addItem("GRCh38 reference (null test)", "GRCh38");
    m_testSample->setToolTip("NA12878 is one openly consented 1000 Genomes sample, extracted at the positions "
                             "the analysis reads with bcftools range requests (a few MB). A reference build is its "
                             "own base at every position, homozygous: every variant the report then claims is a "
                             "statement about the reference genome, and GRCh37's mtDNA is rCRS (haplogroup H).");
    m_testStop = new QPushButton("Stop");
    m_testStop->setEnabled(false);
    m_testStop->setToolTip("Stops the build. What has been fetched so far is kept under "
                           "data/test_genome_cache/, so pressing Get again resumes rather than restarts.");
    m_testStatus = new QLabel;
    m_testStatus->setWordWrap(true);
    testRow->addWidget(m_testGenome);
    testRow->addWidget(m_testSample);
    testRow->addWidget(m_testStop);
    testRow->addWidget(m_testStatus, 1);
    layout->addLayout(testRow);

    auto *hint = new QLabel(
        "The analysis reads array-style genotype files directly. Sequencing reads "
        "(FASTQ) go through the WGS pipeline first, which aligns, calls variants and "
        "converts to the same layout.");
    hint->setWordWrap(true);
    hint->setStyleSheet("color: palette(mid)");
    layout->addWidget(hint);
    layout->addStretch();

    connect(browse, &QPushButton::clicked, this, &InputPage::browse);
    connect(m_genome, &QLineEdit::textChanged, this, &InputPage::onPathChanged);
    connect(m_run, &QPushButton::clicked, this, &InputPage::onRun);
    connect(m_testGenome, &QPushButton::clicked, this, &InputPage::buildTestGenome);
    connect(m_testStop, &QPushButton::clicked, this, [this] {
        m_job.cancel();
        m_testStop->setEnabled(false);
        m_testStatus->setText("Stopping…");
    });
    connect(&m_job, &JobRunner::stage, this, [this](const QString &s, int pct) {
        m_testStatus->setText(QStringLiteral("%1 (%2%)").arg(s).arg(pct));
    });
    connect(&m_job, &JobRunner::finished, this, [this](bool ok, const QString &err) {
        m_testGenome->setEnabled(true);
        m_testSample->setEnabled(true);
        m_testStop->setEnabled(false);
        if (!ok && err == "cancelled") {
            m_testStatus->setText("Stopped. What was fetched is kept; press Get again to resume.");
            return;
        }
        if (ok) {
            const QString sample = m_testSample->currentData().toString();
            const QString path = m_dataDir->text().trimmed()
                + (TestGenome::isReferenceSample(sample) ? "/genome_reference_" : "/genome_1000g_") + sample + ".txt";
            m_testStatus->setText(QStringLiteral("Ready: %1").arg(path));
            setGenomePath(path);
        } else {
            m_testStatus->setText(QStringLiteral("<span style=\"color:palette(highlight)\">%1</span>").arg(err));
        }
    });

    restoreSettings();
}

void InputPage::browse()
{
    const QString start = m_genome->text().isEmpty() ? RepoPaths::dataDir() : m_genome->text();
    const QString path = QFileDialog::getOpenFileName(
        this, "Choose a genome file", start,
        "Genome files (*.txt *.tsv *.csv *.vcf *.fastq *.fq *.gz);;All files (*)");
    if (!path.isEmpty())
        setGenomePath(path);
}

void InputPage::setGenomePath(const QString &path)
{
    m_genome->setText(path);
}

void InputPage::onPathChanged()
{
    const QString path = m_genome->text().trimmed();
    m_isFastq = false;
    if (path.isEmpty() || !QFileInfo(path).isFile()) {
        m_format->setText(path.isEmpty() ? QString() : "file not found");
        m_run->setEnabled(false);
        m_run->setText("Run analysis");
        return;
    }
    const GenomeSniff sniff = sniffGenome(path);
    QString text = QStringLiteral("<b>%1</b> — %2").arg(formatName(sniff.format), sniff.description);
    if (isArrayFormat(sniff.format)) {
        if (path.contains("wgs_work_validate"))
            text += "<br><span style=\"color:palette(highlight)\">This is the chr22 toolchain-validation "
                    "output: synthetic test loci on one chromosome, not a genome. It proves the "
                    "pipeline works; a report on it will be almost entirely \"unknown\".</span>";
        else if (sniff.dataRows < 1000 && sniff.sampledLines < 4000)
            text += "<br><span style=\"color:palette(highlight)\">Very few genotypes — most results "
                    "will come out as unknown. A 23andMe or AncestryDNA export has ~600,000 rows.</span>";
    }
    m_format->setText(text);
    m_isFastq = sniff.format == GenomeFormat::Fastq;
    // Sequencing reads go through the WGS pipeline first; Run does that
    // (and offers to fetch the reference) rather than sending the user to
    // another page.
    m_run->setEnabled(isArrayFormat(sniff.format) || m_isFastq);
    m_run->setText(m_isFastq ? "Run WGS pipeline + analysis" : "Run analysis");
}

void InputPage::onRun()
{
    if (m_isFastq) {
        saveSettings();
        emit wgsRunRequested(m_genome->text().trimmed(), m_name->text().trimmed());
        return;
    }
    AnalysisRunner::Request req;
    req.genomePath = m_genome->text().trimmed();
    req.dataDir = m_dataDir->text().trimmed();
    req.subjectName = m_name->text().trimmed();
    req.outputDir = m_outDir->text().trimmed();
    saveSettings();
    emit runRequested(req);
}

void InputPage::restoreSettings()
{
    QSettings s;
    m_genome->setText(s.value("input/genome").toString());
    m_name->setText(s.value("input/name").toString());
    if (!s.value("input/outDir").toString().isEmpty())
        m_outDir->setText(s.value("input/outDir").toString());
}

void InputPage::saveSettings()
{
    QSettings s;
    s.setValue("input/genome", m_genome->text());
    s.setValue("input/name", m_name->text());
    s.setValue("input/outDir", m_outDir->text());
}

void InputPage::buildTestGenome()
{
    if (m_job.isRunning())
        return;
    const QString sample = m_testSample->currentData().toString();
    if (!TestGenome::isReferenceSample(sample) && ToolLocator::find("bcftools").isEmpty()) {
        m_testStatus->setText(QStringLiteral(
            "<span style=\"color:palette(highlight)\"><b>Requirements not installed.</b></span> "
            "Install them with: <code>%1</code>").arg(ToolLocator::installHint()));
        return;
    }
    TestGenome::Options opt;
    opt.dataDir = m_dataDir->text().trimmed();
    opt.analysisBinary = RepoPaths::analysisBinary();
    opt.sample = sample;
    m_testGenome->setEnabled(false);
    m_testSample->setEnabled(false);
    m_testStop->setEnabled(true);
    m_testStatus->setText("Starting…");
    m_job.start([opt](const Reporter &r, QString *e) {
        TestGenome::Stats stats;
        return TestGenome::build(opt, r, &stats, e);
    });
}
