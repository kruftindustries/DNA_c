#include "DataSourcesPage.h"
#include "TextTable.h"

#include <QFileInfo>
#include <QProcess>
#include <QDirIterator>
#include <QDir>
#include <QFontDatabase>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QJsonObject>
#include <QLabel>
#include <QLocale>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>

#include "ClinVarUpdater.h"
#include "DataVersions.h"
#include "EnsemblLookup.h"
#include "ReferenceGenome.h"
#include "ToolLocator.h"
#include "PharmgkbUpdater.h"
#include "RepoPaths.h"

namespace {

struct SourceRow {
    const char *key;
    const char *label;
    const char *files[2];
};

const SourceRow kSources[] = {
    {"clinvar", "ClinVar", {"clinvar_alleles.tsv", nullptr}},
    {"pharmgkb", "ClinPGx (PharmGKB)", {"clinical_annotations.tsv", "clinical_ann_alleles.tsv"}},
};

} // namespace

DataSourcesPage::DataSourcesPage(QWidget *parent) : QWidget(parent)
{
    auto *layout = new QVBoxLayout(this);

    m_table = new QTableWidget(0, 4);
    m_table->setHorizontalHeaderLabels({"Source", "Release", "Last updated", "Files"});
    m_table->horizontalHeader()->setStretchLastSection(true);
    m_table->verticalHeader()->setVisible(false);
    m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_table->setSelectionMode(QAbstractItemView::NoSelection);
    m_table->setMaximumHeight(210);
    layout->addWidget(m_table);

    auto *bar = new QHBoxLayout;
    m_clinvar = new QPushButton("Update ClinVar");
    m_clinvar->setToolTip("Download variant_summary.txt.gz from NCBI (~430 MB) and reduce it to the GRCh37 variants (a few minutes)");
    m_pharmgkb = new QPushButton("Update ClinPGx");
    m_validate = new QPushButton("Validate ClinPGx");
    m_rsids = new QPushButton("Update rsID positions");
    m_rsids->setToolTip("Look up the GRCh37 position of every rsID the analysis reads (Ensembl); "
                        "used by the WGS pipeline and the test genomes");
    m_reference = new QPushButton("Set up GRCh37 reference…");
    m_reference->setToolTip("Download and index the reference genome the WGS pipeline aligns against "
                            "(~900 MB download, ~12 GB in reference/). Needs samtools and minimap2.");
    m_clearCaches = new QPushButton("Clear caches");
    m_clearCaches->setToolTip("Remove the test-genome, 1000 Genomes index and Ensembl caches under data/. "
                              "They are rebuilt on demand; the annotation data and genomes are kept.");
    auto *refreshButton = new QPushButton("Refresh");
    m_cancel = new QPushButton("Cancel");
    m_cancel->setEnabled(false);
    bar->addWidget(m_clinvar);
    bar->addWidget(m_pharmgkb);
    bar->addWidget(m_validate);
    bar->addWidget(m_rsids);
    bar->addWidget(m_reference);
    bar->addWidget(m_clearCaches);
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
    m_log->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    layout->addWidget(m_log, 1);

    connect(m_clinvar, &QPushButton::clicked, this, [this] {
        const QString dir = RepoPaths::dataDir();
        runJob("Updating ClinVar", [dir](const Reporter &r, QString *e) { return ClinVarUpdater::update(dir, r, e); });
    });
    connect(m_pharmgkb, &QPushButton::clicked, this, [this] {
        const QString dir = RepoPaths::dataDir();
        runJob("Updating ClinPGx", [dir](const Reporter &r, QString *e) { return PharmgkbUpdater::update(dir, r, e); });
    });
    connect(m_validate, &QPushButton::clicked, this, [this] {
        const QString dir = RepoPaths::dataDir();
        runJob("Validating ClinPGx", [dir](const Reporter &r, QString *e) {
            QStringList problems;
            const bool ok = PharmgkbUpdater::validate(dir, r, &problems);
            if (!ok) *e = problems.join("; ");
            return ok;
        });
    });
    connect(refreshButton, &QPushButton::clicked, this, &DataSourcesPage::refresh);
    connect(m_rsids, &QPushButton::clicked, this, [this] {
        const QString dataDir = RepoPaths::dataDir();
        const QString binary = RepoPaths::analysisBinary();
        runJob("Updating rsID positions", [dataDir, binary](const Reporter &r, QString *e) {
            QProcess list;
            list.start(binary, {"--list-rsids"});
            list.waitForFinished(-1);
            if (list.exitCode() != 0) {
                if (e) *e = QStringLiteral("could not run %1 --list-rsids").arg(binary);
                return false;
            }
            const QStringList rsids = outputLines(QString::fromUtf8(list.readAllStandardOutput()));
            QStringList unresolved;
            const bool ok = EnsemblLookup::update(dataDir + "/rsid_positions_grch37.json", rsids, r, &unresolved);
            if (!unresolved.isEmpty())
                r.log(QStringLiteral("  %1 rsID(s) have no GRCh37 position: %2").arg(unresolved.size()).arg(unresolved.join(", ")));
            return ok;
        });
    });
    connect(m_reference, &QPushButton::clicked, this, [this] {
        const QString root = RepoPaths::root();
        runJob("Setting up the GRCh37 reference", [root](const Reporter &r, QString *e) {
            return ReferenceGenome::setup(root, r, e);
        });
    });
    connect(m_clearCaches, &QPushButton::clicked, this, [this] {
        const QString dataDir = RepoPaths::dataDir();
        qint64 freed = 0;
        for (const char *sub : {"test_genome_cache", "1000g_index", "ensembl"}) {
            QDir d(dataDir + '/' + QLatin1String(sub));
            if (!d.exists())
                continue;
            QDirIterator it(d.path(), QDir::Files, QDirIterator::Subdirectories);
            while (it.hasNext()) { it.next(); freed += it.fileInfo().size(); }
            d.removeRecursively();
        }
        m_log->appendPlainText(QStringLiteral("Removed caches (%1)").arg(QLocale().formattedDataSize(freed)));
        refresh();
    });
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
    });

    refresh();
}

void DataSourcesPage::refresh()
{
    const QString dataDir = RepoPaths::dataDir();
    const QJsonObject versions = DataVersions::load(dataDir);
    m_table->setRowCount(0);
    for (const SourceRow &src : kSources) {
        const int row = m_table->rowCount();
        m_table->insertRow(row);
        const QJsonObject v = versions.value(QLatin1String(src.key)).toObject();
        QString release = v.value("release").toString();
        if (release.isEmpty() && v.contains("grch37_variants_written"))
            release = QStringLiteral("%L1 GRCh37 variants").arg(v.value("grch37_variants_written").toDouble(), 0, 'f', 0);
        const QString updated = v.value("updated").toString().left(16).replace('T', ' ');
        QStringList files;
        for (const char *name : src.files) {
            if (!name)
                continue;
            const QFileInfo fi(dataDir + '/' + name);
            files << (fi.exists() ? QStringLiteral("%1 (%2)").arg(name, QLocale().formattedDataSize(fi.size()))
                                  : QStringLiteral("%1 — missing").arg(name));
        }
        m_table->setItem(row, 0, new QTableWidgetItem(QLatin1String(src.label)));
        m_table->setItem(row, 1, new QTableWidgetItem(release.isEmpty() ? "—" : release));
        m_table->setItem(row, 2, new QTableWidgetItem(updated.isEmpty() ? "never" : updated));
        m_table->setItem(row, 3, new QTableWidgetItem(files.join(", ")));
    }

    // The derived data the app fetches on demand.
    auto addRow = [this](const QString &label, const QString &release, const QString &updated, const QString &files) {
        const int row = m_table->rowCount();
        m_table->insertRow(row);
        m_table->setItem(row, 0, new QTableWidgetItem(label));
        m_table->setItem(row, 1, new QTableWidgetItem(release));
        m_table->setItem(row, 2, new QTableWidgetItem(updated));
        m_table->setItem(row, 3, new QTableWidgetItem(files));
    };
    auto stamp = [](const QFileInfo &fi) {
        return fi.exists() ? fi.lastModified().toString("yyyy-MM-dd hh:mm") : QStringLiteral("never");
    };
    {
        const QFileInfo fi(dataDir + "/rsid_positions_grch37.json");
        const QJsonObject lookup = EnsemblLookup::loadLookup(fi.filePath());
        addRow("rsID positions (GRCh37)",
               fi.exists() ? QStringLiteral("%1 rsIDs").arg(lookup.size()) : "—", stamp(fi),
               fi.exists() ? QStringLiteral("rsid_positions_grch37.json (%1)").arg(QLocale().formattedDataSize(fi.size()))
                           : "rsid_positions_grch37.json — missing");
    }
    {
        const ReferenceGenome::Status ref = ReferenceGenome::status(RepoPaths::root());
        QStringList parts;
        parts << (ref.hasFasta ? QStringLiteral("FASTA %1").arg(QLocale().formattedDataSize(ref.fastaBytes)) : "FASTA missing")
              << (ref.hasFai ? "faidx" : "faidx missing")
              << (ref.hasMmi ? QStringLiteral("minimap2 index %1").arg(QLocale().formattedDataSize(ref.mmiBytes)) : "minimap2 index missing");
        addRow("GRCh37 reference (WGS pipeline)", ref.ready() ? "GRCh37 / hs37d5-free b37" : (ref.hasFasta ? "incomplete" : "—"),
               stamp(QFileInfo(ref.fasta)), parts.join(", "));
    }
    {
        qint64 bytes = 0;
        QStringList present;
        for (const char *sub : {"test_genome_cache", "1000g_index", "ensembl"}) {
            QDir d(dataDir + '/' + QLatin1String(sub));
            if (!d.exists())
                continue;
            qint64 b = 0;
            QDirIterator it(d.path(), QDir::Files, QDirIterator::Subdirectories);
            while (it.hasNext()) { it.next(); b += it.fileInfo().size(); }
            bytes += b;
            present << QStringLiteral("%1/ (%2)").arg(QLatin1String(sub), QLocale().formattedDataSize(b));
        }
        addRow("Caches (test genomes, indexes, Ensembl)", present.isEmpty() ? "—" : QLocale().formattedDataSize(bytes),
               "", present.isEmpty() ? "none" : present.join(", "));
    }
    {
        QStringList found, absent;
        for (const QString &t : ToolLocator::wgsTools()) {
            const QString p = ToolLocator::find(t);
            if (p.isEmpty())
                absent << (t == "fastp" ? QStringLiteral("fastp (optional)") : t);
            else
                found << QStringLiteral("%1: %2").arg(t, QDir::toNativeSeparators(p));
        }
        const QStringList required = ToolLocator::missing(ToolLocator::requiredWgsTools());
        addRow("Sequencing tools (FASTQ pipeline, 1000 Genomes sample)",
               required.isEmpty() ? QStringLiteral("ready") : QStringLiteral("missing: %1").arg(required.join(", ")),
               "",
               (found + absent).join("; ") + (required.isEmpty() ? QString() : ". " + ToolLocator::installHint()));
    }
    m_table->resizeColumnsToContents();
    m_reference->setEnabled(!m_job.isRunning() && !ReferenceGenome::status(RepoPaths::root()).ready());
}

void DataSourcesPage::runJob(const QString &title, const JobRunner::Job &job)
{
    if (m_job.isRunning())
        return;
    m_log->clear();
    m_stage->setText(title);
    m_progress->setValue(0);
    setBusy(true);
    m_job.start(job);
}

void DataSourcesPage::setBusy(bool busy)
{
    m_clinvar->setEnabled(!busy);
    m_pharmgkb->setEnabled(!busy);
    m_validate->setEnabled(!busy);
    m_rsids->setEnabled(!busy);
    m_reference->setEnabled(!busy);
    m_clearCaches->setEnabled(!busy);
    m_cancel->setEnabled(busy);
}
