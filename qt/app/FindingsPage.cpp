#include "FindingsPage.h"

#include <QApplication>
#include <QClipboard>
#include <QComboBox>
#include <QFile>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QJsonArray>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QTableView>
#include <QVBoxLayout>

FindingsPage::FindingsPage(QWidget *parent) : QWidget(parent)
{
    auto *layout = new QVBoxLayout(this);

    m_summary = new QLabel("No results yet.");
    layout->addWidget(m_summary);

    m_empty = new QLabel;
    m_empty->setWordWrap(true);
    m_empty->setStyleSheet("color: palette(mid); padding: 1em");
    layout->addWidget(m_empty);

    auto *bar = new QHBoxLayout;
    m_search = new QLineEdit;
    m_search->setPlaceholderText("Search gene, rsID, status, description…");
    m_search->setClearButtonEnabled(true);
    m_category = new QComboBox;
    m_category->addItem("All categories", QString());
    m_minMag = new QSpinBox;
    m_minMag->setRange(0, 6);
    m_minMag->setPrefix("magnitude ≥ ");
    m_count = new QLabel;
    auto *copy = new QPushButton("Copy TSV");
    auto *save = new QPushButton("Export TSV…");
    bar->addWidget(m_search, 1);
    bar->addWidget(m_category);
    bar->addWidget(m_minMag);
    bar->addWidget(m_count);
    bar->addWidget(copy);
    bar->addWidget(save);
    layout->addLayout(bar);

    m_model = new FindingsModel(this);
    m_filter = new FindingsFilter(this);
    m_filter->setSourceModel(m_model);
    m_table = new QTableView;
    m_table->setModel(m_filter);
    m_table->setSortingEnabled(true);
    m_table->sortByColumn(FindingsModel::Magnitude, Qt::DescendingOrder);
    m_table->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_table->setAlternatingRowColors(true);
    m_table->horizontalHeader()->setStretchLastSection(true);
    m_table->verticalHeader()->setVisible(false);
    m_table->setWordWrap(false);
    layout->addWidget(m_table, 1);

    connect(m_search, &QLineEdit::textChanged, this, [this](const QString &t) {
        m_filter->setSearchText(t); updateCount(); });
    connect(m_category, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int i) {
        m_filter->setCategory(m_category->itemData(i).toString()); updateCount(); });
    connect(m_minMag, QOverload<int>::of(&QSpinBox::valueChanged), this, [this](int v) {
        m_filter->setMinimumMagnitude(v); updateCount(); });
    connect(copy, &QPushButton::clicked, this, &FindingsPage::copyTsv);
    connect(save, &QPushButton::clicked, this, &FindingsPage::exportTsv);
    updateCount();
}

void FindingsPage::setResults(const QJsonObject &results)
{
    m_model->setFindings(results.value("findings").toArray());
    const QString current = m_category->currentData().toString();
    m_category->blockSignals(true);
    m_category->clear();
    m_category->addItem("All categories", QString());
    for (const QString &c : m_model->categories())
        m_category->addItem(c, c);
    const int idx = m_category->findData(current);
    m_category->setCurrentIndex(idx < 0 ? 0 : idx);
    m_category->blockSignals(false);
    m_filter->setCategory(m_category->currentData().toString());

    const QJsonObject s = results.value("summary").toObject();
    m_summary->setText(QStringLiteral("%1 variants in genome · %2 database SNPs matched · "
                                      "%3 high, %4 moderate, %5 low impact · %6 drug-gene interactions")
                           .arg(s.value("total_snps").toInt())
                           .arg(s.value("analyzed_snps").toInt())
                           .arg(s.value("high_impact").toInt())
                           .arg(s.value("moderate_impact").toInt())
                           .arg(s.value("low_impact").toInt())
                           .arg(results.value("drug_findings").toArray().size()));
    m_table->resizeColumnsToContents();
    m_haveResults = true;
    updateCount();
}

QVector<Finding> FindingsPage::visibleFindings() const
{
    QVector<Finding> out;
    for (int r = 0; r < m_filter->rowCount(); ++r) {
        const int src = m_filter->mapToSource(m_filter->index(r, 0)).row();
        out << m_model->findings()[src];
    }
    return out;
}

void FindingsPage::updateCount()
{
    m_count->setText(QStringLiteral("%1 of %2").arg(m_filter->rowCount()).arg(m_model->rowCount()));

    // The table explains itself when it has nothing to show.
    QString why;
    if (!m_haveResults)
        why = "This page fills in when an analysis finishes: choose a genome on the "
              "<b>Input</b> page and press <i>Run analysis</i>.";
    else if (m_model->rowCount() == 0)
        why = "The analysis ran, but none of this genome's positions matched the curated SNP "
              "database. That is what to expect from a very small or single-chromosome file "
              "(for example the chr22 toolchain-validation output, which is synthetic test "
              "data, not a genome). A 23andMe or AncestryDNA export, or a full WGS run, "
              "typically matches around 200 database SNPs.";
    else if (m_filter->rowCount() == 0)
        why = "No findings match the current filters.";
    m_empty->setText(why);
    m_empty->setVisible(!why.isEmpty());
}

void FindingsPage::copyTsv()
{
    QApplication::clipboard()->setText(FindingsModel::toTsv(visibleFindings()));
}

void FindingsPage::exportTsv()
{
    const QString to = QFileDialog::getSaveFileName(this, "Export findings", "findings.tsv", "TSV (*.tsv)");
    if (to.isEmpty())
        return;
    QFile f(to);
    if (f.open(QIODevice::WriteOnly | QIODevice::Truncate))
        f.write(FindingsModel::toTsv(visibleFindings()).toUtf8());
}
