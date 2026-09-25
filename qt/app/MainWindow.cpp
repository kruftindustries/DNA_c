#include "MainWindow.h"

#include <QAction>
#include <QApplication>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QJsonDocument>
#include <QFormLayout>
#include <QLineEdit>
#include <QListWidget>
#include <QMenuBar>
#include <QMessageBox>
#include <QSplitter>
#include <QStackedWidget>
#include <QStatusBar>

#include "DataSourcesPage.h"
#include "FindingsPage.h"
#include "InputPage.h"
#include "RepoPaths.h"
#include "ReportPage.h"
#include "RunPage.h"
#include "WgsPage.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
    setWindowTitle("Genetic Health");
    resize(1100, 720);

    m_nav = new QListWidget;
    m_nav->addItems({"Input", "Run", "Report", "Findings", "Data Sources", "WGS Pipeline"});
    m_nav->setFixedWidth(150);

    m_stack = new QStackedWidget;
    m_input = new InputPage;
    m_run = new RunPage;
    m_report = new ReportPage;
    m_findings = new FindingsPage;
    m_data = new DataSourcesPage;
    m_wgs = new WgsPage;
    for (QWidget *w : {static_cast<QWidget *>(m_input), static_cast<QWidget *>(m_run),
                       static_cast<QWidget *>(m_report), static_cast<QWidget *>(m_findings),
                       static_cast<QWidget *>(m_data), static_cast<QWidget *>(m_wgs)})
        m_stack->addWidget(w);

    auto *split = new QSplitter;
    split->addWidget(m_nav);
    split->addWidget(m_stack);
    split->setStretchFactor(1, 1);
    setCentralWidget(split);

    connect(m_nav, &QListWidget::currentRowChanged, m_stack, &QStackedWidget::setCurrentIndex);
    m_nav->setCurrentRow(Input);

    connect(m_input, &InputPage::runRequested, this, &MainWindow::startAnalysis);
    connect(m_input, &InputPage::wgsRunRequested, this, [this](const QString &p, const QString &name) {
        showPage(Wgs);
        m_wgs->startForFile(p, name);
    });
    connect(m_wgs, &WgsPage::reportReady, this, [this](const QString &html, const QString &json) {
        QFile jf(json);
        if (jf.open(QIODevice::ReadOnly))
            m_findings->setResults(QJsonDocument::fromJson(jf.readAll()).object());
        m_report->load(html);
        showPage(Report);
    });
    connect(m_run, &RunPage::cancelRequested, &m_runner, &AnalysisRunner::cancel);
    connect(&m_runner, &AnalysisRunner::logLine, m_run, &RunPage::appendLine);
    connect(&m_runner, &AnalysisRunner::stageChanged, m_run, &RunPage::setStage);
    connect(&m_runner, &AnalysisRunner::finished, this, [this](bool ok, const QString &msg) {
        m_run->setRunning(false);
        statusBar()->showMessage(msg, 10000);
        m_run->appendLine(msg);
        if (!ok)
            return;
        m_findings->setResults(m_runner.results());
        m_report->load(m_runner.reportPath());
        showPage(Report);
    });

    QMenu *file = menuBar()->addMenu("&File");
    file->addAction("&Open genome…", m_input, &InputPage::browse, QKeySequence::Open);
    file->addAction("&Settings…", this, &MainWindow::openSettings);
    file->addSeparator();
    file->addAction("&Quit", qApp, &QApplication::quit, QKeySequence::Quit);
    menuBar()->addMenu("&Help")->addAction("&About", this, &MainWindow::about);

    statusBar()->showMessage(RepoPaths::analysisBinary().isEmpty()
                                 ? "Analysis binary not found — build it with `make -C c` or set it in Settings"
                                 : "Ready");
}

void MainWindow::showPage(Page page)
{
    m_nav->setCurrentRow(page);
}

void MainWindow::startAnalysis(const AnalysisRunner::Request &req)
{
    if (m_runner.isRunning())
        return;
    m_run->clear();
    m_run->setRunning(true);
    showPage(Run);
    m_runner.start(req);
}

void MainWindow::openSettings()
{
    QDialog dlg(this);
    dlg.setWindowTitle("Settings");
    auto *form = new QFormLayout(&dlg);
    auto *data = new QLineEdit(RepoPaths::dataDir());
    auto *binary = new QLineEdit(RepoPaths::analysisBinary());
    form->addRow("Data directory", data);
    form->addRow("Analysis binary", binary);
    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    form->addRow(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);
    if (dlg.exec() != QDialog::Accepted)
        return;
    RepoPaths::setDataDir(data->text().trimmed());
    RepoPaths::setAnalysisBinary(binary->text().trimmed());
    m_data->refresh();
    m_wgs->refresh();
}

void MainWindow::about()
{
    QMessageBox::about(this, "About Genetic Health",
        QStringLiteral("<b>Genetic Health</b><br>Qt %1<br>Report viewer: %2<br>"
                       "Repository: %3<br>Analysis binary: %4")
            .arg(QT_VERSION_STR,
#ifdef HAVE_WEBENGINE
                 "Qt WebEngine",
#else
                 "basic HTML (Qt WebEngine not installed)",
#endif
                 RepoPaths::root().isEmpty() ? "not found" : RepoPaths::root(),
                 RepoPaths::analysisBinary().isEmpty() ? "not found" : RepoPaths::analysisBinary()));
}
