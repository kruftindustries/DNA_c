#include "ReportPage.h"

#include <QDesktopServices>
#include <QFile>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QUrl>
#include <QVBoxLayout>
#ifdef HAVE_WEBENGINE
#include <QWebEngineView>
#else
#include <QTextBrowser>
#endif

ReportPage::ReportPage(QWidget *parent) : QWidget(parent)
{
    auto *layout = new QVBoxLayout(this);
    auto *bar = new QHBoxLayout;
    m_status = new QLabel("No report yet — run an analysis first.");
    auto *open = new QPushButton("Open in browser");
    auto *save = new QPushButton("Save copy…");
    bar->addWidget(m_status, 1);
    bar->addWidget(open);
    bar->addWidget(save);
    layout->addLayout(bar);

#ifdef HAVE_WEBENGINE
    m_view = new QWebEngineView;
#else
    m_view = new QTextBrowser;
    m_view->setOpenExternalLinks(true);
    auto *note = new QLabel(
        "Qt WebEngine is not installed, so this is Qt's basic HTML rendering: "
        "the content is complete but charts, collapsing sections and search need "
        "a real browser. Use <i>Open in browser</i> for the full report.");
    note->setWordWrap(true);
    note->setStyleSheet("color: palette(mid)");
    layout->addWidget(note);
#endif
    layout->addWidget(m_view, 1);

    connect(open, &QPushButton::clicked, this, &ReportPage::openExternally);
    connect(save, &QPushButton::clicked, this, &ReportPage::saveCopy);
}

void ReportPage::load(const QString &htmlPath)
{
    m_path = htmlPath;
    m_status->setText(htmlPath);
#ifdef HAVE_WEBENGINE
    m_view->load(QUrl::fromLocalFile(htmlPath));
#else
    QFile f(htmlPath);
    if (f.open(QIODevice::ReadOnly))
        m_view->setHtml(QString::fromUtf8(f.readAll()));
#endif
}

void ReportPage::openExternally()
{
    if (!m_path.isEmpty())
        QDesktopServices::openUrl(QUrl::fromLocalFile(m_path));
}

void ReportPage::saveCopy()
{
    if (m_path.isEmpty())
        return;
    const QString to = QFileDialog::getSaveFileName(this, "Save report as", m_path, "HTML (*.html)");
    if (to.isEmpty() || to == m_path)
        return;
    QFile::remove(to);
    if (QFile::copy(m_path, to))
        m_status->setText(QStringLiteral("Saved a copy to %1").arg(to));
}
