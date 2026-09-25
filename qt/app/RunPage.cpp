#include "RunPage.h"

#include <QFontDatabase>
#include <QHBoxLayout>
#include <QLabel>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QVBoxLayout>

RunPage::RunPage(QWidget *parent) : QWidget(parent)
{
    auto *layout = new QVBoxLayout(this);
    auto *top = new QHBoxLayout;
    m_stage = new QLabel("Idle");
    m_progress = new QProgressBar;
    m_progress->setRange(0, 100);
    m_cancel = new QPushButton("Cancel");
    m_cancel->setEnabled(false);
    top->addWidget(m_stage);
    top->addWidget(m_progress, 1);
    top->addWidget(m_cancel);
    layout->addLayout(top);

    m_log = new QPlainTextEdit;
    m_log->setReadOnly(true);
    m_log->setMaximumBlockCount(20000);
    m_log->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    layout->addWidget(m_log, 1);

    connect(m_cancel, &QPushButton::clicked, this, &RunPage::cancelRequested);
}

void RunPage::clear() { m_log->clear(); m_progress->setValue(0); m_stage->setText("Idle"); }
void RunPage::appendLine(const QString &line) { m_log->appendPlainText(line); }
void RunPage::setStage(const QString &stage, int percent) { m_stage->setText(stage); m_progress->setValue(percent); }
void RunPage::setRunning(bool running) { m_cancel->setEnabled(running); }
