#pragma once

#include <QWidget>

#include "JobRunner.h"

class QLabel;
class QPlainTextEdit;
class QProgressBar;
class QPushButton;
class QTableWidget;

// The annotation sources the analysis reads: what is present, which release
// it is, and buttons to update them. The updaters are native (core/), run on
// a worker thread with progress here.
class DataSourcesPage : public QWidget {
    Q_OBJECT
public:
    explicit DataSourcesPage(QWidget *parent = nullptr);

    void refresh();

private:
    void runJob(const QString &title, const JobRunner::Job &job);
    void setBusy(bool busy);

    QTableWidget *m_table;
    QPlainTextEdit *m_log;
    QProgressBar *m_progress;
    QLabel *m_stage;
    QPushButton *m_clinvar, *m_pharmgkb, *m_validate, *m_rsids, *m_reference, *m_clearCaches, *m_cancel;
    JobRunner m_job;
};
