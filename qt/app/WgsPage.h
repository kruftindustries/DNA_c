#pragma once

#include <QWidget>

#include "JobRunner.h"

class QCheckBox;
class QLabel;
class QLineEdit;
class QPlainTextEdit;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QTableWidget;

// The FASTQ pipeline, natively: tool discovery, reference status, the chr22
// toolchain validation and the full run, all through core/WgsPipeline on a
// worker thread. No Python involved.
class WgsPage : public QWidget {
    Q_OBJECT
public:
    explicit WgsPage(QWidget *parent = nullptr);

    void setFastqPath(const QString &path);
    void refresh();

    // Run the pipeline on a file chosen elsewhere (the Input page's Run on a
    // FASTQ). Offers to fetch the reference first when it is not set up.
    void startForFile(const QString &fastq, const QString &name);

signals:
    void reportReady(const QString &htmlPath, const QString &jsonPath);

private:
    void validateToolchain();
    void runPipeline();
    void setupReference();
    void startPipeline(bool setupReferenceFirst);
    void runJob(const QString &title, const JobRunner::Job &job);
    void setBusy(bool busy);

    QTableWidget *m_tools;
    QLabel *m_requirements, *m_reference, *m_stage;
    QLineEdit *m_fastq, *m_name;
    QSpinBox *m_threads;
    QCheckBox *m_skipQc, *m_full, *m_keep;
    QPushButton *m_validate, *m_run, *m_setupRef, *m_cancel;
    QProgressBar *m_progress;
    QPlainTextEdit *m_log;
    JobRunner m_job;
    bool m_lastWasRun = false;
};
