#pragma once

#include <QWidget>

#include "AnalysisRunner.h"
#include "JobRunner.h"

class QLineEdit;
class QLabel;
class QPushButton;
class QComboBox;

// Choose the genome file and the run parameters. The file's layout is
// detected from its contents as soon as it is chosen; a FASTQ is routed to
// the WGS page instead of the array analysis.
class InputPage : public QWidget {
    Q_OBJECT
public:
    explicit InputPage(QWidget *parent = nullptr);

    void browse();
    void setGenomePath(const QString &path);

signals:
    void runRequested(const AnalysisRunner::Request &request);
    void wgsRunRequested(const QString &fastqPath, const QString &subjectName);

private:
    void onPathChanged();
    void onRun();
    void restoreSettings();
    void saveSettings();

    void buildTestGenome();

    QLineEdit *m_genome, *m_name, *m_dataDir, *m_outDir;
    QLabel *m_format, *m_testStatus;
    QPushButton *m_run, *m_testGenome, *m_testStop;
    QComboBox *m_testSample;
    JobRunner m_job;
    bool m_isFastq = false;
};
