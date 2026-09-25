// Runs the C analysis binary on a genome: once for the JSON results the
// findings table reads, once for the HTML report. Two runs rather than one
// because the binary's --json and --html are separate outputs; the second
// costs about two seconds, dominated by the ClinVar scan.
#pragma once

#include <QJsonObject>
#include <QObject>

#include "ProcessRunner.h"

class AnalysisRunner : public QObject {
    Q_OBJECT
public:
    explicit AnalysisRunner(QObject *parent = nullptr);

    struct Request {
        QString genomePath;
        QString dataDir;
        QString subjectName;
        QString outputDir;
    };

    void start(const Request &req);
    void cancel();
    bool isRunning() const { return m_stage != Stage::Idle; }

    QJsonObject results() const { return m_results; }
    QString reportPath() const { return m_reportPath; }

signals:
    void logLine(const QString &line);
    void stageChanged(const QString &stage, int percent);
    void finished(bool ok, const QString &message);

private:
    enum class Stage { Idle, Json, Html };
    void onFinished(int code, QProcess::ExitStatus status);

    ProcessRunner m_runner;
    Request m_req;
    Stage m_stage = Stage::Idle;
    QJsonObject m_results;
    QString m_reportPath;
};
