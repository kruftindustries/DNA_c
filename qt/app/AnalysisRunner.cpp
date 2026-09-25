#include "AnalysisRunner.h"

#include <QDir>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonParseError>

#include "RepoPaths.h"

AnalysisRunner::AnalysisRunner(QObject *parent) : QObject(parent)
{
    connect(&m_runner, &ProcessRunner::logLine, this, &AnalysisRunner::logLine);
    connect(&m_runner, &ProcessRunner::finished, this, &AnalysisRunner::onFinished);
    connect(&m_runner, &ProcessRunner::failedToStart, this, [this](const QString &why) {
        m_stage = Stage::Idle;
        emit finished(false, QStringLiteral("could not start %1: %2")
                                 .arg(RepoPaths::analysisBinary(), why));
    });
}

void AnalysisRunner::start(const Request &req)
{
    if (isRunning())
        return;
    m_req = req;
    m_results = QJsonObject();
    m_reportPath.clear();

    const QString binary = RepoPaths::analysisBinary();
    if (binary.isEmpty()) {
        emit finished(false, QStringLiteral("analysis binary not found; build it with "
                                            "`make -C c` or set it in Settings"));
        return;
    }
    QDir().mkpath(req.outputDir);

    m_stage = Stage::Json;
    emit stageChanged(QStringLiteral("Analysing"), 5);
    m_runner.start(binary, {"--data", req.dataDir, "--json", "--quiet", req.genomePath},
                   /*capture=*/true);
}

void AnalysisRunner::cancel()
{
    m_runner.cancel();
}

void AnalysisRunner::onFinished(int code, QProcess::ExitStatus status)
{
    if (status != QProcess::NormalExit || code != 0) {
        m_stage = Stage::Idle;
        emit finished(false, QStringLiteral("analysis exited with code %1").arg(code));
        return;
    }

    if (m_stage == Stage::Json) {
        QJsonParseError err;
        const QJsonDocument doc = QJsonDocument::fromJson(m_runner.stdoutText().toUtf8(), &err);
        if (err.error != QJsonParseError::NoError || !doc.isObject()) {
            m_stage = Stage::Idle;
            emit finished(false, QStringLiteral("could not parse results: %1").arg(err.errorString()));
            return;
        }
        m_results = doc.object();
        emit logLine(QStringLiteral("%1 findings, %2 drug-gene interactions")
                         .arg(m_results.value("findings").toArray().size())
                         .arg(m_results.value("drug_findings").toArray().size()));

        m_stage = Stage::Html;
        emit stageChanged(QStringLiteral("Rendering report"), 60);
        m_reportPath = QDir(m_req.outputDir).filePath("GENETIC_HEALTH_REPORT.html");
        QStringList args{"--data", m_req.dataDir, "--html", m_reportPath, "--quiet"};
        if (!m_req.subjectName.isEmpty())
            args << "--name" << m_req.subjectName;
        args << m_req.genomePath;
        m_runner.start(RepoPaths::analysisBinary(), args);
        return;
    }

    m_stage = Stage::Idle;
    emit stageChanged(QStringLiteral("Done"), 100);
    emit finished(true, QStringLiteral("report written to %1").arg(m_reportPath));
}
