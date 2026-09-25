// One external command with live output, used for everything the GUI
// delegates: the analysis binary, the data downloaders, the WGS pipeline.
#pragma once

#include <QObject>
#include <QProcess>
#include <QProcessEnvironment>
#include <QStringList>

class ProcessRunner : public QObject {
    Q_OBJECT
public:
    explicit ProcessRunner(QObject *parent = nullptr);

    void setEnvironment(const QProcessEnvironment &env);
    void setWorkingDirectory(const QString &dir);

    // Starts the command. `capture` keeps stdout for stdoutText() instead of
    // logging it, for commands whose output is data (e.g. --json).
    void start(const QString &program, const QStringList &args, bool capture = false);
    void cancel();

    bool isRunning() const;
    QString stdoutText() const { return m_stdout; }
    QString commandLine() const { return m_commandLine; }

signals:
    void logLine(const QString &line);
    void finished(int exitCode, QProcess::ExitStatus status);
    void failedToStart(const QString &reason);

private:
    void flushLines(QByteArray &buffer, bool final);

    QProcess m_proc;
    QByteArray m_outBuf, m_errBuf;
    QString m_stdout, m_commandLine;
    bool m_capture = false;
};
