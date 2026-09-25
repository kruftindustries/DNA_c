#include "ProcessRunner.h"

ProcessRunner::ProcessRunner(QObject *parent) : QObject(parent)
{
    connect(&m_proc, &QProcess::readyReadStandardOutput, this, [this] {
        const QByteArray chunk = m_proc.readAllStandardOutput();
        if (m_capture) {
            m_stdout += QString::fromUtf8(chunk);
            return;
        }
        m_outBuf += chunk;
        flushLines(m_outBuf, false);
    });
    connect(&m_proc, &QProcess::readyReadStandardError, this, [this] {
        m_errBuf += m_proc.readAllStandardError();
        flushLines(m_errBuf, false);
    });
    connect(&m_proc, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished), this,
            [this](int code, QProcess::ExitStatus status) {
                if (m_capture)
                    m_stdout += QString::fromUtf8(m_proc.readAllStandardOutput());
                else
                    m_outBuf += m_proc.readAllStandardOutput();
                m_errBuf += m_proc.readAllStandardError();
                flushLines(m_outBuf, true);
                flushLines(m_errBuf, true);
                emit finished(code, status);
            });
    connect(&m_proc, &QProcess::errorOccurred, this, [this](QProcess::ProcessError e) {
        if (e == QProcess::FailedToStart)
            emit failedToStart(m_proc.errorString());
    });
}

void ProcessRunner::setEnvironment(const QProcessEnvironment &env) { m_proc.setProcessEnvironment(env); }
void ProcessRunner::setWorkingDirectory(const QString &dir) { m_proc.setWorkingDirectory(dir); }

void ProcessRunner::start(const QString &program, const QStringList &args, bool capture)
{
    m_capture = capture;
    m_stdout.clear();
    m_outBuf.clear();
    m_errBuf.clear();
    m_commandLine = program + ' ' + args.join(' ');
    emit logLine(QStringLiteral("$ %1").arg(m_commandLine));
    m_proc.start(program, args);
}

void ProcessRunner::cancel()
{
    if (m_proc.state() == QProcess::NotRunning)
        return;
    emit logLine(QStringLiteral("(cancelled)"));
    m_proc.kill();
}

bool ProcessRunner::isRunning() const { return m_proc.state() != QProcess::NotRunning; }

void ProcessRunner::flushLines(QByteArray &buffer, bool final)
{
    int at;
    while ((at = buffer.indexOf('\n')) >= 0) {
        emit logLine(QString::fromUtf8(buffer.left(at)).trimmed());
        buffer.remove(0, at + 1);
    }
    if (final && !buffer.isEmpty()) {
        emit logLine(QString::fromUtf8(buffer).trimmed());
        buffer.clear();
    }
}
