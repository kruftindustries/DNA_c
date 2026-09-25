// Runs one core job on a worker thread and forwards its Reporter callbacks
// to the GUI thread as signals. Cancellation is a flag the job polls.
#pragma once

#include <QObject>
#include <QThread>
#include <atomic>
#include <functional>

#include "Reporter.h"

class JobRunner : public QObject {
    Q_OBJECT
public:
    using Job = std::function<bool(const Reporter &, QString *error)>;

    explicit JobRunner(QObject *parent = nullptr);
    ~JobRunner() override;

    void start(const Job &job);
    void cancel() { m_cancel = true; }
    bool isRunning() const { return m_thread && m_thread->isRunning(); }

signals:
    void logLine(const QString &line);
    void stage(const QString &name, int percent);
    void finished(bool ok, const QString &error);

private:
    QThread *m_thread = nullptr;
    std::atomic<bool> m_cancel{false};
};
