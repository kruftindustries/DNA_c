#include "JobRunner.h"

#include <QMetaObject>

JobRunner::JobRunner(QObject *parent) : QObject(parent) {}

JobRunner::~JobRunner()
{
    if (m_thread) {
        m_cancel = true;
        m_thread->wait();
    }
}

void JobRunner::start(const Job &job)
{
    if (isRunning())
        return;
    m_cancel = false;
    if (m_thread) {
        m_thread->deleteLater();
        m_thread = nullptr;
    }
    m_thread = QThread::create([this, job] {
        Reporter rep;
        rep.log = [this](const QString &s) {
            QMetaObject::invokeMethod(this, [this, s] { emit logLine(s); }, Qt::QueuedConnection);
        };
        rep.stage = [this](const QString &name, int pct) {
            QMetaObject::invokeMethod(this, [this, name, pct] { emit stage(name, pct); }, Qt::QueuedConnection);
        };
        rep.cancelled = [this] { return m_cancel.load(); };
        QString error;
        const bool ok = job(rep, &error);
        QMetaObject::invokeMethod(this, [this, ok, error] { emit finished(ok, error); }, Qt::QueuedConnection);
    });
    m_thread->start();
}
