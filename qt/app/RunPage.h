#pragma once

#include <QWidget>

class QLabel;
class QPlainTextEdit;
class QProgressBar;
class QPushButton;

// Live log and progress for whatever is currently running.
class RunPage : public QWidget {
    Q_OBJECT
public:
    explicit RunPage(QWidget *parent = nullptr);

    void clear();
    void appendLine(const QString &line);
    void setStage(const QString &stage, int percent);
    void setRunning(bool running);

signals:
    void cancelRequested();

private:
    QLabel *m_stage;
    QProgressBar *m_progress;
    QPlainTextEdit *m_log;
    QPushButton *m_cancel;
};
