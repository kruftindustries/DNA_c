#pragma once

#include <QWidget>

class QLabel;
#ifdef HAVE_WEBENGINE
class QWebEngineView;
#else
class QTextBrowser;
#endif

// Shows the generated report. With Qt WebEngine the interactive report
// renders in-app; without it Qt's basic HTML engine shows the content and
// the full report opens in the system browser.
class ReportPage : public QWidget {
    Q_OBJECT
public:
    explicit ReportPage(QWidget *parent = nullptr);

    void load(const QString &htmlPath);
    QString path() const { return m_path; }

private:
    void openExternally();
    void saveCopy();

    QString m_path;
    QLabel *m_status;
#ifdef HAVE_WEBENGINE
    QWebEngineView *m_view;
#else
    QTextBrowser *m_view;
#endif
};
