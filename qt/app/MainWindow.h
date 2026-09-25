#pragma once

#include <QMainWindow>

#include "AnalysisRunner.h"

class QListWidget;
class QStackedWidget;
class InputPage;
class RunPage;
class ReportPage;
class FindingsPage;
class DataSourcesPage;
class WgsPage;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = nullptr);

private:
    enum Page { Input, Run, Report, Findings, DataSources, Wgs };
    void showPage(Page page);
    void startAnalysis(const AnalysisRunner::Request &req);
    void openSettings();
    void about();

    QListWidget *m_nav;
    QStackedWidget *m_stack;
    InputPage *m_input;
    RunPage *m_run;
    ReportPage *m_report;
    FindingsPage *m_findings;
    DataSourcesPage *m_data;
    WgsPage *m_wgs;
    AnalysisRunner m_runner;
};
