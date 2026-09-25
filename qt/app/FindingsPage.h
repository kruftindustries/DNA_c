#pragma once

#include <QJsonObject>
#include <QWidget>

#include "FindingsModel.h"

class QComboBox;
class QLabel;
class QLineEdit;
class QSpinBox;
class QTableView;

// The findings table with search, category and magnitude filters, and
// clipboard/file export of whatever is currently visible.
class FindingsPage : public QWidget {
    Q_OBJECT
public:
    explicit FindingsPage(QWidget *parent = nullptr);

    void setResults(const QJsonObject &results);
    QVector<Finding> visibleFindings() const;

private:
    void updateCount();
    void copyTsv();
    void exportTsv();

    FindingsModel *m_model;
    FindingsFilter *m_filter;
    QTableView *m_table;
    QLineEdit *m_search;
    QComboBox *m_category;
    QSpinBox *m_minMag;
    QLabel *m_count, *m_summary, *m_empty;
    bool m_haveResults = false;
};
