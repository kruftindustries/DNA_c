// The lifestyle findings from the analysis JSON as a sortable, filterable
// table.
#pragma once

#include <QAbstractTableModel>
#include <QJsonArray>
#include <QSortFilterProxyModel>
#include <QVector>

struct Finding {
    QString rsid, gene, category, genotype, status, description, note;
    int magnitude = 0;
};

class FindingsModel : public QAbstractTableModel {
    Q_OBJECT
public:
    enum Column { Rsid, Gene, Category, Genotype, Status, Magnitude, Description, ColumnCount };

    explicit FindingsModel(QObject *parent = nullptr);

    void setFindings(const QJsonArray &findings);
    const QVector<Finding> &findings() const { return m_rows; }
    QStringList categories() const;   // distinct, sorted

    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const override;

    // Tab-separated export of the rows given, header first.
    static QString toTsv(const QVector<Finding> &rows);

private:
    QVector<Finding> m_rows;
};

class FindingsFilter : public QSortFilterProxyModel {
    Q_OBJECT
public:
    explicit FindingsFilter(QObject *parent = nullptr);

    void setMinimumMagnitude(int m);
    void setCategory(const QString &category);   // empty = all
    void setSearchText(const QString &text);

protected:
    bool filterAcceptsRow(int sourceRow, const QModelIndex &sourceParent) const override;
    bool lessThan(const QModelIndex &left, const QModelIndex &right) const override;

private:
    int m_minMagnitude = 0;
    QString m_category, m_search;
};
