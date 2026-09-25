#include "FindingsModel.h"

#include <QJsonObject>
#include <QSet>

FindingsModel::FindingsModel(QObject *parent) : QAbstractTableModel(parent) {}

void FindingsModel::setFindings(const QJsonArray &findings)
{
    beginResetModel();
    m_rows.clear();
    for (const QJsonValue &v : findings) {
        const QJsonObject o = v.toObject();
        Finding f;
        f.rsid = o.value("rsid").toString();
        f.gene = o.value("gene").toString();
        f.category = o.value("category").toString();
        f.genotype = o.value("genotype").toString();
        f.status = o.value("status").toString();
        f.description = o.value("description").toString();
        f.note = o.value("note").toString();
        f.magnitude = o.value("magnitude").toInt();
        m_rows << f;
    }
    endResetModel();
}

QStringList FindingsModel::categories() const
{
    QSet<QString> set;
    for (const Finding &f : m_rows)
        set.insert(f.category);
    QStringList out(set.begin(), set.end());
    out.sort();
    return out;
}

int FindingsModel::rowCount(const QModelIndex &parent) const
{
    return parent.isValid() ? 0 : m_rows.size();
}

int FindingsModel::columnCount(const QModelIndex &parent) const
{
    return parent.isValid() ? 0 : ColumnCount;
}

QVariant FindingsModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid() || index.row() >= m_rows.size())
        return {};
    const Finding &f = m_rows[index.row()];
    if (role == Qt::DisplayRole || role == Qt::EditRole) {
        switch (index.column()) {
        case Rsid: return f.rsid;
        case Gene: return f.gene;
        case Category: return f.category;
        case Genotype: return f.genotype;
        case Status: return f.status;
        case Magnitude: return f.magnitude;
        case Description: return f.description;
        }
    }
    if (role == Qt::ToolTipRole && !f.note.isEmpty())
        return f.note;
    return {};
}

QVariant FindingsModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (orientation != Qt::Horizontal || role != Qt::DisplayRole)
        return {};
    static const char *names[ColumnCount] = {
        "rsID", "Gene", "Category", "Genotype", "Status", "Magnitude", "Description"};
    return section >= 0 && section < ColumnCount ? QString::fromLatin1(names[section]) : QVariant();
}

QString FindingsModel::toTsv(const QVector<Finding> &rows)
{
    QString out = QStringLiteral("rsid\tgene\tcategory\tgenotype\tstatus\tmagnitude\tdescription\n");
    for (const Finding &f : rows)
        out += QStringLiteral("%1\t%2\t%3\t%4\t%5\t%6\t%7\n")
                   .arg(f.rsid, f.gene, f.category, f.genotype, f.status)
                   .arg(f.magnitude)
                   .arg(QString(f.description).replace('\t', ' ').replace('\n', ' '));
    return out;
}

FindingsFilter::FindingsFilter(QObject *parent) : QSortFilterProxyModel(parent) {}

void FindingsFilter::setMinimumMagnitude(int m) { m_minMagnitude = m; invalidateFilter(); }
void FindingsFilter::setCategory(const QString &category) { m_category = category; invalidateFilter(); }
void FindingsFilter::setSearchText(const QString &text) { m_search = text.trimmed(); invalidateFilter(); }

bool FindingsFilter::filterAcceptsRow(int sourceRow, const QModelIndex &sourceParent) const
{
    const auto *model = static_cast<const FindingsModel *>(sourceModel());
    if (!model || sourceRow >= model->findings().size())
        return false;
    const Finding &f = model->findings()[sourceRow];
    if (f.magnitude < m_minMagnitude)
        return false;
    if (!m_category.isEmpty() && f.category != m_category)
        return false;
    if (!m_search.isEmpty()) {
        const Qt::CaseSensitivity cs = Qt::CaseInsensitive;
        if (!f.rsid.contains(m_search, cs) && !f.gene.contains(m_search, cs) &&
            !f.status.contains(m_search, cs) && !f.description.contains(m_search, cs) &&
            !f.category.contains(m_search, cs))
            return false;
    }
    Q_UNUSED(sourceParent);
    return true;
}

bool FindingsFilter::lessThan(const QModelIndex &left, const QModelIndex &right) const
{
    if (left.column() == FindingsModel::Magnitude)
        return left.data().toInt() < right.data().toInt();
    return QSortFilterProxyModel::lessThan(left, right);
}
