#include "ChainLift.h"

#include <QSet>
#include <QTextStream>

#include <algorithm>

#include "GzipStream.h"

namespace {

const QSet<QString> &primaryChromosomes()
{
    static const QSet<QString> s = [] {
        QSet<QString> out;
        for (int i = 1; i <= 22; ++i)
            out.insert(QString::number(i));
        out.insert("X"); out.insert("Y"); out.insert("MT");
        return out;
    }();
    return s;
}

} // namespace

bool ChainLift::parse(QIODevice &in, QString *error)
{
    // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    // size dt dq          (one per aligned block; the last line has only size)
    bool keep = false;
    qint64 t = 0, q = 0, qSize = 0;
    int target = -1;
    bool reverse = false;
    QString source;
    int lineNo = 0;
    while (!in.atEnd()) {
        const QString line = QString::fromLatin1(in.readLine()).trimmed();
        lineNo++;
        if (line.isEmpty())
            continue;
        const QStringList f = line.split(' ', Qt::SkipEmptyParts);
        if (f[0] == "chain") {
            if (f.size() < 12) {
                if (error) *error = QStringLiteral("malformed chain header at line %1").arg(lineNo);
                return false;
            }
            source = f[2];
            const QString qName = f[7];
            keep = primaryChromosomes().contains(source) && primaryChromosomes().contains(qName)
                   && f[4] == "+";
            if (!keep)
                continue;
            t = f[5].toLongLong();
            qSize = f[8].toLongLong();
            reverse = f[9] == "-";
            q = f[10].toLongLong();
            target = m_targets.indexOf(qName);
            if (target < 0) {
                m_targets.append(qName);
                target = m_targets.size() - 1;
            }
            continue;
        }
        if (!keep)
            continue;
        const qint64 size = f[0].toLongLong();
        const qint64 dt = f.size() > 1 ? f[1].toLongLong() : 0;
        const qint64 dq = f.size() > 2 ? f[2].toLongLong() : 0;
        m_blocks[source].append({t, t + size, q, qSize, target, reverse});
        t += size + dt;
        q += size + dq;
    }
    for (auto &blocks : m_blocks)
        std::sort(blocks.begin(), blocks.end(), [](const Block &a, const Block &b) { return a.tStart < b.tStart; });
    if (m_blocks.isEmpty()) {
        if (error) *error = QStringLiteral("no chains on primary chromosomes");
        return false;
    }
    return true;
}

bool ChainLift::load(const QString &gzPath, QString *error)
{
    GzipStream in(gzPath);
    if (!in.open(QIODevice::ReadOnly)) {
        if (error) *error = QStringLiteral("cannot read %1: %2").arg(gzPath, in.errorText());
        return false;
    }
    return parse(in, error);
}

bool ChainLift::lift(const QString &chrom, qint64 pos, Lifted *out) const
{
    const auto it = m_blocks.find(chrom);
    if (it == m_blocks.end())
        return false;
    const qint64 z = pos - 1;
    const QVector<Block> &blocks = it.value();
    // Last block starting at or before z.
    auto b = std::upper_bound(blocks.begin(), blocks.end(), z,
                              [](qint64 v, const Block &blk) { return v < blk.tStart; });
    if (b == blocks.begin())
        return false;
    --b;
    if (z >= b->tEnd)
        return false;
    qint64 qz = b->qStart + (z - b->tStart);
    if (b->reverse)
        qz = b->qSize - 1 - qz;
    out->chrom = m_targets[b->target];
    out->pos = qz + 1;
    out->reverse = b->reverse;
    return true;
}

int ChainLift::blockCount() const
{
    int n = 0;
    for (const auto &blocks : m_blocks)
        n += blocks.size();
    return n;
}
