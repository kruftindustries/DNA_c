#include "TargetRegions.h"

#include <QFile>
#include <QJsonObject>

#include <algorithm>
#include <cstring>
#include <vector>

#include "EnsemblLookup.h"
#include "TextTable.h"

namespace TargetRegions {

namespace {

struct Site {
    char chrom[8];
    qint64 pos;
};

bool allDigits(const QString &s)
{
    if (s.isEmpty())
        return false;
    for (QChar c : s)
        if (!c.isDigit())
            return false;
    return true;
}

// Python sorts by (chrom.zfill(2), pos).
QByteArray zfill2(const char *chrom)
{
    QByteArray c(chrom);
    while (c.size() < 2)
        c.prepend('0');
    return c;
}

} // namespace

qint64 build(const QString &clinvarTsv, const QString &lookupJson, const QString &bedPath,
             const Reporter &reporter, QString *error)
{
    std::vector<Site> sites;
    sites.reserve(5000000);
    auto add = [&](const QString &chrom, const QString &pos) {
        if (chrom.isEmpty() || !allDigits(pos))
            return;
        const QByteArray c = chrom.toUtf8();
        if (c.size() >= int(sizeof(Site::chrom)))
            return;
        Site s;
        std::memset(s.chrom, 0, sizeof s.chrom);
        std::memcpy(s.chrom, c.constData(), size_t(c.size()));
        s.pos = pos.toLongLong();
        sites.push_back(s);
    };

    QFile cv(clinvarTsv);
    if (cv.open(QIODevice::ReadOnly)) {
        TsvTable table(&cv);
        const int cChrom = table.column("chrom"), cPos = table.column("pos");
        qint64 rows = 0;
        while (table.next()) {
            add(table.field(cChrom), table.field(cPos));
            if ((++rows & 0xFFFFF) == 0 && reporter.cancelled()) {
                if (error)
                    *error = QStringLiteral("cancelled");
                return -1;
            }
        }
    }

    const QJsonObject lookup = EnsemblLookup::loadLookup(lookupJson);
    for (auto it = lookup.begin(); it != lookup.end(); ++it) {
        const QJsonObject info = it.value().toObject();
        add(info.value("chrom").toString(), info.value("pos").toString());
    }

    if (sites.empty())
        return 0;

    std::sort(sites.begin(), sites.end(), [](const Site &a, const Site &b) {
        const int c = std::strcmp(zfill2(a.chrom).constData(), zfill2(b.chrom).constData());
        return c != 0 ? c < 0 : a.pos < b.pos;
    });
    sites.erase(std::unique(sites.begin(), sites.end(), [](const Site &a, const Site &b) {
        return a.pos == b.pos && std::strcmp(a.chrom, b.chrom) == 0;
    }), sites.end());

    QFile bed(bedPath);
    if (!bed.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (error)
            *error = QStringLiteral("cannot write %1").arg(bedPath);
        return -1;
    }
    QByteArray out;
    out.reserve(1 << 20);
    for (const Site &s : sites) {
        out += s.chrom;
        out += '\t';
        out += QByteArray::number(s.pos - 1);
        out += '\t';
        out += QByteArray::number(s.pos);
        out += '\n';
        if (out.size() > (1 << 20)) {
            bed.write(out);
            out.clear();
        }
    }
    bed.write(out);
    reporter.log(QStringLiteral("  Built target regions: %L1 positions").arg(qint64(sites.size())));
    return qint64(sites.size());
}

} // namespace TargetRegions
