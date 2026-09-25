// UCSC chain-file liftover: a (chromosome, position) on one assembly to the
// other, with strand. Ensembl publishes GRCh37_to_GRCh38.chain.gz (160 KB)
// beside each genome, which is all the reference-as-sample builder needs to
// find a GRCh37 site in GRCh38 without asking the REST API. The file lists
// every chain twice, once with "chr" names and once without; only the bare
// names are kept. Chains on patch scaffolds are ignored.
#pragma once

#include <QHash>
#include <QIODevice>
#include <QString>
#include <QVector>

class ChainLift {
public:
    struct Lifted {
        QString chrom;
        qint64 pos = 0;        // 1-based
        bool reverse = false;  // the target strand is the reverse of the source
    };

    // Parse chain text from an open device (plain or GzipStream).
    bool parse(QIODevice &in, QString *error);
    bool load(const QString &gzPath, QString *error);

    // False when no chain block covers the position.
    bool lift(const QString &chrom, qint64 pos, Lifted *out) const;

    int blockCount() const;

private:
    struct Block {
        qint64 tStart, tEnd;   // 0-based, half-open, on the source
        qint64 qStart;         // 0-based on the target, in strand coordinates
        qint64 qSize;
        int target;            // index into m_targets
        bool reverse;
    };
    QVector<QString> m_targets;
    QHash<QString, QVector<Block>> m_blocks;   // source chromosome -> blocks sorted by tStart
};
