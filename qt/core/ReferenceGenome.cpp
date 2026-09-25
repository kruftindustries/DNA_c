#include "ReferenceGenome.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QLocale>
#include <QProcess>
#include <QSaveFile>
#include <QUrl>

#include "Downloader.h"
#include "GzipStream.h"
#include "ToolLocator.h"
#include "WgsPipeline.h"

namespace ReferenceGenome {

const char *const kUrl =
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz";

Status status(const QString &root)
{
    const WgsConfig cfg = WgsConfig::defaults(root);
    Status s;
    s.fasta = cfg.refGenome;
    s.fai = cfg.refGenome + ".fai";
    s.mmi = cfg.refMmi;
    s.dir = QFileInfo(s.fasta).path();
    const QFileInfo fa(s.fasta), fai(s.fai), mmi(s.mmi);
    s.hasFasta = fa.exists() && fa.size() > 0;
    s.hasFai = fai.exists() && fai.size() > 0;
    s.hasMmi = mmi.exists() && mmi.size() > 0;
    s.fastaBytes = s.hasFasta ? fa.size() : -1;
    s.mmiBytes = s.hasMmi ? mmi.size() : -1;
    return s;
}

namespace {

// Run a tool, polling so a cancellation kills it.
bool runTool(const QString &program, const QStringList &args, const Reporter &reporter, QString *error)
{
    QProcess p;
    p.setProcessChannelMode(QProcess::MergedChannels);
    p.start(program, args);
    if (!p.waitForStarted(10000)) {
        if (error) *error = QStringLiteral("cannot start %1: %2").arg(program, p.errorString());
        return false;
    }
    while (!p.waitForFinished(500)) {
        if (reporter.cancelled()) {
            p.kill();
            p.waitForFinished(2000);
            if (error) *error = "cancelled";
            return false;
        }
        const QByteArray out = p.readAll();
        if (!out.trimmed().isEmpty())
            reporter.log("    " + QString::fromUtf8(out).trimmed());
    }
    if (p.exitStatus() != QProcess::NormalExit || p.exitCode() != 0) {
        if (error)
            *error = QStringLiteral("%1 failed: %2").arg(QFileInfo(program).fileName(),
                                                        QString::fromUtf8(p.readAll()).trimmed());
        return false;
    }
    return true;
}

} // namespace

bool setup(const QString &root, const Reporter &reporter, QString *error)
{
    const QString samtools = ToolLocator::find("samtools");
    const QString minimap2 = ToolLocator::find("minimap2");
    if (samtools.isEmpty() || minimap2.isEmpty()) {
        if (error)
            *error = QStringLiteral("%1 not installed. %2")
                         .arg(ToolLocator::missing({"samtools", "minimap2"}).join(", "),
                              ToolLocator::installHint());
        return false;
    }
    Status s = status(root);
    QDir().mkpath(s.dir);
    const QString gz = s.fasta + ".gz";

    // 1. Download (skipped when the FASTA already exists).
    if (!s.hasFasta) {
        if (!QFileInfo(gz).exists() || QFileInfo(gz).size() < (800LL << 20)) {
            reporter.stage(QStringLiteral("Downloading the GRCh37 reference (~900 MB)"), 0);
            reporter.log(QStringLiteral(">>> Downloading %1").arg(QLatin1String(kUrl)));
            if (!Downloader::download(QUrl(QLatin1String(kUrl)), gz, reporter, error))
                return false;
        } else {
            reporter.log(QStringLiteral("    %1 already downloaded").arg(gz));
        }

        // 2. Decompress in-process to a plain FASTA (~3 GB), which is what
        //    samtools faidx and bcftools mpileup read.
        reporter.stage(QStringLiteral("Decompressing the reference (~3 GB)"), 40);
        const qint64 total = QFileInfo(gz).size();
        const auto progress = [&](qint64 written) {
            // The plain FASTA is ~3.4x the gzip; report against that.
            reporter.stage(QStringLiteral("Decompressing the reference (%1)")
                               .arg(QLocale().formattedDataSize(written)),
                           40 + int(qMin<qint64>(25, written * 25 / qMax<qint64>(1, total * 34 / 10))));
            return !reporter.cancelled();
        };
        if (!GzipStream::decompressFile(gz, s.fasta, progress, error)) {
            return false;
        }
        QFile::remove(gz);
        reporter.log(QStringLiteral("    Wrote %1").arg(s.fasta));
        s = status(root);
    }

    // 3. samtools index (fast).
    if (!s.hasFai) {
        reporter.stage(QStringLiteral("Indexing with samtools faidx"), 66);
        if (!runTool(samtools, {"faidx", s.fasta}, reporter, error))
            return false;
    }

    // 4. minimap2 index (~5 minutes, ~3.5 GB).
    if (!s.hasMmi) {
        reporter.stage(QStringLiteral("Building the minimap2 index (~5 min)"), 70);
        const QString tmp = s.mmi + ".part";
        if (!runTool(minimap2, {"-x", "sr", "-d", tmp, s.fasta}, reporter, error)) {
            QFile::remove(tmp);
            return false;
        }
        QFile::remove(s.mmi);
        if (!QFile::rename(tmp, s.mmi)) {
            if (error) *error = QStringLiteral("cannot move %1 into place").arg(tmp);
            return false;
        }
    }

    s = status(root);
    reporter.log(QStringLiteral("    Reference ready: %1 (%2), index %3")
                     .arg(s.fasta, QLocale().formattedDataSize(s.fastaBytes),
                          QLocale().formattedDataSize(s.mmiBytes)));
    reporter.stage(QStringLiteral("Reference ready"), 100);
    return true;
}

} // namespace ReferenceGenome
