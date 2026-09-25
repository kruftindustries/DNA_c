// A gzip file as a sequential QIODevice, decompressed in-process with the
// vendored miniz (third_party/miniz), so reading ClinVar's variant_summary
// or Ensembl's chain file needs no gzip binary on any platform. Handles
// multi-member files (bgzip output is a run of gzip members).
#pragma once

#include <QByteArray>
#include <QFile>
#include <QIODevice>

#include <functional>

#include "third_party/miniz/miniz.h"

class GzipStream : public QIODevice {
    Q_OBJECT
public:
    explicit GzipStream(const QString &path, QObject *parent = nullptr);
    ~GzipStream() override;

    bool open(OpenMode mode) override;
    void close() override;
    bool atEnd() const override;
    bool isSequential() const override { return true; }
    bool waitForReadyRead(int msecs) override;
    QString errorText() const { return m_error; }

    // Decompress `gz` to `dest` in 4 MiB chunks, through a QSaveFile so a
    // failure or a cancellation leaves no partial output behind. `progress`
    // is called with the bytes written so far and returns false to cancel
    // (then `error` is "cancelled"). Replaces every `gzip -dc` shell-out:
    // there is no gzip on Windows, and a process that never started reports
    // exit code 0, which once handed samtools an empty FASTA.
    static bool decompressFile(const QString &gz, const QString &dest,
                               const std::function<bool(qint64)> &progress, QString *error);

protected:
    qint64 readData(char *data, qint64 maxSize) override;
    qint64 writeData(const char *, qint64) override { return -1; }

private:
    bool fill();                 // more compressed input; false at end of file
    bool takeByte(quint8 *b);
    bool startMember();          // parse a gzip header, set up raw inflate
    void endMember();

    QFile m_file;
    QByteArray m_in;
    int m_inPos = 0;
    mz_stream m_strm{};
    bool m_inflating = false;
    bool m_finished = false;     // no more members
    int m_members = 0;           // complete members decoded so far
    QString m_error;
};
