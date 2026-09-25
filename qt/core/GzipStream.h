// A gzip file as a sequential QIODevice, decompressed in-process with the
// vendored miniz (third_party/miniz), so reading ClinVar's variant_summary
// or Ensembl's chain file needs no gzip binary on any platform. Handles
// multi-member files (bgzip output is a run of gzip members).
#pragma once

#include <QByteArray>
#include <QFile>
#include <QIODevice>

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
