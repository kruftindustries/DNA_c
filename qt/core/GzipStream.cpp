#include "GzipStream.h"

#include <QSaveFile>

GzipStream::GzipStream(const QString &path, QObject *parent) : QIODevice(parent), m_file(path) {}

GzipStream::~GzipStream()
{
    close();
}

bool GzipStream::decompressFile(const QString &gz, const QString &dest,
                                const std::function<bool(qint64)> &progress, QString *error)
{
    GzipStream in(gz);
    if (!in.open(QIODevice::ReadOnly)) {
        if (error) *error = in.errorText();
        return false;
    }
    QSaveFile out(dest);
    if (!out.open(QIODevice::WriteOnly)) {
        if (error) *error = QStringLiteral("cannot write %1: %2").arg(dest, out.errorString());
        return false;
    }
    qint64 written = 0;
    QByteArray chunk;
    while (!(chunk = in.read(4 << 20)).isEmpty()) {
        if (out.write(chunk) != chunk.size()) {
            if (error) *error = QStringLiteral("cannot write %1: %2").arg(dest, out.errorString());
            return false;
        }
        written += chunk.size();
        if (progress && !progress(written)) {
            if (error) *error = QStringLiteral("cancelled");
            return false;   // QSaveFile discards the partial output
        }
    }
    if (!in.errorText().isEmpty()) {
        if (error) *error = in.errorText();
        return false;
    }
    if (!out.commit()) {
        if (error) *error = QStringLiteral("cannot finish %1: %2").arg(dest, out.errorString());
        return false;
    }
    return true;
}

bool GzipStream::open(OpenMode mode)
{
    if (!(mode & ReadOnly) || (mode & WriteOnly))
        return false;
    if (!m_file.open(QIODevice::ReadOnly)) {
        m_error = m_file.errorString();
        return false;
    }
    m_in.clear();
    m_inPos = 0;
    m_finished = false;
    m_members = 0;
    m_error.clear();
    return QIODevice::open(mode | Unbuffered);
}

void GzipStream::close()
{
    endMember();
    m_file.close();
    QIODevice::close();
}

bool GzipStream::atEnd() const
{
    return m_finished && QIODevice::atEnd();
}

bool GzipStream::waitForReadyRead(int)
{
    return !m_finished;
}

bool GzipStream::fill()
{
    if (m_inPos < m_in.size())
        return true;
    m_in = m_file.read(1 << 16);
    m_inPos = 0;
    return !m_in.isEmpty();
}

bool GzipStream::takeByte(quint8 *b)
{
    if (!fill())
        return false;
    *b = quint8(m_in.at(m_inPos++));
    return true;
}

// RFC 1952: ID1 ID2 CM FLG MTIME(4) XFL OS, then the optional fields FLG
// announces, then a raw deflate stream, then CRC32 and ISIZE.
bool GzipStream::startMember()
{
    quint8 h[10];
    for (int i = 0; i < 10; i++) {
        if (!takeByte(&h[i])) {
            if (i == 0) {
                m_finished = true;   // clean end of file between members
                return false;
            }
            m_error = QStringLiteral("truncated gzip header");
            return false;
        }
    }
    if (h[0] != 0x1f || h[1] != 0x8b || h[2] != 8) {
        // After a complete member, anything that is not another gzip header
        // is trailing data, not an error: samtools' RAZF files (the 1000
        // Genomes reference is one) append a block index after the stream,
        // and gzip itself reads them with "trailing garbage ignored".
        if (m_members > 0) {
            m_finished = true;
            return false;
        }
        m_error = QStringLiteral("not a gzip stream");
        return false;
    }
    const quint8 flags = h[3];
    quint8 b;
    if (flags & 0x04) {   // FEXTRA
        quint8 lo, hi;
        if (!takeByte(&lo) || !takeByte(&hi)) { m_error = "truncated gzip header"; return false; }
        for (int n = lo | (hi << 8); n > 0; n--)
            if (!takeByte(&b)) { m_error = "truncated gzip header"; return false; }
    }
    if (flags & 0x08)     // FNAME
        do { if (!takeByte(&b)) { m_error = "truncated gzip header"; return false; } } while (b);
    if (flags & 0x10)     // FCOMMENT
        do { if (!takeByte(&b)) { m_error = "truncated gzip header"; return false; } } while (b);
    if (flags & 0x02)     // FHCRC
        if (!takeByte(&b) || !takeByte(&b)) { m_error = "truncated gzip header"; return false; }

    memset(&m_strm, 0, sizeof(m_strm));
    if (mz_inflateInit2(&m_strm, -MZ_DEFAULT_WINDOW_BITS) != MZ_OK) {
        m_error = QStringLiteral("inflate init failed");
        return false;
    }
    m_inflating = true;
    m_members++;
    return true;
}

void GzipStream::endMember()
{
    if (m_inflating) {
        mz_inflateEnd(&m_strm);
        m_inflating = false;
    }
}

qint64 GzipStream::readData(char *data, qint64 maxSize)
{
    if (maxSize <= 0)
        return 0;
    for (;;) {
        if (m_finished)
            return 0;
        if (!m_inflating) {
            if (!startMember())
                return m_finished ? 0 : -1;
        }
        if (!fill()) {
            m_error = QStringLiteral("truncated gzip data");
            return -1;
        }
        m_strm.next_in = reinterpret_cast<const unsigned char *>(m_in.constData() + m_inPos);
        m_strm.avail_in = static_cast<unsigned>(m_in.size() - m_inPos);
        m_strm.next_out = reinterpret_cast<unsigned char *>(data);
        m_strm.avail_out = static_cast<unsigned>(qMin<qint64>(maxSize, 1u << 30));
        const int status = mz_inflate(&m_strm, MZ_NO_FLUSH);
        const qint64 produced = qint64(maxSize > (1 << 30) ? (1u << 30) : maxSize) - m_strm.avail_out;
        m_inPos = m_in.size() - int(m_strm.avail_in);
        if (status == MZ_STREAM_END) {
            endMember();
            // CRC32 and ISIZE follow the deflate stream.
            quint8 b;
            for (int i = 0; i < 8; i++)
                if (!takeByte(&b))
                    break;
            if (produced > 0)
                return produced;
            continue;   // next member, if any
        }
        if (status != MZ_OK && status != MZ_BUF_ERROR) {
            m_error = QStringLiteral("corrupt gzip data (%1)").arg(mz_error(status));
            return -1;
        }
        if (produced > 0)
            return produced;
        // No output yet: needs more input (the loop refills).
    }
}
