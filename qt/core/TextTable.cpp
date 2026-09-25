#include "TextTable.h"

TsvReader::TsvReader(QIODevice *device) : m_device(device) {}

int TsvReader::nextByte()
{
    if (m_pos >= m_buffer.size()) {
        if (m_eof)
            return -1;
        m_buffer = m_device->read(1 << 16);
        m_pos = 0;
        if (m_buffer.isEmpty()) {
            // A blocking device (a process) may have nothing yet.
            if (m_device->atEnd() || !m_device->waitForReadyRead(-1)) {
                m_buffer = m_device->readAll();
                if (m_buffer.isEmpty()) {
                    m_eof = true;
                    return -1;
                }
            } else {
                m_buffer = m_device->read(1 << 16);
                if (m_buffer.isEmpty()) {
                    m_eof = true;
                    return -1;
                }
            }
        }
        m_bytesRead += m_buffer.size();
    }
    return static_cast<unsigned char>(m_buffer[m_pos++]);
}

bool TsvReader::readRecord(QStringList &fields)
{
    // Python's csv state machine, non-strict, doublequote=True, no escapechar:
    // a quote only opens a quoted field at the very start of a field.
    for (;;) {
        fields.clear();
        QByteArray field;
        bool sawAnything = false, quoted = false, inQuotes = false;
        for (;;) {
            int c = nextByte();
            if (c < 0) {
                if (!sawAnything)
                    return false;
                fields << QString::fromUtf8(field);
                return true;
            }
            sawAnything = true;
            if (inQuotes) {
                if (c == '"') {
                    int d = nextByte();
                    if (d == '"') {
                        field += '"';
                    } else {
                        inQuotes = false;
                        if (d < 0) {
                            fields << QString::fromUtf8(field);
                            return true;
                        }
                        // Re-process d as an unquoted character.
                        m_pos--;
                    }
                } else if (c == '\r') {
                    // Universal newlines inside a quoted field become "\n".
                    int d = nextByte();
                    if (d != '\n' && d >= 0)
                        m_pos--;
                    field += '\n';
                } else {
                    field += char(c);
                }
                continue;
            }
            if (c == '\t') {
                fields << QString::fromUtf8(field);
                field.clear();
                quoted = false;
                continue;
            }
            if (c == '\n' || c == '\r') {
                if (c == '\r') {
                    int d = nextByte();
                    if (d != '\n' && d >= 0)
                        m_pos--;
                }
                fields << QString::fromUtf8(field);
                break;
            }
            if (c == '"' && field.isEmpty() && !quoted) {
                inQuotes = true;
                quoted = true;
                continue;
            }
            field += char(c);
        }
        // csv.reader yields [] for a blank line and DictReader skips it.
        if (fields.size() == 1 && fields[0].isEmpty())
            continue;
        return true;
    }
}

TsvTable::TsvTable(QIODevice *device) : m_reader(device)
{
    m_ok = m_reader.readRecord(m_header);
}

bool TsvTable::next()
{
    return m_reader.readRecord(m_row);
}

QString TsvTable::field(int column) const
{
    return column >= 0 && column < m_row.size() ? m_row[column] : QString();
}

QByteArray tsvRow(const QStringList &fields)
{
    QByteArray out;
    for (int i = 0; i < fields.size(); ++i) {
        if (i)
            out += '\t';
        const QByteArray f = fields[i].toUtf8();
        const bool needsQuote = f.contains('\t') || f.contains('"') ||
                                f.contains('\r') || f.contains('\n');
        if (!needsQuote) {
            out += f;
            continue;
        }
        out += '"';
        for (char c : f) {
            if (c == '"')
                out += '"';
            out += c;
        }
        out += '"';
    }
    out += "\r\n";
    return out;
}

LineReader::LineReader(QIODevice *device) : m_device(device) {}

bool LineReader::readLine(QByteArray &line)
{
    line.clear();
    for (;;) {
        int nl = m_buffer.indexOf('\n', m_pos);
        if (nl >= 0) {
            line += m_buffer.mid(m_pos, nl - m_pos);
            m_pos = nl + 1;
            if (line.endsWith('\r'))
                line.chop(1);
            return true;
        }
        line += m_buffer.mid(m_pos);
        m_buffer.clear();
        m_pos = 0;
        if (m_eof)
            return !line.isEmpty();
        QByteArray chunk = m_device->read(1 << 16);
        if (chunk.isEmpty()) {
            if (!m_device->atEnd() && m_device->waitForReadyRead(-1))
                chunk = m_device->read(1 << 16);
            if (chunk.isEmpty())
                chunk = m_device->readAll();
        }
        if (chunk.isEmpty()) {
            m_eof = true;
            if (line.endsWith('\r'))
                line.chop(1);
            return !line.isEmpty();
        }
        m_buffer = chunk;
    }
}
