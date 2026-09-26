// Tab-separated tables with Python csv-module semantics.
//
// The ClinVar extract is written by csv.DictWriter and read by csv.DictReader
// (and by the C port's gh_tsv, which matches them). Both sides have rules
// that a plain split on tabs gets wrong: a field that starts with a quote
// runs until the closing quote and may contain tabs and newlines, a doubled
// quote inside it is a literal quote, and on output a field is quoted only
// when it needs to be. This reproduces those rules so the Qt writer's output
// is byte-identical to the Python's.
#pragma once

#include <QByteArray>
#include <QIODevice>
#include <QStringList>

class TsvReader {
public:
    explicit TsvReader(QIODevice *device);

    // Next record as decoded fields; false at end of input. Blank lines are
    // skipped, as csv.DictReader skips them.
    bool readRecord(QStringList &fields);

    qint64 bytesRead() const { return m_bytesRead; }

private:
    int nextByte();
    QIODevice *m_device;
    QByteArray m_buffer;
    int m_pos = 0;
    qint64 m_bytesRead = 0;
    bool m_eof = false;
};

// A table with a header row: fields by column name.
class TsvTable {
public:
    explicit TsvTable(QIODevice *device);

    bool ok() const { return m_ok; }
    const QStringList &header() const { return m_header; }
    int column(const QString &name) const { return m_header.indexOf(name); }
    bool hasColumn(const QString &name) const { return column(name) >= 0; }

    bool next();                                 // advance to the next row
    QString field(int column) const;             // "" when out of range
    QString field(const QString &name) const { return field(column(name)); }
    const QStringList &row() const { return m_row; }
    qint64 bytesRead() const { return m_reader.bytesRead(); }

private:
    TsvReader m_reader;
    QStringList m_header, m_row;
    bool m_ok = false;
};

// Plain lines from a (possibly blocking) device, terminator stripped.
class LineReader {
public:
    explicit LineReader(QIODevice *device);
    bool readLine(QByteArray &line);

private:
    QIODevice *m_device;
    QByteArray m_buffer;
    int m_pos = 0;
    bool m_eof = false;
};

// csv.writer(delimiter="\t") with QUOTE_MINIMAL and the default "\r\n"
// terminator: a field is quoted when it contains a tab, a quote, a CR or an
// LF, with quotes doubled.
QByteArray tsvRow(const QStringList &fields);

// The non-empty lines of a tool's stdout, without line terminators. On
// Windows a MinGW-built program's stdout is in text mode and every line
// ends in "\r\n"; splitting on '\n' alone left a '\r' on each rsID the
// engine listed, so nothing ever matched the position lookup.
QStringList outputLines(const QString &text);
QList<QByteArray> outputLines(const QByteArray &text);
