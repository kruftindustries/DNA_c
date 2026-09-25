#include "ZipArchive.h"

#include <QFile>

#include "third_party/miniz/miniz.h"

namespace ZipArchive {

namespace {

struct Reader {
    mz_zip_archive zip{};
    bool ok = false;
    explicit Reader(const QString &path, QString *error)
    {
        ok = mz_zip_reader_init_file(&zip, QFile::encodeName(path).constData(), 0);
        if (!ok && error)
            *error = QStringLiteral("cannot read %1 as a zip archive (%2)")
                         .arg(path, mz_zip_get_error_string(mz_zip_get_last_error(&zip)));
    }
    ~Reader() { if (ok) mz_zip_reader_end(&zip); }
};

} // namespace

QStringList names(const QString &archive, QString *error)
{
    QStringList out;
    Reader r(archive, error);
    if (!r.ok)
        return out;
    const mz_uint n = mz_zip_reader_get_num_files(&r.zip);
    for (mz_uint i = 0; i < n; i++) {
        mz_zip_archive_file_stat st;
        if (mz_zip_reader_file_stat(&r.zip, i, &st) && !mz_zip_reader_is_file_a_directory(&r.zip, i))
            out << QString::fromUtf8(st.m_filename);
    }
    return out;
}

bool extract(const QString &archive, const QString &member, const QString &dest, QString *error)
{
    Reader r(archive, error);
    if (!r.ok)
        return false;
    const int index = mz_zip_reader_locate_file(&r.zip, member.toUtf8().constData(), nullptr, 0);
    if (index < 0) {
        if (error) *error = QStringLiteral("%1 has no member %2").arg(archive, member);
        return false;
    }
    if (!mz_zip_reader_extract_to_file(&r.zip, mz_uint(index), QFile::encodeName(dest).constData(), 0)) {
        if (error)
            *error = QStringLiteral("extracting %1 failed (%2)")
                         .arg(member, mz_zip_get_error_string(mz_zip_get_last_error(&r.zip)));
        return false;
    }
    return true;
}

} // namespace ZipArchive
