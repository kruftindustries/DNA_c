#include "RepoPaths.h"
#include "ToolLocator.h"

#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QSettings>
#include <QStandardPaths>

namespace {

QString findRootFrom(QString dir)
{
    for (int depth = 0; depth < 8; ++depth) {
        if (QFileInfo(dir + "/c/src/main.c").exists())
            return QDir(dir).absolutePath();
        QDir d(dir);
        if (!d.cdUp())
            break;
        dir = d.absolutePath();
    }
    return QString();
}

QString setting(const char *key)
{
    return QSettings().value(QString::fromLatin1(key)).toString();
}

} // namespace

namespace RepoPaths {

QString root()
{
    const QString override = setting("paths/root");
    if (!override.isEmpty())
        return override;
    QString found = findRootFrom(QCoreApplication::applicationDirPath());
    if (found.isEmpty())
        found = findRootFrom(QDir::currentPath());
    return found;
}

QString dataDir()
{
    const QString override = setting("paths/data");
    if (!override.isEmpty())
        return override;
    return root().isEmpty() ? QString() : root() + "/data";
}

QString analysisBinary()
{
    const QString override = setting("paths/binary");
    if (!override.isEmpty())
        return override;
    return ToolLocator::analysisBinary(root());
}

void setDataDir(const QString &path)        { QSettings().setValue("paths/data", path); }
void setAnalysisBinary(const QString &path) { QSettings().setValue("paths/binary", path); }

} // namespace RepoPaths
