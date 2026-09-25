#include "Locations.h"

#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QStandardPaths>

namespace Locations {

QString findCheckout(const QString &startDir)
{
    QString dir = startDir;
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

namespace {

QString appDataRoot()
{
    // GenericDataLocation rather than AppDataLocation so the app and
    // gh-data (different application names) share one directory.
    const QString base = QStandardPaths::writableLocation(QStandardPaths::GenericDataLocation);
    const QString root = base + "/genetic-health";
    QDir().mkpath(root + "/data");
    return root;
}

} // namespace

QString defaultRoot()
{
    QString found = findCheckout(QCoreApplication::applicationDirPath());
    if (found.isEmpty())
        found = findCheckout(QDir::currentPath());
    return found.isEmpty() ? appDataRoot() : found;
}

bool isPackaged()
{
    return findCheckout(QCoreApplication::applicationDirPath()).isEmpty()
        && findCheckout(QDir::currentPath()).isEmpty();
}

} // namespace Locations
