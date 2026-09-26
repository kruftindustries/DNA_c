#include "ToolLocator.h"

#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QStandardPaths>

namespace ToolLocator {

namespace {

QString exeName(const QString &tool)
{
#ifdef Q_OS_WIN
    return tool + QStringLiteral(".exe");
#else
    return tool;
#endif
}

// Where package managers put binaries that a GUI launched from the desktop
// does not have on PATH (macOS in particular hands Finder-launched apps
// only the system directories).
QStringList packagePrefixDirs()
{
    QStringList dirs;
#ifndef Q_OS_WIN
    dirs << QStringLiteral("/opt/homebrew/bin") << QStringLiteral("/usr/local/bin")
         << QStringLiteral("/opt/local/bin") << QDir::homePath() + QStringLiteral("/.local/bin");
#endif
    return dirs;
}

} // namespace

QString bundledDir()
{
    return QCoreApplication::applicationDirPath() + QStringLiteral("/tools");
}

QString find(const QString &tool)
{
    const QFileInfo bundled(bundledDir() + '/' + exeName(tool));
    if (bundled.isExecutable())
        return bundled.absoluteFilePath();
    const QString onPath = QStandardPaths::findExecutable(tool);
    if (!onPath.isEmpty())
        return onPath;
    return QStandardPaths::findExecutable(tool, packagePrefixDirs());
}

QStringList extraPathDirs(const QStringList &tools)
{
    QStringList dirs;
    for (const QString &t : tools) {
        const QString path = find(t);
        if (path.isEmpty())
            continue;
        const QString dir = QFileInfo(path).absolutePath();
        if (!dirs.contains(dir))
            dirs << dir;
    }
    return dirs;
}

QProcessEnvironment environmentFor(const QStringList &tools)
{
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    const QStringList extra = extraPathDirs(tools);
    if (!extra.isEmpty()) {
        const QString sep = QDir::listSeparator();
        env.insert(QStringLiteral("PATH"),
                   QDir::toNativeSeparators(extra.join(sep)) + sep + env.value(QStringLiteral("PATH")));
    }
    return env;
}

QStringList wgsTools()
{
    return {QStringLiteral("minimap2"), QStringLiteral("samtools"),
            QStringLiteral("bcftools"), QStringLiteral("fastp")};
}

QStringList requiredWgsTools()
{
    return {QStringLiteral("minimap2"), QStringLiteral("samtools"), QStringLiteral("bcftools")};
}

QStringList missing(const QStringList &tools)
{
    QStringList out;
    for (const QString &t : tools)
        if (find(t).isEmpty())
            out << t;
    return out;
}

QString installHint()
{
#if defined(Q_OS_MACOS)
    return QStringLiteral("Install them with Homebrew: brew install minimap2 samtools bcftools "
                          "(fastp is optional: brew install fastp)");
#elif defined(Q_OS_WIN)
    return QStringLiteral("The Windows release ships them in the tools folder beside "
                          "genetic-health-qt.exe (%1); extract the whole zip, or copy that folder "
                          "from it. Without fastp, read QC is skipped.")
        .arg(QDir::toNativeSeparators(bundledDir()));
#else
    return QStringLiteral("Install them with: sudo apt install minimap2 samtools bcftools "
                          "(fastp is optional: sudo apt install fastp)");
#endif
}

QString analysisBinary(const QString &root)
{
    const QString exe = exeName(QStringLiteral("genetic-health"));
    QStringList candidates{QCoreApplication::applicationDirPath() + '/' + exe};
    if (!root.isEmpty())
        candidates << root + '/' + exe << root + "/c/bin/" + exe;
    for (const QString &c : candidates)
        if (QFileInfo(c).isExecutable())
            return c;
    return QStandardPaths::findExecutable(QStringLiteral("genetic-health"));
}

} // namespace ToolLocator
