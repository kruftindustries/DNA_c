#include "ToolLocator.h"

#include <QCoreApplication>

#include <QFileInfo>
#include <QStandardPaths>

namespace ToolLocator {

QString find(const QString &tool)
{
    return QStandardPaths::findExecutable(tool);
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
    if (!extra.isEmpty())
        env.insert(QStringLiteral("PATH"),
                   extra.join(QLatin1Char(':')) + ':' + env.value(QStringLiteral("PATH")));
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
    return QStringLiteral("brew install minimap2 samtools bcftools fastp");
#elif defined(Q_OS_WIN)
    return QStringLiteral("in WSL2 (Ubuntu): sudo apt install minimap2 samtools bcftools fastp");
#else
    return QStringLiteral("sudo apt install minimap2 samtools bcftools fastp");
#endif
}

QString analysisBinary(const QString &root)
{
#ifdef Q_OS_WIN
    const QString exe = QStringLiteral("genetic-health.exe");
#else
    const QString exe = QStringLiteral("genetic-health");
#endif
    QStringList candidates{QCoreApplication::applicationDirPath() + '/' + exe};
    if (!root.isEmpty())
        candidates << root + '/' + exe << root + "/c/bin/" + exe;
    for (const QString &c : candidates)
        if (QFileInfo(c).isExecutable())
            return c;
    return QStandardPaths::findExecutable(QStringLiteral("genetic-health"));
}

} // namespace ToolLocator
