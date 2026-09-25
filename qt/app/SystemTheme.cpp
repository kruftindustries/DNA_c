#include "SystemTheme.h"

#include <QApplication>
#include <QProcess>
#include <QProcessEnvironment>
#include <QSettings>
#include <QStyle>
#include <QStyleFactory>
#include <QStyleHints>
#include <QTimer>

#ifdef QT_DBUS_LIB
#include <QDBusConnection>
#include <QDBusInterface>
#include <QDBusReply>
#include <QDBusVariant>
#endif

namespace {

#ifdef QT_DBUS_LIB
// The portal's answer, or Unknown when there is no portal or no setting.
SystemTheme::Scheme fromPortal()
{
    QDBusInterface portal("org.freedesktop.portal.Desktop", "/org/freedesktop/portal/desktop",
                          "org.freedesktop.portal.Settings", QDBusConnection::sessionBus());
    if (!portal.isValid())
        return SystemTheme::Scheme::Unknown;
    QDBusReply<QDBusVariant> reply = portal.call("Read", "org.freedesktop.appearance", "color-scheme");
    if (!reply.isValid())
        return SystemTheme::Scheme::Unknown;
    QVariant v = reply.value().variant();
    while (v.canConvert<QDBusVariant>())          // Read wraps the value twice
        v = v.value<QDBusVariant>().variant();
    switch (v.toUInt()) {
    case 1:  return SystemTheme::Scheme::Dark;
    case 2:  return SystemTheme::Scheme::Light;
    default: return SystemTheme::Scheme::Unknown;   // 0: no preference
    }
}
#endif

QString gsetting(const QString &key)
{
    QProcess p;
    p.start("gsettings", {"get", "org.gnome.desktop.interface", key});
    if (!p.waitForFinished(2000) || p.exitCode() != 0)
        return QString();
    return QString::fromUtf8(p.readAllStandardOutput()).trimmed().remove('\'');
}

} // namespace

SystemTheme::Scheme SystemTheme::preferred()
{
#if QT_VERSION >= QT_VERSION_CHECK(6, 5, 0)
    switch (QGuiApplication::styleHints()->colorScheme()) {
    case Qt::ColorScheme::Dark:  return Scheme::Dark;
    case Qt::ColorScheme::Light: return Scheme::Light;
    default: break;
    }
#endif
#if defined(Q_OS_WIN)
    QSettings reg("HKEY_CURRENT_USER\\Software\\Microsoft\\Windows\\CurrentVersion\\Themes\\Personalize",
                  QSettings::NativeFormat);
    const QVariant light = reg.value("AppsUseLightTheme");
    if (light.isValid())
        return light.toInt() == 0 ? Scheme::Dark : Scheme::Light;
    return Scheme::Unknown;
#elif defined(Q_OS_MACOS)
    return Scheme::Unknown;   // Qt follows macOS appearance on its own
#else
#ifdef QT_DBUS_LIB
    const Scheme portal = fromPortal();
    if (portal != Scheme::Unknown)
        return portal;
#endif
    const QString scheme = gsetting("color-scheme");
    if (scheme == "prefer-dark")
        return Scheme::Dark;
    if (scheme == "prefer-light")
        return Scheme::Light;
    const QString theme = QProcessEnvironment::systemEnvironment().value("GTK_THEME", gsetting("gtk-theme"));
    if (!theme.isEmpty())
        return theme.contains("dark", Qt::CaseInsensitive) ? Scheme::Dark : Scheme::Light;
    return Scheme::Unknown;
#endif
}

bool SystemTheme::isDark(const QPalette &palette)
{
    return palette.color(QPalette::Window).lightness() < 128;
}

QPalette SystemTheme::darkPalette()
{
    QPalette p;
    const QColor window(0x2b, 0x2b, 0x2b), base(0x1e, 0x1e, 0x1e), alt(0x26, 0x26, 0x26);
    const QColor text(0xe6, 0xe6, 0xe6), disabled(0x80, 0x80, 0x80), highlight(0x3d, 0x8e, 0xc9);
    p.setColor(QPalette::Window, window);
    p.setColor(QPalette::WindowText, text);
    p.setColor(QPalette::Base, base);
    p.setColor(QPalette::AlternateBase, alt);
    p.setColor(QPalette::ToolTipBase, window);
    p.setColor(QPalette::ToolTipText, text);
    p.setColor(QPalette::Text, text);
    p.setColor(QPalette::Button, window);
    p.setColor(QPalette::ButtonText, text);
    p.setColor(QPalette::BrightText, Qt::red);
    p.setColor(QPalette::Link, highlight.lighter(120));
    p.setColor(QPalette::LinkVisited, highlight.lighter(140));
    p.setColor(QPalette::Highlight, highlight);
    p.setColor(QPalette::HighlightedText, Qt::white);
    p.setColor(QPalette::PlaceholderText, disabled);
    p.setColor(QPalette::Mid, QColor(0x9a, 0x9a, 0x9a));   // "color: palette(mid)" hints stay readable
    p.setColor(QPalette::Dark, QColor(0x15, 0x15, 0x15));
    p.setColor(QPalette::Light, QColor(0x40, 0x40, 0x40));
    for (QPalette::ColorRole role : {QPalette::WindowText, QPalette::Text, QPalette::ButtonText}) {
        p.setColor(QPalette::Disabled, role, disabled);
    }
    p.setColor(QPalette::Disabled, QPalette::Highlight, QColor(0x50, 0x50, 0x50));
    return p;
}

SystemTheme::SystemTheme(QApplication &app) : QObject(&app), m_app(app)
{
    const QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    m_managed = env.value("QT_QPA_PLATFORMTHEME").isEmpty() && env.value("QT_STYLE_OVERRIDE").isEmpty();
    m_initial = app.palette();
    m_initialStyle = app.style()->objectName();
#if QT_VERSION >= QT_VERSION_CHECK(6, 5, 0)
    // Qt propagates the scheme into the palette itself from 6.5 on.
    m_managed = false;
#endif
    if (!m_managed)
        return;
    apply();
#ifdef QT_DBUS_LIB
    QDBusConnection::sessionBus().connect(
        "org.freedesktop.portal.Desktop", "/org/freedesktop/portal/desktop",
        "org.freedesktop.portal.Settings", "SettingChanged", this, SLOT(apply()));
#endif
    // Desktops without a portal: notice a change within a few seconds.
    auto *poll = new QTimer(this);
    connect(poll, &QTimer::timeout, this, &SystemTheme::apply);
    poll->start(5000);
}

void SystemTheme::apply()
{
    const Scheme want = preferred();
    if (want == m_applied)
        return;
    if (want == Scheme::Dark && !isDark(m_initial)) {
        // The platform theme did not deliver the dark colours; Fusion is the
        // one style that draws every widget from the palette.
        if (m_app.style()->objectName().compare("fusion", Qt::CaseInsensitive) != 0)
            m_app.setStyle(QStyleFactory::create("Fusion"));
        m_app.setPalette(darkPalette());
    } else if (want != Scheme::Dark) {
        // Back to what the platform gave us. The style stays Fusion once
        // switched: swapping styles under live widgets is not worth it.
        if (m_applied == Scheme::Dark)
            m_app.setPalette(m_initial);
    }
    m_applied = want;
}
