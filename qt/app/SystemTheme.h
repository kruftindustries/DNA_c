// Follow the desktop's light/dark preference.
//
// Qt 5 picks up fonts and dialogs from the platform theme but, on most
// Linux desktops, not the colours: a dark GTK theme leaves a Qt 5 app light.
// This asks the desktop what it prefers -- the XDG settings portal
// (org.freedesktop.appearance color-scheme), then gsettings, then the GTK
// theme name; the registry on Windows -- and, when it prefers dark and the
// palette Qt arrived with is not, switches to the Fusion style with a dark
// palette. It keeps listening for changes so the app flips with the desktop.
// Qt 6.5+ reports the scheme itself and this defers to it. A user who
// configures Qt theming explicitly (QT_QPA_PLATFORMTHEME, QT_STYLE_OVERRIDE)
// is left alone.
#pragma once

#include <QObject>
#include <QPalette>

class QApplication;

class SystemTheme : public QObject {
    Q_OBJECT
public:
    enum class Scheme { Unknown, Light, Dark };

    explicit SystemTheme(QApplication &app);

    // What the desktop prefers right now.
    static Scheme preferred();
    static bool isDark(const QPalette &palette);
    static QPalette darkPalette();

private slots:
    void apply();

private:

    QApplication &m_app;
    QPalette m_initial;
    QString m_initialStyle;
    bool m_managed = true;   // false when the user configures Qt theming themselves
    Scheme m_applied = Scheme::Unknown;
};
