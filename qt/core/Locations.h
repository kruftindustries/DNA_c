// Where the app keeps its data. In a source checkout -- the executable
// sits in it, or the current directory is in it -- data/, reference/ and
// reports/ live at the checkout root. A packaged build (AppImage, Windows
// zip, macOS bundle) has no checkout, so the same layout goes under the
// per-user application data directory: ~/.local/share/genetic-health on
// Linux, %LOCALAPPDATA%\genetic-health on Windows, ~/Library/Application
// Support/genetic-health on macOS. The app, gh-data and the Qt tests all
// resolve it here so they can never disagree.
#pragma once

#include <QString>

namespace Locations {

// The checkout containing `startDir` (c/src/main.c marks it), or empty.
QString findCheckout(const QString &startDir);

// The working root as described above; the data directory is created when
// it has to fall back to application data.
QString defaultRoot();

// True when defaultRoot() is the application data directory, not a checkout.
bool isPackaged();

} // namespace Locations
