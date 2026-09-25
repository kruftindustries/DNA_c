#!/usr/bin/env bash
# Build the macOS release: Genetic Health.app with gh-data and the engine
# inside its Contents/MacOS, wrapped in a dmg by macdeployqt. Run from the
# repository root with Qt's bin on PATH (qmake, macdeployqt), e.g. after
# `brew install qt@5` or from a Qt online-installer kit.
#
#   packaging/macos-dmg.sh                 -> dist/genetic-health-macos-<arch>.dmg
#
# The result is unsigned: a first launch needs right-click -> Open (or
# `xattr -dr com.apple.quarantine "Genetic Health.app"`).
set -euo pipefail
cd "$(dirname "$0")/.."
ARCH=$(uname -m)
DIST=dist
mkdir -p build "$DIST"
( cd build && qmake ../genetic-health.pro && make -j"$(sysctl -n hw.ncpu)" )

APP=genetic-health-qt.app
[ -d "$APP" ] || { echo "expected $APP in the repository root"; exit 1; }
cp gh-data genetic-health "$APP/Contents/MacOS/"
macdeployqt "$APP" -executable="$APP/Contents/MacOS/gh-data" -executable="$APP/Contents/MacOS/genetic-health-qt"
rm -rf "$DIST/Genetic Health.app"
cp -R "$APP" "$DIST/Genetic Health.app"
rm -f "$DIST/genetic-health-macos-$ARCH.dmg"
hdiutil create -volname "Genetic Health" -srcfolder "$DIST/Genetic Health.app" -ov -format UDZO "$DIST/genetic-health-macos-$ARCH.dmg"
rm -rf "$DIST/Genetic Health.app"
ls -la "$DIST"
