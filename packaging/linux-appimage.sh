#!/usr/bin/env bash
# Build the Linux release: an AppImage of the desktop app (with gh-data and the
# engine inside it) and a portable tar.gz of the same tree. Run from the
# repository root after `make build`, or let it build. Needs linuxdeploy and
# its Qt plugin (downloaded into packaging/tools/ if absent) and, unless FUSE
# is available, nothing else -- the tools are run extracted.
#
#   packaging/linux-appimage.sh            -> dist/Genetic_Health-x86_64.AppImage
#                                             dist/genetic-health-linux-x86_64.tar.gz
set -euo pipefail
cd "$(dirname "$0")/.."
ARCH=$(uname -m)
TOOLS=packaging/tools
DIST=dist
mkdir -p "$TOOLS" "$DIST"

[ -x genetic-health-qt ] && [ -x gh-data ] && [ -x genetic-health ] || make build

for t in linuxdeploy-$ARCH.AppImage linuxdeploy-plugin-qt-$ARCH.AppImage; do
    if [ ! -x "$TOOLS/$t" ]; then
        repo=linuxdeploy; [[ $t == *plugin-qt* ]] && repo=linuxdeploy-plugin-qt
        curl -sL -o "$TOOLS/$t" "https://github.com/linuxdeploy/$repo/releases/download/continuous/$t"
        chmod +x "$TOOLS/$t"
    fi
done

rm -rf AppDir
export APPIMAGE_EXTRACT_AND_RUN=1
export OUTPUT="$DIST/Genetic_Health-$ARCH.AppImage"
# gh-data and the engine go in usr/bin beside the app, which is where the app
# looks for the engine first; linuxdeploy bundles the Qt libraries and
# platform plugins they need.
"$TOOLS/linuxdeploy-$ARCH.AppImage" --appdir AppDir \
    --executable genetic-health-qt --executable gh-data --executable genetic-health \
    --desktop-file packaging/genetic-health.desktop --icon-file packaging/genetic-health.png \
    --plugin qt --output appimage

# The same tree as a folder: usr/bin/{genetic-health-qt,gh-data,genetic-health}
# with the libraries in usr/lib, for running gh-data directly or without FUSE.
tar czf "$DIST/genetic-health-linux-$ARCH.tar.gz" -C AppDir usr
rm -rf AppDir
ls -la "$DIST"
