# Builds everything: the C analysis engine, the Qt desktop app, the gh-data
# console tool and the Qt unit tests. The three executables land in this
# directory.
#
#   mkdir build && cd build && qmake .. && make        (Linux, macOS)
#   mkdir build && cd build && qmake .. && mingw32-make (Windows, Qt MinGW kit)
#
# Or open this file in Qt Creator and press Build.
TEMPLATE = subdirs
SUBDIRS = engine app cli tests
engine.file = c/engine.pro
app.file = qt/app/app.pro
cli.file = qt/cli/cli.pro
tests.file = qt/tests/tests.pro
app.depends = engine
cli.depends = engine
tests.depends = app cli
