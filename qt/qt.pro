# Qt front end for the genetic health pipeline.
#
#   mkdir -p build && cd build && qmake .. && make -j
#   QT_QPA_PLATFORM=offscreen ./tests/tst_genomeformat   # etc.
#
# Qt 5.15 or Qt 6. WebEngine is optional: with it the report renders in-app,
# without it the report pane shows a text rendering and opens the full report
# in the system browser.
TEMPLATE = subdirs
SUBDIRS = app cli tests
tests.depends = app cli
