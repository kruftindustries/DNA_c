# The analysis engine as a qmake project, so one build (the top-level
# genetic-health.pro) produces the engine, the desktop app and gh-data on
# every platform. c/Makefile remains the Unix build for the engine and its
# unit tests.
TEMPLATE = app
TARGET = genetic-health
CONFIG += console warn_on
CONFIG -= qt app_bundle
DESTDIR = $$PWD/..

INCLUDEPATH += $$PWD/include
SOURCES += $$files($$PWD/src/*.c)

*-g++*|*-clang*: QMAKE_CFLAGS += -std=c11 -O2
win32-g++: DEFINES += __USE_MINGW_ANSI_STDIO=1
msvc: QMAKE_CFLAGS += /std:c11
unix: LIBS += -lm
