QT += core gui widgets network
CONFIG += c++17 warn_on
DEFINES += QT_DEPRECATED_WARNINGS
INCLUDEPATH += $$PWD/app
include($$PWD/core/core.pri)

qtHaveModule(webenginewidgets) {
    QT += webenginewidgets
    DEFINES += HAVE_WEBENGINE
}

# Sources shared by the app and the tests (everything but main.cpp).
HEADERS += \
    $$PWD/app/GenomeFormat.h \
    $$PWD/app/RepoPaths.h \
    $$PWD/app/ProcessRunner.h \
    $$PWD/app/AnalysisRunner.h \
    $$PWD/app/FindingsModel.h

SOURCES += \
    $$PWD/app/GenomeFormat.cpp \
    $$PWD/app/RepoPaths.cpp \
    $$PWD/app/ProcessRunner.cpp \
    $$PWD/app/AnalysisRunner.cpp \
    $$PWD/app/FindingsModel.cpp
