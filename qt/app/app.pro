include(../common.pri)
TEMPLATE = app
# The desktop's colour-scheme preference comes over D-Bus (XDG portal) on Linux.
unix:!macx: qtHaveModule(dbus): QT += dbus
TARGET = genetic-health-qt
DESTDIR = $$PWD/../..

HEADERS += \
    MainWindow.h InputPage.h RunPage.h ReportPage.h FindingsPage.h \
    DataSourcesPage.h WgsPage.h JobRunner.h SystemTheme.h

SOURCES += \
    main.cpp MainWindow.cpp InputPage.cpp RunPage.cpp ReportPage.cpp \
    FindingsPage.cpp DataSourcesPage.cpp WgsPage.cpp JobRunner.cpp SystemTheme.cpp
