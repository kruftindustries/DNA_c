include(../core/core.pri)
QT -= gui
CONFIG += c++17 warn_on console
CONFIG -= app_bundle
TEMPLATE = app
TARGET = gh-data
VERSION = 0.2.3
DESTDIR = $$PWD/../..
SOURCES += main.cpp
