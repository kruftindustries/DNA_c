# Shared by every tst_*.pro, after the test has included core.pri or
# common.pri. All the tests live in this one directory and each compiles the
# core sources itself, so without private object and moc directories a
# parallel build (make -j) has several of them writing the same .o files at
# once and fails at link time.
QT += testlib
QT -= gui
CONFIG += testcase console
TEMPLATE = app
OBJECTS_DIR = .obj/$$TARGET
MOC_DIR = .moc/$$TARGET
RCC_DIR = .rcc/$$TARGET
