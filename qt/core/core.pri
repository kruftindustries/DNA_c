# Non-GUI core: data updates and the WGS pipeline. QtCore + QtNetwork only,
# shared by the GUI, the gh-data console tool and the tests.
QT += core network
INCLUDEPATH += $$PWD

RESOURCES += $$PWD/core.qrc

HEADERS += \
    $$PWD/Reporter.h \
    $$PWD/ToolLocator.h \
    $$PWD/TextTable.h \
    $$PWD/GzipStream.h \
    $$PWD/DataVersions.h \
    $$PWD/Downloader.h \
    $$PWD/ClinVarUpdater.h \
    $$PWD/PharmgkbUpdater.h \
    $$PWD/EnsemblLookup.h \
    $$PWD/TargetRegions.h \
    $$PWD/VcfConverter.h \
    $$PWD/WgsPipeline.h \
    $$PWD/WgsValidator.h \
    $$PWD/TestGenome.h \
    $$PWD/ChainLift.h \
    $$PWD/RemoteFasta.h \
    $$PWD/ZipArchive.h \
    $$PWD/ReferenceGenome.h \
    $$PWD/Locations.h

SOURCES += \
    $$PWD/Reporter.cpp \
    $$PWD/ToolLocator.cpp \
    $$PWD/TextTable.cpp \
    $$PWD/GzipStream.cpp \
    $$PWD/DataVersions.cpp \
    $$PWD/Downloader.cpp \
    $$PWD/ClinVarUpdater.cpp \
    $$PWD/PharmgkbUpdater.cpp \
    $$PWD/EnsemblLookup.cpp \
    $$PWD/TargetRegions.cpp \
    $$PWD/VcfConverter.cpp \
    $$PWD/WgsPipeline.cpp \
    $$PWD/WgsValidator.cpp \
    $$PWD/TestGenome.cpp \
    $$PWD/ChainLift.cpp \
    $$PWD/RemoteFasta.cpp \
    $$PWD/ZipArchive.cpp \
    $$PWD/ReferenceGenome.cpp \
    $$PWD/Locations.cpp \
    $$PWD/third_party/miniz/miniz.c

# miniz (third_party/miniz, MIT) is plain C compiled with the C compiler;
# the project's C++ warning flags do not apply to it.
DEFINES += MINIZ_NO_DEFLATE_APIS   # miniz.h derives MINIZ_NO_ARCHIVE_WRITING_APIS from it
*-g++*|*-clang*: QMAKE_CFLAGS += -w
