@echo off
rem Build the Windows release: a folder with genetic-health-qt.exe, gh-data.exe,
rem genetic-health.exe and the Qt runtime, zipped. Run from the Qt MinGW
rem command prompt (Start menu -> Qt -> "Qt 5.15.x (MinGW 64-bit)") in the
rem repository root:
rem
rem   packaging\windows-package.bat        -> dist\genetic-health-windows-x64.zip
setlocal
cd /d "%~dp0\.."
if not exist build mkdir build
pushd build
qmake ..\genetic-health.pro || exit /b 1
mingw32-make -j%NUMBER_OF_PROCESSORS% || exit /b 1
popd

set OUT=dist\genetic-health-windows-x64
if exist "%OUT%" rmdir /s /q "%OUT%"
mkdir "%OUT%"
copy /y genetic-health-qt.exe "%OUT%\" >nul
copy /y gh-data.exe "%OUT%\" >nul
copy /y genetic-health.exe "%OUT%\" >nul
copy /y README.md "%OUT%\" >nul
rem The Qt DLLs, platform plugins and the MinGW runtime, for both Qt executables.
rem windeployqt locates the plugins through `qmake -query`; point it at the
rem kit's own qmake explicitly (QT_ROOT_DIR is set by the Qt online installer
rem prompt and by install-qt-action) so a stray qmake on PATH cannot mislead it.
set QMAKE_OPT=
if defined QT_ROOT_DIR set QMAKE_OPT=--qmake "%QT_ROOT_DIR%\bin\qmake.exe"
where qmake
qmake -query QT_INSTALL_PLUGINS
if defined QT_ROOT_DIR dir "%QT_ROOT_DIR%\plugins\platforms"
windeployqt --verbose 2 %QMAKE_OPT% --release --compiler-runtime --no-translations --no-opengl-sw "%OUT%\genetic-health-qt.exe" || exit /b 1
windeployqt %QMAKE_OPT% --release --compiler-runtime --no-translations "%OUT%\gh-data.exe" || exit /b 1
if exist "%OUT%.zip" del "%OUT%.zip"
powershell -NoProfile -Command "Compress-Archive -Path '%OUT%\*' -DestinationPath '%OUT%.zip'" || exit /b 1
dir dist
endlocal
