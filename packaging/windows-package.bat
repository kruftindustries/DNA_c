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
rem The Qt DLLs and the MinGW runtime for both Qt executables. Plugins are
rem copied by hand: windeployqt's release/debug detection rejects the MinGW
rem kit's platform plugin ("Unable to find the platform plugin"), and the app
rem needs exactly two of them.
windeployqt --release --compiler-runtime --no-translations --no-opengl-sw --no-plugins "%OUT%\genetic-health-qt.exe" || exit /b 1
windeployqt --release --compiler-runtime --no-translations --no-plugins "%OUT%\gh-data.exe" || exit /b 1
for /f "delims=" %%p in ('qmake -query QT_INSTALL_PLUGINS') do set QT_PLUGINS=%%p
mkdir "%OUT%\platforms" "%OUT%\styles" 2>nul
copy /y "%QT_PLUGINS%\platforms\qwindows.dll" "%OUT%\platforms\" >nul || exit /b 1
copy /y "%QT_PLUGINS%\styles\qwindowsvistastyle.dll" "%OUT%\styles\" >nul

rem HTTPS (ClinVar, ClinPGx, Ensembl, 1000 Genomes) needs OpenSSL 1.1, which
rem Qt 5 loads at run time and does not ship. The Qt online installer's
rem "OpenSSL 1.1.1 Toolkit" (and install-qt-action's tools_openssl_x64) puts
rem it under Qt\Tools\OpenSSL\Win_x64\bin.
for /f "delims=" %%p in ('qmake -query QT_INSTALL_PREFIX') do set QT_PREFIX=%%p
set OPENSSL_BIN=%QT_PREFIX%\..\..\Tools\OpenSSL\Win_x64\bin
if defined IQTA_TOOLS set OPENSSL_BIN=%IQTA_TOOLS%\OpenSSL\Win_x64\bin
if exist "%OPENSSL_BIN%\libssl-1_1-x64.dll" (
    copy /y "%OPENSSL_BIN%\libssl-1_1-x64.dll" "%OUT%\" >nul
    copy /y "%OPENSSL_BIN%\libcrypto-1_1-x64.dll" "%OUT%\" >nul
) else (
    echo WARNING: OpenSSL 1.1 DLLs not found at %OPENSSL_BIN%; HTTPS downloads will not work without them.
)
dir "%OUT%"
if exist "%OUT%.zip" del "%OUT%.zip"
powershell -NoProfile -Command "Compress-Archive -Path '%OUT%\*' -DestinationPath '%OUT%.zip'" || exit /b 1
dir dist
endlocal
