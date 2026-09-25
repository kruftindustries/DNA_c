@echo off
rem Build the Windows release: a folder with genetic-health-qt.exe, gh-data.exe,
rem genetic-health.exe and the Qt runtime, zipped. Run from the Qt MinGW
rem command prompt (Start menu -> Qt -> "Qt 6.8.x (MinGW 64-bit)") in the
rem repository root. Use a Qt 6 kit (6.8 MinGW 64-bit): its Schannel TLS
rem backend gives the app HTTPS without OpenSSL, which Qt 5 would need and
rem no longer distributes.
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
rem needs a known handful of them.
windeployqt --release --compiler-runtime --no-translations --no-opengl-sw --no-plugins "%OUT%\genetic-health-qt.exe" || exit /b 1
windeployqt --release --compiler-runtime --no-translations --no-plugins "%OUT%\gh-data.exe" || exit /b 1
for /f "delims=" %%p in ('qmake -query QT_INSTALL_PLUGINS') do set QT_PLUGINS=%%p
for /f "delims=" %%v in ('qmake -query QT_VERSION') do set QT_VERSION=%%v
mkdir "%OUT%\platforms" "%OUT%\styles" 2>nul
copy /y "%QT_PLUGINS%\platforms\qwindows.dll" "%OUT%\platforms\" >nul || exit /b 1
rem Native look: qwindowsvistastyle (Qt 5) or qmodernwindowsstyle (Qt 6.7+).
for %%f in ("%QT_PLUGINS%\styles\q*.dll") do copy /y "%%f" "%OUT%\styles\" >nul
rem Qt 6 does TLS through plugins: qschannelbackend needs nothing else.
if exist "%QT_PLUGINS%\tls" (
    mkdir "%OUT%\tls" 2>nul
    for %%f in ("%QT_PLUGINS%\tls\q*.dll") do copy /y "%%f" "%OUT%\tls\" >nul
)

rem Qt 5 instead loads OpenSSL 1.1 at run time and does not ship it. If this
rem is a Qt 5 kit, look for the DLLs where the Qt installer's OpenSSL toolkit
rem put them, or in OPENSSL_BIN.
if "%QT_VERSION:~0,2%"=="5." (
    for /f "delims=" %%p in ('qmake -query QT_INSTALL_PREFIX') do set QT_PREFIX=%%p
    if not defined OPENSSL_BIN set OPENSSL_BIN=%QT_PREFIX%\..\..\Tools\OpenSSL\Win_x64\bin
    if exist "%OPENSSL_BIN%\libssl-1_1-x64.dll" (
        copy /y "%OPENSSL_BIN%\libssl-1_1-x64.dll" "%OUT%\" >nul
        copy /y "%OPENSSL_BIN%\libcrypto-1_1-x64.dll" "%OUT%\" >nul
    ) else (
        echo WARNING: Qt 5 kit and no OpenSSL 1.1 DLLs at %OPENSSL_BIN%; HTTPS downloads will not work. Build with a Qt 6 kit instead.
    )
)
dir "%OUT%"
if exist "%OUT%.zip" del "%OUT%.zip"
powershell -NoProfile -Command "Compress-Archive -Path '%OUT%\*' -DestinationPath '%OUT%.zip'" || exit /b 1
dir dist
endlocal
