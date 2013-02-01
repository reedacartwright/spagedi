@echo off

call "C:\Program Files\Microsoft Visual Studio 9.0\VC\vcvarsall.bat" x86

set MAKE=nmake
set CMAKE=cmake
set SVN=svn
set REPOS=svn://scit.us/documents/projects/spagedi

set RELENG_DIR="%TEMP%\spagedi-releng.%RANDOM%"
mkdir %RELENG_DIR% || exit /B 1

set DEST_DIR="%CD%"
set SOURCE_DIR="%RELENG_DIR%\source"
set BUILD_DIR="%RELENG_DIR%\build"

%SVN% co -q %REPOS% %SOURCE_DIR% || exit /B 1

mkdir %BUILD_DIR% || exit /B 1
cd %BUILD_DIR% || exit /B 1

%CMAKE% -G "NMake Makefiles" %SOURCE_DIR% -DCMAKE_BUILD_TYPE=Release
%MAKE%
%MAKE% package
%MAKE% package_source

xcopy /Q /Y SPAGeDi-1* %DEST_DIR%

cd %DEST_DIR%
rd /S /Q %RELENG_DIR%
