@echo off
set ProductDir=
set toolchain=
set arch=
set archivetag=HEAD

if "%~1"=="-t" goto TOOLCHAIN
goto BEGINVCVAR
:TOOLCHAIN
if "%~2"=="M32" goto M32
if "%~2"=="m32" goto M32
if "%~2"=="M64" goto M64
if "%~2"=="m64" goto M64
set toolchain="%~f2"
goto ENDTOOLCHAIN
:M32
set arch=x86
goto ENDTOOLCHAIN
:M64
set arch=amd64
:ENDTOOLCHAIN
shift
shift

rem Find Visual Studio Install
:BEGINVCVAR
if DEFINED VCINSTALLDIR goto ENDVCVAR
if NOT DEFINED PROCESSOR_ARCHITECTURE goto X86
if /I %PROCESSOR_ARCHITECTURE% == x86 goto X86
:AMD64
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\WOW6432Node\Microsoft\VCExpress\8.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\WOW6432Node\Microsoft\VisualStudio\8.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\WOW6432Node\Microsoft\VCExpress\9.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\WOW6432Node\Microsoft\VisualStudio\9.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\WOW6432Node\Microsoft\VisualStudio\10.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\WOW6432Node\Microsoft\VisualStudio\11.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
goto REST
:X86
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\Microsoft\VCExpress\8.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\Microsoft\VisualStudio\8.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\Microsoft\VCExpress\9.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\Microsoft\VisualStudio\9.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\Microsoft\VisualStudio\10.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
FOR /F "tokens=2*" %%A IN ('REG.EXE QUERY "HKLM\Software\Microsoft\VisualStudio\11.0\Setup\VC" /V "ProductDir" 2^>NUL ^| FIND "REG_SZ"') DO SET ProductDir=%%B
:REST
if NOT DEFINED ProductDir goto ENDVCVAR
:VCVAR
call "%ProductDir%\vcvarsall.bat" %arch%
:ENDVCVAR

if NOT [%1]==[] set archivetag=%1
set buildargs=-DRELENG_TAG=%archivetag%
if DEFINED toolchain set buildargs=%buildargs% -DRELENG_TOOLCHAIN=%toolchain%

cmake %buildargs% -P releng.cmake

:EOF
