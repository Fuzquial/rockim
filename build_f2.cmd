@echo off
rem build_f2.cmd — chantier f2 : preBrokenJoints + instrumentation microsismique.
rem Chaine identique a celle de la campagne p1 (MSVC + Eigen en-tetes du dossier rockim/).
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
if "%~1"=="" (set OUT=rockim_f2.exe) else (set OUT=%~1)
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\rockim\eigen-3.4.0 src\*.cpp /Fe:%OUT%
