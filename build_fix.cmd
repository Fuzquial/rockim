@echo off
rem Build de travail pour le lot de corrections du guide du 2026-09-06.
rem Eigen : la copie locale est dans ..\rockim\eigen-3.4.0 (celle de CMakeLists),
rem PAS dans ..\eigen-3.4.0 comme l'ecrivent build_dev.cmd / build_chk.cmd.
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
if not exist obj_fix mkdir obj_fix
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX ^
   /I include /I ..\rockim\eigen-3.4.0 /Fo:obj_fix\ src\*.cpp /Fe:rockim_fix.exe
