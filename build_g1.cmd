@echo off
rem Build de rockim_g1 — lot « loi de la note FDEM hybride a insertion adaptative ».
rem Eigen : la copie locale est dans ..\rockim\eigen-3.4.0 (celle de CMakeLists).
rem Usage : build_g1.cmd            -> rockim_g1.exe  (objets dans obj_g1)
rem         build_g1.cmd ref        -> rockim_g1ref.exe (objets dans obj_g1ref)
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
set TAG=%1
if "%TAG%"=="" set TAG=g1
if not exist obj_%TAG% mkdir obj_%TAG%
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX ^
   /I include /I ..\rockim\eigen-3.4.0 /Fo:obj_%TAG%\ src\*.cpp /Fe:rockim_%TAG%.exe
