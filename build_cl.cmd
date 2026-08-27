@echo off
rem build du correctif jointShearRange = coulomb (2026-08-27, commits
rem 4cbd74f..3ab1748). Nom d exe DEDIE : rockim_cl.exe — l ancien binaire
rem du run du 26/08 reste disponible pour comparaison.
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim_cl.exe
