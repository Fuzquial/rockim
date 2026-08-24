@echo off
rem Build de la capacite DP-DFH en FEM pur : sortie de l endommagement (2026-08-24).
rem Nom d'exe DEDIE (rockim_dfh.exe) : le rockim_p1 et ses runs en cours ne
rem sont jamais touches, et l'exe d'un calcul actif n'est pas verrouille.
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim_dfh.exe
