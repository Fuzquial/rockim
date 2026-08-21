@echo off
rem build du correctif du signe hydro (AbuAisha 2017) + controle de fermeture.
rem Nom d'exe DEDIE : l'ancien binaire reste disponible pour comparaison, et
rem l'exe d'un run en cours n'est jamais verrouille par la compilation.
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim_hy.exe
