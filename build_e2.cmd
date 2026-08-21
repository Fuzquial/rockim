@echo off
rem build essai 2 : mouillage des joints en adoucissement (cle hydroWetDamage).
rem Nom d exe DEDIE : rockim_e1.exe reste intact pour comparaison.
rem l'exe d'un run en cours n'est jamais verrouille par la compilation.
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim_e2.exe
