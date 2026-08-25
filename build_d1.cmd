@echo off
rem ---------------------------------------------------------------------------
rem Chantier "DIF intrinseque" (point 1 du tableau de comparaison aux articles).
rem
rem   build_d0.cmd -> rockim_d0.exe : la REFERENCE, code non modifie (main).
rem   build_d1.cmd -> rockim_d1.exe : APRES l ajout strainRateDIFArm = envelope.
rem
rem La preuve de bit-neutralite (principe I) se fait en comparant les sorties
rem des deux exe sur la suite fast : 19/19 aux MEMES chiffres.
rem ---------------------------------------------------------------------------
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim_d1.exe
