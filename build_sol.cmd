@echo off
rem Build de la branche joint-handoff : les conventions RELEVEES DANS LE CODE
rem SOURCE de Solidity (ImperialCollegeLondon/solidity-solver-open, LGPL-3.0).
rem Exe DISTINCT : rockim_dev.exe et les autres restent intacts.
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim_sol.exe
