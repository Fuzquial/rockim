@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1
if not exist obj_ENER mkdir obj_ENER
cl /nologo /std:c++17 /EHsc /O2 /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\rockim\eigen-3.4.0 /Fo:obj_ENER\ /Fe:obj_ENER\bench_energie%1.exe tests_f2\bench_energie_matrice.cpp src\MatLaw.cpp src\MatLawDfhPlus.cpp src\Config.cpp
