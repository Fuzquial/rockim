@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1
if not exist obj_LOTA mkdir obj_LOTA
cl /nologo /std:c++17 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /c /I include /I ..\rockim\eigen-3.4.0 /Fo:obj_LOTA\ src\MatLaw.cpp
