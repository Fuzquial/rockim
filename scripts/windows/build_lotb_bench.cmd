@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
set SP=C:\Users\FUZQUI~1\AppData\Local\Temp\claude\C--Users-fuzquianoalricabi-simulations\060e0aa8-e6b8-455c-8565-429881efb048\scratchpad
cd /d C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1
cl /nologo /std:c++17 /EHsc /O2 /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\rockim\eigen-3.4.0 /Fo:obj_LOTB\bench_ /Fe:%SP%\bench_lotb.exe %SP%\bench_lotb.cpp
%SP%\bench_lotb.exe
