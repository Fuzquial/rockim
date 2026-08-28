@echo off
rem ---------------------------------------------------------------------------
rem BUILD de l etat du 2026-08-28 (WP6 + les cinq reparations + reprises).
rem
rem Nom d exe DEDIE : rockim_wp6.exe. rockim_cl.exe (correctif coulomb du
rem 27/08) reste INTACT — un run en cours qui l utilise n est pas touche, et
rem la comparaison A/B reste possible.
rem
rem Contenu par rapport a rockim_cl.exe :
rem   - contactResidualMu (WP6, spec 005) : mu de contact residuel sur
rem     matiere pulverisee, opt-in, defauts bit-identiques ;
rem   - toolStop arrete VRAIMENT l outil (v = 0) ;
rem   - toolImpulseCap enfin lu en percussion 2D ;
rem   - meanTensionCapFactor garde !law_ en 3D (symetrie 2D) ;
rem   - bannieres gcRestitution / jointDeath / avertissements E3-E6 ;
rem   - garde bulkDamage durcie (toute cle law jette).
rem ---------------------------------------------------------------------------
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d %~dp0
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim_wp6.exe
