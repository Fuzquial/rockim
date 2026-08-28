@echo off
rem ---------------------------------------------------------------------------
rem CAMPAGNE DE COUPLAGE (2026-08-28, WP6) — la physique complete du papier
rem 2026 pour la premiere fois, avec son controle.
rem
rem   run_pulv_campagne.cmd            lance P1 puis CTRL, l un apres l autre
rem   run_pulv_campagne.cmd p1         lance seulement P1
rem   run_pulv_campagne.cmd ctrl       lance seulement le controle
rem
rem ECHELLE A TROIS BARREAUX (le run A du 22/08 sert de troisieme point) :
rem   impact_pulv_a        (deja fait)  bulkDamage seul
rem   impact_pulv_coulomb_ctrl          + jointShearRange = coulomb
rem   impact_pulv_coulomb               + contactResidualMu = 0.18
rem CTRL vs A isole le correctif coulomb ; P1 vs CTRL isole WP6.
rem
rem BINAIRE : rockim_wp6.exe (build_wp6.cmd). rockim_cl.exe reste intact.
rem
rem A VERIFIER AU LOG, sans quoi le run ne teste pas ce qu on croit :
rem   [FDEM3D] jointShearRange = coulomb
rem   [FDEM3D] contactResidualMu = 0.18            (P1 seulement)
rem   [FDEM3D] contact residuel : N evaluations    (fin de run, N > 0)
rem Une cle mal orthographiee est IGNOREE EN SILENCE par le lecteur de config.
rem ---------------------------------------------------------------------------
cd /d %~dp0

if not exist rockim_wp6.exe (
  echo ERREUR : rockim_wp6.exe absent. Lance d abord build_wp6.cmd
  exit /b 1
)
if not exist meshes\impact_pulv.msh (
  echo ERREUR : meshes\impact_pulv.msh absent.
  echo Il est genere localement, pas versionne. Regenere-le avec
  echo   python tools\make_impact_mesh.py meshes\impact_pulv.msh ^<echelle^> 2e-4 ^<echelle_roche^> leger
  echo et dis-moi les valeurs employees le 22/08 si tu les as.
  exit /b 1
)

set CIBLE=%1
if "%CIBLE%"=="" set CIBLE=tout

if "%CIBLE%"=="p1"   goto P1
if "%CIBLE%"=="ctrl" goto CTRL

:P1
echo === P1 : physique complete (coulomb + bulkDamage + mu residuel) ===
rockim_wp6.exe bench_impact\configs\impact_pulv_coulomb.cfg out_pulv_coulomb
if "%CIBLE%"=="p1" goto FIN

:CTRL
echo === CTRL : identique SANS contactResidualMu ===
rockim_wp6.exe bench_impact\configs\impact_pulv_coulomb_ctrl.cfg out_pulv_coulomb_ctrl

:FIN
echo === termine. Envoie a Claude les history.csv et les logs des deux runs ===
