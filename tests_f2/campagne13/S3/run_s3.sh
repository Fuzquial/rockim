#!/bin/bash
# S3 (campagne du 13/09) : banc de joint cinematique, 4 fils, sequentiel.
#   1. neuf decks jointbench (secondes chacun) -> logs + jointbench.csv dans le scratchpad
#   2. deux variantes qui DOIVENT etre refusees (cle jb* hors scenario ; mesh pose)
#   3. ancres bitid : fdem3d_kuru9_court (refs par defaut) et fdem3d_yang_v2_court (refs jointlaw)
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/s3"
D="$R/tests_f2/campagne13/S3"
EXE="$S/rockim_s3test.exe"
cd "$R" || exit 2
mkdir -p "$S"
for d in ten_solidity ten_plastic ten_origin sh_solidity sh_plastic sh_origin tilt_sol tilt_pl_maj tilt_pl_any; do
  echo "=== $d ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$D/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep "wall time" "$S/$d.log"
done
for d in bad_key_percussion bad_mesh; do
  echo "=== $d (DOIT echouer) ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$D/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep -i "error" "$S/$d.log" | tail -1
done
echo "=== bitid fdem3d_kuru9_court (refs defaut) ==="; date
python tools/bitid.py --exe "$EXE" --only fdem3d_kuru9_court --outroot "$S/bitid_k" --json "$S/bitid_k.json" > "$S/bitid_k.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_k.log" | tail -3
echo "=== bitid fdem3d_yang_v2_court (refs jointlaw) ==="; date
python tools/bitid.py --exe "$EXE" --refs tools/bitid_refs_jointlaw.json --only fdem3d_yang_v2_court --outroot "$S/bitid_y" --json "$S/bitid_y.json" > "$S/bitid_y.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_y.log" | tail -3
echo "=== FIN ==="; date
