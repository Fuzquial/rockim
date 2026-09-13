#!/bin/bash
# S2 (campagne du 13/09) : banc court de contactForcePairs, 4 fils, sequentiel.
#   1. s25_fc  : deck s = 2,5, T = 30 us, cle posee (8 paires)
#   2. s25_ref : meme deck sans la cle (temoin bit-identique attendu)
#   3. s25_bad_self / s25_bad_name : DOIVENT echouer a l initialisation
#   4. ancre bitid --refs tools/bitid_refs_jointlaw.json --only fdem3d_yang_v2_court
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/s2"
D="$R/tests_f2/campagne13/S2"
EXE="$S/rockim_s2test.exe"
cd "$R" || exit 2
for d in s25_fc s25_ref; do
  echo "=== $d ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$D/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep "wall time" "$S/$d.log"
done
for d in s25_bad_self s25_bad_name; do
  echo "=== $d (DOIT echouer) ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$D/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep -i "contactForcePairs" "$S/$d.log" | tail -2
done
echo "=== bitid fdem3d_yang_v2_court ==="; date
python tools/bitid.py --exe "$EXE" --refs tools/bitid_refs_jointlaw.json --only fdem3d_yang_v2_court --keep --outroot "$S/bitid_y" --json "$S/bitid_y.json" > "$S/bitid_y.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_y.log" | tail -3
echo "=== FIN ==="; date
