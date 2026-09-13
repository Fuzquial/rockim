#!/bin/bash
# S1 mini-tests (2e lancement : l exe dans Temp est refuse par Apex One, on
# lance la copie tests_f2/campagne13/rockim_s1test.exe = build/rockim.exe
# sha256 cc80c1ef...). OMP_NUM_THREADS=4, cwd = racine, sorties scratchpad.
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/s1b"
EXE="$R/tests_f2/campagne13/rockim_s1test.exe"
cd "$R" || exit 2
cp "$S/../s1/"*.cfg "$S/"
echo "=== bitid 2D ucs ==="; date
python tools/bitid.py --exe "$EXE" --only fdem_ucs_yan_adaptive_court --keep --outroot "$S/bitid_u" --json "$S/bitid_u.json" > "$S/bitid_u.log" 2>&1
echo "rc=$?"; tail -4 "$S/bitid_u.log"
for d in u_rupt_slipRef u_cap u_co_slipF u_co_slipRef; do
  echo "=== $d ==="; date
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"
done
echo "=== bitid 3D kuru9 ==="; date
python tools/bitid.py --exe "$EXE" --only fdem3d_kuru9_court --keep --outroot "$S/bitid_k" --json "$S/bitid_k.json" > "$S/bitid_k.log" 2>&1
echo "rc=$?"; tail -4 "$S/bitid_k.log"
for d in k9_rupt k9_cap k9_pl_slipF k9_pl_slipRef; do
  echo "=== $d ==="; date
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"
done
echo "=== FIN ==="; date
