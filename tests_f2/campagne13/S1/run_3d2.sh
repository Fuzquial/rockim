#!/bin/bash
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/s1b"
EXE="$R/tests_f2/campagne13/rockim_s1test.exe"
cd "$R" || exit 2
until grep -q "=== FIN ===" "$S/run_all.log"; do sleep 10; done
echo "=== bitid heilman ==="; date
python tools/bitid.py --exe "$EXE" --only fdem3d_cut3d_heilman_court --keep --outroot "$S/bitid_h" --json "$S/bitid_h.json" > "$S/bitid_h.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_h.log" | tail -2
for d in h_pl_slipF h_pl_slipRef h_ori_slipRef; do
  echo "=== $d ==="; date
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"
done
echo "=== bitid visc ==="; date
python tools/bitid.py --exe "$EXE" --only fdem3d_visc_yan_3d --keep --outroot "$S/bitid_v" --json "$S/bitid_v.json" > "$S/bitid_v.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_v.log" | tail -2
for d in v_slipRef v_co_slipF v_co_slipRef k9_cap2; do
  echo "=== $d ==="; date
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"
done
echo "=== FIN ==="; date
