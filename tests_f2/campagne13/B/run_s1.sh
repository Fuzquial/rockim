#!/bin/bash
# B (campagne du 13/09) : rejeu SEQUENTIEL des mini-tests S1 avec rockim_g1y17.exe, OMP_NUM_THREADS = 4.
# Les 4 decks bitid sont rejoues avec --keep --outroot pour servir de reference a check_s1.py
# (comparaison out_* contre bitid_*/<deck>), comme dans run_all.sh / run_3d2.sh de S1.
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/b17/s1"
EXE="$R/rockim_g1y17.exe"
cd "$R" || exit 2
echo "=== bitid 2D ucs ==="; date
python tools/bitid.py --exe "$EXE" --only fdem_ucs_yan_adaptive_court --keep --outroot "$S/bitid_u" --json "$S/bitid_u.json" > "$S/bitid_u.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_u.log" | tail -2
for d in u_rupt_slipRef u_cap u_co_slipF u_co_slipRef u2_ref u2_slipRef u2_co_slipF u2_co_slipRef; do
  echo "=== $d ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep "wall time" "$S/$d.log"
done
echo "=== bitid 3D kuru9 ==="; date
python tools/bitid.py --exe "$EXE" --only fdem3d_kuru9_court --keep --outroot "$S/bitid_k" --json "$S/bitid_k.json" > "$S/bitid_k.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_k.log" | tail -2
for d in k9_rupt k9_cap k9_cap2 k9_pl_slipF k9_pl_slipRef; do
  echo "=== $d ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep "wall time" "$S/$d.log"
done
echo "=== bitid heilman ==="; date
python tools/bitid.py --exe "$EXE" --only fdem3d_cut3d_heilman_court --keep --outroot "$S/bitid_h" --json "$S/bitid_h.json" > "$S/bitid_h.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_h.log" | tail -2
for d in h_pl_slipF h_pl_slipRef h_ori_slipRef; do
  echo "=== $d ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep "wall time" "$S/$d.log"
done
echo "=== bitid visc ==="; date
python tools/bitid.py --exe "$EXE" --only fdem3d_visc_yan_3d --keep --outroot "$S/bitid_v" --json "$S/bitid_v.json" > "$S/bitid_v.log" 2>&1
echo "rc=$?"; grep "IDENTIQUE\|ECHEC\|DIFF" "$S/bitid_v.log" | tail -2
for d in v_slipRef v_co_slipF v_co_slipRef; do
  echo "=== $d ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep "wall time" "$S/$d.log"
done
echo "=== FIN ==="; date
