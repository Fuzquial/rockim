#!/bin/bash
# S4 (campagne du 13/09) : essais elementaires de la pulverisation, OMP_NUM_THREADS = 4.
# 1) tous les decks de make_decks.py avec l exe S4 ; 2) les temoins *_ref avec rockim_g1y16.exe
# (bit-identite du chemin par defaut sur ces decks) ; 3) ancre bitid fdem3d_kuru9_court (bulkDamage = yang).
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/s4"
EXE="$R/tests_f2/campagne13/S4/rockim_s4test.exe"
OLD="$R/rockim_g1y16.exe"
cd "$R" || exit 2
for d in iso_ref iso_dev iso_tot iso_pri iso_edge iso_edge_tot \
         biax_ref biax_dev biax_tot biax_pri biax_edge biax_edge_tot \
         uni_ref uni_dev uni_tot uni_pri uni_edge uni_edge_tot \
         uni2d_ref uni2d_dev uni2d_tot uni2d_pri uni2d_edge \
         bad_nokey bad_value bad_len; do
  echo "=== $d ==="; date
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"
  tail -1 "$S/$d.log"
done
for d in iso_ref biax_ref uni_ref uni2d_ref; do
  echo "=== $d (rockim_g1y16.exe) ==="; date
  "$OLD" "$S/$d.cfg" "$S/old_$d" > "$S/old_$d.log" 2>&1
  echo "rc=$?"
done
echo "=== bitid fdem3d_kuru9_court ==="; date
python tools/bitid.py --exe "$EXE" --only fdem3d_kuru9_court --json "$S/bitid_kuru9.json" > "$S/bitid_kuru9.log" 2>&1
echo "rc=$?"
grep -i "identique\|differ\|\[bitid\]" "$S/bitid_kuru9.log"
echo "=== FIN ==="; date
