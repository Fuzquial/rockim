#!/bin/bash
# B (campagne du 13/09) : rejeu des mini-tests S4 avec rockim_g1y17.exe, OMP_NUM_THREADS = 4.
# Temoins *_ref rejoues aussi avec rockim_g1y16.exe (bit-identite du chemin par defaut).
# Pas d ancre bitid ici : l ancre complete tourne dans results/bitid_g1y17.log.
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/b17/s4"
EXE="$R/rockim_g1y17.exe"
OLD="$R/rockim_g1y16.exe"
cd "$R" || exit 2
for d in iso_ref iso_dev iso_tot iso_pri iso_edge iso_edge_tot \
         biax_ref biax_dev biax_tot biax_pri biax_edge biax_edge_tot \
         uni_ref uni_dev uni_tot uni_pri uni_edge uni_edge_tot \
         uni2d_ref uni2d_dev uni2d_tot uni2d_pri uni2d_edge \
         bad_nokey bad_value bad_len; do
  echo "=== $d ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"
  tail -1 "$S/$d.log"
done
for d in iso_ref biax_ref uni_ref uni2d_ref; do
  echo "=== $d (rockim_g1y16.exe) ==="; date
  rm -rf "$S/old_$d"
  "$OLD" "$S/$d.cfg" "$S/old_$d" > "$S/old_$d.log" 2>&1
  echo "rc=$?"
done
echo "=== FIN ==="; date
