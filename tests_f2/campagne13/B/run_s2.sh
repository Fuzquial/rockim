#!/bin/bash
# B (campagne du 13/09) : rejeu des mini-tests S2 avec rockim_g1y17.exe, OMP_NUM_THREADS = 4.
# Pas d ancre bitid ici : l ancre 9 decks tourne dans results/bitid_g1y17.log.
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/b17/s2"
D="$R/tests_f2/campagne13/S2"
EXE="$R/rockim_g1y17.exe"
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
echo "=== FIN ==="; date
