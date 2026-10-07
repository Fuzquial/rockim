#!/bin/bash
# S2bis (campagne du 13/09) : le canal insert/roche de contactForcePairs.
#   s25_fcrock : train lance a -20 m/s, T = 4 us, trackGroup = insert
#     -> Fc_insert_rock DOIT devenir non nulle (C1 du banc 30 us dit 0 exact :
#        la variante falsifiante, c est le banc 30 us lui-meme)
# 4 fils, exe = rockim_g1y17.exe (binaire de la campagne).
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/s2bis"
D="$R/tests_f2/campagne13/S2"
EXE="$R/rockim_g1y17.exe"
mkdir -p "$S"
cd "$R" || exit 2
for d in s25_fcrock; do
  echo "=== $d ==="; date
  rm -rf "$S/out_$d"
  "$EXE" "$D/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"; grep "wall time" "$S/$d.log"
done
echo "=== FIN ==="; date
