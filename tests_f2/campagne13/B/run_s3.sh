#!/bin/bash
# B (campagne du 13/09) : rejeu des mini-tests S3 avec rockim_g1y17.exe, OMP_NUM_THREADS = 4.
# Pas d ancre bitid ici : l ancre complete 8/8 + 1/1 tourne dans results/bitid_g1y17.log.
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/b17/s3"
D="$R/tests_f2/campagne13/S3"
EXE="$R/rockim_g1y17.exe"
cd "$R" || exit 2
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
echo "=== FIN ==="; date
