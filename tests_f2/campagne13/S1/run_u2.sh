#!/bin/bash
export OMP_NUM_THREADS=4
R=C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_g1
S="C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/s1b"
EXE="$R/tests_f2/campagne13/rockim_s1test.exe"
cd "$R" || exit 2
for d in u2_ref u2_slipRef u2_co_slipF u2_co_slipRef; do
  echo "=== $d ==="; date
  "$EXE" "$S/$d.cfg" "$S/out_$d" > "$S/$d.log" 2>&1
  echo "rc=$?"
done
echo "=== FIN ==="; date
