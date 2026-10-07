#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
for k in smoke_yy smoke_yy_iso; do
  echo "=== $k : debut $(date +%H:%M:%S) ==="
  ./rockim_yy.exe "configs/$k.cfg" "out_$k" > "run_$k.log" 2>&1
  echo "=== $k : fin $(date +%H:%M:%S) rc=$? ==="
done
echo "=== SMOKES TERMINES ==="
