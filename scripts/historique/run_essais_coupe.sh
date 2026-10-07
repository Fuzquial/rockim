#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
EXE=./rockim_pot.exe
for e in pot001 potxi dt005 slow; do
  echo "=== $e : debut $(date +%H:%M:%S) ==="
  $EXE "tunnel_edz/configs/cut2d_$e.cfg" "out_cut_$e" > "tunnel_edz/run_cut_$e.log" 2>&1
  echo "=== $e : fin $(date +%H:%M:%S) rc=$? ==="
done
echo "=== CHAINE TERMINEE ==="
