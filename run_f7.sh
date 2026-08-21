#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
for S in aniso iso; do
  t0=$(date +%s)
  echo "=== F7 $S : debut $(date +%H:%M) ==="
  ./rockim_hydro.exe bench_abuaisha/configs/f7_${S}_c.cfg out_f7_${S}_c \
      > bench_abuaisha/run_f7_${S}_c.log 2>&1
  echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> bench_abuaisha/run_f7_${S}_c.log
  echo "=== F7 $S : fin, $(( ($(date +%s) - t0) / 60 )) min ==="
done
