#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
for S in aniso iso; do
  t0=$(date +%s)
  echo "=== HYDRO PRODUCTION $S : debut $(date +%H:%M) ==="
  ./rockim_hydro.exe bench_abuaisha/configs/hf_${S}_hydro.cfg out_hfp_${S} \
      > bench_abuaisha/run_hfp_${S}.log 2>&1
  echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> bench_abuaisha/run_hfp_${S}.log
  echo "=== $S fini en $(( ($(date +%s) - t0) / 60 )) min ==="
done
