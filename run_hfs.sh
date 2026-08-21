#!/bin/bash
# B2 — pression de rupture au forage, SIGNE CORRIGE (2026-08-20).
# Un job a la fois : aniso puis iso, en sequence.
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
for S in aniso iso; do
  t0=$(date +%s)
  echo "=== $S : debut $(date +%H:%M:%S) ===" >> bench_abuaisha/run_hfs.log
  ./rockim_hy.exe bench_abuaisha/configs/hf_${S}_hydro_s.cfg out_hfs_${S} \
      > bench_abuaisha/run_hfs_${S}.log 2>&1
  echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> bench_abuaisha/run_hfs_${S}.log
  echo "=== $S fini en $(( ($(date +%s) - t0) / 60 )) min ===" >> bench_abuaisha/run_hfs.log
done
echo "=== CHAINE TERMINEE $(date +%H:%M:%S) ===" >> bench_abuaisha/run_hfs.log
