#!/bin/bash
# IMPACT St Anne 10,66 m/s — spec 005, variante s = 1,5 (nuit du 2026-08-22)
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
t0=$(date +%s)
echo "=== impact stanne s15 : debut $(date +%H:%M:%S) ===" > bench_impact/run_stanne_s15.log
./rockim_i3.exe bench_impact/configs/impact_stanne_s15.cfg out_imp_stanne \
    >> bench_impact/run_stanne_s15.log 2>&1
echo "=== TERMINE (rc=$?) en $(( $(date +%s) - t0 )) s ===" >> bench_impact/run_stanne_s15.log
