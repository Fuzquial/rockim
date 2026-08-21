#!/bin/bash
# ESSAI 1 — isotrope a cible analytique egale (12,00 MPa)
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
t0=$(date +%s)
echo "=== essai 1 : debut $(date +%H:%M:%S) ===" > bench_abuaisha/run_e1.log
./rockim_e1.exe bench_abuaisha/configs/e1_iso_cible12.cfg out_e1_iso12 \
    >> bench_abuaisha/run_e1.log 2>&1
echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> bench_abuaisha/run_e1.log
