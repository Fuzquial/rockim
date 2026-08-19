#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
t0=$(date +%s)
./rockim_yy.exe configs/impact_uni.cfg out_impact_uni > run_impact_uni.log 2>&1
echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> run_impact_uni.log
