#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
t0=$(date +%s)
./rockim_yy.exe configs/indent2d_yan.cfg out_indent2d_yan > run_indent2d_yan.log 2>&1
echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> run_indent2d_yan.log
