#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
t0=$(date +%s)
./rockim_yy.exe configs/indent3d_grad.cfg out_indent3d_grad > run_indent3d_grad.log 2>&1
echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> run_indent3d_grad.log
