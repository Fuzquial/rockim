#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
t0=$(date +%s)
./rockim_yy.exe configs/yang_equiv.cfg out_yang_equiv > run_yang_equiv.log 2>&1
echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> run_yang_equiv.log
