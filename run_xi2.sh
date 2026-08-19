#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
t0=$(date +%s)
./rockim_yy.exe configs/impact3d_yang_court_xi2.cfg out_yang_court_xi2 > run_yang_court_xi2.log 2>&1
echo "=== xi=2 : $(( $(date +%s) - t0 )) s, rc=$? ===" >> run_yang_court_xi2.log
