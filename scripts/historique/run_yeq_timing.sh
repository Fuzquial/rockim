#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
t0=$(date +%s)
./rockim_yy.exe configs/yang_equiv_smoke.cfg out_yeq_smoke > run_yeq_smoke.log 2>&1
echo "=== 18522 pas en $(( $(date +%s) - t0 )) s ===" >> run_yeq_smoke.log
