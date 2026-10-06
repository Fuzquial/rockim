#!/bin/bash
# usage: lance.sh cfg out
cd /tmp/claude-0/-home-user/95ba17a5-6321-5243-8157-c0b3e0616268/scratchpad/runs/coupe
export OMP_NUM_THREADS=2
t0=$(date +%s)
echo "debut $(date -Is) : nice -n 10 timeout 1800 ../../build/rockim $1 $2" > $2.time
nice -n 10 timeout 1800 ../../build/rockim $1 $2 > $2.log 2>&1
rc=$?
t1=$(date +%s)
echo "fin $(date -Is) rc=$rc duree=$((t1-t0)) s" >> $2.time
