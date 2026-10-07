#!/bin/bash
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"
# attendre la fin du 2D
while tasklist 2>/dev/null | grep -qi rockim_yy; do sleep 20; done
echo "=== 2D fini a $(date +%H:%M:%S), lancement du 3D ==="
t0=$(date +%s)
./rockim_yy.exe configs/indent3d_yan.cfg out_indent3d_yan > run_indent3d_yan.log 2>&1
echo "=== 3D TERMINE en $(( $(date +%s) - t0 )) s ==="
