#!/usr/bin/env bash
# Lance les simulations du rapport-guide, un calcul a la fois, depuis la racine du depot.
# Usage : bash docs/rapport_guide/simulations_a_lancer/lancer.sh <groupe> [fils]
#   groupe : yan | calib | diag_contact | abuaisha | tunnel | impact | tout
#   fils   : OMP_NUM_THREADS (defaut 1 ; 1 = sorties identiques au bit au calcul serie)
# Sorties : sim_out/<groupe>/<nom>/ (hors git) et journal sim_out/<groupe>/<nom>.log
set -u
cd "$(git rev-parse --show-toplevel)"
EXE=./build/rockim
[ -x "$EXE" ] || { echo "compiler d'abord : mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j"; exit 1; }
GRP=${1:-tout}; NT=${2:-1}
DIR=docs/rapport_guide/simulations_a_lancer/configs
if [ "$GRP" = tout ]; then GROUPES="yan calib diag_contact abuaisha tunnel impact"; else GROUPES="$GRP"; fi
for g in $GROUPES; do
  mkdir -p sim_out/$g
  for c in $DIR/$g/*.cfg; do
    n=$(basename "$c" .cfg)
    if grep -q "wall time" "sim_out/$g/$n.log" 2>/dev/null; then echo "deja fait : $g/$n"; continue; fi
    echo "=== $g/$n  debut $(date '+%F %T')"
    OMP_NUM_THREADS=$NT nice -n 10 "$EXE" "$c" "sim_out/$g/$n" > "sim_out/$g/$n.log" 2>&1
    echo "    fin $(date '+%F %T')  code $?"
  done
done
