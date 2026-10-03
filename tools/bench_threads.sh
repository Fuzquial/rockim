#!/usr/bin/env bash
# bench_threads.sh — gain de la campagne d'optimisation selon le nombre de fils.
#
# Compare, sur le deck court Kuru9 (tests_f2/bitid/fdem3d_kuru9_court.cfg,
# maillage versionne), quatre variantes a 1, 2, 4 et 10 fils :
#   avant        : code d'avant la campagne (commit aac12f8, meme fonctionnalites)
#   apres        : code actuel, reglages par defaut (resultats bit-identiques a avant)
#   apres_dt     : + dtUpdate = inserted            (opt-in, Wu et al. 2024)
#   apres_dt_vol : + potForce = volume, facteur 5   (opt-in, Liu et al. 2022)
#
# Usage, depuis la racine du depot, apres avoir compile build/rockim :
#   bash tools/bench_threads.sh                 # 1 2 4 10 fils
#   THREADS="4 10" REPS=2 bash tools/bench_threads.sh
# Sortie : results/bench_threads_<date>.csv et un tableau a l'ecran.
# Le binaire « avant » est compile une fois dans ../rockim_avant (git worktree).
set -euo pipefail

THREADS=${THREADS:-"1 2 4 10"}
REPS=${REPS:-1}
REF=aac12f8
DECK=tests_f2/bitid/fdem3d_kuru9_court.cfg
ROOT=$(pwd)
AVANT=$ROOT/../rockim_avant
WORK=$(mktemp -d "${TMPDIR:-/tmp}/bench_threads.XXXXXX")
mkdir -p results
CSV=results/bench_threads_$(date +%Y-%m-%d_%H%M).csv

[ -x build/rockim ] || { echo "build/rockim absent : compiler d'abord"; exit 1; }
[ -f "$DECK" ] || { echo "$DECK absent : lancer depuis la racine du depot"; exit 1; }

# Memes drapeaux OpenMP que le build courant (libomp de Homebrew sur macOS).
CMAKE_OMP=()
if [ "$(uname)" = Darwin ]; then
    OMP=$(brew --prefix libomp)
    CMAKE_OMP=(-DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$OMP/include"
               -DOpenMP_CXX_LIB_NAMES=omp
               -DOpenMP_omp_LIBRARY="$OMP/lib/libomp.dylib")
    otool -L build/rockim | grep -q omp \
        || echo "ATTENTION : build/rockim n'est pas lie a libomp (build serie) : les fils n'auront aucun effet"
fi

if [ ! -x "$AVANT/build/rockim" ]; then
    echo "== compilation de la version d'avant ($REF) dans $AVANT"
    [ -d "$AVANT" ] || git worktree add --detach "$AVANT" "$REF"
    cmake -S "$AVANT" -B "$AVANT/build" -DCMAKE_BUILD_TYPE=Release "${CMAKE_OMP[@]}" >/dev/null 2>&1 \
        && cmake --build "$AVANT/build" -j >/dev/null 2>&1 \
        || { echo "echec de compilation de $REF dans $AVANT"; exit 1; }
fi

cp "$DECK" "$WORK/base.cfg"
{ cat "$DECK"; echo "dtUpdate = inserted"; } > "$WORK/dt.cfg"
{ cat "$DECK"; echo "dtUpdate = inserted"; echo "potForce = volume";
  echo "potVolumeFactor = 5"; } > "$WORK/dt_vol.cfg"

run() {  # run <exe> <cfg> <fils> -> temps mur en s
    local out="$WORK/out_$RANDOM"
    OMP_NUM_THREADS=$3 "$1" "$2" "$out" > "$out.log" 2>&1 \
        || { echo "ECHEC" ; tail -3 "$out.log" >&2; return; }
    sed -n 's/.*wall time: \([0-9.]*\) s.*/\1/p' "$out.log" | tail -1
    rm -rf "$out"
}

echo "fils,rep,avant,apres,apres_dt,apres_dt_vol" > "$CSV"
for t in $THREADS; do
    for r in $(seq 1 "$REPS"); do
        a=$(run "$AVANT/build/rockim" "$WORK/base.cfg" "$t")
        b=$(run "$ROOT/build/rockim" "$WORK/base.cfg" "$t")
        c=$(run "$ROOT/build/rockim" "$WORK/dt.cfg" "$t")
        d=$(run "$ROOT/build/rockim" "$WORK/dt_vol.cfg" "$t")
        echo "$t,$r,$a,$b,$c,$d" >> "$CSV"
        echo "fils=$t rep=$r : avant $a s | apres $b s | +dt $c s | +dt+vol $d s"
    done
done

python3 - "$CSV" <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
cols = ["avant", "apres", "apres_dt", "apres_dt_vol"]
best = {}
for r in rows:
    for c in cols:
        try: v = float(r[c])
        except ValueError: continue
        k = (int(r["fils"]), c)
        best[k] = min(best.get(k, v), v)
ts = sorted({k[0] for k in best})
print("\nTemps mur (s, meilleur des repetitions) et gain par rapport a « avant » au meme nombre de fils")
print("%5s %9s %16s %16s %16s" % ("fils", "avant", "apres", "+dtUpdate", "+dt+volume"))
for t in ts:
    a = best.get((t, "avant"))
    cells = []
    for c in cols[1:]:
        v = best.get((t, c))
        cells.append("%7.1f (x%.2f)" % (v, a / v) if v and a else "%16s" % "-")
    print("%5d %9s %s" % (t, "%.1f" % a if a else "-", " ".join(cells)))
a1 = best.get((1, "avant"))
if a1:
    print("\nAcceleration par les fils (meme version, par rapport a 1 fil) :")
    for c in cols:
        s1 = best.get((1, c))
        if s1:
            print("  %-13s" % c + "  ".join("%d fils x%.2f" % (t, s1 / best[(t, c)])
                                          for t in ts if (t, c) in best))
EOF
echo "CSV : $CSV"
rm -rf "$WORK"
