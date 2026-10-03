#!/usr/bin/env bash
# prof_threads.sh — quelle phase du pas de temps passe mal a l'echelle ?
#
# Lance le deck court Kuru9 (reglages par defaut, pas fixe) avec ROCKIM_PROF=1
# a 1, 4 et 10 fils, lit la ligne « [prof3d] per step (ms) » (elements, insertion,
# joints, contact, outil) et en deduit le RESTE (integration, sorties, bilan :
# temps mur / pas − somme des phases). Pour chaque phase : ms par pas a chaque
# nombre de fils et acceleration par rapport a 1 fil. La phase dont
# l'acceleration plafonne est la prochaine cible.
#
# Usage, depuis la racine du depot, apres avoir compile build/rockim :
#   bash tools/prof_threads.sh
#   THREADS="1 4 6 10" DECK=mon_deck.cfg bash tools/prof_threads.sh
set -euo pipefail

THREADS=${THREADS:-"1 4 10"}
DECK=${DECK:-tests_f2/bitid/fdem3d_kuru9_court.cfg}
EXE=${EXE:-build/rockim}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/prof_threads.XXXXXX")
mkdir -p results
OUT=results/prof_threads_$(date +%Y-%m-%d_%H%M).txt

[ -x "$EXE" ] || { echo "$EXE absent : compiler d'abord"; exit 1; }
[ -f "$DECK" ] || { echo "$DECK absent : lancer depuis la racine du depot"; exit 1; }

for t in $THREADS; do
    echo "== $t fils"
    OMP_NUM_THREADS=$t ROCKIM_PROF=1 "$EXE" "$DECK" "$WORK/o$t" > "$WORK/log$t" 2>&1 \
        || { echo "ECHEC a $t fils"; tail -5 "$WORK/log$t"; exit 1; }
    grep -E "per step \(ms\)|wall time" "$WORK/log$t" | sed 's/^ *[0-9]*%//'
done

python3 - "$WORK" $THREADS <<'EOF' | tee "$OUT"
import re, sys
work, threads = sys.argv[1], [int(t) for t in sys.argv[2:]]
data = {}
for t in threads:
    log = open(f"{work}/log{t}", errors="replace").read()
    m = re.findall(r"per step \(ms\):(.*?)\((\d+) steps\)", log)
    if not m:
        sys.exit(f"pas de ligne [prof3d] a {t} fils : ROCKIM_PROF ignore par ce mode ?")
    phases = dict((k, float(v)) for k, v in re.findall(r"(\w+) ([0-9.]+)", m[-1][0]))
    steps = int(m[-1][1])
    wall = float(re.findall(r"wall time: ([0-9.]+)", log)[-1])
    phases["reste"] = max(0.0, 1e3 * wall / steps - sum(phases.values()))
    phases["TOTAL"] = 1e3 * wall / steps
    data[t] = phases
t0 = threads[0]
names = [k for k in data[t0] if data[t0][k] > 0.005 or k in ("reste", "TOTAL")]
print("\nms par pas (et acceleration par rapport a %d fil%s)" % (t0, "s" if t0 > 1 else ""))
print("%-10s" % "phase" + "".join("%20s" % ("%d fils" % t) for t in threads))
for k in names:
    row = "%-10s" % k
    for t in threads:
        v = data[t].get(k, 0.0)
        row += "%20s" % ("%.3f (x%.2f)" % (v, data[t0][k] / v if v > 0 else 0))
    print(row)
print("\nPart de chaque phase dans le pas a %d fils :" % threads[-1])
tot = data[threads[-1]]["TOTAL"]
for k in names[:-1]:
    print("  %-10s %5.1f %%" % (k, 100 * data[threads[-1]].get(k, 0) / tot))
EOF
echo "Resultat : $OUT"
rm -rf "$WORK"
