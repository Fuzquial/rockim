#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_bench_F.sh — variantes SANS GARDE-FOU (budgetAbortPct = 0) des bancs B, B1 et C,
# lancees DERRIERE la file E (tools/queue_bench_E.sh) : attend « file E terminee »
# dans results/queue_benches_s25.log. UN gros job a la fois.
# Motif (13/09, 14 h 40) : B et C se sont arretes par ENERGY ABORT a 83 et 92 us
# (exces 22 J et 4,4 J) : la loi de Solidity transcrite mot a mot cree de l energie
# en mode mixte. Ici on la laisse faire, comme leur code qui n a pas de bilan, et
# on MESURE l energie creee a cote de la reaction, du bit et des fissures.
#   bash tools/queue_bench_F.sh          (lance en arriere-plan le 13/09 14 h 45)
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y16.exe
QLOG=results/queue_benches_s25.log
log() { echo "[queue F $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }
log "bancs B'/B1'/C' (sans garde-fou) en attente de la fin de la file E"
while ! grep -q "file E terminee" "$QLOG" 2>/dev/null; do sleep 300; done
run() {  # deck out tag
  local out="$2" name
  name="$(basename "$out" | sed 's/^out_//')"
  log "banc $3 : $EXE $1 $out"
  OMP_NUM_THREADS=14 "$EXE" "$1" "$out" > "results/$name.log" 2>&1
  log "banc $3 termine (code $?) : $(tail -1 "$out/history.csv" 2>/dev/null | cut -d, -f1) s"
}
run configs/yang2026_bench_s25_solidity_noabort.cfg out_yang_bench_s25_solidity_noabort "B' (B sans garde-fou)"
run configs/yang2026_bench_s25_law_noabort.cfg out_yang_bench_s25_law_noabort "B1' (loi seule sans garde-fou)"
run configs/yang2026_bench_s25_solidity_visc_noabort.cfg out_yang_bench_s25_solidity_visc_noabort "C' (C sans garde-fou)"
log "file F terminee"
