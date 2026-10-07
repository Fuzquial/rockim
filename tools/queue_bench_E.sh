#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_bench_E.sh — bancs B1 (loi seule) et B2 (penalite seule), lances DERRIERE
# le banc D (tools/queue_bench_D.sh) : attend « banc D termine » dans
# results/queue_benches_s25.log. UN gros job a la fois. Suite au diagnostic
# independant du 12/09 (§6.3 : separer les sensibilites du banc B).
#   bash tools/queue_bench_E.sh          (lance en arriere-plan le 13/09 15 h)
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y16.exe
QLOG=results/queue_benches_s25.log
log() { echo "[queue E $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }
log "bancs B1/B2 en attente de la fin du banc D"
while ! grep -q "banc D termine" "$QLOG" 2>/dev/null; do sleep 300; done
run() {  # deck out tag
  local out="$2" name
  name="$(basename "$out" | sed 's/^out_//')"
  log "banc $3 : $EXE $1 $out"
  OMP_NUM_THREADS=14 "$EXE" "$1" "$out" > "results/$name.log" 2>&1
  log "banc $3 termine (code $?) : $(tail -1 "$out/history.csv" 2>/dev/null | cut -d, -f1) s"
}
run configs/yang2026_bench_s25_law.cfg out_yang_bench_s25_law B1
run configs/yang2026_bench_s25_pen.cfg out_yang_bench_s25_pen B2
log "file E terminee"
