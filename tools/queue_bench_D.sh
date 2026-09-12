#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_bench_D.sh — banc D (temoin v3-P + le seul contact de Solidity), lance
# DERRIERE la file A/B/C (tools/queue_benches_s25.sh) : attend « file terminee »
# dans results/queue_benches_s25.log. UN gros job a la fois.
#   bash tools/queue_bench_D.sh          (lance en arriere-plan le 13/09 14 h 30)
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y16.exe
QLOG=results/queue_benches_s25.log
log() { echo "[queue D $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }
log "banc D en attente de la fin de la file A/B/C"
while ! grep -q "file terminee" "$QLOG" 2>/dev/null; do sleep 300; done
OUT=out_yang_bench_s25_v3P_contact
log "banc D : $EXE configs/yang2026_bench_s25_v3P_contact.cfg $OUT"
OMP_NUM_THREADS=14 "$EXE" configs/yang2026_bench_s25_v3P_contact.cfg "$OUT" > results/yang_bench_s25_v3P_contact.log 2>&1
log "banc D termine (code $?) : $(tail -1 "$OUT/history.csv" 2>/dev/null | cut -d, -f1) s"
