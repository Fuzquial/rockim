#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_benches_s25.sh — file d'attente des trois bancs s = 2,5 du 13/09
# (A temoin v3-P 300 us, B loi de Solidity + penalite edge/25, C = B + visc 4000)
# DERRIERE le run s = 1 (rockim_g1y15.exe) : UN gros job a la fois (regle de
# Fernando). Chaque banc part quand (1) le run s = 1 est termine (processus
# absent) et (2) l'ancre de bit-identite de rockim_g1y16 est ecrite
# (results/bitid_g1y16.done) ; B et C ne partent que si l'ancre est IDENTIQUE.
#   bash tools/queue_benches_s25.sh          (lance en arriere-plan le 13/09 14 h)
# Journal : results/queue_benches_s25.log ; journaux des bancs : results/<out>.log
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y16.exe
QLOG=results/queue_benches_s25.log
log() { echo "[queue $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }
log "en attente : fin du run s = 1 (rockim_g1y15.exe) et ancre bit-identite (results/bitid_g1y16.done)"
while tasklist 2>/dev/null | grep -qi "rockim_g1y15.exe"; do sleep 300; done
log "run s = 1 termine"
while [ ! -f results/bitid_g1y16.done ]; do sleep 120; done
log "ancre ecrite : $(grep -c IDENTIQUE results/bitid_g1y16.log) ligne(s) IDENTIQUE"
run() {  # deck out tag
  local out="$2" name
  name="$(basename "$out" | sed 's/^out_//')"
  log "banc $3 : $EXE $1 $out"
  OMP_NUM_THREADS=14 "$EXE" "$1" "$out" > "results/$name.log" 2>&1
  log "banc $3 termine (code $?) : $(tail -1 "$out/history.csv" 2>/dev/null | cut -d, -f1) s"
}
run configs/yang2026_bench_s25_v3P_300.cfg out_yang_bench_s25_v3P_300 A
if grep -q "8/8 IDENTIQUE" results/bitid_g1y16.log && grep -q "1/1 IDENTIQUE" results/bitid_g1y16.log; then
  run configs/yang2026_bench_s25_solidity.cfg out_yang_bench_s25_solidity B
  run configs/yang2026_bench_s25_solidity_visc.cfg out_yang_bench_s25_solidity_visc C
else
  log "ancre bit-identite NON identique (voir results/bitid_g1y16.log) : bancs B et C NON lances"
fi
log "file terminee"
