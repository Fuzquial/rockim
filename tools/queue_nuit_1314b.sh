#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_nuit_1314b.sh — PROLONGATION du St Anne fin (rock137) SANS LE
# REDEMARRER, puis le Kuru. Remplace la surveillance de queue_nuit_1314.sh,
# dont la borne de 3 h aurait coupe le run a ~125 us, c est-a-dire AVANT que
# les radiales n apparaissent (sur le maillage grossier elles ne se forment
# qu apres 150 us et ne se detachent qu entre 200 et 300 us).
#
# Demande de Fernando le 14/09 a 09 h 10 : « rajoute 2 h dans ce cas oui sans
# re initialiser ». Le processus rockim_g1y19 lance a 06 h 34 CONTINUE : seule
# la surveillance a ete remplacee (l ancienne file a ete arretee, le fils a
# survecu). Nouvelle echeance = 06 h 34 + 5 h = 11 h 34, soit ~210 us au
# rythme mesure de 0,705 us/min.
#
#   bash tools/queue_nuit_1314b.sh       (lance le 14/09 09 h 12)
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y19.exe
QLOG=results/queue_benches_s25.log
PID_STANNE="${1:-16880}"      # PID du St Anne fin deja en cours
FIN_STANNE=$(( $(date +%s) + 8700 ))   # ~2 h 25 de plus a partir de maintenant
log() { echo "[nuit $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }

log "PROLONGATION du St Anne fin (PID $PID_STANNE) : nouvelle echeance dans 2 h 25, sans redemarrage"
while kill -0 "$PID_STANNE" 2>/dev/null; do
  sleep 120
  if [ "$(date +%s)" -ge "$FIN_STANNE" ]; then
    log "St Anne fin : echeance prolongee atteinte, arret propre a $(tail -1 out_stanne2025_rock137/history.csv 2>/dev/null | cut -d, -f1) s"
    kill "$PID_STANNE" 2>/dev/null; sleep 5; break
  fi
done
log "St Anne fin termine : $(tail -1 out_stanne2025_rock137/history.csv 2>/dev/null | cut -d, -f1) s, $(ls out_stanne2025_rock137 2>/dev/null | grep -c 'fdem3d_joints_') trames"

CAP=10800
log "Kuru train1 v5 (borne 3 h) : $EXE configs/yang2026_kuru_train1_v5.cfg out_yang_kuru_train1_v5"
OMP_NUM_THREADS=14 "$EXE" configs/yang2026_kuru_train1_v5.cfg out_yang_kuru_train1_v5 \
    > results/yang_kuru_train1_v5.log 2>&1 &
PID=$!
t0=$(date +%s)
while kill -0 "$PID" 2>/dev/null; do
  sleep 120
  if [ $(( $(date +%s) - t0 )) -ge "$CAP" ]; then
    log "Kuru : borne de 3 h atteinte, arret propre a $(tail -1 out_yang_kuru_train1_v5/history.csv 2>/dev/null | cut -d, -f1) s"
    kill "$PID" 2>/dev/null; sleep 5; break
  fi
done
wait "$PID" 2>/dev/null
log "Kuru train1 v5 termine : $(tail -1 out_yang_kuru_train1_v5/history.csv 2>/dev/null | cut -d, -f1) s, $(ls out_yang_kuru_train1_v5 2>/dev/null | grep -c 'fdem3d_joints_') trames"
log "file de la nuit terminee"
