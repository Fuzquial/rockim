#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_nuit_1314b.sh — PROLONGATION du St Anne fin (rock137) SANS LE
# REDEMARRER, puis le Kuru. Remplace la surveillance de queue_nuit_1314.sh,
# dont la borne de 3 h aurait coupe le run a ~125 us, AVANT l apparition des
# radiales (sur le maillage grossier elles ne se forment qu apres 150 us).
# Demande de Fernando le 14/09 a 09 h 10 : « rajoute 2 h sans re initialiser ».
#
# PIEGE CORRIGE (09 h 27) : `kill -0 <PID>` ne voit PAS un PID Windows depuis
# Git bash — la premiere version a cru le run termine et a lance le Kuru EN
# PARALLELE. La presence du processus se teste par `tasklist`, et l arret par
# `taskkill /PID`. Meme lecon que le 12/09 sur les files de bancs.
#
#   bash tools/queue_nuit_1314b.sh <PID_du_St_Anne_en_cours>
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y19.exe
QLOG=results/queue_benches_s25.log
PID_STANNE="${1:?usage: queue_nuit_1314b.sh <PID>}"
FIN_STANNE=$(( $(date +%s) + 8400 ))          # ~2 h 20 de plus
log() { echo "[nuit $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }
vivant() { tasklist 2>/dev/null | grep -qE "rockim_g1y19.exe +$1 "; }

log "PROLONGATION du St Anne fin (PID $PID_STANNE), echeance dans 2 h 20, SANS redemarrage"
while vivant "$PID_STANNE"; do
  sleep 120
  if [ "$(date +%s)" -ge "$FIN_STANNE" ]; then
    log "St Anne fin : echeance prolongee atteinte, arret a $(tail -1 out_stanne2025_rock137/history.csv 2>/dev/null | cut -d, -f1) s"
    powershell -NoProfile -Command "Stop-Process -Id $PID_STANNE -Force -ErrorAction SilentlyContinue" > /dev/null 2>&1
    sleep 5; break
  fi
done
log "St Anne fin termine : $(tail -1 out_stanne2025_rock137/history.csv 2>/dev/null | cut -d, -f1) s, $(ls out_stanne2025_rock137 2>/dev/null | grep -c 'fdem3d_joints_') trames"

CAP=10800
log "Kuru train1 v5 (borne 3 h)"
OMP_NUM_THREADS=14 "$EXE" configs/yang2026_kuru_train1_v5.cfg out_yang_kuru_train1_v5 \
    > results/yang_kuru_train1_v5.log 2>&1 &
t0=$(date +%s)
sleep 20
KPID=$(tasklist 2>/dev/null | grep -E "^rockim_g1y19.exe" | awk '{print $2}' | tail -1)
log "Kuru lance (PID $KPID)"
while vivant "$KPID"; do
  sleep 120
  if [ $(( $(date +%s) - t0 )) -ge "$CAP" ]; then
    log "Kuru : borne de 3 h atteinte, arret a $(tail -1 out_yang_kuru_train1_v5/history.csv 2>/dev/null | cut -d, -f1) s"
    powershell -NoProfile -Command "Stop-Process -Id $KPID -Force -ErrorAction SilentlyContinue" > /dev/null 2>&1
    sleep 5; break
  fi
done
log "Kuru train1 v5 termine : $(tail -1 out_yang_kuru_train1_v5/history.csv 2>/dev/null | cut -d, -f1) s, $(ls out_yang_kuru_train1_v5 2>/dev/null | grep -c 'fdem3d_joints_') trames"
log "file de la nuit terminee"
