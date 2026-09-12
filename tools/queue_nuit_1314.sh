#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_nuit_1314.sh — le KURU de la nuit du 13 au 14/09, DERRIERE St Anne.
# Demande de Fernando : « ensuite lance st anne 2025, et ensuite Kuru mais
# essaye de faire un run de 2-3 h max ». UN gros job a la fois.
#
# Le budget est tenu par une BORNE DE TEMPS MUR (3 h) et non par un T reduit :
# les trames et history.csv sont ecrits au fil de l eau, donc l arret a 3 h
# garde tout ce qui a ete calcule. dt mesure a la fumee : 1,671 ns (contre
# 1,0 ns du run v3 s = 1 : c est jointPenaltyLength = edge qui le double),
# 119 665 pas pour les 200 us du deck.
#   bash tools/queue_nuit_1314.sh        (lance en arriere-plan le 14/09 01 h 35)
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y18.exe
QLOG=results/queue_benches_s25.log
CAP=10800                     # 3 h de temps mur
log() { echo "[nuit $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }
log "Kuru train1 v5 en attente de la fin de St Anne"
while tasklist 2>/dev/null | grep -qi "rockim_g1y18.exe"; do sleep 300; done
log "machine libre : Kuru train1 v5 (borne de temps mur $((CAP/3600)) h)"
OMP_NUM_THREADS=14 "$EXE" configs/yang2026_kuru_train1_v5.cfg out_yang_kuru_train1_v5 \
    > results/yang_kuru_train1_v5.log 2>&1 &
PID=$!
t0=$(date +%s)
while kill -0 "$PID" 2>/dev/null; do
  sleep 120
  now=$(date +%s)
  if [ $((now - t0)) -ge "$CAP" ]; then
    log "borne de 3 h atteinte : arret propre du Kuru a $(tail -1 out_yang_kuru_train1_v5/history.csv 2>/dev/null | cut -d, -f1) s"
    kill "$PID" 2>/dev/null
    sleep 5
    break
  fi
done
wait "$PID" 2>/dev/null
log "Kuru train1 v5 termine : $(tail -1 out_yang_kuru_train1_v5/history.csv 2>/dev/null | cut -d, -f1) s, $(ls out_yang_kuru_train1_v5 2>/dev/null | grep -c 'fdem3d_joints_') trames"
log "file de la nuit terminee"
