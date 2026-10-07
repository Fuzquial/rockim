#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_nuit_1314c.sh — SANS BORNE DE TEMPS pour le St Anne fin (rock137).
# Demande de Fernando le 14/09 a 10 h 45 : « au pire enleve la limite et on
# verra quand on aura ». Le run va donc jusqu'au bout de son T = 300 us, puis
# le Kuru s'enchaine (lui, borne a 3 h : c'est la demande initiale).
#
# Le processus rockim_g1y19 (St Anne fin) lance a 06 h 34 CONTINUE : seule la
# surveillance est remplacee, comme a 09 h 28. Rythme mesure a 10 h 41 :
# 0,54 us/min (il ralentit avec le nombre de facettes rompues) -> les 166 us
# restants demandent ~5 h, fin vers 16 h.
#
# PIEGES DEJA PAYES (ne pas les reintroduire) :
#   - `kill -0 <PID>` ne voit pas un PID Windows depuis Git bash : la boucle
#     sort aussitot et le job suivant demarre EN PARALLELE (14/09, 09 h 27).
#   - `tasklist /FI ...` : Git bash transforme /FI en chemin.
#   - ne PAS tester la fin par une ligne du journal : les files defectueuses
#     du 14/09 y ont laisse deux « St Anne fin termine » perimees.
#   Test de presence : `tasklist | grep -E "exe +PID "`. Arret : Stop-Process.
#
#   bash tools/queue_nuit_1314c.sh <PID_du_St_Anne_en_cours>
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y19.exe
QLOG=results/queue_benches_s25.log
PID_STANNE="${1:?usage: queue_nuit_1314c.sh <PID>}"
log() { echo "[nuit $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }
vivant() { tasklist 2>/dev/null | grep -qE "rockim_g1y19.exe +$1 "; }

log "St Anne fin (PID $PID_STANNE) : BORNE DE TEMPS RETIREE, il ira au bout des 300 us"
while vivant "$PID_STANNE"; do sleep 300; done
log "St Anne fin acheve de lui-meme : $(tail -1 out_stanne2025_rock137/history.csv 2>/dev/null | cut -d, -f1) s, $(ls out_stanne2025_rock137 2>/dev/null | grep -c 'fdem3d_joints_') trames"

CAP=10800
log "Kuru train1 v5 (borne 3 h, demande initiale de Fernando)"
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
