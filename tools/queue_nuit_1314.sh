#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# queue_nuit_1314.sh — la nuit du 13 au 14/09, TROIS runs a la suite, un seul
# gros job a la fois, chacun borne en TEMPS MUR (les trames et history.csv sont
# ecrits au fil de l eau : l arret a la borne garde tout ce qui est calcule).
#
#   1. (deja lance) St Anne sur `rock25`, 450 us — la MECANIQUE : reaction,
#      penetration, retournement, rebond. Coeur de 839 tetraedres a 3,59 mm :
#      trop grossier pour juger la fissuration (relecture du 14/09).
#   2. St Anne sur `rock137`, 300 us, borne 3 h — la LOCALISATION. Seule
#      variable contre le run 1 : la resolution de la roche (14 030 tetraedres
#      dans le coeur, arete mediane 1,42 mm) ; train FIGE identique.
#   3. Kuru `train1_v5` sur `rock137`, 200 us, borne 3 h — la demande de
#      Fernando (« un run de 2-3 h max »).
#
#   bash tools/queue_nuit_1314.sh        (relance le 14/09 02 h 25)
# ---------------------------------------------------------------------------
cd "$(dirname "$0")/.." || exit 1
EXE=./rockim_g1y19.exe
QLOG=results/queue_benches_s25.log
log() { echo "[nuit $(date '+%d/%m %H:%M')] $*" | tee -a "$QLOG"; }

run_capped() {   # deck out journal borne_s etiquette
  local deck="$1" out="$2" jr="$3" cap="$4" tag="$5"
  log "$tag : $EXE $deck $out (borne $((cap / 3600)) h)"
  OMP_NUM_THREADS=14 "$EXE" "$deck" "$out" > "results/$jr.log" 2>&1 &
  local pid=$! t0
  t0=$(date +%s)
  while kill -0 "$pid" 2>/dev/null; do
    sleep 120
    if [ $(( $(date +%s) - t0 )) -ge "$cap" ]; then
      log "$tag : borne de temps atteinte, arret propre a $(tail -1 "$out/history.csv" 2>/dev/null | cut -d, -f1) s"
      kill "$pid" 2>/dev/null; sleep 5; break
    fi
  done
  wait "$pid" 2>/dev/null
  log "$tag termine : $(tail -1 "$out/history.csv" 2>/dev/null | cut -d, -f1) s, $(ls "$out" 2>/dev/null | grep -c 'fdem3d_joints_') trames"
}

log "file de la nuit : attente de la fin du St Anne grossier (rock25)"
while tasklist 2>/dev/null | grep -qiE "rockim_g1y1[89].exe"; do sleep 300; done
log "machine libre"
run_capped configs/stanne2025_rock137_visc0.cfg out_stanne2025_rock137 stanne2025_rock137 10800 "St Anne rock137 (localisation)"
run_capped configs/yang2026_kuru_train1_v5.cfg  out_yang_kuru_train1_v5 yang_kuru_train1_v5 10800 "Kuru train1 v5"
log "file de la nuit terminee"
