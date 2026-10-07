#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# morning_s1.sh — le dépouillement du run s = 1 (out_yang2026_v3) en une commande.
#   bash tools/morning_s1.sh [out_dir] [tag]
# Produit : results/fig/coupes_<tag>*.png, joints_<tag>*.png, joints_aretes_<tag>*.png,
# le rapport des sept critères et la ligne de la matrice fracture x conservation.
# Se rejoue à tout instant sur les trames déjà écrites (le run peut continuer).
# ---------------------------------------------------------------------------
set -u
OUT="${1:-out_yang2026_v3}"
TAG="${2:-s1_v3}"
LOG="results/$(basename "$OUT" | sed 's/^out_//').log"
cd "$(dirname "$0")/.." || exit 1
echo "== $OUT  ($(date '+%H:%M'))"
echo "-- dernier instant : $(tail -1 "$OUT/history.csv" | cut -d, -f1) s ; trames : $(ls "$OUT" | grep -c 'fdem3d_joints_')"
grep -oE "ENERGY ABORT.{0,140}" "$LOG" | head -1
python tools/law_matrix.py "$OUT"
python tools/yang_report.py "$OUT" "$LOG" 2>&1 | sed -n 1,12p
python tools/fig_joints_cuts.py "$OUT" --title "s = 1, deck v3-P (plastic + coulomb + ratchet)" --stem "results/fig/coupes_${TAG}" 2>&1 | grep -vi warning
python tools/fig_joints_only.py "$OUT" --title "s = 1, deck v3-P" --stem "results/fig/joints_aretes_${TAG}" 2>&1 | grep -vi warning
python tools/fig_joints_only.py "$OUT" --faces --title "s = 1, deck v3-P" --stem "results/fig/joints_${TAG}" 2>&1 | grep -vi warning
echo "== figures : results/fig/coupes_${TAG}.png, joints_aretes_${TAG}.png, joints_${TAG}.png"
