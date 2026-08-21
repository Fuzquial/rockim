#!/bin/bash
# ---------------------------------------------------------------------------
# B1 — balayage dampingLocal sur le cas de reference indent3d_grad.
#
# QUESTION POSEE : l'amortissement local de Cundall ne fait-il que MANGER de
# l'energie, ou le calcul s'appuie-t-il dessus pour tenir ?
#
# Mesure de reference (out_indent3d_grad, dampingLocal = 0,05) :
#   Cundall 0,0914 J  contre  joints 0,0035 J  ->  26 fois la physique.
#
# Les deux runs s'enchainent SEQUENTIELLEMENT (rockim prend tous les coeurs).
# Compter ~5 h chacun, possiblement plus a faible amortissement : moins
# d'amortissement = plus de debris = plus de paires de contact (lecon du
# 17/08, ou 124 fragments ont fait exploser le contact a 1,13 milliard de
# paires, 56 % du temps de calcul).
#
# Le moniteur d'energie est ARME (budgetAbortPct = 5) : un run qui divergerait
# s'arrete PROPREMENT en laissant son autopsie, au lieu de bruler 5 h. Il ne
# touche a aucun flottant de la physique — la comparaison avec la reference
# reste valide.
#
#   bash run_b1_damping.sh
# ---------------------------------------------------------------------------
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"

for TAG in d001 d000; do
  t0=$(date +%s)
  echo "=== DEBUT $TAG a $(date +%H:%M) ==="
  ./rockim_yy.exe configs/indent3d_grad_$TAG.cfg out_indent3d_grad_$TAG \
      > run_indent3d_grad_$TAG.log 2>&1
  echo "=== TERMINE en $(( $(date +%s) - t0 )) s ===" >> run_indent3d_grad_$TAG.log
  echo "=== FIN $TAG a $(date +%H:%M), $(( ($(date +%s) - t0) / 60 )) min ==="
done
