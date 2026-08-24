#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_metriques_compare.py — les quatre observables de l'etude, pour les trois
# schemas, a TEMPS DEPUIS LE RELACHEMENT EGAL (0,39 s).
#
#   python tunnel_edz/tools/fig_metriques_compare.py
# ---------------------------------------------------------------------------
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})

NOMS = ["adaptatif\n(reference)", "adaptatif\n+ pointe /1,6",
        "intrinsique"]
COUL = ["#7a7a7a", "#b3202f", "#1f4e79"]

# mesures faites sur les runs, a temps depuis le relachement egal
JOINTS = [21378, 39500, 35906]          # joints rompus
PROP = [43.7, 57.4, 58.9]               # % de propagation (post-traitement)
MONO = [75.7, 88.3, 84.9]               # % de blocs mono-element
GROS = [37.1, 492.6, 318.3]             # plus gros bloc detache [m2]

fig, ax = plt.subplots(1, 4, figsize=(13.6, 3.9))
X = np.arange(3)

for A, val, titre, unite, log in (
        (ax[0], JOINTS, "Joints rompus", "", False),
        (ax[1], PROP, "Propagation", "%", False),
        (ax[2], MONO, "Blocs d'un seul element", "%", False),
        (ax[3], GROS, "Plus gros bloc detache", r"m$^2$", True)):
    A.bar(X, val, color=COUL, width=0.62)
    for i, v in enumerate(val):
        A.text(i, v * (1.06 if log else 1.0) + (0 if log else max(val) * 0.02),
               ("%.0f" if v >= 100 else "%.1f") % v, ha="center",
               fontsize=9)
    A.set_xticks(X)
    A.set_xticklabels(NOMS, fontsize=8.5)
    A.set_title(titre + ("  [%s]" % unite if unite else ""), fontsize=10.5)
    if log:
        A.set_yscale("log")
        A.set_ylim(20, 1200)
    else:
        A.set_ylim(0, max(val) * 1.22)
    A.grid(axis="y", alpha=0.25)

ax[1].axhline(58.9, color="#1f4e79", ls="--", lw=0.9)
fig.suptitle("Les trois schemas d'insertion a temps depuis le relachement "
             "EGAL (0,39 s) — la pointe relachee amene l'adaptatif "
             "au comportement de l'intrinsique", fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.90])
for ext in ("pdf", "png"):
    fig.savefig("tunnel_edz/fig_metriques_compare." + ext, dpi=165)
print("ecrit : tunnel_edz/fig_metriques_compare")
