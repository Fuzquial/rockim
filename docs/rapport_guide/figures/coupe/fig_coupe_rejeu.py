#!/usr/bin/env python3
"""Rejeu du 6 octobre 2026 (binaire f0209ef) du deck v3 de coupe 2D : penalite
contre contact de Signorini avec plancher plat (etape 6).

Lit donnees/out_v3_rejeu/history.csv et donnees/out_v3sig_rejeu/history.csv
(colonnes t, toolFx, toolFy, toolX, nBroken, nFrag).
(a) force de coupe horizontale |F_x| en fonction de la course, cible Heilman ;
(b) joints rompus cumules.

    python3 fig_coupe_rejeu.py  ->  fig_coupe_rejeu.pdf
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 9})
# --- Style commun du rapport : Computer Modern (CMU Serif, a defaut Latin Modern
# Roman 10), virgule decimale sur toutes les graduations, signe moins ASCII.
import glob as _glob
import matplotlib.ticker as _mticker
from matplotlib import font_manager as _fm
for _f in (_glob.glob("/usr/share/fonts/**/cmun*.[ot]tf", recursive=True)
           + _glob.glob("/usr/share/texmf/fonts/opentype/public/lm/lmroman10-*.otf")
           + _glob.glob("/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-*.otf")):
    _fm.fontManager.addfont(_f)
_fm.fontManager.ttflist = [_e for _e in _fm.fontManager.ttflist
                           if _e.name != "Latin Modern Roman" or "lmroman10" in _e.fname]
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "Latin Modern Roman", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "axes.formatter.use_locale": False,
    "pdf.fonttype": 3,
})
_sf_call = _mticker.ScalarFormatter.__call__


def _sf_virgule(self, x, pos=None):
    s = _sf_call(self, x, pos)
    return s.replace(".", "{,}") if "$" in s else s.replace(".", ",")


_mticker.ScalarFormatter.__call__ = _sf_virgule
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DON = os.path.join(HERE, "donnees")
RUNS = [("out_v3_rejeu", "pénalité (deck v3)", "#b2182b"),
        ("out_v3sig_rejeu", "Signorini, plancher plat", "#2166ac"),
        ("out_v3sig_ratchet", "idem, cliquet des sécantes", "#4dac26")]
FACE = 1.5   # mm : face verticale de l'entaille du maillage cut_step

fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8))
for run, lab, c in RUNS:
    if not os.path.exists(os.path.join(DON, run, "history.csv")):
        continue
    h = np.genfromtxt(os.path.join(DON, run, "history.csv"), delimiter=",", names=True)
    course = h["toolX"] * 1e3 - FACE
    axs[0].plot(course, np.abs(h["toolFx"]) / 1e6, color=c, lw=0.6, label=lab)
    axs[1].plot(course, h["nBroken"], color=c, lw=1.0, label=lab)
    axs[1].plot(course[-1], h["nBroken"][-1], "x", color=c, ms=6)   # arret du garde-fou
axs[0].axhline(3.08, color="k", ls="--", lw=0.8)
axs[0].text(-0.08, 3.2, "pic de Heilman, 3,08 MN/m", fontsize=7)
axs[0].set_xlabel("avance dans la marche (mm)")
axs[0].set_ylabel(r"$|F_x|$ (MN/m)")
axs[0].set_title("(a) force de coupe", fontsize=9)
axs[1].set_xlabel("avance dans la marche (mm)")
axs[1].set_ylabel("joints rompus")
axs[1].set_title("(b) joints rompus (× : arrêt du garde-fou)", fontsize=9)
axs[1].legend(fontsize=7, frameon=False)
for ax in axs:
    ax.set_xlim(-0.1, 0.9)
fig.tight_layout(pad=0.4)
fig.savefig(os.path.join(HERE, "fig_coupe_rejeu.pdf"))
print("ecrit fig_coupe_rejeu.pdf")
