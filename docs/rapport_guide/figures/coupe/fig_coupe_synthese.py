#!/usr/bin/env python3
"""Synthese du banc de coupe 2D (Heilman et al. 2024) : rapport d'injection,
vitesse nodale maximale et pic de force de coupe, archives et rejeu.

Sources (aucune valeur saisie sans provenance) :
- injection = « outil->solide » / « tool work output », relus dans les journaux
  d'archive phd_geothermie/FDEM/rockim/coupe_pdc/methode/run_cut_*.log et dans
  les journaux du rejeu copies dans donnees/ ;
- v max des runs d'archive : non imprimee par les journaux d'aout, reprise de
  RESULTATS_2026-08-18.md (v3 l. 45 ; nogc l. 180 ; a1, epfl, a2 l. 268-270) ;
- pic de F_x des runs d'archive : RESULTATS_2026-08-18.md l. 266-270 ;
- rejeu : journaux et history.csv de donnees/ (binaire f0209ef, 6 octobre 2026) ;
- cible : Heilman et al., fig. 3B, 21,5 kN sur 6,98 mm de largeur = 3,08 MN/m.

    python3 fig_coupe_synthese.py  ->  fig_coupe_synthese.pdf
"""
import os
import re
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
ARCH = "/home/user/phd_geothermie/FDEM/rockim/coupe_pdc/methode"
DON = os.path.join(HERE, "donnees")


def lit(log, cle):
    t = open(log, encoding="utf-8", errors="replace").read()
    m = re.search(cle + r"\s*:\s*([-\d.e+]+)", t)
    return float(m.group(1))


def injection(log):
    return lit(log, r"outil->solide") / lit(log, r"tool work output")


def pic_fx(run):
    h = np.genfromtxt(os.path.join(DON, run, "history.csv"), delimiter=",", names=True)
    return np.abs(h["toolFx"]).max() / 1e6


# (etiquette, journal, v max m/s, pic Fx MN/m)
RUNS = [
    ("v3", f"{ARCH}/run_cut_v3.log", 2544.0, 9.01),
    ("nogc", f"{ARCH}/run_cut_nogc.log", 314.0, 4.41),
    ("a1", f"{ARCH}/run_cut_a1.log", 7153.0, 6.44),
    ("epfl", f"{ARCH}/run_cut_epfl.log", 229.0, 4.33),
    ("a2", f"{ARCH}/run_cut_a2.log", 311.0, 0.798),
]
for lab, run in (("v3\nrejeu", "out_v3_rejeu"), ("Signorini\nrejeu", "out_v3sig_rejeu")):
    log = os.path.join(DON, run + ".log")
    RUNS.append((lab, log, lit(log, r"v nodale max"), pic_fx(run)))

lab = [r[0] for r in RUNS]
inj = np.array([injection(r[1]) for r in RUNS])
vb = np.array([r[2] for r in RUNS]) / 20.0
fx = np.array([r[3] for r in RUNS])
x = np.arange(len(RUNS))
col = ["#b2182b"] * 5 + ["#b2182b", "#2166ac"]
hat = [""] * 5 + ["//", "//"]

fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9))
for ax, y, tit, ref, reflab in (
        (axs[0], inj, "(a) injection outil / travail de corps rigide", 1.0, "1"),
        (axs[1], vb, r"(b) $v_\mathrm{max} / 2v_\mathrm{outil}$", 1.0, "borne 1")):
    ax.bar(x, y, color=col, hatch=hat, edgecolor="k", linewidth=0.4)
    ax.set_yscale("log")
    ax.axhline(ref, color="k", ls="--", lw=0.8)
    ax.axhline(2.0, color="0.5", ls=":", lw=0.8)
    ax.set_title(tit, fontsize=9)
    for xi, yi in zip(x, y):
        ax.text(xi, yi * 1.15, f"{yi:.3g}".replace(".", ","), ha="center", fontsize=6.5)
ax = axs[2]
ax.bar(x, fx, color=col, hatch=hat, edgecolor="k", linewidth=0.4)
ax.axhline(3.08, color="k", ls="--", lw=0.8)
ax.text(len(x) - 0.5, 3.25, "Heilman 3,08", ha="right", fontsize=7)
ax.set_title(r"(c) pic de $|F_x|$ (MN/m)", fontsize=9)
for xi, yi in zip(x, fx):
    ax.text(xi, yi + 0.15, f"{yi:.3g}".replace(".", ","), ha="center", fontsize=6.5)
for ax in axs:
    ax.set_xticks(x)
    ax.set_xticklabels(lab, fontsize=6.5, rotation=0)
    ax.tick_params(axis="y", labelsize=7)
fig.tight_layout(pad=0.4)
fig.savefig(os.path.join(HERE, "fig_coupe_synthese.pdf"))
for r, a, b, c in zip(lab, inj, vb, fx):
    print(r.replace("\n", " "), round(a, 3), round(b, 3), round(c, 3))
