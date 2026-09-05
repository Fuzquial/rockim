#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_bilan.py — LA planche du run de reference : ou va l energie, et ou en
# sont les sept criteres de Yang et al.
#
#   python bench_impact/tools/fig_bilan.py out_temoin_adap --stem fig_bilan
#
# Deux panneaux, une idee chacun :
#   (a) la PARTITION DE L ENERGIE, rockim contre les chiffres publies
#       (ARMA 24-0952 : fissuration 2,6 %, frottement 64,9 % de 49,3 J) ;
#   (b) les SEPT CRITERES de leur Table 3, chacun rapporte a sa fourchette.
#
# Les postes du panneau (a) sont lus dans le RESUME du run (stdout), pas
# recalcules : ce sont ceux que le solveur imprime lui-meme.
# ---------------------------------------------------------------------------
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "font.size": 10,
})

BLEU, ROUGE, GRIS, VERT = "#1f4e79", "#b22222", "#9a9a9a", "#2e7d32"

# Partition publiee (ARMA 24-0952, St Anne a 9,41 m/s, 49,3 J entrants).
PUB = {"fissuration": 0.026, "frottement": 0.649}

# Les sept criteres : (etiquette, unite, borne basse, borne haute)
CRIT = [
    ("contrainte max au bit",    "MPa", 200.0, 260.0),
    ("vitesse d'indentation",    "m/s",   9.40,  9.85),
    ("vitesse de rebond",        "m/s",   6.87,  7.10),
    ("rebond / indentation",     "",      0.72,  0.73),
    ("profondeur d'indentation", "mm",    1.45,  1.60),
    ("fissure radiale max",      "mm",   20.2,  24.5),
    ("rayon de cratere",         "mm",   10.0,  12.1),
]


def lire_resume(log):
    """Les postes du bilan, lus dans le resume imprime par le solveur."""
    t = open(log, encoding="utf-8", errors="replace").read()

    def f(rx, d=float("nan")):
        m = re.search(rx, t)
        return float(m.group(1)) if m else d

    return dict(
        ke0=f(r"energy budget \(V2/B4\): KE ([\d.eE+-]+)"),
        ke1=f(r"energy budget \(V2/B4\): KE [\d.eE+-]+ -> ([\d.eE+-]+)"),
        visq=f(r"dont visqueux \(2 mu D\) : ([\d.eE+-]+) J"),
        fric=f(r"contact *: [\d.eE+-]+ J \(dont frottement ([\d.eE+-]+) J"),
        contact=f(r"contact *: ([\d.eE+-]+) J"),
        joints=f(r"joints *: (-?[\d.eE+-]+) J cohesif"),
        resid=f(r"residu *: (-?[\d.eE+-]+) J"),
        nbrok=f(r"broken joints *: (\d+)"),
        shear=f(r"breakage mode: \d+ tensile, \d+ shear \(([\d.eE+-]+) %"),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--log", default=None)
    ap.add_argument("--stem", default="fig_bilan")
    ap.add_argument("--mes", nargs=7, type=float, required=True,
                    help="les 7 criteres mesures, dans l'ordre de CRIT")
    a = ap.parse_args()

    log = a.log or ("log_" + os.path.basename(a.run).replace("out_", "") + ".txt")
    b = lire_resume(log)
    ke0 = b["ke0"]

    fig, (A, B) = plt.subplots(1, 2, figsize=(11.6, 4.6))
    fig.suptitle("Impact a insert unique, calcaire St Anne — run de reference "
                 "rockim contre Yang $et\\ al.$", fontsize=12.5)

    # ---- (a) ou va l'energie ------------------------------------------------
    postes = ["viscosite", "frottement\nde contact", "joints\n(fissuration)"]
    mes = [b["visq"], b["fric"], abs(b["joints"])]
    pub = [0.0, PUB["frottement"] * ke0, PUB["fissuration"] * ke0]
    x = np.arange(3)
    w = 0.36
    A.bar(x - w / 2, mes, w, color=BLEU, label="rockim (reference)")
    A.bar(x + w / 2, pub, w, color=ROUGE, alpha=0.85,
          label="Yang $et\\ al.$ (ARMA 2024)")
    for xi, v in zip(x - w / 2, mes):
        A.annotate("%.1f" % v, (xi, v), ha="center", va="bottom", fontsize=9)
    for xi, v in zip(x + w / 2, pub):
        if v > 0:
            A.annotate("%.1f" % v, (xi, v), ha="center", va="bottom",
                       fontsize=9, color=ROUGE)
    A.annotate("non publie\n(~0 net)", (x[0] + w / 2, 1.0), ha="center",
               va="bottom", fontsize=8, color=ROUGE, style="italic")
    A.set_xticks(x)
    A.set_xticklabels(postes)
    A.set_ylabel("energie [J]")
    A.set_title("(a)  Ou passe l'energie d'impact  ($KE_0$ = %.1f J)" % ke0,
                loc="left", fontsize=11)
    A.legend(frameon=False, fontsize=9)
    A.grid(axis="y", lw=0.4, alpha=0.4)
    A.set_axisbelow(True)

    # ---- (b) les sept criteres ---------------------------------------------
    noms = [c[0] for c in CRIT]
    y = np.arange(len(CRIT))[::-1]
    for yi, (nom, u, lo, hi), m in zip(y, CRIT, a.mes):
        mid = 0.5 * (lo + hi)
        rl, rh, rm = lo / mid, hi / mid, m / mid
        B.barh(yi, rh - rl, left=rl, height=0.5, color=ROUGE, alpha=0.28)
        ok = rl <= rm <= rh
        B.plot([rm], [yi], "o", ms=7, color=BLEU if ok else "#c25e00",
               zorder=3)
        B.annotate("%.3g %s" % (m, u), (rm, yi), xytext=(0, 9),
                   textcoords="offset points", ha="center", fontsize=8.5,
                   color=BLEU if ok else "#c25e00")
    B.axvline(1.0, color="k", lw=0.7, ls="--")
    B.set_yticks(y)
    B.set_yticklabels(noms, fontsize=9)
    B.set_xlabel("mesure / centre de la fourchette publiee")
    B.set_title("(b)  Les sept criteres de leur Table 3", loc="left",
                fontsize=11)
    B.grid(axis="x", lw=0.4, alpha=0.4)
    B.set_axisbelow(True)
    B.set_xlim(0, 1.6)

    fig.text(0.5, 0.005,
             "rockim : %d joints rompus dont %.2f %% en cisaillement  |  "
             "residu du bilan %.3g J  |  la viscosite dissipe %.0f %% de "
             "l'energie perdue"
             % (b["nbrok"], b["shear"], b["resid"],
                100 * b["visq"] / (b["ke0"] - b["ke1"])),
             ha="center", fontsize=8.5, color="#444444")

    fig.tight_layout(rect=[0, 0.03, 1, 0.94])
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=180)
    print("ecrit : %s.pdf et .png" % a.stem)
    print("  viscosite      %8.2f J   (%.0f %% de l'energie perdue)"
          % (b["visq"], 100 * b["visq"] / (b["ke0"] - b["ke1"])))
    print("  frottement     %8.2f J   contre %.1f J attendus"
          % (b["fric"], PUB["frottement"] * ke0))
    print("  joints         %8.2f J   contre %.1f J attendus"
          % (abs(b["joints"]), PUB["fissuration"] * ke0))


if __name__ == "__main__":
    main()
