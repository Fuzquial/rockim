#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_imperial_live.py — l'etat A MI-VOL de la replique Imperial College.
#
#   python bench_impact/tools/fig_imperial_live.py out_imperial \
#          --stem bench_impact/fig_imperial_live
#
# Trois panneaux, une idee chacun, et TOUS lisent l'INSERT et non le bit :
# la lecon du 2026-08-26 est que le deplacement du bit n'est PAS la
# penetration — le bit se fait comprimer par le piston a son sommet pendant
# que l'insert, 24 cm plus bas, n'a pas encore bouge.
#
#   (a) CHRONOLOGIE : penetration de l'insert et sa vitesse, avec les deux
#       instants PREDITS marques — impact du piston (jeu/vitesse) et arrivee
#       de l'onde a l'insert (longueur/celerite acier). Le test est qu'ils
#       tombent sur les ruptures de pente observees.
#   (b) LA JAUGE a mi-bit contre leur fourchette publiee 200-260 MPa.
#   (c) LES SEPT CRITERES et lesquels sont DEJA acquis.
#
# PAS de courbe force-penetration ici, et c'est deliberé : une reconstruction
# par Newton sur le seul insert ignore la BRASURE qui le lie au bit
# (groupBond.bit.insert), donc elle ne mesure pas la force insert-roche. La
# version fausse, essayee le 2026-08-26, donnait -150 kN. Le F-delta correct
# se fait sur le corps bit+insert+piston une fois le run fini (fig_fp.py).
#
# Le run n'est PAS fini : la figure porte son avancement en clair. C'est un
# constat d'etape, pas un livrable.
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import history

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "font.size": 10,
})

BLEU, ROUGE, GRIS, VERT = "#1f4e79", "#b22222", "#9a9a9a", "#2e7d32"
ORANGE = "#c25e00"

E_ST, RHO_ST = 200e9, 7850.0
GAP_PISTON, V_PISTON, L_BIT, H_INS = 2.0e-4, 10.66, 0.265, 0.0232

# leur Table 3 : (etiquette, unite, borne basse, borne haute, mesurable ?)
CRIT = [
    ("contrainte max au bit",    "MPa", 200.0, 260.0),
    ("vitesse d'indentation",    "m/s",   9.40,  9.85),
    ("vitesse de rebond",        "m/s",   6.87,  7.10),
    ("rebond / indentation",     "",      0.72,  0.73),
    ("profondeur d'indentation", "mm",    1.45,  1.60),
    ("fissure radiale max",      "mm",   20.2,  24.5),
    ("rayon de cratere",         "mm",   10.0,  12.1),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_imperial_live")
    ap.add_argument("--tfin", type=float, default=5.5e-4)
    a = ap.parse_args()

    h = history(a.run)
    t = h["t"]
    tus = t * 1e6
    pen = (h["z_insert"][0] - h["z_insert"]) * 1e3          # mm
    vin = -h["vz_insert"]                                   # m/s, > 0 = enfonce
    sig = -h["szz_bit"] / 1e6                               # MPa, compression > 0
    nb = int(h["nBroken"][-1]) if "nBroken" in h else 0

    # instants PREDITS, independamment du calcul
    t_imp = GAP_PISTON / V_PISTON * 1e6                     # 18,8 us
    c_st = np.sqrt(E_ST / RHO_ST)                           # 5048 m/s
    t_onde = t_imp + (L_BIT - H_INS) / c_st * 1e6           # 66,7 us

    sigpic = float(np.max(sig))
    vmax = float(np.max(vin))

    fig = plt.figure(figsize=(13.6, 5.6))
    gs = fig.add_gridspec(1, 3, wspace=0.42, left=0.06, right=0.975,
                          top=0.80, bottom=0.16)
    fig.suptitle("Replique Imperial College — impact a insert unique, "
                 "calcaire St Anne", fontsize=13, y=0.965)
    fig.text(0.5, 0.905,
             "ETAT A MI-VOL : $t$ = %.1f $\\mu$s sur %.0f  (%.1f %%) — "
             "l'indentation vient de commencer, %d joint(s) rompu(s)"
             % (tus[-1], a.tfin * 1e6, 100 * t[-1] / a.tfin, nb),
             ha="center", fontsize=10, color="#444444", style="italic")

    # ---- (a) chronologie ---------------------------------------------------
    A = fig.add_subplot(gs[0, 0])
    A.plot(tus, pen, color=BLEU, lw=2.0, zorder=3)
    A.set_xlabel(r"temps  [$\mu$s]")
    A.set_ylabel("penetration de l'insert  [mm]", color=BLEU)
    A.tick_params(axis="y", labelcolor=BLEU)
    A.set_ylim(-0.004, max(pen.max() * 1.35, 0.02))
    A2 = A.twinx()
    A2.axhspan(9.40, 9.85, color=ROUGE, alpha=0.15, zorder=0)
    A2.plot(tus, vin, color=ROUGE, lw=1.3, alpha=0.9, zorder=2)
    A2.set_ylabel("vitesse de l'insert  [m/s]", color=ROUGE)
    A2.tick_params(axis="y", labelcolor=ROUGE)
    A2.set_ylim(-0.6, 11.5)
    A2.annotate("leur fourchette\n9,40 - 9,85 m/s", (2, 9.62),
                fontsize=7.5, color=ROUGE, va="center")
    ytop = A.get_ylim()[1]
    for x, lab, dx in ((t_imp, "impact piston", 4),
                       (t_onde, "onde a l'insert", -4)):
        A.axvline(x, color=GRIS, lw=1.0, ls="--", zorder=1)
        A.annotate(lab, (x, ytop), xytext=(dx, -4),
                   textcoords="offset points", fontsize=8, color="#555555",
                   va="top", ha="left" if dx > 0 else "right", rotation=90)
    A.set_title("(a)  Les deux instants PREDITS tombent sur\n"
                "les ruptures de pente observees", loc="left", fontsize=10.5)

    # ---- (b) la jauge contre leur fourchette -------------------------------
    B = fig.add_subplot(gs[0, 1])
    B.axhspan(200.0, 260.0, color=VERT, alpha=0.18, zorder=0)
    B.annotate("leur Table 3 : 200 - 260 MPa", (2, 230), fontsize=8,
               color=VERT, va="center")
    B.plot(tus, sig, color=BLEU, lw=1.5)
    k = int(np.argmax(sig))
    B.plot([tus[k]], [sig[k]], "o", ms=7, color=ROUGE, zorder=3)
    B.annotate("pic %.0f MPa" % sig[k], (tus[k], sig[k]),
               xytext=(8, 8), textcoords="offset points", fontsize=9.5,
               color=ROUGE)
    B.axhline(0, color="k", lw=0.5)
    B.set_xlabel(r"temps  [$\mu$s]")
    B.set_ylabel("contrainte de section a mi-bit  [MPa]")
    B.set_ylim(-20, 290)
    B.set_title("(b)  L'amplitude de l'onde tombe DANS\nleur fourchette "
                "publiee", loc="left", fontsize=10.5)

    # ---- (c) les sept criteres, et ce qui est deja acquis -------------------
    C = fig.add_subplot(gs[0, 2])
    mes = {0: sigpic, 1: vmax}          # les deux seuls deja mesurables
    y = np.arange(len(CRIT))[::-1]
    for yi, (nom, u, lo, hi) in zip(y, CRIT):
        mid = 0.5 * (lo + hi)
        C.barh(yi, hi / mid - lo / mid, left=lo / mid, height=0.5,
               color=ROUGE, alpha=0.28)
    for i, v in mes.items():
        lo, hi = CRIT[i][2], CRIT[i][3]
        mid = 0.5 * (lo + hi)
        rm = v / mid
        ok = lo <= v <= hi
        C.plot([rm], [y[i]], "o", ms=9, color=BLEU if ok else ORANGE, zorder=3)
        C.annotate("%.4g %s" % (v, CRIT[i][1]), (rm, y[i]), xytext=(0, 11),
                   textcoords="offset points", ha="center", fontsize=9,
                   color=BLEU if ok else ORANGE, fontweight="bold")
    for i in range(len(CRIT)):
        if i not in mes:
            C.annotate("pas encore mesurable", (1.0, y[i]), xytext=(0, -2),
                       textcoords="offset points", ha="center", fontsize=7.5,
                       color="#999999", style="italic", va="top")
    C.axvline(1.0, color="k", lw=0.8, ls="--")
    C.set_yticks(y)
    C.set_yticklabels([c[0] for c in CRIT], fontsize=9)
    C.set_xlabel("mesure / centre de la fourchette publiee")
    C.set_xlim(0.55, 1.45)
    C.grid(axis="x", lw=0.4, alpha=0.4)
    C.set_axisbelow(True)
    # Le titre est DYNAMIQUE : il a menti une fois (2026-08-26, la vitesse
    # d indentation a traverse la fourchette puis l a depassee entre deux
    # tirages de la figure). Un titre en dur sur un run VIVANT est un piege.
    nok = sum(1 for i, v in mes.items() if CRIT[i][2] <= v <= CRIT[i][3])
    C.set_title("(c)  Deux criteres deja mesurables :\n%d dans la "
                "fourchette, %d hors" % (nok, len(mes) - nok),
                loc="left", fontsize=10.5)

    fig.text(0.5, 0.035,
             "insert : %.4f mm enfonces a %.2f m/s  |  joints rompus : %d  |  "
             "instants predits : impact %.1f $\\mu$s, onde %.1f $\\mu$s "
             "($c_{acier}$ = %.0f m/s)  |  les 5 autres criteres exigent "
             "l'indentation complete"
             % (pen[-1], vin[-1], nb, t_imp, t_onde, c_st),
             ha="center", fontsize=8.5, color="#444444")

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=180)
    print("ecrit : %s.pdf et .png" % a.stem)
    print("  t = %.1f us (%.1f %%)   penetration insert = %.4f mm"
          % (tus[-1], 100 * t[-1] / a.tfin, pen[-1]))
    print("  v_insert max = %.2f m/s   (leur fourchette 9,40 - 9,85)  %s"
          % (vmax, "DANS" if 9.40 <= vmax <= 9.85 else "hors"))
    print("  pic de jauge = %.0f MPa   (leur fourchette 200 - 260)     %s"
          % (sigpic, "DANS" if 200 <= sigpic <= 260 else "hors"))
    print("  joints rompus = %d" % nb)


if __name__ == "__main__":
    main()
