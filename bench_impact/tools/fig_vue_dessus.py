# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_vue_dessus.py — LA VRAIE vue de dessus : ce qu'on verrait en regardant
# le cratere d'en haut, sans voir au travers.
#
#   python bench_impact/tools/fig_vue_dessus.py out_imperial_coulomb 9 \
#          --stem bench_impact/figures/haut_B9 --t-us 123.7
#
# POURQUOI CE SCRIPT. La vue en plan de fig_fissure.py empile les facettes par
# CATEGORIE (maillage intact < zone de processus < facettes rompues, zorder
# fixe) : une facette profonde peut donc se dessiner PAR-DESSUS une facette de
# surface, et le regard traverse le cratere. C'est lisible comme carte de
# densite, ca ne l'est pas comme relief.
#
# Ici, algorithme du PEINTRE : les facettes sont triees par profondeur (la
# plus profonde dessinee en premier) et remplies de facon OPAQUE, sur un fond
# de roche intacte. Ce qui reste visible est donc exactement ce qui serait vu
# d'en haut — le sommet de chaque colonne de matiere, jamais ce qu'il y a
# dessous.
#
# La profondeur d'un triangle est prise a son sommet le PLUS HAUT (z max) :
# c'est ce point qui intercepte le regard en premier.
#
#   (a) le FACIES : rouge = cisaillement, bleu = traction, opaque ;
#   (b) le RELIEF : la meme scene coloree par la profondeur de ce qu'on voit —
#       lecture directe de la forme du cratere.
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import read_vtu, joints_frame, broken, frame_times, CX, CY, Z_SURF

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42})

BLEU, ROUGE = "#1f4e79", "#b22222"
ROCHE = "#d9ccb4"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("frame", type=int)
    ap.add_argument("--stem", default="bench_impact/figures/vue_dessus")
    ap.add_argument("--t-us", type=float, default=None)
    ap.add_argument("--half", type=float, default=12.0)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)

    p, con, f = read_vtu(joints_frame(a.run, a.frame))
    ctr, nrm, mode, P = broken(p, con, f)
    if len(P) == 0:
        sys.exit("aucune facette rompue dans cette frame")

    # --- projection (x, y) en mm et profondeur du SOMMET LE PLUS HAUT -------
    XY = np.stack([(P[:, :, 0] - CX) * 1e3, (P[:, :, 1] - CY) * 1e3], axis=2)
    z_haut = P[:, :, 2].max(axis=1)                 # le point qui est vu
    prof = (Z_SURF - z_haut) * 1e3                  # mm sous la surface

    # --- ALGORITHME DU PEINTRE : le plus profond d'abord -------------------
    ordre = np.argsort(-prof)                       # prof decroissante
    XY, prof, mode = XY[ordre], prof[ordre], mode[ordre]

    fig, ax = plt.subplots(1, 2, figsize=(9.8, 4.9))
    th = np.linspace(0, 2 * np.pi, 240)

    for k, A in enumerate(ax):
        # fond : la roche intacte, opaque — sans lui on verrait "a travers"
        # les interstices entre facettes.
        A.add_patch(plt.Circle((0, 0), a.half * 1.6, fc=ROCHE, ec="none",
                               zorder=0))
        A.set_xlim(-a.half, a.half)
        A.set_ylim(-a.half, a.half)
        A.set_aspect("equal")
        A.set_xlabel(r"$x$  [mm]")
        A.set_ylabel(r"$y$  [mm]")

    # (a) facies, opaque -----------------------------------------------------
    col = np.where(mode == 2, ROUGE, BLEU)
    ax[0].add_collection(PolyCollection(list(XY), facecolors=list(col),
                                        edgecolors="k", linewidths=0.25,
                                        zorder=2))
    n2 = int((mode == 2).sum())
    ax[0].set_title("(a) faciès vu de dessus — rouge = cisaillement",
                    fontsize=9.4, loc="left")
    ax[0].plot([], [], "s", color=ROUGE, ms=7,
               label="cisaillement — %d" % n2)
    ax[0].plot([], [], "s", color=BLEU, ms=7,
               label="traction — %d" % (len(mode) - n2))
    ax[0].legend(loc="lower right", fontsize=7.8, frameon=True,
                 framealpha=0.95, edgecolor="0.75")

    # (b) relief : profondeur de CE QU ON VOIT -------------------------------
    pc = PolyCollection(list(XY), array=prof, cmap="YlOrBr",
                        norm=plt.Normalize(0, np.percentile(prof, 97)),
                        edgecolors="k", linewidths=0.25, zorder=2)
    ax[1].add_collection(pc)
    cb = fig.colorbar(pc, ax=ax[1], fraction=0.046, pad=0.03)
    cb.set_label("profondeur de la facette vue  [mm]", fontsize=8.2)
    cb.ax.tick_params(labelsize=7.5)
    ax[1].set_title("(b) relief — la forme du cratère", fontsize=9.4,
                    loc="left")

    r_ext = np.hypot(XY[:, :, 0], XY[:, :, 1]).max()
    for A in ax:
        A.plot(r_ext * np.cos(th), r_ext * np.sin(th), ls=(0, (5, 3)),
               lw=1.0, color="#5b2d8e", zorder=4)
    ax[1].annotate("étendue : %.2f mm" % r_ext,
                   xy=(r_ext * 0.71, r_ext * 0.71),
                   xytext=(a.half * 0.16, a.half * 0.83), fontsize=8,
                   color="#5b2d8e",
                   arrowprops=dict(arrowstyle="->", lw=0.7, color="#5b2d8e"))

    ttl = "Le cratère vu d'en haut"
    if a.t_us:
        ttl += r"   ($t = %.1f\ \mu$s)" % a.t_us
    fig.suptitle(ttl, fontsize=11, y=0.99)
    fig.text(0.5, 0.005, "algorithme du peintre : les %d facettes rompues "
             "sont triées par profondeur et remplies de façon OPAQUE — on ne "
             "voit jamais ce qui est dessous" % len(XY),
             ha="center", fontsize=7.6, color="0.35", style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=200, bbox_inches="tight")
    print("ecrit : %s.pdf et .png" % a.stem)
    print("  %d facettes, %d cisaillement (%.0f %%)"
          % (len(mode), n2, 100 * n2 / len(mode)))
    print("  profondeur vue : %.2f a %.2f mm (mediane %.2f)"
          % (prof.min(), prof.max(), np.median(prof)))
    print("  etendue radiale : %.2f mm" % r_ext)
