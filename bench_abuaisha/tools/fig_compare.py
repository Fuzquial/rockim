#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_compare.py — COMPARAISON DES DEUX ETATS DE CONTRAINTE.
#
#   python bench_abuaisha/tools/fig_compare.py out_hfs_aniso out_hfs_iso
#
# Les deux calculs partagent le maillage (189 487 triangles), la roche, le
# schema, la pompe et l horloge (T = 4 ms, 60 trames). La SEULE difference
# est l etat de contrainte lointain : sigma'_H = 6,8 MPa en anisotrope,
# 4,6 en isotrope, sigma'_h = 4,6 dans les deux cas. Toute difference de
# resultat est donc imputable a cette seule cause.
#
# Quatre lectures :
#   (a) pression de puits, les deux courbes et les deux seuils analytiques ;
#   (b) DISTRIBUTION ANGULAIRE des joints rompus — la difference de facies,
#       rendue quantitative : deux lobes sur sigma'_H contre une repartition
#       etalee ;
#   (c) et (d) les facies, a la meme echelle et au meme instant final.
# ---------------------------------------------------------------------------
import argparse
import csv
import glob
import io
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
from plot_tunnel_fields import complete  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 10.5,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "savefig.dpi": 200,
})

CX, CY, RB = 4.0, 4.0, 0.05
FT = 5.0e6
BLEU, ROUGE = "#1f4e79", "#a11122"


def vg(x, n=2):
    return (("%." + str(n) + "f") % x).replace(".", ",")


def kof(f):
    return int(os.path.basename(f).rsplit("_", 1)[1].split(".")[0])


def grab(txt, nm):
    m = re.search('Name="' + nm + '"[^>]*>(.*?)</DataArray>', txt, re.S)
    return np.fromstring(m.group(1), sep=" ")


def lire(run):
    """Historique + derniere trame de joints."""
    d = {}
    h = list(csv.DictReader(open(os.path.join(run, "history.csv"))))
    d["t"] = np.array([float(r["t"]) for r in h]) * 1e3
    d["p"] = np.array([float(r["hydroP"]) for r in h]) / 1e6
    d["V"] = np.array([float(r["hydroVol"]) for r in h])
    d["nb"] = np.array([int(r["nBroken"]) for r in h])
    d["nw"] = np.array([int(float(r["hydroNWet"])) for r in h])
    js = sorted([f for f in glob.glob(
        os.path.join(run, "fdem_joints_[0-9]*.vtu")) if complete(f)], key=kof)
    txt = io.open(js[-1], errors="ignore").read()
    P = np.fromstring(re.search(
        r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", txt, re.S).group(1),
        sep=" ").reshape(-1, 3)[:, :2]
    C = grab(txt, "connectivity").astype(int).reshape(-1, 2)
    D = grab(txt, "damage")
    m = D > 0.5
    d["seg"] = (P[C[m]] - [CX, CY]) * 1e3
    d["D"] = D[m]
    mid = 0.5 * (P[C[m, 0]] + P[C[m, 1]]) - [CX, CY]
    d["r"] = np.hypot(mid[:, 0], mid[:, 1]) * 1e3
    # angle du joint rompu par rapport a l axe x, replie sur [0, 90] deg :
    # une fissure a +170 deg et une a -10 deg pointent dans la meme direction
    a = np.degrees(np.arctan2(mid[:, 1], mid[:, 0]))
    d["ang"] = np.minimum(np.abs(a), 180.0 - np.abs(a))
    d["ipic"] = int(np.argmax(d["p"]))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aniso")
    ap.add_argument("iso")
    ap.add_argument("--stem", default=None)
    a = ap.parse_args()
    A = lire(a.aniso if os.path.isabs(a.aniso) else os.path.join(ROOT, a.aniso))
    I = lire(a.iso if os.path.isabs(a.iso) else os.path.join(ROOT, a.iso))

    fig = plt.figure(figsize=(12.4, 9.6))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.24)

    # --- (a) pression de puits --------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    for d, lab, col, cib in ((A, "anisotrope", BLEU, 12.0),
                             (I, "isotrope", ROUGE, 14.2)):
        ax.plot(d["t"], d["p"], lw=1.8, color=col, label=lab)
        ax.axhline(cib, color=col, ls=(0, (4, 3)), lw=0.9, alpha=0.65)
        ax.plot(d["t"][d["ipic"]], d["p"][d["ipic"]], "o", ms=6,
                mfc="none", mec=col, mew=1.5)
    ax.text(0.12, 12.0, "seuil analytique 12,0", fontsize=8.5, color=BLEU,
            va="bottom")
    ax.text(0.12, 14.2, "seuil analytique 14,2", fontsize=8.5, color=ROUGE,
            va="bottom")
    ax.annotate(vg(A["p"][A["ipic"]]), (A["t"][A["ipic"]], A["p"][A["ipic"]]),
                textcoords="offset points", xytext=(-30, 4), fontsize=9.5,
                color=BLEU)
    ax.annotate(vg(I["p"][I["ipic"]]), (I["t"][I["ipic"]], I["p"][I["ipic"]]),
                textcoords="offset points", xytext=(6, 2), fontsize=9.5,
                color=ROUGE)
    ax.set_xlabel("temps  [ms]")
    ax.set_ylabel("pression de puits  [MPa]")
    ax.set_title("(a)  Pression de puits", loc="left")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.22, lw=0.6)

    # --- (b) distribution angulaire des joints rompus ---------------------
    ax = fig.add_subplot(gs[0, 1])
    bins = np.linspace(0, 90, 19)
    for d, lab, col in ((A, "anisotrope", BLEU), (I, "isotrope", ROUGE)):
        w = np.ones(len(d["ang"])) / max(len(d["ang"]), 1) * 100.0
        ax.hist(d["ang"], bins=bins, weights=w, histtype="step", lw=1.9,
                color=col, label=lab + " (%d joints)" % len(d["ang"]))
    ax.axhline(100.0 / (len(bins) - 1), color="0.45", ls=":", lw=1.0)
    ax.text(21, 100.0 / (len(bins) - 1) + 1.2, "répartition uniforme",
            fontsize=8.5, color="0.35")
    ax.set_xlabel(r"angle du joint rompu à l'axe de $\sigma^{\prime}_H$  [deg]")
    ax.set_ylabel("part des joints rompus  [%]")
    ax.set_title("(b)  Distribution angulaire de la fissuration", loc="left")
    ax.set_xlim(0, 90)
    ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.22, lw=0.6)

    # --- (c) et (d) les deux facies, MEME echelle -------------------------
    z = 1.06 * max(np.abs(A["seg"]).max(), np.abs(I["seg"]).max())
    for k, (d, lab, sh) in enumerate(((A, "anisotrope", 6.8),
                                      (I, "isotrope", 4.6))):
        ax = fig.add_subplot(gs[1, k])
        ax.add_patch(plt.Circle((0.0, 0.0), RB * 1e3, fc="0.94", ec="0.55",
                                lw=0.7, zorder=1))
        lc = LineCollection(d["seg"], array=d["D"], cmap="Greys", lw=1.5,
                            clim=(-0.25, 1.0), zorder=2)
        ax.add_collection(lc)
        ax.set_xlim(-z, z); ax.set_ylim(-z, z)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$   [mm]   ($\sigma^{\prime}_H$ = " + vg(sh, 1)
                      + " MPa)")
        ax.set_ylabel(r"$y$   [mm]")
        lettre = "(c)" if k == 0 else "(d)"
        ax.set_title(lettre + "  Faciès " + lab + r",  $t$ = "
                     + vg(d["t"][-1]) + " ms", loc="left")
        ax.grid(alpha=0.18, lw=0.5)
        ax.text(0.03, 0.03, "%d joints rompus" % d["nb"][-1],
                transform=ax.transAxes, fontsize=9.5, va="bottom",
                bbox=dict(fc="w", ec="0.75", alpha=0.9, pad=3))

    fig.suptitle("Fracturation hydraulique au forage : effet de l'état de "
                 "contrainte lointain", fontsize=12.5, y=0.965)
    fig.text(0.5, 0.932, "mêmes maillage, roche, schéma et pompe ; seule "
             r"$\sigma^{\prime}_H$ diffère (6,8 contre 4,6 MPa)",
             ha="center", fontsize=9.5, color="0.35")

    stem = a.stem or os.path.join(ROOT, "bench_abuaisha", "fig_compare_hf")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    fig.savefig(stem + ".png", bbox_inches="tight")
    print("écrit :", stem + ".pdf  et  .png")
    for d, lab, cib in ((A, "anisotrope", 12.0), (I, "isotrope", 14.2)):
        print("  %-11s pic %s MPa (cible %s, %+.1f %%) | %d joints | "
              "portée max %s mm | angle médian %s deg"
              % (lab, vg(d["p"][d["ipic"]]), vg(cib, 1),
                 (d["p"][d["ipic"]] - cib) / cib * 100, d["nb"][-1],
                 vg(d["r"].max(), 0), vg(float(np.median(d["ang"])), 0)))


if __name__ == "__main__":
    main()
