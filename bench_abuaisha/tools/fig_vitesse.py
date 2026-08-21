#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_vitesse.py — TRAJECTOIRES DE FISSURE SUR LE CHAMP DE VITESSE.
#
#   python bench_abuaisha/tools/fig_vitesse.py out_hfs_aniso --sh 6.8 --sv 4.6
#
# Reproduit la presentation de leurs figures 7, 10, 12, 13 et 15 : la norme de
# la vitesse en fond, les trajectoires de fissure par-dessus, une barre
# d'echelle de 10 cm, l'instant en haut a gauche. Leur echelle de couleur va
# de 0 a 0,12 m/s ; on la garde par defaut pour que les planches soient
# superposables, et le maximum reellement atteint est imprime.
#
# INSTANTS. Leur figure 12 montre 1,26 / 1,38 / 1,68 ms pour une rupture a
# 1,08 ms, soit +0,18 / +0,30 / +0,60 ms apres le pic. On prend les memes
# DECALAGES par rapport a notre propre pic : les deux campagnes n'ont pas la
# meme horloge (leur module de fluide n'est pas donne dans l'article, et le
# notre — eau a 2,2 GPa — fait monter la pression 2,7 fois moins vite).
#
# La vitesse est NODALE ; on la ramene a l'element par moyenne des trois
# sommets, ce qui est aussi ce que fait un rendu par cellules.
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
from matplotlib.collections import PolyCollection, LineCollection

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


def vg(x, n=2):
    return (("%." + str(n) + "f") % x).replace(".", ",")


def grab(txt, nm):
    m = re.search('Name="' + nm + '"[^>]*>(.*?)</DataArray>', txt, re.S)
    return np.fromstring(m.group(1), sep=" ")


def points(txt):
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--sh", type=float, default=6.8)
    ap.add_argument("--sv", type=float, default=4.6)
    ap.add_argument("--vmax", type=float, default=0.12,
                    help="haut de leur echelle de couleur [m/s]")
    ap.add_argument("--dt", type=float, nargs="+", default=[0.18, 0.30, 0.60],
                    help="decalages apres le pic [ms], comme leur fig. 12")
    ap.add_argument("--xmax", type=float, default=0.48)
    ap.add_argument("--ymax", type=float, default=0.30)
    ap.add_argument("--stem", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(ROOT, a.run)

    h = list(csv.DictReader(open(os.path.join(run, "history.csv"))))
    th = np.array([float(r["t"]) for r in h]) * 1e3
    ph = np.array([float(r["hydroP"]) for r in h]) / 1e6
    tpic = th[int(np.argmax(ph))]

    frw = list(csv.DictReader(open(os.path.join(run, "frames.csv"))))
    tfr = np.array([float(r["t"]) for r in frw]) * 1e3

    def kof(f):
        return int(os.path.basename(f).rsplit("_", 1)[1].split(".")[0])

    fs = sorted([f for f in glob.glob(os.path.join(run, "fdem_[0-9]*.vtu"))
                 if complete(f)], key=kof)
    js = {kof(f): f for f in glob.glob(
        os.path.join(run, "fdem_joints_[0-9]*.vtu")) if complete(f)}
    ks = [kof(f) for f in fs if kof(f) in js]
    cible = [tpic + d for d in a.dt]
    pris = [min(ks, key=lambda q: abs(tfr[q] - c)) for c in cible]
    if len(set(pris)) < len(pris):
        print("  (attention : deux instants tombent sur la meme trame)")

    # --- geometrie, lue une seule fois -------------------------------------
    t0 = io.open(fs[0], errors="ignore").read()
    P = points(t0)
    C = grab(t0, "connectivity").astype(int).reshape(-1, 3)
    cen = P[C].mean(axis=1)
    sel = np.where((np.abs(cen[:, 0] - CX) < a.xmax) &
                   (np.abs(cen[:, 1] - CY) < a.ymax))[0]
    verts = [(P[C[i]] - [CX, CY]) * 1e3 for i in sel]
    tj0 = io.open(js[ks[0]], errors="ignore").read()
    JP = points(tj0)
    JC = grab(tj0, "connectivity").astype(int).reshape(-1, 2)
    JSEG = (JP[JC] - [CX, CY]) * 1e3
    del t0, tj0

    n = len(pris)
    fig, axes = plt.subplots(n, 1, figsize=(8.4, 3.05 * n))
    if n == 1:
        axes = [axes]
    vmaxi = 0.0
    for ax, k in zip(axes, pris):
        txt = io.open([f for f in fs if kof(f) == k][0], errors="ignore").read()
        V = grab(txt, "velocity").reshape(-1, 3)[:, :2]
        mag = np.hypot(V[:, 0], V[:, 1])
        el = mag[C[sel]].mean(axis=1)          # nodal -> element
        vmaxi = max(vmaxi, float(el.max()))
        pc = PolyCollection(verts, array=el, cmap="coolwarm", lw=0,
                            clim=(0.0, a.vmax), rasterized=True)
        ax.add_collection(pc)

        d = grab(io.open(js[k], errors="ignore").read(), "damage")
        msk = d > 0.5                          # la TRAJECTOIRE, pas le halo
        ax.add_collection(LineCollection(JSEG[msk], colors="#d0006f", lw=0.9,
                                         zorder=3))
        ax.add_patch(plt.Circle((0.0, 0.0), RB * 1e3, fc="0.35", ec="0.2",
                                lw=0.6, zorder=4))

        j = int(np.argmin(np.abs(th - tfr[k])))
        ax.text(0.015, 0.93, "$t$ = " + vg(tfr[k]) + " ms   ($t_{pic}$ + "
                + vg(tfr[k] - tpic) + ")", transform=ax.transAxes,
                fontsize=10, va="top", color="0.15")
        ax.text(0.985, 0.93, "$p$ = " + vg(ph[j]) + " MPa",
                transform=ax.transAxes, fontsize=10, va="top", ha="right",
                color="0.15")
        # barre d'echelle de 10 cm, comme leurs figures
        x0, y0 = -a.xmax * 1e3 + 30, -a.ymax * 1e3 + 40
        ax.plot([x0, x0 + 100], [y0, y0], "k-", lw=1.6)
        ax.text(x0 + 50, y0 + 10, "10 cm", ha="center", fontsize=9)
        ax.set_xlim(-a.xmax * 1e3, a.xmax * 1e3)
        ax.set_ylim(-a.ymax * 1e3, a.ymax * 1e3)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.6)

    cb = fig.colorbar(pc, ax=axes, orientation="horizontal", fraction=0.055,
                      pad=0.035, aspect=46)
    cb.set_label("norme de la vitesse   [m/s]")
    cb.outline.set_linewidth(0.6)

    reg = "anisotrope" if abs(a.sh - a.sv) > 1e-9 else "isotrope"
    fig.suptitle("Trajectoires de fissure dans le champ de vitesse — état "
                 "de contrainte " + reg, fontsize=12.5, y=0.955)

    stem = a.stem or os.path.join(
        ROOT, "bench_abuaisha",
        "fig_vitesse_" + os.path.basename(run).replace("out_hfs_", ""))
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    fig.savefig(stem + ".png", bbox_inches="tight")
    print("écrit :", stem + ".pdf  et  .png")
    print("  pic à %s ms ; trames prises : %s"
          % (vg(tpic), ", ".join(vg(tfr[k]) + " ms" for k in pris)))
    print("  échelle affichée 0 – %s m/s (la leur) ; maximum réellement "
          "atteint %s m/s" % (vg(a.vmax), vg(vmaxi, 3)))


if __name__ == "__main__":
    main()
