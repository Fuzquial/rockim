#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_impact3d_oblique.py — vue oblique de haut, ENVELOPPES SEULES : la
# surface de rupture (facettes des joints rompus) et le contour du bloc.
# Rien d'autre : ni maillage, ni champ, ni element intact.
#
#   python tunnel_edz/plot_impact3d_oblique.py out_indent3d_ye [--elev 28]
#                                              [--azim -58] [--toolr 0.004]
#
# Les facettes de joint rompues forment, par leur union, l'enveloppe de la
# zone fracturee — c'est la surface que la fissuration a creee. On les trace
# telles quelles, coloriees par mode (traction / cisaillement), avec les douze
# aretes du bloc en fil de fer pour l'echelle et la silhouette de l'outil.
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_impact3d import read3d  # noqa: E402
from plot_tunnel_fields import complete  # noqa: E402

C_TEN, C_SHR = "#1B8A3A", "#C8342B"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_indent3d_ye")
    ap.add_argument("--elev", type=float, default=28.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--toolr", type=float, default=0.004)
    ap.add_argument("--zoom", type=float, default=0.010)
    ap.add_argument("--titre", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)
    out = os.path.join(HERE, os.path.basename(run) + "_enveloppe.png")

    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem3d_[0-9]*.vtu")))
          if complete(f)]
    jn = [f for f in sorted(glob.glob(os.path.join(run, "fdem3d_joints_[0-9]*.vtu")))
          if complete(f)]
    k = min(len(el), len(jn)) - 1
    P0, _, _ = read3d(el[0], [])
    W, D, H = P0[:, 0].max(), P0[:, 1].max(), P0[:, 2].max()
    cx, cy = 0.5 * W, 0.5 * D

    PJ, CJ, SJ = read3d(jn[k], ["breakMode"], ncell=3)
    bm = SJ["breakMode"]
    brk = bm > 0
    tri = PJ[CJ[brk]]
    mode = bm[brk]
    print(f"trame {k} : {int(brk.sum())} facettes rompues "
          f"({int((mode == 1).sum())} traction, {int((mode == 2).sum())} cisaillement)")

    fig = plt.figure(figsize=(15.5, 7.4))
    for i, (only, ttl) in enumerate([(None, "les deux modes"),
                                     (1, "traction seule"),
                                     (2, "cisaillement seul")]):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")
        sel = np.ones(len(mode), bool) if only is None else (mode == only)
        if sel.any():
            cols = [C_TEN if m == 1 else C_SHR for m in mode[sel]]
            pc = Poly3DCollection(tri[sel], facecolors=cols, edgecolors="none",
                                  alpha=0.85)
            ax.add_collection3d(pc)
        # les douze aretes du bloc, en fil de fer
        c = [(0, 0, 0), (W, 0, 0), (W, D, 0), (0, D, 0),
             (0, 0, H), (W, 0, H), (W, D, H), (0, D, H)]
        ed = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
              (0, 4), (1, 5), (2, 6), (3, 7)]
        ax.add_collection3d(Line3DCollection([[c[i0], c[i1]] for i0, i1 in ed],
                                             colors="0.55", lw=0.8))
        # silhouette de l'outil : cercle au niveau de la surface
        if a.toolr > 0:
            th = np.linspace(0, 2 * np.pi, 80)
            ax.plot(cx + a.toolr * np.cos(th), cy + a.toolr * np.sin(th),
                    np.full_like(th, H), color="#0B4F9E", lw=1.8)
        ax.set_xlim(cx - a.zoom, cx + a.zoom)
        ax.set_ylim(cy - a.zoom, cy + a.zoom)
        ax.set_zlim(H - 1.6 * a.zoom, H)
        ax.set_box_aspect((1, 1, 0.8))
        ax.view_init(elev=a.elev, azim=a.azim)
        ax.set_title(f"({chr(97+i)}) {ttl} — {int(sel.sum())} facettes",
                     fontsize=11)
        ax.set_xlabel("x [mm]", fontsize=8, labelpad=-4)
        ax.set_ylabel("y [mm]", fontsize=8, labelpad=-4)
        ax.tick_params(labelsize=6, pad=-2)
        ax.locator_params(nbins=4)
        # panneaux de fond retires : ils n'apportent rien et alourdissent
        for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
            pane.fill = False
            pane.set_edgecolor("0.9")
        ax.grid(False)
        ax.set_xticklabels([f"{v*1e3:.0f}" for v in ax.get_xticks()])
        ax.set_yticklabels([f"{v*1e3:.0f}" for v in ax.get_yticks()])
        ax.set_zticklabels([f"{v*1e3:.0f}" for v in ax.get_zticks()])

    fig.suptitle(a.titre or f"{os.path.basename(run)} — enveloppe de rupture, "
                 f"vue oblique de haut", fontsize=12.5)
    fig.tight_layout()
    fig.savefig(out, dpi=165)
    print("ecrit :", out)


if __name__ == "__main__":
    main()
