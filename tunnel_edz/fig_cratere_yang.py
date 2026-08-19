#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_cratere_yang.py — rendu du cratere et des fissures dans l esthetique des
# planches de Yang et al. (IJRMMS 191, 2025) : facettes de joint REMPLIES sur
# fond clair, vue de dessus et coupe axiale, pas de grille, pas de nuage de
# points.
#
#   python tunnel_edz/fig_cratere_yang.py out_imp3d_ultra [--frame 17]
#                                         [--out fig_cratere_ref_yang]
#
# Les joints sont des TRIANGLES en 3D : on les trace tels quels (PolyCollection)
# au lieu de leurs centroides. C est ce qui donne la lecture de facies — un
# nuage de points ne distingue pas une fissure plane d une zone diffuse.
# ---------------------------------------------------------------------------
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.labelsize": 10, "axes.titlesize": 10.5,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8.5,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# palette proche de leurs planches : brique sombre sur fond clair
C_BROKEN = "#8C2F27"
C_DAMAGE = "#E0B7A0"
C_ROCK   = "#EDEAE6"
C_ROCKED = "#C9C3BB"
C_INSERT = "#4A6FA5"


def read_joints(path):
    d = open(path, "r", errors="ignore").read()
    g = lambda p: re.search(p, d, re.S).group(1)
    P = np.fromstring(g(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>"),
                      sep=" ").reshape(-1, 3)
    C = np.fromstring(g(r'Name="connectivity"[^>]*>(.*?)</DataArray>'),
                      sep=" ", dtype=np.int64).reshape(-1, 3)
    dm = np.fromstring(g(r'Name="damage"[^>]*>(.*?)</DataArray>'), sep=" ")
    return P, C, dm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_imp3d_ultra")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--out", default=None)
    ap.add_argument("--half", type=float, default=6.0, help="demi-fenetre [mm]")
    ap.add_argument("--slab", type=float, default=1.2, help="epaisseur coupe [mm]")
    a = ap.parse_args()

    run = os.path.join(ROOT, a.run)
    js = sorted(f for f in os.listdir(run) if f.startswith("fdem3d_joints_"))
    fj = js[a.frame] if a.frame < 0 else js[a.frame]
    P, C, dm = read_joints(os.path.join(run, fj))

    tri = P[C]                                     # (nJ, 3, 3)
    cen = tri.mean(axis=1)
    brk = dm >= 1.0
    dmg = (dm > 0.05) & ~brk
    if not brk.any():
        print("aucun joint casse dans", fj); return

    # centre d impact = barycentre des joints casses ; surface = z max
    x0, y0 = cen[brk, 0].mean(), cen[brk, 1].mean()
    # Surface libre = z max des joints ROMPUS. Prendre le z max de TOUS les
    # points serait faux : le VTU porte aussi les joints de l INSERT, 22 mm
    # plus haut, et la coupe sortirait vide (erreur payee le 2026-08-19).
    zs = tri[brk][:, :, 2].max()
    H = a.half * 1e-3

    fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.6),
                           gridspec_kw=dict(width_ratios=[1.0, 1.15]))

    # ================= (a) vue de dessus =================================
    p = ax[0]
    p.set_facecolor(C_ROCK)
    sel = (np.abs(cen[:, 0] - x0) < H) & (np.abs(cen[:, 1] - y0) < H)
    for mask, col, z, al in ((dmg & sel, C_DAMAGE, 1, 0.75),
                             (brk & sel, C_BROKEN, 2, 0.95)):
        if not mask.any():
            continue
        pol = (tri[mask][:, :, :2] - [x0, y0]) * 1e3
        p.add_collection(PolyCollection(pol, facecolors=col, edgecolors="none",
                                        alpha=al, zorder=z))
    p.set_xlim(-a.half, a.half)
    p.set_ylim(-a.half, a.half)
    p.set_aspect("equal")
    p.set_xlabel("x [mm]")
    p.set_ylabel("y [mm]")
    p.set_title("(a) vue de dessus", fontsize=10.5)

    # ================= (b) coupe axiale ==================================
    p = ax[1]
    p.set_facecolor("white")
    S = a.slab * 1e-3
    sl = (np.abs(cen[:, 1] - y0) < S) & (np.abs(cen[:, 0] - x0) < H)
    # silhouette de roche
    p.add_patch(plt.Rectangle((-a.half, -1.35 * a.half), 2 * a.half,
                              1.35 * a.half, fc=C_ROCK, ec=C_ROCKED, lw=0.8,
                              zorder=0))
    for mask, col, z, al in ((dmg & sl, C_DAMAGE, 1, 0.8),
                             (brk & sl, C_BROKEN, 2, 0.95)):
        if not mask.any():
            continue
        pol = np.stack([(tri[mask][:, :, 0] - x0) * 1e3,
                        (tri[mask][:, :, 2] - zs) * 1e3], axis=-1)
        p.add_collection(PolyCollection(pol, facecolors=col, edgecolors="none",
                                        alpha=al, zorder=z))
    # insert (R = 11 mm) au repos sur la surface
    th = np.linspace(np.pi, 2 * np.pi, 200)
    R = 11.0
    p.plot(R * np.cos(th), R * np.sin(th) + R, color=C_INSERT, lw=2.0,
           zorder=3, solid_capstyle="round")
    p.axhline(0.0, color=C_ROCKED, lw=0.9, zorder=1)
    p.set_xlim(-a.half, a.half)
    p.set_ylim(-1.35 * a.half, 0.42 * a.half)
    p.set_aspect("equal")
    p.set_xlabel("x [mm]")
    p.set_ylabel("profondeur sous la surface  [mm]")
    p.set_title("(b) coupe axiale, tranche |y| < %.1f mm" % a.slab, fontsize=10.5)

    for p in ax:
        for s in ("top", "right"):
            p.spines[s].set_visible(False)
        p.tick_params(length=3)

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    ax[1].legend(handles=[Patch(fc=C_BROKEN, label="joint rompu (D = 1)"),
                          Patch(fc=C_DAMAGE, label="endommagé (D > 0,05)"),
                          Line2D([], [], color=C_INSERT, lw=2, label="insert")],
                 loc="lower right", framealpha=0.95)

    nb, nd = int(brk.sum()), int(dmg.sum())
    fig.suptitle("%s — trame %s : %d joints rompus, %d endommagés"
                 % (a.run, re.search(r"(\d+)\.vtu", fj).group(1), nb, nd),
                 fontsize=11, y=1.0)
    fig.tight_layout()
    out = a.out or ("fig_cratere_%s_yang" % a.run.replace("out_", ""))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, out + "." + ext), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print("ecrit :", out, "|", nb, "rompus,", nd, "endommages")


if __name__ == "__main__":
    main()
