#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_impact_planche.py — planche complete d un impact 3D, dans l esthetique
# des figures de Yang, Xiang, Naderi, Wang, Aising, Ugarte & Latham
# (IJRMMS 191, 2025, 106125) : facettes de joint REMPLIES, palette brique sur
# fond clair, perspective + vue de dessus + coupe axiale + portees radiales.
#
#   python tunnel_edz/fig_impact_planche.py out_impact_uni [--frame -1]
#          [--half 22] [--slab 1.5] [--out fig_impact]
#
# POURQUOI DES FACETTES ET PAS DES POINTS : un nuage de centroides ne
# distingue pas une fissure PLANE d une zone d endommagement diffuse. Les
# fissures radiales sont precisement des surfaces quasi-planes rayonnant du
# cratere : elles n apparaissent qu en trace remplie.
#
# Les quatre panneaux :
#   (a) perspective 3D des joints rompus, couleur = profondeur — leur fig. 3
#   (b) vue de dessus, rompus sur fond d endommages — le bol du cratere
#   (c) coupe axiale — mediane sous l insert, laterales, profondeur
#   (d) portee radiale par secteur angulaire — la longueur de fissure radiale
#       qui est un de leurs sept criteres (leur Table 3 : 20 a 24,5 mm)
# ---------------------------------------------------------------------------
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.labelsize": 9.5, "axes.titlesize": 10,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.5,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

C_BROKEN, C_DAMAGE = "#8C2F27", "#E0B7A0"
C_ROCK, C_EDGE, C_INS = "#EDEAE6", "#C9C3BB", "#4A6FA5"


def read_joints(path):
    d = open(path, "r", errors="ignore").read()
    g = lambda p: re.search(p, d, re.S).group(1)
    P = np.fromstring(g(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>"),
                      sep=" ").reshape(-1, 3)
    C = np.fromstring(g(r'Name="connectivity"[^>]*>(.*?)</DataArray>'),
                      sep=" ", dtype=np.int64).reshape(-1, 3)
    dm = np.fromstring(g(r'Name="damage"[^>]*>(.*?)</DataArray>'), sep=" ")
    bm = np.fromstring(g(r'Name="breakMode"[^>]*>(.*?)</DataArray>'), sep=" ") \
        if 'Name="breakMode"' in d else np.zeros(len(C))
    return P, C, dm, bm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_impact_uni")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--half", type=float, default=22.0, help="demi-fenetre [mm]")
    ap.add_argument("--slab", type=float, default=1.5, help="demi-tranche [mm]")
    ap.add_argument("--sectors", type=int, default=24)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    run = os.path.join(ROOT, a.run)
    js = sorted(f for f in os.listdir(run) if f.startswith("fdem3d_joints_"))
    if not js:
        print("aucun VTU de joints dans", run); return
    fj = js[a.frame]
    P, C, dm, bm = read_joints(os.path.join(run, fj))

    tri = P[C]
    cen = tri.mean(axis=1)
    brk = dm >= 1.0
    dmg = (dm > 0.05) & ~brk
    if not brk.any():
        print("aucun joint rompu dans", fj, "— rien a tracer"); return

    x0, y0 = cen[brk, 0].mean(), cen[brk, 1].mean()
    zs = tri[brk][:, :, 2].max()          # surface libre = joints rompus les
    H = a.half * 1e-3                     # plus hauts (jamais le z max global,
    S = a.slab * 1e-3                     # qui appartient a l INSERT)

    fig = plt.figure(figsize=(12.4, 9.6))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.22)

    # ================= (a) perspective 3D ================================
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    sel = brk & (np.abs(cen[:, 0] - x0) < H) & (np.abs(cen[:, 1] - y0) < H)
    T = (tri[sel] - [x0, y0, zs]) * 1e3
    if len(T):
        d = T[:, :, 2].mean(axis=1)
        norm = plt.Normalize(d.min(), 0.0)
        cols = plt.cm.YlOrRd_r(norm(d))
        pc = Poly3DCollection(T, facecolors=cols, edgecolors="none", alpha=0.92)
        ax.add_collection3d(pc)
    ax.set_xlim(-a.half, a.half); ax.set_ylim(-a.half, a.half)
    ax.set_zlim(-1.6 * a.half, 0.25 * a.half)
    ax.set_box_aspect((1, 1, 0.9))
    ax.view_init(elev=22, azim=-58)
    ax.set_xlabel("x [mm]", labelpad=-4); ax.set_ylabel("y [mm]", labelpad=-4)
    ax.set_zlabel("z [mm]", labelpad=-4)
    ax.tick_params(labelsize=7, pad=-2)
    ax.set_title("(a) joints rompus en perspective", fontsize=10)

    # ================= (b) vue de dessus =================================
    bx = fig.add_subplot(gs[0, 1])
    bx.set_facecolor(C_ROCK)
    win = (np.abs(cen[:, 0] - x0) < H) & (np.abs(cen[:, 1] - y0) < H)
    for m, c, z, al in ((dmg & win, C_DAMAGE, 1, 0.7), (brk & win, C_BROKEN, 2, 0.95)):
        if m.any():
            bx.add_collection(PolyCollection(
                (tri[m][:, :, :2] - [x0, y0]) * 1e3,
                facecolors=c, edgecolors="none", alpha=al, zorder=z))
    bx.set_xlim(-a.half, a.half); bx.set_ylim(-a.half, a.half)
    bx.set_aspect("equal"); bx.set_xlabel("x [mm]"); bx.set_ylabel("y [mm]")
    bx.set_title("(b) vue de dessus", fontsize=10)

    # ================= (c) coupe axiale ==================================
    cx = fig.add_subplot(gs[1, 0])
    sl = (np.abs(cen[:, 1] - y0) < S) & (np.abs(cen[:, 0] - x0) < H)
    cx.add_patch(plt.Rectangle((-a.half, -1.6 * a.half), 2 * a.half,
                               1.6 * a.half, fc=C_ROCK, ec=C_EDGE, lw=0.8, zorder=0))
    for m, c, z, al in ((dmg & sl, C_DAMAGE, 1, 0.8), (brk & sl, C_BROKEN, 2, 0.95)):
        if m.any():
            pol = np.stack([(tri[m][:, :, 0] - x0) * 1e3,
                            (tri[m][:, :, 2] - zs) * 1e3], axis=-1)
            cx.add_collection(PolyCollection(pol, facecolors=c,
                                             edgecolors="none", alpha=al, zorder=z))
    th = np.linspace(np.pi, 2 * np.pi, 200)
    R = 8.51
    cx.plot(R * np.cos(th), R * np.sin(th) + R, color=C_INS, lw=2.0, zorder=3)
    cx.axhline(0.0, color=C_EDGE, lw=0.9, zorder=1)
    cx.set_xlim(-a.half, a.half); cx.set_ylim(-1.6 * a.half, 0.35 * a.half)
    cx.set_aspect("equal"); cx.set_xlabel("x [mm]")
    cx.set_ylabel("profondeur sous la surface [mm]")
    cx.set_title("(c) coupe axiale, tranche |y| < %.1f mm" % a.slab, fontsize=10)

    # ================= (d) portee radiale par secteur ====================
    dx = fig.add_subplot(gs[1, 1], projection="polar")
    rr = np.hypot(cen[brk, 0] - x0, cen[brk, 1] - y0) * 1e3
    tt = np.arctan2(cen[brk, 1] - y0, cen[brk, 0] - x0)
    bins = np.linspace(-np.pi, np.pi, a.sectors + 1)
    idx = np.digitize(tt, bins) - 1
    port = np.array([rr[idx == k].max() if (idx == k).any() else 0.0
                     for k in range(a.sectors)])
    ctr = 0.5 * (bins[:-1] + bins[1:])
    dx.bar(ctr, port, width=2 * np.pi / a.sectors, color=C_BROKEN, alpha=0.85,
           edgecolor="white", lw=0.6)
    rcr = np.percentile(rr, 95)
    dx.plot(np.linspace(-np.pi, np.pi, 200), [rcr] * 200, color=C_INS,
            lw=1.6, ls="--")
    dx.set_title("(d) portée radiale par secteur [mm]\n"
                 "tireté bleu : rayon de cratère (p95) = %.1f mm" % rcr,
                 fontsize=10, pad=14)
    dx.set_theta_zero_location("E")
    dx.tick_params(labelsize=7.5)

    from matplotlib.patches import Patch
    bx.legend(handles=[Patch(fc=C_BROKEN, label="joint rompu (D = 1)"),
                       Patch(fc=C_DAMAGE, label="endommagé (D > 0,05)")],
              loc="upper right", framealpha=0.95)

    nb, nd = int(brk.sum()), int(dmg.sum())
    fig.suptitle("%s — trame %s : %d joints rompus, %d endommagés  |  "
                 "portée radiale max %.1f mm, médiane des secteurs %.1f mm"
                 % (a.run, re.search(r"(\d+)\.vtu", fj).group(1), nb, nd,
                    port.max(), np.median(port)), fontsize=11, y=0.97)
    out = a.out or ("fig_impact_%s" % a.run.replace("out_", ""))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, out + "." + ext), dpi=190,
                    bbox_inches="tight")
    plt.close(fig)
    print("ecrit : %s | %d rompus, %d endommages | R_crat %.1f mm, "
          "portee max %.1f mm" % (out, nb, nd, rcr, port.max()))


if __name__ == "__main__":
    main()
