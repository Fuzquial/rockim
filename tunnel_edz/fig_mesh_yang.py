#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_mesh_yang.py — controle visuel du maillage equivalent au banc de Yang
# et al. (IJRMMS 191, 2025), dans l esprit de leur fig. 6 : coupe axiale
# graduee, vue de dessus de la surface, zoom sur la zone de contact.
#
#   python tunnel_edz/fig_mesh_yang.py [out_yang_equiv] [--slab 1.5]
#
# Lit la trame 0 du run (geometrie non deformee). Trace les ARETES des
# tetraedres de la tranche, pas leurs centroides : c est la seule facon de
# voir si la gradation est reguliere ou si elle saute.
# ---------------------------------------------------------------------------
import argparse
import itertools
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.labelsize": 9.5, "axes.titlesize": 10,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
C_ROCK, C_INS = "#6E6A66", "#2E5E8C"


def read(path):
    d = open(path, "r", errors="ignore").read()
    g = lambda p: re.search(p, d, re.S).group(1)
    P = np.fromstring(g(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>"),
                      sep=" ").reshape(-1, 3)
    C = np.fromstring(g(r'Name="connectivity"[^>]*>(.*?)</DataArray>'),
                      sep=" ", dtype=np.int64).reshape(-1, 4)
    ph = np.fromstring(g(r'Name="phase"[^>]*>(.*?)</DataArray>'), sep=" ")
    return P, C, ph


def edges(T, ax_i, ax_j):
    """6 aretes de chaque tetra, projetees sur deux axes."""
    seg = []
    for i, j in itertools.combinations(range(4), 2):
        seg.append(np.stack([T[:, i, [ax_i, ax_j]], T[:, j, [ax_i, ax_j]]], 1))
    return np.concatenate(seg, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_yang_equiv")
    ap.add_argument("--slab", type=float, default=1.5, help="demi-tranche [mm]")
    ap.add_argument("--half", type=float, default=125.0, help="demi-fenetre [mm]")
    a = ap.parse_args()

    P, C, ph = read(os.path.join(ROOT, a.run, "fdem3d_0000.vtu"))
    T = P[C]
    cen = T.mean(axis=1)
    ins = ph > 0.5                                  # phase 1 = insert
    x0, y0 = cen[ins, 0].mean(), cen[ins, 1].mean()
    # phase est un champ de CELLULE : indexer P par ins serait faux. On prend
    # les noeuds des tetras de ROCHE (erreur payee le 2026-08-19).
    zs = P[np.unique(C[~ins].ravel())][:, 2].max()
    S = a.slab * 1e-3

    fig, ax = plt.subplots(1, 3, figsize=(12.6, 4.2))

    # ---- (a) coupe axiale, vue large ---------------------------------
    for p, half, ttl in ((ax[0], a.half, "(a) coupe axiale — bloc entier"),
                         (ax[2], 0.22 * a.half, "(c) zoom sur la zone de contact")):
        H = half * 1e-3
        sl = (np.abs(cen[:, 1] - y0) < S) & (np.abs(cen[:, 0] - x0) < H) \
             & (cen[:, 2] > zs - 2.2 * H) & (cen[:, 2] < zs + 0.35 * H)
        for m, c, lw in ((sl & ~ins, C_ROCK, 0.18), (sl & ins, C_INS, 0.30)):
            if not m.any():
                continue
            e = edges(T[m], 0, 2)
            e[:, :, 0] = (e[:, :, 0] - x0) * 1e3
            e[:, :, 1] = (e[:, :, 1] - zs) * 1e3
            p.add_collection(LineCollection(e, colors=c, lw=lw, alpha=0.85))
        p.set_xlim(-half, half)
        p.set_ylim(-2.2 * half, 0.35 * half)
        p.set_aspect("equal")
        p.set_xlabel("x [mm]")
        p.set_ylabel("z − z$_{surface}$ [mm]")
        p.set_title(ttl, fontsize=10)

    # ---- (b) vue de dessus de la surface -----------------------------
    p = ax[1]
    near = (~ins) & (cen[:, 2] > zs - 3e-3)
    e = edges(T[near], 0, 1)
    e[:, :, 0] = (e[:, :, 0] - x0) * 1e3
    e[:, :, 1] = (e[:, :, 1] - y0) * 1e3
    p.add_collection(LineCollection(e, colors=C_ROCK, lw=0.12, alpha=0.8))
    p.set_xlim(-a.half, a.half)
    p.set_ylim(-a.half, a.half)
    p.set_aspect("equal")
    p.set_xlabel("x [mm]")
    p.set_ylabel("y [mm]")
    p.set_title("(b) surface vue de dessus (3 mm sous la peau)", fontsize=10)

    for p in ax:
        for s in ("top", "right"):
            p.spines[s].set_visible(False)

    # ---- statistiques de gradation, imprimees ------------------------
    ed = np.zeros(len(T))
    for i, j in itertools.combinations(range(4), 2):
        ed += np.linalg.norm(T[:, i] - T[:, j], axis=1)
    ed /= 6.0
    r = np.hypot(cen[:, 0] - x0, cen[:, 1] - y0)
    lig = []
    for lo, hi in [(0, 12.5), (12.5, 25), (25, 50), (50, 125)]:
        m = (~ins) & (r >= lo * 1e-3) & (r < hi * 1e-3)
        if m.sum():
            lig.append("r %3d-%3d mm : %6d tet, arete med %5.2f mm"
                       % (lo, hi, m.sum(), np.median(ed[m]) * 1e3))
    lig.append("insert        : %6d tet, arete med %5.2f mm"
               % (ins.sum(), np.median(ed[ins]) * 1e3))
    print("\n".join("  " + l for l in lig))

    fig.suptitle("Maillage équivalent au banc de Yang — %d tétraèdres "
                 "(roche %d, insert %d)"
                 % (len(T), int((~ins).sum()), int(ins.sum())),
                 fontsize=11, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_mesh_yang." + ext), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print("  ecrit : fig_mesh_yang.pdf / .png")


if __name__ == "__main__":
    main()
