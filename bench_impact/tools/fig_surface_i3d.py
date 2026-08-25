#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_surface_i3d.py — l'endommagement VU DE HAUT, sur la seule face libre.
#
#   python bench_impact/tools/fig_surface_i3d.py out_imp3d_gros \
#          --stem bench_impact/fig_surf --demi 40
#
# On ne garde que les faces de tetraedres POSEES SUR LA SURFACE (leurs trois
# sommets a z = z_max) et on les colore par l'endommagement du tetraedre
# parent : c'est ce que verrait une camera au-dessus de l'eprouvette, sans le
# brouillard des elements profonds qu'une projection de volume superpose.
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})


def read(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    P = np.fromstring(re.search(
        r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").reshape(-1, 3)
    con = np.fromstring(re.search(
        r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").astype(int).reshape(-1, 4)
    m = re.search(r'Name="(?:dfhD|damage)"[^>]*>\s*(.*?)\s*</DataArray>',
                  s, re.S)
    D = np.fromstring(m.group(1), sep=" ") if m else None
    return P, con, D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_surf")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--demi", type=float, default=40.0, help="demi-fenetre mm")
    ap.add_argument("--tol", type=float, default=1e-4, help="tolerance z [m]")
    a = ap.parse_args()

    fs = [f for f in sorted(glob.glob(a.run + "/fem3d_[0-9]*.vtu"))
          + sorted(glob.glob(a.run + "/fdem3d_[0-9]*.vtu"))
          if "joints" not in f]
    P0, con, _ = read(fs[0])
    _, _, D = read(fs[a.frame])

    zs = P0[:, 2].max()
    haut = P0[:, 2] > zs - a.tol                    # noeuds de la face libre
    # une face de tetraedre est SUR la surface si 3 de ses 4 noeuds y sont
    faces, val = [], []
    combos = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    for e in range(len(con)):
        n4 = con[e]
        if haut[n4].sum() < 3:
            continue
        for c in combos:
            f3 = n4[list(c)]
            if haut[f3].all():
                faces.append(f3)
                val.append(D[e])
                break
    faces = np.array(faces)
    val = np.array(val)

    cx = 0.5 * (P0[:, 0].min() + P0[:, 0].max())
    cy = 0.5 * (P0[:, 1].min() + P0[:, 1].max())
    x = (P0[:, 0] - cx) * 1e3
    y = (P0[:, 1] - cy) * 1e3

    fig, A = plt.subplots(figsize=(6.9, 6.2))
    tp = A.tripcolor(Triangulation(x, y, faces), facecolors=val,
                     cmap="inferno", vmin=0, vmax=1, rasterized=True)
    cb = fig.colorbar(tp, ax=A, pad=0.02)
    cb.set_label(r"$D_{\max}$ (DP-DFH) en surface")
    th = np.linspace(0, 2 * np.pi, 200)
    for R in (10, 20, 30):
        A.plot(R * np.cos(th), R * np.sin(th), color="w", lw=0.6, alpha=0.5)
    A.set_xlim(-a.demi, a.demi)
    A.set_ylim(-a.demi, a.demi)
    A.set_aspect("equal")
    A.set_xlabel("x [mm]")
    A.set_ylabel("y [mm]")
    n = (val > 0.9).sum()
    A.set_title("Impact 3D DP-DFH — endommagement de la FACE LIBRE\n"
                "%d faces de surface, %d a D > 0,9 (cercles blancs : "
                "10, 20, 30 mm)" % (len(val), n), fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)
    print("écrit : %s | %d faces de surface, %d a D > 0,9"
          % (a.stem, len(val), n))


if __name__ == "__main__":
    main()
