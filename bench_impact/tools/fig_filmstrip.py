#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_filmstrip.py — la bande de trames de l'article (leur fig. 14) : une
# colonne par instant, vue de dessus au-dessus, coupe verticale en dessous.
#
#   python bench_impact/tools/fig_filmstrip.py out_imp_fidele 14 15 16 17 18 \
#          --stem bench_impact/fig_strip
#
# Rendu par FACES, avec la meme discrimination geometrique que fig_impact :
# sub-vertical = fissure (rouge traction / jaune cisaillement), sub-horizontal
# = zone broyee (rose pale). C'est ce qui fait apparaitre l'etoile radiale.
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
from imp_lib import (CX, CY, Z_SURF, broken, frame_times, joints_frame,
                     read_vtu)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})
ROUGE, JAUNE = "#b22222", "#e6a817"


def faces(AX, P, mode, n, a0, a1, c0, c1):
    vert = np.abs(n[:, 2]) <= 0.6
    for mm, col, al, z in ((~vert, "#e8b7b7", 0.30, 1),
                           (vert & (mode < 1.5), ROUGE, 0.85, 3),
                           (vert & (mode >= 1.5), JAUNE, 0.9, 4)):
        if mm.any():
            po = (P[mm][:, :, (a0, a1)] - np.array([c0, c1])) * 1e3
            AX.add_collection(PolyCollection(
                po, facecolors=col, edgecolors=col, linewidths=0.2,
                alpha=al, zorder=z))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("frames", nargs="+", type=int)
    ap.add_argument("--stem", default="fig_strip")
    ap.add_argument("--half", type=float, default=35.0, help="demi-fenetre mm")
    a = ap.parse_args()

    tf = frame_times(a.run)
    N = len(a.frames)
    fig, ax = plt.subplots(2, N, figsize=(2.9 * N, 6.4))
    if N == 1:
        ax = ax.reshape(2, 1)
    fig.suptitle("Impact St Anne 10,66 m/s — propagation des fissures "
                 "(run « fidèle », roche au 1 mm)", fontsize=12)

    for j, k in enumerate(a.frames):
        pts, con, f = read_vtu(joints_frame(a.run, k))
        c, n, mode, P = broken(pts, con, f)
        tk = tf.get(k, 0.0) * 1e6
        A, B = ax[0, j], ax[1, j]
        if len(c):
            faces(A, P, mode, n, 0, 1, CX, CY)
            s5 = np.abs(c[:, 1] - CY) < 0.005
            faces(B, P[s5], mode[s5], n[s5], 0, 2, CX, Z_SURF)
        for X in (A, B):
            X.set_xlim(-a.half, a.half)
            X.set_aspect("equal")
            X.set_xticks([-30, 0, 30])
        A.set_ylim(-a.half, a.half)
        A.set_yticks([-30, 0, 30])
        B.set_ylim(-a.half, 5)
        B.set_yticks([-30, -15, 0])
        B.axhline(0, color="#333", lw=0.7)
        A.set_title("%.0f $\\mu$s   (%d rompus)" % (tk, len(c)), fontsize=10)
        if j:
            A.set_yticklabels([])
            B.set_yticklabels([])
        else:
            A.set_ylabel("vue de dessus\ny [mm]")
            B.set_ylabel("coupe $|y|<5$ mm\nz [mm]")
        B.set_xlabel("x [mm]")
        print("  trame %2d : t = %6.1f us, %5d joints rompus" % (k, tk, len(c)))

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)
    print("écrit : %s.pdf et .png" % a.stem)


if __name__ == "__main__":
    main()
