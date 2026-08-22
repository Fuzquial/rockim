#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_wang11.py — reproduction de la fig. 11 de Wang et al. (2024) :
# champ de deplacement U autour de la cavite (configuration DEFORMEE) et
# trois zooms sur les modes de ruine — serrage aux reins, rupture de
# compression-cisaillement, soulevement du radier.
#
#   python tunnel_edz/tools/fig_wang11.py out_tun_ref_iso \
#          --stem tunnel_edz/fig_wang11
#
# Les VTU de rockim portent les positions COURANTES : U = P(t_fin) - P(0),
# exactement la methode d'edz_metrics. Echelle arc-en-ciel comme l'article
# (leur bornage 0 - 0,37 m).
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})

CX = CY = 50.0


def read_tri_vtu(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    pts = np.fromstring(
        re.search(r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S).group(1), sep=" ").reshape(-1, 3)[:, :2]
    con = np.fromstring(
        re.search(r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S).group(1), sep=" ").astype(int)
    ncell = con.size // 3
    return pts, con.reshape(ncell, 3)


def read_lines_vtu(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    pts = np.fromstring(
        re.search(r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S).group(1), sep=" ").reshape(-1, 3)[:, :2]
    con = np.fromstring(
        re.search(r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S).group(1), sep=" ").astype(int).reshape(-1, 2)
    dmg = np.fromstring(
        re.search(r'Name="damage"[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S).group(1), sep=" ")
    return pts, con, dmg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_wang11")
    ap.add_argument("--umax", type=float, default=0.0,
                    help="borne de l'echelle [m] (0 = max mesure)")
    a = ap.parse_args()

    frames = sorted(glob.glob(a.run + "/fdem_[0-9]*.vtu"))
    P0, _ = read_tri_vtu(frames[0])
    P1, tri = read_tri_vtu(frames[-1])
    U = np.linalg.norm(P1 - P0, axis=1)
    umax = a.umax or float(U.max())
    x, y = P1[:, 0] - CX, P1[:, 1] - CY
    jf = sorted(glob.glob(a.run + "/fdem_joints_[0-9]*.vtu"))
    jp, jc, jd = read_lines_vtu(jf[-1])
    seg = (jp - [CX, CY])[jc[jd >= 0.999]]

    fig = plt.figure(figsize=(12.6, 9.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[2.1, 1.0], hspace=0.28,
                          wspace=0.25)
    A = fig.add_subplot(gs[0, :])
    tp = A.tripcolor(x, y, tri, U, cmap="jet", vmin=0, vmax=umax,
                     shading="gouraud", rasterized=True)
    A.add_collection(LineCollection(seg, colors="white", linewidths=0.25,
                                    alpha=0.55))
    cb = fig.colorbar(tp, ax=A, pad=0.015)
    cb.set_label("U [m]")
    A.set_xlim(-28, 28)
    A.set_ylim(-24, 24)
    A.set_aspect("equal")
    A.set_xlabel("x [m]")
    A.set_ylabel("y [m]")
    A.set_title("Champ de déplacement (configuration déformée) et réseau "
                "de fissures", fontsize=12)
    A.annotate("clé de voûte", (0, 6.2), ha="center", fontsize=9,
               color="white")
    A.annotate("rein", (7.2, 0), fontsize=9, color="white")
    A.annotate("radier", (0, -6.6), ha="center", fontsize=9, color="white")

    zooms = ((-9.5, -0.5, -4.0, 5.0, "Serrage aux reins (squeezing)"),
             (0.5, 9.5, -4.0, 5.0, "Compression–cisaillement"),
             (-4.5, 4.5, -9.5, -0.5, "Soulèvement du radier"))
    for j, (x0, x1, y0, y1, ttl) in enumerate(zooms):
        B = fig.add_subplot(gs[1, j])
        B.tripcolor(x, y, tri, U, cmap="jet", vmin=0, vmax=umax,
                    shading="gouraud", rasterized=True)
        B.add_collection(LineCollection(seg, colors="white",
                                        linewidths=0.45, alpha=0.75))
        B.set_xlim(x0, x1)
        B.set_ylim(y0, y1)
        B.set_aspect("equal")
        B.set_title(ttl, fontsize=10)
        B.set_xticks([])
        B.set_yticks([])

    fig.suptitle("Tunnel de Hutou Beishan — leur fig. 11 "
                 "(σ₀ = 5 MPa, λ = 1) : U max = %.3f m (publié : 0,347)"
                 % U.max(), fontsize=13)
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)
    print("écrit : %s | U max %.3f m" % (a.stem, U.max()))


if __name__ == "__main__":
    main()
