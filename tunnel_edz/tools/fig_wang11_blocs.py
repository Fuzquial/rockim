#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_wang11_blocs.py — les zooms de la fig. 11 de Wang et al. (2024) en
# RENDU ELEMENTS : chaque triangle est dessine a sa position COURANTE
# (deformee), en bloc rouge plein ; les fissures ouvertes et les blocs qui
# se decollent apparaissent en blanc (le fond passe a travers), comme dans
# leurs panneaux du bas. Complement de fig_wang11.py (champ U lisse).
#
#   python tunnel_edz/tools/fig_wang11_blocs.py out_tun_ref_stab \
#          --stem tunnel_edz/fig_wang11_blocs
#
# Modes de ruine et emplacements d'apres leur §4 : extrusion (squeezing) a la
# voute, compression-cisaillement au rein, soulevement au radier.
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})

CX = CY = 50.0
ZOOMS = ((-3.5, 3.5, 1.5, 8.5, "Extrusion en voûte (squeezing)"),
         (1.0, 8.0, -3.5, 3.5, "Compression–cisaillement au rein"),
         (-3.5, 3.5, -8.5, -1.5, "Soulèvement du radier"))


def read_tri_vtu(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    pts = np.fromstring(
        re.search(r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S).group(1), sep=" ").reshape(-1, 3)[:, :2]
    con = np.fromstring(
        re.search(r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S).group(1), sep=" ").astype(int)
    return pts, con.reshape(-1, 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_wang11_blocs")
    ap.add_argument("--frame", type=int, default=-1,
                    help="indice de frame (defaut : la derniere)")
    a = ap.parse_args()

    frames = [f for f in sorted(glob.glob(a.run + "/fdem_[0-9]*.vtu"))
              if "joints" not in os.path.basename(f)]
    f = frames[a.frame]
    t = None
    fcsv = os.path.join(a.run, "frames.csv")
    if os.path.exists(fcsv):
        rows = open(fcsv).read().splitlines()[1:]
        idx = a.frame if a.frame >= 0 else len(frames) + a.frame
        for r in rows:
            p = r.split(",")
            if int(p[0]) == idx:
                t = float(p[1])
    P, tri = read_tri_vtu(f)
    xy = P - [CX, CY]
    ctr = xy[tri].mean(axis=1)

    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.8))
    for A, (x0, x1, y0, y1, ttl) in zip(ax, ZOOMS):
        m = ((ctr[:, 0] > x0 - 1) & (ctr[:, 0] < x1 + 1) &
             (ctr[:, 1] > y0 - 1) & (ctr[:, 1] < y1 + 1))
        A.add_collection(PolyCollection(
            xy[tri[m]], facecolors="#cd3a2b", edgecolors="#7e1f14",
            linewidths=0.06))
        A.set_xlim(x0, x1)
        A.set_ylim(y0, y1)
        A.set_aspect("equal")
        A.set_title(ttl, fontsize=11)
        A.set_xticks([])
        A.set_yticks([])
        A.set_facecolor("white")
    lbl = "t = %.3f s" % t if t is not None else os.path.basename(f)
    fig.suptitle("Modes de ruine en rendu éléments — leurs zooms de "
                 "fig. 11 (%s ; blanc = fissures ouvertes et vides)" % lbl,
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=170)
    print("écrit : %s (frame %s, %s)" % (a.stem, a.frame, lbl))


if __name__ == "__main__":
    main()
