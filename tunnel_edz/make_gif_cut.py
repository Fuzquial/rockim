#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_gif_cut.py — animation d'un essai de coupe : maillage DEFORME, cutter a
# sa position courante, fissures. Pensee pour voir CE QUE FAIT L'OUTIL, y
# compris quand le run part en vrille.
#
#   python tunnel_edz/make_gif_cut.py out_cut_heilman [--depth 0.001016]
#          [--rake 20] [--x0 0 --x1 0.014 --y0 0.014 --y1 0.021]
#
# Deux precautions pour rester lisible sur un run divergent :
#   * la fenetre est FIXE et cadree sur la zone de coupe ;
#   * les triangles dont un sommet a fui au-dela de la boite d'affichage
#     elargie sont ECARTES du trace — sans quoi un seul noeud ejecte a
#     2 000 km tire des traits en travers de toute l'image.
# Le compteur d'elements ecartes est affiche : c'est lui qui dit quand le
# calcul cesse d'etre physique.
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

C_TEN, C_SHR = "#1B8A3A", "#C8342B"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_cut_heilman")
    ap.add_argument("--depth", type=float, default=0.001016)
    ap.add_argument("--rake", type=float, default=20.0)
    ap.add_argument("--face", type=float, default=0.013)
    ap.add_argument("--x0", type=float, default=0.0)
    ap.add_argument("--x1", type=float, default=0.014)
    ap.add_argument("--y0", type=float, default=0.0145)
    ap.add_argument("--y1", type=float, default=0.0215)
    ap.add_argument("--duree", type=int, default=650)
    ap.add_argument("--dpi", type=int, default=105)
    # Run EN COURS : history.csv s arrete au temps courant, pas a T. Sans ces
    # deux options le script deduirait tf de la derniere ligne et placerait le
    # cutter jusqu a 8 % trop en arriere — une trace qui ment.
    ap.add_argument("--T", type=float, default=0.0, help="duree physique visee")
    ap.add_argument("--frames", type=int, default=0, help="nb de trames visees")
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)
    out = os.path.join(HERE, os.path.basename(run) + "_coupe.gif")

    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    jn = [f for f in sorted(glob.glob(os.path.join(run, "fdem_joints_[0-9]*.vtu")))
          if complete(f)]
    n = min(len(el), len(jn))
    P0, C, _ = read_vtu(el[0], [])
    W, H = P0[:, 0].max(), P0[:, 1].max()
    d = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",",
                      names=True, invalid_raise=False)
    b = np.radians(a.rake)

    imgs = []
    for k in range(n):
        P, _, S = read_vtu(el[k], ["velocity"])
        # triangles encore dans la boite elargie : les autres ont fui
        box = ((P[:, 0] > -0.5 * W) & (P[:, 0] < 1.5 * W) &
               (P[:, 1] > -0.5 * H) & (P[:, 1] < 1.5 * H))
        keep = box[C].all(axis=1)
        nfly = int((~keep).sum())
        PJ, CJ, SJ = read_vtu(jn[k], ["breakMode"], ncell=2)
        bm = SJ["breakMode"]
        seg = PJ[CJ[bm > 0]]
        segok = np.array([s for s in seg
                          if np.all(np.abs(s[:, 0] - 0.5 * W) < W)
                          and np.all(np.abs(s[:, 1] - 0.5 * H) < H)]) \
            if len(seg) else np.zeros((0, 2, 2))

        tf = (a.T * k / a.frames) if (a.T > 0.0 and a.frames > 0)             else d["t"][-1] * k / max(n - 1, 1)
        xc = np.interp(tf, d["t"], d["toolX"])
        fx = np.interp(tf, d["t"], np.abs(d["toolFx"])) / 1e6
        edge = np.array([xc, H - a.depth])
        face = edge + a.face * np.array([-np.sin(b), np.cos(b)])

        fig, ax = plt.subplots(figsize=(11, 5.6))
        if keep.any():
            ax.triplot(mtri.Triangulation(P[:, 0], P[:, 1], C[keep]),
                       lw=0.25, color="0.55")
        if len(segok):
            ax.add_collection(LineCollection(segok, colors=C_SHR, lw=1.0))
        ax.plot([edge[0], face[0]], [edge[1], face[1]], color="#0B4F9E", lw=3)
        ax.plot(edge[0], edge[1], "o", color="#0B4F9E", ms=7)
        ax.set_aspect("equal")
        ax.set_xlim(a.x0, a.x1)
        ax.set_ylim(a.y0, a.y1)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(f"trame {k}/{n-1}   t = {tf*1e6:.0f} us   arete a "
                     f"x = {xc*1e3:.2f} mm   |Fx| = {fx:.2f} MN/m   "
                     f"elements enfuis : {nfly}", fontsize=11)
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=a.dpi)
        plt.close(fig)
        buf.seek(0)
        imgs.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))
        print(f"  trame {k}: arete {xc*1e3:6.2f} mm, {nfly:5d} elements enfuis, "
              f"{len(segok):5d} fissures")

    dur = [a.duree] * (len(imgs) - 1) + [3 * a.duree]
    imgs[0].save(out, save_all=True, append_images=imgs[1:], duration=dur,
                 loop=0, optimize=True)
    print(f"ecrit : {out}  ({len(imgs)} images, "
          f"{os.path.getsize(out)/1e6:.1f} Mo)")


if __name__ == "__main__":
    main()
