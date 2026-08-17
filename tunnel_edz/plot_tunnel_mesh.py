#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Figure de CONTROLE du maillage du tunnel en fer a cheval (Wang et al. 2024,
# Front. Earth Sci. 12:1517816) genere par
#   tools/make_unstructured_mesh.py tunnelhs 100 100 0.22 18 2.0 \
#                                   meshes/tunnel_hs.msh 1
#
# Trois panneaux :
#   (a) le massif entier avec ses trois zones de taille (leur fig. 6b) ;
#   (b) un zoom de 30 m sur la cavite (la zone qui portera l'EDZ) ;
#   (c) le PROFIL seul, cote, a comparer trait pour trait a leur fig. 6a —
#       c'est le panneau qui sert a valider (ou corriger) ma lecture du dessin.
#
# Lancement sans argument depuis ce dossier :  python plot_tunnel_mesh.py
# ---------------------------------------------------------------------------
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from make_unstructured_mesh import TUNNEL_HS, tunnel_hs_profile  # noqa: E402

MESH = os.path.join(HERE, "..", "meshes", "tunnel_hs.msh")
OUT = os.path.join(HERE, "tunnel_mesh_check.png")


def read_msh22(path):
    """Noeuds + triangles d'un MSH 2.2 ASCII (le format que rockim relit)."""
    xy, tri = [], []
    idx = {}
    with open(path) as f:
        it = iter(f)
        for line in it:
            if line.startswith("$Nodes"):
                for _ in range(int(next(it))):
                    p = next(it).split()
                    idx[int(p[0])] = len(xy)
                    xy.append((float(p[1]), float(p[2])))
            elif line.startswith("$Elements"):
                for _ in range(int(next(it))):
                    p = next(it).split()
                    if int(p[1]) != 2:            # 2 = triangle
                        continue
                    n = int(p[2])
                    tri.append([idx[int(v)] for v in p[3 + n:6 + n]])
    return np.array(xy), np.array(tri, dtype=int)


def profile_xy(cx, cy0, n=400):
    """Contour analytique du tunnel (memes formules que le generateur)."""
    g = tunnel_hs_profile()
    hs, ht, ys, yw = g["halfSpan"], g["height"], g["ySpring"], g["yWall"]
    tx, ty = g["tang"]
    seg = []
    # radier : arc de centre (0, r4) du point bas a la tangence
    a0 = np.arctan2(0.0 - g["r4"], 0.0)
    a1 = np.arctan2(ty - g["r4"], tx)
    t = np.linspace(a0, a1, n // 4)
    seg.append(np.c_[g["r4"] * np.cos(t), g["r4"] + g["r4"] * np.sin(t)])
    # conge : arc de centre c3 de la tangence au pied du piedroit
    b0 = np.arctan2(ty - g["c3"][1], tx - g["c3"][0])
    b1 = np.arctan2(yw - g["c3"][1], hs - g["c3"][0])
    t = np.linspace(b0, b1, n // 4)
    seg.append(np.c_[g["c3"][0] + g["r3"] * np.cos(t),
                     g["c3"][1] + g["r3"] * np.sin(t)])
    # piedroit droit
    seg.append(np.c_[[hs, hs], [yw, ys]])
    # voute : demi-arc droit puis gauche
    t = np.linspace(0.0, np.pi, n)
    seg.append(np.c_[g["rCrown"] * np.cos(t), ys + g["rCrown"] * np.sin(t)])
    right = np.vstack(seg[:3])
    crown = seg[3]
    left = right[::-1] * np.array([-1.0, 1.0])
    P = np.vstack([right, crown, left])
    return P[:, 0] + cx, P[:, 1] + cy0, g


def main():
    xy, tri = read_msh22(MESH)
    T = mtri.Triangulation(xy[:, 0], xy[:, 1], tri)
    W, H = xy[:, 0].max(), xy[:, 1].max()
    cx, cy0 = 0.5 * W, 0.5 * H - 0.5 * TUNNEL_HS["height"]
    px, py, g = profile_xy(cx, cy0)

    fig = plt.figure(figsize=(15.5, 5.4))
    ax = fig.add_subplot(1, 3, 1)
    ax.triplot(T, lw=0.08, color="0.35")
    ax.plot(px, py, color="#1FA7B5", lw=1.4)
    ax.set_aspect("equal")
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_title(f"(a) massif {W:.0f} x {H:.0f} m — {len(tri)} triangles")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")

    ax = fig.add_subplot(1, 3, 2)
    ax.triplot(T, lw=0.25, color="0.35")
    ax.plot(px, py, color="#1FA7B5", lw=1.6)
    ax.set_aspect("equal")
    ax.set_xlim(cx - 15, cx + 15)
    ax.set_ylim(0.5 * H - 15, 0.5 * H + 15)
    ax.set_title("(b) zoom 30 m — la zone fine porte l'EDZ")
    ax.set_xlabel("x [m]")

    ax = fig.add_subplot(1, 3, 3)
    ax.plot(px - cx, py - cy0, color="#1FA7B5", lw=2.0)
    ax.set_aspect("equal")
    hs, ht, ys, yw = g["halfSpan"], g["height"], g["ySpring"], g["yWall"]
    # cotes de l'article, reportees sur la reconstruction
    def cote(x0, y0, x1, y1, txt, dx=0.0, dy=0.25):
        ax.annotate("", (x1, y1), (x0, y0),
                    arrowprops=dict(arrowstyle="<->", lw=0.8, color="0.25"))
        ax.text(0.5 * (x0 + x1) + dx, 0.5 * (y0 + y1) + dy, txt, ha="center",
                va="bottom", fontsize=8, color="0.15")
    cote(-hs, -1.1, hs, -1.1, f"{2*hs:.2f} m  (cote 11 m)")
    cote(hs + 1.6, 0.0, hs + 1.6, ht, f"{ht:.2f} m", dx=0.7, dy=-0.2)
    cote(-hs - 1.0, ys, -hs - 1.0, ht, f"{ht-ys:.2f} m", dx=-1.5, dy=-0.3)
    cote(-hs - 1.0, yw, -hs - 1.0, ys, f"{ys-yw:.2f}\n(cote 1,6)",
         dx=-1.5, dy=-0.6)
    cote(-hs - 1.0, 0.0, -hs - 1.0, yw, f"{yw:.2f} m", dx=-1.5, dy=-0.3)
    ax.plot([0], [ys], "k+", ms=8)
    ax.text(0.2, ys + 0.15, f"centre voute R{g['rCrown']:.2f}", fontsize=8)
    ax.plot(*g["tang"], "o", ms=3.5, color="#C8553D")
    ax.annotate(f"tangence radier/conge\nR{g['r4']:.1f} et R{g['r3']:.2f}",
                xy=g["tang"], xytext=(7.8, -1.7), fontsize=8, color="#C8553D",
                ha="center", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.7, color="#C8553D"))
    ax.set_xlim(-11, 11)
    ax.set_ylim(-2.2, 10.5)
    ax.set_title("(c) profil reconstruit — a comparer a leur fig. 6a")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")

    fig.tight_layout()
    fig.savefig(OUT, dpi=160)
    print("ecrit :", OUT)


if __name__ == "__main__":
    main()
