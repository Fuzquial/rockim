#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_montage_impact.py — MONTRER LE MONTAGE d'impact avant de le lancer.
#
#   python tools/fig_montage_impact.py meshes/impact_yang_s1.msh \
#          [--out results/montage.png] [--v0 9.0]
#
# Complement de fig_mesh3d.py, qui montre la GRADATION du maillage. Ici on
# montre l'ASSEMBLAGE : quel corps est ou, quels jeux les separent, et en
# combien de temps l'onde va de l'un a l'autre. C'est ce qui manquait pour
# comprendre pourquoi un run court ne montre rien.
#
# Les trois panneaux :
#   (a) la chaine complete, chaque corps d'une couleur, cotee ;
#   (b) le zoom sur l'insert et la surface de la roche, le seul endroit ou la
#       physique interessante se passe, avec le jeu initial ;
#   (c) la CHRONOLOGIE : fermeture du jeu piston-bit a v0, traversee du bit a
#       la celerite de l'acier, arrivee sur la roche. En nombre de PAS de
#       temps aussi, parce que c'est ce qui se paie.
#
# La celerite prise pour le bit est celle d'une onde de barre, c = sqrt(E/rho),
# avec les proprietes acier de la Table 1 de Yang et al. 2026 (200 GPa, 7850).
# C'est la bonne pour un barreau mince devant la longueur d'onde ; elle
# sous-estime de quelques pour cent la celerite 3D de dilatation.
#
# Post-traitement PUR : lit un MSH 2.2 ASCII, n'ecrit qu'un PNG.
# ---------------------------------------------------------------------------
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
COL = {"rock": "#c8b89a", "insert": "#c0392b", "bit": "#4a6fa5",
       "piston": "#2e4a6b", "circlip": "#7b9e3f", "plate": "#a5824a"}


def read_msh22(path):
    names, nodes, tets = {}, {}, {}
    with open(path) as f:
        it = iter(f)
        for line in it:
            L = line.strip()
            if L == "$PhysicalNames":
                for _ in range(int(next(it))):
                    p = next(it).split()
                    names[int(p[1])] = p[2].strip('"')
            elif L == "$Nodes":
                for _ in range(int(next(it))):
                    p = next(it).split()
                    nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            elif L == "$Elements":
                for _ in range(int(next(it))):
                    p = next(it).split()
                    if int(p[1]) != 4:
                        continue
                    nt = int(p[2])
                    tets.setdefault(int(p[3]), []).append(
                        [int(x) for x in p[3 + nt:3 + nt + 4]])
    ids = sorted(nodes)
    idx = {k: i for i, k in enumerate(ids)}
    P = np.array([nodes[k] for k in ids])
    T = {names[ph]: np.array([[idx[a] for a in t] for t in v], dtype=int)
         for ph, v in tets.items()}
    return P, T


def slice_polys(P, T, axis, value):
    """Section exacte des tetraedres par le plan axis = value."""
    s = P[:, axis] - value
    keep = [i for i in range(3) if i != axis]
    out = []
    sv = s[T]
    cross = (sv.min(1) < 0) & (sv.max(1) > 0)
    for t, sgn in zip(T[cross], sv[cross]):
        pts = []
        for a, b in EDGES:
            if (sgn[a] < 0) != (sgn[b] < 0):
                w = sgn[a] / (sgn[a] - sgn[b])
                pts.append(P[t[a]] + w * (P[t[b]] - P[t[a]]))
        if len(pts) < 3:
            continue
        q = np.array(pts)[:, keep]
        if len(q) == 4:
            c = q.mean(0)
            q = q[np.argsort(np.arctan2(q[:, 1] - c[1], q[:, 0] - c[0]))]
        out.append(q)
    return out


def main():
    ap = argparse.ArgumentParser(description="montrer le montage d'impact")
    ap.add_argument("msh")
    ap.add_argument("--out", default="results/montage_impact.png")
    ap.add_argument("--v0", type=float, default=9.0, help="vitesse piston m/s")
    ap.add_argument("--dt", type=float, default=9.42436e-10,
                    help="pas de temps du solveur [s]")
    ap.add_argument("--spp", type=float, default=0.125,
                    help="cout mesure par pas [s]")
    a = ap.parse_args()

    P, T = read_msh22(a.msh)
    polys = {k: slice_polys(P, v, 1, 0.0) for k, v in T.items()}
    zr = {k: (P[np.unique(v)][:, 2].min(), P[np.unique(v)][:, 2].max())
          for k, v in T.items()}

    fig = plt.figure(figsize=(14.5, 8.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.35])
    A, B, C = (fig.add_subplot(gs[0, i]) for i in range(3))

    for ax, lw in ((A, 0.08), (B, 0.5)):
        for k, pl in polys.items():
            ax.add_collection(PolyCollection(
                pl, facecolors=COL.get(k, "0.6"), edgecolors="k",
                linewidths=lw, alpha=0.95))

    # ---- (a) la chaine ---------------------------------------------------
    A.set_xlim(-0.14, 0.14)
    A.set_ylim(-0.16, 0.54)
    A.set_aspect("equal")
    A.set_xlabel("x (m)")
    A.set_ylabel("z (m)")
    A.set_title("(a) La chaine, six corps")
    for k, (z0, z1) in zr.items():
        A.annotate("%s" % k, (0.10, 0.5 * (z0 + z1)), fontsize=9,
                   color=COL.get(k, "k"), fontweight="bold", ha="left")
    A.axhline(0.0, color="crimson", lw=0.8, ls="--")

    # ---- (b) le zoom sur le contact ------------------------------------
    B.set_xlim(-0.016, 0.016)
    B.set_ylim(-0.016, 0.026)
    B.set_aspect("equal")
    B.set_xlabel("x (m)")
    B.set_title("(b) Insert carbure R 8,51 mm sur la roche")
    gap = (zr["insert"][0] - zr["rock"][1]) * 1e3
    B.annotate("jeu initial %.2f mm" % gap, (0.0, 0.0),
               xytext=(0.007, -0.006), fontsize=9, color="crimson",
               arrowprops=dict(arrowstyle="->", color="crimson", lw=1.2))
    B.axhline(0.0, color="crimson", lw=0.8, ls="--")

    # ---- (c) la chronologie ---------------------------------------------
    c = np.sqrt(200e9 / 7850.0)
    gpb = zr["piston"][0] - zr["bit"][1]
    Lbit = zr["bit"][1] - zr["bit"][0]
    t1 = gpb / a.v0
    t2 = t1 + Lbit / c
    ev = [("le piston part\n(jeu %.2f mm)" % (gpb * 1e3), 0.0),
          ("il touche le bit", t1),
          ("l'onde sort du bit\net charge la roche", t2),
          ("pic de force\n(leur fig. 9a)", t2 + 30e-6),
          ("rebond du bit\n(leur fig. 9d)", t2 + 400e-6)]
    # Les trois premiers evenements tiennent dans 70 us sur une echelle de
    # 800 : empiles a la verticale ils se recouvrent. Chaque etiquette est
    # donc posee sur sa PROPRE ligne, reliee a son instant par un trait.
    yb = 0.0
    C.barh(yb, t2 * 1e6, color="0.85", edgecolor="k", height=0.42,
           label="preambule : RIEN n'atteint la roche")
    C.barh(yb, (800e-6 - t2) * 1e6, left=t2 * 1e6, color="#c0392b",
           alpha=0.35, edgecolor="k", height=0.42, label="physique utile")
    for i, (lab, t) in enumerate(ev):
        y = 0.75 + 0.62 * i
        x = t * 1e6
        C.plot([x, x], [yb + 0.21, y - 0.06], color="0.45", lw=0.9)
        C.plot([x], [y - 0.06], "o", color="0.3", ms=3)
        C.annotate("%s  —  %.0f us,  %s pas,  %.1f h"
                   % (lab.replace("\n", " "), x, format(int(t / a.dt), ","),
                      t / a.dt * a.spp / 3600.0),
                   (x, y), fontsize=8.5, ha="left" if x < 450 else "right",
                   va="bottom")
    C.set_xlim(-50, 860)
    C.set_ylim(-0.6, 4.3)
    C.set_yticks([])
    C.set_xlabel("temps simule (us)")
    C.set_title("(c) Chronologie, et ce qu'elle coute")
    C.legend(fontsize=8, loc="lower right")

    fig.suptitle("Montage d'impact Yang et al. 2026 — %s | %d tetraedres, "
                 "231 471 joints, dt = %.2f ns, %.3f s/pas mesure"
                 % (os.path.basename(a.msh), sum(len(v) for v in T.values()),
                    a.dt * 1e9, a.spp), fontsize=12)
    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=130)
    print("wrote " + a.out)
    print("  jeu piston-bit %.3f mm -> fermeture a %.1f us" % (gpb * 1e3, t1 * 1e6))
    print("  bit %.0f mm, c = %.0f m/s -> traversee %.1f us"
          % (Lbit * 1e3, c, Lbit / c * 1e6))
    print("  la roche est chargee a %.0f us = %d pas = %.1f h de calcul"
          % (t2 * 1e6, t2 / a.dt, t2 / a.dt * a.spp / 3600.0))


if __name__ == "__main__":
    main()
