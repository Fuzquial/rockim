#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_mesh3d.py — MONTRER un maillage tetraedrique AVANT de lancer.
#
#   python tools/fig_mesh3d.py meshes/impact3d_yang.msh [--out fig.png]
#                              [--zoom R] [--plane y]
#
# POURQUOI CE SCRIPT EXISTE. Regle du 2026-09-07 : 60 runs ont ete rendus
# ininterpretables parce que le maillage n'avait pas ete regarde avant
# (5,5 elements par grain la ou il en fallait 35 a 90). Un maillage se VOIT ;
# un compte d'elements ne dit rien de la gradation ni de la qualite.
#
# CE QU'IL TRACE.
#   (a) COUPE EXACTE par le plan median. Chaque tetraedre traverse par le plan
#       est reellement DECOUPE (triangle ou quadrilatere selon le nombre
#       d'aretes coupees), pas approxime par une tranche d'elements dont le
#       centre tombe pres du plan — une tranche d'epaisseur fixe melange les
#       tailles et donne une image fausse de la gradation.
#   (b) ZOOM sur la zone d'impact, ou se joue tout le facies.
#   (c) PROFIL de la taille de maille realisee en fonction de la distance au
#       point d'impact, avec les paliers vises en pointilles. C'est le
#       controle falsifiant : un champ de taille peut etre correctement ecrit
#       et mal realise (Gmsh lisse, et MeshSizeExtendFromBoundary peut
#       ecraser le champ). On verifie ce qui est SORTI, pas ce qui a ete
#       demande.
#   (d) QUALITE : distribution du rapport de forme (rayon insphere sur
#       circonsphere, normalise a 1 pour le tetraedre regulier). En explicite,
#       un seul element degenere fixe le pas de temps de tout le calcul.
#
# Post-traitement PUR : lit un MSH 2.2 ASCII, n'ecrit qu'un PNG.
# ---------------------------------------------------------------------------
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection


def read_msh22(path):
    """Noeuds et tetraedres d'un MSH 2.2 ASCII."""
    nodes, tets = {}, []
    with open(path) as f:
        it = iter(f)
        for line in it:
            line = line.strip()
            if line == "$Nodes":
                n = int(next(it))
                for _ in range(n):
                    p = next(it).split()
                    nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            elif line == "$Elements":
                m = int(next(it))
                for _ in range(m):
                    p = next(it).split()
                    if int(p[1]) != 4:            # 4 = tetraedre
                        continue
                    nt = int(p[2])
                    tets.append([int(x) for x in p[3 + nt:3 + nt + 4]])
    ids = sorted(nodes)
    idx = {k: i for i, k in enumerate(ids)}
    P = np.array([nodes[k] for k in ids])
    T = np.array([[idx[a] for a in t] for t in tets], dtype=int)
    return P, T


EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def slice_tets(P, T, axis, value):
    """Section EXACTE des tetraedres par le plan axis = value.

    Chaque tetraedre dont les sommets ne sont pas tous du meme cote est
    decoupe : les aretes changeant de signe donnent 3 ou 4 points, qui
    forment la face de section. Les quadrilateres sont reordonnes par angle
    autour de leur centroide (l'ordre des aretes ne donne pas un polygone
    simple).
    """
    s = P[:, axis] - value
    keep = [i for i in range(3) if i != axis]
    polys = []
    sv = s[T]                                     # (nTet, 4)
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
        polys.append(q)
    return polys


def tet_quality(P, T):
    """Rapport insphere/circonsphere x 3, egal a 1 pour le tetraedre regulier.

    r_in = 3V/A_tot ; R_circ par la formule du determinant. Le facteur 3
    normalise : pour un tetraedre regulier R = 3 r.
    """
    a, b, c, d = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]], P[T[:, 3]]
    V = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0
    def area(p, q, r):
        return 0.5 * np.linalg.norm(np.cross(q - p, r - p), axis=1)
    A = area(a, b, c) + area(a, b, d) + area(a, c, d) + area(b, c, d)
    rin = 3.0 * V / np.maximum(A, 1e-300)
    la, lb, lc = b - a, c - a, d - a
    num = np.linalg.norm(
        np.linalg.norm(la, axis=1, keepdims=True) ** 2 * np.cross(lb, lc)
        + np.linalg.norm(lb, axis=1, keepdims=True) ** 2 * np.cross(lc, la)
        + np.linalg.norm(lc, axis=1, keepdims=True) ** 2 * np.cross(la, lb),
        axis=1)
    R = num / np.maximum(12.0 * V, 1e-300)
    return 3.0 * rin / np.maximum(R, 1e-300), V


def main():
    ap = argparse.ArgumentParser(description="montrer un maillage 3D")
    ap.add_argument("msh")
    ap.add_argument("--out", default=None)
    ap.add_argument("--split", action="store_true",
                    help="une figure PAR PANNEAU au lieu de la planche a trois")
    ap.add_argument("--zoom", type=float, default=0.030,
                    help="demi-largeur du zoom autour de l'impact [m]")
    ap.add_argument("--steps", type=float, nargs="*", default=None,
                    help="paliers vises : h1 r1 h2 r2 ... [m], pour le profil")
    ap.add_argument("--impact", type=float, nargs=3, default=None,
                    metavar=("X", "Y", "Z"),
                    help="point d'impact impose [m]. NECESSAIRE des que le "
                         "maillage contient PLUSIEURS CORPS : le defaut prend "
                         "le centre de la face du haut, qui est alors le "
                         "sommet du piston et non la surface de la roche.")
    a = ap.parse_args()

    P, T = read_msh22(a.msh)
    if a.impact:
        xC, yC, zTop = a.impact
    else:
        # point d'impact : centre de la face superieure
        zTop = P[:, 2].max()
        xC, yC = 0.5 * (P[:, 0].min() + P[:, 0].max()), \
                 0.5 * (P[:, 1].min() + P[:, 1].max())
    CE = P[T].mean(1)
    r = np.linalg.norm(CE - np.array([xC, yC, zTop]), axis=1)
    # taille realisee = arete moyenne du tetraedre
    L = np.mean([np.linalg.norm(P[T[:, i]] - P[T[:, j]], axis=1)
                 for i, j in EDGES], axis=0)
    qual, V = tet_quality(P, T)

    polys = slice_tets(P, T, 1, yC)
    bandeau = ("%s — %d tetraedres, %d noeuds | qualite : mediane %.2f, "
               "min %.3f, %d sous 0,10 | h de %.2f a %.1f mm"
               % (os.path.basename(a.msh), len(T), len(P),
                  np.median(qual), qual.min(), int((qual < 0.10).sum()),
                  L.min() * 1e3, L.max() * 1e3))
    out = a.out or (os.path.splitext(a.msh)[0] + "_vue.png")

    if a.split:
        # Une figure PAR PANNEAU (demande du 15/09) : la planche a trois cases
        # devient illisible reduite a une colonne d'article, et on ne veut
        # souvent qu'un seul des trois panneaux.
        base, ext = os.path.splitext(out)
        figs, A = {}, None
        for suff, taille in (("a_coupe", (9.0, 7.6)),
                             ("b_zoom", (7.2, 7.2)),
                             ("c_gradation", (8.4, 5.4))):
            f, ax = plt.subplots(figsize=taille)
            figs[suff] = (f, ax, base + "_" + suff + ext)
        A = figs["a_coupe"][1]
        Bx = figs["b_zoom"][1]
        C = figs["c_gradation"][1]
        fig = None
    else:
        fig = plt.figure(figsize=(14.0, 9.0))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0])
        A = fig.add_subplot(gs[0, :])
        Bx = fig.add_subplot(gs[1, 0])
        C = fig.add_subplot(gs[1, 1])
        figs = None

    # ---- (a) coupe exacte par y = yC -----------------------------------
    A.add_collection(PolyCollection(polys, facecolors="none",
                                    edgecolors="0.25", linewidths=0.18))
    A.set_xlim(P[:, 0].min(), P[:, 0].max())
    A.set_ylim(P[:, 2].min(), P[:, 2].max())
    A.set_aspect("equal")
    A.set_xlabel("x (m)")
    A.set_ylabel("z (m)")
    A.set_title("(a) Coupe exacte par le plan median — %d faces de section"
                % len(polys))
    A.plot([xC], [zTop], "v", color="crimson", ms=11, zorder=5)

    # ---- (b) zoom sur l'impact -----------------------------------------
    Bx.add_collection(PolyCollection(polys, facecolors="none",
                                     edgecolors="0.15", linewidths=0.45))
    Bx.set_xlim(xC - a.zoom, xC + a.zoom)
    Bx.set_ylim(zTop - 2 * a.zoom, zTop + 0.1 * a.zoom)
    Bx.set_aspect("equal")
    Bx.set_xlabel("x (m)")
    Bx.set_ylabel("z (m)")
    Bx.set_title("(b) Zone d'impact (%.0f mm de part et d'autre)"
                 % (a.zoom * 1e3))
    Bx.plot([xC], [zTop], "v", color="crimson", ms=13, zorder=5)

    # ---- (c) profil de taille realisee ---------------------------------
    rb = np.linspace(0, r.max(), 45)
    ib = np.clip(np.digitize(r, rb) - 1, 0, len(rb) - 2)
    med = np.array([np.median(L[ib == k]) if (ib == k).any() else np.nan
                    for k in range(len(rb) - 1)])
    lo = np.array([np.percentile(L[ib == k], 5) if (ib == k).any() else np.nan
                   for k in range(len(rb) - 1)])
    hi = np.array([np.percentile(L[ib == k], 95) if (ib == k).any() else np.nan
                   for k in range(len(rb) - 1)])
    rc = 0.5 * (rb[1:] + rb[:-1])
    C.fill_between(rc * 1e3, lo * 1e3, hi * 1e3, color="0.8",
                   label="5-95 % des elements")
    C.plot(rc * 1e3, med * 1e3, "-", color="k", lw=1.8, label="mediane")
    if a.steps:
        for k in range(0, len(a.steps) - 1, 2):
            h, rr = a.steps[k], a.steps[k + 1]
            C.plot([0, rr * 1e3], [h * 1e3, h * 1e3], "--", color="crimson",
                   lw=1.1)
            C.plot([rr * 1e3, rr * 1e3], [0, h * 1e3], ":", color="crimson",
                   lw=1.0)
        C.plot([], [], "--", color="crimson", label="paliers vises")
    C.set_xlabel("distance au point d'impact (mm)")
    C.set_ylabel("arete moyenne realisee (mm)")
    C.set_title("(c) Gradation REALISEE (controle falsifiant)")
    C.legend(fontsize=8)
    C.grid(alpha=0.3)

    if figs is not None:
        for suff in ("a_coupe", "b_zoom", "c_gradation"):
            f, _, chemin = figs[suff]
            f.suptitle(bandeau, fontsize=9.5)
            f.tight_layout()
            f.savefig(chemin, dpi=150, bbox_inches="tight")
            plt.close(f)
            print("wrote " + chemin)
    else:
        fig.suptitle(bandeau, fontsize=12)
        fig.tight_layout()
        fig.savefig(out, dpi=130)
        print("wrote " + out)
    print("  %d tetraedres, %d noeuds" % (len(T), len(P)))
    print("  arete : min %.3f mm, mediane %.3f mm, max %.2f mm"
          % (L.min() * 1e3, np.median(L) * 1e3, L.max() * 1e3))
    print("  qualite (1 = regulier) : mediane %.3f, min %.4f, "
          "%d elements sous 0,10 (%.3f %%)"
          % (np.median(qual), qual.min(), int((qual < 0.10).sum()),
             100.0 * (qual < 0.10).mean()))
    print("  volume total %.6f m3 (cylindre exact : compare a pi R2 H)"
          % V.sum())


if __name__ == "__main__":
    main()
