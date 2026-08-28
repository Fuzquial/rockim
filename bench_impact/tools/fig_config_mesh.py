# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_config_mesh.py — la CONFIGURATION du modele et son MAILLAGE.
#
#   python bench_impact/tools/fig_config_mesh.py \
#          meshes/impact_fidele_s15.msh --stem bench_impact/figures/fig
#
# Deux figures, une idee chacune :
#   fig_config.pdf  — le montage : 6 corps, 3 materiaux, conditions aux
#                     limites et les deux JEUX (piston/bit 0,2 mm et
#                     insert/roche 0,02 mm) qui fixent la chronologie.
#   fig_mesh.pdf    — le maillage : gradation du bord vers l'insert, avec un
#                     zoom sur la zone d'impact ou tout se joue.
#
# La coupe est une INTERSECTION EXACTE des tetraedres avec le plan y = 0
# (marching tetrahedra), et non une tranche epaisse : les aretes tracees sont
# donc les vraies traces des faces sur le plan, sans recouvrement.
# ---------------------------------------------------------------------------
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42})

NAMES = {1: "rock", 2: "insert", 3: "bit", 4: "piston", 5: "circlip",
         6: "plate"}
# couleur PAR MATERIAU (Table 4 de l'article) : la lecture doit etre
# « qui est en acier, qui est en carbure », pas « qui est quel corps ».
MAT = {1: "rock", 2: "carbide", 3: "steel", 4: "steel", 5: "carbide",
       6: "steel"}
COL = {"rock": "#c9b79c", "steel": "#8fa8bf", "carbide": "#d98c45"}
EDGE = {"rock": "#6b5f4c", "steel": "#4a6076", "carbide": "#8a5320"}


def read_msh(path):
    """Noeuds et tetraedres par groupe physique d'un .msh gmsh 2.2."""
    L = open(path).read().splitlines()
    i = L.index("$Nodes")
    nn = int(L[i + 1])
    nid = {}
    P = np.zeros((nn, 3))
    for k in range(nn):
        p = L[i + 2 + k].split()
        nid[int(p[0])] = k
        P[k] = [float(x) for x in p[1:4]]
    j = L.index("$Elements")
    ne = int(L[j + 1])
    tets, phys = [], []
    for k in range(ne):
        p = L[j + 2 + k].split()
        if int(p[1]) != 4:                      # 4 = tetraedre
            continue
        tets.append([nid[int(x)] for x in p[-4:]])
        phys.append(int(p[3]))
    return P, np.array(tets), np.array(phys)


def cut_y0(P, tets, phys, y0=0.0):
    """Intersection exacte des tets avec le plan y = y0.

    Renvoie (polygones, groupe). Chaque tet coupe donne un triangle (1 noeud
    d'un cote) ou un quadrilatere (2-2). Les sommets sont ordonnes par angle
    autour du centroide — licite, le polygone etant convexe et plan.
    """
    polys, grp = [], []
    E = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    d = P[:, 1] - y0
    for T, ph in zip(tets, phys):
        dd = d[T]
        if np.all(dd > 0) or np.all(dd < 0):
            continue
        pts = []
        for a, b in E:
            da, db = dd[a], dd[b]
            if da * db < 0:                     # l'arete traverse le plan
                w = da / (da - db)
                pts.append(P[T[a]] + w * (P[T[b]] - P[T[a]]))
        if len(pts) < 3:
            continue
        Q = np.array(pts)[:, [0, 2]] * 1e3      # projection (x, z) en mm
        c = Q.mean(axis=0)
        Q = Q[np.argsort(np.arctan2(Q[:, 1] - c[1], Q[:, 0] - c[0]))]
        polys.append(Q)
        grp.append(ph)
    return polys, np.array(grp)


def draw(ax, polys, grp, lw=0.0, alpha=1.0):
    for ph in sorted(set(grp)):
        sel = [p for p, g in zip(polys, grp) if g == ph]
        m = MAT[ph]
        ax.add_collection(PolyCollection(
            sel, facecolors=COL[m], edgecolors=EDGE[m] if lw else "none",
            linewidths=lw, alpha=alpha, zorder=2))


# ---------------------------------------------------------------------------
def fig_config(polys, grp, stem):
    fig, (ax, axz) = plt.subplots(
        1, 2, figsize=(7.0, 6.2), gridspec_kw={"width_ratios": [1.55, 1]})

    draw(ax, polys, grp, lw=0.0)
    ax.set_xlim(-135, 135)
    ax.set_ylim(-200, 545)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ (mm)")
    ax.set_ylabel(r"$z$ (mm)")

    # --- le chargement et les appuis ---------------------------------------
    ax.annotate("", xy=(0, 452), xytext=(0, 512),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#b3242b"))
    ax.annotate("piston lancé à\n" r"$v_0 = 10{,}66$ m/s", xy=(-2, 470),
                xytext=(-18, 420), color="#b3242b", fontsize=8.5,
                va="center", ha="right",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#b3242b"))
    for xx in np.linspace(-125, 125, 21):       # encastrement de la base
        ax.plot([xx, xx - 9], [-150, -161], lw=0.6, color="0.25",
                solid_capstyle="butt", zorder=3)
    ax.plot([-125, 125], [-150, -150], lw=1.1, color="0.25", zorder=3)
    ax.text(0, -185, "base encastrée", ha="center", va="top", fontsize=8.5,
            color="0.25")

    # --- qui est quoi : SEULS les grands corps ici, les petites pieces sont
    #     etiquetees dans le zoom (b), ou il y a la place de les distinguer --
    for nm, xt, zt, xa, za in [("piston", 46, 400, 12, 400),
                               ("taillant", 46, 170, 14, 170),
                               ("roche", 60, -80, 20, -80)]:
        ax.annotate(nm, xy=(xa, za), xytext=(xt, zt), fontsize=8.5,
                    ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
    ax.annotate("outil (détail en b)", xy=(22, 33), xytext=(52, 90),
                fontsize=8, ha="left", va="center", color="#b3242b",
                arrowprops=dict(arrowstyle="->", lw=0.7, color="#b3242b"))
    # legende des MATERIAUX (la lecture utile : qui est acier, qui carbure)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=COL[m], ec=EDGE[m],
                                     lw=0.5, label=l)
                       for m, l in [("rock", "calcaire St Anne"),
                                    ("steel", "acier"),
                                    ("carbide", "carbure")]],
              loc="upper right", fontsize=7.8, frameon=True, framealpha=1,
              edgecolor="0.7", handlelength=1.2, borderpad=0.5)
    ax.set_title("(a) le montage — 6 corps, 3 matériaux", fontsize=9.5,
                 loc="left")

    # --- zoom sur le jeu, qui commande la chronologie ----------------------
    draw(axz, polys, grp, lw=0.25)
    axz.set_xlim(-30, 30)
    axz.set_ylim(-9, 46)
    axz.set_aspect("equal")
    axz.set_xlabel(r"$x$ (mm)")
    axz.set_ylabel(r"$z$ (mm)")
    axz.axhline(0, lw=0.7, ls=(0, (4, 3)), color="#b3242b", zorder=5)
    axz.annotate(r"jeu $= 20\ \mu$m" "\n(invisible ici)",
                 xy=(-7.5, 0.2), xytext=(-29, 9), fontsize=7.8,
                 color="#b3242b", va="center",
                 arrowprops=dict(arrowstyle="->", lw=0.7, color="#b3242b"))
    for nm, xt, zt, xa, za in [("insert", 17, 12, 7, 12),
                               ("circlip", 20, 27, 17.5, 31),
                               ("plaque", 20, 42, 24, 36),
                               ("taillant", -29, 27, -13, 26)]:
        axz.annotate(nm, xy=(xa, za), xytext=(xt, zt), fontsize=8,
                     ha="left", va="center",
                     arrowprops=dict(arrowstyle="-", lw=0.5, color="0.35"))
    axz.set_title("(b) zoom : l'insert est POSÉ sur la roche", fontsize=9.5,
                  loc="left")

    for a in (ax, axz):
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig("%s_config.%s" % (stem, ext), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)


def fig_mesh(polys, grp, P, tets, phys, stem):
    fig, (ax, axz) = plt.subplots(
        1, 2, figsize=(7.0, 6.2), gridspec_kw={"width_ratios": [1.55, 1]})

    draw(ax, polys, grp, lw=0.16, alpha=0.95)
    ax.set_xlim(-135, 135)
    ax.set_ylim(-200, 545)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ (mm)")
    ax.set_ylabel(r"$z$ (mm)")
    ax.add_patch(plt.Rectangle((-30, -32), 60, 62, fill=False, lw=0.9,
                               ec="#b3242b", zorder=6))
    ax.set_title("(a) la coupe $y = 0$ — %d tétraèdres au total"
                 % len(tets), fontsize=9.5, loc="left")

    draw(axz, polys, grp, lw=0.35, alpha=0.95)
    axz.set_xlim(-30, 30)
    axz.set_ylim(-32, 30)
    axz.set_aspect("equal")
    axz.set_xlabel(r"$x$ (mm)")
    axz.set_ylabel(r"$z$ (mm)")
    axz.set_title("(b) sous l'insert : la zone qui décide",
                  fontsize=9.5, loc="left")
    for s in axz.spines.values():
        s.set_edgecolor("#b3242b")
        s.set_linewidth(0.9)

    # --- gradation reelle DE LA ROCHE : arete moyenne par distance a l'axe -
    #     (la roche seule : c'est elle qui est graduee, l'acier ne porte que
    #      l'onde et garde une taille uniforme)
    v = P[tets[phys == 1]]
    h = np.mean([np.linalg.norm(v[:, a] - v[:, b], axis=1)
                 for a, b in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3),
                              (2, 3)]], axis=0) * 1e3
    ctr = v.mean(axis=1)
    r = np.hypot(ctr[:, 0], ctr[:, 1]) * 1e3
    txt = ("roche : arête moyenne\n"
           r"$r < 25$ mm : %.1f" "\n"
           r"$25$–$50$ : %.1f" "\n"
           r"$r > 50$ : %.1f mm" % (
               h[r < 25].mean(), h[(r >= 25) & (r < 50)].mean(),
               h[r >= 50].mean()))
    ax.text(-131, 535, txt, fontsize=7.6, va="top", color="0.25",
            bbox=dict(fc="white", ec="0.7", lw=0.5, pad=3.0))

    for a in (ax, axz):
        a.spines[["top", "right"]].set_visible(False) if a is ax else None
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig("%s_mesh.%s" % (stem, ext), dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("msh")
    ap.add_argument("--stem", default="bench_impact/figures/fig_imperial")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)
    P, tets, phys = read_msh(a.msh)
    polys, grp = cut_y0(P, tets, phys)
    print("coupe y = 0 : %d polygones sur %d tets" % (len(polys), len(tets)))
    fig_config(polys, grp, a.stem)
    fig_mesh(polys, grp, P, tets, phys, a.stem)
    print("ecrit : %s_config.pdf / %s_mesh.pdf" % (a.stem, a.stem))
