# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_cratere.py — LE CRATERE VU D'EN HAUT, EN ELEMENTS (leur convention).
#
#   python bench_impact/tools/fig_cratere.py out_imperial_coulomb 9 \
#          --stem bench_impact/figures/crat_B9 --t-us 123.7
#
# POURQUOI CE SCRIPT, ET EN QUOI IL DIFFERE DE fig_vue_dessus.py.
# Yang et al. ne dessinent PAS les joints. Leur fig. 8 (bloc « rock cracking »),
# leur fig. 11 (crateres simules) et leurs fig. 14-15 montrent des ELEMENTS :
# la roche telle qu'on la voit, coloree par l'endommagement ou par le fragment
# auquel elle appartient. Un reseau de facettes de joints, meme trie par
# profondeur, reste un objet de modelisateur : il montre les INTERFACES, pas
# la matiere. D'ou cette vue-ci, qui rend les TETRAEDRES.
#
# Rendu : algorithme du peintre sur les tetraedres de ROCHE, les plus profonds
# d'abord, remplis de facon OPAQUE. Chaque tetraedre est projete sur (x, y) et
# dessine par son ENVELOPPE CONVEXE 2D (triangle ou quadrilatere selon
# l'orientation) — pas par une face choisie arbitrairement.
#
#   (a) RELIEF : couleur = profondeur du sommet vu. C'est la forme du cratere.
#   (b) FRAGMENTS : la roche encore solidaire en gris, les fragments detaches
#       en couleur (un ton par fragment) — l'equivalent de leur fig. 11.
#
# Le champ `fragment` du VTU vaut 0 pour le bloc intact ; les identifiants 1 a
# 3 sont l'OUTIL (insert+bit+circlip), le piston et la plaque — ecartes ici
# par le filtre phase = roche. Les fragments de roche commencent donc a 4.
# ---------------------------------------------------------------------------
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import CX, CY, Z_SURF, frame_times

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42})

PH_ROCK = 0


def read_elems(path):
    s = open(path, encoding="utf-8", errors="ignore").read()

    def arr(n):
        m = re.search(r'Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % n, s, re.S)
        return None if m is None else np.fromstring(m.group(1), sep=" ")

    m = re.search(r"<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>", s,
                  re.S)
    P = np.fromstring(m.group(1), sep=" ").reshape(-1, 3)
    con = arr("connectivity").astype(int).reshape(-1, 4)
    return P, con, arr("phase"), arr("fragment")


def hull4(Q):
    """Enveloppe convexe de 4 points 2D (monotone chain). Rend 3 ou 4 points.

    Une projection de tetraedre donne soit un quadrilatere, soit un triangle
    avec un sommet interieur : un tri par angle inclurait ce sommet interieur
    et creuserait un faux creneau. D'ou le vrai calcul d'enveloppe.
    """
    p = Q[np.lexsort((Q[:, 1], Q[:, 0]))]

    def half(pts):
        h = []
        for q in pts:
            while len(h) >= 2 and (h[-1][0] - h[-2][0]) * (q[1] - h[-2][1]) - \
                    (h[-1][1] - h[-2][1]) * (q[0] - h[-2][0]) <= 0:
                h.pop()
            h.append(q)
        return h[:-1]

    return np.array(half(p) + half(p[::-1]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("frame", type=int)
    ap.add_argument("--stem", default="bench_impact/figures/cratere")
    ap.add_argument("--t-us", type=float, default=None)
    ap.add_argument("--half", type=float, default=14.0)
    ap.add_argument("--vis-mm", type=float, default=3.0,
                    help="profondeur au-dela de laquelle un tet est considere"
                         " comme jamais visible (reglage de l echelle)")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)

    P, con, ph, frag = read_elems(
        os.path.join(a.run, "fdem3d_%04d.vtu" % a.frame))
    V = P[con]                                        # (n, 4, 3)
    ctr = V.mean(axis=1)
    rock = ph == PH_ROCK
    win = ((np.abs(ctr[:, 0] - CX) * 1e3 < a.half + 3) &
           (np.abs(ctr[:, 1] - CY) * 1e3 < a.half + 3) & rock)
    idx = np.nonzero(win)[0]
    z_haut = V[idx, :, 2].max(axis=1)                 # sommet vu
    prof = (Z_SURF - z_haut) * 1e3                    # mm sous la surface
    ordre = np.argsort(-prof)                         # peintre : profond d'abord
    idx, prof = idx[ordre], prof[ordre]

    polys = [hull4(np.c_[(V[i, :, 0] - CX) * 1e3, (V[i, :, 1] - CY) * 1e3])
             for i in idx]
    fr = frag[idx]
    detache = fr >= 4                                 # 0 = bloc, 1-3 = outil

    print("tetraedres de roche dessines : %d (fenetre %.0f mm)"
          % (len(polys), a.half))
    print("  profondeur vue : %.2f a %.2f mm" % (prof.min(), prof.max()))
    print("  dont detaches  : %d tets, %d fragments distincts"
          % (detache.sum(), len(np.unique(fr[detache]))))

    fig, ax = plt.subplots(1, 2, figsize=(10.0, 5.0))

    # --- (a) relief ---------------------------------------------------------
    # L'ECHELLE se regle sur la surface REELLEMENT VUE. Normaliser sur tous les
    # tetraedres de la fenetre inclut ceux enfouis a plusieurs centimetres, qui
    # sont recouverts et ne s'affichent jamais : la dynamique du cratere s'y
    # ecrase a rien (constate : P92 = 28 mm pour un cratere de 0,5 mm).
    vis = prof < a.vis_mm
    vmax = float(np.percentile(prof[vis], 98)) if vis.any() else 1.0
    print("  echelle du relief : 0 a %.2f mm (sur les %d tets a moins de "
          "%.0f mm)" % (vmax, vis.sum(), a.vis_mm))
    pc = PolyCollection(polys, array=prof, cmap="YlOrBr",
                        norm=plt.Normalize(max(prof.min(), -0.5), vmax),
                        edgecolors="none")
    ax[0].add_collection(pc)
    cb = fig.colorbar(pc, ax=ax[0], fraction=0.046, pad=0.03)
    cb.set_label("profondeur de la surface vue  [mm]", fontsize=8.2)
    cb.ax.tick_params(labelsize=7.5)
    ax[0].set_title("(a) le relief : la forme du cratère", fontsize=9.4,
                    loc="left")

    # --- (b) fragments ------------------------------------------------------
    rng = np.random.default_rng(1)
    tons = plt.get_cmap("turbo")(rng.random(int(frag.max()) + 2))
    cols = [("0.78" if not d else tons[int(k)])
            for d, k in zip(detache, fr)]
    ax[1].add_collection(PolyCollection(polys, facecolors=cols,
                                        edgecolors="none"))
    ax[1].set_title("(b) les fragments détachés (gris = roche solidaire)",
                    fontsize=9.4, loc="left")
    ax[1].text(0.03, 0.045, "%d tétraèdres détachés\nen %d fragments"
               % (detache.sum(), len(np.unique(fr[detache]))),
               transform=ax[1].transAxes, fontsize=8, va="bottom",
               bbox=dict(fc="white", ec="0.75", lw=0.5, pad=3))

    for A in ax:
        A.set_xlim(-a.half, a.half)
        A.set_ylim(-a.half, a.half)
        A.set_aspect("equal")
        A.set_xlabel(r"$x$  [mm]")
        A.set_ylabel(r"$y$  [mm]")

    ttl = "Le cratère en ÉLÉMENTS, vu d'en haut"
    if a.t_us:
        ttl += r"   ($t = %.1f\ \mu$s)" % a.t_us
    fig.suptitle(ttl, fontsize=11, y=0.99)
    fig.text(0.5, 0.005, "rendu ÉLÉMENTS (leur convention : ils ne dessinent "
             "pas les joints) — tétraèdres de roche triés par profondeur et "
             "remplis opaques, projetés par leur enveloppe convexe",
             ha="center", fontsize=7.6, color="0.35", style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=200, bbox_inches="tight")
    print("ecrit : %s.pdf et .png" % a.stem)
