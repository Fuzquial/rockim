# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_deformation.py — LE DEPLACEMENT de la roche sous l'insert.
#
#   python bench_impact/tools/fig_deformation.py out_imperial_coulomb 8 \
#          --stem bench_impact/figures/def_B8 --t-us 110
#
# Les VTU de rockim n'ecrivent AUCUN champ de deformation. Mais ils ecrivent
# les positions, et les points se correspondent d'une frame a l'autre
# (connectivite identique, verifie) : le deplacement u = P(t) - P(0) est donc
# une mesure DIRECTE, pas une reconstruction.
#
#   (a) DEPLACEMENT VERTICAL u_z : l'enfoncement, et jusqu'ou il se propage.
#   (b) MAILLAGE DEFORME, amplifie : le cratere se voit en forme, la ou un
#       champ colore ne donne qu'une intensite. L'amplification est ECRITE
#       sur la figure — sans elle, 0,4 mm sur 36 mm de fenetre est invisible.
#
# Coupe exacte y = y_axe (marching tetrahedra), comme fig_contraintes.py.
# ---------------------------------------------------------------------------
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42})

CX = CY = 0.125
Z_SURF = 0.15
PH_ROCK = 0


def read(path, want_fields=True):
    s = open(path, encoding="utf-8", errors="ignore").read()
    m = re.search(r"<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>", s,
                  re.S)
    P = np.fromstring(m.group(1), sep=" ").reshape(-1, 3)
    c = re.search(r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S)
    con = np.fromstring(c.group(1), sep=" ").astype(int).reshape(-1, 4)
    ph = None
    if want_fields:
        m = re.search(r'Name="phase"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S)
        ph = np.fromstring(m.group(1), sep=" ")
    return P, con, ph


def cut(P, con, y0, amp=1.0, P0=None):
    """Polygones de l'intersection avec y = y0, dans la config DEFORMEE.

    Si amp != 1, la position tracee est P0 + amp (P - P0) : le plan de coupe
    reste celui de la config NON deformee, sans quoi l'amplification ferait
    entrer et sortir des tetraedres de la coupe.
    """
    E = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    ref = P if P0 is None else P0
    d = ref[:, 1] - y0
    dd = d[con]
    hit = ~(np.all(dd > 0, axis=1) | np.all(dd < 0, axis=1))
    Q = P if amp == 1.0 else P0 + amp * (P - P0)
    polys, keep = [], []
    for i in np.nonzero(hit)[0]:
        T, e = con[i], dd[i]
        pts = []
        for a, b in E:
            if e[a] * e[b] < 0:
                w = e[a] / (e[a] - e[b])
                pts.append(Q[T[a]] + w * (Q[T[b]] - Q[T[a]]))
        if len(pts) < 3:
            continue
        R = np.array(pts)
        R = np.c_[(R[:, 0] - CX) * 1e3, (R[:, 2] - Z_SURF) * 1e3]
        cc = R.mean(axis=0)
        polys.append(R[np.argsort(np.arctan2(R[:, 1] - cc[1],
                                             R[:, 0] - cc[0]))])
        keep.append(i)
    return polys, np.array(keep)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("frame", type=int)
    ap.add_argument("--stem", default="bench_impact/figures/deformation")
    ap.add_argument("--t-us", type=float, default=None)
    ap.add_argument("--half", type=float, default=18.0)
    ap.add_argument("--deep", type=float, default=18.0)
    ap.add_argument("--amp", type=float, default=8.0)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)

    P0, con, ph = read(os.path.join(a.run, "fdem3d_0000.vtu"))
    Pn, con2, _ = read(os.path.join(a.run, "fdem3d_%04d.vtu" % a.frame),
                       False)
    assert np.array_equal(con, con2), "connectivite differente entre frames"
    u = Pn - P0
    uz = -u[:, 2]                      # enfoncement compte POSITIF vers le bas
    uz_t = uz[con].mean(axis=1) * 1e6  # par tetraedre, en microns

    polys, keep = cut(Pn, con, CY, 1.0, P0)
    rock = ph[keep] == PH_ROCK
    v = uz_t[keep]
    zone = np.array([np.abs(p[:, 0]).min() < a.half and
                     p[:, 1].min() > -a.deep for p in polys])
    m = rock & zone
    print("coupe %d polygones, %d de roche dans la fenetre" % (len(polys),
                                                               m.sum()))
    print("  u_z roche : %.1f a %.1f um" % (v[m].min(), v[m].max()))

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.5))

    rp = [p for p, k in zip(polys, rock) if k]
    op = [p for p, k in zip(polys, rock) if not k]
    ax = axes[0]
    ax.add_collection(PolyCollection(op, facecolors="0.82",
                                     edgecolors="0.62", linewidths=0.15))
    pc = PolyCollection(rp, array=v[rock], cmap="YlGnBu",
                        norm=plt.Normalize(0, np.percentile(v[m], 99.5)),
                        edgecolors="none")
    ax.add_collection(pc)
    cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"$u_z$ (enfoncement)  [$\mu$m]", fontsize=8.2)
    cb.ax.tick_params(labelsize=7.5)
    ax.set_title(r"(a) déplacement vertical $u_z$", fontsize=9.2, loc="left")

    # --- (b) le maillage DEFORME, amplifie ---------------------------------
    # SEULE LA ROCHE est amplifiee : amplifier l'outil le ferait plonger huit
    # fois trop profond, ce qui donnerait a lire un enfoncement imaginaire.
    # L'outil est trace a sa position VRAIE, en gris, comme repere.
    pa, ka = cut(Pn, con, CY, a.amp, P0)
    ra = ph[ka] == PH_ROCK
    ax2 = axes[1]
    ax2.add_collection(PolyCollection(op, facecolors="none",
                                      edgecolors="0.55", linewidths=0.3,
                                      zorder=3))
    ax2.add_collection(PolyCollection(
        [p for p, k in zip(pa, ra) if k], facecolors="#e8dcc6",
        edgecolors="#6b5f4c", linewidths=0.22, zorder=2))
    ax2.set_title(r"(b) roche déformée, amplifiée $\times %g$" % a.amp,
                  fontsize=9.2, loc="left")
    ax2.text(-a.half + 0.6, -a.deep + 0.8,
             "l'outil (gris) est à sa position VRAIE,\nnon amplifiée",
             fontsize=7.4, color="0.35", va="bottom")

    for x in (ax, ax2):
        x.set_xlim(-a.half, a.half)
        x.set_ylim(-a.deep, 8)
        x.set_aspect("equal")
        x.axhline(0, lw=0.6, ls=(0, (4, 3)), color="0.35", zorder=5)
        x.set_xlabel(r"$x$ depuis l'axe de l'insert  [mm]")
        x.set_ylabel(r"$z$  [mm]")

    ttl = "La déformation de la roche"
    if a.t_us:
        ttl += r"   ($t = %.1f\ \mu$s)" % a.t_us
    fig.suptitle(ttl, fontsize=11, y=0.99)
    fig.text(0.5, 0.005, "$u = P(t) - P(0)$, mesure DIRECTE sur les positions "
             "(la connectivité est identique d'une frame à l'autre) ; "
             "le trait tireté est la surface initiale",
             ha="center", fontsize=7.6, color="0.35", style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=200, bbox_inches="tight")
    print("ecrit : %s.pdf et .png" % a.stem)
