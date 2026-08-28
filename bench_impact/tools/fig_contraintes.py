# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_contraintes.py — LE CHAMP DE CONTRAINTE sous l'insert, en coupe exacte.
#
#   python bench_impact/tools/fig_contraintes.py out_imperial/fdem3d_0007.vtu \
#          --stem bench_impact/figures/sig_A --t-us 96.2
#
# Pourquoi cette figure : le correctif `jointShearRange = coulomb` porte sur
# la rupture en CISAILLEMENT des joints COMPRIMES. Ces deux mots se lisent
# ici — (a) montre que sous l'insert le champ est en compression (Hertz),
# (b) montre que le cisaillement y est pourtant intense. C'est exactement la
# combinaison que la plage cohesion-seule rendait incassable.
#
# La coupe est une INTERSECTION EXACTE des tetraedres avec le plan y = y_axe
# (marching tetrahedra), pas une tranche epaisse : chaque polygone porte la
# contrainte de SON tetraedre, sans moyenne ni recouvrement.
#
# CONVENTION DE SIGNE : sigma1 est la contrainte principale MAJEURE telle que
# rockim l'ecrit (traction positive). Sous l'insert elle est negative.
# ---------------------------------------------------------------------------
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import TwoSlopeNorm

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42})

CX = CY = 0.125          # axe de l'impact (le solveur translate a lo = 0)
Z_SURF = 0.15            # sommet de la roche
PH_ROCK = 0              # phase 0 = rock (ordre de la cle `phases`)


def read_vtu(path):
    s = open(path, encoding="utf-8", errors="ignore").read()

    def arr(name):
        m = re.search(r'Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % name, s,
                      re.S)
        return None if m is None else np.fromstring(m.group(1), sep=" ")

    m = re.search(r"<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>", s,
                  re.S)
    P = np.fromstring(m.group(1), sep=" ").reshape(-1, 3)
    con = arr("connectivity").astype(int).reshape(-1, 4)
    f = {n: arr(n) for n in ("sigma1", "tauMax", "vonMises", "phase")}
    return P, con, f


def cut(P, con, vals, y0):
    """Polygones de l'intersection des tets avec y = y0, et valeur par tet."""
    E = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    d = P[:, 1] - y0
    dd = d[con]
    hit = ~(np.all(dd > 0, axis=1) | np.all(dd < 0, axis=1))
    polys, keep = [], []
    for i in np.nonzero(hit)[0]:
        T, e = con[i], dd[i]
        pts = []
        for a, b in E:
            if e[a] * e[b] < 0:
                w = e[a] / (e[a] - e[b])
                pts.append(P[T[a]] + w * (P[T[b]] - P[T[a]]))
        if len(pts) < 3:
            continue
        Q = np.array(pts)
        Q = np.c_[(Q[:, 0] - CX) * 1e3, (Q[:, 2] - Z_SURF) * 1e3]
        c = Q.mean(axis=0)
        polys.append(Q[np.argsort(np.arctan2(Q[:, 1] - c[1],
                                             Q[:, 0] - c[0]))])
        keep.append(i)
    return polys, {k: v[keep] for k, v in vals.items()}


def panel(ax, polys, v, cmap, norm, title, cbl, half, deep, other=None):
    # l'outil (carbure) en gris uni : il porte le GPa, le colorer sur
    # l'echelle de la roche ne dirait rien et detournerait l'oeil.
    if other:
        ax.add_collection(PolyCollection(other, facecolors="0.82",
                                         edgecolors="0.62", linewidths=0.15,
                                         zorder=1))
    pc = PolyCollection(polys, array=v, cmap=cmap, norm=norm,
                        edgecolors="none", zorder=2)
    ax.add_collection(pc)
    ax.set_xlim(-half, half)
    ax.set_ylim(-deep, 6)
    ax.set_aspect("equal")
    ax.axhline(0, lw=0.6, color="0.35", zorder=4)
    ax.set_xlabel(r"$x$ depuis l'axe de l'insert  [mm]")
    ax.set_ylabel(r"$z$  [mm]")
    ax.set_title(title, fontsize=9.2, loc="left")
    cb = ax.figure.colorbar(pc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(cbl, fontsize=8.2)
    cb.ax.tick_params(labelsize=7.5)
    return pc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("vtu")
    ap.add_argument("--stem", default="bench_impact/figures/contraintes")
    ap.add_argument("--t-us", type=float, default=None)
    ap.add_argument("--half", type=float, default=18.0)
    ap.add_argument("--deep", type=float, default=18.0)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)

    P, con, f = read_vtu(a.vtu)
    polys, v = cut(P, con, f, CY)
    # l'ECHELLE se regle sur la ROCHE : l'insert en carbure porte le GPa et
    # ecraserait toute la dynamique du champ qu'on veut lire.
    rock = v["phase"] == PH_ROCK
    s1 = v["sigma1"] / 1e6
    tau = v["tauMax"] / 1e6
    zone = np.array([np.abs(p[:, 0]).min() < a.half and
                     p[:, 1].min() > -a.deep for p in polys])
    m = rock & zone
    print("coupe : %d polygones, dont %d de roche dans la fenetre"
          % (len(polys), m.sum()))
    print("  sigma1 roche : %.1f a %.1f MPa" % (s1[m].min(), s1[m].max()))
    print("  tauMax roche : %.1f a %.1f MPa" % (tau[m].min(), tau[m].max()))

    rp = [p for p, k in zip(polys, rock) if k]
    op = [p for p, k in zip(polys, rock) if not k]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.5))
    lim = np.percentile(np.abs(s1[m]), 99.5)
    panel(axes[0], rp, s1[rock], "RdBu_r",
          TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim),
          r"(a) $\sigma_1$ : bleu = compression",
          r"$\sigma_1$  [MPa]", a.half, a.deep, op)
    tl = np.percentile(tau[m], 99.5)
    panel(axes[1], rp, tau[rock], "magma_r", plt.Normalize(0, tl),
          r"(b) cisaillement $\tau_{\max}$",
          r"$\tau_{\max}$  [MPa]", a.half, a.deep, op)
    # le repere qui donne son sens au panneau (b)
    axes[1].text(-a.half + 0.7, -a.deep + 1.0,
                 r"cohésion $c = 18{,}8$ MPa" "\n"
                 r"$\tau_{\max}$ culmine à %.0f MPa (%.0f$\times$ $c$)"
                 % (tau[m].max(), tau[m].max() / 18.8),
                 fontsize=7.6, va="bottom", color="0.2",
                 bbox=dict(fc="white", ec="0.7", lw=0.5, pad=2.5))

    ttl = "Le champ sous l'insert"
    if a.t_us:
        ttl += r"   ($t = %.1f\ \mu$s)" % a.t_us
    fig.suptitle(ttl, fontsize=11, y=0.99)
    fig.text(0.5, 0.005, "coupe exacte $y = y_{\\mathrm{axe}}$ (marching "
             "tetrahedra) — chaque polygone porte la contrainte de SON "
             "tétraèdre ; échelle réglée sur la roche seule",
             ha="center", fontsize=7.6, color="0.35", style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=200, bbox_inches="tight")
    print("ecrit : %s.pdf et .png" % a.stem)
