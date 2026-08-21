#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# parker_compare.py — H3, le PONT DE NON-REGRESSION du module hydro.
#
#   python bench_abuaisha/tools/parker_compare.py out_parker_c out_parker_hydro_c
#
# Le cas de Parker (annexe A d'AbuAisha et al. 2017) est le meme probleme
# physique par deux chemins numeriques differents :
#
#   * confiningPressure — une pression suiveuse imposee sur des faces
#     EXTERIEURES D'ORIGINE. Deja valide le 2026-08-18. Ne teste PAS le fluide.
#   * hydro = on        — la meme pression, mais portee par la CAVITE FLUIDE :
#     connexite depuis la source, volume par le theoreme de Green, chargement
#     des levres. C'est le module hydro qui est en cause.
#
# Si la cavite, la connexite et le SIGNE du chargement sont justes, les deux
# courbes doivent se superposer. Tout ecart est un bug, pas une physique.

# MESURE SIGNEE depuis le 2026-08-20 : une ouverture NEGATIVE est une
# interpenetration, donc un signe de chargement inverse. Avant cette date la
# mesure etait une valeur absolue et ce test ne pouvait pas voir ce bug-la.
#
# La solution analytique (Parker 1981, p. 33), en deformation plane :
#     w(x) = 2 sigma' (1 - nu^2) / E * sqrt(l^2 - x^2)
#
# CONVENTION, la meme que parker_check.py : l'eq. A.1 est la DEMI-ouverture
# classique de Sneddon pour une fissure de Griffith sous pression interne.
# L'ouverture TOTALE levre a levre vaut le DOUBLE, 4 sigma'(1-nu^2)/E sqrt(...).
# Comme on mesure l'ecart complet entre les deux levres, c'est a 2 x A.1 qu'il
# faut le comparer. (La figure A.21 de l'article nomme A.1 « aperture », ce qui
# preterait a confusion : elle trace la demi-ouverture.)
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tunnel_edz"))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

plt.rcParams.update({"font.size": 11, "figure.dpi": 110})


def aperture(run, half):
    """Ouverture SIGNEE levre a levre. Positive = la fissure S'OUVRE.

    (!) CORRECTIF DU 2026-08-20. Cette fonction rendait `g.max() - g.min()`,
    une valeur ABSOLUE : elle etait aveugle au signe et a valide une
    INTERPENETRATION comme une ouverture. C'est ce qui a fait passer H3 sur un
    module hydro dont le chargement etait de signe inverse (cf. le correctif
    de hydroForces(), src/FdemSolver.cpp).

    Les deux levres etant CONFONDUES a t = 0, leur position initiale ne peut
    pas les distinguer. Le discriminant est donc la GEOMETRIE DE L'ELEMENT
    PORTEUR : en FDEM chaque triangle possede ses propres noeuds (3 par
    element), et le centroide du triangle est d'un cote ou de l'autre de la
    ligne de discontinuite. On mesure alors

        w = <y des noeuds de la levre HAUTE> - <y des noeuds de la levre BASSE>

    qui change de signe si les levres s'interpenetrent.
    """
    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    if not el:
        raise SystemExit(f"aucune trame complete dans {run}")
    P0, C, _ = read_vtu(el[0], [])
    P, _, _ = read_vtu(el[-1], [])
    cy = 0.5 * (P0[:, 1].max() + P0[:, 1].min())
    cx = 0.5 * (P0[:, 0].max() + P0[:, 0].min())

    # cote de chaque noeud = signe de (y du centroide de SON triangle - cy).
    # Les noeuds etant dedoubles par element, chacun est ecrit exactement une
    # fois par cette affectation.
    side = np.zeros(len(P0))
    side[C.ravel()] = np.repeat(P0[C, 1].mean(axis=1), C.shape[1]) - cy

    on = (np.abs(P0[:, 1] - cy) < 1e-9) & (np.abs(P0[:, 0] - cx) <= half + 1e-9)
    xs, ys, sd = np.round(P0[on, 0], 9), P[on, 1], side[on]
    xa, wa = [], []
    for xv in np.unique(xs):
        g = np.abs(xs - xv) < 1e-12
        up, lo = g & (sd > 0), g & (sd < 0)
        if not up.any() or not lo.any():
            continue          # pointe de fissure : une seule levre presente
        xa.append(xv - cx)
        wa.append(ys[up].mean() - ys[lo].mean())
    o = np.argsort(xa)
    return np.array(xa)[o], np.array(wa)[o], len(el) - 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("conf", help="run par confiningPressure")
    ap.add_argument("hydro", help="run par le module hydro")
    ap.add_argument("--sig", type=float, default=2.0e6)
    ap.add_argument("--l", type=float, default=0.75)
    ap.add_argument("--E", type=float, default=45.0e9)
    ap.add_argument("--nu", type=float, default=0.2)
    ap.add_argument("--out", default="parker_H3_pont_hydro")
    a = ap.parse_args()

    root = os.path.join(HERE, "..", "..")
    runs = {}
    for tag, r in (("conf", a.conf), ("hydro", a.hydro)):
        p = r if os.path.isabs(r) else os.path.join(root, r)
        runs[tag] = aperture(p, a.l)

    xs = np.linspace(-a.l, a.l, 400)
    wth = 2.0 * (2.0 * a.sig * (1.0 - a.nu ** 2) / a.E) * np.sqrt(
        np.maximum(a.l ** 2 - xs ** 2, 0.0))          # x2 : ouverture TOTALE
    w0 = 2.0 * (2.0 * a.sig * (1.0 - a.nu ** 2) / a.E) * a.l

    fig, (ax, axr) = plt.subplots(
        1, 2, figsize=(13.2, 5.0), gridspec_kw={"width_ratios": [1.6, 1.0]})

    ax.plot(xs, wth * 1e3, "k-", lw=2.2, zorder=3,
            label="Sneddon / Parker — ouverture TOTALE (2 $\times$ éq. A.1)")
    style = {"conf": ("o", "#c8342b", "pression de confinement"),
             "hydro": ("s", "#1e4b8c", "module HYDRO (cavité fluide)")}
    txt = []
    for tag in ("conf", "hydro"):
        x, w, kf = runs[tag]
        m, c, lab = style[tag]
        wc = np.interp(0.0, x, w)
        e = 100.0 * (wc - w0) / w0
        ax.plot(x, w * 1e3, m, ms=4.0, mfc="none", mew=1.1, color=c,
                label=f"{lab} — centre {e:+.2f} %")
        txt.append((lab, wc, e))
    ax.set_xlabel("abscisse depuis le centre de la discontinuité  [m]")
    ax.set_ylabel("ouverture lèvre à lèvre  [mm]")
    ax.set_title("H3 — le même problème par deux chemins numériques")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower center", framealpha=0.95)

    # panneau de droite : l'ecart entre les DEUX CHEMINS, qui est la question
    xc, wc_, _ = runs["conf"]
    xh, wh_, _ = runs["hydro"]
    xi = np.linspace(-0.9 * a.l, 0.9 * a.l, 200)
    d = 100.0 * (np.interp(xi, xh, wh_) - np.interp(xi, xc, wc_)) \
        / np.maximum(np.interp(xi, xc, wc_), 1e-30)
    axr.axhline(0.0, color="k", lw=1.0)
    axr.plot(xi, d, "-", color="#1b8a3a", lw=1.8)
    axr.fill_between(xi, d, 0.0, color="#1b8a3a", alpha=0.15)
    axr.set_xlabel("abscisse  [m]")
    axr.set_ylabel("écart hydro − confinement  [%]")
    axr.set_title(f"écart entre les deux chemins\nmédian {np.median(np.abs(d)):.2f} %"
                  f", max {np.max(np.abs(d)):.2f} %")
    axr.grid(alpha=0.3)

    fig.suptitle("Ouverture d'une discontinuité sous pression uniforme de 2 MPa  —  "
                 "AbuAisha et al. (2017), annexe A", y=0.99, fontsize=12)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.text(0.5, 0.005,
             "Le confinement charge des faces extérieures d'origine : il ne teste pas le fluide. "
             "Le module hydro charge une cavité dont la connexité, le volume et le signe sont calculés.",
             ha="center", fontsize=9, style="italic", color="0.35")
    out = os.path.join(HERE, "..", a.out + ".png")
    fig.savefig(out, dpi=155)
    fig.savefig(out.replace(".png", ".pdf"))
    print(f"ecrit : {out}")
    print(f"  theorie w(0)          = {w0*1e3:.4f} mm")
    for lab, wc, e in txt:
        print(f"  {lab:34s} w(0) = {wc*1e3:.4f} mm  ({e:+.2f} %)")
    print(f"  ecart entre chemins   : median {np.median(np.abs(d)):.2f} %, "
          f"max {np.max(np.abs(d)):.2f} %")


if __name__ == "__main__":
    main()
