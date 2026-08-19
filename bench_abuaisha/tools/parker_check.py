#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# parker_check.py — LE controle analytique du benchmark AbuAisha et al. (2017).
#
#   python bench_abuaisha/tools/parker_check.py out_parker [--sig 2e6]
#          [--l 0.75] [--E 45e9] [--nu 0.2]
#
# Leur annexe A compare l'ouverture d'une discontinuite sous pression uniforme
# a la solution fermee de Parker (1981, p. 33), en deformation plane :
#
#     w(x) = 2 sigma' (1 - nu^2) / E * sqrt(l^2 - x^2)                  (A.1)
#
# ou sigma' est la contrainte effective d'ouverture (p moins la contrainte
# normale a la discontinuite), l la DEMI-longueur et w l'ouverture TOTALE.
#
# Verification du jeu de parametres de l'article (leur figure A.20) :
#   sigma' = 12 - 10 = 2 MPa, l = 0,75 m, nu = 0,2
#   w(0) = 2 x 2e6 x 0,96 x 0,75 / E = 2,88e6 / E
#   leur figure A.21 lit w(0) ~ 0,065 mm  ->  E = 44,3 GPa
# Le texte annonce « E=45 MPa », ce qui est une coquille : c'est 45 GPa, et
# l'arithmetique ci-dessus le confirme a 1,5 % pres. On retient 45 GPa.
#
# CE QUI EST MESURE : les levres etant des noeuds DEDOUBLES a la meme abscisse
# (greffon Crack de gmsh), l'ouverture est simplement l'ecart vertical entre
# les deux copies. Aucune interpolation, aucun post-traitement discutable.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_parker")
    ap.add_argument("--sig", type=float, default=2.0e6, help="sigma' [Pa]")
    ap.add_argument("--l", type=float, default=0.75, help="demi-longueur [m]")
    ap.add_argument("--E", type=float, default=45.0e9)
    ap.add_argument("--nu", type=float, default=0.2)
    ap.add_argument("--frame", type=int, default=-1)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", "..", a.run)
    out = os.path.join(HERE, "..", os.path.basename(run) + "_parker.png")

    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    if not el:
        raise SystemExit(f"aucune trame complete dans {run}")
    k = len(el) - 1 if a.frame < 0 else a.frame
    P0, C, _ = read_vtu(el[0], [])
    P, _, _ = read_vtu(el[k], [])
    cy = 0.5 * (P0[:, 1].max() + P0[:, 1].min())
    cx = 0.5 * (P0[:, 0].max() + P0[:, 0].min())

    # noeuds initialement SUR la ligne de fissure
    on = (np.abs(P0[:, 1] - cy) < 1e-9) & (np.abs(P0[:, 0] - cx) <= a.l + 1e-9)
    if on.sum() == 0:
        raise SystemExit("aucun noeud sur la ligne de fissure — le greffon "
                         "Crack a-t-il bien tourne au maillage ?")
    # EN FDEM CHAQUE ELEMENT A SES PROPRES NOEUDS : rockim annonce 3 x n_elem
    # noeuds. Sur la ligne de fissure il y a donc ~6 noeuds par abscisse (un par
    # element incident), pas 2. L'ouverture se mesure comme l'ECART TOTAL entre
    # le groupe du haut et celui du bas a chaque abscisse — robuste au nombre de
    # copies, contrairement a un appariement deux a deux.
    xs, ys = P0[on, 0], P[on, 1]
    xa, wa = [], []
    for xv in np.unique(np.round(xs, 9)):
        g = ys[np.abs(np.round(xs, 9) - xv) < 1e-12]
        if len(g) < 2:
            continue                        # pointe non dedoublee
        xa.append(xv - cx)
        wa.append(g.max() - g.min())        # ecart levre a levre COMPLET
    xa, wa = np.array(xa), np.array(wa)
    o = np.argsort(xa)
    xa, wa = xa[o], wa[o]
    if len(xa) < 5:
        raise SystemExit(f"seulement {len(xa)} paires de levres trouvees")

    # solution de Parker (1981), eq. A.1
    def parker(x):
        return (2.0 * a.sig * (1.0 - a.nu ** 2) / a.E) * np.sqrt(
            np.maximum(a.l ** 2 - x ** 2, 0.0))

    # CONVENTION. L'eq. A.1 est la DEMI-ouverture classique d'une fissure de
    # Griffith sous pression interne (Sneddon) : l'ouverture totale vaut
    # 4 sigma'(1-nu^2)/E sqrt(l^2-x^2), soit le DOUBLE. Leur figure A.21 la
    # nomme « aperture » mais trace bien A.1. On mesure l'ecart levre a levre
    # COMPLET, donc on le compare a 2 x A.1 — et on affiche les deux lectures.
    wth = 2.0 * parker(xa)
    interior = np.abs(xa) < 0.9 * a.l          # les pointes sont singulieres
    err = 100.0 * (wa[interior] - wth[interior]) / np.maximum(wth[interior], 1e-30)
    w0num = np.interp(0.0, xa, wa)
    w0th = 2.0 * parker(np.array([0.0]))[0]

    print(f"trame {k}/{len(el)-1}, {len(xa)} paires de levres")
    print(f"  ouverture au centre : numerique {w0num*1e3:.4f} mm | "
          f"Parker {w0th*1e3:.4f} mm | ecart {100*(w0num-w0th)/w0th:+.2f} %")
    print(f"  sur |x| < 0,9 l : ecart moyen {err.mean():+.2f} %, "
          f"ecart absolu median {np.median(np.abs(err)):.2f} %, "
          f"max {np.abs(err).max():.2f} %")

    xf = np.linspace(-a.l, a.l, 400)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xf, 2.0 * parker(xf) * 1e3, "-", color="0.15", lw=1.8,
            label="2 x Parker (1981), eq. A.1 = ouverture totale")
    ax.plot(xa, wa * 1e3, "o", mfc="none", mec="#C8342B", ms=5, mew=1.2,
            label="rockim FDEM")
    ax.set_xlabel("abscisse depuis le centre de la discontinuite [m]")
    ax.set_ylabel("ouverture totale levre a levre [mm]")
    ax.set_title(f"Ouverture sous pression uniforme — $\\sigma'$ = "
                 f"{a.sig/1e6:g} MPa, $l$ = {a.l:g} m, E = {a.E/1e9:g} GPa",
                 fontsize=10.5)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=155)
    print("ecrit :", os.path.normpath(out))


if __name__ == "__main__":
    main()
