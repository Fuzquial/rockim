#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# crack_orientation.py — le X conjugué est-il là, ou est-ce que je l'ai vu
# parce que je le cherchais ?
#
#   python tunnel_edz/tools/crack_orientation.py out_tun_ref_iso [out_tun_ref]
#
# THEORIE. Autour d'une cavité sous contrainte hydrostatique, la contrainte
# principale la plus compressive est TANGENTIELLE. Dans un matériau de
# Mohr-Coulomb, les surfaces de rupture par cisaillement font un angle
# +-(45 - phi/2) avec cette direction, soit +-28 deg pour phi = 34 deg. Un
# faciès conjugué doit donc montrer DEUX pics symétriques a +-28 deg de la
# tangente locale — et rien de tel si les fissures ne font que suivre le
# maillage.
#
# MESURE. Pour chaque fissure de cisaillement (breakMode = 2) on calcule
# l'angle entre sa direction propre et la tangente locale (perpendiculaire au
# rayon depuis le centre du tunnel), replié dans [-90, 90].
#
# CONTROLE. Le meme calcul est fait sur les ARETES NON ROMPUES de la meme
# zone : c'est l'offre du maillage. Si les fissures rompues montrent des pics
# que les arêtes disponibles n'ont pas, la sélection est PHYSIQUE. Si les deux
# histogrammes se ressemblent, on ne mesure que le mailleur.
# ---------------------------------------------------------------------------
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

PHI = 34.0
THEO = 45.0 - 0.5 * PHI          # 28 deg


def profile(run, cx=50.0, cy=50.0, rmin=6.0, rmax=18.0):
    jn = [f for f in sorted(glob.glob(os.path.join(run, "fdem_joints_[0-9]*.vtu")))
          if complete(f)]
    P, C, S = read_vtu(jn[-1], ["damage", "breakMode"], ncell=2)
    bm = S["breakMode"]
    seg = P[C]
    xm = seg[:, :, 0].mean(axis=1) - cx
    ym = seg[:, :, 1].mean(axis=1) - cy
    r = np.hypot(xm, ym)
    d = seg[:, 1] - seg[:, 0]
    ang = np.degrees(np.arctan2(d[:, 1], d[:, 0]))
    tang = np.degrees(np.arctan2(xm, -ym))        # tangente locale
    rel = (ang - tang + 90.0) % 180.0 - 90.0      # dans [-90, 90]
    zone = (r > rmin) & (r < rmax)
    L = np.linalg.norm(d, axis=1)
    return rel, zone, bm, L


def hist(rel, w, lab):
    h, edges = np.histogram(rel, bins=18, range=(-90, 90), weights=w)
    h = 100 * h / h.sum()
    ctr = 0.5 * (edges[1:] + edges[:-1])
    print(f"  {lab:26s} " + " ".join(f"{v:4.1f}" for v in h))
    return h, ctr


def main():
    runs = sys.argv[1:] or ["out_tun_ref_iso"]
    print(f"Angle des fissures de CISAILLEMENT avec la tangente locale, "
          f"secteurs de 10 deg de -90 a +90.")
    print(f"Attendu pour un facies conjugue de Mohr-Coulomb (phi = {PHI:g} deg) :"
          f" deux pics a +-{THEO:.0f} deg.\n")
    for rn in runs:
        run = rn if os.path.isabs(rn) else os.path.join(HERE, "..", "..", rn)
        rel, zone, bm, L = profile(run)
        print(f"--- {os.path.basename(rn)} ---")
        shear = zone & (bm == 2)
        offer = zone & (bm == 0)
        hs, ctr = hist(rel[shear], L[shear], "fissures cisaillement")
        ho, _ = hist(rel[offer], L[offer], "aretes NON rompues (offre)")
        sel = hs / np.maximum(ho, 1e-9)
        print("  " + " " * 26 + " " + " ".join(f"{v:4.2f}" for v in sel)
              + "   <- selection (rompu/offre)")
        # les deux secteurs les plus selectionnes
        k = np.argsort(sel)[::-1][:4]
        pk = sorted(ctr[k])
        print(f"  secteurs les plus selectionnes : "
              + ", ".join(f"{c:+.0f} deg (x{sel[list(ctr).index(c)]:.2f})"
                          for c in pk))
        # test symetrique : moyenne de la selection dans +-[20,40] contre
        # les zones proches de 0 et de +-90 (pur cisaillement tangentiel/radial)
        conj = (np.abs(ctr) > 18) & (np.abs(ctr) < 42)
        flat = (np.abs(ctr) < 10) | (np.abs(ctr) > 80)
        print(f"  selection moyenne dans la bande conjuguee +-[18,42] deg : "
              f"{sel[conj].mean():.2f}")
        print(f"  selection moyenne hors bande (0 et +-90)               : "
              f"{sel[flat].mean():.2f}")
        verdict = ("CONJUGUE" if sel[conj].mean() > 1.15 * sel[flat].mean()
                   else "pas de selection conjuguee nette")
        print(f"  -> {verdict}\n")


if __name__ == "__main__":
    main()
