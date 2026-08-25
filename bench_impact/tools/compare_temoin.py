#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# compare_temoin.py — le depouillement du TEMOIN du point 1 : adaptatif contre
# intrinseque, sur deux runs APPARIES qui ne different que par le schema.
#
#   python bench_impact/tools/compare_temoin.py out_temoin_adap out_temoin_intr
#
# Reutilise imp_lib (les 7 criteres de leur Table 3) et AJOUTE ce que
# fig_impact.py ne fait pas : le BILAN D ENERGIE confronte aux chiffres publies.
#
# LA question de ce temoin n est pas la fissuration — elle est deja bonne. Le
# papier ARMA 2024 (Yang, Xiang, Naderi, Wang, Latham, Aising, Gerbaud,
# Ugarte, ARMA 24-0952) mesure, sur St Anne a 9,41 m/s pour 49,3 J entrants :
#
#     energie de FISSURATION   1,3 a 1,66 J   =  2,6 %
#     energie de FROTTEMENT   32,0 J          = 64,9 %
#     amortissement + erreur   ~ +7 J transitoire, ~0 en fin
#
# Le run adaptatif du 2026-08-22 donnait eJnt = 1,10 J — la BONNE valeur — et
# eFric = 0,66 J contre 32. C est le frottement entre fragments qui manque,
# d un facteur ~50, pas la surface de fissure. Ce script met ce rapport au
# centre du tableau.
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import (broken, frames_of, history, joints_frame, metrics,
                     read_vtu)

# Cibles publiees. Les 7 criteres viennent de la Table 3 / fig. 9-10 de
# IJRMMS 191 (2025) a 10,66 m/s ; les energies de l ARMA 2024 a 9,41 m/s.
CIBLES = {
    "szz":    ("contrainte max au bit",      "MPa", "200 - 260"),
    "vind":   ("vitesse d indentation",      "m/s", "9,40 - 9,85"),
    "vreb":   ("vitesse de rebond",          "m/s", "6,87 - 7,10"),
    "ratio":  ("rebond / indentation",       "-",   "0,72 - 0,73"),
    "depth":  ("profondeur d indentation",   "mm",  "~1,53"),
    "radial": ("fissure radiale max",        "mm",  "20,2 - 24,5"),
    "crater": ("rayon de cratere",           "mm",  "10,0 - 12,1"),
}
E_FISSURATION_CIBLE = 1.3      # J, ARMA 2024 fig. 4, fin de run
E_FROTTEMENT_CIBLE = 32.0      # J, idem


def depouille(run):
    """Les 7 criteres + le bilan d energie d un run."""
    h = history(run)
    ks = frames_of(run)
    pts, con, f = read_vtu(joints_frame(run, ks[-1]))
    c, n, mode, P = broken(pts, con, f)
    m = metrics(c)

    vz = h["vz_bit"]
    vind = -vz.min()
    apres = int(np.argmin(vz))
    vreb = vz[apres:].max()
    zi = h["z_insert"]

    d = dict(
        szz=(np.abs(h["szz_bit"]).max() / 1e6) if "szz_bit" in h else np.nan,
        vind=vind,
        vreb=vreb,
        ratio=(vreb / vind) if vind > 0 else np.nan,
        depth=(zi[0] - zi.min()) * 1e3,
        radial=m["radial"] * 1e3,
        crater=m["crater"] * 1e3,
        nbroken=m["n"],
        t_fin=h["t"][-1] * 1e6,
    )
    # Bilan d energie V2/B4 : les colonnes sont des travaux CUMULES signes,
    # negatif = preleve au solide. On les rend positifs pour la lecture.
    for k in ("eEl", "eJnt", "eGc", "eFric", "eCund", "eLys"):
        d[k] = -h[k][-1] if k in h else np.nan
    return d


def ligne(nom, unite, a, b, cible):
    """Une ligne du tableau : adaptatif | intrinseque | cible."""
    def fmt(x):
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return "     —"
        return "%9.3g" % x if abs(x) < 1e-2 else "%9.2f" % x
    return "  %-26s %-5s %s %s   %s" % (nom, unite, fmt(a), fmt(b), cible)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("adaptatif")
    ap.add_argument("intrinseque")
    a = ap.parse_args()

    A = depouille(a.adaptatif)
    I = depouille(a.intrinseque)

    print()
    print("=" * 78)
    print("  TEMOIN DU POINT 1 — schema d insertion, tout le reste egal")
    print("  adaptatif : %s   (t = %.0f us)" % (a.adaptatif, A["t_fin"]))
    print("  intrinseque: %s  (t = %.0f us)" % (a.intrinseque, I["t_fin"]))
    print("=" * 78)
    print()
    print("  %-26s %-5s %9s %9s   %s"
          % ("critere", "unite", "ADAPT", "INTRIN", "publie"))
    print("  " + "-" * 74)
    for k, (nom, unite, cible) in CIBLES.items():
        print(ligne(nom, unite, A[k], I[k], cible))
    print(ligne("joints rompus", "-", A["nbroken"], I["nbroken"], ""))

    print()
    print("  BILAN D ENERGIE — la vraie question du temoin")
    print("  " + "-" * 74)
    print(ligne("fissuration (eJnt)", "J", A["eJnt"], I["eJnt"],
                "%.1f  (2,6 %% de 49,3)" % E_FISSURATION_CIBLE))
    print(ligne("frottement (eFric)", "J", A["eFric"], I["eFric"],
                "%.1f  (64,9 %%)" % E_FROTTEMENT_CIBLE))
    print(ligne("contact general (eGc)", "J", A["eGc"], I["eGc"], ""))
    print(ligne("elastique stocke (eEl)", "J", A["eEl"], I["eEl"], ""))
    print(ligne("Cundall (eCund)", "J", A["eCund"], I["eCund"], "0 (coupe)"))
    print()
    for nom, d in (("adaptatif", A), ("intrinseque", I)):
        if d["eFric"] and not np.isnan(d["eFric"]) and d["eFric"] > 0:
            print("    %-12s deficit de frottement : facteur %.0f"
                  % (nom, E_FROTTEMENT_CIBLE / d["eFric"]))
    print()
    print("  LECTURE. Si l intrinseque comble le deficit de frottement, le")
    print("  schema d insertion EST la cause et le point 1 se ferme. S il ne le")
    print("  comble pas, la cause est ailleurs — et le tableau designe alors")
    print("  les points 2 et 3 (le joint rompu ne meurt jamais en compression,")
    print("  donc le contact ne prend jamais le relais avec son mu = 0,6).")
    print()


if __name__ == "__main__":
    main()
