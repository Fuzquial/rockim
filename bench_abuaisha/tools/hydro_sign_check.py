#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# hydro_sign_check.py — LE CONTROLE DE SIGNE DU MODULE HYDRO.
#
#   python bench_abuaisha/tools/hydro_sign_check.py out_signe_conf out_signe_hydro
#
# H3 transpose au forage (cf. bench_abuaisha/configs/signe_*.cfg). Les deux
# chemins chargent la MEME liste de faces avec la MEME rampe : sans rupture,
# ils doivent donner le meme deplacement de paroi, et ce deplacement doit etre
# VERS L'EXTERIEUR. La cible est analytique — trou pressurise en deformation
# plane, u_r = p a / 2G.
#
# Ce controle existe parce que H3 (parker_compare.py) mesurait une valeur
# ABSOLUE et etait aveugle au signe. Ici le signe du rayon ne se cache pas.
# ---------------------------------------------------------------------------
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

CX, CY, RB = 4.0, 4.0, 0.05
E, NU, P = 35.0e9, 0.27, 12.0e6


def wall_dr(run):
    fs = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    if not fs:
        raise SystemExit("aucune trame complete dans " + run)
    P0, _, _ = read_vtu(fs[0], [])
    Pf, _, _ = read_vtu(fs[-1], [])
    r0 = np.hypot(P0[:, 0] - CX, P0[:, 1] - CY)
    w = np.where(np.abs(r0 - RB) < 2.5e-3)[0]
    rf = np.hypot(Pf[w, 0] - CX, Pf[w, 1] - CY)
    return (rf - r0[w]).mean(), len(w), len(fs)


def main():
    a, b = sys.argv[1], sys.argv[2]
    G = E / (2.0 * (1.0 + NU))
    ur = P * RB / (2.0 * G)
    da, na, fa = wall_dr(a)
    db, nb, fb = wall_dr(b)
    print("cible analytique (Lame, u_r = p a / 2G) : %+9.3f um" % (ur * 1e6))
    print("  %-18s %+9.3f um   (%d noeuds de paroi, %d trames)"
          % ("confinement", da * 1e6, na, fa))
    print("  %-18s %+9.3f um   (%d noeuds de paroi, %d trames)"
          % ("hydro", db * 1e6, nb, fb))
    ok = True
    if db <= 0.0:
        print("ECHEC : la paroi RENTRE — le fluide serre la cavite au lieu "
              "de l'ouvrir. Le signe du chargement est inverse.")
        ok = False
    ecart = abs(db - da) / max(abs(da), 1e-30) * 100.0
    print("  ecart entre chemins : %.3e %%" % ecart)
    if ecart > 1e-6:
        print("ECHEC : les deux chemins divergent alors qu'ils chargent les "
              "memes faces avec la meme rampe.")
        ok = False
    err = (db - ur) / ur * 100.0
    print("  ecart a l analytique : %+.2f %%" % err)
    if abs(err) > 15.0:
        print("ECHEC : trop loin de Lame.")
        ok = False
    print("[PASS]" if ok else "[FAIL]")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
