#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# wall_convergence.py — la convergence de la paroi, mesuree ROBUSTEMENT.
#
#   python tunnel_edz/tools/wall_convergence.py out_tun_s3 out_tun_s4 ...
#
# POURQUOI. `edz_metrics` rapporte le deplacement MAXIMAL sur tous les noeuds.
# Un seul bloc detache qui part en vol suffit a le fixer : la mesure n'est
# alors plus la convergence du tunnel mais la trajectoire d'un debris. Le
# balayage du 2026-08-17 a produit un resultat non monotone (0,389 m a 5 MPa
# contre 0,204 m a 6 MPa) qui sent exactement ce piege.
#
# CE QU'ON MESURE ICI, sur les seuls noeuds de la PAROI (a moins de 0,4 m du
# contour initial) et en projetant sur la direction RADIALE (positif = vers
# l'interieur du tunnel) :
#   moyenne   : la convergence d'ensemble, insensible a un bloc isole
#   p90       : la convergence des zones les plus sollicitees
#   max       : pour comparaison avec la mesure fragile
#   fraction  : part des noeuds de paroi qui depassent 2 x la moyenne
#               (un chiffre eleve = deformation localisee ; tres eleve = debris)
# ---------------------------------------------------------------------------
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402
from plot_tunnel_mesh import profile_xy  # noqa: E402
from make_unstructured_mesh import TUNNEL_HS  # noqa: E402


def wall_nodes(P0, cx, cy, tol=0.4):
    """Noeuds initialement sur le contour de la cavite."""
    px, py, _ = profile_xy(cx, cy - 0.5 * TUNNEL_HS["height"], n=1200)
    C = np.stack([px, py], axis=1)
    # distance au contour par recherche du point de contour le plus proche
    sel = np.hypot(P0[:, 0] - cx, P0[:, 1] - cy) < 9.0        # pre-filtre
    idx = np.where(sel)[0]
    d = np.full(len(idx), 1e9)
    for k in range(0, len(C), 4):                             # 1 point sur 4
        d = np.minimum(d, np.hypot(P0[idx, 0] - C[k, 0], P0[idx, 1] - C[k, 1]))
    return idx[d < tol]


def convergence(run):
    """(moyenne, p90, max, % de noeuds au-dela de 2x la moyenne), en metres.

    Utilisable comme fonction : c'est ce que trace plot_sweep.py.
    """
    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    if len(el) < 2:
        return None
    P0, _, _ = read_vtu(el[0], [])
    P, _, _ = read_vtu(el[-1], [])
    cx, cy = 0.5 * P0[:, 0].max(), 0.5 * P0[:, 1].max()
    w = wall_nodes(P0, cx, cy)
    u = P[w] - P0[w]
    r = np.stack([P0[w, 0] - cx, P0[w, 1] - cy], axis=1)
    r /= np.linalg.norm(r, axis=1)[:, None]
    conv = -(u * r).sum(axis=1)
    return (conv.mean(), np.percentile(conv, 90), conv.max(),
            100.0 * np.mean(conv > 2.0 * conv.mean()))


def main():
    runs = sys.argv[1:]
    print(f"{'run':18s} {'moyenne':>9s} {'p90':>8s} {'max':>8s} "
          f"{'noeuds':>7s} {'>2x moy':>8s}")
    for rn in runs:
        run = rn if os.path.isabs(rn) else os.path.join(HERE, "..", "..", rn)
        el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
              if complete(f)]
        if len(el) < 2:
            print(f"{os.path.basename(rn):18s}  (pas assez de trames)")
            continue
        P0, _, _ = read_vtu(el[0], [])
        P, _, _ = read_vtu(el[-1], [])
        W, H = P0[:, 0].max(), P0[:, 1].max()
        cx, cy = 0.5 * W, 0.5 * H
        w = wall_nodes(P0, cx, cy)
        u = P[w] - P0[w]
        r = np.stack([P0[w, 0] - cx, P0[w, 1] - cy], axis=1)
        r /= np.linalg.norm(r, axis=1)[:, None]
        conv = -(u * r).sum(axis=1)          # positif = vers l'interieur
        frac = 100.0 * np.mean(conv > 2.0 * conv.mean())
        print(f"{os.path.basename(rn):18s} {conv.mean():8.4f} m "
              f"{np.percentile(conv, 90):7.4f} {conv.max():7.4f} "
              f"{len(w):7d} {frac:7.1f} %")


if __name__ == "__main__":
    main()
