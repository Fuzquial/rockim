#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# mesh_isotropy.py — le maillage est-il VRAIMENT non structuré ?
#
#   python tunnel_edz/tools/mesh_isotropy.py meshes/tunnel_hs.msh [--png sortie]
#
# Trois mesures, parce que « non structuré » ne se juge pas à l'oeil :
#   1. HISTOGRAMME D'ORIENTATION des arêtes (secteurs de 10 deg, pondéré par la
#      longueur). Un maillage isotrope donne 5,6 % par secteur. Un réseau
#      hexagonal parfait concentre tout sur 0/60/120 deg. L'indicateur retenu
#      est le rapport pic/creux et l'anisotropie R = |somme(exp(6 i theta))|/N,
#      qui vaut 0 pour l'isotropie et 1 pour un réseau hexagonal parfait
#      (l'ordre 6 est celui d'un pavage de triangles équilatéraux).
#   2. HISTOGRAMME DU PLUS PETIT ANGLE de chaque triangle. Un Delaunay
#      aléatoire étale de 20 a 60 deg ; un maillage quasi structuré se masse
#      contre 60 deg (triangles tous équilatéraux).
#   3. Un ZOOM d'image, pour voir.
#
# Motivation : sur le run tunnel du 2026-08-17, les fissures suivent des lignes
# droites sur plusieurs mètres. Il faut savoir si c'est la physique ou le
# mailleur.
# ---------------------------------------------------------------------------
import argparse
import os

import numpy as np


def read_msh22(path):
    xy, tri, idx = [], [], {}
    with open(path) as f:
        it = iter(f)
        for line in it:
            if line.startswith("$Nodes"):
                for _ in range(int(next(it))):
                    p = next(it).split()
                    idx[int(p[0])] = len(xy)
                    xy.append((float(p[1]), float(p[2])))
            elif line.startswith("$Elements"):
                for _ in range(int(next(it))):
                    p = next(it).split()
                    if int(p[1]) != 2:
                        continue
                    n = int(p[2])
                    tri.append([idx[int(v)] for v in p[3 + n:6 + n]])
    return np.array(xy), np.array(tri, dtype=int)


def stats(P, T, cx, cy, rmax, label):
    c = P[T].mean(axis=1)
    keep = np.hypot(c[:, 0] - cx, c[:, 1] - cy) < rmax        # zone fine
    T = T[keep]
    a = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    e = [a[1] - a[0], a[2] - a[1], a[0] - a[2]]
    ang = np.concatenate([np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 180
                          for v in e])
    ln = np.concatenate([np.linalg.norm(v, axis=1) for v in e])
    h, _ = np.histogram(ang, bins=18, range=(0, 180), weights=ln)
    h = 100 * h / h.sum()
    # ordre 6 : un pavage équilatéral a ses arêtes a 0, 60, 120 deg
    R = np.abs(np.sum(ln * np.exp(6j * np.radians(ang)))) / ln.sum()
    L = np.stack([np.linalg.norm(v, axis=1) for v in e], axis=1)
    cosA = [(L[:, (i + 1) % 3] ** 2 + L[:, (i + 2) % 3] ** 2 - L[:, i] ** 2)
            / (2 * L[:, (i + 1) % 3] * L[:, (i + 2) % 3]) for i in range(3)]
    A = np.degrees(np.arccos(np.clip(np.stack(cosA, axis=1), -1, 1)))
    amin = A.min(axis=1)
    print(f"--- {label} : {len(T)} triangles dans la zone fine ---")
    print("  orientation des aretes, secteurs de 10 deg (isotrope = 5,6 %) :")
    print("   " + " ".join(f"{v:4.1f}" for v in h))
    print(f"  pic/creux = {h.max() / h.min():.2f}   "
          f"anisotropie hexagonale R6 = {R:.3f}  "
          f"(0 = isotrope, 1 = reseau parfait)")
    print(f"  plus petit angle : median {np.median(amin):.1f} deg, "
          f"p5 {np.percentile(amin, 5):.1f}, min {amin.min():.1f} ; "
          f"{100 * np.mean(amin > 50):.1f} % des triangles au-dessus de 50 deg")
    return R, h.max() / h.min(), np.median(amin)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("--cx", type=float, default=50.0)
    ap.add_argument("--cy", type=float, default=50.0)
    ap.add_argument("--rmax", type=float, default=20.0)
    ap.add_argument("--png", default=None)
    a = ap.parse_args()
    P, T = read_msh22(a.mesh)
    stats(P, T, a.cx, a.cy, a.rmax, os.path.basename(a.mesh))
    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
        fig, ax = plt.subplots(1, 2, figsize=(13, 6.4))
        for k, w in enumerate((8.0, 3.0)):
            ax[k].triplot(mtri.Triangulation(P[:, 0], P[:, 1], T), lw=0.45,
                          color="0.25")
            ax[k].set_xlim(a.cx - w, a.cx + w)
            ax[k].set_ylim(a.cy + 4, a.cy + 4 + 2 * w)
            ax[k].set_aspect("equal")
            ax[k].set_title(f"{os.path.basename(a.mesh)} — fenetre {2*w:.0f} m "
                            "au-dessus de la voute")
        fig.tight_layout()
        fig.savefig(a.png, dpi=170)
        print("ecrit :", a.png)


if __name__ == "__main__":
    main()
