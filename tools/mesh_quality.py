#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# mesh_quality.py — les SLIVERS d'un maillage .msh v2.2 : diametre inscrit
# h = 6 V / (somme des aires des faces) par tetra (la mesure que rockim
# appelle « h inscrit » et qui commande le pas de temps), par corps physique.
#   python tools/mesh_quality.py meshes/a.msh [meshes/b.msh ...] [--worst 5]
# Imprime : N, h min, quantiles 0,1 % / 1 %, nombre de tetras sous 0,3 / 0,5 mm,
# et les pires elements (position, h, arete moyenne, corps).
# ---------------------------------------------------------------------------
import sys
import numpy as np


def read_msh(path):
    names, nodes, tets, phys = {}, None, [], []
    with open(path, errors="replace") as f:
        lines = f.read().split(chr(10))
    i = 0
    while i < len(lines):
        L = lines[i].strip()
        if L == "$PhysicalNames":
            n = int(lines[i + 1])
            for k in range(n):
                d, tag, nm = lines[i + 2 + k].split(maxsplit=2)
                names[int(tag)] = nm.strip().strip('"')
            i += n + 3
        elif L == "$Nodes":
            n = int(lines[i + 1])
            arr = np.array(" ".join(lines[i + 2:i + 2 + n]).split(), float).reshape(n, 4)
            nodes = np.zeros((int(arr[:, 0].max()) + 1, 3))
            nodes[arr[:, 0].astype(int)] = arr[:, 1:]
            i += n + 3
        elif L == "$Elements":
            n = int(lines[i + 1])
            for k in range(n):
                p = lines[i + 2 + k].split()
                if p[1] == "4":
                    nt = int(p[2])
                    phys.append(int(p[3]))
                    tets.append([int(x) for x in p[3 + nt:7 + nt]])
            i += n + 3
        else:
            i += 1
    return nodes, np.array(tets), np.array(phys), names


def inscribed(nodes, tets):
    P = nodes[tets]                                   # (n, 4, 3)
    a, b, c, d = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    V = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0
    A = (np.linalg.norm(np.cross(b - a, c - a), axis=1) + np.linalg.norm(np.cross(b - a, d - a), axis=1)
         + np.linalg.norm(np.cross(c - a, d - a), axis=1) + np.linalg.norm(np.cross(c - b, d - b), axis=1)) / 2.0
    h = 6.0 * V / A
    E = np.stack([np.linalg.norm(b - a, axis=1), np.linalg.norm(c - a, axis=1), np.linalg.norm(d - a, axis=1),
                  np.linalg.norm(c - b, axis=1), np.linalg.norm(d - b, axis=1), np.linalg.norm(d - c, axis=1)], 1)
    return h, E.mean(1), P.mean(1), V


def report(path, worst=5):
    nodes, tets, phys, names = read_msh(path)
    h, em, cen, V = inscribed(nodes, tets)
    print("== %s : %d tetras, h min %.4f mm, q0.1%% %.4f mm, q1%% %.4f mm, mediane %.3f mm ; < 0,3 mm : %d ; < 0,5 mm : %d"
          % (path, len(h), h.min() * 1e3, np.quantile(h, 1e-3) * 1e3, np.quantile(h, 1e-2) * 1e3,
             np.median(h) * 1e3, (h < 3e-4).sum(), (h < 5e-4).sum()))
    for tag in sorted(set(phys.tolist())):
        m = phys == tag
        hh = h[m]
        print("   %-8s N %7d  h min %.4f mm  q1%% %.4f  mediane %.3f  arete med. %.3f mm  h/arete min %.3f  < 0,3 mm : %d"
              % (names.get(tag, str(tag)), m.sum(), hh.min() * 1e3, np.quantile(hh, 1e-2) * 1e3,
                 np.median(hh) * 1e3, np.median(em[m]) * 1e3, (hh / em[m]).min(), (hh < 3e-4).sum()))
    idx = np.argsort(h)[:worst]
    for j in idx:
        print("   pire : h %.4f mm  arete moy. %.3f mm  h/arete %.3f  corps %-7s  x %+.4f y %+.4f z %+.4f  (r %.4f)"
              % (h[j] * 1e3, em[j] * 1e3, h[j] / em[j], names.get(int(phys[j]), "?"),
                 cen[j, 0], cen[j, 1], cen[j, 2], np.hypot(cen[j, 0], cen[j, 1])))
    return h


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    worst = 5
    if "--worst" in sys.argv:
        worst = int(sys.argv[sys.argv.index("--worst") + 1])
        args = [a for a in args if a != str(worst)]
    for p in args:
        report(p, worst)
