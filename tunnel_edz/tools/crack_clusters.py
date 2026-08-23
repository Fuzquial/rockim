#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# crack_clusters.py — la COALESCENCE, mesuree.
#
#   python tunnel_edz/tools/crack_clusters.py out_A out_B ... [--t 0.41]
#
# Compter les joints rompus ne dit rien du MOTIF : 28 000 fissures isolees et
# 28 000 fissures chainees donnent le meme nombre. Ici on agrege les joints
# rompus en COMPOSANTES CONNEXES (deux joints rompus qui partagent un noeud
# appartiennent a la meme fissure) et on regarde la distribution des tailles.
#
# Localisation = peu de composantes, tres longues.  Diffusion = beaucoup de
# composantes courtes.  C'est l'observable qui distingue le nuage des blocs.
#
# Topologie prise sur la trame 0 (configuration NON deformee, ou les copies
# co-localisees d'un meme noeud coincident encore) ; endommagement pris a la
# trame demandee. L'ordre des joints est stable d'une trame a l'autre.
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import os
import re

import numpy as np


def read_joints(path, want_damage=True):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    pts = np.fromstring(re.search(
        r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").reshape(-1, 3)[:, :2]
    con = np.fromstring(re.search(
        r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").astype(int).reshape(-1, 2)
    dmg = None
    if want_damage:
        dmg = np.fromstring(re.search(
            r'Name="damage"[^>]*>\s*(.*?)\s*</DataArray>', s,
            re.S).group(1), sep=" ")
    return pts, con, dmg


class UF:
    def __init__(self, n):
        self.p = np.arange(n)

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def analyse(run, tcible):
    fs = sorted(glob.glob(run + "/fdem_joints_[0-9]*.vtu"))
    if not fs:
        print("%-20s : aucun VTU joints" % run)
        return
    # trame la plus proche du temps cible
    ft = {}
    fcsv = os.path.join(run, "frames.csv")
    if os.path.exists(fcsv):
        for line in open(fcsv).read().splitlines()[1:]:
            p = line.split(",")
            ft[int(p[0])] = float(p[1])
    idx = len(fs) - 1
    if tcible and ft:
        cand = [(abs(v - tcible), k) for k, v in ft.items() if k < len(fs)]
        idx = min(cand)[1] if cand else idx
    P0, con, _ = read_joints(fs[0], want_damage=False)
    _, _, D = read_joints(fs[idx])

    # noeuds uniques par coordonnees NON deformees
    key = np.round(P0, 6)
    _, node = np.unique(key, axis=0, return_inverse=True)
    brk = np.where(D >= 0.999)[0]
    uf = UF(node.max() + 1)
    for j in brk:
        uf.union(node[con[j, 0]], node[con[j, 1]])
    root = np.array([uf.find(node[con[j, 0]]) for j in brk])
    _, lab = np.unique(root, return_inverse=True)
    ncl = lab.max() + 1 if len(lab) else 0
    taille = np.bincount(lab) if len(lab) else np.array([0])

    # longueur geometrique de chaque joint rompu (config non deformee)
    L = np.linalg.norm(P0[con[brk, 0]] - P0[con[brk, 1]], axis=1)
    lcl = np.bincount(lab, weights=L) if len(lab) else np.array([0.0])

    print("%-20s trame %2d (t = %.3f s)" % (run, idx, ft.get(idx, float("nan"))))
    print("   joints rompus      : %6d" % len(brk))
    print("   fissures (composantes connexes) : %5d" % ncl)
    print("   joints par fissure : moyenne %5.1f | mediane %3.0f | max %5d"
          % (taille.mean(), np.median(taille), taille.max()))
    print("   longueur de fissure [m] : moyenne %5.2f | p95 %5.2f | max %6.2f"
          % (lcl.mean(), np.percentile(lcl, 95), lcl.max()))
    gros = np.sort(lcl)[::-1][:10].sum()
    print("   part des 10 plus grandes fissures : %4.1f %% de la longueur"
          % (100 * gros / lcl.sum()))
    iso = (taille <= 2).sum()
    print("   fissures de 1-2 joints (isolees)  : %5d  (%4.1f %%)"
          % (iso, 100 * iso / max(ncl, 1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", type=float, default=0.0,
                    help="temps cible [s] (0 = derniere trame)")
    a = ap.parse_args()
    for r in a.runs:
        analyse(r.rstrip("/\\"), a.t)


if __name__ == "__main__":
    main()
