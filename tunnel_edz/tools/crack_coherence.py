#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# crack_coherence.py — les fissures vont-elles DROIT ?
#
#   python tunnel_edz/tools/crack_coherence.py out_A out_B ... [--t 0.41]
#
# Une macro-fracture est une chaine d aretes rompues QUASI COLINEAIRES. Un
# nuage est une collection d aretes rompues d orientations quelconques. On
# mesure donc, pour chaque paire d aretes rompues partageant un sommet,
# l angle entre leurs directions (0 deg = parfaitement alignees, 90 = en T).
#
# On mesure aussi le BRANCHEMENT : le nombre d aretes rompues par sommet.
# Une fissure propre en porte 2 (elle traverse) ; 3 et plus = bifurcation.
#
# Topologie et directions prises sur la configuration NON deformee.
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import os
import re

import numpy as np


def jvtu(path, dmg=True):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    pts = np.fromstring(re.search(
        r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").reshape(-1, 3)[:, :2]
    con = np.fromstring(re.search(
        r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").astype(int).reshape(-1, 2)
    D = None
    if dmg:
        D = np.fromstring(re.search(
            r'Name="damage"[^>]*>\s*(.*?)\s*</DataArray>', s,
            re.S).group(1), sep=" ")
    return pts, con, D


def analyse(run, tcible):
    fs = sorted(glob.glob(run + "/fdem_joints_[0-9]*.vtu"))
    if not fs:
        print("%-20s : aucun VTU joints" % run)
        return
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

    P0, con, _ = jvtu(fs[0], dmg=False)
    _, _, D = jvtu(fs[idx])
    _, node = np.unique(np.round(P0, 6), axis=0, return_inverse=True)
    brk = np.where(D >= 0.999)[0]
    if len(brk) < 10:
        print("%-20s : trop peu de rompus" % run)
        return

    # direction unitaire de chaque arete rompue (non orientee)
    v = P0[con[brk, 1]] - P0[con[brk, 0]]
    v /= np.linalg.norm(v, axis=1)[:, None]

    # aretes rompues incidentes a chaque sommet
    inc = {}
    for k, j in enumerate(brk):
        for nd in (node[con[j, 0]], node[con[j, 1]]):
            inc.setdefault(nd, []).append(k)

    ang = []
    deg = []
    for nd, ks in inc.items():
        deg.append(len(ks))
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                c = abs(float(np.dot(v[ks[i]], v[ks[j]])))
                ang.append(np.degrees(np.arccos(min(1.0, c))))
    ang = np.array(ang)
    deg = np.array(deg)

    print("%-20s trame %2d (t = %.3f s) — %d aretes rompues"
          % (run, idx, ft.get(idx, float("nan")), len(brk)))
    print("   angle entre aretes voisines : median %5.1f deg | "
          "moyenne %5.1f | p25 %5.1f" % (np.median(ang), ang.mean(),
                                         np.percentile(ang, 25)))
    print("   paires QUASI ALIGNEES (< 30 deg) : %5.1f %%"
          % (100.0 * (ang < 30).mean()))
    print("   paires en T (> 60 deg)          : %5.1f %%"
          % (100.0 * (ang > 60).mean()))
    print("   sommets : degre moyen %.2f | %5.1f %% de degre >= 3 "
          "(bifurcations)" % (deg.mean(), 100.0 * (deg >= 3).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", type=float, default=0.0)
    a = ap.parse_args()
    for r in a.runs:
        analyse(r.rstrip("/\\"), a.t)


if __name__ == "__main__":
    main()
