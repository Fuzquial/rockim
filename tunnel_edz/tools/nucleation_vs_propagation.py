#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# nucleation_vs_propagation.py — le mecanisme, mesure.
#
#   python tunnel_edz/tools/nucleation_vs_propagation.py out_A out_B ...
#
# Une fissure GRANDIT si le joint qui casse touche un joint deja casse
# (propagation en pointe) ; elle NUCLEE s il casse en terrain vierge. Le
# rapport des deux distingue un schema qui propage d un schema qui essaime.
#
# Methode : entre deux trames consecutives, on releve les joints qui
# franchissent D = 0,999 ; pour chacun on regarde s il partage un noeud
# (numerotation de la configuration NON deformee) avec un joint deja rompu a
# la trame precedente. Le premier front (trame ou tout est nouveau) est
# exclu du ratio.
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


def analyse(run):
    fs = sorted(glob.glob(run + "/fdem_joints_[0-9]*.vtu"))
    if len(fs) < 3:
        print("%-20s : trop peu de trames" % run)
        return
    P0, con, _ = jvtu(fs[0], dmg=False)
    _, node = np.unique(np.round(P0, 6), axis=0, return_inverse=True)
    n1, n2 = node[con[:, 0]], node[con[:, 1]]

    ft = {}
    fcsv = os.path.join(run, "frames.csv")
    if os.path.exists(fcsv):
        for line in open(fcsv).read().splitlines()[1:]:
            p = line.split(",")
            ft[int(p[0])] = float(p[1])

    prev = np.zeros(len(con), bool)
    touche = np.zeros(node.max() + 1, bool)      # noeuds portes par du rompu
    nprop = nnuc = 0
    print("%-20s (%d trames)" % (run, len(fs)))
    print("   trame   t [s]   nouveaux   propagation   nucleation   %prop")
    for i, f in enumerate(fs):
        _, _, D = jvtu(f)
        now = D >= 0.999
        new = now & ~prev
        k = np.where(new)[0]
        if i > 1 and len(k):                     # trame 0/1 = amorcage
            adj = touche[n1[k]] | touche[n2[k]]
            p, n = int(adj.sum()), int((~adj).sum())
            nprop += p
            nnuc += n
            print("   %5d  %6.3f   %7d   %11d   %10d   %5.1f %%"
                  % (i, ft.get(i, float("nan")), len(k), p, n,
                     100.0 * p / max(len(k), 1)))
        touche[n1[now]] = True
        touche[n2[now]] = True
        prev = now
    tot = nprop + nnuc
    print("   ---> BILAN : %d propagations / %d nucleations "
          "-> %.1f %% de propagation" % (nprop, nnuc,
                                         100.0 * nprop / max(tot, 1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    for r in a.runs:
        analyse(r.rstrip("/\\"))


if __name__ == "__main__":
    main()
