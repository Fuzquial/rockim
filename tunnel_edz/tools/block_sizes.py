#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# block_sizes.py — la taille des BLOCS, mesuree.
#
#   python tunnel_edz/tools/block_sizes.py out_A out_B ... [--t 0.41] [--r 25]
#
# Pourquoi pas crack_clusters : dans tous nos runs le reseau de fissures est
# deja PERCOLANT (une composante connexe contient 96-99 % des joints rompus).
# Ce qui distingue le nuage des blocs n'est donc pas la connexite des
# fissures, c'est la taille des morceaux INTACTS qu'elles decoupent.
#
# Methode : graphe des elements (2 triangles voisins = une arete), on RETIRE
# les aretes dont le joint est rompu (D >= 0,999), et on prend les composantes
# connexes. Un bloc = un paquet d'elements encore solidaires. On pondere par
# l'aire pour parler en m2, et on restreint a un disque de rayon --r autour du
# tunnel (sinon le massif lointain intact ecrase la statistique).
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import os
import re

import numpy as np

CX = CY = 50.0


def arr(s, name, n=None):
    a = np.fromstring(re.search(
        r'Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % name, s, re.S).group(1),
        sep=" ")
    return a


def points(s):
    return np.fromstring(re.search(
        r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").reshape(-1, 3)[:, :2]


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


def analyse(run, tcible, rmax):
    mesh = sorted(f for f in glob.glob(run + "/fdem_[0-9]*.vtu")
                  if "joints" not in os.path.basename(f))
    jts = sorted(glob.glob(run + "/fdem_joints_[0-9]*.vtu"))
    if not mesh or not jts:
        print("%-20s : trames manquantes" % run)
        return
    ft = {}
    fcsv = os.path.join(run, "frames.csv")
    if os.path.exists(fcsv):
        for line in open(fcsv).read().splitlines()[1:]:
            p = line.split(",")
            ft[int(p[0])] = float(p[1])
    idx = len(jts) - 1
    if tcible and ft:
        cand = [(abs(v - tcible), k) for k, v in ft.items() if k < len(jts)]
        idx = min(cand)[1] if cand else idx

    s0 = io.open(mesh[0], encoding="utf-8", errors="ignore").read()
    P0 = points(s0)
    tri = arr(s0, "connectivity").astype(int).reshape(-1, 3)
    # noeuds uniques par coordonnees non deformees
    _, node = np.unique(np.round(P0, 6), axis=0, return_inverse=True)
    a, b, c = (P0[tri[:, i]] for i in range(3))
    area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                        - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
    ctr = (a + b + c) / 3.0
    rad = np.hypot(ctr[:, 0] - CX, ctr[:, 1] - CY)

    # arete (paire de noeuds uniques) -> elements adjacents
    edge = {}
    for e in range(len(tri)):
        n3 = node[tri[e]]
        for i in range(3):
            k = (min(n3[i], n3[(i + 1) % 3]), max(n3[i], n3[(i + 1) % 3]))
            edge.setdefault(k, []).append(e)

    sj0 = io.open(jts[0], encoding="utf-8", errors="ignore").read()
    Pj = points(sj0)
    conj = arr(sj0, "connectivity").astype(int).reshape(-1, 2)
    sjt = io.open(jts[idx], encoding="utf-8", errors="ignore").read()
    D = arr(sjt, "damage")
    nodej = np.array([np.argmin(np.sum((np.round(Pj[i], 6) - 0) ** 2))
                      for i in range(0)])  # placeholder (non utilise)
    # projection des extremites de joints sur la numerotation des noeuds :
    # meme arrondi, on retrouve l'indice unique par recherche dans un dict
    keymap = {}
    for i, k in enumerate(map(tuple, np.round(P0, 6))):
        keymap.setdefault(k, node[i])
    broken = set()
    miss = 0
    for j in np.where(D >= 0.999)[0]:
        k1 = keymap.get(tuple(np.round(Pj[conj[j, 0]], 6)))
        k2 = keymap.get(tuple(np.round(Pj[conj[j, 1]], 6)))
        if k1 is None or k2 is None:
            miss += 1
            continue
        broken.add((min(k1, k2), max(k1, k2)))

    uf = UF(len(tri))
    for k, els in edge.items():
        if len(els) == 2 and k not in broken:
            uf.union(els[0], els[1])
    root = np.array([uf.find(e) for e in range(len(tri))])

    sel = rad <= rmax
    _, lab = np.unique(root[sel], return_inverse=True)
    A = np.bincount(lab, weights=area[sel])
    A = np.sort(A)[::-1]
    n1 = (A < 2 * np.median(area[sel])).sum()      # blocs mono-element
    print("%-20s trame %2d (t = %.3f s) — disque r <= %.0f m"
          % (run, idx, ft.get(idx, float("nan")), rmax))
    print("   aretes rompues appariees : %d (%d non trouvees)"
          % (len(broken), miss))
    print("   blocs : %5d   | aire mediane %6.3f m2 | moyenne %6.3f m2"
          % (len(A), np.median(A), A.mean()))
    print("   5 plus gros blocs [m2] : " + " ".join("%.1f" % x for x in A[:5]))
    print("   blocs mono-element : %5d  (%4.1f %%)"
          % (n1, 100 * n1 / max(len(A), 1)))
    # taille caracteristique hors massif : on ecarte le plus gros (le massif)
    if len(A) > 1:
        Ab = A[1:]
        print("   hors massif : aire moyenne %6.3f m2, p95 %6.3f, max %6.2f"
              % (Ab.mean(), np.percentile(Ab, 95), Ab.max()))
        print("   -> cote equivalent moyen : %.2f m" % np.sqrt(Ab.mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", type=float, default=0.0)
    ap.add_argument("--r", type=float, default=25.0)
    a = ap.parse_args()
    for r in a.runs:
        analyse(r.rstrip("/\\"), a.t, a.r)


if __name__ == "__main__":
    main()
