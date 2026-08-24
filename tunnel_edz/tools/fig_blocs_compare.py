#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_blocs_compare.py — les trois schemas d'insertion, cote a cote, chaque
# element COLORE PAR LA TAILLE DU BLOC auquel il appartient.
#
#   python tunnel_edz/tools/fig_blocs_compare.py
#
# C'est le resultat central de l'etude rendu lisible : la ou l'anneau est
# GRANULE, tous les elements sont bleu fonce (bloc = 1 element) ; la ou de
# VRAIS BLOCS se detachent, ils apparaissent en clair, d'un seul tenant.
# ---------------------------------------------------------------------------
import glob
import io
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})

CX = CY = 50.0
RUNS = [("out_tun_ref_stab", 0.41, "adaptatif (reference)"),
        ("out_tun_tip16", 0.41, r"adaptatif + pointe relachee /1,6"),
        ("out_tun_intr", 0.54, "intrinsique (tous joints des t = 0)")]


def arr(s, name):
    return np.fromstring(re.search(
        r'Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % name, s, re.S).group(1),
        sep=" ")


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


def blocs(run, tcible):
    mesh = sorted(f for f in glob.glob(run + "/fdem_[0-9]*.vtu")
                  if "joints" not in os.path.basename(f))
    jts = sorted(glob.glob(run + "/fdem_joints_[0-9]*.vtu"))
    ft = {}
    for line in open(os.path.join(run, "frames.csv")).read().splitlines()[1:]:
        p = line.split(",")
        ft[int(p[0])] = float(p[1])
    idx = min([(abs(v - tcible), k) for k, v in ft.items()
               if k < len(jts)])[1]

    s0 = io.open(mesh[0], encoding="utf-8", errors="ignore").read()
    P0 = points(s0)
    tri = arr(s0, "connectivity").astype(int).reshape(-1, 3)
    _, node = np.unique(np.round(P0, 6), axis=0, return_inverse=True)
    a, b, c = (P0[tri[:, i]] for i in range(3))
    area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                        - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))

    edge = {}
    for e in range(len(tri)):
        n3 = node[tri[e]]
        for i in range(3):
            k = (min(n3[i], n3[(i + 1) % 3]), max(n3[i], n3[(i + 1) % 3]))
            edge.setdefault(k, []).append(e)

    sj0 = io.open(jts[0], encoding="utf-8", errors="ignore").read()
    Pj, conj = points(sj0), arr(sj0, "connectivity").astype(int).reshape(-1, 2)
    D = arr(io.open(jts[idx], encoding="utf-8", errors="ignore").read(),
            "damage")
    keymap = {}
    for i, k in enumerate(map(tuple, np.round(P0, 6))):
        keymap.setdefault(k, node[i])
    broken, seg = set(), []
    for j in np.where(D >= 0.999)[0]:
        k1 = keymap.get(tuple(np.round(Pj[conj[j, 0]], 6)))
        k2 = keymap.get(tuple(np.round(Pj[conj[j, 1]], 6)))
        if k1 is None or k2 is None:
            continue
        broken.add((min(k1, k2), max(k1, k2)))
        seg.append([Pj[conj[j, 0]] - [CX, CY], Pj[conj[j, 1]] - [CX, CY]])

    uf = UF(len(tri))
    for k, els in edge.items():
        if len(els) == 2 and k not in broken:
            uf.union(els[0], els[1])
    root = np.array([uf.find(e) for e in range(len(tri))])
    _, lab = np.unique(root, return_inverse=True)
    A = np.bincount(lab, weights=area)
    return P0, tri, A[lab], np.array(seg), ft[idx], len(broken)


fig, ax = plt.subplots(1, 3, figsize=(14.4, 5.3))
for A, (run, tc, titre) in zip(ax, RUNS):
    P0, tri, ablk, seg, t, nb = blocs(run, tc)
    x, y = P0[:, 0] - CX, P0[:, 1] - CY
    v = np.log10(np.clip(ablk, 1e-3, 1e3))
    tp = A.tripcolor(x, y, tri, v, cmap="viridis", vmin=-2, vmax=2,
                     shading="flat", rasterized=True)
    if len(seg):
        A.add_collection(LineCollection(seg, colors="k", linewidths=0.15,
                                        alpha=0.35))
    A.set_xlim(-26, 26)
    A.set_ylim(-24, 24)
    A.set_aspect("equal")
    A.set_xlabel("x [m]")
    A.set_title("%s\nt = %.2f s, %d joints rompus" % (titre, t, nb),
                fontsize=10)
    if A is ax[0]:
        A.set_ylabel("y [m]")
ax[0].figure.colorbar(tp, ax=ax, pad=0.015, shrink=0.85,
                      label=r"$\log_{10}$(aire du bloc [m$^2$])")
fig.suptitle("Taille des blocs intacts decoupes par la fissuration — "
             "bleu = un seul element (anneau granule), jaune = gros bloc",
             fontsize=12)
for ext in ("pdf", "png"):
    fig.savefig("tunnel_edz/fig_blocs_compare." + ext, dpi=165,
                bbox_inches="tight")
print("ecrit : tunnel_edz/fig_blocs_compare")
