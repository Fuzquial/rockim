#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_all_meshes.py — planche des maillages des simulations preparees.
# Vue d'ensemble + zoom sur la zone raffinee pour chacun.
# ---------------------------------------------------------------------------
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "cm", "font.size": 8.5})

P1 = "../rockim/rockim_p1/meshes/"
MESHES = [
    ("meshes/tunnel_hs_iso.msh", "tunnel_hs_iso",
     "3 decks schistosite\n(tunnel_wp00/45/90)", (50, 50), 11.0),
    (P1 + "hf_bore.msh", "hf_bore",
     "9 decks banc AbuAisha\n(e1-e9, e11, e13)", None, 0.11),
    (P1 + "hf_bore_w.msh", "hf_bore_w",
     "1 deck (e10_plateau_aniso)\nzone raffinee 1,0 m", None, 0.11),
    (P1 + "montney.msh", "montney",
     "1 deck (e12_montney)\nessai de terrain", None, 0.16),
]


def read_msh(path):
    lines = open(path, errors="ignore").read().split("\n")
    i = lines.index("$Nodes")
    nn = int(lines[i + 1])
    xy = np.empty((nn, 2))
    for k in range(nn):
        p = lines[i + 2 + k].split()
        xy[k] = (float(p[1]), float(p[2]))
    j = lines.index("$Elements")
    ne = int(lines[j + 1])
    tris = []
    for k in range(ne):
        p = lines[j + 2 + k].split()
        if p[1] == "2":
            nt = int(p[2])
            tris.append([int(v) - 1 for v in p[3 + nt:6 + nt]])
    return xy, np.array(tris)


fig, axes = plt.subplots(2, 4, figsize=(15.0, 7.6))

for col, (path, name, use, ctr, half) in enumerate(MESHES):
    try:
        xy, tris = read_msh(path)
    except Exception as e:
        for r in (0, 1):
            axes[r][col].text(0.5, 0.5, "ABSENT", ha="center", va="center")
            axes[r][col].axis("off")
        print(f"  {name}: ABSENT ({e})")
        continue
    # arete mediane, pour le rapport
    e0 = np.linalg.norm(xy[tris[:, 1]] - xy[tris[:, 0]], axis=1)
    hmed = float(np.median(e0))
    if ctr is None:
        ctr = (float(xy[:, 0].mean()), float(xy[:, 1].mean()))
    W = xy[:, 0].max() - xy[:, 0].min()
    H = xy[:, 1].max() - xy[:, 1].min()
    print(f"  {name}: {len(xy)} noeuds, {len(tris)} triangles, "
          f"{W:.1f} x {H:.1f} m, h_med = {hmed*1000:.1f} mm")

    for row, hf in enumerate([None, half]):
        ax = axes[row][col]
        if hf is None:
            sel = np.ones(len(tris), bool)
            x0, x1 = xy[:, 0].min(), xy[:, 0].max()
            y0, y1 = xy[:, 1].min(), xy[:, 1].max()
            lw, ttl = 0.05, f"{name}\n{len(tris):,} triangles".replace(",", " ")
        else:
            c = xy[tris].mean(axis=1)
            sel = ((np.abs(c[:, 0] - ctr[0]) < hf) &
                   (np.abs(c[:, 1] - ctr[1]) < hf))
            x0, x1 = ctr[0] - hf, ctr[0] + hf
            y0, y1 = ctr[1] - hf, ctr[1] + hf
            lw, ttl = 0.22, f"zoom  ($h_{{med}}$ = {hmed*1000:.1f} mm)"
        ax.add_collection(PolyCollection(
            xy[tris[sel]], facecolors="none", edgecolors="0.55",
            linewidths=lw, rasterized=True))
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_title(ttl, fontsize=8.5, pad=5)
        ax.tick_params(labelsize=7)
        if row == 1:
            ax.set_xlabel(use, fontsize=8, labelpad=6)

fig.suptitle("Maillages des simulations preparees — aucune n'est lancee",
             fontsize=12, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("tunnel_schisto/meshes_planche.png", dpi=155)
fig.savefig("tunnel_schisto/meshes_planche.pdf", dpi=155)
print("  ecrit : tunnel_schisto/meshes_planche.png / .pdf")
