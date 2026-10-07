#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_mesh.py — figure du maillage tunnel + des plans de schistosite REELS.
#
# Les plans ne sont PAS recalcules ici : ils sont lus dans le champ
# `weakPlane` que le solveur ecrit dans fdem_joints_0000.vtu. La figure montre
# donc ce que rockim a effectivement selectionne, et non ce qu'un script
# independant croit qu'il aurait du selectionner.
#
#   python tunnel_schisto/plot_mesh.py meshes/tunnel_hs_iso.msh out_dry45
# ---------------------------------------------------------------------------
import sys, re, base64, zlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

plt.rcParams.update({"font.family": "serif", "font.serif": ["CMU Serif", "DejaVu Serif"],
                     "mathtext.fontset": "cm", "font.size": 9})

msh, outdir = sys.argv[1], sys.argv[2]

# ---- maillage Gmsh 2.2 ASCII ----------------------------------------------
lines = open(msh, errors="ignore").read().split("\n")
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
    if p[1] == "2":                       # 2 = triangle 3 noeuds
        ntag = int(p[2])
        tris.append([int(v) - 1 for v in p[3 + ntag:6 + ntag]])
tris = np.array(tris)
print(f"  maillage : {nn} noeuds, {len(tris)} triangles")

# ---- joints + champ weakPlane, lus dans le VTU du solveur -----------------
raw = open(f"{outdir}/fdem_joints_0000.vtu", "rb").read().decode("utf8", "ignore")

def darr(name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  raw, re.S)
    return np.fromstring(m.group(1).strip(), sep=" ") if m else None

pts = darr("Points") if darr("Points") is not None else None
m = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', raw, re.S)
pts = np.fromstring(m.group(1).strip(), sep=" ").reshape(-1, 3)[:, :2]
m = re.search(r'<DataArray[^>]*Name="connectivity"[^>]*>(.*?)</DataArray>', raw, re.S)
conn = np.fromstring(m.group(1).strip(), sep=" ", dtype=int).reshape(-1, 2)
wp = darr("weakPlane")
print(f"  joints   : {len(conn)}, dont {int(wp.sum())} sur plans de schistosite")

segs = pts[conn]
weak = segs[wp > 0.5]

# ---- figure ---------------------------------------------------------------
fig = plt.figure(figsize=(10.5, 5.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.16)

for ax, (half, ttl) in zip(
        [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])],
        [(50.0, "(a)  modele complet $100\\times100$ m"),
         (11.0, "(b)  zoom sur la galerie")]):
    x0, x1 = 50 - half, 50 + half
    keep = np.all((xy[tris][:, :, 0] > x0 - 3) & (xy[tris][:, :, 0] < x1 + 3) &
                  (xy[tris][:, :, 1] > 50 - half - 3) &
                  (xy[tris][:, :, 1] < 50 + half + 3), axis=1)
    ax.add_collection(PolyCollection(
        xy[tris[keep]], facecolors="none", edgecolors="0.75",
        linewidths=0.10 if half > 20 else 0.30))
    kw = np.all((weak[:, :, 0] > x0 - 2) & (weak[:, :, 0] < x1 + 2) &
                (weak[:, :, 1] > 50 - half - 2) &
                (weak[:, :, 1] < 50 + half + 2), axis=1)
    ax.add_collection(LineCollection(
        weak[kw], colors="#c0392b", linewidths=0.45 if half > 20 else 1.1))
    ax.set_xlim(x0, x1); ax.set_ylim(50 - half, 50 + half)
    ax.set_aspect("equal"); ax.set_title(ttl, fontsize=10, pad=7)
    ax.set_xlabel("$x$ [m]"); ax.set_ylabel("$y$ [m]")
    ax.tick_params(labelsize=8)

T1 = ("Maillage du tunnel de Hutou Beishan (106 298 triangles) "
      "et schistosite a $45^\\circ$, espacement $1{,}5$ m")
T2 = ("en rouge : les 8782 joints selectionnes par la cle weakPlanes, "
      "lus dans la sortie du solveur")
fig.suptitle(T1 + "\n" + T2, fontsize=10.5, y=0.985)
fig.subplots_adjust(top=0.86, bottom=0.10)
fig.savefig("tunnel_schisto/mesh_schisto45.png", dpi=190)
fig.savefig("tunnel_schisto/mesh_schisto45.pdf")
print("  ecrit : tunnel_schisto/mesh_schisto45.png / .pdf")
