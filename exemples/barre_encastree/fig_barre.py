#!/usr/bin/env python3
"""Deformee de la barre encastree (amplifiee) et histoire de l appui.
   python3 exemples/barre_encastree/fig_barre.py <dossier_de_sortie_du_run>"""
import os, re, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9})
run = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))

def vtu(p):
    t = open(p).read(); out = {}
    for m in re.finditer(r'<DataArray([^>]*)>(.*?)</DataArray>', t, re.S):
        nm = re.search(r'Name="([^"]*)"', m.group(1))
        nc = re.search(r'NumberOfComponents="(\d+)"', m.group(1))
        a = np.array(m.group(2).split(), dtype=float)
        if nc and int(nc.group(1)) > 1: a = a.reshape(-1, int(nc.group(1)))
        out.setdefault(nm.group(1) if nm else "Points", a)
    return out
fr = sorted(f for f in os.listdir(run) if re.match(r"fdem3d_\d+\.vtu", f))
v0, v1 = vtu(os.path.join(run, fr[0])), vtu(os.path.join(run, fr[-1]))
X0, X1 = v0["Points"], v1["Points"]
U = X1 - X0
amp = 1000.0
conn = v1["connectivity"].astype(int).reshape(-1, 4)
faces = np.concatenate([conn[:, [1, 2, 3]], conn[:, [0, 3, 2]],
                        conn[:, [0, 1, 3]], conn[:, [0, 2, 1]]])
key = np.sort(X0[faces].round(9).reshape(len(faces), -1), axis=1)
# faces exterieures : celles dont les trois sommets sont sur la boite
lo, hi = X0.min(0), X0.max(0)
on = lambda P: (np.isclose(P, lo, atol=1e-9) | np.isclose(P, hi, atol=1e-9))
ext = np.array([any(on(X0[f]).all(axis=0)) for f in faces])
F = faces[ext]
P = (X0 + amp * U)[F] * 1e3
ux = U[F][:, :, 0].mean(1) * 1e6
fig = plt.figure(figsize=(7.2, 2.9))
ax = fig.add_subplot(1, 2, 1, projection="3d")
norm = plt.Normalize(ux.min(), 0)
pc = Poly3DCollection(P, facecolor=plt.cm.viridis(norm(ux)), edgecolor="k",
                      linewidths=0.08)
ax.add_collection3d(pc)
ax.set_xlim(-12, 100); ax.set_ylim(0, 20); ax.set_zlim(0, 20)
ax.set_box_aspect((5.6, 1, 1)); ax.view_init(20, -60)
ax.set_xlabel("$x$ [mm]"); ax.tick_params(labelsize=7)
sm = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
fig.colorbar(sm, ax=ax, orientation="horizontal", shrink=0.7, pad=0.02,
             label=r"$u_x$ [$\mu$m]")
ax.set_yticks([]); ax.set_zticks([])
ax.set_title("deformee x 1000 (encastree en x = 0, pression 5 MPa en x = 100 mm)",
             fontsize=8)
h = list(csv.DictReader(open(os.path.join(run, "history.csv"))))
t = np.array([float(r["t"]) for r in h]) * 1e6
a2 = fig.add_subplot(1, 2, 2)
a2.plot(t, [-float(r["F_xmax_x"]) for r in h], "r-", lw=1, label="charge (xmax)")
a2.plot(t, [float(r["RF_xmin_x"]) for r in h], "k-", lw=0.8,
        label="reaction de l encastrement (xmin)")
a2.set_xlabel(r"$t$ [$\mu$s]"); a2.set_ylabel("force selon $x$ [N]")
a2.legend(frameon=False, loc="lower right")
a2.set_title("fin : 2 000,1 N appliques, 2 000,7 N repris", fontsize=8)
fig.tight_layout()
for e in ("pdf", "png"):
    fig.savefig(os.path.join(here, "barre_resultat." + e), dpi=200,
                bbox_inches="tight")
