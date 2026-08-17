# -*- coding: utf-8 -*-
"""Les fragments sont-ils EJECTES ? Champ de vitesse dans la zone de coupe,
fragments detaches identifies, et vitesse d'ejection comparee a la vitesse
de coupe."""
import csv, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUN = os.path.join(BASE, "runs", "demo_cutting3")
H, DOC, RAKE, VCUT = 20.0, 0.5, 20.0, 2.0


def arr(t, n, ncomp=1):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    if not m:
        return None
    v = np.fromstring(m.group(1), sep=" ")
    return v.reshape(-1, ncomp) if ncomp > 1 else v


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def cutter(ax, xt, shown=4.0):
    a = np.deg2rad(RAKE)
    tip = np.array([xt, H - DOC])
    top = tip + shown * np.array([-np.sin(a), np.cos(a)])
    ax.add_patch(Polygon([tip, top, top + [-2.5, 0], tip + [-2.5, 0]],
                         closed=True, facecolor="0.3", edgecolor="k",
                         lw=1.0, zorder=6))


fs = sorted(glob.glob(os.path.join(RUN, "fdem_[0-9][0-9][0-9][0-9].vtu")))
t0 = open(fs[0]).read()
conn = arr(t0, "connectivity").astype(int).reshape(-1, 3)
tl = open(fs[-1]).read()
P = pts(tl) * 1e3
frag = arr(tl, "fragment")
vel = arr(tl, "velocity", 3)
rows = [r for r in csv.DictReader(open(os.path.join(RUN, "history.csv")))
        if all(v not in (None, "") for v in r.values())]
x = np.array([float(r["toolX"]) for r in rows]) * 1e3
xt = x[-1]

# vitesse par element (moyenne des 3 noeuds)
vmag = np.linalg.norm(vel[:, :2], axis=1)
vel_el = vmag[conn].mean(axis=1)
# fragment principal = le plus nombreux
ids, cnt = np.unique(frag, return_counts=True)
main = ids[np.argmax(cnt)]
det = frag != main

fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.0),
                       gridspec_kw={"width_ratios": [1, 1, 1.1]})

# (a) champ de vitesse
a = ax[0]
pc = a.add_collection(PolyCollection([P[c] for c in conn], array=vel_el,
                                     cmap="inferno", edgecolors="none"))
pc.set_clim(0, max(VCUT, vel_el.max()))
cb = fig.colorbar(pc, ax=a, shrink=0.8); cb.set_label("|v| (m/s)")
cutter(a, xt)
a.set_xlim(xt - 4, xt + 9); a.set_ylim(H - 6, H + 4); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
a.set_title("Champ de vitesse (outil à %.0f m/s)" % VCUT, fontsize=10)

# (b) fragments detaches
a = ax[1]
cols = ["crimson" if d else "0.85" for d in det]
a.add_collection(PolyCollection([P[c] for c in conn], facecolors=cols,
                                edgecolors="0.6", linewidths=0.12))
cutter(a, xt)
a.set_xlim(xt - 4, xt + 9); a.set_ylim(H - 6, H + 4); a.set_aspect("equal")
a.set_xlabel("x (mm)")
a.set_title("Fragments détachés (rouge) — %d élément(s) sur %d"
            % (det.sum(), len(conn)), fontsize=10)

# (c) histogramme des vitesses
a = ax[2]
if det.any():
    a.hist(vel_el[det], bins=20, color="crimson", alpha=0.85,
           label="éléments détachés")
a.hist(vel_el[~det], bins=40, color="0.6", alpha=0.6, label="massif")
a.axvline(VCUT, color="C0", lw=2, ls="--", label="vitesse de coupe 2 m/s")
a.set_yscale("log")
a.set_xlabel("|v| (m/s)"); a.set_ylabel("nombre d'éléments")
a.set_title("Vitesses : les débris sont-ils éjectés ?", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

fig.suptitle("Éjection des débris — coupe PDC à %.2f mm de course" % xt,
             fontsize=12)
fig.tight_layout()
out = os.path.join(BASE, "figures", "cutting_ejecta.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
print("elements detaches : %d | v max detaches : %.2f m/s | v max massif : %.2f m/s"
      % (det.sum(), vel_el[det].max() if det.any() else 0, vel_el[~det].max()))
print("v mediane des detaches : %.2f m/s (vitesse de coupe %.1f)"
      % (np.median(vel_el[det]) if det.any() else 0, VCUT))
