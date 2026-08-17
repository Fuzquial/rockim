# -*- coding: utf-8 -*-
"""Coupe PDC v3 : zoom sur la zone de coupe + force. Le domaine fait 80 fois
la passe : on ne montre que les 15 premiers mm, la ou tout se joue."""
import csv, glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUN = os.path.join(BASE, "runs", "demo_cutting3")
H, DOC, RAKE = 20.0, 0.5, 20.0


def arr(t, n):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    return np.fromstring(m.group(1), sep=" ") if m else None


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
rows = [r for r in csv.DictReader(open(os.path.join(RUN, "history.csv")))
        if all(v not in (None, "") for v in r.values())]
g = lambda k: np.array([float(r[k]) for r in rows])
x, fx, t = g("toolX") * 1e3, np.abs(g("toolFx")) * 1e-3, g("t") * 1e3
nb = g("nBroken")

sel = [max(1, len(fs) - 3), max(1, len(fs) - 2), len(fs) - 1]
fig = plt.figure(figsize=(15, 7.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1])

for j, idx in enumerate(sel):
    a = fig.add_subplot(gs[0, j])
    P = pts(open(fs[idx]).read()) * 1e3
    a.add_collection(PolyCollection([P[c] for c in conn], facecolors="#ddd9d0",
                                    edgecolors="0.55", linewidths=0.12))
    jf = fs[idx].replace("fdem_", "fdem_joints_")
    n = 0
    if os.path.exists(jf):
        jt = open(jf).read()
        jP, jc = pts(jt) * 1e3, arr(jt, "connectivity").astype(int)
        tb, off = arr(jt, "tBreak"), arr(jt, "offsets").astype(int)
        segs, st = [], 0
        for i, en in enumerate(off):
            ii = jc[st:en]; st = en
            if tb is not None and tb[i] >= 0 and len(ii) >= 2:
                segs.append(jP[ii[:2]])
        n = len(segs)
        if segs:
            a.add_collection(LineCollection(segs, colors="crimson", lw=1.6))
    tf = t[-1] * idx / max(1, len(fs) - 1)
    xt = float(np.interp(tf, t, x))
    cutter(a, xt)
    a.set_xlim(xt - 4, xt + 9); a.set_ylim(H - 6, H + 4)
    a.set_aspect("equal"); a.set_xlabel("x (mm)")
    if j == 0:
        a.set_ylabel("y (mm)")
    a.set_title("outil à %.2f mm — %d joints rompus" % (xt, n), fontsize=10)

# --- force ---------------------------------------------------------------
a = fig.add_subplot(gs[1, :])
m = x > -0.1
a.plot(x[m], fx[m], "C0", lw=1.1)
a.axhline(np.median(fx[x > 0]) if (x > 0).any() else 0, color="C3", lw=1.3,
          ls="--", label="médiane %.0f kN/m" % np.median(fx[x > 0]))
a2 = a.twinx()
a2.plot(x[m], nb[m], "C7", lw=1.2, alpha=0.7)
a2.set_ylabel("joints rompus", color="C7")
a.set_xlabel("position de l'arête (mm)")
a.set_ylabel("force de coupe $F_x$ (kN/m)")
a.set_title("Force de coupe — la signature cyclique de la coupe fragile",
            fontsize=11)
a.legend(fontsize=9); a.grid(alpha=0.3)

fig.suptitle("Coupe PDC — domaine 40 × 20 mm, passe 0,5 mm (80 × la passe en "
             "longueur), 2893 éléments", fontsize=12)
fig.tight_layout()
out = os.path.join(BASE, "figures", "cutting_v3.png")
fig.savefig(out, dpi=140)
print("ecrit", out, "| frames", len(fs))
print("force en matiere : mediane %.0f, max %.0f kN/m | joints %d"
      % (np.median(fx[x > 0]), fx[x > 0].max(), nb[-1]))
