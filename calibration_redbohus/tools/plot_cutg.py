# -*- coding: utf-8 -*-
"""Coupe PDC sur maillage GRADUE homogene : vue d'ensemble du maillage,
zoom sur la zone de coupe a trois instants, et force."""
import csv, glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUN = os.path.join(BASE, "runs", "demo_cut_graded")
H, DOC, RAKE = 15.0, 0.5, 20.0


def arr(t, n, nc=1):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    if not m:
        return None
    v = np.fromstring(m.group(1), sep=" ")
    return v.reshape(-1, nc) if nc > 1 else v


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def cutter(ax, xt, shown=3.0):
    a = np.deg2rad(RAKE)
    tip = np.array([xt, H - DOC])
    top = tip + shown * np.array([-np.sin(a), np.cos(a)])
    ax.add_patch(Polygon([tip, top, top + [-2.0, 0], tip + [-2.0, 0]],
                         closed=True, facecolor="0.3", edgecolor="k",
                         lw=1.0, zorder=7))


fs = sorted(glob.glob(os.path.join(RUN, "fdem_[0-9][0-9][0-9][0-9].vtu")))
t0 = open(fs[0]).read()
conn = arr(t0, "connectivity").astype(int).reshape(-1, 3)
P0 = pts(t0) * 1e3
rows = [r for r in csv.DictReader(open(os.path.join(RUN, "history.csv")))
        if all(v not in (None, "") for v in r.values())]
g = lambda k: np.array([float(r[k]) for r in rows])
x, fx, t = g("toolX") * 1e3, np.abs(g("toolFx")) * 1e-3, g("t") * 1e3

fig = plt.figure(figsize=(15.5, 8.6))
gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.25, 0.9])

# --- (a) le maillage gradue ------------------------------------------------
a = fig.add_subplot(gs[0, :])
a.add_collection(PolyCollection([P0[c] for c in conn], facecolors="#e6e2d8",
                                edgecolors="0.45", linewidths=0.18))
a.axhline(H - 2.5, color="C0", lw=1.2, ls="--")
a.text(0.5, H - 2.2, "bande raffinée (2,5 mm) — élément 0,20 mm",
       color="C0", fontsize=9)
a.text(0.5, 1.0, "fond grossier — élément 0,44 mm", color="0.35", fontsize=9)
cutter(a, -2.0)
a.set_xlim(-4, 31); a.set_ylim(-1, H + 3); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
a.set_title("Maillage gradué homogène — 2546 triangles, 57 %% dans la bande "
            "de coupe (pas de grains)", fontsize=11)

# --- (b) zoom a trois instants ---------------------------------------------
sel = [len(fs) // 3, 2 * len(fs) // 3, len(fs) - 1]
for j, idx in enumerate(sel):
    a = fig.add_subplot(gs[1, j])
    P = pts(open(fs[idx]).read()) * 1e3
    a.add_collection(PolyCollection([P[c] for c in conn], facecolors="#e6e2d8",
                                    edgecolors="0.55", linewidths=0.15))
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
            a.add_collection(LineCollection(segs, colors="crimson", lw=1.3))
    tf = t[-1] * idx / max(1, len(fs) - 1)
    xt = float(np.interp(tf, t, x))
    cutter(a, xt)
    a.set_xlim(xt - 2.5, xt + 5); a.set_ylim(H - 3.5, H + 2.5)
    a.set_aspect("equal"); a.set_xlabel("x (mm)")
    if j == 0:
        a.set_ylabel("y (mm)")
    a.set_title("outil à %.2f mm — %d joints rompus" % (xt, n), fontsize=10)

# --- (c) force -------------------------------------------------------------
a = fig.add_subplot(gs[2, :])
m = x > -0.2
a.plot(x[m], fx[m], "C0", lw=1.0)
med = np.median(fx[x > 0])
a.axhline(med, color="C3", ls="--", lw=1.3, label="médiane %.0f kN/m" % med)
a.set_xlabel("position de l'arête (mm)")
a.set_ylabel("$F_x$ (kN/m)")
a.set_title("Force de coupe", fontsize=10)
a.legend(fontsize=9); a.grid(alpha=0.3)

fig.suptitle("Coupe PDC — maillage gradué, matériau homogène + Weibull",
             fontsize=13)
fig.tight_layout()
out = os.path.join(BASE, "figures", "cutting_graded.png")
fig.savefig(out, dpi=140)
print("ecrit", out, "| frames", len(fs))
