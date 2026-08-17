# -*- coding: utf-8 -*-
"""Qui fait monter la force ? Les eclats qui s'envolent, ou le tas de debris
qui ne peut pas s'evacuer ? Reponse : le tas."""
import csv, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
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


fs = sorted(glob.glob(os.path.join(RUN, "fdem_[0-9][0-9][0-9][0-9].vtu")))
conn = arr(open(fs[0]).read(), "connectivity").astype(int).reshape(-1, 3)
rows = [r for r in csv.DictReader(open(os.path.join(RUN, "history.csv")))
        if all(v not in (None, "") for v in r.values())]
x = np.array([float(r["toolX"]) for r in rows]) * 1e3
fx = np.abs(np.array([float(r["toolFx"]) for r in rows])) * 1e-3

xt, pile, hot, base = [], [], [], []
for f in fs[16:]:
    t = open(f).read()
    P = pts(t) * 1e3
    cy = P[conn][:, :, 1].mean(axis=1)
    vm = np.linalg.norm(arr(t, "velocity", 3)[:, :2], axis=1)
    ve = vm[conn].mean(axis=1)
    i = int(os.path.basename(f)[5:9])
    xx = x[min(len(x) - 1, int(i / (len(fs) - 1) * (len(x) - 1)))]
    s = np.abs(x - xx) < 0.08
    if s.sum() < 4:
        continue
    xt.append(xx); pile.append((cy > 15.02).sum())
    hot.append((ve > 20).sum()); base.append(np.percentile(fx[s], 10))
xt, pile, hot, base = map(np.array, (xt, pile, hot, base))

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))

# (a) les deux candidats au cours du temps
a = ax[0]
a.plot(xt, pile, "C0-o", ms=3.5, lw=1.6, label="tas de débris (élém. au-dessus\nde la surface initiale)")
a.plot(xt, hot, "crimson", marker="s", ms=3.5, lw=1.6, label="éclats > 20 m/s")
a.set_xlabel("position de l'arête (mm)"); a.set_ylabel("nombre d'éléments")
a.set_title("(a) Le tas grossit ; les éclats, non", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# (b) correlation force / tas
a = ax[1]
a.plot(pile, base, "C0-o", ms=5, lw=1.4)
for i in range(0, len(pile), max(1, len(pile) // 5)):
    a.annotate("%.1f mm" % xt[i], (pile[i], base[i]), fontsize=7,
               xytext=(4, -8), textcoords="offset points", color="0.35")
a.set_xlabel("taille du tas (éléments)")
a.set_ylabel("force de base $F_x$ (kN/m, hors pics)")
a.set_title("(b) La force suit le tas, pas les éclats", fontsize=10)
a.grid(alpha=0.3)

# (c) le tas, en image
a = ax[2]
tl = open(fs[-1]).read()
P = pts(tl) * 1e3
cy = P[conn][:, :, 1].mean(axis=1)
inpile = cy > 15.02
cols = ["crimson" if p else "#ddd9d0" for p in inpile]
a.add_collection(PolyCollection([P[c] for c in conn], facecolors=cols,
                                edgecolors="0.55", linewidths=0.15))
a.axhline(15.0, color="C0", lw=1.2, ls="--")
a.text(x[-1] + 2.5, 15.3, "surface initiale", color="C0", fontsize=8)
aa = np.deg2rad(RAKE); tip = np.array([x[-1], H - DOC])
top = tip + 3.5 * np.array([-np.sin(aa), np.cos(aa)])
a.add_patch(Polygon([tip, top, top + [-2, 0], tip + [-2, 0]], closed=True,
                    facecolor="0.3", edgecolor="k", lw=1, zorder=7))
a.set_xlim(x[-1] - 4, x[-1] + 8); a.set_ylim(H - 5, H + 5); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
a.set_title("(c) %d éléments empilés devant la face\nd'attaque, sans issue en 2D"
            % inpile.sum(), fontsize=10)

fig.suptitle("La force ne monte pas à cause des éclats éjectés, mais du tas de "
             "débris qui ne peut pas s'évacuer", fontsize=12.5)
fig.tight_layout()
out = os.path.join(BASE, "figures", "cutting_graded_cause.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
print("tas %d -> %d elements | force base %.0f -> %.0f kN/m | eclats chauds %d -> %d"
      % (pile[0], pile[-1], base[0], base[-1], hot[0], hot[-1]))
