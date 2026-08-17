# -*- coding: utf-8 -*-
"""Diagnostic du run de coupe gradue : il DIVERGE. Trois preuves cote a cote —
force qui croit sans borne, energie creee, noeuds a vitesse aberrante."""
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
H, DOC, RAKE, VCUT = 15.0, 0.5, 20.0, 2.0


def arr(t, n, nc=1):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    if not m:
        return None
    v = np.fromstring(m.group(1), sep=" ")
    return v.reshape(-1, nc) if nc > 1 else v


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


rows = [r for r in csv.DictReader(open(os.path.join(RUN, "history.csv")))
        if all(v not in (None, "") for v in r.values())]
g = lambda k: np.array([float(r[k]) for r in rows])
x, fx = g("toolX") * 1e3, np.abs(g("toolFx")) * 1e-3
W = g("work")
diss = np.abs(g("eCund")) + np.abs(g("eLys")) + np.abs(g("eGc")) + np.abs(g("eFric"))

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.7))

# (a) force : croissance sans borne
a = ax[0]
m = x > -0.2
a.plot(x[m], fx[m], "C0", lw=0.9)
for lo, hi in [(0, .5), (.5, 1.), (1., 1.5), (1.5, 2.2)]:
    s = (x >= lo) & (x < hi)
    if s.sum():
        a.hlines(np.median(fx[s]), lo, hi, color="crimson", lw=2.4)
        a.text((lo + hi) / 2, np.median(fx[s]) * 1.5, "%.0f" % np.median(fx[s]),
               color="crimson", ha="center", fontsize=8)
a.axhline(68, color="C2", ls="--", lw=1.4, label="v3 (GBM) : 68 kN/m")
a.set_yscale("log")
a.set_xlabel("position de l'arête (mm)"); a.set_ylabel("$F_x$ (kN/m)")
a.set_title("(a) La force croît sans borne\n(médianes rouges par tranche)",
            fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, which="both")

# (b) energie : creee, pas dissipee
a = ax[1]
a.plot(x[m], W[m] * 1e-3, "C0", lw=1.6, label="travail de l'outil")
a.plot(x[m], diss[m] * 1e-3, "crimson", lw=1.6,
       label="dissipé (amort. + Lysmer + $G_c$ + frott.)")
a.set_xlabel("position de l'arête (mm)"); a.set_ylabel("énergie (kJ/m)")
a.set_title("(b) 30 × plus dissipé qu'injecté\n→ de l'énergie est CRÉÉE",
            fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# (c) champ de vitesse final
a = ax[2]
tl = open(sorted(glob.glob(os.path.join(RUN, "fdem_[0-9][0-9][0-9][0-9].vtu")))[-1]).read()
conn = arr(tl, "connectivity").astype(int).reshape(-1, 3)
P = pts(tl) * 1e3
vm = np.linalg.norm(arr(tl, "velocity", 3)[:, :2], axis=1)
vel_el = vm[conn].mean(axis=1)
pc = a.add_collection(PolyCollection([P[c] for c in conn], array=np.log10(np.maximum(vel_el, 1e-3)),
                                     cmap="inferno", edgecolors="none"))
pc.set_clim(-2, 3)
cb = fig.colorbar(pc, ax=a, shrink=0.85)
cb.set_label(r"$\log_{10}|v|$  (m/s)")
xt = x[-1]
aa = np.deg2rad(RAKE); tip = np.array([xt, H - DOC])
top = tip + 3.0 * np.array([-np.sin(aa), np.cos(aa)])
a.add_patch(Polygon([tip, top, top + [-2, 0], tip + [-2, 0]], closed=True,
                    facecolor="0.3", edgecolor="k", lw=1, zorder=7))
a.set_xlim(xt - 4, xt + 8); a.set_ylim(H - 7, H + 4); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
a.set_title("(c) $v_{max}$ = %.0f m/s pour une coupe à 2 m/s\n(%d nœuds > 20 m/s ; pic de 1725 m/s frame 28)"
            % (vm.max(), int((vm > 20).sum())), fontsize=10)

fig.suptitle("Coupe sur maillage gradué — le calcul DIVERGE (instabilité de "
             "contact), run arrêté", fontsize=12.5)
fig.tight_layout()
out = os.path.join(BASE, "figures", "cutting_graded_diagnostic.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
