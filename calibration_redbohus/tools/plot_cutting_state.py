# -*- coding: utf-8 -*-
"""Etat courant de la simulation de coupe : ou casse-t-il, et la force
est-elle localisee devant l'outil ou diffuse dans tout le bloc ?"""
import csv, glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUN = os.path.join(BASE, "runs", sys.argv[1] if len(sys.argv) > 1
                   else "demo_cutting2")


def arr(t, n):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    return np.fromstring(m.group(1), sep=" ") if m else None


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


fs = sorted(glob.glob(os.path.join(RUN, "fdem_[0-9][0-9][0-9][0-9].vtu")))
t0 = open(fs[0]).read()
conn = arr(t0, "connectivity").astype(int).reshape(-1, 3)
rows = [r for r in csv.DictReader(open(os.path.join(RUN, "history.csv")))
        if all(v not in (None, "") for v in r.values())]
g = lambda k: np.array([float(r[k]) for r in rows])
x, fx, t = g("toolX") * 1e3, np.abs(g("toolFx")) * 1e-3, g("t") * 1e3

fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.0),
                       gridspec_kw={"width_ratios": [1, 1, 1.15]})

for j, idx in enumerate((max(0, len(fs) // 2), len(fs) - 1)):
    a = ax[j]
    P = pts(open(fs[idx]).read()) * 1e3
    a.add_collection(PolyCollection([P[c] for c in conn], facecolors="0.88",
                                    edgecolors="0.6", linewidths=0.1))
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
            a.add_collection(LineCollection(segs, colors="crimson", lw=1.2))
    # position de l'outil a cette frame
    tf = t[-1] * idx / max(1, len(fs) - 1)
    xt = np.interp(tf, t, x)
    a.plot([xt, xt], [0, 18], "C0--", lw=1.2)
    a.annotate("outil", (xt, 19), color="C0", fontsize=9, ha="center")
    P0 = pts(t0) * 1e3
    a.set_xlim(-3, P0[:, 0].max() + 2); a.set_ylim(-1, 22)
    a.set_aspect("equal"); a.set_xlabel("x (mm)")
    a.set_title("frame %d — %d joints rompus" % (idx, n), fontsize=10)

a = ax[2]
a.plot(x, fx, "C0", lw=1.0)
a.axhline(300, color="C2", lw=1.6, ls="--")
a.text(x.min() + 0.3, 380, "ordre attendu ≈ 300 kN/m\n(MSE granite × 1,5 mm)",
       color="C2", fontsize=9)
a.set_yscale("log")
a.set_xlabel("position de l'outil (mm)")
a.set_ylabel("force de coupe $F_x$ (kN/m)")
a.set_title("Force de coupe — échelle log")
a.grid(alpha=0.3, which="both")

fig.suptitle("Coupe PDC à 2 m/s : la force reste 10× trop élevée — "
             "diagnostic en cours", fontsize=12)
fig.tight_layout()
out = os.path.join(BASE, "figures", "cutting_state.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
print("force mediane (en matiere) %.0f kN/m | attendu ~300"
      % np.median(fx[x > 0]))
