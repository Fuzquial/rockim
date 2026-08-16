# -*- coding: utf-8 -*-
"""Demo : rupture nette + resistance residuelle, FDEM 2D homogene."""
import csv, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUN = os.path.join(BASE, "runs", "demo_rupture")


def arr(t, n):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    return np.fromstring(m.group(1), sep=" ") if m else None


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


rows = [r for r in csv.DictReader(open(os.path.join(RUN, "history.csv")))
        if all(v not in (None, "") for v in r.values())]
h = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
s = np.abs(h["sigma"]) * 1e-6
e = np.abs(h["epsGauge"]) * 100
nb = h["nBroken"]
ip = int(np.argmax(s))

fig = plt.figure(figsize=(14.5, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1, 1])

# --- courbe -----------------------------------------------------------
a = fig.add_subplot(gs[0])
a.plot(e, s, "C0", lw=2)
a.plot(e[ip], s[ip], "C3o", ms=8)
res = np.median(s[e > e[ip] * 2.5]) if (e > e[ip] * 2.5).any() else s[-1]
a.axhline(res, color="C2", lw=1.4, ls="--")
a.annotate("PIC %.1f MPa\nà ε = %.3f %%" % (s[ip], e[ip]), (e[ip], s[ip]),
           textcoords="offset points", xytext=(14, -6), fontsize=10, color="C3")
a.annotate("résistance résiduelle ≈ %.1f MPa\n(%.0f %% du pic)"
           % (res, 100 * res / s[ip]), (e.max() * 0.55, res),
           textcoords="offset points", xytext=(0, 14), fontsize=10, color="C2")
a2 = a.twinx()
a2.plot(e, nb, "C7", lw=1.1, alpha=0.7)
a2.set_ylabel("joints rompus", color="C7")
a.set_xlabel("déformation axiale (%)"); a.set_ylabel(r"$\sigma_1$ (MPa)")
a.set_title("Compression simple — rupture nette puis plateau résiduel")
a.grid(alpha=0.3)

# --- facies : deux instants -------------------------------------------
fs = sorted(glob.glob(os.path.join(RUN, "fdem_[0-9][0-9][0-9][0-9].vtu")))
t0 = open(fs[0]).read()
conn = arr(t0, "connectivity").astype(int).reshape(-1, 3)
stat = arr(t0, "ftScale")
for j, (idx, ttl) in enumerate(((min(4, len(fs) - 1), "juste après le pic"),
                                (len(fs) - 1, "état résiduel"))):
    ax = fig.add_subplot(gs[j + 1])
    P = pts(open(fs[idx]).read())
    ax.add_collection(PolyCollection([P[c] * 1e3 for c in conn],
                                     facecolors="0.9", edgecolors="0.65",
                                     linewidths=0.1))
    jf = fs[idx].replace("fdem_", "fdem_joints_")
    n = 0
    if os.path.exists(jf):
        jt = open(jf).read()
        jP, jc = pts(jt), arr(jt, "connectivity").astype(int)
        tb, off = arr(jt, "tBreak"), arr(jt, "offsets").astype(int)
        segs, st = [], 0
        for i, en in enumerate(off):
            ii = jc[st:en]; st = en
            if tb is not None and tb[i] >= 0 and len(ii) >= 2:
                segs.append(jP[ii[:2]] * 1e3)
        n = len(segs)
        if segs:
            ax.add_collection(LineCollection(segs, colors="crimson", lw=1.3))
    P0 = pts(t0)
    ax.set_xlim(P0[:, 0].min() * 1e3 - 2, P0[:, 0].max() * 1e3 + 2)
    ax.set_ylim(P0[:, 1].min() * 1e3 - 2, P0[:, 1].max() * 1e3 + 2)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("%s — %d joints rompus" % (ttl, n), fontsize=10)

fig.suptitle("FDEM 2D homogène, panoplie complète + hétérogénéité de Weibull "
             "(m = 6) : la rupture LOCALISE et la charge s'effondre", fontsize=12)
fig.tight_layout()
out = os.path.join(BASE, "figures", "demo_rupture.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
print("pic %.1f MPa a eps %.3f %% | residuel %.1f MPa (%.0f %%) | joints %d"
      % (s[ip], e[ip], res, 100 * res / s[ip], nb[-1]))
