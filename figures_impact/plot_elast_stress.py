# -*- coding: utf-8 -*-
"""Balayage elastique : champs de von Mises A PENETRATION EGALE (seule
comparaison licite, l'instant du meme etat differant d'un maillage a l'autre)
+ branches de charge.

Tailles d'element MESUREES dans le maillage, pas les valeurs nominales du
generateur (elles differaient d'un facteur 1,7 — corrige le 2026-08-17).

  python plot_elast_stress.py [penetration_cible_mm]
"""
import csv, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
TWIN, PTARGET = 4e-5, float(sys.argv[1]) if len(sys.argv) > 1 else 0.20
FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))

# (run, libelle, h roche MESURE, h insert MESURE, couleur)
RUNS = [("p1_banc_mid", "82k uniforme",  4.74, 2.82, "C0"),
        ("p1_banc_fin", "259k uniforme", 3.29, 1.88, "C2"),
        ("p1_grad_15",  "gradué 1,5 mm", 1.69, 1.66, "C1"),
        ("p1_grad_7",   "gradué 0,7 mm", 0.80, 0.80, "C3")]


def txt(p):
    return open(p, encoding="utf-8", errors="ignore").read()


def arr(s, n, d=float):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, s, re.S)
    return None if not m else np.fromstring(m.group(1), sep=" ", dtype=d)


def pts_of(s):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", s, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)


def load(run):
    D = os.path.join(ROOT, "out_elast_" + run)
    p = os.path.join(D, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if len(rows) < 50 or float(rows[-1]["t"]) < 0.95 * TWIN:
        return None
    g = lambda k: np.array([float(r[k]) for r in rows])
    z = g("grpZ") * 1e3
    pen, fz, t = z[0] - z, np.abs(g("grpFz")) * 1e-3, g("t")
    # frame dont la PENETRATION est la plus proche de la cible
    fr = {int(r["frame"]): float(r["t"])
          for r in csv.DictReader(open(os.path.join(D, "frames.csv")))}
    best, bd = None, 1e9
    for k, tk in fr.items():
        pk = float(np.interp(tk, t, pen))
        if abs(pk - PTARGET) < bd:
            best, bd, bp = k, abs(pk - PTARGET), pk
    return dict(pen=pen, fz=fz, D=D, frame=best, penf=bp,
                nb=int(rows[-1]["nBroken"]))


data = []
for run, lab, hr, hi, c in RUNS:
    d = load(run)
    if d is None:
        print("  (en cours ou absent : %s)" % run)
        continue
    d.update(run=run, lab=lab, hr=hr, hi=hi, c=c)
    data.append(d)
if not data:
    raise SystemExit("rien a tracer")

n = len(data)
fig = plt.figure(figsize=(4.0 * n, 8.6))
gs = fig.add_gridspec(2, n, height_ratios=[1.2, 1])

for j, d in enumerate(data):
    a = fig.add_subplot(gs[0, j])
    f = os.path.join(d["D"], "fdem3d_%04d.vtu" % d["frame"])
    s = txt(f)
    P = pts_of(s) * 1e3
    T = arr(s, "connectivity", int).reshape(-1, 4)
    vm = arr(s, "vonMises") * 1e-6
    ph = arr(s, "phase", int)
    cen = P[T].mean(axis=1)
    band = (ph == 0) & (np.abs(cen[:, 1] - 60) < 1.5)
    polys, vals = [], []
    for k in np.where(band)[0]:
        q = P[T[k]][:, [0, 2]]
        for fa in FACES:
            polys.append(q[list(fa)]); vals.append(vm[k])
    pc = a.add_collection(PolyCollection(polys, array=np.array(vals),
                                         cmap="inferno", edgecolors="none"))
    pc.set_clim(0, 300)
    if j == n - 1:
        cb = fig.colorbar(pc, ax=a, shrink=0.85); cb.set_label("von Mises (MPa)")
    a.axhline(120, color="C0", lw=0.9, ls="--")
    a.set_xlim(40, 80); a.set_ylim(96, 124); a.set_aspect("equal")
    a.set_xlabel("x (mm)")
    if j == 0:
        a.set_ylabel("z (mm)")
    a.set_title("%s\nroche %.2f mm — δ = %.3f mm"
                % (d["lab"], d["hr"], d["penf"]), fontsize=9.5, color=d["c"])

a = fig.add_subplot(gs[1, :max(1, n // 2)])
for d in data:
    m = d["pen"] > 0
    a.plot(d["pen"][m], d["fz"][m], color=d["c"], lw=1.8,
           label="%s — roche %.2f mm" % (d["lab"], d["hr"]))
a.axvline(PTARGET, color="0.6", lw=1.0, ls=":")
a.set_xlabel("pénétration (mm)"); a.set_ylabel(r"$|F_z|$ (kN)")
a.set_title("Branches de charge (élastique pur, 0 rupture)", fontsize=10)
a.legend(fontsize=8.5); a.grid(alpha=0.3)

a = fig.add_subplot(gs[1, max(1, n // 2):])
for x, mk in ((0.10, "o"), (0.15, "s"), (0.20, "^")):
    hs, fs = [], []
    for d in sorted(data, key=lambda q: -q["hr"]):
        if d["pen"].max() < x:
            continue
        hs.append(d["hr"]); fs.append(float(np.interp(x, d["pen"], d["fz"])))
    if len(hs) >= 2:
        a.plot(hs, fs, mk + "-", lw=1.6, ms=7, label=r"$\delta$ = %.2f mm" % x)
a.set_xlabel("taille d'élément MESURÉE dans la roche (mm)")
a.set_ylabel(r"$|F_z|$ à pénétration fixée (kN)")
a.set_title("Convergence : y a-t-il une asymptote ?", fontsize=10)
a.invert_xaxis(); a.legend(fontsize=8.5); a.grid(alpha=0.3)

fig.suptitle("Convergence du contact — champs comparés à pénétration égale "
             "(δ ≈ %.2f mm)" % PTARGET, fontsize=12.5)
fig.tight_layout()
out = os.path.join(HERE, "elast_stress.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
for d in sorted(data, key=lambda q: -q["hr"]):
    print("  %-15s roche %.2f mm  frame %d a delta %.3f mm  F = %.2f kN  "
          "joints rompus %d" % (d["lab"], d["hr"], d["frame"], d["penf"],
                                float(np.interp(PTARGET, d["pen"], d["fz"])),
                                d["nb"]))
