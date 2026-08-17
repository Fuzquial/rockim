# -*- coding: utf-8 -*-
"""Les maillages du balayage de convergence, en coupe, avec l'insert — et les
branches de charge correspondantes. Objectif : voir CE QUI EST RAFFINE et ce
qui ne l'est pas.

  python plot_meshes_contact.py
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
TWIN = 4e-5

MESHES = [("p1_banc_mid", "82k uniforme\nroche 2,76 mm", "C0"),
          ("p1_grad_15",  "gradué\nroche 1,50 mm",       "C1"),
          ("p1_banc_fin", "259k uniforme\nroche 1,88 mm", "C2"),
          ("p1_grad_7",   "gradué\nroche 0,70 mm",       "C3")]
FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def read_msh(path):
    """MSH 2.2 ASCII -> (points [mm], tets, tag physique par tet)."""
    nodes, tets, phys = {}, [], []
    with open(path) as fh:
        it = iter(fh)
        for line in it:
            if line.startswith("$Nodes"):
                for _ in range(int(next(it))):
                    t = next(it).split()
                    nodes[int(t[0])] = (float(t[1]), float(t[2]), float(t[3]))
            elif line.startswith("$Elements"):
                for _ in range(int(next(it))):
                    t = next(it).split()
                    if t[1] == "4":                     # tetraedre
                        ntag = int(t[2])
                        phys.append(int(t[3]) if ntag >= 1 else 0)
                        tets.append([int(x) for x in t[3 + ntag:3 + ntag + 4]])
    P = np.array([nodes[i] for i in range(1, len(nodes) + 1)]) * 1e3
    return P, np.array(tets) - 1, np.array(phys)


def load_hist(sub):
    p = os.path.join(ROOT, "out_elast_" + sub, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if len(rows) < 50 or float(rows[-1]["t"]) < 0.95 * TWIN:
        return None
    g = lambda k: np.array([float(r[k]) for r in rows])
    z = g("grpZ") * 1e3
    return z[0] - z, np.abs(g("grpFz")) * 1e-3


fig = plt.figure(figsize=(16, 8.6))
gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 1])

for j, (name, lab, col) in enumerate(MESHES):
    a = fig.add_subplot(gs[0, j])
    p = os.path.join(ROOT, "meshes", name + ".msh")
    if not os.path.exists(p):
        a.text(.5, .5, "maillage absent", ha="center"); a.axis("off"); continue
    P, T, ph = read_msh(p)
    cen = P[T].mean(axis=1)
    band = np.abs(cen[:, 1] - 60.0) < 1.5
    tags = np.unique(ph)
    ins_tag = tags[np.argmin([(ph == t).sum() for t in tags])]  # l'insert = le petit
    for sel, fc, ec in ((band & (ph != ins_tag), "#e8e4da", "0.55"),
                        (band & (ph == ins_tag), "#9fb6c4", "0.25")):
        polys = []
        for k in np.where(sel)[0]:
            q = P[T[k]][:, [0, 2]]
            polys += [q[list(f)] for f in FACES]
        if polys:
            a.add_collection(PolyCollection(polys, facecolors=fc,
                                            edgecolors=ec, linewidths=0.25))
    a.axhline(120, color="C3", lw=1.0, ls="--")
    n_ins = (ph == ins_tag).sum()
    a.set_xlim(38, 82); a.set_ylim(102, 136); a.set_aspect("equal")
    a.set_xlabel("x (mm)")
    if j == 0:
        a.set_ylabel("z (mm)")
    a.set_title("%s\n%d k tets — insert : %d tets" % (lab, len(T) / 1000, n_ins),
                fontsize=9.5, color=col)

# --- branches de charge ---------------------------------------------------
a = fig.add_subplot(gs[1, :2])
got = []
for name, lab, col in MESHES:
    h = load_hist(name)
    if h is None:
        continue
    pen, fz = h
    m = pen > 0
    a.plot(pen[m], fz[m], color=col, lw=1.8,
           label=lab.replace("\n", " — "))
    got.append((name, lab, col, pen, fz))
for x in (0.10, 0.15, 0.20):
    a.axvline(x, color="0.85", lw=0.8, zorder=0)
a.set_xlabel("pénétration de l'insert (mm)")
a.set_ylabel(r"force de contact $|F_z|$ (kN)")
a.set_title("Branche de charge — élastique pur, aucune rupture", fontsize=10)
a.legend(fontsize=8.5); a.grid(alpha=0.3)

# --- force a penetration fixee -------------------------------------------
a = fig.add_subplot(gs[1, 2:])
for x, mk in ((0.10, "o"), (0.15, "s"), (0.20, "^")):
    hs, fs, cs = [], [], []
    for name, lab, col, pen, fz in got:
        if pen.max() < x:
            continue
        hs.append(float(lab.split()[-2].replace(",", ".")))
        fs.append(float(np.interp(x, pen, fz))); cs.append(col)
    if len(hs) >= 2:
        o = np.argsort(hs)[::-1]
        a.plot(np.array(hs)[o], np.array(fs)[o], mk + "-", lw=1.5, ms=7,
               label=r"$\delta$ = %.2f mm" % x)
a.set_xlabel("taille d'élément dans la ROCHE (mm)  —  l'insert reste à 2,8 mm")
a.set_ylabel(r"$|F_z|$ à pénétration fixée (kN)")
a.set_title("Convergence : la courbe s'aplatit-elle ?", fontsize=10)
a.invert_xaxis(); a.legend(fontsize=8.5); a.grid(alpha=0.3)

fig.suptitle("Ce qui est raffiné, et ce qui ne l'est pas — bleu = insert "
             "(inchangé dans les 4 cas), beige = roche", fontsize=12.5)
fig.tight_layout()
out = os.path.join(HERE, "meshes_contact.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
