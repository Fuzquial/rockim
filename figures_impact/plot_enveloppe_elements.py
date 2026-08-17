# -*- coding: utf-8 -*-
"""ENVELOPPE PAR ELEMENTS : on ne trace pas les facettes de joint (qui
s'interpenetrent et que matplotlib trie mal en 3D) mais les TETRAEDRES
adjacents aux joints casses. On obtient un volume solide, une frontiere nette,
et un volume endommage mesurable.

Un tetraedre est classe :
  mode I  seul  -> BLEU     (tous ses joints casses sont en ouverture)
  mode II seul  -> ROUGE    (tous en glissement)
  mixte         -> VIOLET   (les deux)

Une figure par mode, chacune avec la surface du bloc et le cercle de Hertz.

  python plot_enveloppe_elements.py [run]      (defaut out_imp3d_ultra)
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
RUN = sys.argv[1] if len(sys.argv) > 1 else "out_imp3d_ultra"
A_HERTZ = 1.45
BLEU, ROUGE, VIOLET = "#1f6fb4", "#c0392b", "#7d3c98"
FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))

D = os.path.join(ROOT, RUN)
fj = sorted(x for x in os.listdir(D)
            if re.fullmatch(r"fdem3d_joints_\d{4}\.vtu", x))[-1]
fe = fj.replace("_joints", "")


def rd(p, n, d=float):
    s = open(os.path.join(D, p), encoding="utf-8", errors="ignore").read()
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, s, re.S)
    return None if not m else np.fromstring(m.group(1), sep=" ", dtype=d)


def pts_of(p):
    s = open(os.path.join(D, p), encoding="utf-8", errors="ignore").read()
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", s, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3) * 1e3


# --- elements ------------------------------------------------------------
P = pts_of(fe)
T = rd(fe, "connectivity", int).reshape(-1, 4)
Q = P[T]
cen = Q.mean(axis=1)
vol = np.abs(np.einsum("ij,ij->i", Q[:, 1] - Q[:, 0],
                       np.cross(Q[:, 2] - Q[:, 0], Q[:, 3] - Q[:, 0]))) / 6

# --- joints casses -> tetraedres adjacents -------------------------------
# La correspondance joint -> elements n'est PAS dans le VTU : on la retrouve
# geometriquement, chaque facette de joint etant portee par deux tetraedres.
# On associe chaque facette aux tets dont le centroide est le plus proche.
Pj = pts_of(fj)
conn, off = rd(fj, "connectivity", int), rd(fj, "offsets", int)
dmg, bm = rd(fj, "damage"), rd(fj, "breakMode")
idx = np.where(dmg >= 1.0)[0]
jcen = np.array([Pj[conn[(0 if i == 0 else off[i - 1]):off[i]][:3]].mean(axis=0)
                 for i in idx])
jmode = bm[idx]

# rayon de recherche : la taille locale d'element
h = (6 * vol) ** (1 / 3)
hloc = np.median(h[np.hypot(cen[:, 0] - 60, cen[:, 1] - 60) < 5])
flag = np.zeros(len(T), dtype=int)          # bit 1 = mode I, bit 2 = mode II
for c, m in zip(jcen, jmode):
    d2 = np.sum((cen - c) ** 2, axis=1)
    near = np.where(d2 < (0.9 * hloc) ** 2)[0]
    if len(near) == 0:
        near = [int(np.argmin(d2))]
    for k in near:
        flag[k] |= (1 if m == 1 else 2)

onlyI, onlyII, mixte = flag == 1, flag == 2, flag == 3
print("  tetras adjacents a un joint casse : %d" % (flag > 0).sum())
print("    mode I seul  : %5d   volume %7.2f mm3" % (onlyI.sum(), vol[onlyI].sum()))
print("    mode II seul : %5d   volume %7.2f mm3" % (onlyII.sum(), vol[onlyII].sum()))
print("    mixte        : %5d   volume %7.2f mm3" % (mixte.sum(), vol[mixte].sum()))
print("    VOLUME ENDOMMAGE TOTAL : %.2f mm3" % vol[flag > 0].sum())

sel_all = flag > 0
rad = max(3.0, 1.25 * np.hypot(cen[sel_all, 0] - 60, cen[sel_all, 1] - 60).max())
dep = max(3.0, 1.25 * (120 - cen[sel_all, 2]).max())


def tri_proj(mask, cols):
    out = []
    for k in np.where(mask)[0]:
        q = Q[k][:, cols]
        out += [q[list(f)] for f in FACES]
    return out


for tag, principal, col, titre in (
        ("modeI", onlyI | mixte, BLEU, "MODE I — ouverture / traction"),
        ("modeII", onlyII | mixte, ROUGE, "MODE II — glissement / cisaillement")):
    autre = sel_all & ~principal
    fig, ax = plt.subplots(1, 2, figsize=(13.2, 5.6))
    for k, (cols, xl, yl, xlim, ylim, t) in enumerate((
            ([0, 1], "x (mm)", "y (mm)", (60 - rad, 60 + rad),
             (60 - rad, 60 + rad), "(a) Vue de dessus"),
            ([0, 2], "x (mm)", "z (mm)", (60 - rad, 60 + rad),
             (120 - dep, 121), "(b) Coupe — enveloppe complète"))):
        a = ax[k]
        p_autre = tri_proj(autre, cols)
        if p_autre:
            a.add_collection(PolyCollection(p_autre, facecolors="0.85",
                                            edgecolors="0.7", linewidths=0.15))
        pm = tri_proj(mixte & principal, cols)
        if pm:
            a.add_collection(PolyCollection(pm, facecolors=VIOLET,
                                            edgecolors="0.35", linewidths=0.15))
        pp = tri_proj(principal & ~mixte, cols)
        if pp:
            a.add_collection(PolyCollection(pp, facecolors=col,
                                            edgecolors="0.35", linewidths=0.15))
        if k == 0:
            a.add_patch(plt.Circle((60, 60), A_HERTZ, fill=False, ec="k",
                                   lw=1.7, ls="--"))
            a.text(60, 60 + A_HERTZ + 0.1 * rad,
                   "contact Hertz  a = %.2f mm" % A_HERTZ, fontsize=8.5,
                   ha="center")
            a.legend(handles=[
                Line2D([], [], marker="s", ls="", ms=10, color=col,
                       label="ce mode seul"),
                Line2D([], [], marker="s", ls="", ms=10, color=VIOLET,
                       label="mixte (les deux modes)"),
                Line2D([], [], marker="s", ls="", ms=10, color="0.85",
                       label="l'autre mode seul")], fontsize=8, loc="upper right")
        else:
            a.axhline(120, color="0.35", lw=1.2)
            a.text(60 - 0.95 * rad, 120.3, "surface du bloc", fontsize=8,
                   color="0.4")
        a.set_xlim(*xlim); a.set_ylim(*ylim); a.set_aspect("equal")
        a.set_xlabel(xl); a.set_ylabel(yl); a.set_title(t, fontsize=10.5)

    nprop = int((principal & ~mixte).sum())
    fig.suptitle("%s — enveloppe par ÉLÉMENTS : %d tétraèdres propres + %d "
                 "mixtes, volume %.2f mm³ (total endommagé %.2f mm³)"
                 % (titre, nprop, int(mixte.sum()), vol[principal].sum(),
                    vol[sel_all].sum()), fontsize=11.5)
    fig.tight_layout()
    out = os.path.join(HERE, "env_elem_%s_%s.png"
                       % (tag, RUN.replace("out_imp3d_", "")))
    fig.savefig(out, dpi=140)
    print("  ecrit", os.path.basename(out))
