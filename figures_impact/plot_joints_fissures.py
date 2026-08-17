# -*- coding: utf-8 -*-
"""Les joints FISSURES des trois impacts 3D : ou l'endommagement s'installe,
et jusqu'ou il va. Vue de dessus + coupe, une colonne par run.

Un joint est trace des que D > 0 ; la couleur est D (0 = intact, 1 = rompu).
Les joints pleinement rompus sont cercles en noir.

  python plot_joints_fissures.py
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

DEFAUT = [("out_imp3d_homog", "ft = 10 MPa, maillage 4,7 mm\n$\\ell_{cz}$ = 35 mm — HORS fenêtre"),
          ("out_banc_mid",    "ft = 10 MPa, insertion intrinsèque\n$\\ell_{cz}$ = 35 mm — hors fenêtre"),
          ("out_imp3d_ultra", "ft = 87 MPa, maillage 0,46 mm\n$\\ell_{cz}$ = 1,0 mm — DANS la fenêtre")]
# un ou plusieurs noms de run en argument remplacent la selection par defaut
RUNS = [(r, r.replace("out_", "")) for r in sys.argv[1:]] if len(sys.argv) > 1 else DEFAUT


def txt(p):
    return open(p, encoding="utf-8", errors="ignore").read()


def arr(s, n, d=float):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, s, re.S)
    return None if not m else np.fromstring(m.group(1), sep=" ", dtype=d)


def pts_of(s):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", s, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)


fig, ax = plt.subplots(2, 3, figsize=(15.5, 9.2))
for j, (sub, lab) in enumerate(RUNS):
    D = os.path.join(ROOT, sub)
    fr = sorted(f for f in os.listdir(D)
                if re.fullmatch(r"fdem3d_joints_\d{4}\.vtu", f))[-1]
    s = txt(os.path.join(D, fr))
    P = pts_of(s) * 1e3
    conn = arr(s, "connectivity", int)
    off = arr(s, "offsets", int)
    dmg = arr(s, "damage")
    tbr = arr(s, "tBreak")

    sel = np.where(dmg > 0)[0]
    tri = []
    for i in sel:
        a0 = 0 if i == 0 else off[i - 1]
        tri.append(P[conn[a0:off[i]][:3]])
    tri = np.array(tri) if len(tri) else np.zeros((0, 3, 3))
    brk = dmg[sel] >= 1.0
    cen = tri.mean(axis=1) if len(tri) else np.zeros((0, 3))

    # CADRAGE PROPRE A CHAQUE RUN : la fissuration d'un run bien resolu tient
    # dans quelques millimetres (le contact fait 1,5 mm de rayon) et
    # apparaissait comme un point dans une fenetre fixe de 70 mm.
    if len(cen):
        rad = max(4.0, 1.35 * np.hypot(cen[:, 0] - 60, cen[:, 1] - 60).max())
        dep = max(4.0, 1.35 * (120 - cen[:, 2].min()))
    else:
        rad, dep = 35.0, 25.0
    for row, cols, xl, yl, xlim, ylim in (
            (0, [0, 1], "x (mm)", "y (mm)", (60 - rad, 60 + rad), (60 - rad, 60 + rad)),
            (1, [0, 2], "x (mm)", "z (mm)", (60 - rad, 60 + rad), (120 - dep, 122))):
        a = ax[row, j]
        if len(tri):
            pc = a.add_collection(PolyCollection(
                tri[:, :, cols], array=dmg[sel], cmap="YlOrRd",
                edgecolors="0.4", linewidths=0.3))
            pc.set_clim(0, 1)
            if brk.any():
                # taille du marqueur ADAPTEE au nombre : dimensionne pour 4
                # joints, il fusionnait en tache noire illisible a 347
                n = int(brk.sum())
                ms = 90.0 if n <= 20 else (25.0 if n <= 100 else 6.0)
                lw = 1.4 if n <= 20 else 0.5
                a.scatter(cen[brk][:, cols[0]], cen[brk][:, cols[1]], s=ms,
                          facecolors="none", edgecolors="k", lw=lw, zorder=6,
                          label="rompus (D = 1) — %d" % n)
        if row == 0:
            if rad > 12:
                a.add_patch(plt.Circle((60, 60), 11, fill=False, ec="C0",
                                       lw=1.6, ls="--"))
                a.text(60, 73, "insert R = 11 mm", color="C0", fontsize=8,
                       ha="center")
            else:   # zoom serre : on repere le rayon de contact de Hertz
                a.add_patch(plt.Circle((60, 60), 1.45, fill=False, ec="C0",
                                       lw=1.6, ls="--"))
                a.text(60, 60 + 0.22 * rad, "contact Hertz a = 1,45 mm",
                       color="C0", fontsize=8, ha="center")
            a.set_title("%s\n%d joints endommagés, %d rompus"
                        % (lab, len(sel), int(brk.sum())), fontsize=9.5)
        else:
            a.axhline(120, color="0.5", lw=1.0)
            a.text(60 - 0.95 * rad, 120.6, "surface du bloc", color="0.4", fontsize=7.5)
            a.set_title("coupe verticale — profondeur atteinte : %.1f mm"
                        % (120 - cen[:, 2].min() if len(cen) else 0),
                        fontsize=9.5)
        a.set_xlim(*xlim); a.set_ylim(*ylim); a.set_aspect("equal")
        a.set_xlabel(xl); a.set_ylabel(yl)
        if row == 0 and j == 0 and brk.any():
            a.legend(fontsize=8, loc="upper right")

cb = fig.colorbar(pc, ax=ax.ravel().tolist(), shrink=0.55, pad=0.02)
cb.set_label("endommagement D du joint  (1 = rompu)")
fig.suptitle("Les joints fissurés — impact insert/granite à 8 m/s, "
             "état final (l'endommagement est irréversible)", fontsize=13)
out = os.path.join(HERE, "joints_fissures.png")
fig.savefig(out, dpi=140, bbox_inches="tight")
print("ecrit", out)
