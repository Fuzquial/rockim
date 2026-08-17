# -*- coding: utf-8 -*-
"""Les joints FISSURES des trois impacts 3D : ou l'endommagement s'installe,
et jusqu'ou il va. Vue de dessus + coupe, une colonne par run.

Un joint est trace des que D > 0 ; la couleur est D (0 = intact, 1 = rompu).
Les joints pleinement rompus sont cercles en noir.

  python plot_joints_fissures.py
"""
import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

RUNS = [("out_banc_mid",    "14/08 — insertion intrinsèque\n(tous les joints existent)"),
        ("out_imp3d_homog", "17/08 — adaptative, homogène\n(82 joints insérés)"),
        ("out_imp3d_weib",  "17/08 — adaptative + Weibull\n(132 joints insérés)")]


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

    for row, cols, xl, yl, xlim, ylim in (
            (0, [0, 1], "x (mm)", "y (mm)", (25, 95), (25, 95)),
            (1, [0, 2], "x (mm)", "z (mm)", (25, 95), (95, 122))):
        a = ax[row, j]
        if len(tri):
            pc = a.add_collection(PolyCollection(
                tri[:, :, cols], array=dmg[sel], cmap="YlOrRd",
                edgecolors="0.4", linewidths=0.3))
            pc.set_clim(0, 1)
            if brk.any():
                a.scatter(cen[brk][:, cols[0]], cen[brk][:, cols[1]], s=90,
                          facecolors="none", edgecolors="k", lw=1.4, zorder=6,
                          label="rompus (D = 1)")
        if row == 0:
            a.add_patch(plt.Circle((60, 60), 11, fill=False, ec="C0", lw=1.6,
                                   ls="--"))
            a.text(60, 73, "insert R = 11 mm", color="C0", fontsize=8,
                   ha="center")
            a.set_title("%s\n%d joints endommagés, %d rompus"
                        % (lab, len(sel), int(brk.sum())), fontsize=9.5)
        else:
            a.axhline(120, color="0.5", lw=1.0)
            a.text(27, 120.6, "surface du bloc", color="0.4", fontsize=7.5)
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
