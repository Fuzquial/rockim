# -*- coding: utf-8 -*-
"""UNE FIGURE PAR MODE de rupture, en 3D, avec l'enveloppe du bloc.

  mode I  (ouverture / traction)     -> BLEU
  mode II (glissement / cisaillement) -> ROUGE

Chaque figure porte deux vues : a gauche le mode seul, zoome a l'echelle des
fissures ; a droite le meme objet dans l'enveloppe du bloc 120^3, pour situer.
L'autre mode est rappele en gris translucide, pour que la comparaison soit
possible sans superposer les deux couleurs.

  python plot_mode_enveloppe.py [run]      (defaut out_imp3d_ultra)
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
RUN = sys.argv[1] if len(sys.argv) > 1 else "out_imp3d_ultra"
A_HERTZ, W = 1.45, 120.0
BLEU, ROUGE = "#1f6fb4", "#c0392b"

D = os.path.join(ROOT, RUN)
f = sorted(x for x in os.listdir(D)
           if re.fullmatch(r"fdem3d_joints_\d{4}\.vtu", x))[-1]
s = open(os.path.join(D, f), encoding="utf-8", errors="ignore").read()


def arr(n, d=float):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, s, re.S)
    return None if not m else np.fromstring(m.group(1), sep=" ", dtype=d)


P = np.fromstring(re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>",
                            s, re.S).group(1), sep=" ").reshape(-1, 3) * 1e3
conn, off = arr("connectivity", int), arr("offsets", int)
dmg, bm = arr("damage"), arr("breakMode")
idx = np.where(dmg >= 1.0)[0]
tri = np.array([P[conn[(0 if i == 0 else off[i - 1]):off[i]][:3]] for i in idx])
mode = bm[idx]
cen = tri.mean(axis=1)


def cube_edges(w):
    """aretes du bloc [0,w]^3 pour l'enveloppe."""
    c = np.array([[x, y, z] for x in (0, w) for y in (0, w) for z in (0, w)])
    seg = []
    for i in range(8):
        for j in range(i + 1, 8):
            if np.sum(np.abs(c[i] - c[j]) > 1e-9) == 1:
                seg.append([c[i], c[j]])
    return seg


def cercle(rad, z, n=80):
    th = np.linspace(0, 2 * np.pi, n)
    return np.column_stack([60 + rad * np.cos(th), 60 + rad * np.sin(th),
                            np.full(n, z)])


for v, col, nom, titre in ((1, BLEU, "modeI_traction",
                            "MODE I — ouverture / traction"),
                           (2, ROUGE, "modeII_cisaillement",
                            "MODE II — glissement / cisaillement")):
    sel, autre = mode == v, mode != v
    n, nA = int(sel.sum()), int(autre.sum())
    fig = plt.figure(figsize=(13.6, 6.4))

    for k in (0, 1):
        a = fig.add_subplot(1, 2, k + 1, projection="3d")
        # l'AUTRE mode, en gris translucide, pour situer
        if nA:
            pc = Poly3DCollection(tri[autre], facecolors="0.72", alpha=0.18,
                                  edgecolors="none")
            a.add_collection3d(pc)
        pc = Poly3DCollection(tri[sel], facecolors=col, alpha=0.9,
                              edgecolors="k", linewidths=0.15)
        a.add_collection3d(pc)

        if k == 1:                              # enveloppe du bloc entier
            a.add_collection3d(Line3DCollection(cube_edges(W), colors="0.55",
                                               linewidths=0.9))
            a.set_xlim(0, W); a.set_ylim(0, W); a.set_zlim(0, W)
            a.set_title("dans l'enveloppe du bloc 120³ mm", fontsize=10)
        else:                                   # zoom a l'echelle des fissures
            rad = max(3.0, 1.3 * np.hypot(cen[:, 0] - 60, cen[:, 1] - 60).max())
            dep = max(3.0, 1.3 * (120 - cen[:, 2]).max())
            a.plot(*cercle(A_HERTZ, 120.05).T, color="k", lw=1.8, ls="--")
            a.text(60, 60 + A_HERTZ, 120.4, "contact Hertz\na = %.2f mm" % A_HERTZ,
                   fontsize=8, ha="center")
            a.set_xlim(60 - rad, 60 + rad); a.set_ylim(60 - rad, 60 + rad)
            a.set_zlim(120 - dep, 120.6)
            a.set_title("zoom — %d facettes de ce mode" % n, fontsize=10)
        a.set_xlabel("x (mm)", fontsize=8.5); a.set_ylabel("y (mm)", fontsize=8.5)
        a.set_zlabel("z (mm)", fontsize=8.5)
        a.tick_params(labelsize=7.5)
        a.view_init(elev=22, azim=-58)
        try:
            a.set_box_aspect((1, 1, 0.8))
        except Exception:
            pass

    r = np.hypot(cen[sel][:, 0] - 60, cen[sel][:, 1] - 60)
    prof = 120 - cen[sel][:, 2]
    fig.suptitle("%s — %d facettes sur %d (%.0f %%)   |   rayon médian %.2f mm, "
                 "profondeur médiane %.2f mm, max %.2f mm\n(l'autre mode en gris "
                 "translucide)" % (titre, n, len(mode), 100 * n / len(mode),
                                   np.median(r), np.median(prof), prof.max()),
                 fontsize=11.5)
    fig.tight_layout()
    out = os.path.join(HERE, "mode_%s_%s.png"
                       % (nom, RUN.replace("out_imp3d_", "")))
    fig.savefig(out, dpi=140)
    print("ecrit %s  (%d facettes, rayon med %.2f mm, prof med %.2f, max %.2f)"
          % (os.path.basename(out), n, np.median(r), np.median(prof), prof.max()))
