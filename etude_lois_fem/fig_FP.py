# -*- coding: utf-8 -*-
"""Courbes force-penetration des percussions de la matrice (force filtree 10 us, penetration = descente
de l'insert depuis l'instant de contact toolDelay). 4 panneaux : R vs P ; ablations a P = 0 ; P = 100 R vs GIIc1 ;
convergence en maillage a P = 0."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_fem3d import read_history
plt.rcParams.update({"font.family": "serif", "font.serif": ["CMU Serif", "Latin Modern Roman", "DejaVu Serif"],
                     "mathtext.fontset": "cm", "font.size": 9, "axes.titlesize": 9, "legend.fontsize": 7.5})
HERE = os.path.dirname(os.path.abspath(__file__)); M = os.path.join(HERE, "matrice"); OUT = os.path.join(HERE, "figures")

def fp(name, t0, win=10e-6):
    head, h = read_history(os.path.join(M, "out_" + name, "history.csv")); col = {n: i for i, n in enumerate(head)}
    t = h[:, 0]; dt = np.median(np.diff(t)); n = max(1, int(round(win / dt)))
    F = np.convolve(h[:, col["toolFz"]] * 1e-3, np.ones(n) / n, mode="same")     # force de compression sur l'insert, positive
    ic = np.argmax((t >= t0) & (np.abs(h[:, col["toolFz"]]) > 200.0))              # premier contact (200 N brut)
    d = (h[ic, col["toolZ"]] - h[:, col["toolZ"]]) * 1e3                               # penetration depuis le contact
    sel = t >= t[ic]
    return d[sel], F[sel]

def draw(a, runs, title):
    for name, lab, c, t0, ls in runs:
        try: d, F = fp(name, t0)
        except Exception as e: print("manque", name, e); continue
        a.plot(d, F, color=c, lw=1.1, ls=ls, label=lab)
    a.set_xlabel("pénétration [mm]"); a.set_ylabel("force de contact filtrée 10 µs [kN]"); a.set_title(title); a.grid(alpha=0.25); a.legend(frameon=False)
    a.set_xlim(left=-0.02); a.set_ylim(bottom=-2); a.axhline(0, color="k", lw=0.4)

fig, ax = plt.subplots(2, 2, figsize=(8.6, 6.4))
draw(ax[0, 0], [("C_T1_R_P000", "P = 0", "C0", 0, "-"), ("D_T1_R_P050", "P = 50 MPa", "C2", 60e-6, "-"), ("C_T1_R_P100", "P = 100 MPa", "C3", 60e-6, "-")],
     "référence R : effet du confinement (aire = travail absorbé)")
draw(ax[0, 1], [("C_T1_R_P000", "R", "C0", 0, "-"), ("D_T1_lin_P000", "méridien linéaire", "C7", 0, "--"), ("D_T1_apex_P000", "cut-off à l'apex (ft 23,1)", "C1", 0, "-")],
     "ablations à P = 0")
draw(ax[1, 0], [("C_T1_R_P100", "R, G_IIc = 10 N/mm", "C3", 60e-6, "-"), ("C_T1_R_P100_GIIc1", "G_IIc = 1 N/mm", "C1", 60e-6, "-")],
     "P = 100 MPa : énergie de rupture en compression")
draw(ax[1, 1], [("C_T1_R_P000_h15", "cœur 1,5 mm", "0.8", 0, "-"), ("C_T1_R_P000_h10", "1,0 mm", "0.6", 0, "-"), ("C_T1_R_P000_h075", "0,75 mm", "0.4", 0, "-"),
                ("C_T1_R_P000", "0,5 mm (R)", "C0", 0, "-"), ("C_T1_R_P000_c035", "0,35 mm", "k", 0, "-")],
     "convergence en maillage, R à P = 0")
fig.suptitle("Percussion 16 J, insert R = 7,94 mm, Red Bohus : force–pénétration", y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig7_FP_percussion.pdf")); fig.savefig(os.path.join(OUT, "fig7_FP_percussion.png"), dpi=130)
print("-> fig7")
