# -*- coding: utf-8 -*-
"""Trois figures d'aperçu (PDF vectoriel + PNG) depuis les sorties existantes :
  1. F(t) filtré et enfoncement, phase A (élastique, dpr ref, apex, dpdfh) + référence R de la matrice
  2. partition de la dissipation de R (wPlas, wDamT, wDamC) et érosion
  3. coupe du cratère de R (rendu éléments, tranche |y - D/2| < 0,8 mm) : D, omega_c, érodés"""
import io, os, sys, re, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_fem3d import read_history, read_vtu_cells

plt.rcParams.update({"font.family": "serif", "font.serif": ["CMU Serif", "Latin Modern Roman", "DejaVu Serif"],
                     "mathtext.fontset": "cm", "font.size": 9, "axes.titlesize": 9, "legend.fontsize": 8})
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)

def smooth(t, y, win=10e-6):
    if len(t) < 3: return y
    dt = np.median(np.diff(t)); n = max(1, int(round(win / dt)))
    k = np.ones(n) / n
    return np.convolve(y, k, mode="same")

runs = [("phaseA/out_A1_elastic", "élastique (A1)", "0.5"), ("phaseA/out_A2_dpr_ref", "dpr, ft 9, érosion héritée (A2)", "C3"),
        ("phaseA/out_A3_dpr_apex", "dpr, cut-off à l'apex 23 MPa (A3)", "C1"), ("phaseA/out_A7_dpdfh", "dpdfh (A7)", "C2"),
        ("matrice/out_C_T1_R_P000", "référence R de la matrice (P = 0)", "k")]
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))
for d, lab, c in runs:
    p = os.path.join(HERE, d, "history.csv")
    if not os.path.exists(p): continue
    head, h = read_history(p); col = {n: i for i, n in enumerate(head)}
    t = h[:, 0] * 1e6; fz = -h[:, col["toolFz"]] * 1e-3; z = (h[0, col["toolZ"]] - h[:, col["toolZ"]]) * 1e3
    ax[0].plot(t, smooth(h[:, 0], fz), color=c, lw=1.1, label=lab)
    ax[1].plot(t, z, color=c, lw=1.1)
ax[0].set_xlabel("t [µs]"); ax[0].set_ylabel("force sur l'insert, filtrée 10 µs [kN]"); ax[0].set_xlim(0, 300)
ax[1].set_xlabel("t [µs]"); ax[1].set_ylabel("enfoncement de l'insert [mm]"); ax[1].set_xlim(0, 300)
ax[0].legend(frameon=False, loc="upper right");
for a in ax: a.grid(alpha=0.25)
fig.suptitle("Percussion R = 7,94 mm, 16 J, Red Bohus harmonisé, P = 0 : quatre lois, une référence", y=1.02)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1_Ft_phaseA_R.pdf")); fig.savefig(os.path.join(OUT, "fig1_Ft_phaseA_R.png"), dpi=130)

# --- figure 2 : partition de la dissipation de R ---
p = os.path.join(HERE, "matrice/out_C_T1_R_P000/history.csv")
head, h = read_history(p); col = {n: i for i, n in enumerate(head)}
t = h[:, 0] * 1e6
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))
E0 = 0.5 * 0.5 * 8.0 ** 2
ax[0].plot(t, h[:, col["work"]], "k", lw=1.2, label="travail de l'outil")
ax[0].plot(t, h[:, col["wPlas"]], "C3", lw=1.1, label="plastique (cône DP puissance)")
ax[0].plot(t, h[:, col["wDamC"]], "C1", lw=1.1, label="endommagement compressif ω_c")
ax[0].plot(t, h[:, col["wDamT"]], "C0", lw=1.1, label="endommagement tractif D")
ax[0].plot(t, h[:, col["toolKE"]], "0.5", lw=1.0, ls="--", label="énergie cinétique de l'insert")
ax[0].axhline(E0, color="0.8", lw=0.8); ax[0].text(2, E0 * 0.97, "16 J", va="top", fontsize=7, color="0.4")
ax[0].set_xlabel("t [µs]"); ax[0].set_ylabel("énergie [J]"); ax[0].legend(frameon=False, fontsize=7); ax[0].grid(alpha=0.25)
ax[1].plot(t, h[:, col["nEroSpall"]], "C0", lw=1.1, label="érodés : spall (98 % de Gf dissipés)")
ax[1].plot(t, h[:, col["nEroGeo"]], "C3", lw=1.1, label="érodés : soupape géométrique")
ax[1].plot(t, h[:, col["V_D09"]] * 1e9 / 10, "C0", ls=":", lw=1.0, label="V(D ≥ 0,9) / 10 [mm³]")
ax[1].plot(t, h[:, col["V_wc05"]] * 1e9 / 10, "C1", ls=":", lw=1.0, label="V(ω_c ≥ 0,5) / 10 [mm³]")
ax[1].set_xlabel("t [µs]"); ax[1].set_ylabel("éléments / volumes"); ax[1].legend(frameon=False, fontsize=7); ax[1].grid(alpha=0.25)
fig.suptitle("Référence R à P = 0 : où va l'énergie, et ce que la suppression retire", y=1.02)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_partition_R_P0.pdf")); fig.savefig(os.path.join(OUT, "fig2_partition_R_P0.png"), dpi=130)

# --- figure 3 : coupe du cratère en rendu éléments ---
def slab_polys(vtu, ycut, halfw):
    txt = io.open(vtu, encoding="utf-8", errors="ignore").read()
    m = re.search(r"<Points>\s*<DataArray[^>]*>\s*(.*?)\s*</DataArray>", txt, re.S)
    pts = np.fromstring(m.group(1), sep=" ").reshape(-1, 3)
    con, vol, f = read_vtu_cells(vtu)
    cen = pts[con].mean(axis=1)
    sel = np.abs(cen[:, 1] - ycut) < halfw
    polys = []
    for tet in con[sel]:
        P = pts[tet][:, [0, 2]] * 1e3
        # enveloppe convexe des 4 points projetés (ordre angulaire)
        c = P.mean(axis=0); ang = np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0])
        polys.append(P[np.argsort(ang)])
    return polys, {k: v[sel] for k, v in f.items()}

vtus = sorted(glob.glob(os.path.join(HERE, "matrice/out_C_T1_R_P000/fem3d_*.vtu")))
if vtus:
    polys, f = slab_polys(vtus[-1], 0.024, 0.0008)
    fig, ax = plt.subplots(1, 3, figsize=(8.4, 2.6), sharey=True)
    for a, key, title, cmap, vmax in ((ax[0], "damage", "endommagement tractif D", "Blues", 1.0),
                                      (ax[1], "omegaC", "endommagement compressif ω_c", "Oranges", 1.0),
                                      (ax[2], "epvEq", "déformation plastique équivalente", "Reds", 0.3)):
        vals = f[key].copy(); ero = f["eroded"] > 0.5
        pc = PolyCollection(polys, array=np.where(ero, np.nan, vals), cmap=cmap, edgecolors="none", clim=(0, vmax))
        a.add_collection(pc)
        if ero.any():
            a.add_collection(PolyCollection([p for p, e in zip(polys, ero) if e], facecolors="k", edgecolors="none"))
        a.set_xlim(8, 40); a.set_ylim(14, 33); a.set_aspect("equal"); a.set_title(title)
        a.set_xlabel("x [mm]"); fig.colorbar(pc, ax=a, fraction=0.046, pad=0.02)
    ax[0].set_ylabel("z [mm]")
    fig.suptitle("Référence R, P = 0, t = 300 µs : coupe |y − 24| < 0,8 mm sous l'insert (noir = éléments supprimés)", y=1.02)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_coupe_R_P0.pdf")); fig.savefig(os.path.join(OUT, "fig3_coupe_R_P0.png"), dpi=130)
print("figures ->", OUT)
