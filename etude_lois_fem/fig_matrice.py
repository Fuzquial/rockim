# -*- coding: utf-8 -*-
"""Figures de la phase C (percussion, maillages Delaunay, référence R) :
  4. convergence en maillage de la référence (P = 0 et 100) : W absorbé, F pic, restitution, enfoncement
  5. effet du confinement : F(t), enfoncement, partition de l'énergie — R à P = 0 / P = 100 / P = 100 avec G_IIc = 1
  6. coupes sous l'insert (rendu éléments) : ω_c et déformation plastique, R P = 0 / R P = 100 / G_IIc = 1 P = 100"""
import io, os, sys, re, glob, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_fem3d import read_history, read_vtu_cells

plt.rcParams.update({"font.family": "serif", "font.serif": ["CMU Serif", "Latin Modern Roman", "DejaVu Serif"],
                     "mathtext.fontset": "cm", "font.size": 9, "axes.titlesize": 9, "legend.fontsize": 7.5})
HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.join(HERE, "matrice"); OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
rows = {r["run"].replace("out_", ""): r for r in csv.DictReader(io.open(os.path.join(M, "resultats.csv"), encoding="utf-8"))}
def g(name, key):
    try: return float(rows[name][key])
    except Exception: return float("nan")

# ---- figure 4 : convergence ---------------------------------------------------
cores = [("h15", 1.5), ("h10", 1.0), ("h075", 0.75), ("", 0.5), ("c035", 0.35)]
fig, ax = plt.subplots(1, 4, figsize=(9.2, 2.5))
for P, c in (("P000", "C0"), ("P100", "C3")):
    xs, W, F, e, d = [], [], [], [], []
    for tag, h in cores:
        name = "C_T1_R_%s%s" % (P, ("_" + tag) if tag else "")
        if name not in rows: continue
        xs.append(h); W.append(g(name, "W_abs_J")); F.append(g(name, "F_pic_kN")); e.append(g(name, "e_r")); d.append(g(name, "dmax_mm"))
    for a, y, lab in zip(ax, (W, F, e, d), ("travail absorbé [J]", "force pic [kN]", "restitution", "enfoncement max [mm]")):
        a.plot(xs, y, "o-", color=c, lw=1.1, ms=4, label="P = %d MPa" % int(P[1:]))
        a.set_xlabel("taille de maille sous l'insert [mm]"); a.set_ylabel(lab); a.invert_xaxis(); a.grid(alpha=0.25)
    # réalisations à 0,5 mm
    for k, (key, a) in enumerate(zip(("W_abs_J", "F_pic_kN", "e_r", "dmax_mm"), ax)):
        for sd in (2, 3):
            v = g("C_T1_R_%s_s%d" % (P, sd), key)
            if not np.isnan(v): a.plot([0.5], [v], "x", color=c, ms=5)
ax[0].legend(frameon=False)
fig.suptitle("Convergence en maillage de la référence R (Delaunay gradué ; croix = réalisations à 0,5 mm)", y=1.03)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig4_convergence_R.pdf")); fig.savefig(os.path.join(OUT, "fig4_convergence_R.png"), dpi=130)

# ---- figure 5 : P = 0 vs 100 vs GIIc1 -------------------------------------------
def smooth(t, y, win=10e-6):
    dt = np.median(np.diff(t)); n = max(1, int(round(win / dt))); return np.convolve(y, np.ones(n) / n, mode="same")
runs = [("C_T1_R_P000", "R, P = 0", "C0", 0), ("C_T1_R_P100", "R, P = 100 MPa", "C3", 60e-6), ("C_T1_R_P100_GIIc1", "G_IIc = 1 N/mm, P = 100", "C1", 60e-6)]
fig, ax = plt.subplots(1, 3, figsize=(9.2, 2.6))
for name, lab, c, t0 in runs:
    head, h = read_history(os.path.join(M, "out_" + name, "history.csv")); col = {n: i for i, n in enumerate(head)}
    t = (h[:, 0] - t0) * 1e6; sel = t >= -5
    fz = -h[:, col["toolFz"]] * 1e-3
    ax[0].plot(t[sel], smooth(h[:, 0], fz)[sel], color=c, lw=1.1, label=lab)
    i0 = np.argmax(h[:, 0] >= t0)
    ax[1].plot(t[sel], (h[i0, col["toolZ"]] - h[:, col["toolZ"]])[sel] * 1e3, color=c, lw=1.1)
    for key, ls in (("wPlas", "-"), ("wDamC", "--"), ("wDamT", ":")):
        ax[2].plot(t[sel], h[sel, col[key]], color=c, lw=1.1, ls=ls)
ax[0].set_xlabel("t après contact [µs]"); ax[0].set_ylabel("force filtrée 10 µs [kN]"); ax[0].legend(frameon=False); ax[0].set_xlim(0, 300)
ax[1].set_xlabel("t après contact [µs]"); ax[1].set_ylabel("enfoncement [mm]"); ax[1].set_xlim(0, 300)
ax[2].set_xlabel("t après contact [µs]"); ax[2].set_ylabel("énergie dissipée [J]"); ax[2].set_xlim(0, 300)
ax[2].plot([], [], "k-", label="plastique"); ax[2].plot([], [], "k--", label="ω_c (compression)"); ax[2].plot([], [], "k:", label="D (traction)"); ax[2].legend(frameon=False)
for a in ax: a.grid(alpha=0.25)
fig.suptitle("Percussion 16 J, R = 7,94 mm, cœur 0,5 mm : le confinement rend la roche presque élastique, sauf si G_IIc est faible", y=1.03)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig5_confinement_R.pdf")); fig.savefig(os.path.join(OUT, "fig5_confinement_R.png"), dpi=130)

# ---- figure 6 : coupes ---------------------------------------------------------
def slab(vtu, ycut, halfw):
    txt = io.open(vtu, encoding="utf-8", errors="ignore").read()
    m = re.search(r"<Points>\s*<DataArray[^>]*>\s*(.*?)\s*</DataArray>", txt, re.S)
    pts = np.fromstring(m.group(1), sep=" ").reshape(-1, 3)
    con, vol, f = read_vtu_cells(vtu)
    cen = pts[con].mean(axis=1); sel = np.abs(cen[:, 1] - ycut) < halfw
    polys, big = [], []
    for tet in con[sel]:
        Q = pts[tet] * 1e3; big.append(max(np.linalg.norm(Q[i] - Q[j]) for i in range(4) for j in range(i)) > 2.0)   # noeuds envoles (element supprime)
        P = Q[:, [0, 2]]; c = P.mean(axis=0); ang = np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]); polys.append(P[np.argsort(ang)])
    fs = {k: v[sel] for k, v in f.items()}; fs["big"] = np.array(big)
    return polys, fs
fig, ax = plt.subplots(2, 3, figsize=(9.2, 5.0), sharex=True, sharey=True)
for j, (name, lab, _, _) in enumerate(runs):
    vtus = sorted(glob.glob(os.path.join(M, "out_" + name, "fem3d_*.vtu")))
    if not vtus: continue
    polys, f = slab(vtus[-1], 0.024, 0.0004)
    ero = f["eroded"] > 0.5; big = f["big"]
    for i, (key, title, cmap, vmax) in enumerate((("omegaC", "ω_c", "Oranges", 1.0), ("epvEq", "ε̄ᵖ", "Reds", 0.3))):
        a = ax[i, j]
        pc = PolyCollection([p for p, b in zip(polys, big) if not b], array=np.where(ero, np.nan, f[key])[~big], cmap=cmap, edgecolors="none", clim=(0, vmax)); a.add_collection(pc)
        if ero.any(): a.add_collection(PolyCollection([p for p, e, b in zip(polys, ero, big) if e and not b], facecolors="k", edgecolors="none"))
        a.set_xlim(14, 34); a.set_ylim(22, 33); a.set_aspect("equal"); a.set_title("%s — %s" % (lab, title)); fig.colorbar(pc, ax=a, fraction=0.04, pad=0.02)
        if i == 1: a.set_xlabel("x [mm]")
        if j == 0: a.set_ylabel("z [mm]")
fig.suptitle("Coupes |y − 24| < 0,4 mm sous l'insert à la fin du run (noir = éléments supprimés)", y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig6_coupes_R.pdf")); fig.savefig(os.path.join(OUT, "fig6_coupes_R.png"), dpi=130)
print("figures 4-6 ->", OUT)
