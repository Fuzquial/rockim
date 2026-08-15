# -*- coding: utf-8 -*-
"""Figures de la campagne : sensibilite (phase B), maillage GBM et facies."""
import csv, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)
PARAMS = ["ft", "cohesion", "frictionDeg", "Gf", "gfShearFactor", "crushCap"]
TARGET = {"ucs": 126.6, "bts": 10.27, "tx20": 424.8}


def arr(t, n):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    return np.fromstring(m.group(1), sep=" ") if m else None


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


# ---------------------------------------------------------------- sensibilite
def fig_sensitivity():
    rows = {r["tag"]: r for r in csv.DictReader(open(
        os.path.join(BASE, "screen_results.csv")))}
    C = rows["C"]
    cols = {"ucs": "ucs_peak_MPa", "bts": "bts_sigma_t_MPa",
            "tx20": "tx20_peak_MPa"}
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.6))
    for a, (k, col) in zip(ax, cols.items()):
        c = float(C[col])
        y = np.arange(len(PARAMS))
        lo = [float(rows[p + "_lo"][col]) for p in PARAMS]
        hi = [float(rows[p + "_hi"][col]) for p in PARAMS]
        for i, (l, h) in enumerate(zip(lo, hi)):
            a.plot([l, h], [i, i], "-", color="0.75", lw=6, solid_capstyle="butt")
            a.plot(l, i, "o", color="C0", ms=7)
            a.plot(h, i, "o", color="C3", ms=7)
        a.axvline(c, color="k", lw=1, ls="--", label="centre")
        a.axvline(TARGET[k], color="C2", lw=2, label="cible exp.")
        a.set_yticks(y); a.set_yticklabels(PARAMS if k == "ucs" else [])
        a.set_xlabel({"ucs": "UCS (MPa)", "bts": r"$\sigma_t$ (MPa)",
                      "tx20": r"$\sigma_1$ à $\sigma_3$ = 20 MPa (MPa)"}[k])
        a.set_title(k.upper())
        a.grid(alpha=0.3, axis="x")
        if k == "ucs":
            a.legend(fontsize=8, loc="lower right")
    fig.suptitle("Phase B — criblage : effet de chaque paramètre porté à ses bornes "
                 "(bleu = borne basse, rouge = borne haute)")
    fig.tight_layout()
    p = os.path.join(FIG, "phaseB_sensibilite.png")
    fig.savefig(p, dpi=140)
    print("ecrit", p)


# ------------------------------------------------------------- maillage + facies
def panel(ax, run, frame, title, show_grain=True):
    d = os.path.join(BASE, "runs", run)
    t0 = open(os.path.join(d, "fdem_0000.vtu")).read()
    P0, conn = pts(t0), arr(t0, "connectivity").astype(int).reshape(-1, 3)
    grain = arr(t0, "grain")
    f = os.path.join(d, "fdem_%04d.vtu" % frame)
    if not os.path.exists(f):
        f = os.path.join(d, "fdem_0000.vtu")
    P = pts(open(f).read())
    polys = [P[c] * 1e3 for c in conn]
    if show_grain and grain is not None:
        rng = np.random.default_rng(3)
        lut = rng.random((int(grain.max()) + 2, 3)) * 0.45 + 0.5
        cols = [lut[int(g)] for g in grain]
    else:
        cols = ["0.85"] * len(polys)
    ax.add_collection(PolyCollection(polys, facecolors=cols,
                                     edgecolors="0.35", linewidths=0.15))
    # joints rompus (tBreak >= 0) sur la meme frame
    jf = os.path.join(d, "fdem_joints_%04d.vtu" % frame)
    if os.path.exists(jf):
        jt = open(jf).read()
        jP, jc = pts(jt), arr(jt, "connectivity").astype(int)
        tb = arr(jt, "tBreak")
        off = arr(jt, "offsets").astype(int)
        segs, s = [], 0
        for i, e in enumerate(off):
            idx = jc[s:e]; s = e
            if tb is not None and tb[i] >= 0 and len(idx) >= 2:
                segs.append(jP[idx[:2]] * 1e3)
        if segs:
            ax.add_collection(LineCollection(segs, colors="crimson", lw=1.4))
            title += "  (%d joints rompus)" % len(segs)
    ax.set_xlim(P0[:, 0].min() * 1e3 - 1, P0[:, 0].max() * 1e3 + 1)
    ax.set_ylim(P0[:, 1].min() * 1e3 - 1, P0[:, 1].max() * 1e3 + 1)
    ax.set_aspect("equal"); ax.set_title(title, fontsize=10)
    ax.set_xlabel("x (mm)")


def fig_mesh_facies():
    fig, ax = plt.subplots(1, 3, figsize=(13, 6.2))
    panel(ax[0], "C_ucs_s4211", 0, "Maillage GBM Voronoï\n(grains 5 mm, Delaunay 2 mm)")
    panel(ax[1], "C_ucs_s4211", 3, "UCS — faciès de rupture")
    panel(ax[2], "C_bts_s4211", 3, "Brésilien — faciès")
    ax[0].set_ylabel("y (mm)")
    fig.suptitle("Modèle de calibration Red Bohus — maillage à grains de Voronoï "
                 "et faciès simulés (jeu central)")
    fig.tight_layout()
    p = os.path.join(FIG, "modele_et_facies.png")
    fig.savefig(p, dpi=140)
    print("ecrit", p)


if __name__ == "__main__":
    fig_sensitivity()
    fig_mesh_facies()
