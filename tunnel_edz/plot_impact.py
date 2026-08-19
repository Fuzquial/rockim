#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_impact.py — planche pour un run d'impact sur bloc plein (pas de contour
# de tunnel a tracer, contrairement a plot_tunnel_fields.py).
#
#   python tunnel_edz/plot_impact.py out_impact_tun [--zoom 12]
#
# (a) fissures par mode           (b) sigma_1        (c) tau_max
# (d) deformation de cisaillement (e) deplacement    (f) force outil vs temps
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

C_TEN, C_SHR = "#1B8A3A", "#C8342B"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_impact_tun")
    ap.add_argument("--zoom", type=float, default=12.0)
    ap.add_argument("--frame", type=int, default=-1)
    # rayon de l'outil [m] : trace un DEMI-CERCLE tangent a la surface, a
    # l'echelle du cas. (La premiere version dessinait un segment de +-1 m code
    # en dur pour le cas metrique — sur un bloc de 40 mm il barrait toute la
    # figure.) 0 = ne rien tracer.
    ap.add_argument("--toolr", type=float, default=0.0)
    ap.add_argument("--titre", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)
    out = os.path.join(HERE, os.path.basename(run) + "_planche.png")

    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    jn = [f for f in sorted(glob.glob(os.path.join(run, "fdem_joints_[0-9]*.vtu")))
          if complete(f)]
    k = min(len(el), len(jn)) - 1 if a.frame < 0 else a.frame
    P0, C, _ = read_vtu(el[0], [])
    P, _, S = read_vtu(el[k], ["sigmaXX", "sigmaYY", "sigmaXY"])
    W, H = P0[:, 0].max(), P0[:, 1].max()
    cx, cy = 0.5 * W, H                       # point d'impact : haut, centre

    sxx, syy, sxy = S["sigmaXX"], S["sigmaYY"], S["sigmaXY"]
    tr, dif = 0.5 * (sxx + syy), 0.5 * (sxx - syy)
    rad = np.hypot(dif, sxy)
    s1, tmax = (tr + rad) / 1e6, rad / 1e6
    X, x = P0[C], P[C]
    A = np.stack([X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]], axis=2)
    B = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]], axis=2)
    F = B @ np.linalg.inv(A)
    E = 0.5 * (np.transpose(F, (0, 2, 1)) @ F - np.eye(2))
    gmax = 200.0 * np.hypot(0.5 * (E[:, 0, 0] - E[:, 1, 1]), E[:, 0, 1])
    umag = np.linalg.norm(P - P0, axis=1)
    uel = umag[C].mean(axis=1)
    T0 = mtri.Triangulation(P0[:, 0], P0[:, 1], C)

    PJ, CJ, SJ = read_vtu(jn[k], ["damage", "breakMode"], ncell=2)
    bm = SJ["breakMode"]
    brk = (bm > 0) | (SJ["damage"] >= 0.999)
    seg, mode = PJ[CJ[brk]], bm[brk]

    Z = a.zoom
    # Echelles calees sur la ZONE D'IMPACT et non sur tout le domaine : le pic
    # est tres localise sous le poincon, et un percentile global affiche le
    # champ lointain en saturant l'interet (defaut de la premiere version).
    ctr = P0[C].mean(axis=1)
    near = np.hypot(ctr[:, 0] - cx, ctr[:, 1] - cy) < 0.6 * Z
    def lim(f, q=98.0):
        return np.percentile(np.abs(f[near]), q) if near.any() else np.abs(f).max()

    fig = plt.figure(figsize=(16.5, 9.2))
    panels = [(None, None, None, "(a) fissures"),
              (s1, "RdBu_r", dict(vmin=-lim(s1), vmax=lim(s1)),
               r"(b) $\sigma_1$ [MPa]"),
              (tmax, "inferno", dict(vmax=lim(tmax)),
               r"(c) $\tau_{max}$ [MPa]"),
              (gmax, "magma", dict(vmax=lim(gmax, 99.0)),
               r"(d) $\gamma_{max}$ [%]"),
              (uel, "turbo", {}, f"(e) |u| [m] — max {umag.max():.4f} m")]
    for i, (fld, cmap, kw, ttl) in enumerate(panels):
        ax = fig.add_subplot(2, 3, i + 1)
        if fld is None:
            ax.add_collection(LineCollection(seg[mode == 2], colors=C_SHR, lw=1.2))
            ax.add_collection(LineCollection(seg[mode == 1], colors=C_TEN, lw=1.2))
            nt, ns = int((mode == 1).sum()), int((mode == 2).sum())
            ttl += f" : {nt} traction, {ns} cisaillement"
            ax.plot([], [], color=C_TEN, lw=3, label="traction")
            ax.plot([], [], color=C_SHR, lw=3, label="cisaillement")
            ax.legend(loc="lower left", fontsize=8)
        else:
            h = ax.tripcolor(T0, facecolors=fld, cmap=cmap, **kw)
            fig.colorbar(h, ax=ax, shrink=0.82)
        if a.toolr > 0.0:
            th = np.linspace(np.pi, 2 * np.pi, 120)
            ax.plot(cx + a.toolr * np.cos(th), cy + a.toolr + a.toolr * np.sin(th),
                    color="#0B4F9E", lw=2.2)
        ax.set_aspect("equal")
        ax.set_xlim(cx - Z, cx + Z)
        ax.set_ylim(cy - 1.6 * Z, cy + 0.15 * Z)
        ax.set_title(ttl, fontsize=10.5)
        ax.set_xlabel("x [m]")
        if i % 3 == 0:
            ax.set_ylabel("y [m]")

    ax = fig.add_subplot(2, 3, 6)
    d = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",",
                      names=True, invalid_raise=False)
    ax.plot(d["t"] * 1e3, np.abs(d["toolFy"]) / 1e6, color="#0B4F9E", lw=1.8)
    ax2 = ax.twinx()
    ax2.plot(d["t"] * 1e3, d["nBroken"], color=C_SHR, lw=1.5)
    ax2.set_ylabel("joints rompus", color=C_SHR)
    ax.set_xlabel("temps [ms]")
    ax.set_ylabel("force outil |Fy| [MN/m]", color="#0B4F9E")
    ax.set_title("(f) force d'impact et fissuration", fontsize=10.5)
    ax.grid(alpha=0.3)

    fig.suptitle(a.titre or os.path.basename(run), fontsize=12.5)
    fig.tight_layout()
    fig.savefig(out, dpi=155)
    print("ecrit :", out)


if __name__ == "__main__":
    main()
