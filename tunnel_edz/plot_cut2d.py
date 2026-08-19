#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_cut2d.py — planche d'un essai de coupe 2D au cutter PDC.
#
#   python tunnel_edz/plot_cut2d.py out_cut2d_ye [--rake 20] [--depth 0.002]
#
# (a) fissures par mode + silhouette du cutter a sa position finale
# (b) cisaillement maximal    (c) deformation de cisaillement (zone plastique)
# (d) FORCES DE COUPE contre la distance parcourue — la courbe qui dit s'il y
#     a un REGIME DE COUPE : des cycles montee/effondrement signent le
#     detachement de copeaux successifs, une force plate signe du broyage.
# (e) rapport tangentiel/normal, l'observable de l'essai de rayure
# (f) fissuration cumulee
#
# Le cutter est trace a sa position courante : arete de coupe + face inclinee
# de l'angle de coupe arriere.
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
    ap.add_argument("run", nargs="?", default="out_cut2d_ye")
    ap.add_argument("--rake", type=float, default=20.0)
    ap.add_argument("--depth", type=float, default=0.002)
    ap.add_argument("--face", type=float, default=0.013)
    ap.add_argument("--frame", type=int, default=-1)
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

    sxx, syy, sxy = S["sigmaXX"], S["sigmaYY"], S["sigmaXY"]
    rad = np.hypot(0.5 * (sxx - syy), sxy)
    tmax = rad / 1e6
    X, x = P0[C], P[C]
    A = np.stack([X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]], axis=2)
    B = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]], axis=2)
    F = B @ np.linalg.inv(A)
    E = 0.5 * (np.transpose(F, (0, 2, 1)) @ F - np.eye(2))
    gmax = 200.0 * np.hypot(0.5 * (E[:, 0, 0] - E[:, 1, 1]), E[:, 0, 1])
    T0 = mtri.Triangulation(P0[:, 0], P0[:, 1], C)

    PJ, CJ, SJ = read_vtu(jn[k], ["damage", "breakMode"], ncell=2)
    bm = SJ["breakMode"]
    brk = (bm > 0) | (SJ["damage"] >= 0.999)
    seg, mode = PJ[CJ[brk]], bm[brk]

    d = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",",
                      names=True, invalid_raise=False)
    xs = (d["toolX"] - d["toolX"][0]) * 1e3          # course [mm]
    fx, fy = np.abs(d["toolFx"]) / 1e6, np.abs(d["toolFy"]) / 1e6
    # position du cutter a la trame tracee
    tf = d["t"][-1] * k / max(len(el) - 1, 1)
    xc = np.interp(tf, d["t"], d["toolX"])
    b = np.radians(a.rake)
    edge = np.array([xc, H - a.depth])
    face = edge + a.face * np.array([-np.sin(b), np.cos(b)])

    fig = plt.figure(figsize=(16.5, 9.0))
    for i, (fld, cmap, ttl) in enumerate(
            [(None, None, "(a) fissures"),
             (tmax, "inferno", r"(b) $\tau_{max}$ [MPa]"),
             (gmax, "magma", r"(c) $\gamma_{max}$ [%] — zone plastique")]):
        ax = fig.add_subplot(2, 3, i + 1)
        if fld is None:
            ax.add_collection(LineCollection(seg[mode == 2], colors=C_SHR, lw=0.7))
            ax.add_collection(LineCollection(seg[mode == 1], colors=C_TEN, lw=0.7))
            nt, ns = int((mode == 1).sum()), int((mode == 2).sum())
            ttl += f" : {nt} traction, {ns} cisaillement"
        else:
            near = P0[C].mean(axis=1)[:, 1] > H - 0.004
            h = ax.tripcolor(T0, facecolors=fld, cmap=cmap,
                             vmax=np.percentile(fld[near], 98))
            fig.colorbar(h, ax=ax, shrink=0.8)
        ax.plot([edge[0], face[0]], [edge[1], face[1]], color="#0B4F9E", lw=2.5)
        ax.plot(edge[0], edge[1], "o", color="#0B4F9E", ms=5)
        ax.plot([0, W], [H, H], color="0.6", lw=0.8)
        ax.set_aspect("equal")
        ax.set_xlim(0, W)
        ax.set_ylim(H - 0.006, H + 0.004)
        ax.set_title(ttl, fontsize=10.5)
        ax.set_xlabel("x [m]")
        if i == 0:
            ax.set_ylabel("y [m]")

    ax = fig.add_subplot(2, 3, 4)
    ax.plot(xs, fx, color="#C8342B", lw=1.0, label="tangentielle $F_x$")
    ax.plot(xs, fy, color="#0B4F9E", lw=1.0, label="normale $F_y$")
    ax.set_xlabel("course du cutter [mm]")
    ax.set_ylabel("force [MN/m]")
    ax.set_title("(d) forces de coupe", fontsize=10.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(2, 3, 5)
    ok = fy > 0.02 * max(fy.max(), 1e-12)
    ax.plot(xs[ok], (fx / np.maximum(fy, 1e-12))[ok], color="#6A3D9A", lw=1.0)
    ax.axhline(np.median((fx / np.maximum(fy, 1e-12))[ok]), color="0.4", ls="--",
               lw=1.2, label=f"mediane {np.median((fx/np.maximum(fy,1e-12))[ok]):.2f}")
    ax.set_xlabel("course du cutter [mm]")
    ax.set_ylabel(r"$F_x / F_y$")
    ax.set_title("(e) rapport tangentiel / normal", fontsize=10.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 5)

    ax = fig.add_subplot(2, 3, 6)
    ax.plot(xs, d["nBroken"], color=C_SHR, lw=1.6)
    ax.set_xlabel("course du cutter [mm]")
    ax.set_ylabel("joints rompus")
    ax.set_title("(f) fissuration cumulee", fontsize=10.5)
    ax.grid(alpha=0.3)

    fig.suptitle(f"Coupe 2D au cutter PDC — rake {a.rake:g} deg, passe "
                 f"{a.depth*1e3:g} mm, loi de Ye  ({os.path.basename(run)}, "
                 f"trame {k})", fontsize=12.5)
    fig.tight_layout()
    fig.savefig(out, dpi=155)
    print(f"course {xs.max():.2f} mm, Fx pic {fx.max():.3f} / median "
          f"{np.median(fx[ok]):.3f} MN/m, Fy pic {fy.max():.3f}, "
          f"{int(d['nBroken'][-1])} fissures")
    print("ecrit :", out)


if __name__ == "__main__":
    main()
