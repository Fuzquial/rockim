#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_joints_cuts.py — COUPES EXACTES du reseau de joints rompus (fdem3d).
#
#   python tools/fig_joints_cuts.py out_dir [--frame k] [--stem results/fig/x]
#                                  [--title "..."] [--lim 30] [--depth 25]
#
# Chaque facette rompue est INTERSECTEE avec le plan de coupe : la trace
# d'une fissure est alors un segment (pas un triangle projete), trace en
# trait epais, rouge = traction, orange = cisaillement. Huit coupes :
#   rangee 1 — verticales x-z a y = 0, 5, 10 mm et y-z a x = 0 ;
#   rangee 2 — horizontales x-y a z = -1, -3, -6, -10 mm sous la surface.
# En tete de chaque coupe : nombre de traces et longueur cumulee (mm).
# Reperes : insert R 8,51 mm ; Yang 9 m/s : cratere ~7 mm, radiales ~10 mm.
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench_impact", "tools"))
from imp_lib import (CX, CY, Z_SURF, broken, frame_times, frames_of,  # noqa
                     joints_frame, read_vtu)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})

ROUGE, ORANGE = "#b22222", "#e08a00"
R_INSERT, R_CRATERE, R_RADIAL = 8.51, 7.0, 10.0
LW = 2.0


def traces(P, k, v):
    """Segments d'intersection des triangles P (n,3,3) avec le plan x_k = v.
    Retourne (m,2,3) et l'index des triangles coupes."""
    d = P[:, :, k] - v
    s = np.sign(d)
    cross = (s.max(axis=1) > 0) & (s.min(axis=1) < 0)
    idx = np.where(cross)[0]
    segs, keep = [], []
    for i in idx:
        pts = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            da, db = d[i, a], d[i, b]
            if da * db < 0:
                t = da / (da - db)
                pts.append(P[i, a] + t * (P[i, b] - P[i, a]))
        if len(pts) == 2:
            segs.append(pts)
            keep.append(i)
    if not segs:
        return np.zeros((0, 2, 3)), np.zeros(0, int)
    return np.array(segs), np.array(keep)


def coupe(AX, P, mode, k, v, ax0, ax1, titre):
    segs, keep = traces(P, k, v)
    ltot = 0.0
    if len(segs):
        S = segs[:, :, (ax0, ax1)]
        ltot = float(np.linalg.norm(S[:, 1] - S[:, 0], axis=1).sum())
        mk = mode[keep]
        for mm, col, z in ((mk >= 1.5, ORANGE, 3), (mk < 1.5, ROUGE, 4)):
            if mm.any():
                AX.add_collection(LineCollection(S[mm], colors=col, linewidths=LW,
                                                 alpha=0.9, zorder=z,
                                                 capstyle="round"))
    AX.set_title("%s  —  %d traces, %.0f mm" % (titre, len(segs), ltot),
                 loc="left", fontsize=10)
    return len(segs), ltot


def reperes_xy(AX, lim):
    th = np.linspace(0, 2 * np.pi, 200)
    for r, col, ls in ((R_INSERT, "#333", "--"), (R_CRATERE, "#1f4e79", ":"),
                       (R_RADIAL, "#1f4e79", "-.")):
        AX.plot(r * np.cos(th), r * np.sin(th), ls, color=col, lw=0.9, zorder=6)
    AX.set_xlim(-lim, lim); AX.set_ylim(-lim, lim); AX.set_aspect("equal")
    AX.set_xlabel("x [mm]"); AX.set_ylabel("y [mm]")


def reperes_xz(AX, lim, depth, y0, xlabel="x [mm]"):
    AX.axhline(0, color="#333", lw=1.0, zorder=6)
    # section de l'insert (sphere R 8,51 centree en surface, trace a y = y0)
    if abs(y0) < R_INSERT:
        rr = np.sqrt(R_INSERT ** 2 - y0 ** 2)
        AX.axvline(-rr, color="#333", ls="--", lw=0.8, zorder=6)
        AX.axvline(rr, color="#333", ls="--", lw=0.8, zorder=6)
    for r, col, ls in ((R_CRATERE, "#1f4e79", ":"), (R_RADIAL, "#1f4e79", "-.")):
        if abs(y0) < r:
            rr = np.sqrt(r ** 2 - y0 ** 2)
            AX.axvline(-rr, color=col, ls=ls, lw=0.8, zorder=6)
            AX.axvline(rr, color=col, ls=ls, lw=0.8, zorder=6)
    AX.set_xlim(-lim, lim); AX.set_ylim(-depth, 4); AX.set_aspect("equal")
    AX.set_xlabel(xlabel); AX.set_ylabel("z sous la surface [mm]")


def main():
    global LW
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default=None)
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--title", default=None)
    ap.add_argument("--lim", type=float, default=30.0)
    ap.add_argument("--depth", type=float, default=25.0)
    ap.add_argument("--ys", default="0,5,10", help="coupes verticales y = ... (mm)")
    ap.add_argument("--zs", default="-1,-3,-6,-10", help="coupes horizontales z = ... (mm)")
    ap.add_argument("--lw", type=float, default=LW)
    a = ap.parse_args()
    LW = a.lw

    ks = frames_of(a.run)
    k = ks[-1] if a.frame < 0 else a.frame
    tk = frame_times(a.run).get(k, float("nan"))
    pts, con, f = read_vtu(joints_frame(a.run, k))
    c, n, mode, P = broken(pts, con, f)
    nb = len(c)
    # en mm, repere de l'impact (axe en 0,0 ; surface z = 0, profondeur < 0)
    Pm = (P - np.array([CX, CY, Z_SURF])) * 1e3 if nb else np.zeros((0, 3, 3))

    ys = [float(x) for x in a.ys.split(",")]
    zs = [float(x) for x in a.zs.split(",")]
    ncol = max(len(ys) + 1, len(zs))
    fig, AX = plt.subplots(2, ncol, figsize=(4.6 * ncol, 8.4),
                           gridspec_kw=dict(height_ratios=(1.0, 2.1), hspace=0.45))
    fig.suptitle((a.title or a.run) + "  —  coupes du reseau de joints rompus a "
                 "$t$ = %.0f $\\mu$s (%d facettes) ; rouge traction, orange cisaillement"
                 % (tk * 1e6, nb), fontsize=13)
    print("== %s  frame %d  t = %.1f us  %d facettes rompues" % (a.run, k, tk * 1e6, nb))

    # ---- rangee 1 : coupes verticales --------------------------------------
    for j, y0 in enumerate(ys):
        A = AX[0, j]
        nt, lt = coupe(A, Pm, mode, 1, y0, 0, 2, "coupe $y$ = %g mm" % y0)
        reperes_xz(A, a.lim, a.depth, y0)
        print("   y = %5.1f mm : %4d traces, %7.1f mm" % (y0, nt, lt))
    A = AX[0, len(ys)]
    nt, lt = coupe(A, Pm, mode, 0, 0.0, 1, 2, "coupe $x$ = 0 mm")
    reperes_xz(A, a.lim, a.depth, 0.0, xlabel="y [mm]")
    print("   x = %5.1f mm : %4d traces, %7.1f mm" % (0.0, nt, lt))
    for j in range(len(ys) + 1, ncol):
        AX[0, j].axis("off")

    # ---- rangee 2 : coupes horizontales ------------------------------------
    for j, z0 in enumerate(zs):
        A = AX[1, j]
        nt, lt = coupe(A, Pm, mode, 2, z0, 0, 1, "coupe $z$ = %g mm" % z0)
        reperes_xy(A, a.lim)
        print("   z = %5.1f mm : %4d traces, %7.1f mm" % (z0, nt, lt))
    for j in range(len(zs), ncol):
        AX[1, j].axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.95], h_pad=2.5)
    stem = a.stem or os.path.join("results", "fig", "coupes_" + os.path.basename(a.run))
    for ext in ("pdf", "png"):
        fig.savefig(stem + "." + ext, dpi=170)
    print("   ecrit : %s.pdf / .png" % stem)


if __name__ == "__main__":
    main()
