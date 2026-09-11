#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_joints_only.py — RIEN QUE LES JOINTS ROMPUS d'un impact 3D (fdem3d).
#
#   python tools/fig_joints_only.py out_dir [--frame k] [--stem results/fig/x]
#                                  [--title "..."] [--lim 30] [--depth 25]
#
# Six vues, deux lectures du meme reseau :
#   rangee 1 — MODE de rupture : rouge = traction, jaune = cisaillement ;
#              (a) 3D, (b) vue de dessus, (c) coupe |y| < 5 mm ;
# Rendu par ARETES (facettes non remplies) par defaut, --faces pour le plein.
#   rangee 2 — ORIENTATION : faces sub-verticales (|n_z| <= 0,6) = FISSURES
#              (radiales, cone) en rouge sombre ombre par la profondeur,
#              faces sub-horizontales = zone BROYEE en rose pale ;
#              (d) 3D, (e) vue de dessus, (f) coupe.
# Cercles de reference (Yang 2026, 9 m/s) : cratere ~7 mm, radiales ~10 mm,
# et le rayon de l'insert 8,51 mm. Reutilise imp_lib (bench_impact/tools).
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench_impact", "tools"))
from imp_lib import (CX, CY, Z_SURF, broken, frame_times, frames_of,  # noqa
                     joints_frame, metrics, read_vtu)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})

ROUGE, JAUNE = "#b22222", "#e6a817"
ROUGE_SOMBRE = np.array([0.72, 0.13, 0.13])
ROSE = "#e8b7b7"
R_INSERT = 8.51          # mm, rayon de l'insert (deck 2D, meme insert)
R_CRATERE_YANG = 7.0     # mm, leur cratere a 9 m/s
R_RADIAL_YANG = 10.0     # mm, leurs radiales a 9 m/s


def cercles(AX, plein=True):
    th = np.linspace(0, 2 * np.pi, 200)
    for r, col, ls, lab in ((R_INSERT, "#333", "--", "insert R 8,5"),
                            (R_CRATERE_YANG, "#1f4e79", ":", "Yang cratere 7"),
                            (R_RADIAL_YANG, "#1f4e79", "-.", "Yang radiales 10")):
        AX.plot(r * np.cos(th), r * np.sin(th), ls, color=col, lw=0.9,
                zorder=6, label=lab)


def coupe_refs(AX, lim):
    AX.axhline(0, color="#333", lw=0.9, zorder=6)
    for r, col, ls in ((R_INSERT, "#333", "--"), (R_CRATERE_YANG, "#1f4e79", ":"),
                       (R_RADIAL_YANG, "#1f4e79", "-.")):
        AX.axvline(-r, color=col, ls=ls, lw=0.8, zorder=6)
        AX.axvline(r, color=col, ls=ls, lw=0.8, zorder=6)


EDGES = True             # --faces pour revenir au rendu plein


def poly2d(AX, P, ax0, ax1, c0, c1, col, al, z, lw=0.25):
    if len(P) == 0:
        return
    po = (P[:, :, (ax0, ax1)] - np.array([c0, c1])) * 1e3
    if EDGES:            # rien que les ARETES des facettes rompues
        AX.add_collection(PolyCollection(po, facecolors="none", edgecolors=col,
                                         linewidths=max(lw, 1.1), alpha=min(al + 0.1, 1.0),
                                         zorder=z))
    else:
        AX.add_collection(PolyCollection(po, facecolors=col, edgecolors=col,
                                         linewidths=lw, alpha=al, zorder=z))


def vue_mode_2d(AX, P, mode, ax0, ax1, c0, c1):
    poly2d(AX, P[mode < 1.5], ax0, ax1, c0, c1, ROUGE, 0.75, 3)
    poly2d(AX, P[mode >= 1.5], ax0, ax1, c0, c1, JAUNE, 0.85, 4)


def vue_orient_2d(AX, P, n, ax0, ax1, c0, c1):
    vert = np.abs(n[:, 2]) <= 0.6
    poly2d(AX, P[~vert], ax0, ax1, c0, c1, ROSE if not EDGES else "#c8c8c8", 0.35 if not EDGES else 0.6, 1)
    poly2d(AX, P[vert], ax0, ax1, c0, c1, ROUGE, 0.9, 3, lw=0.4)


def vue_3d(B, P, n, mode, lim, depth, par_mode):
    if len(P):
        P3 = (P - np.array([CX, CY, Z_SURF])) * 1e3
        if par_mode:
            for mm, col, al in ((mode < 1.5, ROUGE, 0.8), (mode >= 1.5, JAUNE, 0.9)):
                if mm.any():
                    if EDGES:
                        B.add_collection3d(Poly3DCollection(
                            P3[mm], facecolors=(1, 1, 1, 0), edgecolors=col,
                            linewidths=0.8, alpha=al))
                    else:
                        B.add_collection3d(Poly3DCollection(
                            P3[mm], facecolors=col, edgecolors="none", alpha=al))
        else:
            vert = np.abs(n[:, 2]) <= 0.6
            if (~vert).any():
                if EDGES:
                    B.add_collection3d(Poly3DCollection(
                        P3[~vert], facecolors=(1, 1, 1, 0), edgecolors="#b0b0b0",
                        linewidths=0.5, alpha=0.5))
                else:
                    B.add_collection3d(Poly3DCollection(
                        P3[~vert], facecolors="#c9c9c9", edgecolors="none", alpha=0.25))
            if vert.any():
                zc = P3[vert].mean(axis=1)[:, 2]
                lo, hi = float(zc.min()), float(zc.max())
                sh = 0.35 + 0.65 * (zc - lo) / max(hi - lo, 1e-9)
                cols = np.clip(np.outer(sh, ROUGE_SOMBRE), 0, 1)
                if EDGES:
                    B.add_collection3d(Poly3DCollection(
                        P3[vert], facecolors=(1, 1, 1, 0), edgecolors=cols,
                        linewidths=0.8, alpha=0.9))
                else:
                    B.add_collection3d(Poly3DCollection(
                        P3[vert], facecolors=cols, edgecolors="none", alpha=0.9))
    th = np.linspace(0, 2 * np.pi, 120)
    B.plot(R_INSERT * np.cos(th), R_INSERT * np.sin(th), 0 * th, "--",
           color="#333", lw=0.8)
    B.set_xlim(-lim, lim); B.set_ylim(-lim, lim); B.set_zlim(-depth, 4)
    B.view_init(elev=30, azim=-55)
    B.set_xlabel("x [mm]"); B.set_ylabel("y [mm]"); B.set_zlabel("z [mm]")
    B.set_box_aspect((1, 1, 0.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default=None)
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--title", default=None)
    ap.add_argument("--lim", type=float, default=30.0, help="demi-largeur, mm")
    ap.add_argument("--depth", type=float, default=25.0, help="profondeur, mm")
    ap.add_argument("--ysec", type=float, default=5.0, help="demi-epaisseur de coupe, mm")
    ap.add_argument("--faces", action="store_true", help="faces pleines (defaut : aretes)")
    a = ap.parse_args()
    global EDGES
    EDGES = not a.faces

    ks = frames_of(a.run)
    k = ks[-1] if a.frame < 0 else a.frame
    tk = frame_times(a.run).get(k, float("nan"))
    pts, con, f = read_vtu(joints_frame(a.run, k))
    c, n, mode, P = broken(pts, con, f)
    m = metrics(c)
    nb = len(c)
    nT = int((mode < 1.5).sum()) if nb else 0
    vert = (np.abs(n[:, 2]) <= 0.6) if nb else np.zeros(0, bool)
    nV = int(vert.sum())
    r = np.hypot(c[:, 0] - CX, c[:, 1] - CY) * 1e3 if nb else np.zeros(0)
    zr = (Z_SURF - c[:, 2]) * 1e3 if nb else np.zeros(0)
    rV = float(r[vert].max()) if nV else 0.0
    zV = float(zr[vert].max()) if nV else 0.0

    print("== %s  frame %d  t = %.1f us" % (a.run, k, tk * 1e6))
    print("   joints rompus        : %d  (traction %d = %.0f %%, cisaillement %d)"
          % (nb, nT, 100.0 * nT / max(nb, 1), nb - nT))
    print("   faces sub-verticales : %d = %.0f %%  (fissures) ; sub-horizontales %d (broye)"
          % (nV, 100.0 * nV / max(nb, 1), nb - nV))
    print("   extension radiale    : toutes %.1f mm ; fissures %.1f mm   [Yang ~10]"
          % (m["radial"] * 1e3, rV))
    print("   rayon de cratere     : %.1f mm (faces < 3 mm sous la surface) [Yang ~7]"
          % (m["crater"] * 1e3))
    print("   profondeur           : toutes %.1f mm ; fissures %.1f mm" % (m["depth"] * 1e3, zV))

    lim, depth = a.lim, a.depth
    fig = plt.figure(figsize=(17.0, 10.6))
    fig.suptitle((a.title or a.run) + "  —  joints rompus a $t$ = %.0f $\\mu$s : "
                 "%d faces, %.0f %% traction, %.0f %% sub-verticales"
                 % (tk * 1e6, nb, 100.0 * nT / max(nb, 1), 100.0 * nV / max(nb, 1)),
                 fontsize=13)

    # ---- rangee 1 : mode ---------------------------------------------------
    A = fig.add_subplot(2, 3, 1, projection="3d")
    vue_3d(A, P, n, mode, lim, depth, par_mode=True)
    A.set_title("(a)  3D — rouge traction, jaune cisaillement", loc="left", fontsize=11)

    B = fig.add_subplot(2, 3, 2)
    if nb:
        vue_mode_2d(B, P, mode, 0, 1, CX, CY)
    cercles(B)
    B.set_xlim(-lim, lim); B.set_ylim(-lim, lim); B.set_aspect("equal")
    B.set_xlabel("x [mm]"); B.set_ylabel("y [mm]")
    B.set_title("(b)  Vue de dessus, par mode", loc="left", fontsize=11)
    B.legend(frameon=False, fontsize=8, loc="upper right")

    C = fig.add_subplot(2, 3, 3)
    if nb:
        s5 = np.abs(c[:, 1] - CY) < a.ysec * 1e-3
        vue_mode_2d(C, P[s5], mode[s5], 0, 2, CX, Z_SURF)
    coupe_refs(C, lim)
    C.set_xlim(-lim, lim); C.set_ylim(-depth, 5); C.set_aspect("equal")
    C.set_xlabel("x [mm]"); C.set_ylabel("z sous la surface [mm]")
    C.set_title("(c)  Coupe $|y| < %g$ mm, par mode" % a.ysec, loc="left", fontsize=11)

    # ---- rangee 2 : orientation -------------------------------------------
    D = fig.add_subplot(2, 3, 4, projection="3d")
    vue_3d(D, P, n, mode, lim, depth, par_mode=False)
    D.set_title("(d)  3D — fissures (sub-verticales) ombrees, broye en gris",
                loc="left", fontsize=11)

    E = fig.add_subplot(2, 3, 5)
    if nb:
        vue_orient_2d(E, P, n, 0, 1, CX, CY)
    cercles(E)
    E.set_xlim(-lim, lim); E.set_ylim(-lim, lim); E.set_aspect("equal")
    E.set_xlabel("x [mm]"); E.set_ylabel("y [mm]")
    E.set_title("(e)  Vue de dessus — fissures rouge, broye rose", loc="left", fontsize=11)

    F = fig.add_subplot(2, 3, 6)
    if nb:
        s5 = np.abs(c[:, 1] - CY) < a.ysec * 1e-3
        vue_orient_2d(F, P[s5], n[s5], 0, 2, CX, Z_SURF)
    coupe_refs(F, lim)
    F.set_xlim(-lim, lim); F.set_ylim(-depth, 5); F.set_aspect("equal")
    F.set_xlabel("x [mm]"); F.set_ylabel("z sous la surface [mm]")
    F.set_title("(f)  Coupe $|y| < %g$ mm — fissures / broye" % a.ysec,
                loc="left", fontsize=11)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    stem = a.stem or os.path.join("results", "fig", "joints_" + os.path.basename(a.run))
    for ext in ("pdf", "png"):
        fig.savefig(stem + "." + ext, dpi=170)
    print("   ecrit : %s.pdf / .png" % stem)


if __name__ == "__main__":
    main()
