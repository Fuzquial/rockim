#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_impact.py — LA planche de l'impact a insert unique (spec 005, WP5).
#
#   python bench_impact/tools/fig_impact.py out_imp_stanne --stem bench_impact/fig_stanne
#
# RENDU PAR FACES, comme les figures de l'article : une fissure radiale est
# un plan vertical dont la projection en vue de dessus est un TRAIT — c'est
# ce trait qui dessine l'etoile, la ou des centroides ne font qu'un nuage.
#   (a) vitesses du bit et du piston (leur fig. 8) ;
#   (b) le reseau de fissures en 3D, cratere compris (leurs fig. 11 et 16) ;
#   (c) fissures vues de DESSUS, faces projetees, rouge = traction,
#       jaune = cisaillement (leurs fig. 7 et 14, rangee du haut) ;
#   (d) coupe verticale |y| < 5 mm (leur fig. 14, rangee du bas).
# Les 7 metriques de leur Table 3 sont imprimees en tete.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import (CX, CY, Z_SURF, broken, frame_times, frames_of, history,
                     joints_frame, metrics, read_vtu)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})

ROUGE, JAUNE = "#b22222", "#e6a817"


def faces_2d(AX, P, mode, n, ax0, ax1, c0, c1):
    """Faces projetees dans le plan (ax0, ax1), en mm.

    La lisibilite des figures de l'article vient de la GEOMETRIE : une
    fissure radiale/mediane est un plan SUB-VERTICAL, la zone broyee est
    faite de faces de toutes orientations. On rend donc :
      - les faces sub-horizontales (|n_z| > 0,6) en rose pale : le broye ;
      - les faces sub-verticales en rouge/jaune francs : les FISSURES,
        qui se projettent en traits et dessinent l'etoile.
    """
    vert = np.abs(n[:, 2]) <= 0.6
    lots = ((~vert, "#e8b7b7", 0.30, 1),
            (vert & (mode < 1.5), ROUGE, 0.85, 3),
            (vert & (mode >= 1.5), JAUNE, 0.9, 4))
    for mm, col, al, z in lots:
        if not mm.any():
            continue
        po = (P[mm][:, :, (ax0, ax1)] - np.array([c0, c1])) * 1e3
        AX.add_collection(PolyCollection(
            po, facecolors=col, edgecolors=col, linewidths=0.2,
            alpha=al, zorder=z))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_impact")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--title", default=None,
                    help="titre de la planche (defaut : le cas St Anne 10,66 m/s)")
    a = ap.parse_args()

    h = history(a.run)
    t = h["t"] * 1e6
    ks = frames_of(a.run)
    k = ks[-1] if a.frame < 0 else a.frame
    tk = frame_times(a.run).get(k, h["t"][-1])
    pts, con, f = read_vtu(joints_frame(a.run, k))
    c, n, mode, P = broken(pts, con, f)
    m = metrics(c)

    vz = h["vz_bit"]
    vInd = -vz.min()
    after = np.argmin(vz)
    vReb = vz[after:].max()
    zi = h["z_insert"]
    depth = (zi[0] - zi.min()) * 1e3
    szz = np.abs(h["szz_bit"]).max() / 1e6 if "szz_bit" in h else float("nan")

    print("  == 7 criteres (leur Table 3, 10,66 m/s, fourchettes maillages) ==")
    print("  contrainte max au bit  : %7.1f MPa   (leur fig. 9a : ~200-260)" % szz)
    print("  vitesse d indentation  : %7.2f m/s   (9,40 - 9,85)" % vInd)
    print("  vitesse de rebond      : %7.2f m/s   (6,87 - 7,10)" % vReb)
    print("  profondeur d indentation: %6.2f mm    (~1,53)" % depth)
    print("  fissure radiale max    : %7.1f mm    (20,2 - 24,5)" % (m["radial"] * 1e3))
    print("  rayon de cratere       : %7.1f mm    (10,0 - 12,1)" % (m["crater"] * 1e3))
    print("  profondeur fissuree    : %7.1f mm" % (m["depth"] * 1e3))
    print("  joints rompus          : %7d" % m["n"])

    fig = plt.figure(figsize=(12.4, 10.2))
    fig.suptitle(a.title or ("Impact à insert unique — calcaire St Anne, piston à "
                 "10,66 m/s  (schéma adaptatif, DIF Yang fig. 2)"), fontsize=13)

    A = fig.add_subplot(2, 2, 1)
    if "vz_piston" in h:                # montage complet ; absent en allege
        A.plot(t, h["vz_piston"], color="#888", lw=1.2, label="piston")
    A.plot(t, vz, color="#1f4e79", lw=1.6, label="bit")
    A.axhline(0, color="k", lw=0.5)
    A.annotate("indentation %.2f m/s" % vInd, (t[after], vz.min()),
               textcoords="offset points", xytext=(8, -2), fontsize=9)
    A.set_xlabel(r"temps [$\mu$s]")
    A.set_ylabel(r"$v_z$  [m/s]")
    A.set_title("(a)  Vitesses des corps", loc="left", fontsize=11)
    A.legend(frameon=False, fontsize=9)

    # (b) le reseau en 3D — leurs fig. 11/16 : faces rompues, ombrees par la
    # profondeur, vue plongeante
    B = fig.add_subplot(2, 2, 2, projection="3d")
    if len(c):
        r45 = np.hypot(c[:, 0] - CX, c[:, 1] - CY) < 0.045
        P3 = (P[r45] - np.array([CX, CY, Z_SURF])) * 1e3
        n3 = n[r45]
        vert = np.abs(n3[:, 2]) <= 0.6
        # le broye en gris clair (leur cratere), les fissures en rouge ombre
        if (~vert).any():
            B.add_collection3d(Poly3DCollection(
                P3[~vert], facecolors="#c9c9c9", edgecolors="none", alpha=0.25))
        if vert.any():
            zc = P3[vert].mean(axis=1)[:, 2]
            lo, hi = float(zc.min()), float(zc.max())
            sh = 0.35 + 0.65 * (zc - lo) / max(hi - lo, 1e-9)
            cols = np.outer(sh, np.array([0.72, 0.13, 0.13]))
            B.add_collection3d(Poly3DCollection(
                P3[vert], facecolors=np.clip(cols, 0, 1),
                edgecolors="none", alpha=0.9))
        B.set_xlim(-40, 40); B.set_ylim(-40, 40); B.set_zlim(-30, 4)
    B.view_init(elev=32, azim=-55)
    B.set_xlabel("x [mm]"); B.set_ylabel("y [mm]"); B.set_zlabel("z [mm]")
    B.set_title("(b)  Réseau de fissures, 3D", loc="left", fontsize=11)
    B.set_box_aspect((1, 1, 0.45))

    C = fig.add_subplot(2, 2, 3)
    if len(c):
        faces_2d(C, P, mode, n, 0, 1, CX, CY)
    th = np.linspace(0, 2 * np.pi, 100)
    C.plot(m["crater"] * 1e3 * np.cos(th), m["crater"] * 1e3 * np.sin(th),
           ":", color="#666", lw=0.8, zorder=4)
    C.set_xlim(-45, 45); C.set_ylim(-45, 45)
    C.set_aspect("equal")
    C.set_xlabel("x [mm]"); C.set_ylabel("y [mm]")
    C.set_title("(c)  Fissures, vue de dessus  ($t$ = %.0f $\\mu$s)"
                % (tk * 1e6), loc="left", fontsize=11)

    D = fig.add_subplot(2, 2, 4)
    if len(c):
        s5 = np.abs(c[:, 1] - CY) < 0.005
        faces_2d(D, P[s5], mode[s5], n[s5], 0, 2, CX, Z_SURF)
    D.axhline(0, color="#333", lw=0.8)
    D.set_xlim(-45, 45); D.set_ylim(-30, 6)
    D.set_aspect("equal")
    D.set_xlabel("x [mm]"); D.set_ylabel("z sous la surface [mm]")
    D.set_title("(d)  Coupe verticale $|y| < 5$ mm", loc="left", fontsize=11)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=170)
    print("écrit : %s.pdf  et  .png" % a.stem)


if __name__ == "__main__":
    main()
