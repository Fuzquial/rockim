#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_impact.py — LA planche de l'impact a insert unique (spec 005, WP5).
#
#   python bench_impact/tools/fig_impact.py out_imp_stanne --stem bench_impact/fig_stanne
#
# Quatre panneaux, au format des figures de l'article :
#   (a) vitesses du bit et du piston — indentation, rebond (leur fig. 8) ;
#   (b) la jauge sigma_zz a mi-bit — l'onde de frappe (leur fig. 9a) ;
#   (c) fissures vues de DESSUS, rouge = traction, jaune = cisaillement
#       (leurs fig. 7 et 14, rangee du haut) ;
#   (d) coupe verticale |y| < 5 mm (leur fig. 14, rangee du bas).
# Les 7 metriques mesurees sont imprimees face aux fourchettes de leur
# Table 3 (maillages non structures, 10,66 m/s).
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import (Z_SURF, broken, frame_times, frames_of, history,
                     joints_frame, metrics, read_vtu)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})

ROUGE, JAUNE = "#b22222", "#e6a817"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_impact")
    ap.add_argument("--frame", type=int, default=-1)
    a = ap.parse_args()

    h = history(a.run)
    t = h["t"] * 1e6                        # us
    ks = frames_of(a.run)
    k = ks[-1] if a.frame < 0 else a.frame
    tk = frame_times(a.run).get(k, h["t"][-1])
    pts, con, f = read_vtu(joints_frame(a.run, k))
    c, n, mode, _ = broken(pts, con, f)
    m = metrics(c)

    # metriques d energie (leur fig. 8)
    vz = h["vz_bit"]
    vInd = -vz.min()
    after = np.argmin(vz)
    vReb = vz[after:].max()
    zi = h["z_insert"]
    depth = (zi[0] - zi.min()) * 1e3        # mm
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

    fig, ax = plt.subplots(2, 2, figsize=(11.5, 9.6))
    fig.suptitle("Impact à insert unique — calcaire St Anne, piston à "
                 "10,66 m/s  (schéma adaptatif, DIF Yang fig. 2)",
                 fontsize=13)

    A = ax[0, 0]
    A.plot(t, h["vz_piston"], color="#888", lw=1.2, label="piston")
    A.plot(t, vz, color="#1f4e79", lw=1.6, label="bit")
    A.axhline(0, color="k", lw=0.5)
    A.annotate("indentation %.2f m/s" % vInd, (t[after], vz.min()),
               textcoords="offset points", xytext=(8, -2), fontsize=9)
    A.set_xlabel(r"temps [$\mu$s]")
    A.set_ylabel(r"$v_z$  [m/s]")
    A.set_title("(a)  Vitesses des corps", loc="left", fontsize=11)
    A.legend(frameon=False, fontsize=9)

    B = ax[0, 1]
    if "szz_bit" in h:
        B.plot(t, -h["szz_bit"] / 1e6, color="#1f4e79", lw=1.4)
    B.set_xlabel(r"temps [$\mu$s]")
    B.set_ylabel(r"$-\sigma_{zz}$ à mi-bit  [MPa]")
    B.set_title("(b)  Jauge du bit (compression $>0$)", loc="left", fontsize=11)

    def scat(AX, x, y):
        s = mode >= 1.5                     # II = cisaillement
        AX.scatter(x[~s] * 1e3, y[~s] * 1e3, s=2.0, c=ROUGE, lw=0,
                   label="traction")
        AX.scatter(x[s] * 1e3, y[s] * 1e3, s=2.0, c=JAUNE, lw=0,
                   label="cisaillement")
        AX.set_aspect("equal")

    C = ax[1, 0]
    if len(c):
        scat(C, c[:, 0], c[:, 1])
    th = np.linspace(0, 2 * np.pi, 100)
    C.plot(m["crater"] * 1e3 * np.cos(th), m["crater"] * 1e3 * np.sin(th),
           ":", color="#666", lw=0.8)
    C.set_xlabel("x [mm]")
    C.set_ylabel("y [mm]")
    C.set_title("(c)  Fissures, vue de dessus  ($t$ = %.0f $\\mu$s)" % (tk * 1e6),
                loc="left", fontsize=11)
    C.legend(frameon=False, fontsize=9, markerscale=4)

    D = ax[1, 1]
    if len(c):
        s5 = np.abs(c[:, 1]) < 0.005
        cc, mm5 = c[s5], mode[s5]
        s = mm5 >= 1.5
        D.scatter(cc[~s][:, 0] * 1e3, (cc[~s][:, 2] - Z_SURF) * 1e3, s=2.5,
                  c=ROUGE, lw=0)
        D.scatter(cc[s][:, 0] * 1e3, (cc[s][:, 2] - Z_SURF) * 1e3, s=2.5,
                  c=JAUNE, lw=0)
    D.axhline(0, color="#333", lw=0.8)
    D.set_xlabel("x [mm]")
    D.set_ylabel("z sous la surface [mm]")
    D.set_aspect("equal")
    D.set_title("(d)  Coupe verticale $|y| < 5$ mm", loc="left", fontsize=11)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=170)
    print("écrit : %s.pdf  et  .png" % a.stem)


if __name__ == "__main__":
    main()
