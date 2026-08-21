#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_postpic.py — LA PHASE POST-PIC du benchmark B2 (AbuAisha et al. 2017).
#
#   python bench_abuaisha/tools/fig_postpic.py out_hfs_aniso --sh 6.8 --sv 4.6
#
# Ce que la figure principale ne montre pas : ce qui se passe APRES la rupture.
# C'est la partie que le benchmark declarait hors de portee tant que la
# pression ne suivait pas la fissure — leur figure 11 y lit un plateau de
# propagation vers 5,5 MPa, la contrainte lointaine effective.
#
#   (a) la chute de pression, avec le plateau vise en reference
#   (b) la propagation : joints rompus, faces mouillees, volume de la cavite
#   (c) le facies a la derniere trame
#
# Lecture bon marche : (a) et (b) ne lisent que history.csv, (c) une seule
# trame de joints.
# ---------------------------------------------------------------------------
import argparse
import csv
import glob
import io
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
from plot_tunnel_fields import complete  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 10.5,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "savefig.dpi": 200,
})

CX, CY, RB = 4.0, 4.0, 0.05
FT = 5.0e6


def vg(x, n=2):
    return (("%." + str(n) + "f") % x).replace(".", ",")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--sh", type=float, default=6.8)
    ap.add_argument("--sv", type=float, default=4.6)
    ap.add_argument("--plateau", type=float, default=5.5,
                    help="plateau de propagation lu dans leur figure 11")
    ap.add_argument("--rfine", type=float, default=0.46,
                    help="rayon de la zone raffinee [m], mesure sur le maillage")
    ap.add_argument("--stem", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(ROOT, a.run)

    h = list(csv.DictReader(open(os.path.join(run, "history.csv"))))
    t = np.array([float(r["t"]) for r in h]) * 1e3
    p = np.array([float(r["hydroP"]) for r in h]) / 1e6
    V = np.array([float(r["hydroVol"]) for r in h])
    nb = np.array([int(r["nBroken"]) for r in h])
    nw = np.array([int(float(r["hydroNWet"])) for r in h])
    V0, ip = V[0], int(np.argmax(p))
    t0 = t[ip] - 0.15                       # on cadre juste avant le pic

    fig = plt.figure(figsize=(14.8, 4.9))
    gs = fig.add_gridspec(1, 3, wspace=0.44)

    # --- (a) la chute de pression ------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    m = t >= t0
    ax.plot(t[m], p[m], lw=1.9, color="#1f4e79")
    ax.plot(t[ip], p[ip], "o", ms=6.5, mfc="none", mec="#a11122", mew=1.5)
    ax.annotate("pic  " + vg(p[ip]) + " MPa", (t[ip], p[ip]),
                textcoords="offset points", xytext=(8, 4), fontsize=10,
                color="#a11122")
    ax.axhline(a.plateau, color="#0a7", ls=(0, (4, 3)), lw=1.1)
    ax.annotate("plateau de propagation visé  " + vg(a.plateau, 1) + " MPa",
                (t[m][-1], a.plateau), textcoords="offset points",
                xytext=(-6, 6), ha="right", fontsize=9, color="#076")
    ax.axhline(a.sv, color="0.55", ls=":", lw=1.0)
    ax.annotate(r"$\sigma^{\prime}_h$ = " + vg(a.sv, 1) + " MPa",
                (t[m][0], a.sv), textcoords="offset points", xytext=(4, 4),
                fontsize=9, color="0.4")
    ax.set_xlabel("temps  [ms]")
    ax.set_ylabel("pression de puits  [MPa]")
    ax.set_title("(a)  Chute post-pic", loc="left")
    ax.grid(alpha=0.22, lw=0.6)

    # --- (b) la propagation -------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(t[m], nb[m], lw=1.9, color="#7a2048", label="joints rompus")
    ax.plot(t[m], nw[m] - nw[0], lw=1.5, ls="--", color="#c07",
            label="faces mouillées créées")
    ax.set_xlabel("temps  [ms]")
    ax.set_ylabel("nombre")
    ax.set_title("(b)  Propagation", loc="left")
    ax.grid(alpha=0.22, lw=0.6)
    ax.legend(fontsize=9, loc="upper left")
    ax2 = ax.twinx()
    ax2.plot(t[m], (V[m] / V0 - 1.0) * 100.0, lw=1.7, color="#12796f")
    ax2.set_ylabel(r"$V/V_0 - 1$   [%]", color="#12796f")
    ax2.tick_params(axis="y", colors="#12796f")

    # --- (c) le facies ------------------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    jf = [f for f in sorted(glob.glob(
        os.path.join(run, "fdem_joints_[0-9]*.vtu"))) if complete(f)]
    txt = io.open(jf[-1], errors="ignore").read()

    def grab(nm):
        mm = re.search('Name="' + nm + '"[^>]*>(.*?)</DataArray>', txt, re.S)
        return np.fromstring(mm.group(1), sep=" ")

    JP = np.fromstring(re.search(
        r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", txt, re.S).group(1),
        sep=" ").reshape(-1, 3)[:, :2]
    JC = grab("connectivity").astype(int).reshape(-1, 2)
    JD = grab("damage")
    msk = JD > 0.01
    seg = (JP[JC[msk]] - [CX, CY]) * 1e3
    ax.add_patch(plt.Circle((0.0, 0.0), RB * 1e3, fc="0.94", ec="0.55",
                            lw=0.7, zorder=1))
    # LA LIMITE DE LA ZONE RAFFINEE, mesuree sur le maillage : au-dela, les
    # elements passent de 3 mm a 0,3 m et la fissure ne peut plus se propager
    # de facon signifiante. Une aile qui touche ce cercle est bornee par le
    # MAILLAGE, plus par la physique.
    ax.add_patch(plt.Circle((0.0, 0.0), a.rfine * 1e3, fc="none", ec="#c60",
                            ls=(0, (5, 4)), lw=1.1, zorder=3))
    ax.annotate("limite de la zone raffinée", (0.0, a.rfine * 1e3),
                textcoords="offset points", xytext=(0, 5), ha="center",
                fontsize=8.5, color="#c60")
    lc = LineCollection(seg, array=JD[msk], cmap="Greys", lw=1.6,
                        clim=(-0.25, 1.0), zorder=2)
    ax.add_collection(lc)
    cb = fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"endommagement du joint   $D$")
    cb.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    cb.outline.set_linewidth(0.6)
    z = max(0.08, 1.14 * a.rfine, 1.10 * float(np.abs(seg).max()) / 1e3)
    ax.set_xlim(-z * 1e3, z * 1e3)
    ax.set_ylim(-z * 1e3, z * 1e3)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$   [mm]")
    ax.set_ylabel(r"$y$   [mm]")
    ax.set_title("(c)  Faciès  ($t$ = " + vg(t[-1]) + " ms)", loc="left")

    # demi-longueurs des deux ailes
    mid = 0.5 * (seg[:, 0, :] + seg[:, 1, :])
    lr = mid[mid[:, 0] > 0, 0].max() if (mid[:, 0] > 0).any() else 0.0
    lg = -mid[mid[:, 0] < 0, 0].min() if (mid[:, 0] < 0).any() else 0.0
    ax.text(0.03, 0.03, "aile droite  " + vg(lr, 0) + " mm\naile gauche  "
            + vg(lg, 0) + " mm", transform=ax.transAxes, fontsize=9.5,
            va="bottom", bbox=dict(fc="w", ec="0.75", alpha=0.9, pad=3))

    reg = "anisotrope" if abs(a.sh - a.sv) > 1e-9 else "isotrope"
    fig.suptitle("Phase post-pic — état de contrainte " + reg
                 + r"   ($\sigma^{\prime}_H$ = " + vg(a.sh, 1)
                 + r" MPa,  $\sigma^{\prime}_h$ = " + vg(a.sv, 1) + " MPa)",
                 fontsize=12.5, y=0.99)

    stem = a.stem or os.path.join(
        ROOT, "bench_abuaisha",
        "fig_postpic_" + os.path.basename(run).replace("out_hfs_", ""))
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    fig.savefig(stem + ".png", bbox_inches="tight")
    print("écrit :", stem + ".pdf  et  .png")
    print("  pic %s MPa -> %s MPa à t = %s ms  (%.0f %% de chute)"
          % (vg(p[ip]), vg(p[-1]), vg(t[-1]), (1 - p[-1] / p[ip]) * 100))
    print("  ailes : droite %s mm, gauche %s mm | %d joints rompus, %d faces"
          % (vg(lr, 0), vg(lg, 0), nb[-1], nw[-1]))


if __name__ == "__main__":
    main()
