#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_couplage.py — L'EFFET DU COUPLAGE HYDRO, A MEME INSTANT.
#
#   python bench_abuaisha/tools/fig_couplage.py out_f7_aniso out_hfp_aniso
#
# Les deux runs partagent le maillage (189 487 elements), la roche, l'etat de
# contrainte, le schema et l'horloge (T = 3 ms, 60 trames). La seule difference
# est le CHEMIN DE LA PRESSION :
#
#   * out_f7_aniso  : `confiningPressure`, rampe imposee sur les faces
#     d'ORIGINE du forage. Les faces nees de la fissuration ne recoivent rien.
#   * out_hfp_aniso : `hydro = on`, pompe a debit constant. La pression est
#     portee par une CAVITE FLUIDE dont la connexite suit la fissuration.
#
# HONNETETE DE LA COMPARAISON. Le run couple est en cours : on ne compare donc
# pas les etats FINAUX (le temoin finit a 18 MPa de paroi, le couple n'en est
# pas la) mais la MEME TRAME, c'est-a-dire le meme instant physique. La
# pression de paroi du temoin y vaut 24 MPa x t / 4 ms ; elle est affichee.
# ---------------------------------------------------------------------------
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tunnel_edz"))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

plt.rcParams.update({"font.size": 10.5, "figure.dpi": 110})
CX, CY, RB = 4.0, 4.0, 0.05
TFIN, NFRAME = 3.0e-3, 60


def frames(run):
    """{index -> chemin} des trames de joints completes."""
    out = {}
    for f in sorted(glob.glob(os.path.join(run, "fdem_joints_0*.vtu"))):
        if complete(f):
            out[int(re.search(r"_(\d+)\.vtu$", f).group(1))] = f
    if not out:
        raise SystemExit(f"aucune trame complete dans {run}")
    return out


def cracks(path):
    # ncell=2 : les cellules d'un VTU de joints sont des SEGMENTS a 2 noeuds.
    P, C, D = read_vtu(path, ["damage"], ncell=2)
    d = np.asarray(D["damage"])
    return P, np.asarray(C), np.where(d >= 1.0)[0]


def panel(ax, path, titre, half):
    P, seg, br = cracks(path)
    port = 0.0
    if len(br):
        L = np.stack([P[seg[br, 0], :2], P[seg[br, 1], :2]], axis=1)
        ax.add_collection(LineCollection(L - np.array([CX, CY]),
                                         colors="#8c2318", linewidths=1.6,
                                         zorder=3))
        mid = 0.5 * (L[:, 0] + L[:, 1]) - np.array([CX, CY])
        port = np.hypot(mid[:, 0], mid[:, 1]).max()
    th = np.linspace(0, 2 * np.pi, 200)
    ax.fill(RB * np.cos(th), RB * np.sin(th), color="0.35", zorder=4)
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.set_xlabel("x  [m]")
    ax.set_title(f"{titre}\n{len(br)} joints rompus — portée max "
                 f"{port*100:.1f} cm  (paroi : 5,0 cm)", fontsize=10.5)
    return len(br), port


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sans")
    ap.add_argument("avec")
    ap.add_argument("--half", type=float, default=0.10)
    ap.add_argument("--out", default="f7_effet_couplage")
    a = ap.parse_args()
    root = os.path.join(HERE, "..", "..")
    rs = a.sans if os.path.isabs(a.sans) else os.path.join(root, a.sans)
    ra = a.avec if os.path.isabs(a.avec) else os.path.join(root, a.avec)

    fs, fa = frames(rs), frames(ra)
    k = max(fa)                                  # derniere trame du run couple
    if k not in fs:
        k = max(i for i in fs if i <= k)
    t = k * TFIN / NFRAME
    import math
    pw = 24.0 * 0.5 * (1.0 - math.cos(math.pi * min(t / 4.0e-3, 1.0)))

    fig = plt.figure(figsize=(15.4, 6.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.2], wspace=0.30)
    ax0, ax1, ax2 = (fig.add_subplot(gs[i]) for i in range(3))

    n0, p0 = panel(ax0, fs[k], "SANS couplage — pression imposée, signe JUSTE\n"
                   f"paroi à {pw:.1f} MPa à cet instant", a.half)
    ax0.set_ylabel("y  [m]")

    h = os.path.join(ra, "history.csv")
    rows = list(csv.DictReader(open(h))) if os.path.exists(h) else []
    tt = np.array([float(x["t"]) for x in rows]) * 1e3
    pp = np.array([float(x["hydroP"]) for x in rows]) / 1e6
    nb = np.array([float(x["nBroken"]) for x in rows])
    nw = np.array([float(x["hydroNWet"]) for x in rows])
    p_now = float(np.interp(t * 1e3, tt, pp)) if len(tt) else float("nan")

    n1, p1 = panel(ax1, fa[k], "AVEC couplage — SIGNE INVERSÉ (bug)\n"
                   f"cavité à {p_now:.1f} MPa à cet instant", a.half)
    ax1.legend(handles=[
        Line2D([0], [0], color="#8c2318", lw=2, label="joint rompu (D = 1)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="0.35",
               markersize=9, label="forage $\\varnothing$ 0,1 m")],
        loc="lower right", fontsize=8.5, framealpha=0.95)

    if len(tt):
        ax2.plot(tt, pp, "-", color="#1e4b8c", lw=1.9)
        j = int(np.argmax(pp))
        ax2.plot(tt[j], pp[j], "o", color="#c8342b", ms=7, zorder=6)
        ax2.annotate(f"pic {pp[j]:.2f} MPa", (tt[j], pp[j]),
                     textcoords="offset points", xytext=(-8, 12),
                     ha="right", color="#c8342b", fontsize=9.5)
        ax2.axhline(12.0, color="k", ls="--", lw=1.2)
        ax2.text(tt[-1], 11.4, "seuil analytique 12 MPa (éq. 10)",
                 ha="right", va="top", fontsize=9)
        ax2b = ax2.twinx()
        ax2b.plot(tt, nw, "-", color="#1b8a3a", lw=1.3, alpha=0.85)
        ax2b.plot(tt, nb, "--", color="#b07d1a", lw=1.3, alpha=0.9)
        ax2b.set_ylabel("nombre de joints", color="0.3")
        ax2b.legend(handles=[
            Line2D([0], [0], color="#1e4b8c", lw=2, label="pression de cavité"),
            Line2D([0], [0], color="#1b8a3a", lw=1.6, label="joints mouillés"),
            Line2D([0], [0], color="#b07d1a", lw=1.6, ls="--",
                   label="joints rompus")],
            loc="center left", fontsize=8.5, framealpha=0.95)
        ax2.set_xlabel("temps  [ms]")
        ax2.set_ylabel("pression de cavité  [MPa]", color="#1e4b8c")
        ax2.set_title("pression de puits et front mouillé\n"
                      "(leur figure 11) — run EN COURS", fontsize=10.5)
        ax2.grid(alpha=0.3)
        ax2.set_ylim(0, 13.5)

    fig.suptitle(f"Effet du couplage hydro-mécanique au même instant "
                 f"t = {t*1e3:.2f} ms  —  état anisotrope "
                 "$\\sigma'_H = -6{,}8$, $\\sigma'_h = -4{,}6$ MPa, "
                 "189 487 éléments", fontsize=12, y=0.99)
    fig.tight_layout(rect=[0, 0.115, 1, 0.94])
    fig.text(0.5, 0.028,
             "hydroForces() applique la pression SUIVANT la normale sortante au lieu de son OPPOSÉ : "
             "elle serre le forage au lieu de l'ouvrir.\nD'où une compression orthoradiale, et un écaillage "
             "de paroi aligné sur $\\sigma'_h$ (breakout) au lieu d'ailes de fracturation alignées sur "
             "$\\sigma'_H$. Rayon de paroi : $-15 \\rightarrow -21\\ \\mu$m quand $p$ monte (le trou se ferme), "
             "contre $-10 \\rightarrow -6\\ \\mu$m à gauche (il s'ouvre).",
             ha="center", fontsize=9, style="italic", color="0.3")
    out = os.path.join(HERE, "..", a.out + ".png")
    fig.savefig(out, dpi=155)
    fig.savefig(out.replace(".png", ".pdf"))
    print(f"ecrit : {out}   (trame {k}, t = {t*1e3:.3f} ms)")
    print(f"  sans couplage : {n0} joints, portee {p0*100:.1f} cm, paroi {pw:.2f} MPa")
    print(f"  avec couplage : {n1} joints, portee {p1*100:.1f} cm, cavite {p_now:.2f} MPa")


if __name__ == "__main__":
    main()
