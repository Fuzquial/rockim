#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_b2.py — PLANCHE LIVRABLE du benchmark B2 (AbuAisha et al. 2017).
#
#   python bench_abuaisha/tools/fig_b2.py out_hfs_aniso --sh 6.8 --sv 4.6
#   python bench_abuaisha/tools/fig_b2.py out_hfs_iso   --sh 4.6 --sv 4.6
#
# Sort un PDF vectoriel (les champs sont rasterises pour la taille du fichier,
# le texte et les axes restent vectoriels) et un PNG.
#
# Destinee au PARTAGE : pas de jargon interne, aucun renvoi aux numeros
# d'equation de l'article. Les references vont dans la legende du document.
# Relancer la commande apres la fin du run met la planche a jour.
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
from matplotlib.collections import PolyCollection, LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

# Computer Modern n'est pas installe sur cette machine : le rendu mathematique
# reste CM (fontset matplotlib), le texte courant prend STIXGeneral, le serif
# le plus proche disponible qui porte les accents.
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
    """Nombre a la francaise : virgule decimale."""
    return (("%." + str(n) + "f") % x).replace(".", ",")


def repere(ax):
    """Les deux fleches qui rappellent l'orientation des contraintes."""
    ax.annotate("", xy=(0.97, 0.07), xytext=(0.79, 0.07),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="<->", lw=0.9, color="0.35"))
    ax.text(0.88, 0.10, r"$\sigma^{\prime}_H$", transform=ax.transAxes,
            ha="center", fontsize=10, color="0.35")
    ax.annotate("", xy=(0.07, 0.94), xytext=(0.07, 0.76),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="<->", lw=0.9, color="0.35"))
    ax.text(0.125, 0.85, r"$\sigma^{\prime}_h$", transform=ax.transAxes,
            va="center", fontsize=10, color="0.35")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--sh", type=float, default=6.8)
    ap.add_argument("--sv", type=float, default=4.6)
    ap.add_argument("--zoom", type=float, default=0.30)
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
    pcib = -a.sh + 3.0 * a.sv + FT / 1e6

    fig = plt.figure(figsize=(10.6, 8.6))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.24)

    # --- (a) pression de puits ---------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(t, p, lw=1.8, color="#1f4e79")
    ax.axhline(pcib, color="0.45", ls=(0, (4, 3)), lw=0.9)
    ax.text(t[1], pcib, "  seuil analytique  " + vg(pcib, 1) + " MPa",
            color="0.35", fontsize=9, va="bottom")
    ax.plot(t[ip], p[ip], "o", ms=6.5, mfc="none", mec="#a11122", mew=1.5)
    ax.annotate(vg(p[ip]) + " MPa", (t[ip], p[ip]),
                textcoords="offset points", xytext=(-12, 9), ha="right",
                fontsize=10, color="#a11122")
    ax.set_xlabel("temps  [ms]")
    ax.set_ylabel("pression de puits  [MPa]")
    ax.set_title("(a)  Pression de puits", loc="left")
    ax.grid(alpha=0.22, lw=0.6)

    # --- (b) volume de la cavite -------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(t, (V / V0 - 1.0) * 100.0, lw=1.8, color="#12796f")
    ax.axhline(0.0, color="0.45", lw=0.8)
    ax.axvline(t[ip], color="#a11122", ls=(0, (4, 3)), lw=0.9)
    ax.set_xlabel("temps  [ms]")
    ax.set_ylabel(r"$V/V_0 - 1$   [%]")
    ax.set_title("(b)  Volume de la cavité", loc="left")
    ax.grid(alpha=0.22, lw=0.6)

    # --- champs -------------------------------------------------------------
    # (c) est pris A L'INSTANT DU PIC : c'est ce champ-la qui dicte ou la
    # fissure s'amorce. Le prendre a la fin le noierait sous la compression
    # post-rupture et les ondes reflechies.
    fs = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    frw = list(csv.DictReader(open(os.path.join(run, "frames.csv"))))
    tfr = np.array([float(r["t"]) for r in frw]) * 1e3

    def kof(f):
        return int(os.path.basename(f).split("_")[1].split(".")[0])

    kpic = min((kof(f) for f in fs), key=lambda q: abs(tfr[q] - t[ip]))
    fpic = [f for f in fs if kof(f) == kpic][0]
    P, C, F = read_vtu(fpic, ["sigmaXX", "sigmaYY", "sigmaXY"])
    tf = tfr[kof(fs[-1])]
    sxx, syy, sxy = F["sigmaXX"], F["sigmaYY"], F["sigmaXY"]
    s1 = 0.5 * (sxx + syy) + np.hypot(0.5 * (sxx - syy), sxy)
    cen = P[C].mean(axis=1)
    zc = min(a.zoom, 0.15)
    sel = np.where((np.abs(cen[:, 0] - CX) < zc) &
                   (np.abs(cen[:, 1] - CY) < zc))[0]

    # --- (c) contrainte principale majeure, AU PIC -------------------------
    axc = fig.add_subplot(gs[1, 0])
    val = s1[sel] / 1e6
    v = float(np.percentile(np.abs(val), 98.0))
    pc = PolyCollection([(P[C[i]] - [CX, CY]) * 1e3 for i in sel], array=val,
                        cmap="RdBu_r", lw=0, clim=(-v, v), rasterized=True)
    axc.add_collection(pc)
    th = np.linspace(0.0, 2.0 * np.pi, 240)
    axc.plot(RB * 1e3 * np.cos(th), RB * 1e3 * np.sin(th), "k-", lw=0.7)
    cb = fig.colorbar(pc, ax=axc, fraction=0.046, pad=0.03)
    cb.set_label(r"$\sigma_1$   [MPa]")
    cb.outline.set_linewidth(0.6)
    axc.set_title("(c)  Contrainte principale majeure, au pic  "
                  + "($t$ = " + vg(tfr[kpic]) + " ms)", loc="left")

    # --- (d) fissuration ----------------------------------------------------
    axd = fig.add_subplot(gs[1, 1])
    jf = [f for f in sorted(glob.glob(
        os.path.join(run, "fdem_joints_[0-9]*.vtu"))) if complete(f)]
    txt = io.open(jf[-1], errors="ignore").read()

    def jarr(nm):
        m = re.search('Name="' + nm + '"[^>]*>(.*?)</DataArray>', txt, re.S)
        return np.fromstring(m.group(1), sep=" ")

    JP = np.fromstring(re.search(
        r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", txt, re.S).group(1),
        sep=" ").reshape(-1, 3)[:, :2]
    JC = jarr("connectivity").astype(int).reshape(-1, 2)
    JD = jarr("damage")
    m = JD > 0.01
    axd.add_patch(plt.Circle((0.0, 0.0), RB * 1e3, fc="0.94", ec="0.55",
                             lw=0.7, zorder=1))
    lc = LineCollection((JP[JC[m]] - [CX, CY]) * 1e3, array=JD[m],
                        cmap="Greys", lw=1.7, clim=(-0.25, 1.0), zorder=2)
    axd.add_collection(lc)
    cb = fig.colorbar(lc, ax=axd, fraction=0.046, pad=0.03)
    cb.set_label(r"endommagement du joint   $D$")
    cb.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    cb.outline.set_linewidth(0.6)
    axd.set_title("(d)  Fissuration  ($t$ = " + vg(tf) + " ms)", loc="left")

    # la fenetre de (d) suit la fissure ; celle de (c) reste sur le forage
    zd = max(0.08, 1.12 * float(np.abs((JP[JC[m]] - [CX, CY])).max()))
    for ax, z in ((axc, zc), (axd, zd)):
        ax.set_xlim(-z * 1e3, z * 1e3)
        ax.set_ylim(-z * 1e3, z * 1e3)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$   [mm]")
        ax.set_ylabel(r"$y$   [mm]")
        repere(ax)

    reg = "anisotrope" if abs(a.sh - a.sv) > 1e-9 else "isotrope"
    fig.suptitle("Fracturation hydraulique au forage — état de contrainte "
                 + reg + r"   ($\sigma^{\prime}_H$ = " + vg(a.sh, 1)
                 + r" MPa,  $\sigma^{\prime}_h$ = " + vg(a.sv, 1) + " MPa)",
                 fontsize=12.5, y=0.972)
    fig.text(0.5, 0.933, str(nb[-1]) + " joints rompus,  " + str(nw[-1])
             + " faces mouillées  (105 au départ)", ha="center",
             fontsize=9.5, color="0.35")

    stem = a.stem or os.path.join(
        ROOT, "bench_abuaisha",
        "fig_b2_" + os.path.basename(run).replace("out_hfs_", ""))
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    fig.savefig(stem + ".png", bbox_inches="tight")
    print("écrit :", stem + ".pdf  et  .png")
    print("  pic %s MPa à t = %s ms | seuil analytique %s MPa | écart %+.1f %%"
          % (vg(p[ip]), vg(t[ip], 3), vg(pcib, 1),
             (p[ip] - pcib) / pcib * 100.0))


if __name__ == "__main__":
    main()
