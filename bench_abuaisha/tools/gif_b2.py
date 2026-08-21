#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# gif_b2.py — ANIMATION du benchmark B2 (AbuAisha et al. 2017), trois panneaux.
#
#   python bench_abuaisha/tools/gif_b2.py out_hfs_aniso --sh 6.8 --sv 4.6
#
#   1. la pression de puits, avec le curseur du temps
#   2. la contrainte principale majeure
#   3. la fissuration (endommagement des joints)
#
# ECONOMIE DE LECTURE. Les VTU font 30 Mo piece et il y en a deux par trame.
# La geometrie (points, connectivite) est IDENTIQUE d'une trame a l'autre —
# les deplacements valent quelques micrometres pour une fenetre de 600 mm,
# donc invisibles. On ne la lit donc QU'UNE FOIS, et on ne relit ensuite que
# les champs. C'est ce qui fait tenir l'animation en quelques minutes au lieu
# d'une heure.
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
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import PolyCollection, LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
from plot_tunnel_fields import complete  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 10,
    "axes.titlesize": 10.5,
    "axes.linewidth": 0.8,
})

CX, CY, RB = 4.0, 4.0, 0.05
FT = 5.0e6


def vg(x, n=2):
    return (("%." + str(n) + "f") % x).replace(".", ",")


def grab(txt, name):
    m = re.search('Name="' + name + '"[^>]*>(.*?)</DataArray>', txt, re.S)
    return np.fromstring(m.group(1), sep=" ")


def points(txt):
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--sh", type=float, default=6.8)
    ap.add_argument("--sv", type=float, default=4.6)
    ap.add_argument("--zoom", type=float, default=0.30)
    ap.add_argument("--fps", type=float, default=8.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(ROOT, a.run)

    h = list(csv.DictReader(open(os.path.join(run, "history.csv"))))
    th = np.array([float(r["t"]) for r in h]) * 1e3
    ph = np.array([float(r["hydroP"]) for r in h]) / 1e6
    nbh = np.array([int(r["nBroken"]) for r in h])
    nwh = np.array([int(float(r["hydroNWet"])) for r in h])
    ipk = int(np.argmax(ph))
    pcib = -a.sh + 3.0 * a.sv + FT / 1e6

    frw = list(csv.DictReader(open(os.path.join(run, "frames.csv"))))
    tfr = np.array([float(r["t"]) for r in frw]) * 1e3

    def kof(f):
        # rsplit : "fdem_0043.vtu" ET "fdem_joints_0043.vtu" donnent 43
        return int(os.path.basename(f).rsplit("_", 1)[1].split(".")[0])

    fs = sorted([f for f in glob.glob(os.path.join(run, "fdem_[0-9]*.vtu"))
                 if complete(f)], key=kof)
    js = {kof(f): f for f in glob.glob(
        os.path.join(run, "fdem_joints_[0-9]*.vtu")) if complete(f)}
    ks = [kof(f) for f in fs if kof(f) in js]
    print("trames exploitables :", len(ks))

    # --- geometrie, lue UNE FOIS -------------------------------------------
    t0 = io.open(fs[0], errors="ignore").read()
    P = points(t0)
    C = grab(t0, "connectivity").astype(int).reshape(-1, 3)
    cen = P[C].mean(axis=1)
    sel = np.where((np.abs(cen[:, 0] - CX) < a.zoom) &
                   (np.abs(cen[:, 1] - CY) < a.zoom))[0]
    verts = [(P[C[i]] - [CX, CY]) * 1e3 for i in sel]

    tj = io.open(js[ks[0]], errors="ignore").read()
    JP = points(tj)
    JC = grab(tj, "connectivity").astype(int).reshape(-1, 2)
    JSEG = (JP[JC] - [CX, CY]) * 1e3
    del t0, tj

    # --- echelle des couleurs, figee sur la trame du pic --------------------
    kpic = min(ks, key=lambda q: abs(tfr[q] - th[ipk]))
    tx = io.open([f for f in fs if kof(f) == kpic][0], errors="ignore").read()
    sxx, syy, sxy = grab(tx, "sigmaXX"), grab(tx, "sigmaYY"), grab(tx, "sigmaXY")
    s1 = 0.5 * (sxx + syy) + np.hypot(0.5 * (sxx - syy), sxy)
    vmax = float(np.percentile(np.abs(s1[sel] / 1e6), 98.0))
    del tx

    # --- mise en place ------------------------------------------------------
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.4, 5.1))
    fig.subplots_adjust(left=0.045, right=0.955, top=0.855,
                        bottom=0.125, wspace=0.34)

    ax1.plot(th, ph, lw=1.6, color="#1f4e79")
    ax1.axhline(pcib, color="0.45", ls=(0, (4, 3)), lw=0.9)
    ax1.text(th[1], pcib, "  seuil analytique  " + vg(pcib, 1) + " MPa",
             color="0.35", fontsize=8.5, va="bottom")
    ax1.plot(th[ipk], ph[ipk], "o", ms=6, mfc="none", mec="#a11122", mew=1.4)
    ax1.annotate(vg(ph[ipk]) + " MPa", (th[ipk], ph[ipk]),
                 textcoords="offset points", xytext=(-10, 8), ha="right",
                 fontsize=9.5, color="#a11122")
    cur, = ax1.plot([], [], "o", ms=7, color="#1f4e79")
    vline = ax1.axvline(th[0], color="0.6", lw=0.8)
    ax1.set_xlabel("temps  [ms]")
    ax1.set_ylabel("pression de puits  [MPa]")
    ax1.set_title("Pression de puits", loc="left")
    ax1.grid(alpha=0.22, lw=0.6)
    ax1.set_xlim(th[0], th[-1])

    pc = PolyCollection(verts, cmap="RdBu_r", lw=0, clim=(-vmax, vmax),
                        array=np.zeros(len(sel)))
    ax2.add_collection(pc)
    thc = np.linspace(0.0, 2.0 * np.pi, 240)
    ax2.plot(RB * 1e3 * np.cos(thc), RB * 1e3 * np.sin(thc), "k-", lw=0.7)
    cb = fig.colorbar(pc, ax=ax2, fraction=0.046, pad=0.03)
    cb.set_label(r"$\sigma_1$   [MPa]")
    cb.outline.set_linewidth(0.6)
    ax2.set_title("Contrainte principale majeure", loc="left")

    ax3.add_patch(plt.Circle((0.0, 0.0), RB * 1e3, fc="0.94", ec="0.55",
                             lw=0.7, zorder=1))
    lc = LineCollection(np.zeros((0, 2, 2)), cmap="Greys", lw=1.6,
                        clim=(-0.25, 1.0), zorder=2)
    ax3.add_collection(lc)
    cb = fig.colorbar(lc, ax=ax3, fraction=0.046, pad=0.03)
    cb.set_label(r"endommagement du joint   $D$")
    cb.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    cb.outline.set_linewidth(0.6)
    ax3.set_title("Fissuration", loc="left")

    for ax in (ax2, ax3):
        ax.set_xlim(-a.zoom * 1e3, a.zoom * 1e3)
        ax.set_ylim(-a.zoom * 1e3, a.zoom * 1e3)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$   [mm]")
        ax.set_ylabel(r"$y$   [mm]")
        ax.annotate("", xy=(0.97, 0.07), xytext=(0.80, 0.07),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="<->", lw=0.9, color="0.35"))
        ax.text(0.885, 0.10, r"$\sigma^{\prime}_H$", transform=ax.transAxes,
                ha="center", fontsize=9.5, color="0.35")
        ax.annotate("", xy=(0.07, 0.94), xytext=(0.07, 0.77),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="<->", lw=0.9, color="0.35"))
        ax.text(0.125, 0.855, r"$\sigma^{\prime}_h$", transform=ax.transAxes,
                va="center", fontsize=9.5, color="0.35")

    reg = "anisotrope" if abs(a.sh - a.sv) > 1e-9 else "isotrope"
    sup = fig.suptitle("", fontsize=12, y=0.965)

    def draw(i):
        k = ks[i]
        tx = io.open([f for f in fs if kof(f) == k][0], errors="ignore").read()
        sxx = grab(tx, "sigmaXX")
        syy = grab(tx, "sigmaYY")
        sxy = grab(tx, "sigmaXY")
        s1 = 0.5 * (sxx + syy) + np.hypot(0.5 * (sxx - syy), sxy)
        pc.set_array(s1[sel] / 1e6)

        d = grab(io.open(js[k], errors="ignore").read(), "damage")
        msk = d > 0.01
        lc.set_segments(JSEG[msk])
        lc.set_array(d[msk])

        tt = tfr[k]
        j = int(np.argmin(np.abs(th - tt)))
        cur.set_data([th[j]], [ph[j]])
        vline.set_xdata([th[j], th[j]])
        sup.set_text("Fracturation hydraulique au forage, état de contrainte "
                     + reg + "   —   $t$ = " + vg(tt) + " ms,   $p$ = "
                     + vg(ph[j]) + " MPa,   " + str(nbh[j])
                     + " joints rompus,   " + str(nwh[j]) + " faces mouillées")
        if i % 6 == 0:
            print("  trame %d/%d" % (i + 1, len(ks)))
        return pc, lc, cur, vline, sup

    out = a.out or os.path.join(ROOT, "bench_abuaisha", "gif_b2_"
                                + os.path.basename(run).replace("out_hfs_", "")
                                + ".gif")
    ani = FuncAnimation(fig, draw, frames=len(ks), blit=False)
    ani.save(out, writer=PillowWriter(fps=a.fps), dpi=95)
    print("écrit :", out, "(%.1f Mo)" % (os.path.getsize(out) / 1e6))


if __name__ == "__main__":
    main()
