# -*- coding: utf-8 -*-
"""Planche de facies de la base LHS : modes de rupture obtenus selon les
parametres (classification de Jiang et al. 2025 : fissure unique, V conjugue,
fissures multiples non paralleles, fendage vertical)."""
import csv, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUNS = os.path.join(BASE, "runs")


def arr(t, n):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    return np.fromstring(m.group(1), sep=" ") if m else None


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def last_frame(d, pref="fdem"):
    fs = sorted(glob.glob(os.path.join(d, "%s_[0-9][0-9][0-9][0-9].vtu" % pref)))
    return fs[-1] if fs else None


def draw(ax, run, title):
    d = os.path.join(RUNS, run)
    f0 = os.path.join(d, "fdem_0000.vtu")
    fl = last_frame(d)
    if not (os.path.exists(f0) and fl):
        ax.set_axis_off(); ax.set_title(title + "\n(pas de sortie)", fontsize=8)
        return 0
    t0 = open(f0).read()
    conn = arr(t0, "connectivity").astype(int).reshape(-1, 3)
    P = pts(open(fl).read())
    ax.add_collection(PolyCollection([P[c] * 1e3 for c in conn],
                                     facecolors="0.88", edgecolors="0.6",
                                     linewidths=0.08))
    n = 0
    jf = fl.replace("fdem_", "fdem_joints_")
    if os.path.exists(jf):
        jt = open(jf).read()
        jP, jc = pts(jt), arr(jt, "connectivity").astype(int)
        tb, off = arr(jt, "tBreak"), arr(jt, "offsets").astype(int)
        segs, s = [], 0
        for i, e in enumerate(off):
            idx = jc[s:e]; s = e
            if tb is not None and tb[i] >= 0 and len(idx) >= 2:
                segs.append(jP[idx[:2]] * 1e3)
        n = len(segs)
        if segs:
            ax.add_collection(LineCollection(segs, colors="crimson", lw=1.1))
    P0 = pts(t0)
    ax.set_xlim(P0[:, 0].min() * 1e3 - 1, P0[:, 0].max() * 1e3 + 1)
    ax.set_ylim(P0[:, 1].min() * 1e3 - 1, P0[:, 1].max() * 1e3 + 1)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("%s\n%d joints rompus" % (title, n), fontsize=8)
    return n


def main():
    # jeux LHS termines, tries par UCS
    rows = []
    for lg in sorted(glob.glob(os.path.join(RUNS, "L*_ucs_s4211.log"))):
        tag = os.path.basename(lg).split("_")[0]
        txt = open(lg, errors="replace").read()
        m = re.search(r"peak axial stress\s*=\s*([\d.]+)", txt)
        if m:
            rows.append((float(m.group(1)), tag))
    rows.sort()
    if len(rows) < 6:
        print("pas assez de runs finis (%d)" % len(rows)); return
    pick = [rows[0], rows[len(rows) // 5], rows[2 * len(rows) // 5],
            rows[3 * len(rows) // 5], rows[4 * len(rows) // 5], rows[-1]]

    fig, ax = plt.subplots(2, 6, figsize=(15, 7.2))
    for j, (ucs, tag) in enumerate(pick):
        draw(ax[0][j], "%s_ucs_s4211" % tag, "%s — UCS %.0f MPa" % (tag, ucs))
        draw(ax[1][j], "%s_bts_s4211" % tag, "%s — brésilien" % tag)
    fig.suptitle("Faciès de rupture de la base LHS — de l'éprouvette la plus faible "
                 "(gauche) à la plus résistante (droite)\nhaut : compression simple   "
                 "bas : brésilien", fontsize=11)
    fig.tight_layout()
    p = os.path.join(BASE, "figures", "facies_lhs.png")
    fig.savefig(p, dpi=140)
    print("ecrit", p, "|", len(rows), "jeux disponibles")


if __name__ == "__main__":
    main()
