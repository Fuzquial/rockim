#!/usr/bin/env python3
"""Rejeu du 2026-10-06 (binaire f0209ef) du tunnel EDZ de Wang et al. (2024)
sur maillage REDUIT -- depouillement direct des trames VTU.

    python3 fig_tunnel_rejeu.py RUN1:etiquette1 RUN2:etiquette2 ...

Pour chaque run : derniere trame, joints rompus (damage >= 0,999), rayon
d'EDZ (quantile 95 % et maximum des distances des milieux de joints rompus au
centre (50, 50), memes definitions que tunnel_edz/tools/edz_metrics.py),
longueur fissuree, modes (breakMode : 1 traction, 2 cisaillement), et
deplacement de paroi (noeuds des aretes libres de la cavite, trame 0 ->
derniere trame). Ecrit rejeu_metrics.json et fig_tunnel_rejeu.pdf.
"""
import json
import os
import re
import sys
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.tri import Triangulation

plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix",
                     "font.size": 9, "pdf.fonttype": 42})
HERE = os.path.dirname(os.path.abspath(__file__))
CX, CY = 50.0, 50.0


def arr(txt, name, dtype=float):
    m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
    return None if m is None else np.array(m.group(1).split(), dtype=float).astype(dtype)


def pts(txt):
    m = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', txt, re.S)
    return np.array(m.group(1).split(), dtype=float).reshape(-1, 3)[:, :2]


def frames(run, joints):
    pre = "fdem_joints_" if joints else "fdem_"
    fs = sorted(f for f in os.listdir(run) if f.startswith(pre) and f.endswith(".vtu")
                and (joints or "joints" not in f))
    return [os.path.join(run, f) for f in fs]


def analyse(run):
    fe = frames(run, False)
    fj = frames(run, True)
    t0 = open(fe[0]).read()
    t1 = open(fe[-1]).read()
    P0, P1 = pts(t0), pts(t1)
    C = arr(t1, "connectivity", int).reshape(-1, 3)
    U = np.linalg.norm(P1 - P0, axis=1)
    # aretes libres de la cavite : aretes (en coordonnees initiales) portees
    # par un seul triangle, a moins de 10 m du centre
    key = lambda i: (round(P0[i, 0], 6), round(P0[i, 1], 6))
    cnt = Counter()
    for tri in C:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            cnt[frozenset((key(tri[a]), key(tri[b])))] += 1
    wall = []
    for tri in C:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            e = frozenset((key(tri[a]), key(tri[b])))
            if cnt[e] == 1:
                xm = 0.5 * (P0[tri[a]] + P0[tri[b]])
                if np.hypot(xm[0] - CX, xm[1] - CY) < 10.0:
                    wall += [tri[a], tri[b]]
    wall = np.unique(wall)
    tj = open(fj[-1]).read()
    PJ = pts(tj)
    CJ = arr(tj, "connectivity", int).reshape(-1, 2)
    D = arr(tj, "damage")
    mode = arr(tj, "breakMode", int)
    brk = D >= 0.999
    seg = PJ[CJ[brk]]
    mid = seg.mean(axis=1)
    r = np.hypot(mid[:, 0] - CX, mid[:, 1] - CY)
    L = np.linalg.norm(seg[:, 1] - seg[:, 0], axis=1).sum()
    hist = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",", names=True)
    fr = np.genfromtxt(os.path.join(run, "frames.csv"), delimiter=",", names=True)
    res = dict(run=os.path.basename(run), t_end=float(fr["t"][-1]),
               triangles=int(len(C)), joints=int(len(D)), broken=int(brk.sum()),
               tensile=int(((mode == 1) & brk).sum()), shear=int(((mode == 2) & brk).sum()),
               inserted=int(hist["nInserted"][-1]),
               edz_radius_p95_m=float(np.percentile(r, 95)) if brk.any() else 0.0,
               edz_radius_max_m=float(r.max()) if brk.any() else 0.0,
               edz_halfaxis_x_p95_m=float(np.percentile(np.abs(mid[:, 0] - CX), 95)) if brk.any() else 0.0,
               edz_halfaxis_y_p95_m=float(np.percentile(np.abs(mid[:, 1] - CY), 95)) if brk.any() else 0.0,
               crack_length_m=float(L),
               u_wall_mean_m=float(U[wall].mean()), u_wall_median_m=float(np.median(U[wall])),
               u_wall_p95_m=float(np.percentile(U[wall], 95)), u_max_m=float(U.max()))
    geo = dict(P0=P0, P1=P1, C=C, U=U, seg=seg, mode=mode[brk])
    return res, geo


def main():
    items = [a.split(":", 1) for a in sys.argv[1:]]
    out = []
    geos = []
    for run, lab in items:
        res, geo = analyse(run)
        res["label"] = lab
        out.append(res)
        geos.append(geo)
        print(json.dumps(res, ensure_ascii=False))
    with open(os.path.join(HERE, "rejeu_metrics.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    n = len(items)
    fig, ax = plt.subplots(2, n, figsize=(2.6 * n + 0.6, 5.9), squeeze=False)
    vmax = max(np.percentile(g["U"], 99.5) for g in geos)
    AMP = 5.0
    for k, (g, res) in enumerate(zip(geos, out)):
        P, C = g["P1"], g["C"]
        cen = P[C].mean(axis=1)
        keep = (np.abs(cen[:, 0] - CX) < 24) & (np.abs(cen[:, 1] - CY) < 24)
        tri = Triangulation(P[:, 0] - CX, P[:, 1] - CY, C[keep])
        Uel = g["U"][C[keep]].mean(axis=1)
        a = ax[0, k]
        pc = a.tripcolor(tri, facecolors=Uel, cmap="viridis", vmin=0, vmax=vmax,
                         rasterized=True)
        s_ = g["seg"] - np.array([CX, CY])
        col = np.where(g["mode"] == 2, "#ef3b2c", "#ffffff")
        a.add_collection(LineCollection(s_, colors=col, linewidths=0.35))
        a.set_xlim(-22, 22)
        a.set_ylim(-22, 22)
        a.set_aspect("equal")
        a.set_title("%s\n%d joints rompus, EDZ %s m" % (
            res["label"], res["broken"],
            f"{res['edz_radius_max_m']:.1f}".replace(".", ",")), fontsize=8.5)
        # rendu elements, zoom au rein droit, deplacements amplifies x AMP
        a = ax[1, k]
        P0 = g["P0"]
        Pa = P0 + AMP * (P - P0)
        cz = Pa[C].mean(axis=1)
        kz = (cz[:, 0] - CX > 0) & (cz[:, 0] - CX < 14) & (np.abs(cz[:, 1] - CY + 1) < 7)
        triz = Triangulation(Pa[:, 0] - CX, Pa[:, 1] - CY, C[kz])
        a.tripcolor(triz, facecolors=g["U"][C[kz]].mean(axis=1), cmap="viridis",
                    vmin=0, vmax=vmax, edgecolors="none", rasterized=True)
        a.set_xlim(2, 12)
        a.set_ylim(-6, 4)
        a.set_aspect("equal")
        a.set_facecolor("white")
        a.set_xlabel("x [m]")
        a.set_title("rein droit, éléments (u × 5)", fontsize=8.5)
    ax[0, 0].set_ylabel("y [m]")
    vf = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
    for aa in ax.ravel():
        aa.xaxis.set_major_formatter(vf)
        aa.yaxis.set_major_formatter(vf)
    ax[1, 0].set_ylabel("y [m]")
    cb = fig.colorbar(pc, ax=ax.ravel().tolist(), shrink=0.6, pad=0.02)
    cb.set_label(r"$|\mathbf{u}|$ [m]")
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v:g}".replace(".", ",")))
    fig.savefig(os.path.join(HERE, "fig_tunnel_rejeu.pdf"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(HERE, "fig_tunnel_rejeu.png"), dpi=150, bbox_inches="tight")

    # histoire : joints rompus contre temps depuis le debut du relachement
    fig, a = plt.subplots(figsize=(4.6, 3.0))
    for (run, lab), res in zip(items, out):
        h = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",", names=True)
        ts = cfg_get(os.path.join(run, "config_effective.cfg"), "excavStart")
        a.plot(h["t"] - ts, h["nBroken"], label=lab)
    a.axvspan(0, 0.08, color="0.9")
    a.text(0.04, a.get_ylim()[1] * 0.92, "rampe", ha="center", fontsize=8)
    a.set_xlabel("temps depuis le début du relâchement [s]")
    a.set_ylabel("joints rompus")
    a.set_xlim(left=-0.02)
    a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v:g}".replace(".", ",")))
    a.legend(frameon=False, fontsize=8)
    a.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_tunnel_rejeu_histoire.pdf"))
    fig.savefig(os.path.join(HERE, "fig_tunnel_rejeu_histoire.png"), dpi=150)


def cfg_get(path, key):
    for line in open(path):
        line = line.split("#")[0]
        if "=" in line and line.split("=")[0].strip() == key:
            return float(line.split("=")[1])
    raise KeyError(key)


if __name__ == "__main__":
    main()
