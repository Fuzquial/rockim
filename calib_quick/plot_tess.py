#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_tess.py — planche des tessellations GBM (frame 0 des VTU rockim) :
# triangles colores par phase, joints de grain en noir, histogramme des
# tailles de grain. Usage :
#   python calib_quick/plot_tess.py out_dir:titre [out_dir:titre ...] [--out fig.png]
# ---------------------------------------------------------------------------
import math
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection

plt.rcParams.update({"font.family": "serif", "font.size": 9})
PHASE_COL = ["#f1d9a5", "#9ecae1", "#6b4c2a", "#c7e9c0", "#fdae6b"]
PHASE_NAME = ["feldspath", "quartz", "biotite"]


def read_vtu(path):
    txt = open(path, encoding="utf-8", errors="replace").read()

    def arr(name, dtype=float):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        return np.array(m.group(1).split(), dtype=dtype)
    pts = arr("Points" if 'Name="Points"' in txt else "").reshape(-1, 3)[:, :2] if False else None
    m = re.search(r'<Points>\s*<DataArray[^>]*>(.*?)</DataArray>', txt, re.S)
    pts = np.array(m.group(1).split(), dtype=float).reshape(-1, 3)[:, :2]
    conn = arr("connectivity", int); off = arr("offsets", int)
    tri = np.array([conn[o - 3:o] for o in off if True])
    grain = arr("grain", int); phase = arr("phase", int)
    return pts, tri, grain, phase


def grain_stats(pts, tri, grain):
    u = pts[tri[:, 1]] - pts[tri[:, 0]]; v = pts[tri[:, 2]] - pts[tri[:, 0]]
    a = np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]) / 2
    ga = np.zeros(grain.max() + 1)
    np.add.at(ga, grain, a)
    return 2 * np.sqrt(ga[ga > 0] / math.pi)


def boundaries(pts, tri, grain):
    key = {}
    for t, g in zip(tri, grain):
        for k in range(3):
            P, Q = pts[t[k]], pts[t[(k + 1) % 3]]
            kk = tuple(sorted((tuple(np.round(P, 7)), tuple(np.round(Q, 7)))))
            key.setdefault(kk, []).append(g)
    return [np.array(kk) for kk, gs in key.items() if len(set(gs)) > 1 or len(gs) == 1]


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    argv = sys.argv[1:]
    out = "fig_tess.png"
    if "--out" in argv:
        k = argv.index("--out"); out = argv[k + 1]; argv = argv[:k] + argv[k + 2:]
    args = [a for a in argv if not a.startswith("--")]
    n = len(args)
    fig, axs = plt.subplots(2, n, figsize=(3.4 * n, 8.2), gridspec_kw={"height_ratios": [3, 1]})
    axs = np.array(axs).reshape(2, n)
    for c, spec in enumerate(args):
        d, title = spec.split(":", 1)
        vtu = os.path.join(d, "fdem_0000.vtu")
        pts, tri, grain, phase = read_vtu(vtu)
        ax = axs[0, c]
        ax.add_collection(PolyCollection(pts[tri], facecolors=[PHASE_COL[p] for p in phase], edgecolors="none"))
        ax.add_collection(LineCollection(boundaries(pts, tri, grain), colors="k", linewidths=0.6))
        ax.set_aspect("equal"); ax.set_xlim(pts[:, 0].min(), pts[:, 0].max()); ax.set_ylim(pts[:, 1].min(), pts[:, 1].max())
        ax.set_xticks([]); ax.set_yticks([])
        deq = grain_stats(pts, tri, grain)
        sd = np.log(deq).std()
        ax.set_title(f"{title}\n{len(deq)} grains, d_eq {deq.mean()*1e3:.2f} mm, sd ln d = {sd:.2f}", fontsize=9)
        h = axs[1, c]
        h.hist(deq * 1e3, bins=np.linspace(0, max(8, deq.max() * 1e3), 25), color="0.5")
        h.set_xlabel("d_eq [mm]"); h.set_ylabel("grains" if c == 0 else "")
        print(f"{title:40s} {len(deq):4d} grains  d_eq {deq.mean()*1e3:.2f} ± {deq.std()*1e3:.2f} mm  sd ln d {sd:.3f}  min {deq.min()*1e3:.2f} max {deq.max()*1e3:.2f}")
    from matplotlib.patches import Patch
    axs[0, 0].legend(handles=[Patch(color=PHASE_COL[i], label=PHASE_NAME[i]) for i in range(3)], loc="lower left", fontsize=7, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("figure :", out)


if __name__ == "__main__":
    main()
