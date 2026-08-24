#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_i3d_dfh.py — l impact 3D DP-DFH : vue de dessus (l ETOILE) et coupe.
#
#   python bench_impact/tools/fig_i3d_dfh.py out_imp3d_dfh --stem bench_impact/fig_i3d
#
# On ne trace que les tetraedres ROMPUS (D >= seuil), projetes : de dessus
# pour lire les branches radiales, en coupe verticale pour le cratere.
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})


def read(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    P = np.fromstring(re.search(
        r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").reshape(-1, 3)
    con = np.fromstring(re.search(
        r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").astype(int).reshape(-1, 4)
    D = np.fromstring(re.search(
        r'Name="(?:dfhD|damage)"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S).group(1),
        sep=" ")
    return P, con, D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_i3d")
    ap.add_argument("--seuil", type=float, default=0.9)
    ap.add_argument("--frame", type=int, default=-1)
    a = ap.parse_args()

    fs = [f for f in sorted(glob.glob(a.run + "/fem3d_[0-9]*.vtu"))
          + sorted(glob.glob(a.run + "/fdem3d_[0-9]*.vtu"))
          if "joints" not in f]
    P0, con, _ = read(fs[0])
    _, _, D = read(fs[a.frame])
    ctr = P0[con].mean(axis=1)
    zs = P0[:, 2].max()                              # surface libre = z max
    x = (ctr[:, 0] - 0.5 * (P0[:, 0].min() + P0[:, 0].max())) * 1e3
    y = (ctr[:, 1] - 0.5 * (P0[:, 1].min() + P0[:, 1].max())) * 1e3
    z = (ctr[:, 2] - zs) * 1e3                        # 0 = surface, <0 dedans
    br = D >= a.seuil

    fig, ax = plt.subplots(1, 2, figsize=(12.0, 5.3))
    A, B = ax
    A.scatter(x[br], y[br], s=3, c=-z[br], cmap="inferno_r", vmin=0, vmax=25,
              linewidths=0)
    th = np.linspace(0, 2 * np.pi, 200)
    for R in (10, 20, 30):
        A.plot(R * np.cos(th), R * np.sin(th), color="0.75", lw=0.6, zorder=0)
    A.set_xlim(-40, 40)
    A.set_ylim(-40, 40)
    A.set_aspect("equal")
    A.set_xlabel("x [mm]")
    A.set_ylabel("y [mm]")
    A.set_title("(a)  Vue de dessus — cercles a 10, 20, 30 mm",
                loc="left", fontsize=11)

    cut = br & (np.abs(y) < 5.0)
    B.scatter(x[cut], z[cut], s=4, c="#b3202f", linewidths=0)
    B.axhline(0, color="k", lw=0.6)
    B.set_xlim(-40, 40)
    B.set_ylim(-40, 5)
    B.set_aspect("equal")
    B.set_xlabel("x [mm]")
    B.set_ylabel("z sous la surface [mm]")
    B.set_title("(b)  Coupe verticale |y| < 5 mm", loc="left", fontsize=11)

    r = np.hypot(x[br], y[br])
    fig.suptitle("Impact 3D DP-DFH (rockim, FEM pur) — %d tets a D >= %.2f, "
                 "extension radiale max %.1f mm, profondeur %.1f mm"
                 % (int(br.sum()), a.seuil, r.max(), -z[br].min()), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)
    print("ecrit : %s | rayon max %.1f mm | profondeur %.1f mm | %d tets"
          % (a.stem, r.max(), -z[br].min(), int(br.sum())))


if __name__ == "__main__":
    main()
