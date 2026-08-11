#!/usr/bin/env python3
"""Quick-look plots from a rockim output directory (no ParaView needed).

usage: python3 tools/plot_results.py <out_dir> [--title "..."]

Produces, depending on what it finds in <out_dir>:
  <out_dir>/plot_field.png    damage map (FEM) or particles + broken bonds (DEM)
  <out_dir>/plot_history.png  tool force-time history
"""
import argparse
import csv
import os
import sys

import glob
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def vtu_array(text, name, dtype=float):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  text, re.S)
    return np.fromstring(m.group(1), sep=" ", dtype=dtype)


def vtu_points(text):
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def read_csv(path):
    with open(path) as f:
        r = csv.reader(f)
        head = next(r)
        rows = [[float(x) for x in row] for row in r if row]
    data = {h: np.array([row[i] for row in rows]) for i, h in enumerate(head)}
    return data


def plot_fem(out, ax):
    # Render the true triangulation of the last frame (the mesh is crossed
    # CST triangles; scattering element centroids, as an earlier version did,
    # produces a misleading dotted/"hexagonal" texture).
    vtus = sorted(glob.glob(os.path.join(out, "fem_[0-9]*.vtu")))
    if vtus:
        txt = open(vtus[-1]).read()
        pts = vtu_points(txt)
        conn = vtu_array(txt, "connectivity", int).reshape(-1, 3)
        dmg = vtu_array(txt, "damage")
        pc = ax.tripcolor(pts[:, 0], pts[:, 1], conn, facecolors=dmg,
                          cmap="inferno", vmin=0, vmax=1)
        plt.colorbar(pc, ax=ax, label="damage D", shrink=0.8)
        ax.set_title("FEM: damage field, deformed mesh (holes = eroded)")
        return
    d = read_csv(os.path.join(out, "fem_final_elements.csv"))
    live = d["eroded"] < 0.5
    sc = ax.scatter(d["cx"][live], d["cy"][live], c=d["damage"][live],
                    s=4, cmap="inferno", vmin=0, vmax=1, marker="s")
    ax.scatter(d["cx"][~live], d["cy"][~live], color="white", s=4, marker="s")
    plt.colorbar(sc, ax=ax, label="damage D", shrink=0.8)
    ax.set_title("FEM: damage field (white = eroded)")


def plot_dem(out, ax):
    p = read_csv(os.path.join(out, "dem_final_particles.csv"))
    b = read_csv(os.path.join(out, "dem_final_bonds.csv"))
    broken = b["broken"] > 0.5
    for x1, y1, x2, y2 in zip(b["x1"][broken], b["y1"][broken],
                              b["x2"][broken], b["y2"][broken]):
        ax.plot([x1, x2], [y1, y2], color="red", lw=0.5, alpha=0.6, zorder=3)
    frag = p["fragment"]
    main = frag == 0
    ax.scatter(p["x"][main], p["y"][main], s=3, color="0.6", zorder=2)
    sc = ax.scatter(p["x"][~main], p["y"][~main], c=frag[~main], s=3,
                    cmap="turbo", zorder=2)
    ax.set_title("DEM: fragments (grey = main body) + broken bonds (red)")


def plot_dem3d(out, ax):
    """Mid-depth slice (|y - yc| < one particle diameter) of the final state."""
    p = read_csv(os.path.join(out, "dem3d_final_particles.csv"))
    yc = 0.5 * (p["y"].min() + p["y"].max())
    sel = np.abs(p["y"] - yc) < 2.05 * p["r"]
    frag = p["fragment"][sel]
    main = frag == 0
    ax.scatter(p["x"][sel][main], p["z"][sel][main], s=6, color="0.6", zorder=2)
    if (~main).any():
        ax.scatter(p["x"][sel][~main], p["z"][sel][~main], c=frag[~main], s=6,
                   cmap="turbo", zorder=3)
    ax.set_title("DEM3D: mid-depth slice, fragments (grey = main body)")


def plot_fdem(out, ax):
    e = read_csv(os.path.join(out, "fdem_final_elements.csv"))
    j = read_csv(os.path.join(out, "fdem_final_joints.csv"))
    frag = e["fragment"]
    main_b = frag == 0
    ax.scatter(e["cx"][main_b], e["cy"][main_b], s=4, color="0.75", zorder=1)
    if (~main_b).any():
        ax.scatter(e["cx"][~main_b], e["cy"][~main_b], c=frag[~main_b], s=5,
                   cmap="turbo", zorder=3)
    br = j["damage"] >= 0.99
    if br.any():
        from matplotlib.collections import LineCollection
        segs = np.stack([np.stack([j["x1"][br], j["y1"][br]], 1),
                         np.stack([j["x2"][br], j["y2"][br]], 1)], 1)
        ax.add_collection(LineCollection(segs, colors="crimson", lw=0.5, zorder=2))
    ax.set_title("FDEM: fragments (grey = main body) + broken joints (red)")


def plot_fdem3d(out, ax):
    e = read_csv(os.path.join(out, "fdem3d_final_elements.csv"))
    yc = 0.5 * (e["cy"].min() + e["cy"].max())
    hm = np.median(np.diff(np.unique(np.round(e["cy"], 5)))) * 3
    sel = np.abs(e["cy"] - yc) < hm
    frag = e["fragment"][sel]
    mb = frag == 0
    ax.scatter(e["cx"][sel][mb], e["cz"][sel][mb], s=8, color="0.75", zorder=1)
    if (~mb).any():
        ax.scatter(e["cx"][sel][~mb], e["cz"][sel][~mb], c=frag[~mb], s=9,
                   cmap="turbo", zorder=3)
    ax.set_title("FDEM3D: mid-depth slice, fragments (grey = main body)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--title", default="")
    a = ap.parse_args()

    is_fem = os.path.exists(os.path.join(a.out, "fem_final_elements.csv"))
    is_dem = os.path.exists(os.path.join(a.out, "dem_final_particles.csv"))
    is_d3 = os.path.exists(os.path.join(a.out, "dem3d_final_particles.csv"))
    is_fd = os.path.exists(os.path.join(a.out, "fdem_final_elements.csv"))
    is_f3 = os.path.exists(os.path.join(a.out, "fdem3d_final_elements.csv"))
    if not (is_fem or is_dem or is_d3 or is_fd or is_f3):
        sys.exit("no rockim final CSVs found in " + a.out)

    fig, ax = plt.subplots(figsize=(9, 5))
    (plot_fem if is_fem else plot_dem if is_dem else
     plot_dem3d if is_d3 else plot_fdem if is_fd else plot_fdem3d)(a.out, ax)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    if a.title:
        fig.suptitle(a.title)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "plot_field.png"), dpi=150)

    hist = os.path.join(a.out, "history.csv")
    if os.path.exists(hist):
        h = read_csv(hist)
        if "toolFz" in h:
            fig2, ax2 = plt.subplots(figsize=(9, 4))
            ax2.plot(h["t"] * 1e3, h["toolFz"] / 1e3, label="F$_z$ on tool")
            ax2.plot(h["t"] * 1e3, h["toolFx"] / 1e3, label="F$_x$ on tool", alpha=0.8)
            ax2.set_xlabel("t [ms]")
            ax2.set_ylabel("force [kN]")
            ax2.legend()
            ax2.grid(alpha=0.3)
            ax2.set_title("Tool force history" + (" — " + a.title if a.title else ""))
            fig2.tight_layout()
            fig2.savefig(os.path.join(a.out, "plot_history.png"), dpi=150)
        elif "toolFy" in h:
            fig2, ax2 = plt.subplots(figsize=(9, 4))
            ax2.plot(h["t"] * 1e3, h["toolFy"] / 1e6, label="F$_y$ on tool")
            ax2.plot(h["t"] * 1e3, h["toolFx"] / 1e6, label="F$_x$ on tool", alpha=0.8)
            ax2.set_xlabel("t [ms]")
            ax2.set_ylabel("force [MN/m]")
            ax2.legend()
            ax2.grid(alpha=0.3)
            ax2.set_title("Tool force history" + (" — " + a.title if a.title else ""))
            fig2.tight_layout()
            fig2.savefig(os.path.join(a.out, "plot_history.png"), dpi=150)


if __name__ == "__main__":
    main()
