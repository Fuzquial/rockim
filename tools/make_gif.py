#!/usr/bin/env python3
"""Animated GIF of a rockim run, tool included.

usage:  make_gif.py <config.cfg> <run_dir> [out.gif]

Reads the per-frame VTU files plus frames.csv (frame -> time, tool pose) and
history.csv (tool force, interpolated at frame times for the title), renders
each frame with matplotlib and assembles a GIF with Pillow.

DEM runs : particles coloured by speed, broken bonds drawn in red, tool as a
           grey disc.
FEM runs : surviving mesh coloured by damage (eroded elements = holes), tool
           as a grey disc or flat punch.
"""
import io
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, Rectangle
from PIL import Image


def read_cfg(path):
    cfg = {}
    for line in open(path):
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            cfg[k.strip()] = v.strip()
    return cfg


def vtu_array(text, name, dtype=float):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  text, re.S)
    return np.fromstring(m.group(1), sep=" ", dtype=dtype)


def vtu_points(text):
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def vtu_points3(text):
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)


def draw_tool(ax, cfg, mode, scen, tx, ty):
    if mode == "fem" and scen == "percussion":
        shape = cfg.get("toolShape", "flat")
    else:
        shape = cfg.get("toolShape", "disc")
    if shape == "flat":
        w = float(cfg.get("toolWidth", 0.02))
        ax.add_patch(Rectangle((tx - w / 2, ty), w, 0.012,
                               fc="0.35", ec="k", zorder=5))
    else:
        r = float(cfg.get("toolRadius", 0.01))
        ax.add_patch(Circle((tx, ty), r, fc="0.35", ec="k", zorder=5))


def main():
    cfg_path, run = sys.argv[1], sys.argv[2]
    gif = sys.argv[3] if len(sys.argv) > 3 else run + "/animation.gif"
    cfg = read_cfg(cfg_path)
    mode, scen = cfg.get("mode", "fem"), cfg.get("scenario", "percussion")
    W, H = float(cfg.get("W", 0.2)), float(cfg.get("H", 0.1))

    frames = np.genfromtxt(run + "/frames.csv", delimiter=",", names=True)
    hist = np.genfromtxt(run + "/history.csv", delimiter=",", names=True)
    if "toolFz" in hist.dtype.names:
        fmag = np.sqrt(hist["toolFx"]**2 + hist["toolFy"]**2 + hist["toolFz"]**2)
    else:
        fmag = np.hypot(hist["toolFx"], hist["toolFy"])

    # room above the specimen for the tool and the ejecta
    ylim = (-0.02 * H, H * (1.55 if scen == "percussion" else 1.75))
    xlim = (-0.06 * W, 1.06 * W)

    images = []
    for k in range(len(frames["frame"])):
        t, tx, ty = frames["t"][k], frames["toolX"][k], frames["toolY"][k]
        fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=110)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])

        if mode == "fdem3d":
            # mid-depth slice: element centroids coloured by von Mises,
            # broken-joint centroids in red, sliced spherical tool
            D = float(cfg.get("D", 0.08))
            H = float(cfg.get("H", 0.06))
            nyc = int(cfg.get("ny", 20))
            hm = D / nyc
            yc = 0.5 * D
            ty = frames["toolZ"][k]
            txt = open(run + "/fdem3d_%04d.vtu" % k).read()
            pts = vtu_points3(txt)
            conn = vtu_array(txt, "connectivity", int).reshape(-1, 4)
            svm = vtu_array(txt, "vonMises")
            cen = pts[conn].mean(axis=1)
            sel = np.abs(cen[:, 1] - yc) < 0.7 * hm
            ax.scatter(cen[sel, 0], cen[sel, 2], s=14, c=svm[sel],
                       cmap="Blues", vmin=0,
                       vmax=3 * float(cfg.get("cohesion", 25e6)), lw=0, zorder=2)
            jtxt = open(run + "/fdem3d_joints_%04d.vtu" % k).read()
            jpts = vtu_points3(jtxt)
            jconn = vtu_array(jtxt, "connectivity", int).reshape(-1, 3)
            dmg = vtu_array(jtxt, "damage")
            jcen = jpts[jconn].mean(axis=1)
            br = (dmg >= 0.99) & (np.abs(jcen[:, 1] - yc) < 0.7 * hm)
            if br.any():
                ax.scatter(jcen[br, 0], jcen[br, 2], s=9, color="crimson",
                           zorder=3)
            R = float(cfg.get("toolRadius", 0.015))
            dy = frames["toolY"][k] - yc
            if abs(dy) < R:
                ax.add_patch(Circle((tx, ty), np.sqrt(R * R - dy * dy),
                                    fc="0.35", ec="k", zorder=5))
        elif mode == "dem3d":
            # mid-depth slice |y - D/2| < ~one layer, projected on (x, z)
            D = float(cfg.get("D", 0.1))
            pr = float(cfg.get("particleRadius", 1.5e-3))
            yc, half = 0.5 * D, 1.05 * pr
            ty = frames["toolZ"][k]
            txt = open(run + "/dem3d_particles_%04d.vtu" % k).read()
            pts = vtu_points3(txt)
            spd = vtu_array(txt, "speed")
            rad = vtu_array(txt, "radius")
            sel = np.abs(pts[:, 1] - yc) < half
            fig.canvas.draw()
            per = (ax.transData.transform((1, 0))
                   - ax.transData.transform((0, 0)))[0]
            s = (2 * rad[sel] * per * 72.0 / fig.dpi) ** 2
            ax.scatter(pts[sel, 0], pts[sel, 2], s=s, c=spd[sel],
                       cmap="viridis", vmin=0,
                       vmax=max(1.0, np.percentile(spd, 99)), lw=0, zorder=2)
            btxt = open(run + "/dem3d_bonds_%04d.vtu" % k).read()
            bpts = vtu_points3(btxt)
            conn = vtu_array(btxt, "connectivity", int).reshape(-1, 2)
            state = vtu_array(btxt, "state")
            inb = (np.abs(bpts[conn[:, 0], 1] - yc) < half) \
                  & (np.abs(bpts[conn[:, 1], 1] - yc) < half)
            broken = conn[(state > 0.5) & inb]
            if len(broken):
                segs = np.stack([bpts[broken[:, 0]][:, [0, 2]],
                                 bpts[broken[:, 1]][:, [0, 2]]], 1)
                ax.add_collection(LineCollection(segs, colors="crimson",
                                                 lw=0.5, zorder=3))
            R = float(cfg.get("toolRadius", 0.015))
            dy = frames["toolY"][k] - yc
            if abs(dy) < R:
                ax.add_patch(Circle((tx, ty), np.sqrt(R * R - dy * dy),
                                    fc="0.35", ec="k", zorder=5))
        elif mode == "fdem":
            txt = open(run + "/fdem_%04d.vtu" % k).read()
            pts = vtu_points(txt)
            conn = vtu_array(txt, "connectivity", int).reshape(-1, 3)
            svm = vtu_array(txt, "vonMises")
            ax.tripcolor(pts[:, 0], pts[:, 1], conn, facecolors=svm,
                         cmap="Blues", vmin=0, vmax=3 * float(cfg.get("cohesion", 25e6)))
            jtxt = open(run + "/fdem_joints_%04d.vtu" % k).read()
            jpts = vtu_points(jtxt)
            jconn = vtu_array(jtxt, "connectivity", int).reshape(-1, 2)
            dmg = vtu_array(jtxt, "damage")
            br = dmg >= 0.99
            if br.any():
                segs = np.stack([jpts[jconn[br, 0]], jpts[jconn[br, 1]]], 1)
                ax.add_collection(LineCollection(segs, colors="crimson",
                                                 lw=0.6, zorder=3))
        elif mode == "dem":
            txt = open(run + "/dem_particles_%04d.vtu" % k).read()
            pts = vtu_points(txt)
            spd = vtu_array(txt, "speed")
            rad = vtu_array(txt, "radius")
            # marker size: radius in data units -> points^2
            fig.canvas.draw()
            per = (ax.transData.transform((1, 0))
                   - ax.transData.transform((0, 0)))[0]      # px per data unit
            s = (2 * rad * per * 72.0 / fig.dpi) ** 2
            ax.scatter(pts[:, 0], pts[:, 1], s=s, c=spd, cmap="viridis",
                       vmin=0, vmax=max(1.0, np.percentile(spd, 99)),
                       lw=0, zorder=2)
            btxt = open(run + "/dem_bonds_%04d.vtu" % k).read()
            bpts = vtu_points(btxt)
            conn = vtu_array(btxt, "connectivity", int).reshape(-1, 2)
            state = vtu_array(btxt, "state")
            broken = conn[state > 0.5]
            if len(broken):
                segs = np.stack([bpts[broken[:, 0]], bpts[broken[:, 1]]], 1)
                ax.add_collection(LineCollection(segs, colors="crimson",
                                                 lw=0.5, zorder=3))
        else:  # fem
            txt = open(run + "/fem_%04d.vtu" % k).read()
            pts = vtu_points(txt)
            conn = vtu_array(txt, "connectivity", int).reshape(-1, 3)
            dmg = vtu_array(txt, "damage")
            svm = vtu_array(txt, "vonMises")
            # light-grey backdrop: eroded elements leave holes that read as
            # the excavated crater, distinct from unstressed (white) rock
            ax.add_patch(Rectangle((0, 0), W, H, fc="0.88", ec="none",
                                   zorder=1))
            svmax = 3.0 * float(cfg.get("cohesion", 25e6))
            ax.tripcolor(pts[:, 0], pts[:, 1], conn, facecolors=svm,
                         cmap="Blues", vmin=0, vmax=svmax, zorder=2)
            sel = dmg > 0.02          # overlay: the growing damage/crack field
            if sel.any():
                ax.tripcolor(pts[:, 0], pts[:, 1], conn[sel],
                             facecolors=dmg[sel], cmap="inferno",
                             vmin=0, vmax=1, zorder=3)

        if mode != "dem3d":
            draw_tool(ax, cfg, mode, scen, tx, ty)
        F = np.interp(t, hist["t"], fmag)
        if mode == "fdem3d":
            # mid-depth slice: element centroids coloured by von Mises,
            # broken-joint centroids in red, sliced spherical tool
            D = float(cfg.get("D", 0.08))
            H = float(cfg.get("H", 0.06))
            nyc = int(cfg.get("ny", 20))
            hm = D / nyc
            yc = 0.5 * D
            ty = frames["toolZ"][k]
            txt = open(run + "/fdem3d_%04d.vtu" % k).read()
            pts = vtu_points3(txt)
            conn = vtu_array(txt, "connectivity", int).reshape(-1, 4)
            svm = vtu_array(txt, "vonMises")
            cen = pts[conn].mean(axis=1)
            sel = np.abs(cen[:, 1] - yc) < 0.7 * hm
            ax.scatter(cen[sel, 0], cen[sel, 2], s=14, c=svm[sel],
                       cmap="Blues", vmin=0,
                       vmax=3 * float(cfg.get("cohesion", 25e6)), lw=0, zorder=2)
            jtxt = open(run + "/fdem3d_joints_%04d.vtu" % k).read()
            jpts = vtu_points3(jtxt)
            jconn = vtu_array(jtxt, "connectivity", int).reshape(-1, 3)
            dmg = vtu_array(jtxt, "damage")
            jcen = jpts[jconn].mean(axis=1)
            br = (dmg >= 0.99) & (np.abs(jcen[:, 1] - yc) < 0.7 * hm)
            if br.any():
                ax.scatter(jcen[br, 0], jcen[br, 2], s=9, color="crimson",
                           zorder=3)
            R = float(cfg.get("toolRadius", 0.015))
            dy = frames["toolY"][k] - yc
            if abs(dy) < R:
                ax.add_patch(Circle((tx, ty), np.sqrt(R * R - dy * dy),
                                    fc="0.35", ec="k", zorder=5))
        elif mode == "dem3d":
            ax.set_title("%s %s (mid slice)   t = %.0f µs   |F| = %.1f kN"
                         % (mode.upper(), scen, t * 1e6, F / 1e3), fontsize=10)
        else:
            ax.set_title("%s %s   t = %.0f µs   |F| = %.2f MN/m"
                         % (mode.upper(), scen, t * 1e6, F / 1e6), fontsize=10)
        fig.tight_layout(pad=0.4)

        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        plt.close(fig)
        buf.seek(0)
        images.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))

    images[0].save(gif, save_all=True, append_images=images[1:],
                   duration=90, loop=0)
    print("wrote", gif, "(%d frames)" % len(images))


if __name__ == "__main__":
    main()
