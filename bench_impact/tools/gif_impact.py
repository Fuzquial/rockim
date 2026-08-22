#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# gif_impact.py — l'impact anime (spec 005, WP5) : la vitesse du bit avec
# curseur, les fissures vues de dessus, la coupe verticale.
#
#   python bench_impact/tools/gif_impact.py out_imp_stanne --out impact.gif
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import (Z_SURF, broken, frame_times, frames_of, history,
                     joints_frame, read_vtu)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})
ROUGE, JAUNE = "#b22222", "#e6a817"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out", default="impact.gif")
    ap.add_argument("--fps", type=float, default=6.0)
    a = ap.parse_args()

    h = history(a.run)
    th, vz = h["t"] * 1e6, h["vz_bit"]
    ks = frames_of(a.run)
    tf = frame_times(a.run)
    data = []
    for k in ks:
        pts, con, f = read_vtu(joints_frame(a.run, k))
        c, n, mode, _ = broken(pts, con, f)
        data.append((k, tf.get(k, 0.0) * 1e6, c, mode))
        print("  trame %d/%d : %d joints rompus" % (k, ks[-1], len(c)))

    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.6))
    fig.suptitle("Impact à insert unique — St Anne, 10,66 m/s", fontsize=12)

    def draw(i):
        k, tk, c, mode = data[i]
        for A in ax:
            A.clear()
        A, B, C = ax
        A.plot(th, vz, color="#1f4e79", lw=1.4)
        A.axvline(tk, color="#b22222", lw=1.0)
        A.axhline(0, color="k", lw=0.5)
        A.set_xlabel(r"temps [$\mu$s]")
        A.set_ylabel(r"$v_z$ bit [m/s]")
        A.set_title("vitesse du bit", fontsize=10)
        if len(c):
            s = mode >= 1.5
            B.scatter(c[~s][:, 0] * 1e3, c[~s][:, 1] * 1e3, s=1.5, c=ROUGE, lw=0)
            B.scatter(c[s][:, 0] * 1e3, c[s][:, 1] * 1e3, s=1.5, c=JAUNE, lw=0)
            s5 = np.abs(c[:, 1]) < 0.005
            cc, m5 = c[s5], mode[s5]
            ss = m5 >= 1.5
            C.scatter(cc[~ss][:, 0] * 1e3, (cc[~ss][:, 2] - Z_SURF) * 1e3,
                      s=2.0, c=ROUGE, lw=0)
            C.scatter(cc[ss][:, 0] * 1e3, (cc[ss][:, 2] - Z_SURF) * 1e3,
                      s=2.0, c=JAUNE, lw=0)
        for AX, L in ((B, 45), (C, 45)):
            AX.set_xlim(-L, L)
            AX.set_aspect("equal")
        B.set_ylim(-45, 45)
        C.set_ylim(-45, 8)
        C.axhline(0, color="#333", lw=0.7)
        B.set_title("vue de dessus  ($t$ = %.0f $\\mu$s)" % tk, fontsize=10)
        C.set_title("coupe $|y|<5$ mm", fontsize=10)
        B.set_xlabel("x [mm]"); B.set_ylabel("y [mm]")
        C.set_xlabel("x [mm]"); C.set_ylabel("z [mm]")
        fig.tight_layout(rect=[0, 0, 1, 0.93])

    ani = FuncAnimation(fig, draw, frames=len(data))
    ani.save(a.out, writer=PillowWriter(fps=a.fps), dpi=105)
    print("écrit : %s (%.1f Mo)" % (a.out, os.path.getsize(a.out) / 1e6))


if __name__ == "__main__":
    main()
