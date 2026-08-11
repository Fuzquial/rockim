#!/usr/bin/env python3
"""Force-penetration curves from rockim history files.

Usage: plot_force_penetration.py out_dir[:label] [out_dir2[:label2] ...]
                                 [--H 0.2] [--R 0.015] [--out fp.png]

Penetration is the geometric indentation of the disc's lowest point below the
initial surface: delta = (H + R) - toolY(t). The tool starts with a small gap,
so curves begin slightly left of zero; contact starts at delta = 0. The area
enclosed by each curve is the work exchanged with the rock.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help="out_dir or out_dir:label")
    ap.add_argument("--H", type=float, default=0.2)
    ap.add_argument("--R", type=float, default=0.015)
    ap.add_argument("--out", default="force_penetration.png")
    a = ap.parse_args()

    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]

    for k, spec in enumerate(a.runs):
        run, _, label = spec.partition(":")
        label = label or os.path.basename(run.rstrip("/"))
        h = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",",
                          names=True, invalid_raise=False)
        if "toolY" not in h.dtype.names:
            sys.exit(run + ": no toolY in history (tension run?)")
        F = np.hypot(h["toolFx"], h["toolFy"]) / 1e6          # MN/m
        delta = (a.H + a.R - h["toolY"]) * 1e3                # mm
        # align delta = 0 on first actual contact (discrete packings top out
        # slightly below the nominal surface, shifting the raw curve)
        touch = np.nonzero(F > 0.005 * F.max())[0]
        if len(touch):
            delta -= delta[touch[0]]
        work = h["work"][-1]
        v0, vf = h["toolVy"][0], h["toolVy"][-1]
        rest = (" e~%.2f" % abs(vf / v0)) if vf > 0 else " no rebound"
        ax.plot(delta, F, lw=1.4, color=colors[k % 4],
                label="%s  (W = %.0f J/m,%s)" % (label, work, rest))

    ax.axvline(0, color="0.8", lw=0.8, zorder=0)
    ax.set_xlabel("penetration of disc lowest point below initial surface [mm]")
    ax.set_ylabel("|F| on tool [MN/m]")
    ax.set_title("Force-penetration, identical insert (disc R = 15 mm, 5 kg, 8 m/s)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
