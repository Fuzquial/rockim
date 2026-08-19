#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# gif_fissures.py — animation de la fissuration d un run d indentation 2D.
#
#   python tunnel_edz/gif_fissures.py out_indent2d_yan [--half 14] [--fps 6]
#
# CE QUE CHAQUE TRAME MONTRE :
#   fond gris    : les elements de roche DEFORMES (on voit donc le cratere se
#                  creuser, pas seulement les fissures apparaitre)
#   trait rouge  : joint ROMPU en traction     (D = 1, breakMode 1)
#   trait jaune  : joint ROMPU en cisaillement (D = 1, breakMode 2)
#   trait pale   : joint ENDOMMAGE mais pas rompu (0,05 < D < 1)
#   arc bleu     : l outil, a sa position lue dans frames.csv
#
# POURQUOI TRACER AUSSI LES ENDOMMAGES : la zone de process d un joint cohesif
# a une longueur finie (l_cz = E Gf / ft^2 = 1,92 mm ici). Ne montrer que les
# D = 1 fait croire a une fissure qui avance par sauts alors qu elle avance
# continument — l endommage est la pointe de fissure.
#
# Le repere est recentre sur l outil et sur la surface libre, de sorte que
# l origine soit le point d impact.
# ---------------------------------------------------------------------------
import argparse
import os
import re
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.labelsize": 11, "axes.titlesize": 12,
    "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "legend.fontsize": 9.5,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

C_TEN, C_SHE, C_DMG = "#C0392B", "#E8B62C", "#D9B8A8"
C_ROCK, C_EDGE, C_TOOL = "#D8D4CF", "#B8B2AA", "#2E5E8C"


def grab(txt, name, ncomp=1, dtype=float):
    if name == "points":
        m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    else:
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
    if m is None:
        return None
    a = np.fromstring(m.group(1), sep=" ", dtype=dtype)
    return a.reshape(-1, ncomp) if ncomp > 1 else a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_indent2d_yan")
    ap.add_argument("--half", type=float, default=13.0, help="demi-fenetre [mm]")
    ap.add_argument("--fps", type=int, default=6)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = os.path.join(ROOT, a.run)

    els = sorted(f for f in os.listdir(run)
                 if re.match(r"fdem_\d+\.vtu$", f) or re.match(r"fdem3d_\d+\.vtu$", f))
    jts = sorted(f for f in os.listdir(run) if "joints_" in f)
    if not jts:
        sys.exit("aucun VTU de joints dans " + run)
    fr = np.genfromtxt(os.path.join(run, "frames.csv"), delimiter=",",
                       names=True, invalid_raise=False)

    # --- repere : outil a la trame 0, surface libre de la roche -------------
    t0 = open(os.path.join(run, els[0]), errors="ignore").read()
    P0 = grab(t0, "points", 3)
    x0 = float(fr["toolX"][0])
    ys = P0[:, 1].max()

    tmp = os.path.join(HERE, "_gif_frames")
    os.makedirs(tmp, exist_ok=True)
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))

    H = a.half
    nb_hist = []
    for k, (fe, fj) in enumerate(zip(els, jts)):
        te = open(os.path.join(run, fe), errors="ignore").read()
        tj = open(os.path.join(run, fj), errors="ignore").read()
        P = grab(te, "points", 3)
        C = grab(te, "connectivity", 3, dtype=np.int64)
        JP = grab(tj, "points", 3)
        JC = grab(tj, "connectivity", 2, dtype=np.int64)
        dm = grab(tj, "damage")
        bm = grab(tj, "breakMode")
        if bm is None:
            bm = np.ones(len(JC))

        fig, ax = plt.subplots(figsize=(7.6, 6.4))
        ax.set_facecolor("white")

        # --- roche deformee -------------------------------------------------
        T = P[C]
        cen = T.mean(axis=1)
        keep = (np.abs(cen[:, 0] - x0) < 1.4 * H * 1e-3) & \
               (cen[:, 1] > ys - 2.2 * H * 1e-3)
        pol = np.stack([(T[keep][:, :, 0] - x0) * 1e3,
                        (T[keep][:, :, 1] - ys) * 1e3], axis=-1)
        ax.add_collection(PolyCollection(pol, facecolors=C_ROCK,
                                         edgecolors=C_EDGE, lw=0.05, zorder=0))

        # --- joints ---------------------------------------------------------
        S = JP[JC]
        seg = np.stack([(S[:, :, 0] - x0) * 1e3, (S[:, :, 1] - ys) * 1e3], -1)
        brk = dm >= 1.0
        dmg = (dm > 0.05) & ~brk
        for m, c, lw, z in ((dmg, C_DMG, 0.9, 1),
                            (brk & (bm != 2), C_TEN, 1.35, 3),
                            (brk & (bm == 2), C_SHE, 1.35, 2)):
            if m.any():
                ax.add_collection(LineCollection(seg[m], colors=c, lw=lw,
                                                 zorder=z, alpha=0.95))

        # --- outil ----------------------------------------------------------
        R, th = 4.0, np.linspace(0, 2 * np.pi, 240)
        ty = (float(fr["toolY"][k]) - ys) * 1e3
        ax.plot(R * np.cos(th), R * np.sin(th) + ty, color=C_TOOL, lw=2.2,
                zorder=4)

        nT = int((brk & (bm != 2)).sum()); nS = int((brk & (bm == 2)).sum())
        nb_hist.append(nT + nS)
        ax.set_xlim(-H, H); ax.set_ylim(-1.55 * H, 0.42 * H)
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]"); ax.set_ylabel("z − z$_{surface}$ [mm]")
        ax.set_title("t = %6.1f µs      %d joints rompus  "
                     "(%d traction, %d cisaillement)"
                     % (fr["t"][k] * 1e6, nT + nS, nT, nS), fontsize=11.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(tmp, "f%03d.png" % k), dpi=110)
        plt.close(fig)
        print("  trame %2d/%d  t = %6.1f us  %5d rompus"
              % (k + 1, len(jts), fr["t"][k] * 1e6, nT + nS))

    out = os.path.join(HERE, (a.out or "gif_fissures_" + a.run.replace("out_", ""))
                       + ".gif")
    # imagemagick si dispo, sinon Pillow
    try:
        subprocess.run(["magick", "-delay", str(int(100 / a.fps)), "-loop", "0",
                        os.path.join(tmp, "f*.png"), out], check=True)
    except Exception:
        from PIL import Image
        ims = [Image.open(os.path.join(tmp, f))
               for f in sorted(os.listdir(tmp)) if f.endswith(".png")]
        ims[0].save(out, save_all=True, append_images=ims[1:],
                    duration=int(1000 / a.fps), loop=0)
    print("\necrit : %s  (%d trames, %d -> %d joints rompus)"
          % (out, len(jts), nb_hist[0], nb_hist[-1]))


if __name__ == "__main__":
    main()
