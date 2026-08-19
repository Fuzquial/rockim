#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# gif_fissures3d.py — animation de la fissuration d un run d indentation 3D,
# en DEUX panneaux par trame : coupe axiale et vue de dessus.
#
#   python tunnel_edz/gif_fissures3d.py out_indent3d_grad [--half 12] [--slab 2]
#
# POURQUOI DEUX VUES. En 3D une seule projection ment : la coupe axiale ne voit
# pas les fissures radiales (elles sont PERPENDICULAIRES au plan de coupe et
# n y laissent qu une trace), la vue de dessus ne voit pas la profondeur du
# cratere. Les deux ensemble suffisent a lire le faciete.
#
# Les joints sont des FACETTES triangulaires : on les trace remplies en vue de
# dessus (une fissure radiale est une surface, un nuage de points ne la montre
# pas) et en trace projetee dans la tranche.
#
#   rouge  : rompu en traction        jaune : rompu en cisaillement
#   pale   : endommage (0,05 < D < 1) — la pointe de fissure
#
# Se lance en cours de calcul : ne lit que les trames deja ecrites.
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
from matplotlib.collections import PolyCollection

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.labelsize": 10, "axes.titlesize": 10.5,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

C_TEN, C_SHE, C_DMG = "#C0392B", "#E8B62C", "#D9B8A8"
C_ROCK, C_TOOL = "#EDEAE6", "#2E5E8C"


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
    ap.add_argument("run", nargs="?", default="out_indent3d_grad")
    ap.add_argument("--half", type=float, default=12.0)
    ap.add_argument("--slab", type=float, default=2.0, help="demi-tranche [mm]")
    ap.add_argument("--fps", type=int, default=5)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = os.path.join(ROOT, a.run)

    jts = sorted(f for f in os.listdir(run) if "joints_" in f)
    els = sorted(f for f in os.listdir(run) if re.match(r"fdem3d_\d+\.vtu$", f))
    if not jts:
        sys.exit("aucun VTU de joints dans " + run)
    fr = np.genfromtxt(os.path.join(run, "frames.csv"), delimiter=",",
                       names=True, invalid_raise=False)
    n = min(len(jts), len(fr["t"]))

    # --- repere fixe, pris sur la ROCHE de la trame 0 ----------------------
    # phase est un champ de CELLULE : on remonte aux noeuds des tetras de roche
    # (indexer les points par le masque de cellules serait faux).
    t0 = open(os.path.join(run, els[0]), errors="ignore").read()
    P0 = grab(t0, "points", 3)
    C0 = grab(t0, "connectivity", 4, dtype=np.int64)
    ph = grab(t0, "phase")
    if ph is not None and (ph > 0.5).any():
        zs = P0[np.unique(C0[ph <= 0.5].ravel())][:, 2].max()
    else:
        zs = P0[:, 2].max()
    x0 = float(fr["toolX"][0]) if "toolX" in fr.dtype.names else 0.0
    y0 = float(fr["toolY"][0]) if "toolY" in fr.dtype.names else 0.0

    tmp = os.path.join(HERE, "_gif3d_frames")
    os.makedirs(tmp, exist_ok=True)
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))

    H, S = a.half * 1e-3, a.slab * 1e-3
    for k in range(n):
        tj = open(os.path.join(run, jts[k]), errors="ignore").read()
        JP = grab(tj, "points", 3)
        JC = grab(tj, "connectivity", 3, dtype=np.int64)
        dm = grab(tj, "damage")
        bm = grab(tj, "breakMode")
        if bm is None:
            bm = np.ones(len(JC))
        tri = JP[JC]
        cen = tri.mean(axis=1)
        brk = dm >= 1.0
        dmg = (dm > 0.05) & ~brk

        fig, ax = plt.subplots(1, 2, figsize=(11.4, 5.0))
        for p in ax:
            p.set_facecolor(C_ROCK)

        win = (np.abs(cen[:, 0] - x0) < H) & (np.abs(cen[:, 1] - y0) < H)
        sl = win & (np.abs(cen[:, 1] - y0) < S)

        # ---- (a) coupe axiale --------------------------------------------
        for m, c, z, al in ((dmg & sl, C_DMG, 1, 0.8),
                            (brk & sl & (bm == 2), C_SHE, 2, 0.95),
                            (brk & sl & (bm != 2), C_TEN, 3, 0.95)):
            if m.any():
                pol = np.stack([(tri[m][:, :, 0] - x0) * 1e3,
                                (tri[m][:, :, 2] - zs) * 1e3], axis=-1)
                ax[0].add_collection(PolyCollection(pol, facecolors=c,
                                                    edgecolors="none",
                                                    alpha=al, zorder=z))
        R, th = 4.0, np.linspace(np.pi, 2 * np.pi, 200)
        tz = (float(fr["toolZ"][k]) - zs) * 1e3 if "toolZ" in fr.dtype.names else 0.0
        ax[0].plot(R * np.cos(th), R * np.sin(th) + tz + R, color=C_TOOL,
                   lw=2.0, zorder=4)
        ax[0].axhline(0.0, color="#B8B2AA", lw=0.8, zorder=1)
        ax[0].set_xlim(-a.half, a.half)
        ax[0].set_ylim(-1.5 * a.half, 0.4 * a.half)
        ax[0].set_xlabel("x [mm]")
        ax[0].set_ylabel("z − z$_{surface}$ [mm]")
        ax[0].set_title("(a) coupe axiale, |y| < %.1f mm" % a.slab, fontsize=10.5)

        # ---- (b) vue de dessus -------------------------------------------
        for m, c, z, al in ((dmg & win, C_DMG, 1, 0.8),
                            (brk & win & (bm == 2), C_SHE, 2, 0.95),
                            (brk & win & (bm != 2), C_TEN, 3, 0.95)):
            if m.any():
                ax[1].add_collection(PolyCollection(
                    (tri[m][:, :, :2] - [x0, y0]) * 1e3, facecolors=c,
                    edgecolors="none", alpha=al, zorder=z))
        ax[1].set_xlim(-a.half, a.half); ax[1].set_ylim(-a.half, a.half)
        ax[1].set_xlabel("x [mm]"); ax[1].set_ylabel("y [mm]")
        ax[1].set_title("(b) vue de dessus", fontsize=10.5)

        for p in ax:
            p.set_aspect("equal")
            for s in ("top", "right"):
                p.spines[s].set_visible(False)

        nT = int((brk & (bm != 2)).sum()); nS = int((brk & (bm == 2)).sum())
        fig.suptitle("t = %6.1f µs      %d joints rompus (%d traction, "
                     "%d cisaillement),  %d endommagés"
                     % (fr["t"][k] * 1e6, nT + nS, nT, nS, int(dmg.sum())),
                     fontsize=11.5, y=1.00)
        fig.tight_layout()
        fig.savefig(os.path.join(tmp, "f%03d.png" % k), dpi=105,
                    bbox_inches="tight")
        plt.close(fig)
        print("  trame %2d/%d  t = %6.1f us  %5d rompus, %5d endommages"
              % (k + 1, n, fr["t"][k] * 1e6, nT + nS, int(dmg.sum())))

    out = os.path.join(HERE, (a.out or "gif_fissures_" + a.run.replace("out_", ""))
                       + ".gif")
    try:
        subprocess.run(["magick", "-delay", str(int(100 / a.fps)), "-loop", "0",
                        os.path.join(tmp, "f*.png"), out], check=True)
    except Exception:
        from PIL import Image
        ims = [Image.open(os.path.join(tmp, f))
               for f in sorted(os.listdir(tmp)) if f.endswith(".png")]
        ims[0].save(out, save_all=True, append_images=ims[1:],
                    duration=int(1000 / a.fps), loop=0)
    print("\necrit :", out)


if __name__ == "__main__":
    main()
