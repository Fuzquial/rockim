#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# wall_history.py — HISTORIQUE U(t) de la paroi d'un run tunnel et VERDICT
# de stabilisation. COMPLEMENT de wall_convergence.py (convergence RADIALE
# finale projetee sur le vrai contour, multi-runs) : ici on suit le TEMPS.
# Pas le maximum du champ (pollue par les debris volants),
# on suit la MOYENNE par secteur des noeuds de paroi (r0 dans [4,5 ; 7,5] m).
#
#   python tunnel_edz/tools/wall_history.py out_tun_ref_stab \
#          --stem tunnel_edz/fig_wall_stab
#
# Verdict : STABILISE si le gain de U moyen de paroi entre 0,8*t_end et t_end
# est < 5 % de la valeur finale (le run de reference a 0,25 s en gagnait 20 %).
# Lecture VTU identique a fig_wang11.py (positions COURANTES des noeuds).
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})

CX = CY = 50.0
SECTEURS = (("voûte", 80, 100), ("rein g.", 170, 190),
            ("rein d.", -10, 10), ("radier", 260, 280))


def vtu_points(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    a = np.fromstring(
        re.search(r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
                  re.S).group(1), sep=" ")
    return a.reshape(-1, 3)[:, :2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_wall_history")
    a = ap.parse_args()

    frames = [f for f in sorted(glob.glob(a.run + "/fdem_[0-9]*.vtu"))
              if "joints" not in os.path.basename(f)]
    ft = {}
    for line in open(a.run + "/frames.csv").read().splitlines()[1:]:
        p = line.split(",")
        ft[int(p[0])] = float(p[1])
    times = np.array([ft[i] for i in range(len(frames))])

    P0 = vtu_points(frames[0])
    r0 = np.hypot(P0[:, 0] - CX, P0[:, 1] - CY)
    ang = np.degrees(np.arctan2(P0[:, 1] - CY, P0[:, 0] - CX)) % 360.0
    wall = (r0 >= 4.5) & (r0 <= 7.5)
    masks = [(nom, wall & (((ang >= a0 % 360) | (ang <= a1)) if a0 < 0
                           else ((ang >= a0) & (ang <= a1))))
             for nom, a0, a1 in SECTEURS]

    hist = {nom: [] for nom, _ in masks}
    hist["paroi (moy.)"] = []
    for f in frames:
        U = np.linalg.norm(vtu_points(f) - P0, axis=1)
        for nom, m in masks:
            hist[nom].append(float(U[m].mean()))
        hist["paroi (moy.)"].append(float(U[wall].mean()))

    # verdict sur la moyenne de paroi
    um = np.array(hist["paroi (moy.)"])
    i80 = int(np.argmin(np.abs(times - 0.8 * times[-1])))
    gain = um[-1] - um[i80]
    pct = 100.0 * gain / um[-1] if um[-1] > 0 else 0.0
    stable = pct < 5.0
    verdict = "STABILISE" if stable else "CONVERGE ENCORE"

    # distribution finale hors debris (r final < 4 m ou U > 0,3 m)
    P1 = vtu_points(frames[-1])
    U1 = np.linalg.norm(P1 - P0, axis=1)
    r1 = np.hypot(P1[:, 0] - CX, P1[:, 1] - CY)
    keep = wall & (r1 >= 4.0) & (U1 <= 3.0 * um[-1] + 0.1)
    p50, p95 = np.percentile(U1[keep], [50, 95])

    fig, ax = plt.subplots(figsize=(7.6, 4.9))
    for nom, _ in masks:
        ax.plot(times, hist[nom], lw=1.1, label=nom)
    ax.plot(times, um, "k-", lw=2.0, label="paroi (moy.)")
    ax.axvspan(0.02, 0.10, color="0.92", zorder=0)
    ax.text(0.06, ax.get_ylim()[1] * 0.02, "relâchement", ha="center",
            fontsize=8, color="0.4")
    ax.axvline(times[i80], color="0.6", ls="--", lw=0.8)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("U moyen du secteur [m]")
    ax.set_title("Convergence de paroi — %s : %s "
                 "(gain %.1f %% sur le dernier cinquième)"
                 % (os.path.basename(a.run.rstrip("/\\")), verdict, pct),
                 fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)

    print("run           : %s (%d frames, t_end %.4f s)"
          % (a.run, len(frames), times[-1]))
    print("VERDICT       : %s (gain paroi %.4f m = %.1f %% du final "
          "entre 0,8 t_end et t_end ; seuil 5 %%)" % (verdict, gain, pct))
    print("U paroi final : moy %.4f | p50 %.4f | p95 %.4f | "
          "max hors débris %.4f m" % (um[-1], p50, p95, U1[keep].max()))
    for nom, _ in masks:
        print("   %-9s : %.4f m" % (nom, hist[nom][-1]))
    print("écrit : %s.pdf/.png" % a.stem)
    return 0 if stable else 1


if __name__ == "__main__":
    sys.exit(main())
