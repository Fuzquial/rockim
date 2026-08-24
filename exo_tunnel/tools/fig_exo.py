#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_exo.py — le champ d'endommagement DP-DFH de la cavite pressurisee,
# rendu comme leur fig_bench6_montage (Abaqus).
#
#   python exo_tunnel/tools/fig_exo.py out_exo_t003 --stem exo_tunnel/fig_t003
#
# Trace DMAX = max(D1,D2,D3) (le SDV 2 de la VUMAT) sur la configuration
# NON deformee, echelle 0-1, et compte les ARMES RADIALES par la methode
# angulaire du banc 6 : histogramme de l endommagement sur une couronne, un
# pic = une arme.
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import os
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

CX = CY = 0.1          # centre du trou [m]
R = 0.01               # rayon du trou


def read(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    P = np.fromstring(re.search(
        r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").reshape(-1, 3)[:, :2]
    tri = np.fromstring(re.search(
        r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s,
        re.S).group(1), sep=" ").astype(int).reshape(-1, 3)
    D = np.fromstring(re.search(
        r'Name="dfhD"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S).group(1),
        sep=" ")
    return P, tri, D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_exo")
    ap.add_argument("--rlo", type=float, default=0.015)
    ap.add_argument("--rhi", type=float, default=0.045)
    a = ap.parse_args()

    fs = [f for f in sorted(glob.glob(a.run + "/fdem_[0-9]*.vtu"))
          if "joints" not in os.path.basename(f)]
    P0, tri, _ = read(fs[0])
    _, _, D = read(fs[-1])
    ctr = P0[tri].mean(axis=1)
    x, y = (ctr[:, 0] - CX) * 1e3, (ctr[:, 1] - CY) * 1e3      # mm
    r = np.hypot(x, y)
    th = np.degrees(np.arctan2(y, x)) % 360.0

    # --- (a) carte d endommagement (rendu ELEMENTS, comme leurs champs) ----
    fig = plt.figure(figsize=(11.8, 5.2))
    A = fig.add_subplot(1, 2, 1)
    xt, yt = (P0[:, 0] - CX) * 1e3, (P0[:, 1] - CY) * 1e3
    tp = A.tripcolor(xt, yt, tri, D, cmap="inferno", vmin=0, vmax=1,
                     shading="flat", rasterized=True)
    cb = fig.colorbar(tp, ax=A, pad=0.02)
    cb.set_label(r"$D_{\max}$")
    thc = np.linspace(0, 2 * np.pi, 200)
    A.plot(R * 1e3 * np.cos(thc), R * 1e3 * np.sin(thc), color="w", lw=0.8)
    A.set_xlim(-60, 60)
    A.set_ylim(-60, 60)
    A.set_aspect("equal")
    A.set_xlabel("x [mm]")
    A.set_ylabel("y [mm]")
    A.set_title("(a)  Endommagement DP-DFH", loc="left", fontsize=11)

    # --- (b) comptage ANGULAIRE des armes (methode du banc 6) -------------
    B = fig.add_subplot(1, 2, 2)
    m = (r >= a.rlo * 1e3) & (r <= a.rhi * 1e3)
    nb = 180
    bins = np.linspace(0, 360, nb + 1)
    idx = np.clip(np.digitize(th[m], bins) - 1, 0, nb - 1)
    prof = np.zeros(nb)
    for k, dv in zip(idx, D[m]):
        prof[k] = max(prof[k], dv)
    # une arme = un groupe angulaire contigu au-dessus du seuil
    seuil = 0.5
    up = prof > seuil
    narm = 0
    for i in range(nb):
        if up[i] and not up[i - 1]:
            narm += 1
    B.plot(0.5 * (bins[:-1] + bins[1:]), prof, color="#b3202f", lw=1.2)
    B.axhline(seuil, color="0.5", ls="--", lw=0.8)
    B.set_xlabel(r"$\theta$ [deg]")
    B.set_ylabel(r"$D_{\max}$ sur la couronne %.0f-%.0f mm"
                 % (a.rlo * 1e3, a.rhi * 1e3))
    B.set_ylim(0, 1.05)
    B.set_xlim(0, 360)
    B.set_title("(b)  %d armes radiales (leur banc 6 : 4 a 15 selon le taux)"
                % narm, loc="left", fontsize=11)
    B.grid(alpha=0.25)

    frac = float((D[m] > 0.5).mean())
    fig.suptitle("Cavite pressurisee 250 MPa — rockim, FEM pur + law = dpdfh "
                 "(%d elements, fraction endommagee %.3f)"
                 % (len(D), frac), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)
    print("ecrit : %s | armes = %d | D>0,5 sur couronne = %.3f | D max %.4f"
          % (a.stem, narm, frac, D.max()))


if __name__ == "__main__":
    main()
