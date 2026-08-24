#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_fp_i3d.py — la courbe FORCE-PENETRATION de l'impact 3D DP-DFH.
#
#   python bench_impact/tools/fig_fp_i3d.py out_imp3d_dfh --stem bench_impact/fig_fp_i3d
#
# Ici l'outil est RIGIDE : le solveur ecrit directement sa force de contact
# (toolFz) et sa position (toolZ). Pas de reconstruction par quantite de
# mouvement a faire, contrairement aux runs FDEM a corps maille.
# La penetration est mesuree depuis le PREMIER CONTACT (|Fz| > seuil).
# ---------------------------------------------------------------------------
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_fp_i3d")
    a = ap.parse_args()

    h = np.genfromtxt(a.run + "/history.csv", delimiter=",", names=True)
    t, F, z, vz = h["t"], h["toolFz"], h["toolZ"], h["toolVz"]

    k0 = int(np.argmax(np.abs(F) > 1.0e3))            # premier contact
    p = (z[k0] - z) * 1e3                             # penetration [mm]
    Fk = np.abs(F) / 1e3                              # [kN]

    kpk = int(np.argmax(Fk))
    kmax = int(np.argmax(p))                          # penetration maximale
    print("premier contact t = %.1f us | pic %.0f kN a p = %.2f mm"
          % (t[k0] * 1e6, Fk[kpk], p[kpk]))
    print("penetration max %.2f mm a t = %.1f us | v : %.2f -> %.2f m/s"
          % (p[kmax], t[kmax] * 1e6, vz[0], vz[-1]))

    fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.7))
    A, B = ax
    A.plot(t * 1e6, vz, color="#1f4e79", lw=1.5)
    A.axhline(0, color="k", lw=0.5)
    A.axvline(t[k0] * 1e6, color="0.6", ls="--", lw=0.8)
    A.set_xlabel(r"t [$\mu$s]")
    A.set_ylabel(r"$v_z$ de l'outil [m/s]")
    A.set_title("(a)  Vitesse de l'insert", loc="left", fontsize=11)
    A.grid(alpha=0.25)

    seg = np.array([p[k0:], Fk[k0:]]).T.reshape(-1, 1, 2)
    seg = np.concatenate([seg[:-1], seg[1:]], axis=1)
    lc = LineCollection(seg, cmap="viridis", array=t[k0:-1] * 1e6, lw=1.6)
    B.add_collection(lc)
    cb = fig.colorbar(lc, ax=B, pad=0.02)
    cb.set_label(r"t [$\mu$s]")
    B.annotate("pic %.0f kN" % Fk[kpk], (p[kpk], Fk[kpk]),
               textcoords="offset points", xytext=(8, 4), fontsize=9)
    B.set_xlim(min(0, p[k0:].min()) - 0.1, p[k0:].max() * 1.08)
    B.set_ylim(0, Fk[k0:].max() * 1.15)
    B.set_xlabel("pénétration p [mm]")
    B.set_ylabel(r"$|F_z|$ outil-roche [kN]")
    B.set_title("(b)  Force-pénétration (outil rigide, force directe)",
                loc="left", fontsize=11)
    B.grid(alpha=0.25)

    ke0 = 0.5 * 1.18 * vz[0] ** 2
    fig.suptitle("Impact 3D DP-DFH — insert rigide 1,18 kg a %.1f m/s "
                 "(%.0f J) : pic %.0f kN, enfoncement max %.2f mm"
                 % (abs(vz[0]), ke0, Fk[kpk], p[kmax]), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)
    print("écrit :", a.stem)


if __name__ == "__main__":
    main()
