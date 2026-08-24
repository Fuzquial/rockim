#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_vitesse_i3d.py — la vitesse de l'insert au cours de l'impact 3D, avec
# les evenements marques, et le bilan d'energie cinetique associe.
#
#   python bench_impact/tools/fig_vitesse_i3d.py out_imp3d_dfh \
#          --stem bench_impact/fig_vitesse_i3d
# ---------------------------------------------------------------------------
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_vitesse_i3d")
    ap.add_argument("--masse", type=float, default=1.18)
    a = ap.parse_args()

    h = np.genfromtxt(a.run + "/history.csv", delimiter=",", names=True)
    t, vz, F, z = h["t"], h["toolVz"], h["toolFz"], h["toolZ"]
    tus = t * 1e6
    k0 = int(np.argmax(np.abs(F) > 1.0e3))
    ke = 0.5 * a.masse * vz ** 2

    fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.6))
    A, B = ax

    A.plot(tus, vz, color="#1f4e79", lw=1.8)
    A.axhline(0, color="k", lw=0.7)
    A.axvline(tus[k0], color="0.55", ls="--", lw=0.9)
    A.annotate("premier contact\n%.1f µs" % tus[k0], (tus[k0], vz[k0]),
               textcoords="offset points", xytext=(10, -26), fontsize=9,
               color="0.35")
    A.annotate("v initiale %.2f m/s" % abs(vz[0]), (tus[3], vz[3]),
               textcoords="offset points", xytext=(12, 8), fontsize=9)
    A.annotate("v finale %.2f m/s\n(TOUJOURS descendante :\naucun rebond)"
               % abs(vz[-1]), (tus[-1], vz[-1]),
               textcoords="offset points", xytext=(-150, 16), fontsize=9,
               color="#b3202f")
    # rupture de pente = l instant ou la roche cede en masse
    dv = np.gradient(vz, t)
    kb = int(np.argmax(dv[k0:])) + k0
    A.plot(tus[kb], vz[kb], "o", color="#b3202f", ms=6)
    A.annotate("rupture de pente %.0f µs" % tus[kb], (tus[kb], vz[kb]),
               textcoords="offset points", xytext=(-135, -20), fontsize=9,
               color="#b3202f")
    A.set_xlabel(r"t [$\mu$s]")
    A.set_ylabel(r"$v_z$ de l'insert [m/s]")
    A.set_title("(a)  Vitesse", loc="left", fontsize=11)
    A.grid(alpha=0.25)

    B.plot(tus, ke, color="#1f4e79", lw=1.8, label="énergie cinétique")
    B.fill_between(tus, ke, ke[0], color="#b3202f", alpha=0.12,
                   label="cédée à la roche")
    B.axhline(ke[0], color="0.55", ls="--", lw=0.9)
    B.set_xlabel(r"t [$\mu$s]")
    B.set_ylabel("énergie cinétique de l'insert [J]")
    B.set_title("(b)  %.0f J sur %.0f cédés (%.0f %%)"
                % (ke[0] - ke[-1], ke[0], 100 * (1 - ke[-1] / ke[0])),
                loc="left", fontsize=11)
    B.legend(frameon=False, fontsize=9, loc="center right")
    B.grid(alpha=0.25)

    fig.suptitle("Impact 3D DP-DFH — insert rigide de %.2f kg" % a.masse,
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)
    print("v %.2f -> %.2f m/s | KE %.1f -> %.1f J | rupture de pente %.0f us"
          % (vz[0], vz[-1], ke[0], ke[-1], tus[kb]))
    print("écrit :", a.stem)


if __name__ == "__main__":
    main()
