#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_xi_compare.py — les deux tests courts "enveloppe de Yang, machinerie de
# Yan" sur p1_ultra (2026-08-19), a la seule difference de la viscosite :
#
#   out_yang_court      bulkViscosityXi = 0.5   (le quart du critique)
#   out_yang_court_xi2  bulkViscosityXi = 2.0   (la valeur de Yan, Table 1)
#
# CE QUE CETTE PLANCHE MONTRE : que le choix de xi ne se joue PAS sur la
# mecanique a ce stade — les deux courbes se superposent — mais sur le cout,
# mesure au chronometre a un facteur 2,37.
#
# CE QU ELLE NE MONTRE PAS : la fracturation. 5 joints casses sur 247 226 dans
# les deux cas, pas de decharge, pas de rebond. T = T_ref / 6.
#
#   python tunnel_edz/fig_xi_compare.py
# ---------------------------------------------------------------------------
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.labelsize": 9.5, "axes.titlesize": 10,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.5,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

RUNS = [("out_yang_court",     r"$\xi = 0{,}5$  (quart du critique)", "#4A6FA5", "-",  730),
        ("out_yang_court_xi2", r"$\xi = 2$  (valeur de Yan)",         "#8C2F27", "--", 1731)]


def load(run):
    return np.genfromtxt(os.path.join(ROOT, run, "history.csv"), delimiter=",",
                         names=True, invalid_raise=False)


def main():
    fig, ax = plt.subplots(1, 3, figsize=(10.8, 3.5))

    for run, lab, c, ls, wall in RUNS:
        h = load(run)
        t = h["t"] * 1e6
        F = h["grpFz"] / 1e3

        # (a) force et vitesse
        ax[0].plot(t, F, color=c, ls=ls, lw=1.6, label=lab)

        # (b) force-penetration, branche de charge
        i0 = int(np.argmax(np.abs(h["grpFz"]) > 1.0))
        d = (h["grpZ"][i0] - h["grpZ"]) * 1e6
        m = np.arange(len(F)) >= i0
        ax[1].plot(d[m], F[m], color=c, ls=ls, lw=1.6, label=lab)

        # (c) poste elements
        ax[2].plot(t, -h["eEl"] * 1e3, color=c, ls=ls, lw=1.6, label=lab)

    ax[0].set_xlabel("temps [µs]")
    ax[0].set_ylabel("force axiale sur l'insert  [kN]")
    ax[0].set_title("(a) les deux courbes se superposent", fontsize=10)
    ax[0].legend(loc="upper left", framealpha=0.95)

    ax[1].set_xlabel("pénétration  δ  [µm]")
    ax[1].set_ylabel("force axiale  [kN]")
    ax[1].set_title("(b) F–p, branche de charge seule", fontsize=10)

    ax[2].set_xlabel("temps [µs]")
    ax[2].set_ylabel("poste « éléments »  [mJ]")
    ax[2].set_title("(c) énergie prélevée par les éléments", fontsize=10)
    ax[2].text(0.045, 0.94,
               "dissipation visqueuse en fin de test\n"
               r"  $\xi = 2$ :  9,29 mJ  (5,8 % du poste)" "\n"
               "  Cundall : 0 mJ dans les deux",
               transform=ax[2].transAxes, va="top", fontsize=7.8,
               bbox=dict(fc="white", ec="0.82", lw=0.6, pad=3.5))

    for p in ax:
        p.grid(lw=0.4, color="0.93")
        p.set_axisbelow(True)
        for s in ("top", "right"):
            p.spines[s].set_visible(False)

    fig.suptitle("Viscosité de Yan : le choix de $\\xi$ ne change pas la mécanique à ce stade, "
                 "il coûte un facteur 2,37 au chronomètre", fontsize=10.5, y=1.03)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_xi_compare." + ext), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)

    # ecart quantifie, imprime pour ne pas se fier a l oeil
    a, b = load(RUNS[0][0]), load(RUNS[1][0])
    Fa, Fb = a["grpFz"][-1] / 1e3, b["grpFz"][-1] / 1e3
    print("ecrit : fig_xi_compare.pdf / .png")
    print("  force finale : %.4f kN (xi=0,5) contre %.4f kN (xi=2) -> ecart %.2f %%"
          % (Fa, Fb, 100.0 * abs(Fb - Fa) / abs(Fa)))
    print("  duree        : %d s contre %d s -> facteur %.2f"
          % (RUNS[0][4], RUNS[1][4], RUNS[1][4] / RUNS[0][4]))


if __name__ == "__main__":
    main()
