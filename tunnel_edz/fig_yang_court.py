#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_yang_court.py — planche du test court "enveloppe de Yang, machinerie de
# Yan" sur le maillage de production p1_ultra (2026-08-19).
#
#   run : out_yang_court    (impact3d_yang_court.cfg, T = 2e-5 s = T_ref / 6)
#
# CE QUE CETTE PLANCHE N EST PAS : un resultat de percussion. A T/6 l impact
# commence a peine — 5 joints casses sur 247 226, l insert n a perdu que
# 0,27 m/s sur 8. Il n y a NI decharge, NI rebond, NI fissure radiale. Elle
# montre que la mecanique numerique tient, rien d autre. Le titre le dit.
#
#   python tunnel_edz/fig_yang_court.py
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
    "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUN = "out_yang_court"


def main():
    h = np.genfromtxt(os.path.join(ROOT, RUN, "history.csv"), delimiter=",",
                      names=True, invalid_raise=False)
    t = h["t"] * 1e6                                   # us
    F = h["grpFz"] / 1e3                               # kN
    v = -h["grpVz"]                                    # m/s, positif = descend

    fig, ax = plt.subplots(1, 3, figsize=(10.4, 3.4))

    # ---- (a) force et vitesse ------------------------------------------
    p = ax[0]
    p.plot(t, F, color="#0B4F9E", lw=1.4)
    p.set_xlabel("temps [µs]")
    p.set_ylabel("force axiale sur l'insert  [kN]", color="#0B4F9E")
    p.tick_params(axis="y", labelcolor="#0B4F9E")
    q = p.twinx()
    q.plot(t, v, color="#C8342B", lw=1.4, ls="--")
    q.set_ylabel("vitesse de l'insert  [m s$^{-1}$]", color="#C8342B")
    q.tick_params(axis="y", labelcolor="#C8342B")
    p.set_title("(a) charge en cours — 8,00 → 7,73 m s⁻¹", fontsize=9.5)

    # ---- (b) force-penetration, BRANCHE DE CHARGE SEULE -----------------
    p = ax[1]
    i0 = int(np.argmax(np.abs(h["grpFz"]) > 1.0))
    d = (h["grpZ"][i0] - h["grpZ"]) * 1e6               # um
    m = np.arange(len(F)) >= i0
    p.plot(d[m], F[m], color="#0B4F9E", lw=1.6)
    W = np.trapezoid(F[m] * 1e3, d[m] * 1e-6)
    p.set_xlabel("pénétration de l'insert  δ  [µm]")
    p.set_ylabel("force axiale  [kN]")
    p.set_title("(b) F–p : branche de CHARGE seule,\npas de décharge à T/6", fontsize=9.5)
    p.text(0.04, 0.92, "aire = %.3f J\nδ max = %.0f µm" % (W, d[m].max()),
           transform=p.transAxes, va="top", fontsize=8,
           bbox=dict(fc="white", ec="0.8", lw=0.6, pad=3))

    # ---- (c) postes du bilan --------------------------------------------
    p = ax[2]
    for nom, cle, c, ls in (("éléments", "eEl", "#0B4F9E", "-"),
                            ("contact", "eGc", "#8A5A00", "-"),
                            ("joints", "eJnt", "#1B8A3A", "-"),
                            ("frottement", "eFric", "#8A5A00", ":"),
                            ("Cundall", "eCund", "#C8342B", "-")):
        p.plot(t, -h[cle] * 1e3, color=c, ls=ls, lw=1.5, label=nom)
    p.axhline(0.0, color="0.8", lw=0.6)
    p.set_xlabel("temps [µs]")
    p.set_ylabel("énergie prélevée  [mJ]")
    p.legend(loc="upper left", framealpha=0.95, ncol=2)
    p.set_title("(c) Cundall à zéro : la viscosité de Yan\nl'a remplacé", fontsize=9.5)

    for p in ax:
        p.grid(lw=0.4, color="0.92")
        p.set_axisbelow(True)
    fig.suptitle("Test de convergence — enveloppe de Yang (éq. 1) sur la machinerie de Yan, "
                 "maillage p1_ultra, T = T$_{ref}$/6", fontsize=10.5, y=1.03)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_yang_court." + ext), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print("ecrit : fig_yang_court.pdf / .png   (aire de charge %.4f J)" % W)


if __name__ == "__main__":
    main()
