#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_lambda.py — la synthese du balayage LAMBDA (Wang et al. 2024,
# fig. 15f et 16a) : demi-axes de l'EDZ et comptes de fissures contre le
# coefficient de pression laterale.
#
#   python tunnel_edz/tools/fig_lambda.py --stem tunnel_edz/fig_lambda
#
# Entrees : les JSON d'edz_metrics.py, un par cas (voir RUNS ci-dessous).
# Les valeurs publiees sont LUES SUR LEURS FIGURES (±1 m / ±10 %) et
# tracees en gris pour situer la tendance, pas pour un fit.
# ---------------------------------------------------------------------------
import argparse
import io
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})

RUNS = [(0.50, "out_tunnel_lam0p5"), (0.75, "out_tunnel_lam0p75"),
        (1.00, "out_tun_ref_iso"), (1.25, "out_tunnel_lam1p25"),
        (1.50, "out_tunnel_lam1p5")]

# leur fig. 15f, lue a la regle (la numerisation vaut +-1 m)
PUB_LAM = [0.50, 0.75, 1.00, 1.25, 1.50]
PUB_H = [24.3, 21.0, 18.0, 15.0, 12.8]
PUB_V = [16.5, 17.8, 18.0, 15.2, 13.8]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", default="tunnel_edz/fig_lambda")
    a = ap.parse_args()

    lam, hx, hy, ncr = [], [], [], []
    for L, run in RUNS:
        p = os.path.join(run, "edz_metrics.json")
        if not os.path.exists(p):
            print("  (absent : %s — cas saute)" % p)
            continue
        d = json.load(io.open(p))
        lam.append(L)
        hx.append(2 * d["edz_halfaxis_x_p95_m"])       # axes COMPLETS, comme eux
        hy.append(2 * d["edz_halfaxis_y_p95_m"])
        ncr.append(int(d["broken"]))
        print("  lambda %.2f : axe h %.1f m, axe v %.1f m, %s fissures"
              % (L, hx[-1], hy[-1], ncr[-1]))
    if len(lam) < 2:
        print("  moins de 2 cas : rien a tracer")
        return

    fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.9))
    A = ax[0]
    A.plot(PUB_LAM, PUB_H, "s--", color="#aaa", ms=4,
           label="axe horizontal — leur fig. 15f (lu)")
    A.plot(PUB_LAM, PUB_V, "o--", color="#ccc", ms=4,
           label="axe vertical — leur fig. 15f (lu)")
    A.plot(lam, hx, "s-", color="#b22222", ms=5, label="axe horizontal — rockim")
    A.plot(lam, hy, "o-", color="#1f4e79", ms=5, label="axe vertical — rockim")
    A.set_xlabel(r"coefficient de pression latérale $\lambda$")
    A.set_ylabel("axe de l'EDZ  [m]")
    A.set_title("(a)  Ellipse de l'EDZ (leur fig. 15f)", loc="left",
                fontsize=11)
    A.legend(frameon=False, fontsize=8)

    B = ax[1]
    B.plot(lam, ncr, "ks-", ms=5, label="rockim")
    B.set_xlabel(r"$\lambda$")
    B.set_ylabel("fissures (total)")
    B.set_title("(b)  Nombre de fissures (leur fig. 16a : décroissant)",
                loc="left", fontsize=11)
    B.legend(frameon=False, fontsize=9)
    fig.suptitle("Tunnel de Hutou Beishan — balayage du coefficient de "
                 "pression latérale ($\\sigma_h$ = 5 MPa, "
                 "$\\sigma_v = 5/\\lambda$)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=165)
    print("écrit : %s.pdf et .png" % a.stem)


if __name__ == "__main__":
    main()
