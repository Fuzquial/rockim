#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_cible_dumoulin.py - le vecteur de CIBLES de calibration Red Bohus,
# lu dans Dumoulin et al., Geomech. Energy Environ. 40 (2024) 100592.
#
#   python tools/fig_cible_dumoulin.py
#
# Ce que la figure montre, et pourquoi.
#
# Panneau (a) : l'ENVELOPPE. Les 12 essais triaxiaux de la Table 6 (sigma_3 =
# 20 / 50 / 75 / 100 MPa, 3 repetitions chacun, dispersion 2-6 MPa soit moins
# de 1 %) et les 4 UCS de la Table 1 (126 +- 19 MPa, 41 % d'etendue). Deux
# ajustements y sont portes :
#   - Mohr-Hoek  q = B sigma_3^n + sigma_c  a sigma_c LIBRE : il passe a 4,5
#     MPa RMS mais extrapole sigma_c = -15 MPa. Les triaxiaux SEULS ne
#     contraignent pas l'UCS ; c'est exactement pourquoi l'article prend 20 MPa
#     et non 0 comme etat de reference (Table 9).
#   - Mohr-Hoek a sigma_c IMPOSE a l'UCS mesure : B = 56,7, n = 0,538, RMS
#     6,2 MPa. L'exposant retombe sur le alpha = 0,54 de l'eq. 5 de l'article,
#     ajuste independamment sur les trois granites.
# La droite de Mohr-Coulomb tangente a 20 MPa est tracee pour montrer ce qu'un
# frottement CONSTANT donne : elle sur-predit de 100 MPa a 100 MPa de
# confinement. La concavite n'est pas un detail de forme.
#
# Panneau (b) : la PENTE LOCALE dq/dsigma_3, qui est la lecture utile pour un
# calage. Elle vaut 14,0 entre 0 et 20 MPa, 6,5 entre 20 et 50, 3,8 entre 75 et
# 100 : un facteur 3,7 d'un bout a l'autre. Aucun jeu (c, phi) constant ne
# couvre cela, et c'est la mesure qui justifie d'activer une enveloppe
# meridienne en puissance (law = dpr + meridian = power + mhForm = principal)
# plutot que d'ajuster phi confinement par confinement.
#
# Les points de simulation disponibles sont superposes quand ils existent.
# ---------------------------------------------------------------------------
import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --- Table 6 (triaxiaux) et Table 1 (UCS), Red Bohus ----------------------
S3 = np.array([20, 20, 20, 50, 50, 50, 75, 75, 75, 100, 100, 100], float)
Q = np.array([402.4, 407.9, 404.3, 602.0, 597.9, 597.1,
              702.3, 697.3, 712.7, 794.7, 800.0, 803.1])
UCS = np.array([112.1, 158.2, 117.9, 114.8])


def peak(run):
    """Pic FILTRE (mediane glissante) du deviateur d'un run rockim.

    Le pic BRUT est inexploitable des que beaucoup de joints s'inserent : chaque
    insertion est un transitoire de contrainte d'un pas, et sur les cas a faible
    cohesion le maximum de sigma est un de ces pics isoles - mesure du
    2026-09-11 : 405,7 MPa de pic brut pour un plateau reel a 233,7.
    """
    h = os.path.join(HERE, run, "history.csv")
    if not os.path.exists(h):
        return None
    r = list(csv.DictReader(open(h)))
    t = np.array([float(x["t"]) for x in r])
    s = np.array([float(x["sigma"]) for x in r])
    dly = 0.0
    p = os.path.join(HERE, run, "config_effective.cfg")
    if os.path.exists(p):
        for line in open(p, encoding="utf-8", errors="replace"):
            line = line.split("#")[0]
            if "pullDelay" in line and "=" in line:
                try:
                    dly = float(line.split("=")[1].strip())
                except ValueError:
                    pass
    q = (s - np.interp(dly, t, s)) / 1e6
    k = 25
    qm = np.array([np.median(q[max(0, i - k):i + k + 1]) for i in range(len(q))])
    return float(qm.max())


def main():
    mh = lambda x, B, n, sc: B * np.power(np.maximum(x, 0), n) + sc
    pf, _ = curve_fit(mh, S3, Q, p0=[60, 0.6, 150], maxfev=40000)
    uc = UCS.mean()
    pc, _ = curve_fit(lambda x, B, n: mh(x, B, n, uc), S3, Q,
                      p0=[60, 0.6], maxfev=40000)

    fig, (A, B) = plt.subplots(1, 2, figsize=(12.5, 5.0))
    x = np.linspace(0, 105, 400)

    A.plot(x, mh(x, *pf), "-", color="0.35", lw=1.6,
           label="Mohr-Hoek, $\\sigma_c$ libre : $n$ = %.3f (RMS 4,5 MPa)" % pf[1])
    A.plot(x, mh(x, pc[0], pc[1], uc), "-", color="crimson", lw=2.0,
           label="Mohr-Hoek, $\\sigma_c$ = UCS mesure : $B$ = %.1f, $n$ = %.3f"
                 % (pc[0], pc[1]))
    # tangente de Mohr-Coulomb a 20 MPa : ce que donne un frottement constant
    k20 = pc[0] * pc[1] * 20 ** (pc[1] - 1)
    A.plot(x, mh(20, pc[0], pc[1], uc) + k20 * (x - 20), "--", color="crimson",
           lw=1.2, label="Coulomb tangent a 20 MPa (pente %.1f)" % k20)
    A.plot(S3, Q, "o", ms=7, mfc="w", mec="k", mew=1.4,
           label="Table 6 : 12 triaxiaux")
    A.plot(np.zeros_like(UCS), UCS, "s", ms=7, mfc="w", mec="darkblue", mew=1.4,
           label="Table 1 : 4 UCS (126 $\\pm$ 19)")

    for run, s3, lab in [("out_size_L1_s4211", 0, "sim. c = 25"),
                         ("out_cal20_c25", 20, None),
                         ("out_cal20_c65", 20, "sim. c = 65")]:
        v = peak(run)
        if v is None:
            continue
        A.plot(s3, v, "*", ms=15, color="green", zorder=5, label=lab)

    A.set_xlabel("$\\sigma_3$ (MPa)")
    A.set_ylabel("deviateur au pic $q = \\sigma_1 - \\sigma_3$ (MPa)")
    A.set_title("(a) Enveloppe Red Bohus : la concavite est la mesure")
    A.legend(fontsize=8, loc="lower right")
    A.grid(alpha=0.3)
    A.set_xlim(-4, 105)
    A.set_ylim(0, 900)

    seg = [(0, uc), (20, Q[S3 == 20].mean()), (50, Q[S3 == 50].mean()),
           (75, Q[S3 == 75].mean()), (100, Q[S3 == 100].mean())]
    xm = [(seg[i][0] + seg[i + 1][0]) / 2 for i in range(4)]
    pm = [(seg[i + 1][1] - seg[i][1]) / (seg[i + 1][0] - seg[i][0])
          for i in range(4)]
    B.plot(x, pc[0] * pc[1] * np.power(np.maximum(x, 1e-6), pc[1] - 1), "-",
           color="crimson", lw=2.0, label="Mohr-Hoek ajuste")
    B.plot(xm, pm, "o-", color="k", ms=8, lw=1.4,
           label="pente mesuree, segment par segment")
    for a, b in zip(xm, pm):
        B.annotate("%.1f" % b, (a, b), textcoords="offset points",
                   xytext=(6, 7), fontsize=9)
    B.set_xlabel("$\\sigma_3$ (MPa)")
    B.set_ylabel("$dq/d\\sigma_3$ (-)")
    B.set_title("(b) La pente varie d'un facteur 3,7 : aucun $\\varphi$ constant")
    B.legend(fontsize=9)
    B.grid(alpha=0.3)
    B.set_xlim(-4, 105)
    B.set_ylim(0, 16)

    fig.suptitle("Cibles de calibration Red Bohus, Dumoulin et al. "
                 "Geomech. Energy Environ. 40 (2024) 100592", fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "results", "cible_dumoulin.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    print("wrote " + out)
    print("  Mohr-Hoek sigma_c libre   : B = %.2f  n = %.3f  sigma_c = %.1f"
          % tuple(pf))
    print("  Mohr-Hoek sigma_c = %.1f  : B = %.2f  n = %.3f" % (uc, pc[0], pc[1]))


if __name__ == "__main__":
    main()
