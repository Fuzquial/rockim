#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_fp.py — courbes FORCE - PENETRATION (F-p) et F(t), p(t) d'un impact 3D
# a outil MAILLE (toolShape = none : toolFz vaut 0, la force s'estime).
#
#   python tools/fig_fp.py out_A[:label] out_B[:label] ... [--stem results/fig/fp]
#
# Deux estimations de la force sur la roche, tracees toutes les deux :
#   (1) jauge : F_g = sigma_zz(mi-bit) x A_bit, A = pi (15 mm)^2 = 7,07 cm^2 —
#       c'est la force transmise par l'onde a mi-hauteur du bit (leur fig. 9a),
#       retardee et filtree par le trajet dans le bit ;
#   (2) quantite de mouvement : F_r = m_p dv_p/dt + M_train dv_bit/dt (piston +
#       bit + insert + circlip), lissee sur ~2 us — la roche est la seule force
#       exterieure au systeme piston + train (la poussee du piston est interne),
#       donc F_r est la REACTION DE LA ROCHE, positive vers le haut. Les masses
#       sont lues dans le journal du run (lignes « corps '...' : masse ») : elles
#       different d un maillage a l autre (facettisation des cylindres : piston
#       1,057 kg a s = 1, 0,777 kg a s = 2,5).
# Penetration p = z_insert(0) - z_insert(t) (l'insert part POSE sur la roche).
# ---------------------------------------------------------------------------
import argparse
import csv
import os
import sys

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

A_BIT = np.pi * 0.015 ** 2            # section du bit, m^2 (Phi 30)
M_TRAIN = 1.0802 + 0.0631335 + 0.0142247   # bit + insert + circlip du maillage s = 2,5, kg


def masses(run):
    """masses des corps lues dans results/<run>.log (ou le smoke s = 1)"""
    name = os.path.basename(run.rstrip("/\\"))
    cands = [os.path.join("results", name[4:] + ".log") if name.startswith("out_") else None,
             os.path.join("results", "yang_v3_smoke.log")]
    m = {}
    import re
    for c in cands:
        if not c or not os.path.isfile(c):
            continue
        for mm in re.finditer(r"corps '(\w+)':.*?masse = ([\d.eE+-]+) kg", open(c, errors="replace").read()):
            m.setdefault(mm.group(1), float(mm.group(2)))
        if "piston" in m and "bit" in m:
            break
    mp = m.get("piston", 0.7767)
    mt = m.get("bit", 1.0802) + m.get("insert", 0.0631) + m.get("circlip", 0.0142)
    return mp, mt


def load(run):
    with open(os.path.join(run, "history.csv"), newline="") as f:
        r = [x for x in csv.DictReader(f) if x]
    t = np.array([float(x["t"]) for x in r])
    vz = np.array([float(x["vz_bit"]) for x in r])
    vp = np.array([float(x["vz_piston"]) for x in r])
    zi = np.array([float(x["z_insert"]) for x in r])
    szz = np.array([float(x["szz_bit"]) for x in r])
    return t, vz, vp, zi, szz


def smooth(y, t, tau):
    """moyenne glissante de largeur tau (s) sur une grille non uniforme"""
    if len(t) < 5:
        return y
    dt = np.median(np.diff(t))
    n = max(1, int(round(tau / dt)))
    if n <= 1:
        return y
    k = np.ones(n) / n
    return np.convolve(y, k, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--stem", default="results/fig/fp")
    ap.add_argument("--tau", type=float, default=2e-6, help="lissage de F_d, s")
    ap.add_argument("--mass", type=float, default=M_TRAIN, help="masse du train de frappe, kg")
    a = ap.parse_args()

    fig, AX = plt.subplots(2, 2, figsize=(13.5, 9.2))
    cols = ["#1f4e79", "#b22222", "#2e7d32", "#e08a00", "#6a1b9a"]
    for i, spec in enumerate(a.runs):
        run, _, lab = spec.partition(":")
        lab = lab or os.path.basename(run)
        t, vz, vp, zi, szz = load(run)
        mp, mt = masses(run)
        p = (zi[0] - zi) * 1e3                         # mm
        Fg = -szz * A_BIT * 1e-3                        # kN, compression > 0
        # reaction de la roche (vers le haut > 0) : derivee de la quantite de
        # mouvement totale piston + train, lissee deux fois
        dP = mp * np.gradient(vp, t) + mt * np.gradient(vz, t)
        Fr = smooth(smooth(dP, t, a.tau), t, a.tau) * 1e-3   # kN
        c = cols[i % len(cols)]
        tu = t * 1e6
        lab2 = "%s  (piston %.3f kg = %.1f J, train %.3f kg)" % (lab, mp, 0.5 * mp * vp[0] ** 2, mt)
        AX[0, 0].plot(p, Fr, color=c, lw=1.4, label=lab2)
        AX[0, 1].plot(p, Fg, color=c, lw=1.2, label=lab)
        AX[1, 0].plot(tu, Fr, color=c, lw=1.4, label=lab + " (roche, q. de mvt)")
        AX[1, 0].plot(tu, Fg, color=c, lw=0.9, ls="--", label=lab + " (jauge mi-bit)")
        AX[1, 1].plot(tu, p, color=c, lw=1.4, label=lab)
        print("%-34s t_fin %.0f us  p_max %.3f mm  F_roche max %.1f kN  F_jauge max %.1f kN  piston %.3f kg (%.1f J) train %.3f kg"
              % (lab, tu[-1], p.max(), Fr.max(), Fg.max(), mp, 0.5 * mp * vp[0] ** 2, mt))

    AX[0, 0].set_xlabel("pénétration $p$ [mm]"); AX[0, 0].set_ylabel("$F_r = \\dot P_{piston+train}$ [kN]")
    AX[0, 0].set_title("(a)  Réaction de la roche (quantité de mouvement) contre pénétration", loc="left", fontsize=11)
    AX[0, 1].set_xlabel("pénétration $p$ [mm]"); AX[0, 1].set_ylabel("$F_g = \\sigma_{zz}\\,A_{bit}$ [kN]")
    AX[0, 1].set_title("(b)  Force de jauge (mi-bit, onde) contre pénétration", loc="left", fontsize=11)
    AX[1, 0].set_xlabel(r"temps [$\mu$s]"); AX[1, 0].set_ylabel("$F$ [kN]")
    AX[1, 0].set_title("(c)  Forces contre le temps", loc="left", fontsize=11)
    AX[1, 1].set_xlabel(r"temps [$\mu$s]"); AX[1, 1].set_ylabel("$p$ [mm]")
    AX[1, 1].set_title("(d)  Pénétration de l'insert", loc="left", fontsize=11)
    for A in AX.flat:
        A.axhline(0, color="k", lw=0.5)
        A.legend(frameon=False, fontsize=8)
    fig.suptitle("Force – pénétration, insert unique sur Kuru Grey, piston 9 m/s  "
                 "(A bit = %.2f cm², lissage %.0f µs ; masses lues dans les journaux)"
                 % (A_BIT * 1e4, a.tau * 1e6), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=160)
    print("ecrit : %s.pdf / .png" % a.stem)


if __name__ == "__main__":
    main()
