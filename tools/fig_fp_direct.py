#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_fp_direct.py — force - penetration avec la REACTION MESUREE de la roche
# (colonne Fc_rock_insert_z de contactForcePairs, S2 de la campagne du 13/09),
# a cote de l estimateur historique m dv/dt du train supposé RIGIDE.
#
#   python tools/fig_fp_direct.py out_A[:label] ... [--stem results/fig/fpd]
#
# Motif (relecture independante du 13/09, ECARTS §8) : « fig_fp.py assimile la
# vitesse de tout le train a celle du bit ; c est une approximation pendant la
# propagation des ondes. La nouvelle sortie de force insert-roche doit devenir
# la mesure principale, controlee par l impulsion. » Les deux estimateurs sont
# traces ensemble, avec leurs IMPULSIONS : l impulsion est la grandeur qui doit
# coincider (elle integre les oscillations d onde), le maximum non.
# ---------------------------------------------------------------------------
import argparse
import csv
import os
import re

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


def masses(run):
    name = os.path.basename(run.rstrip("/\\"))
    cands = [os.path.join("results", name[4:] + ".log") if name.startswith("out_") else None,
             os.path.join("results", name + ".log")]
    m = {}
    for c in cands:
        if not c or not os.path.isfile(c):
            continue
        for mm in re.finditer(r"corps '(\w+)':.*?masse = ([\d.eE+-]+) kg",
                              open(c, errors="replace").read()):
            m.setdefault(mm.group(1), float(mm.group(2)))
        if "piston" in m:
            break
    return (m.get("piston", 1.057),
            m.get("bit", 1.288) + m.get("insert", 0.0646) + m.get("circlip", 0.0142))


def load(run):
    with open(os.path.join(run, "history.csv"), newline="") as f:
        r = [x for x in csv.DictReader(f) if x]
    col = {k: np.array([float(x[k]) for x in r]) for k in r[0]}
    return col


def smooth(y, t, tau):
    dt = np.median(np.diff(t))
    n = max(1, int(round(tau / dt)))
    if n <= 1:
        return y
    k = np.ones(n) / n
    return np.convolve(y, k, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--stem", default="results/fig/fp_direct")
    ap.add_argument("--tau", type=float, default=2e-6)
    ap.add_argument("--pair", default="rock_insert",
                    help="couple de contactForcePairs porteur de la reaction")
    a = ap.parse_args()

    fig, AX = plt.subplots(2, 2, figsize=(13.5, 9.2))
    cols = ["#1f4e79", "#b22222", "#2e7d32", "#e08a00", "#6a1b9a"]
    for i, spec in enumerate(a.runs):
        run, _, lab = spec.partition(":")
        lab = lab or os.path.basename(run)
        c = load(run)
        t = c["t"]
        tu = t * 1e6
        p = (c["z_insert"][0] - c["z_insert"]) * 1e3
        mp, mt = masses(run)
        dP = mp * np.gradient(-c["vz_piston"], t) + mt * np.gradient(-c["vz_bit"], t)
        Fest = -smooth(smooth(dP, t, a.tau), t, a.tau) * 1e-3
        key = "Fc_" + a.pair + "_z"
        col = cols[i % len(cols)]
        AX[0, 0].plot(p, Fest, color=col, lw=1.0, ls="--",
                      label="%s — estimateur $m\\,\\dot v$ (train rigide)" % lab)
        AX[1, 0].plot(tu, Fest, color=col, lw=1.0, ls="--")
        txt = "%-28s p %.3f mm  v_bit max %.2f m/s" % (lab, p.max(), (-c["vz_bit"]).max())
        if key in c:
            Fr = c[key] * 1e-3
            I1 = np.trapezoid(Fr * 1e3, t)
            I2 = np.trapezoid(Fest * 1e3, t)
            AX[0, 0].plot(p, Fr, color=col, lw=1.6,
                          label="%s — **réaction MESURÉE** (%s)" % (lab, a.pair))
            AX[1, 0].plot(tu, Fr, color=col, lw=1.6)
            AX[0, 1].plot(tu, np.array([np.trapezoid(Fr[:k + 1] * 1e3, t[:k + 1])
                                        for k in range(len(t))]), color=col, lw=1.6,
                          label="%s — mesurée" % lab)
            AX[0, 1].plot(tu, np.array([np.trapezoid(Fest[:k + 1] * 1e3, t[:k + 1])
                                        for k in range(len(t))]), color=col, lw=1.0, ls="--",
                          label="%s — estimateur" % lab)
            txt += ("  |  F max mesurée %.1f kN contre %.1f kN estimée (+%.0f %%)"
                    "  |  impulsion %.3f contre %.3f N.s (%+.1f %%)"
                    % (Fr.max(), Fest.max(), 100 * (Fest.max() - Fr.max()) / max(Fr.max(), 1e-9),
                       I1, I2, 100 * (I2 - I1) / abs(I1) if I1 else 0.0))
        for extra, sty in (("Fc_plate_bit_z", ":"), ("Fc_piston_bit_z", "-.")):
            if extra in c and np.abs(c[extra]).max() > 1.0:
                AX[1, 1].plot(tu, c[extra] * 1e-3, color=col, lw=0.9, ls=sty,
                              label="%s — %s" % (lab, extra.replace("Fc_", "").replace("_z", "")))
        print(txt)

    AX[0, 0].set_xlabel("pénétration $p$ [mm]"); AX[0, 0].set_ylabel("$F$ [kN]")
    AX[0, 0].set_title("(a)  Force – pénétration : réaction mesurée contre estimateur",
                       loc="left", fontsize=11)
    AX[0, 1].set_xlabel(r"temps [$\mu$s]"); AX[0, 1].set_ylabel("impulsion [N·s]")
    AX[0, 1].set_title("(b)  Impulsion cumulée — c'est ELLE qui doit coïncider",
                       loc="left", fontsize=11)
    AX[1, 0].set_xlabel(r"temps [$\mu$s]"); AX[1, 0].set_ylabel("$F$ [kN]")
    AX[1, 0].set_title("(c)  Forces contre le temps", loc="left", fontsize=11)
    AX[1, 1].set_xlabel(r"temps [$\mu$s]"); AX[1, 1].set_ylabel("$F_z$ [kN]")
    AX[1, 1].set_title("(d)  Autres couples de contact (plaque/bit, piston/bit)",
                       loc="left", fontsize=11)
    for A in AX.flat:
        A.axhline(0, color="k", lw=0.5)
        h, l = A.get_legend_handles_labels()
        if h:
            A.legend(frameon=False, fontsize=7.5)
    fig.suptitle("Réaction de la roche MESURÉE dans le contact (`contactForcePairs`) "
                 "contre l'estimateur du train rigide", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=160)
    print("ecrit : %s.pdf / .png" % a.stem)


if __name__ == "__main__":
    main()
