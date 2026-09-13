#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_retournement.py — DATER le point de retour d'un impact dont le calcul
# n'est pas termine, sans extrapoler n'importe quoi.
#
#   python tools/fig_retournement.py out_A[:label] ... [--stem results/fig/ret]
#                                    [--periode 20] [--depuis 170] [--corps insert]
#
# Pourquoi ce script. La vitesse instantanee de l'insert OSCILLE : le train de
# frappe renvoie des ondes, et sur ce cas la periode vaut une vingtaine de
# microsecondes pour une amplitude de plus ou moins 2 m/s. Extrapoler le
# passage a zero sur les derniers pas donne donc un chiffre qui n'a aucun sens
# (le 13/09, trois points donnaient 212 us alors que la moyenne en donne 260).
#
# Ce qu'on trace, en haut, la vitesse instantanee et sa moyenne glissante sur
# une periode exactement ; en bas, la moyenne seule avec la droite ajustee
# depuis --depuis et son intersection avec zero.
#
# Ce que le chiffre vaut. C'est une BORNE BASSE : la deceleration faiblit a
# l'approche du retournement, quand la force de contact retombe. Le lire comme
# « pas avant », jamais comme « a ».
#
# La periode n'est pas devinee : elle est mesuree sur le signal (autocorrelation
# du residu apres retrait d'une tendance lineaire), et --periode ne sert qu'a
# forcer une valeur si la mesure echoue. La valeur retenue est imprimee.
# ---------------------------------------------------------------------------
import argparse
import csv
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

COL = {"insert": "vz_insert", "bit": "vz_bit", "piston": "vz_piston"}


def load(run, corps):
    chemin = run.split(":")[0]
    label = run.split(":", 1)[1] if ":" in run else os.path.basename(chemin)
    t, v = [], []
    with open(os.path.join(chemin, "history.csv")) as fh:
        r = csv.DictReader(fh)
        for row in r:
            t.append(float(row["t"]) * 1e6)
            v.append(float(row[COL[corps]]))
    return label, np.asarray(t), np.asarray(v)


def periode_mesuree(t, v, tmin):
    """Periode dominante du residu, par autocorrelation. None si indecidable."""
    m = t >= tmin
    tt, vv = t[m], v[m]
    if tt.size < 200:
        return None
    a, b = np.polyfit(tt, vv, 1)
    res = vv - (a * tt + b)
    res -= res.mean()
    ac = np.correlate(res, res, mode="full")[res.size - 1:]
    if ac[0] <= 0:
        return None
    ac = ac / ac[0]
    dt = float(np.median(np.diff(tt)))
    # premier maximum local franc apres le premier passage sous zero
    i = 1
    while i < ac.size and ac[i] > 0.0:
        i += 1
    meilleur, val = None, 0.15
    while i < ac.size - 1:
        if ac[i] > ac[i - 1] and ac[i] >= ac[i + 1] and ac[i] > val:
            meilleur, val = i, ac[i]
        i += 1
    return None if meilleur is None else meilleur * dt


def glissante(t, v, largeur):
    """Moyenne sur une fenetre centree de `largeur` us, bords exclus."""
    tc, vm = [], []
    demi = largeur / 2.0
    for i, c in enumerate(t):
        if c - demi < t[0] or c + demi > t[-1]:
            continue
        m = (t >= c - demi) & (t <= c + demi)
        tc.append(c)
        vm.append(v[m].mean())
    return np.asarray(tc), np.asarray(vm)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("runs", nargs="+")
    p.add_argument("--corps", default="insert", choices=sorted(COL))
    p.add_argument("--periode", type=float, default=None,
                   help="force la periode en us au lieu de la mesurer")
    p.add_argument("--depuis", type=float, default=170.0,
                   help="debut de la fenetre d'ajustement, en us")
    p.add_argument("--stem", default="results/fig/retournement")
    a = p.parse_args()

    fig, (h, b) = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
    for k, run in enumerate(a.runs):
        label, t, v = load(run, a.corps)
        T = a.periode or periode_mesuree(t, v, a.depuis)
        source = "imposee" if a.periode else "mesuree"
        if T is None:
            T, source = 20.0, "defaut"
        tc, vm = glissante(t, v, T)
        c = "C%d" % k

        h.plot(t, v, color=c, lw=0.6, alpha=0.45)
        h.plot(tc, vm, color=c, lw=1.9,
               label="%s  (periode %s : %.1f us)" % (label, source, T))
        b.plot(tc, vm, color=c, lw=1.9, label=label)

        m = tc >= a.depuis
        if m.sum() >= 3:
            pente, ord0 = np.polyfit(tc[m], vm[m], 1)
            if pente > 0:
                t0 = -ord0 / pente
                xs = np.linspace(a.depuis, max(t0 * 1.02, tc[-1]), 50)
                b.plot(xs, pente * xs + ord0, color=c, ls="--", lw=1.1)
                b.plot([t0], [0.0], marker="o", ms=7, mfc="none", color=c)
                b.annotate(u"retournement\nau plus tôt %.0f µs" % t0,
                           xy=(t0, 0.0), xytext=(t0 - 4, -2.6),
                           ha="right", va="top", color=c, fontsize=9,
                           arrowprops=dict(arrowstyle="->", color=c, lw=0.8))
                b.set_xlim(right=t0 * 1.06)
                print("%-24s periode %s %.1f us | pente %+.4f (m/s)/us "
                      "| zero a %.0f us (BORNE BASSE)"
                      % (label, source, T, pente, t0))
            else:
                print("%-24s pente non decroissante : rien a extrapoler" % label)

    for ax in (h, b):
        ax.axhline(0.0, color="0.4", lw=0.8)
        ax.grid(alpha=0.3)
    h.set_ylabel(r"$v_z$ de l'%s  [m/s]" % a.corps)
    h.set_title("Le retournement ne se lit pas sur la vitesse instantanee")
    h.legend(fontsize=8.5, loc="lower right")
    b.set_ylabel(r"$v_z$ moyennee sur une periode  [m/s]")
    b.set_xlabel(r"temps  [$\mu$s]")
    b.set_title("Sur la moyenne, il se date -- et le chiffre est une borne basse")
    fig.tight_layout()

    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)
    for e in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, e), dpi=170, bbox_inches="tight")
    print("ecrit : %s.pdf / .png" % a.stem)


main()
