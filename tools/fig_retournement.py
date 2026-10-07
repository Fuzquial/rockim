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
    """Retourne label, t[us], vz[m/s] et la force de contact roche/insert [kN].

    La force est le CRITERE qui distingue une pause d'un retournement : sous une
    pause elle reste plate, a un vrai retournement elle retombe. Sans elle on
    date des croisements de vitesse qui n'en sont pas (erreur du 14/09).
    """
    chemin = run.split(":")[0]
    label = run.split(":", 1)[1] if ":" in run else os.path.basename(chemin)
    t, v, F = [], [], []
    with open(os.path.join(chemin, "history.csv")) as fh:
        r = csv.DictReader(fh)
        for row in r:
            t.append(float(row["t"]) * 1e6)
            v.append(float(row[COL[corps]]))
            F.append(float(row.get("Fc_rock_insert_z", "nan")) / 1000.0)
    return label, np.asarray(t), np.asarray(v), np.asarray(F)


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
    # Une periode doit tenir plusieurs fois dans la fenetre, sinon ce n'est pas
    # une oscillation mais la forme du signal : apres le rebond l'autocorrelation
    # renvoyait 106 us sur une fenetre de 120, ce qui vidait la moyenne de sens.
    imax = int((tt[-1] - tt[0]) / 3.0 / dt)
    meilleur, val = None, 0.15
    while i < min(ac.size - 1, imax):
        if ac[i] > ac[i - 1] and ac[i] >= ac[i + 1] and ac[i] > val:
            meilleur, val = i, ac[i]
        i += 1
    return None if meilleur is None else meilleur * dt


def amplitude(t, v, centre, largeur=10.0):
    """Ecart-type du residu autour de la droite locale = ce qui reste d'oscillation."""
    m = (t >= centre - largeur / 2) & (t <= centre + largeur / 2)
    if m.sum() < 20:
        return None
    a, b = np.polyfit(t[m], v[m], 1)
    return float((v[m] - (a * t[m] + b)).std())


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
    eteinte = False
    mesure = None
    for k, run in enumerate(a.runs):
        label, t, v, Fc = load(run, a.corps)
        if k == 0 and np.isfinite(Fc).any():
            hF = h.twinx()
            hF.plot(t, Fc, color="C2", lw=1.0, alpha=.55)
            hF.set_ylabel(u"force de contact roche/insert  [kN]", color="C2",
                          fontsize=9)
            hF.tick_params(axis="y", labelcolor="C2", labelsize=8)
            hF.set_ylim(0, float(np.nanmax(Fc)) * 2.6)
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

        # --- si le retournement a EU LIEU, on le lit, on ne l'extrapole plus ----
        # --- un croisement de zero n'est PAS un retournement --------------------
        # Erreur commise le 14/09 : la vitesse a croise zero a 260,1 us, j'ai
        # annonce le retournement et un accord a 2 % avec Yang. Elle est
        # redescendue a -1,9 m/s et l'enfoncement a repris. C'etait une PAUSE
        # sous charge maintenue (la force de contact n'est jamais retombee).
        # Un retournement se verifie : apres lui la vitesse reste positive et
        # l'enfoncement cesse de croitre. Tant que la fin de l'enregistrement
        # ne le confirme pas, on marque le croisement comme provisoire.
        # on parcourt TOUS les croisements montants, on retient le premier confirme
        mont = np.where((v[:-1] <= 0.0) & (v[1:] > 0.0) & (t[1:] > a.depuis))[0]
        crois = mont
        for j in mont:
            t0 = (t[j] + (0.0 - v[j]) * (t[j + 1] - t[j]) / (v[j + 1] - v[j]))
            apres = t >= t0
            frac_neg = float((v[apres] < 0.0).mean()) if apres.sum() else 1.0
            if frac_neg < 0.15:
                mesure = t0
                h.axvline(t0, color="0.15", lw=1.1)
                h.plot([t0], [0.0], marker="*", ms=13, color="0.15")
                h.annotate(u"retournement CONFIRMÉ\n%.1f µs" % t0,
                           xy=(t0, 0.0), xytext=(t0 - 8, -3.0),
                           ha="right", va="top", color="0.15", fontsize=9.5,
                           arrowprops=dict(arrowstyle="->", color="0.15", lw=0.9))
                print("%-24s RETOURNEMENT CONFIRME a %.2f us "
                      "(la vitesse reste positive ensuite)" % (label, t0))
                break
            h.axvline(t0, color="C3", lw=0.9, ls="--")
            h.plot([t0], [0.0], marker="x", ms=9, color="C3", mew=2)
            h.annotate(u"pause %.0f µs\n(la force ne retombe pas)" % t0,
                       xy=(t0, 0.0), xytext=(t0 - 6, 1.9),
                       ha="right", va="bottom", color="C3", fontsize=8.5,
                       arrowprops=dict(arrowstyle="->", color="C3", lw=0.8))
            print("%-24s croisement a %.2f us : PAUSE (%.0f %% du temps suivant "
                  "vz < 0), pas un retournement" % (label, t0, 100 * frac_neg))
        if mesure is None and mont.size:
            h.annotate(u"%d pause(s) sous charge,\naucun retournement confirmé"
                       % mont.size,
                       xy=(t[mont[0] + 1], 0.0), xytext=(t[mont[0] + 1] - 8, -3.0),
                       ha="right", va="top", color="C3", fontsize=9,
                       arrowprops=dict(arrowstyle="->", color="C3", lw=0.9))
        # l'enfoncement maximal, lui, ne se discute pas
        zmax = None
        try:
            import csv as _c
            zz, tt2 = [], []
            for r in _c.DictReader(open(os.path.join(a.runs[k].split(":")[0],
                                                     "history.csv"))):
                tt2.append(float(r["t"]) * 1e6)
                zz.append(float(r["z_insert"]) * 1e3)
            zz, tt2 = np.asarray(zz), np.asarray(tt2)
            i = int(np.argmin(zz[tt2 > 100.0]))
            zmax = float(tt2[tt2 > 100.0][i])
            print("%-24s enfoncement MAXIMAL a %.1f us (%.4f mm)"
                  % (label, zmax, zz[tt2 < 40].mean() - zz[tt2 > 100.0][i]))
        except Exception:                                        # noqa: BLE001
            pass

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

        # --- quand l'oscillation est morte, la moyenne glissante MENT -----------
        # Elle est centree : tant que sa fenetre contient encore la partie
        # oscillante, plus negative, elle traine derriere et repousse le zero.
        # On compare donc l'amplitude residuelle du debut de la fenetre
        # d'ajustement a celle de la fin ; si elle s'est effondree, on trace
        # aussi l'extrapolation DIRECTE de la vitesse instantanee, qui est
        # alors le bon estimateur.
        a0 = amplitude(t, v, a.depuis + 10.0)
        a1 = amplitude(t, v, t[-1] - 5.0)
        if a0 and a1 and a1 < 0.30 * a0:
            eteinte = True
            for w, style in ((15.0, "-"), (6.0, ":")):
                m = t >= t[-1] - w
                if m.sum() < 10:
                    continue
                pe, o0 = np.polyfit(t[m], v[m], 1)
                if pe <= 0:
                    continue
                t0 = -o0 / pe
                xs = np.linspace(t[-1] - w, t0 * 1.01, 30)
                h.plot(xs, pe * xs + o0, color="0.25", ls=style, lw=1.2)
                if style == "-":
                    h.plot([t0], [0.0], marker="s", ms=6, mfc="none", color="0.25")
                    h.annotate(u"oscillation éteinte :\nzéro direct à %.0f µs" % t0,
                               xy=(t0, 0.0), xytext=(t0 - 6, -3.4),
                               ha="right", va="top", color="0.25", fontsize=9,
                               arrowprops=dict(arrowstyle="->", color="0.25", lw=0.8))
                print("%-24s amplitude residuelle %.3f -> %.3f m/s : OSCILLATION "
                      "ETEINTE, extrapolation directe sur %.0f us -> zero a %.1f us"
                      % (label, a0, a1, w, t0))
            print("%-24s => la moyenne glissante est alors BIAISEE TARDIVE, "
                  "c'est le zero direct qu'il faut lire" % label)

    for ax in (h, b):
        ax.axhline(0.0, color="0.4", lw=0.8)
        ax.grid(alpha=0.3)
    h.set_ylabel(r"$v_z$ de l'%s  [m/s]" % a.corps)
    h.legend(fontsize=8.5, loc="lower right")
    b.set_ylabel(r"$v_z$ moyennee sur une periode  [m/s]")
    b.set_xlabel(r"temps  [$\mu$s]")
    if mesure is not None:
        h.set_title(u"Le retournement a eu lieu : il est mesuré, plus extrapolé")
        b.set_title(u"Ce que la moyenne glissante prédisait, à comparer")
    elif 'crois' in dir() and crois.size:
        h.set_title(u"La vitesse a croisé zéro, puis est repartie vers le bas")
        b.set_title(u"Le retournement n'a pas eu lieu : rien à dater encore")
    elif eteinte:
        h.set_title(u"L'oscillation s'est éteinte : c'est ICI que le zéro se lit")
        b.set_title(u"La moyenne, elle, traîne derrière et repousse le zéro à tort")
    else:
        h.set_title(u"Tant que ça oscille, le zéro ne se lit pas ici")
        b.set_title(u"Sur la moyenne il se date, et le chiffre est une borne basse")
    fig.tight_layout()

    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)
    for e in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, e), dpi=170, bbox_inches="tight")
    print("ecrit : %s.pdf / .png" % a.stem)


main()
