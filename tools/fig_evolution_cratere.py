#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_evolution_cratere.py — les mesures du cratere CONTRE LE TEMPS.
#
#   python tools/fig_evolution_cratere.py out_A [--frames 10-18] [--sectors 16]
#                                         [--stem results/fig/evo] [--cache f.json]
#
# Pourquoi. crater_metrics.py mesure une trame. On lisait donc l'evolution en
# comparant des nombres a la main d'un message a l'autre, ce qui est exactement
# la ou l'on se trompe : le 14/09 j'ai ecrit que le detachement s'etait arrete
# alors que le compteur de fragments n'avance que par paliers de trame.
#
# Ce script appelle crater_metrics sur chaque trame, met les resultats en cache
# (ils ne changent plus une fois la trame ecrite) et trace quatre panneaux :
# rayon de peau et extension totale, profondeur, volume detache, et portee des
# radiales. Les reperes de Yang (rayon ~7 mm, enfoncement ~1,0 mm) sont poses
# en pointilles quand ils existent.
#
# La distinction qui compte : le RAYON DE PEAU est la grandeur comparable a
# l'experience (ce que la camera voit), l'EXTENSION TOTALE inclut les fissures
# qui rayonnent hors du cratere et n'a pas d'equivalent mesure.
# ---------------------------------------------------------------------------
import argparse
import json
import os
import re
import subprocess
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
YANG = dict(rCrater=7.0, depth=1.0)

MOTIFS = {
    "t":        r"frame (\d+), (\d+) joints casses",
    "rCrater":  r"R_crater \(p95 peau\)\s*=\s*([-\d.]+) mm",
    "rMax":     r"R_max \(tous casses\)\s*=\s*([-\d.]+) mm",
    "depth":    r"profondeur\s*=\s*([-\d.]+) mm",
    "aire":     r"aire cassee\s*=\s*([-\d.]+) mm",
    "vDet":     r"volume detache\s*=\s*([-\d.]+) mm",
    "vEnd":     r"volume endommage\s*=\s*([-\d.]+) mm",
    "radMoy":   r"fissures radiales.*portee moy ([-\d.]+) mm",
    "radMax":   r"fissures radiales.*max ([-\d.]+) mm",
    "brasMoy":  r"bras endommages.*portee moy ([-\d.]+) mm",
    "nBroken":  r"(\d+) joints casses",
}


def temps_des_trames(run):
    """Lit frames.csv : indice -> temps en us."""
    d = {}
    with open(os.path.join(run, "frames.csv")) as fh:
        next(fh)
        for ligne in fh:
            c = ligne.split(",")
            if len(c) >= 2:
                d[int(c[0])] = float(c[1]) * 1e6
    return d


def retournement(run):
    """Premier passage a zero de vz_insert, en us. None s'il n'a pas eu lieu."""
    import csv as _csv
    t, v = [], []
    with open(os.path.join(run, "history.csv")) as fh:
        for r in _csv.DictReader(fh):
            t.append(float(r["t"]) * 1e6)
            v.append(float(r["vz_insert"]))
    t, v = np.asarray(t), np.asarray(v)
    i = np.where((v > 0.0) & (t > 100.0))[0]
    if not i.size:
        return None
    j = i[0]
    return float(t[j - 1] + (0.0 - v[j - 1]) * (t[j] - t[j - 1])
                 / (v[j] - v[j - 1]))


def mesure(run, k, secteurs):
    cmd = [sys.executable, os.path.join(ROOT, "tools", "crater_metrics.py"),
           run, "--frame", str(k), "--sectors", str(secteurs)]
    s = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT).stdout
    d = {}
    for cle, motif in MOTIFS.items():
        if cle == "t":
            continue
        m = re.search(motif, s)
        if m:
            d[cle] = float(m.group(1))
    return d if "rCrater" in d else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run")
    p.add_argument("--frames", default=None, help="ex. 10-18 ; defaut : toutes")
    p.add_argument("--sectors", type=int, default=16)
    p.add_argument("--cache", default=None)
    p.add_argument("--stem", default="results/fig/evolution_cratere")
    a = p.parse_args()

    temps = temps_des_trames(a.run)
    if a.frames:
        i, j = (int(x) for x in a.frames.split("-"))
        ks = [k for k in sorted(temps) if i <= k <= j]
    else:
        ks = sorted(temps)
    ks = [k for k in ks
          if os.path.exists(os.path.join(a.run, "fdem3d_joints_%04d.vtu" % k))]

    cache = {}
    chemin = a.cache or (a.stem + "_cache.json")
    if os.path.exists(chemin):
        cache = json.load(open(chemin))

    lignes = []
    for k in ks:
        cle = str(k)
        if cle not in cache:
            d = mesure(a.run, k, a.sectors)
            if d is None:
                print("  trame %2d : illisible, ignoree" % k)
                continue
            cache[cle] = d
            print("  trame %2d (%6.1f us) mesuree" % (k, temps[k]))
        d = dict(cache[cle])
        d["t"] = temps[k]
        lignes.append(d)
    json.dump(cache, open(chemin, "w"), indent=1)
    if len(lignes) < 2:
        print("pas assez de trames mesurees")
        return

    t = np.array([d["t"] for d in lignes])

    def col(c):
        return np.array([d.get(c, np.nan) for d in lignes])

    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.6), sharex=True)
    (a1, a2), (a3, a4) = axes

    a1.plot(t, col("rCrater"), "o-", color="C0",
            label=u"rayon de peau (comparable à l'essai)")
    a1.plot(t, col("rMax"), "s--", color="C3",
            label=u"extension totale (fissures comprises)")
    a1.axhline(YANG["rCrater"], color="0.35", ls=":", lw=1.2)
    a1.annotate(u"Yang ≈ 7 mm", xy=(t[-1], YANG["rCrater"]),
                xytext=(-4, -12), textcoords="offset points",
                ha="right", va="top", fontsize=8.5, color="0.35")
    a1.set_ylabel("rayon  [mm]")
    a1.legend(fontsize=8, loc="upper left")

    # le retournement, s'il a eu lieu : la ou le chargement cesse
    tr = retournement(a.run)
    if tr and t[0] < tr < t[-1] * 1.25:
        for ax in (a1, a2, a3, a4):
            ax.axvline(tr, color="0.15", lw=1.0, alpha=.65)
        import matplotlib.transforms as mtr
        for ax in (a2, a4):
            ax.text(tr, 0.04, u"retournement %.0f µs " % tr,
                    transform=mtr.blended_transform_factory(ax.transData,
                                                            ax.transAxes),
                    rotation=90, ha="right", va="bottom",
                    fontsize=8.5, color="0.15")

    a2.plot(t, col("depth"), "o-", color="C2")
    a2.set_ylabel(u"profondeur du cratère  [mm]")

    a3.plot(t, col("vDet"), "o-", color="C4", label=u"détaché")
    a3.set_ylabel(u"volume  [mm$^3$]")
    a3.set_xlabel(u"temps  [$\\mu$s]")
    a3.legend(fontsize=8)

    a4.plot(t, col("radMoy"), "o-", color="C1", label=u"radiales, moyenne")
    a4.plot(t, col("radMax"), "^--", color="C1", alpha=.6,
            label=u"radiales, maximum")
    a4.plot(t, col("brasMoy"), "v:", color="C5",
            label=u"bras endommagés, moyenne")
    a4.set_ylabel(u"portée  [mm]")
    a4.set_xlabel(u"temps  [$\\mu$s]")
    a4.legend(fontsize=8)

    for ax in (a1, a2, a3, a4):
        ax.grid(alpha=.3)
    fig.suptitle(u"St Anne rock137 — ce que le cratère fait au cours du temps",
                 fontsize=13, y=.98)
    fig.tight_layout(rect=[0, 0, 1, .96])

    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)
    for e in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, e), dpi=170, bbox_inches="tight")
    print("ecrit : %s.pdf / .png  (%d trames)" % (a.stem, len(lignes)))

    print("\n  trame   t[us]  rPeau  rTot  prof   vDet   radMoy")
    for d in lignes:
        print("  %6.1f  %5.2f %5.2f %5.2f %7.0f %6.2f"
              % (d["t"], d.get("rCrater", np.nan), d.get("rMax", np.nan),
                 d.get("depth", np.nan), d.get("vDet", np.nan),
                 d.get("radMoy", np.nan)))


main()
