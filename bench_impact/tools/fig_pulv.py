#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_pulv.py — l'equivalent de la fig. 18 de Yang et al. 2026 (IJRMMS 206) :
# nombre d'elements pulverises (bulkDamage : D = Dmax) au cours du temps,
# et l'energie dissipee par ce canal (colonne bdWork quand elle existe).
# Ecrit le 2026-08-30 pour la campagne de replique Kuru (comble le seul
# manque de post-traitement releve par AUDIT_actifs_replique.md §4).
#
#   python bench_impact/tools/fig_pulv.py out_kuru9 --stem fig_kuru9
# ---------------------------------------------------------------------------
import argparse
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--stem", default=None, help="prefixe des fichiers de sortie")
    ap.add_argument("--ref", type=float, default=None,
                    help="valeur cible du papier (ex. 360 a 9 m/s) tracee en pointille")
    args = ap.parse_args()
    stem = args.stem or os.path.join(args.outdir, "fig_pulv")

    with open(os.path.join(args.outdir, "history.csv")) as f:
        rows = list(csv.DictReader(f))
    t = np.array([float(r["t"]) for r in rows]) * 1e6            # us
    if "nPulv" not in rows[0]:
        raise SystemExit("history.csv sans colonne nPulv : bulkDamage n'etait pas arme")
    npulv = np.array([float(r["nPulv"]) for r in rows])
    bd = np.array([float(r["bdWork"]) for r in rows]) if "bdWork" in rows[0] else None

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(t, npulv, color="#B3543B", lw=1.4,
            label="elements a D = Dmax (rockim)")
    if args.ref is not None:
        ax.axhline(args.ref, color="0.3", lw=0.9, ls="--",
                   label=f"papier : ~{args.ref:.0f}")
    ax.set_xlabel("temps (us)")
    ax.set_ylabel("elements pulverises")
    ax.set_title("Evolution du comptage pulverise (leur fig. 18)")
    ax.legend(frameon=False, fontsize=7)
    if bd is not None and np.any(bd != 0):
        ax2 = ax.twinx()
        ax2.plot(t, np.abs(bd), color="#4C6A9C", lw=1.0, alpha=0.8)
        ax2.set_ylabel("energie de pulverisation |bdWork| (J)", color="#4C6A9C",
                       fontsize=7)
        ax2.tick_params(axis="y", colors="#4C6A9C", labelsize=6)
    fig.tight_layout()
    fig.savefig(stem + "_fig18.png", dpi=300)
    print(f"nPulv final = {npulv[-1]:.0f} ; max = {npulv.max():.0f}"
          + (f" ; bdWork final = {bd[-1]:.3f} J" if bd is not None else ""))
    print(f"ecrit : {stem}_fig18.png")


if __name__ == "__main__":
    main()
