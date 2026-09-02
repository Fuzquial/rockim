#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_controls.py — courbes q(eps) des runs de controle, groupees par
# confinement, contre la bande experimentale Red Bohus (cibles converties en
# deformation plane : eps x (1 - nu^2) a E physique). Deviateur CORRIGE
# (q = sigma - sigma3, extract.load_run). Marqueurs : premiere insertion (CI),
# premiere rupture (CD), pic.
#
#   python calib_quick/calib/plot_controls.py out_a out_b ... [--out fig.png] [--labels "a" "b" ...]
# ---------------------------------------------------------------------------
import argparse
import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from extract import load_run, observables, TARGETS, seuils_exp  # noqa: E402

plt.rcParams.update({"font.family": "serif", "font.size": 9})


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--out", default="fig_controls.png")
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--title", default="Runs de contrôle — déviateur corrigé, cibles × (1−ν²)")
    a = ap.parse_args()
    targets = json.load(open(TARGETS))["confinements"]
    seuils = seuils_exp()
    runs = []
    for i, d in enumerate(a.dirs):
        if not os.path.exists(os.path.join(d, "history.csv")):
            print("absent :", d); continue
        r = load_run(d)
        r["name"] = a.labels[i] if a.labels and i < len(a.labels) else os.path.basename(d).replace("out_", "")
        r["obs"] = observables(r, {"confinements": targets})
        runs.append(r)
    confs = sorted({int(round(r["s3"])) for r in runs})
    n = len(confs)
    fig, axs = plt.subplots(2, n, figsize=(5.2 * n, 8.0), squeeze=False, gridspec_kw={"height_ratios": [2.2, 1]})
    cmap = plt.get_cmap("tab10")
    for j, c in enumerate(confs):
        ax, ax2 = axs[0, j], axs[1, j]
        tg = targets.get(str(c))
        sub = [r for r in runs if int(round(r["s3"])) == c]
        f = sub[0]["eps_factor"] if sub else 1.0
        if tg:
            eT = np.array(tg["eps_grid_microstrain"]) * 1e-6 * f * 100; qT = np.array(tg["q_mean_MPa"]); sT = np.array(tg["q_std_MPa"])
            ax.fill_between(eT, qT - sT, qT + sT, color="0.85", label=f"exp ±1σ (n = {int(np.max(tg['n_rep'])) if isinstance(tg.get('n_rep'), list) else tg.get('n_rep', '?')})")
            ax.plot(eT, qT, color="0.3", lw=1.2)
            se = seuils.get(str(c))
            if se:
                ax.axhline(se["q_CI"], color="0.5", ls=":", lw=0.8); ax.text(0.01, se["q_CI"] + 4, "CI exp", fontsize=7, color="0.4")
                ax.axhline(se["q_CD"], color="0.5", ls="--", lw=0.8); ax.text(0.01, se["q_CD"] + 4, "CD exp", fontsize=7, color="0.4")
        for k, r in enumerate(sub):
            col = cmap(k % 10)
            o = r["obs"]
            ax.plot(100 * r["eps"], r["q"], color=col, lw=1.3,
                    label=f"{r['name']}: pic {o['q_peak']:.0f}, CI {o['CI_frac']:.2f}, CD {o['CD_frac']:.2f}, chute {o['drop']:.2f}")
            ipk = int(np.nanargmax(r["q"]))
            ax.plot(100 * r["eps"][ipk], r["q"][ipk], "o", color=col, ms=4)
            if np.isfinite(o["q_CI"]):
                i = int(np.argmax(r["ni"] >= 1)); ax.plot(100 * r["eps"][i], r["q"][i], "^", color=col, ms=5, mfc="none")
            if np.isfinite(o["q_CD"]):
                i = int(np.argmax(r["nb"] >= 1)); ax.plot(100 * r["eps"][i], r["q"][i], "s", color=col, ms=4, mfc="none")
            ax2.plot(100 * r["eps"], r["ni"], color=col, lw=1.0)
            ax2.plot(100 * r["eps"], r["nb"], color=col, lw=1.0, ls="--")
        ax.set_title(f"σ₃ = {c} MPa"); ax.set_xlabel("ε axiale [%]"); ax.set_ylabel("q = σ₁ − σ₃ [MPa]")
        ax.set_xlim(0, 1.3); ax.set_ylim(0, None); ax.legend(fontsize=7, loc="lower right")
        ax2.set_xlabel("ε axiale [%]"); ax2.set_ylabel("joints insérés (—) / rompus (- -)"); ax2.set_xlim(0, 1.3)
    fig.suptitle(a.title + "   (▲ première insertion = CI, □ première rupture = CD, ● pic)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(a.out, dpi=140)
    print("figure :", a.out)


if __name__ == "__main__":
    main()
