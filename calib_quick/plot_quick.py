#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_quick.py — les trois cas rapides (homogene / Weibull / GBM) contre la
# bande experimentale Red Bohus a sigma3 = 50 MPa.
#
#   python calib_quick/plot_quick.py [out_q1_homog_P050 out_q2_weibull_P050 out_q3_gbm_P050]
#
# q(eps) : sigma = colonne `sigma` de history.csv (gripFy / (W thk), Pa),
# offset de consolidation = moyenne de sigma juste avant pullDelay ;
# eps = deplacement impose des mors (rampe cosinus analytique, comme
# simcurve.py de la campagne d aout) / H. Cibles : targets_triax_bohus.json
# (moyenne des 3 replicats sur grille commune, bande = ecart-type).
# ---------------------------------------------------------------------------
import json
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS = os.path.join(HERE, "..", "..", "rockim_f1", "calib_triax3d", "targets_triax_bohus.json")
plt.rcParams.update({"font.family": "serif", "font.size": 9})


def parse_cfg(path):
    cfg = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            cfg[k.strip()] = v.strip()
    return cfg


def grip_displacement(t, pullV, ramp, delay):
    tau = np.maximum(t - delay, 0.0)
    V = abs(pullV)
    if ramp <= 0:
        return V * tau
    u = V * (0.5 * np.minimum(tau, ramp) - (ramp / (2 * math.pi))
             * np.sin(math.pi * np.minimum(tau, ramp) / ramp))
    return u + V * np.maximum(tau - ramp, 0.0)


def load_curve(cfg_path, out_dir):
    cfg = parse_cfg(cfg_path)
    H = float(cfg.get("H", 0.04))
    hist = np.genfromtxt(os.path.join(out_dir, "history.csv"), delimiter=",",
                         names=True, invalid_raise=False)
    t = hist["t"]; sig = np.abs(hist["sigma"]); nb = hist["nBroken"]
    ok = np.isfinite(t) & np.isfinite(sig)
    t, sig, nb = t[ok], sig[ok], nb[ok]
    pullV = float(cfg.get("pullV", -0.25)); ramp = float(cfg.get("pullRamp", 0.0))
    delay = float(cfg.get("pullDelay", 0.0))
    pre = (t > 0.7 * delay) & (t <= delay)
    # CORRECTION C0 (2026-09-02, critique A5) : les mors bloques pendant pullDelay
    # laissent sigma_yy = nu/(1-nu) sigma3 (16,7 MPa a 50) et non sigma3 ; le
    # deviateur comparable a l essai est q = sigma - sigma3, avec eps = 0 quand
    # sigma atteint sigma3. L ancien offset (moyenne avant pullDelay) surestimait
    # q de 33 MPa a 50 MPa et de 13 MPa a 20 MPa.
    s3 = float(cfg.get("confiningPressure", 0))
    eps_raw = grip_displacement(t, pullV, ramp, delay) / H
    q = (sig - s3) / 1e6
    load = t > delay
    reach = np.where(load & (sig >= s3))[0]
    i0 = int(reach[0]) if len(reach) else int(np.argmax(load))
    eps = eps_raw - eps_raw[i0]
    load = load & (np.arange(len(t)) >= i0)
    return dict(eps=eps[load], q=q[load], nb=nb[load], s3=s3 / 1e6,
                sig0=s3 / 1e6, nb0=int(nb[pre][-1]) if pre.sum() else 0)


def fit_E(eps, q, frac=(0.2, 0.5)):
    ipk = int(np.nanargmax(q))
    m = (q[:ipk + 1] >= frac[0] * q[ipk]) & (q[:ipk + 1] <= frac[1] * q[ipk])
    return np.polyfit(eps[:ipk + 1][m], q[:ipk + 1][m], 1)[0] / 1e3 if m.sum() >= 5 else float("nan")  # GPa


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    runs = sys.argv[1:] or ["out_q1_homog_P050", "out_q2_weibull_P050", "out_q3_gbm_P050"]
    names = {"out_q1_homog_P050": "1 — homogène (Gmsh, sans grain)",
             "out_q2_weibull_P050": "2 — Weibull m = 8 sur les joints",
             "out_q3_gbm_P050": "3 — GBM équivalent (3 phases, α = 1)"}
    tg = json.load(open(TARGETS))["confinements"]["50"]
    eT = np.array(tg["eps_grid_microstrain"]) * 1e-6; qT = np.array(tg["q_mean_MPa"]); sT = np.array(tg["q_std_MPa"])
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2))
    ax = axs[0]
    ax.fill_between(100 * eT, qT - sT, qT + sT, color="0.8", label="exp ±1σ (Bohus, σ₃ = 50)")
    ax.plot(100 * eT, qT, color="0.35", lw=1.2)
    cols = ["#1f77b4", "#c0392b", "#2ca02c"]
    print(f"{'cas':38s} {'pic MPa':>8s} {'ε_pic %':>8s} {'E GPa':>6s} {'chute':>6s} {'rompus':>7s}")
    print(f"{'cible exp (σ₃ = 50)':38s} {tg['q_peak_mean_MPa']:8.1f} {tg['eps_peak_microstrain']/1e4:8.2f} {'~77':>6s} {100*tg['chute_fraction_moyenne']:5.0f}% {'':>7s}")
    for run, col in zip(runs, cols):
        cfg = os.path.join(HERE, run.replace("out_", "") + ".cfg")
        if not os.path.exists(os.path.join(run, "history.csv")):
            print(f"{names.get(run, run):38s}  (pas de history.csv)"); continue
        c = load_curve(cfg, run)
        ipk = int(np.nanargmax(c["q"])); qpk = c["q"][ipk]; epk = c["eps"][ipk]
        drop = (qpk - c["q"][-1]) / qpk
        E = fit_E(c["eps"], c["q"])
        ax.plot(100 * c["eps"], c["q"], color=col, lw=1.4, label=names.get(run, run))
        axs[1].plot(100 * c["eps"], c["nb"] - c["nb0"], color=col, lw=1.4)
        print(f"{names.get(run, run):38s} {qpk:8.1f} {100*epk:8.2f} {E:6.1f} {100*drop:5.0f}% {int(c['nb'][-1]-c['nb0']):7d}   ({100*(qpk/tg['q_peak_mean_MPa']-1):+.0f} % sur le pic)")
    ax.set_xlabel("déformation axiale ε [%]"); ax.set_ylabel("déviateur q [MPa]"); ax.set_xlim(0, 1.4); ax.set_ylim(0, 750)
    ax.legend(fontsize=8, loc="lower right"); ax.set_title("(a) q(ε) — trois représentations, mêmes moyennes")
    axs[1].set_xlabel("déformation axiale ε [%]"); axs[1].set_ylabel("joints rompus"); axs[1].set_title("(b) rupture des joints"); axs[1].set_xlim(0, 1.4)
    fig.suptitle("Calibration rapide Red Bohus — triaxial σ₃ = 50 MPa, éprouvette 20 × 40 mm, bulk élastique, joints de la sonde 4 (non calibrés)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(HERE, "fig_quick_P050.png"), dpi=140)
    print("figure :", os.path.join(HERE, "fig_quick_P050.png"))


if __name__ == "__main__":
    main()
