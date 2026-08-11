#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# calibrate_bohus.py — calibration bayésienne du GBM rockim sur les cibles
# RED BOHUS de la base (phd/CONTINUUM.md §1-2) :
#
#   * traction directe quasi-statique : sigma_t = 18,3 MPa
#     (échelle de Weibull 18,7 MPa @ V_eff = 189 mm3, m = 23,
#      Saadati/Shariati 2022 -> moyenne 18,7*Gamma(1+1/23))
#   * UCS = 126,6 MPa (triaxial réel sigma3 = 0, Dumoulin 2024)
#
# HISTORIQUE DES PARAMÉTRISATIONS (les leçons de la passe 1 et des sondes) :
#   1. theta = (ft, cohesion), phi = 40 fixé : l'UCS plafonnait à ~180 MPa
#      quel que soit (c, phi) — le CRUSH CAP du bulk (défaut 8*cohesion)
#      passait SOUS le niveau UCS à basse cohésion et contrôlait la
#      compression (plastification en volume). Le cap est un garde-fou
#      anti-instabilité, pas de la physique -> DÉCOUPLÉ (400 MPa) ici.
#   2. Cap découplé, la vérité GBM apparaît : UCS(c=16, phi=40) = 304 MPa —
#      le réseau de joints PLATS est sur-résistant en compression
#      (UCS/sigma_t ~ 15-20 contre 6,9 pour Bohus), le fait connu de la
#      littérature GBM. Le levier dominant est l'ANGLE DE FRICTION DES
#      JOINTS (phi_joint != phi macro ; calibrations DEM/GBM typiques :
#      10-25 deg).
#   -> Paramétrisation finale : theta = (cohesion_joint, phi_joint),
#      ft_joint FIXÉ à 34 MPa (sigma_t ~ 0,55*ft, confirmé deux fois,
#      quasi orthogonal aux deux autres — le GP de la passe 1 mesurait
#      lengthscale ft = inf sur l'UCS).
#
# Machinerie : LHS -> un GP par observable -> MH 4 chaînes (R-hat), corner ;
# deux essais rockim par point (traction + UCS), POOL PARALLÈLE 1 thread/run.
#
# Usage :
#   python tools/calibrate_bohus.py --exe ./rockim.exe [--n 14] [--jobs 12]
# ---------------------------------------------------------------------------
import argparse
import os
import re
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bayes_bench import GP, lhs                                   # noqa: E402

Y_TARGET = np.array([18.3, 126.6])          # MPa : sigma_t, UCS
SIG_LAB = 0.10 * Y_TARGET                   # dispersion labo declaree 10 %
FT_FIXED = 34.0e6                           # Pa, calé par la traction
THETA_LO = np.array([4.0e6, 10.0])          # cohesion [Pa], phi_joint [deg]
THETA_HI = np.array([20.0e6, 35.0])
SCALE = np.array([1e6, 1.0])                # affichage : MPa, deg
NAMES = ["cohesion_joint [MPa]", "phi_joint [deg]"]

CFG = """
mode = fdem
scenario = tension
mesh = voronoi
T = {T}
frames = 2
W = 0.03
H = 0.06
grainSize = 0.003
grainSeeding = random
lloydIters = 2
refineLevels = 0
seed = {seed}
rho = 2620
E = 52e9
nu = 0.25
ft = {ft}
cohesion = {coh}
frictionDeg = {phi}
Gf = 70
gfShearFactor = 10
gbAlphaTen = 0.5
gbAlphaCoh = 0.5
gbAlphaGf = 0.5
crushCap = 400e6
jointPenaltyFactor = 20
pullV = {pullV}
pullRamp = {ramp}
gripLateralFree = true
dampingLocal = 0.7
verifyFt = false
"""

TESTS = {
    "t": dict(T="8e-4", pullV="0.08", ramp="2e-4"),      # traction directe
    "c": dict(T="1.5e-3", pullV="-0.2", ramp="3e-4"),    # UCS
}


def make_jobs(exe, theta, seed, workdir, tag):
    jobs = []
    for tk, tv in TESTS.items():
        cfgp = os.path.join(workdir, f"cb_{tag}_{tk}.cfg")
        outd = os.path.join(workdir, f"cb_{tag}_{tk}")
        with open(cfgp, "w") as f:
            f.write(CFG.format(ft="%.6g" % FT_FIXED,
                               coh="%.6g" % theta[0],
                               phi="%.6g" % theta[1],
                               seed=seed, **tv))
        jobs.append(([exe, cfgp, outd],
                     os.path.join(workdir, f"cb_{tag}_{tk}.log"),
                     (tag, tk)))
    return jobs


def run_pool(jobs, max_par):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    pending, running, done = list(jobs), [], {}
    t0 = time.time()
    while pending or running:
        while pending and len(running) < max_par:
            cmd, log, key = pending.pop(0)
            f = open(log, "w")
            p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT,
                                 env=env)
            running.append((p, f, log, key))
        time.sleep(3)
        for it in running[:]:
            p, f, log, key = it
            if p.poll() is None:
                continue
            f.close()
            running.remove(it)
            txt = open(log).read()
            if p.returncode != 0:
                raise RuntimeError(f"rockim a échoué ({key}):\n"
                                   + txt[-800:])
            m = re.search(r"peak macro stress = ([0-9.eE+-]+) MPa", txt)
            if not m:
                raise RuntimeError(f"pic introuvable ({key}):\n" + txt[-800:])
            done[key] = float(m.group(1))
            print(f"  [pool] {key} -> {done[key]:.2f} MPa "
                  f"({len(done)}/{len(jobs)}, {time.time() - t0:.0f} s)",
                  flush=True)
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="./rockim.exe")
    ap.add_argument("--n", type=int, default=14)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--out", default="calib_bohus_out")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(7)
    lo, hi = THETA_LO, THETA_HI

    th0 = 0.5 * (lo + hi)
    jobs = []
    for s in (101, 202, 303):
        jobs.extend(make_jobs(args.exe, th0, s, args.out, f"rep{s}"))
    X = lo + lhs(args.n, 2, rng) * (hi - lo)
    for i in range(args.n):
        jobs.extend(make_jobs(args.exe, X[i], 12345, args.out, f"tr{i}"))

    print(f"[calib] ft fixé = {FT_FIXED / 1e6:.0f} MPa ; "
          f"{len(jobs)} runs rockim (pool {args.jobs} x 1 thread)", flush=True)
    res = run_pool(jobs, args.jobs)

    reps = np.array([[res[(f"rep{s}", "t")], res[(f"rep{s}", "c")]]
                     for s in (101, 202, 303)])
    sig_seed = reps.std(0, ddof=1)
    sig = np.sqrt(SIG_LAB ** 2 + sig_seed ** 2)
    print(f"[calib] réplicats au centre : {reps.mean(0).round(2)} MPa, "
          f"écart graines {sig_seed.round(2)}, sigma total {sig.round(2)}",
          flush=True)

    Y = np.array([[res[(f"tr{i}", "t")], res[(f"tr{i}", "c")]]
                  for i in range(args.n)])
    np.savetxt(os.path.join(args.out, "design.csv"),
               np.column_stack([X, Y]), delimiter=",",
               header="cohesion,phiDeg,sigma_t,UCS")

    Xu = (X - lo) / (hi - lo)
    gps = [GP(Xu, Y[:, k]) for k in range(2)]
    for k, gp in enumerate(gps):
        print(f"[calib] GP{k}: lengthscales {gp.ls.round(3)}, "
              f"LOO 2-sigma {gp.loo():.2f}", flush=True)

    def logpost(th):
        if np.any(th < lo) or np.any(th > hi):
            return -np.inf
        u = ((th - lo) / (hi - lo))[None, :]
        lp = 0.0
        for k, gp in enumerate(gps):
            mu, var = gp.predict(u)
            v = var[0] + sig[k] ** 2
            lp += -0.5 * ((Y_TARGET[k] - mu[0]) ** 2 / v
                          + np.log(2 * np.pi * v))
        return lp

    chains = []
    for c in range(4):
        th = lo + rng.random(2) * (hi - lo)
        lp = logpost(th)
        acc, S = 0, []
        for it in range(6000):
            prop = th + rng.normal(0, 0.05, 2) * (hi - lo)
            lpp = logpost(prop)
            if np.log(rng.random()) < lpp - lp:
                th, lp = prop, lpp
                acc += 1
            S.append(th.copy())
        chains.append(np.array(S[1500:]))
        print(f"[calib] chaîne {c}: acceptation {acc / 6000:.2f}", flush=True)
    S = np.vstack(chains)

    cm = np.array([c.mean(0) for c in chains])
    cv = np.array([c.var(0, ddof=1) for c in chains])
    W, B = cv.mean(0), len(chains[0]) * cm.var(0, ddof=1)
    rhat = np.sqrt((W * (1 - 1 / len(chains[0])) + B / len(chains[0])) / W)
    print(f"[calib] R-hat = {rhat.round(3)}", flush=True)

    med = np.percentile(S, 50, axis=0)
    print("\n[calib] ---- postérieur ----", flush=True)
    for k in range(2):
        q = np.percentile(S[:, k], [2.5, 50, 97.5]) / SCALE[k]
        pr = (hi[k] - lo[k]) / np.sqrt(12.0)
        ident = S[:, k].std() / pr
        print(f"[calib] {NAMES[k]}: médiane {q[1]:.1f}  CI95 "
              f"[{q[0]:.1f}, {q[2]:.1f}]  (post/prior sd = {ident:.2f})",
              flush=True)

    jobs = []
    for s in (404, 505, 606):
        jobs.extend(make_jobs(args.exe, med, s, args.out, f"cf{s}"))
    res2 = run_pool(jobs, args.jobs)
    conf = np.array([[res2[(f"cf{s}", "t")], res2[(f"cf{s}", "c")]]
                     for s in (404, 505, 606)])
    print(f"\n[calib] confirmation à theta_median "
          f"({med[0] / 1e6:.1f} MPa, {med[1]:.1f} deg), ft = "
          f"{FT_FIXED / 1e6:.0f} MPa :", flush=True)
    print(f"[calib]   sigma_t = {conf[:, 0].mean():.2f} MPa "
          f"(cible {Y_TARGET[0]}), UCS = {conf[:, 1].mean():.2f} MPa "
          f"(cible {Y_TARGET[1]})", flush=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6.5))
    for k in range(2):
        axs[k][k].hist(S[:, k] / SCALE[k], bins=40, color="#1e4b8c",
                       alpha=0.85)
        axs[k][k].set_xlim(lo[k] / SCALE[k], hi[k] / SCALE[k])
        axs[k][k].set_title(NAMES[k], fontsize=9)
    axs[1][0].hist2d(S[:, 0] / SCALE[0], S[:, 1] / SCALE[1], bins=40,
                     cmap="Blues")
    axs[1][0].set_xlabel(NAMES[0])
    axs[1][0].set_ylabel(NAMES[1])
    axs[0][1].axis("off")
    axs[0][1].text(0.02, 0.30,
                   "cibles Red Bohus :\n$\\sigma_t$ = 18,3 MPa "
                   "(Saadati/Shariati QS)\nUCS = 126,6 MPa (Dumoulin)\n"
                   f"ft fixé = {FT_FIXED / 1e6:.0f} MPa\n"
                   f"{args.n} points LHS, GBM 2D 3 mm\n"
                   "alphas GB = 0,5 ; cap découplé", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "posterior_bohus.png"), dpi=150)
    fig.savefig(os.path.join(args.out, "posterior_bohus.pdf"))
    np.savetxt(os.path.join(args.out, "posterior_samples.csv"), S,
               delimiter=",", header="cohesion,phiDeg")
    print(f"[calib] écrit {args.out}/posterior_bohus.png|pdf + design.csv",
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
