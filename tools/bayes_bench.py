#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# bayes_bench.py — end-to-end SYNTHETIC benchmark of the emulator-based
# Bayesian calibration pipeline (the Osthus/Kennedy-O'Hagan strategy) on
# rockim's GBM tension case, with a KNOWN truth:
#
#   1. "experiment": run rockim at a hidden theta* = (gbAlphaTen, gbAlphaCoh)
#      for several mesh seeds — the seed-to-seed scatter IS the experimental
#      uncertainty;
#   2. training design: Latin hypercube over the prior box, one rockim run
#      per point (the only expensive stage);
#   3. emulator: one Gaussian process per observable (peak stress, broken
#      joints), RBF kernel, hyperparameters by marginal likelihood;
#   4. inference: Metropolis-Hastings on the emulated likelihood;
#   5. verdict: does the 95 % credible interval contain theta*? Which
#      component is identified (posterior << prior) and which is not?
#
# The physics makes the expected answer sharp: the tension peak is governed
# by the BOUNDARY TENSILE strength (alphaTen), while alphaCoh (shear
# cohesion) barely matters for a mode-I test — the posterior must be tight
# on alphaTen and close to the prior on alphaCoh. That identifiability
# contrast is the pedagogical point (and the sanity check of the pipeline
# before spending real FDEM budget on it).
#
# Self-contained: numpy + matplotlib (scipy optional, improves GP hyperopt).
# Usage:
#   python tools/bayes_bench.py --exe ./rockim.exe [--n 12] [--out bayes_out]
# ---------------------------------------------------------------------------
import argparse
import os
import re
import subprocess
import sys

import numpy as np

CFG_TEMPLATE = """
mode = fdem
scenario = tension
mesh = voronoi
T = 6e-4
frames = 2
W = 0.03
H = 0.06
grainSize = 0.006
grainSeeding = random
lloydIters = 2
refineLevels = 0
seed = {seed}
rho = 2650
E = 50e9
nu = 0.25
ft = 10e6
cohesion = 25e6
frictionDeg = 40
Gf = 70
gfShearFactor = 10
phases = alpha beta
phase.alpha.fraction = 0.6
phase.beta.fraction = 0.4
phase.beta.E = 70e9
gbAlphaTen = {aten}
gbAlphaCoh = {acoh}
gbAlphaGf = 0.5
jointPenaltyFactor = 20
pullV = 0.05
dampingLocal = 0.7
verifyFt = false
"""


def run_model(exe, theta, seed, workdir, tag):
    """One rockim run -> observables [peak stress MPa, broken joints]."""
    cfg = CFG_TEMPLATE.format(aten=theta[0], acoh=theta[1], seed=seed)
    cfgp = os.path.join(workdir, f"bb_{tag}.cfg")
    outd = os.path.join(workdir, f"bb_{tag}")
    with open(cfgp, "w") as f:
        f.write(cfg)
    r = subprocess.run([exe, cfgp, outd], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"rockim failed for {tag}:\n{r.stdout}\n{r.stderr}")
    m = re.search(r"peak macro stress = ([0-9.eE+-]+) MPa", r.stdout)
    b = re.search(r"broken joints = (\d+)", r.stdout)
    return np.array([float(m.group(1)), float(b.group(1))])


def lhs(n, d, rng):
    """Latin hypercube in [0,1]^d."""
    x = (np.arange(n)[:, None] + rng.random((n, d))) / n
    for k in range(d):
        rng.shuffle(x[:, k])
    return x


class GP:
    """Minimal RBF Gaussian process with per-dim lengthscales + nugget."""

    def __init__(self, X, y):
        self.X, self.ym, self.ys = X, y.mean(), y.std() + 1e-12
        self.y = (y - self.ym) / self.ys
        self.opt()

    def k(self, A, B, ls, amp):
        d2 = ((A[:, None, :] - B[None, :, :]) / ls) ** 2
        return amp * np.exp(-0.5 * d2.sum(-1))

    def nll(self, p):
        ls, amp, ng = np.exp(p[:-2]), np.exp(p[-2]), np.exp(p[-1])
        K = self.k(self.X, self.X, ls, amp) + ng * np.eye(len(self.X))
        try:
            L = np.linalg.cholesky(K)
        except np.linalg.LinAlgError:
            return 1e30
        a = np.linalg.solve(L.T, np.linalg.solve(L, self.y))
        return 0.5 * self.y @ a + np.log(np.diag(L)).sum()

    def opt(self):
        p0 = np.log([0.3, 0.3, 1.0, 1e-3])
        try:
            from scipy.optimize import minimize
            self.p = minimize(self.nll, p0, method="Nelder-Mead",
                              options={"maxiter": 400}).x
        except ImportError:                     # grid fallback
            best, self.p = 1e30, p0
            for l1 in (0.15, 0.3, 0.6):
                for l2 in (0.15, 0.3, 0.6):
                    for ng in (1e-4, 1e-2):
                        p = np.log([l1, l2, 1.0, ng])
                        v = self.nll(p)
                        if v < best:
                            best, self.p = v, p
        ls, amp, ng = np.exp(self.p[:-2]), np.exp(self.p[-2]), np.exp(self.p[-1])
        K = self.k(self.X, self.X, ls, amp) + ng * np.eye(len(self.X))
        self.L = np.linalg.cholesky(K)
        self.a = np.linalg.solve(self.L.T, np.linalg.solve(self.L, self.y))
        self.ls, self.amp, self.ng = ls, amp, ng

    def predict(self, Xs):
        Ks = self.k(Xs, self.X, self.ls, self.amp)
        mu = Ks @ self.a
        v = np.linalg.solve(self.L, Ks.T)
        var = np.maximum(self.amp - (v ** 2).sum(0), 1e-12)
        return mu * self.ys + self.ym, var * self.ys ** 2

    def loo(self):
        """Leave-one-out check: fraction of points inside 2 sigma."""
        Ki = np.linalg.inv(self.k(self.X, self.X, self.ls, self.amp)
                           + self.ng * np.eye(len(self.X)))
        mu = self.y - Ki @ self.y / np.diag(Ki)
        var = 1.0 / np.diag(Ki)
        z = np.abs(self.y - mu) / np.sqrt(var)
        return (z < 2.0).mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="./rockim.exe")
    ap.add_argument("--n", type=int, default=12, help="LHS training runs")
    ap.add_argument("--out", default="bayes_out")
    ap.add_argument("--truth", type=float, nargs=2, default=[0.45, 0.60])
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(0)
    lo, hi = np.array([0.2, 0.2]), np.array([1.0, 1.0])

    # -- 1. synthetic experiment: 3 mesh seeds at the hidden truth ----------
    print(f"[bench] truth theta* = {args.truth} (hidden from the inference)")
    reps = np.array([run_model(args.exe, args.truth, s, args.out, f"exp{s}")
                     for s in (101, 202, 303)])
    y_obs = reps.mean(0)
    sig_exp = reps.std(0, ddof=1) + np.array([0.05, 0.5])   # floor
    print(f"[bench] observables (peak MPa, broken): {y_obs} +- {sig_exp}")

    # -- 2. LHS design + training runs --------------------------------------
    X = lo + lhs(args.n, 2, rng) * (hi - lo)
    Y = np.array([run_model(args.exe, X[i], 12345, args.out, f"tr{i}")
                  for i in range(args.n)])
    print(f"[bench] {args.n} training runs done")

    # -- 3. emulators (input scaled to [0,1]) --------------------------------
    Xu = (X - lo) / (hi - lo)
    gps = [GP(Xu, Y[:, k]) for k in range(Y.shape[1])]
    for k, gp in enumerate(gps):
        print(f"[bench] GP{k}: lengthscales {gp.ls.round(3)}, "
              f"LOO 2-sigma coverage {gp.loo():.2f}")

    # -- 4. Metropolis-Hastings ---------------------------------------------
    def logpost(th):
        if np.any(th < lo) or np.any(th > hi):
            return -np.inf
        u = ((th - lo) / (hi - lo))[None, :]
        lp = 0.0
        for k, gp in enumerate(gps):
            mu, var = gp.predict(u)
            v = var[0] + sig_exp[k] ** 2
            lp += -0.5 * ((y_obs[k] - mu[0]) ** 2 / v + np.log(2 * np.pi * v))
        return lp

    chains = []
    for c in range(4):
        th = lo + rng.random(2) * (hi - lo)
        lp = logpost(th)
        acc, S = 0, []
        for it in range(4000):
            prop = th + rng.normal(0, 0.05, 2)
            lpp = logpost(prop)
            if np.log(rng.random()) < lpp - lp:
                th, lp = prop, lpp
                acc += 1
            S.append(th.copy())
        chains.append(np.array(S[1000:]))
        print(f"[bench] chain {c}: acceptance {acc / 4000:.2f}")
    S = np.vstack(chains)

    # Gelman-Rubin
    cm = np.array([c.mean(0) for c in chains])
    cv = np.array([c.var(0, ddof=1) for c in chains])
    W, B = cv.mean(0), len(chains[0]) * cm.var(0, ddof=1)
    rhat = np.sqrt((W * (1 - 1 / len(chains[0])) + B / len(chains[0])) / W)
    print(f"[bench] R-hat = {rhat.round(3)}")

    # -- 5. verdict ----------------------------------------------------------
    names = ["gbAlphaTen", "gbAlphaCoh"]
    prior_sd = (hi - lo) / np.sqrt(12.0)
    print("\n[bench] ---- posterior ----")
    ok = True
    for k in range(2):
        q = np.percentile(S[:, k], [2.5, 50, 97.5])
        ident = S[:, k].std() / prior_sd[k]
        inside = q[0] <= args.truth[k] <= q[2]
        ok &= inside
        print(f"[bench] {names[k]}: median {q[1]:.3f}  CI95 "
              f"[{q[0]:.3f}, {q[2]:.3f}]  truth {args.truth[k]}"
              f"  {'INSIDE' if inside else 'OUTSIDE'}"
              f"  (posterior/prior sd = {ident:.2f}"
              f" -> {'identified' if ident < 0.5 else 'NOT identified'})")
    print(f"[bench] pipeline verdict: "
          f"{'PASS — truth recovered' if ok else 'FAIL — truth outside CI95'}")

    # -- corner plot ---------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(6, 6))
    for k in range(2):
        axs[k][k].hist(S[:, k], bins=40, color="#1e4b8c", alpha=0.8)
        axs[k][k].axvline(args.truth[k], color="r")
        axs[k][k].set_xlim(lo[k], hi[k])
        axs[k][k].set_title(names[k], fontsize=9)
    axs[1][0].hist2d(S[:, 0], S[:, 1], bins=40, cmap="Blues")
    axs[1][0].plot(*args.truth, "r+", ms=12, mew=2)
    axs[1][0].set_xlabel(names[0])
    axs[1][0].set_ylabel(names[1])
    axs[0][1].axis("off")
    axs[0][1].text(0.05, 0.5, "rouge = vérité cachée\nMH 4 chaînes,\n"
                   f"{args.n} runs FDEM", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "posterior_corner.png"), dpi=150)
    np.savetxt(os.path.join(args.out, "posterior_samples.csv"), S,
               delimiter=",", header="gbAlphaTen,gbAlphaCoh")
    print(f"[bench] wrote {args.out}/posterior_corner.png + samples")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
