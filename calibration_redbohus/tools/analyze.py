# -*- coding: utf-8 -*-
"""Phases B/D — analyse : criblage, emulateur GP, inversion, Pareto, posterieur.

  python analyze.py screen                 # sensibilite (phase B)
  python analyze.py fit lhs_results.csv    # GP + ARD + inversion + posterieur

Cadre : Kennedy & O'Hagan (2001) — emulateur GP + terme de discrepance ;
criblage par longueurs ARD (MUCM ProcAutomaticRelevanceDetermination) ;
inversion multi-objectif (front de Pareto, Ye et al. 2025) ; posterieur par
Metropolis-Hastings sur l'emulateur.
"""
import csv, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
PARAMS = ["ft", "cohesion", "frictionDeg", "Gf", "gfShearFactor", "crushCap"]
SPACE = {"ft": (20e6, 60e6), "cohesion": (10e6, 90e6),
         "frictionDeg": (15.0, 60.0), "Gf": (30.0, 300.0),
         "gfShearFactor": (1.0, 12.0), "crushCap": (200e6, 1500e6)}

# cibles experimentales : (valeur, ecart-type EXPERIMENTAL)
# eps_pic ajoute le 2026-08-16 : a raideur identique les eprouvettes simulees
# cassaient deux fois trop tot — deux jeux peuvent donner le meme pic avec des
# deformations au pic tres differentes (argument de Ye et al. : calibrer sur la
# courbe, pas sur un scalaire). La cible triaxiale est tres serree (+-2 %) ;
# celle de l'UCS repose sur 2 essais locaux exploitables sur 4, d'ou son grand
# ecart-type (les deux autres sont tronques avant le pic par le decrochage).
TARGETS = {"ucs": (126.6, 21.4), "bts": (10.27, 0.98), "tx20": (424.8, 2.8),
           "tx20eps": (0.661, 0.014), "ucseps": (0.234, 0.075)}
OBS = {"ucs": "ucs_peak_MPa", "bts": "bts_sigma_t_MPa", "tx20": "tx20_peak_MPa",
       "tx20eps": "tx20_eps_pk", "ucseps": "ucs_eps_pk"}


def load(path):
    rows = list(csv.DictReader(open(path)))
    X, Y, ok = [], {k: [] for k in OBS}, []
    for r in rows:
        try:
            x = [float(r[p]) for p in PARAMS]
            y = {k: float(r[c]) for k, c in OBS.items()}
        except (ValueError, KeyError):
            continue
        X.append(x)
        for k in OBS:
            Y[k].append(y[k])
        ok.append(r)
    return np.array(X), {k: np.array(v) for k, v in Y.items()}, ok


def norm(X):
    lo = np.array([SPACE[p][0] for p in PARAMS])
    hi = np.array([SPACE[p][1] for p in PARAMS])
    return (X - lo) / (hi - lo), lo, hi


def screen():
    """Phase B : effet de chaque parametre porte a ses deux bornes."""
    path = os.path.join(BASE, "screen_results.csv")
    rows = {r["tag"]: r for r in csv.DictReader(open(path))}
    C = rows["C"]
    print("centre : UCS %s | BTS %s | tx20 %s\n"
          % (C["ucs_peak_MPa"], C["bts_sigma_t_MPa"], C["tx20_peak_MPa"]))
    print("%-15s %22s %22s %22s" % ("parametre", "UCS (lo -> hi)",
                                    "BTS (lo -> hi)", "tx20 (lo -> hi)"))
    sens = {}
    for p in PARAMS:
        line, s = "%-15s" % p, {}
        for k, col in OBS.items():
            try:
                lo = float(rows[p + "_lo"][col]); hi = float(rows[p + "_hi"][col])
                c = float(C[col])
                rel = abs(hi - lo) / max(abs(c), 1e-9) * 100.0
                s[k] = rel
                line += " %10.1f ->%9.1f" % (lo, hi)
            except (ValueError, KeyError):
                line += " %21s" % "-"
                s[k] = 0.0
        sens[p] = s
        print(line)
    print("\namplitude relative (%% du centre) — le criblage de la phase B :")
    print("%-15s %8s %8s %8s   %s" % ("parametre", "UCS", "BTS", "tx20", "verdict"))
    keep = []
    for p in PARAMS:
        s = sens[p]
        m = max(s.values())
        v = "RETENU" if m > 10.0 else "fige"
        if m > 10.0:
            keep.append(p)
        print("%-15s %8.1f %8.1f %8.1f   %s"
              % (p, s["ucs"], s["bts"], s["tx20"], v))
    json.dump(keep, open(os.path.join(BASE, "screen_keep.json"), "w"))
    print("\nretenus :", keep)


def fit(path):
    X, Y, rows = load(path)
    print("base : %d jeux valides" % len(X))
    Xn, lo, hi = norm(X)
    gps, ard = {}, {}
    for k, y in Y.items():
        ker = (ConstantKernel(1.0, (1e-3, 1e3))
               * Matern(length_scale=np.ones(len(PARAMS)),
                        length_scale_bounds=(1e-2, 1e2), nu=2.5)
               + WhiteKernel(1e-2, (1e-6, 1e1)))
        g = GaussianProcessRegressor(kernel=ker, normalize_y=True,
                                     n_restarts_optimizer=6, random_state=0)
        g.fit(Xn, y)
        gps[k] = g
        ls = g.kernel_.k1.k2.length_scale
        ard[k] = dict(zip(PARAMS, np.round(ls, 3)))
        print("  GP %-5s  R2(train) = %.3f   ARD = %s"
              % (k, g.score(Xn, y), ard[k]))
    # --- pertinence ARD : longueur courte = parametre influent -------------
    keys = list(OBS)
    print("\npertinence ARD (1/longueur, normalisee par sortie, %) :")
    print("%-15s" % "parametre" + "".join("%9s" % k for k in keys))
    for i, p in enumerate(PARAMS):
        line = "%-15s" % p
        for k in keys:
            inv = 1.0 / np.array([ard[k][q] for q in PARAMS])
            line += "%9.1f" % (inv[i] / inv.sum() * 100.0)
        print(line)

    # --- objectifs : erreur relative ponderee par l'ecart-type experimental
    def objs(xn):
        o = []
        for k, (t, s) in TARGETS.items():
            m = gps[k].predict(np.atleast_2d(xn))[0]
            o.append(abs(m - t) / t * 100.0)
        return np.array(o)

    # --- exploration dense sur l'emulateur (gratuit) -----------------------
    rng = np.random.default_rng(7)
    S = rng.random((200000, len(PARAMS)))
    P = {k: gps[k].predict(S) for k in OBS}
    err = np.stack([np.abs(P[k] - TARGETS[k][0]) / TARGETS[k][0] * 100.0
                    for k in TARGETS], axis=1)
    # front de Pareto (non domine)
    idx = np.argsort(err.sum(axis=1))[:4000]
    E, Sx = err[idx], S[idx]
    nd = []
    for i in range(len(E)):
        if not np.any(np.all(E <= E[i], axis=1) & np.any(E < E[i], axis=1)):
            nd.append(i)
    print("\nfront de Pareto : %d points non domines" % len(nd))
    best = idx[np.argmin(err[idx].sum(axis=1))]
    xb = S[best] * (hi - lo) + lo
    print("\ncompromis minimal (somme des erreurs relatives) :")
    for p, v in zip(PARAMS, xb):
        print("   %-15s %12.4g" % (p, v))
    for k in TARGETS:
        m = gps[k].predict(np.atleast_2d(S[best]))[0]
        print("   %-5s emulateur %8.2f   cible %8.2f   ecart %+6.1f %%"
              % (k, m, TARGETS[k][0], 100 * (m - TARGETS[k][0]) / TARGETS[k][0]))

    # --- posterieur bayesien (Metropolis sur l'emulateur) ------------------
    # vraisemblance gaussienne ponderee par les ecarts-types EXPERIMENTAUX,
    # plus un terme de discrepance de modele (Kennedy & O'Hagan 2001)
    disc = {"ucs": 0.10, "bts": 0.10, "tx20": 0.15,
            "tx20eps": 0.15, "ucseps": 0.20}            # 10-20 % de biais
    def logpost(xn):
        if np.any(xn < 0) or np.any(xn > 1):
            return -np.inf
        L = 0.0
        for k, (t, s) in TARGETS.items():
            m, sd = gps[k].predict(np.atleast_2d(xn), return_std=True)
            tot = np.hypot(np.hypot(s, sd[0]), disc[k] * t)
            L += -0.5 * ((m[0] - t) / tot) ** 2
        return L
    x = S[best].copy(); lp = logpost(x)
    chain, acc = [], 0
    step = 0.05
    for it in range(20000):
        xp = x + rng.normal(0, step, len(PARAMS))
        lpp = logpost(xp)
        if np.log(rng.random()) < lpp - lp:
            x, lp = xp, lpp
            acc += 1
        if it > 4000 and it % 5 == 0:
            chain.append(x.copy())
    C = np.array(chain) * (hi - lo) + lo
    print("\nposterieur (%d echantillons, acceptation %.0f %%) :"
          % (len(C), 100.0 * acc / 20000))
    print("%-15s %12s %12s %12s" % ("parametre", "median", "CI95 bas", "CI95 haut"))
    post = {}
    for i, p in enumerate(PARAMS):
        q = np.percentile(C[:, i], [50, 2.5, 97.5])
        post[p] = list(np.round(q, 4))
        print("%-15s %12.4g %12.4g %12.4g" % (p, q[0], q[1], q[2]))
    json.dump({"best": dict(zip(PARAMS, xb.tolist())), "posterior": post},
              open(os.path.join(BASE, "calibration_result.json"), "w"), indent=1)

    # --- figures ------------------------------------------------------------
    fig, ax = plt.subplots(1, len(TARGETS), figsize=(4.3 * len(TARGETS), 4))
    for a, k in zip(ax, TARGETS):
        yp = gps[k].predict(Xn)
        a.plot(Y[k], yp, "o", ms=4)
        lim = [min(Y[k].min(), yp.min()), max(Y[k].max(), yp.max())]
        a.plot(lim, lim, "k--", lw=0.8)
        a.axhline(TARGETS[k][0], color="C3", lw=1, ls=":")
        u = "%" if k.endswith("eps") else "MPa"
        a.set_xlabel("rockim (%s)" % u); a.set_ylabel("émulateur GP (%s)" % u)
        a.set_title("%s — R² = %.3f" % (k, gps[k].score(Xn, Y[k])))
    fig.suptitle("Émulateur GP (Matérn 5/2 ARD) — qualité d'ajustement")
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, "figures", "gp_fit.png"), dpi=140)

    fig, ax = plt.subplots(2, 3, figsize=(13, 7))
    for i, p in enumerate(PARAMS):
        a = ax[i // 3][i % 3]
        a.hist(C[:, i], bins=40, color="C0", alpha=0.8)
        a.axvline(xb[i], color="C3", lw=1.5)
        a.set_title(p, fontsize=10)
    fig.suptitle("Postérieur bayésien des paramètres de joints (rouge = compromis retenu)")
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, "figures", "posterior.png"), dpi=140)
    print("\nfigures ecrites dans figures/")


if __name__ == "__main__":
    if sys.argv[1] == "screen":
        screen()
    else:
        fit(sys.argv[2] if len(sys.argv) > 2
            else os.path.join(BASE, "lhs_results.csv"))
