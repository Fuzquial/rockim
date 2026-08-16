# -*- coding: utf-8 -*-
"""ENRICHISSEMENT ADAPTATIF de la base — l'etape qui manquait.

Principe (Jones, Schonlau & Welch 1998 pour l'idee ; arXiv:1809.10784 et
2411.17858 pour la variante CALIBRATION) : au lieu d'un plan statique, on
place les nouveaux calculs la ou l'emulateur est le plus incertain DANS LA
REGION PLAUSIBLE — celle ou le posterieur a du poids. Un plan uniforme
depense la moitie de ses points dans des zones sans interet ; c'est ce qui a
donne un R2 croise de 0,44 sur l'UCS avec 44 points en 6D.

Selection par lot avec DIVERSITE : on prend les points de forte variance
ponderee par la vraisemblance, en imposant une distance minimale entre eux
pour ne pas les empiler au meme endroit.

  python enrich.py [N]     -> ecrit enrich_points.json (N jeux, defaut 16)
"""
import json, os, sys
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import load, norm, PARAMS, OBS, TARGETS, SPACE

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
DMIN = 0.18                      # distance minimale entre nouveaux points
DISC = {"ucs": 0.10, "bts": 0.10, "tx20": 0.15, "tx20eps": 0.15,
        "ucseps": 0.20}


def main():
    X, Y, _ = load(os.path.join(BASE, "lhs_results.csv"))
    Xn, lo, hi = norm(X)
    print("base actuelle : %d jeux" % len(Xn))

    gps = {}
    for k, y in Y.items():
        ker = (ConstantKernel(1.0, (1e-3, 1e3))
               * Matern(length_scale=np.ones(len(PARAMS)),
                        length_scale_bounds=(1e-2, 1e2), nu=2.5)
               + WhiteKernel(1e-2, (1e-6, 1e1)))
        g = GaussianProcessRegressor(kernel=ker, normalize_y=True,
                                     n_restarts_optimizer=6, random_state=0)
        g.fit(Xn, y)
        gps[k] = g

    rng = np.random.default_rng(2026)
    S = rng.random((120000, len(PARAMS)))
    logL = np.zeros(len(S))
    sig2 = np.zeros(len(S))
    for k, (t, s) in TARGETS.items():
        m, sd = gps[k].predict(S, return_std=True)
        tot = np.sqrt(s ** 2 + sd ** 2 + (DISC[k] * t) ** 2)
        logL += -0.5 * ((m - t) / tot) ** 2
        sig2 += (sd / max(abs(t), 1e-9)) ** 2      # incertitude relative
    w = np.exp(logL - logL.max())                  # vraisemblance normalisee

    # score d'acquisition : incertitude x plausibilite
    score = np.sqrt(sig2) * w
    order = np.argsort(-score)

    chosen = []
    for i in order:
        if len(chosen) >= N:
            break
        x = S[i]
        if chosen and min(np.linalg.norm(x - c) for c in chosen) < DMIN:
            continue
        if min(np.linalg.norm(x - z) for z in Xn) < DMIN * 0.6:
            continue                               # deja couvert par la base
        chosen.append(x)
    print("retenus : %d points (distance min imposee %.2f)" % (len(chosen), DMIN))

    sets = {}
    for j, xn in enumerate(chosen):
        v = xn * (hi - lo) + lo
        p = {k: float(v[i]) for i, k in enumerate(PARAMS)}
        sets["E%03d" % j] = p
    json.dump(sets, open(os.path.join(BASE, "enrich_points.json"), "w"), indent=1)

    print("\naperçu des jeux proposes (predictions de l'emulateur) :")
    print("%-6s %8s %8s %8s %8s   %7s %7s" % ("tag", "cohesion", "phi", "ft",
                                              "Gf", "UCS~", "BTS~"))
    for tag, p in list(sets.items()):
        xn = (np.array([p[k] for k in PARAMS]) - lo) / (hi - lo)
        u = gps["ucs"].predict(np.atleast_2d(xn))[0]
        b = gps["bts"].predict(np.atleast_2d(xn))[0]
        print("%-6s %8.3g %8.1f %8.3g %8.1f   %7.0f %7.1f"
              % (tag, p["cohesion"], p["frictionDeg"], p["ft"], p["Gf"], u, b))
    print("\necrit", os.path.join(BASE, "enrich_points.json"))


if __name__ == "__main__":
    main()
