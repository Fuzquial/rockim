#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# apso.py - PSO a coefficients adaptatifs (APSO) pour minimiser une fonction
# sur [0, 1]^d (parametres normalises, design.ParamSpace.to_unit).
#
#   python calib_quick/calib/apso.py --test
# ---------------------------------------------------------------------------
"""Adaptive Particle Swarm Optimization (APSO).

Reference : Z.-H. Zhan, J. Zhang, Y. Li, H. S.-H. Chung, "Adaptive Particle Swarm
Optimization", IEEE Trans. Systems, Man and Cybernetics - Part B, 39(6):1362-1381, 2009.

Algorithme (minimisation)
-------------------------
Essaim de N particules de position x_i, vitesse v_i, meilleur personnel p_i, meilleur
global g. Mise a jour PSO classique (Shi & Eberhart 1998) :
    v_i <- w v_i + c1 r1 (p_i - x_i) + c2 r2 (g - x_i),   x_i <- x_i + v_i
avec vitesse bornee a vmax = 0.2 x (hi - lo) par dimension et REFLEXION aux bornes
(position repliee, vitesse inversee). A chaque generation :

1. Estimation de l etat evolutif (ESE, section III de l article). Distance moyenne de
   chaque particule aux autres d_i = 1/(N-1) sum_j ||x_i - x_j||, puis facteur evolutif
       f = (d_g - d_min) / (d_max - d_min)  dans [0, 1]
   ou d_g est la distance moyenne de la particule globalement meilleure : f grand = le
   meilleur est isole (exploration ou saut), f petit = l essaim s est resserre autour de
   lui (convergence).

2. Classification floue en quatre etats S1 exploration, S2 exploitation, S3 convergence,
   S4 saut hors d un optimum local, avec les fonctions d appartenance trapezoidales de
   l article (eq. 9-12, fig. 5) :
       S1 : 0 (f<=0.4) | 5f-2 (0.4,0.6] | 1 (0.6,0.7] | -10f+8 (0.7,0.8] | 0 (f>0.8)
       S2 : 0 (f<=0.2) | 10f-2 (0.2,0.3] | 1 (0.3,0.4] | -5f+3 (0.4,0.6] | 0 (f>0.6)
       S3 : 1 (f<=0.1) | -5f+1.5 (0.1,0.3] | 0 (f>0.3)
       S4 : 0 (f<=0.7) | 5f-3.5 (0.7,0.9] | 1 (f>0.9)
   Defuzzification par singleton : etat d appartenance maximale ; dans une zone de
   recouvrement on suit la sequence naturelle S1 -> S2 -> S3 -> S4 -> S1 : on prend le
   successeur de l etat precedent s il est candidat, sinon on garde l etat precedent s il
   est candidat, sinon le maximum.

3. Adaptation des coefficients (section IV) :
   - inertie w(f) = 1 / (1 + 1.5 exp(-2.6 f)), bornee a [0.4, 0.9] ;
   - c1, c2 par etat, avec un increment delta tire uniformement dans [0.05, 0.1]
     (increment 'leger' = delta / 2) :
         S1 exploration  : c1 += delta,    c2 -= delta
         S2 exploitation : c1 += delta/2,  c2 -= delta/2
         S3 convergence  : c1 += delta/2,  c2 += delta/2
         S4 saut         : c1 -= delta,    c2 += delta
     chacun borne a [1.5, 2.5] ; si c1 + c2 > 4, renormalisation c_k <- 4 c_k / (c1 + c2).

4. Apprentissage elitiste (ELS, section IV-C) en etat de convergence : une dimension k
   tiree au hasard du meilleur global est perturbee
       g'_k = g_k + (hi_k - lo_k) N(0, sigma^2),  sigma = 1.0 - 0.9 t / T   (de 1.0 a 0.1)
   si f(g') < f(g), g' remplace le meilleur global (et la particule qui le porte) ; sinon
   g' remplace la pire particule courante de l essaim.

Interface
---------
apso(f_obj, bounds_unit, n_particles=30, iters=200, seed=None) -> (x_best, f_best, hist)
    f_obj : x (d,) -> float (ou (N, d) -> (N,) si vectorized=True) ; NaN/inf = +inf.
    bounds_unit : (d, 2) [lo, hi] - normalement des zeros et des uns.
    hist : dict de tableaux par generation (f_best, f_evol, state, w, c1, c2, n_eval,
           spread) + x_best final ; state code 1..4.
    Arret anticipe optionnel : ftol + patience (amelioration relative de f_best inferieure
    a ftol pendant `patience` generations).

Dependances : numpy uniquement.
"""
import sys

import numpy as np

STATE_NAMES = {1: "exploration", 2: "exploitation", 3: "convergence", 4: "jumping_out"}


# ---------------------------------------------------------------------------
# Estimation de l etat evolutif
# ---------------------------------------------------------------------------
def mean_distances(X):
    """d_i = distance euclidienne moyenne de la particule i aux N-1 autres."""
    N = len(X)
    if N < 2:
        return np.zeros(N)
    diff = X[:, None, :] - X[None, :, :]
    D = np.sqrt((diff ** 2).sum(-1))
    return D.sum(axis=1) / (N - 1)


def evolutionary_factor(X, i_best):
    """f = (d_g - d_min) / (d_max - d_min) ; 0 si l essaim est totalement resserre."""
    d = mean_distances(X)
    dmin, dmax = d.min(), d.max()
    if dmax - dmin <= 1e-300:
        return 0.0, d.mean()
    return float((d[i_best] - dmin) / (dmax - dmin)), float(d.mean())


def membership(f):
    """Appartenances (mu_S1, mu_S2, mu_S3, mu_S4) de l article, f dans [0, 1]."""
    f = min(max(float(f), 0.0), 1.0)
    # S1 exploration
    if f <= 0.4:
        s1 = 0.0
    elif f <= 0.6:
        s1 = 5.0 * f - 2.0
    elif f <= 0.7:
        s1 = 1.0
    elif f <= 0.8:
        s1 = -10.0 * f + 8.0
    else:
        s1 = 0.0
    # S2 exploitation
    if f <= 0.2:
        s2 = 0.0
    elif f <= 0.3:
        s2 = 10.0 * f - 2.0
    elif f <= 0.4:
        s2 = 1.0
    elif f <= 0.6:
        s2 = -5.0 * f + 3.0
    else:
        s2 = 0.0
    # S3 convergence
    if f <= 0.1:
        s3 = 1.0
    elif f <= 0.3:
        s3 = -5.0 * f + 1.5
    else:
        s3 = 0.0
    # S4 saut hors d un optimum local
    if f <= 0.7:
        s4 = 0.0
    elif f <= 0.9:
        s4 = 5.0 * f - 3.5
    else:
        s4 = 1.0
    return (s1, s2, s3, s4)


def classify(f, prev_state):
    """Defuzzification singleton avec suivi de la sequence S1->S2->S3->S4->S1."""
    mu = membership(f)
    cand = [k + 1 for k in range(4) if mu[k] > 0.0]
    if len(cand) == 1:
        return cand[0]
    if not cand:
        return int(np.argmax(mu)) + 1
    nxt = prev_state % 4 + 1
    if nxt in cand:
        return nxt
    if prev_state in cand:
        return prev_state
    return int(np.argmax(mu)) + 1


def inertia(f):
    """w(f) = 1 / (1 + 1.5 exp(-2.6 f)), bornee a [0.4, 0.9]."""
    return float(np.clip(1.0 / (1.0 + 1.5 * np.exp(-2.6 * f)), 0.4, 0.9))


def adapt_coefficients(c1, c2, state, rng, c_lo=1.5, c_hi=2.5, c_sum=4.0):
    """Adaptation de c1, c2 selon l etat (table de la section IV-B)."""
    delta = rng.uniform(0.05, 0.1)
    if state == 1:
        c1, c2 = c1 + delta, c2 - delta
    elif state == 2:
        c1, c2 = c1 + 0.5 * delta, c2 - 0.5 * delta
    elif state == 3:
        c1, c2 = c1 + 0.5 * delta, c2 + 0.5 * delta
    else:
        c1, c2 = c1 - delta, c2 + delta
    c1 = float(np.clip(c1, c_lo, c_hi)); c2 = float(np.clip(c2, c_lo, c_hi))
    if c1 + c2 > c_sum:
        s = c1 + c2
        c1, c2 = c_sum * c1 / s, c_sum * c2 / s
    return c1, c2


# ---------------------------------------------------------------------------
# Bornes
# ---------------------------------------------------------------------------
def reflect(X, V, lo, hi, max_pass=8):
    """Reflexion aux bornes : position repliee, vitesse inversee sur la dimension."""
    X = np.array(X, dtype=float, copy=True)
    V = None if V is None else np.array(V, dtype=float, copy=True)
    for _ in range(max_pass):
        below = X < lo
        above = X > hi
        if not (below.any() or above.any()):
            break
        X = np.where(below, 2.0 * lo - X, X)
        X = np.where(above, 2.0 * hi - X, X)
        if V is not None:
            V = np.where(below | above, -V, V)
    X = np.clip(X, lo, hi)
    return X, V


# ---------------------------------------------------------------------------
# APSO
# ---------------------------------------------------------------------------
def apso(f_obj, bounds_unit, n_particles=30, iters=200, seed=None, vmax_frac=0.2,
         x0=None, vectorized=False, ftol=None, patience=None, sigma_max=1.0, sigma_min=0.1,
         verbose=False, callback=None):
    """Minimise f_obj sur les bornes -> (x_best (d,), f_best, hist)."""
    rng = np.random.default_rng(seed)
    b = np.asarray(bounds_unit, dtype=float)
    if b.ndim != 2 or b.shape[1] != 2:
        raise ValueError("bounds_unit doit etre (d, 2)")
    lo, hi = b[:, 0].copy(), b[:, 1].copy()
    if not np.all(hi > lo):
        raise ValueError("bornes : hi > lo requis")
    d = len(lo)
    span = hi - lo
    vmax = vmax_frac * span
    N = int(n_particles)
    if N < 2:
        raise ValueError("n_particles >= 2")

    n_eval = [0]

    def evaluate(P):
        P = np.atleast_2d(P)
        n_eval[0] += len(P)
        if vectorized:
            F = np.asarray(f_obj(P), dtype=float).ravel()
        else:
            F = np.array([float(f_obj(p)) for p in P], dtype=float)
        F[~np.isfinite(F)] = np.inf
        return F

    # initialisation
    X = lo + rng.random((N, d)) * span
    if x0 is not None:
        X[0] = np.clip(np.asarray(x0, dtype=float), lo, hi)
    V = rng.uniform(-vmax, vmax, (N, d))
    F = evaluate(X)
    P = X.copy(); Fp = F.copy()
    ib = int(np.argmin(Fp)); G = P[ib].copy(); Fg = float(Fp[ib])
    w, c1, c2 = 0.9, 2.0, 2.0
    state = 1
    keys = ("f_best", "f_evol", "state", "w", "c1", "c2", "n_eval", "spread")
    hist = {k: [] for k in keys}
    stall = 0
    f_ref = Fg

    for it in range(int(iters)):
        # 1-2. etat evolutif
        f_evol, spread = evolutionary_factor(X, ib)
        state = classify(f_evol, state)
        # 3. coefficients
        w = inertia(f_evol)
        c1, c2 = adapt_coefficients(c1, c2, state, rng)
        # PSO
        r1 = rng.random((N, d)); r2 = rng.random((N, d))
        V = w * V + c1 * r1 * (P - X) + c2 * r2 * (G - X)
        V = np.clip(V, -vmax, vmax)
        X, V = reflect(X + V, V, lo, hi)
        F = evaluate(X)
        better = F < Fp
        P[better] = X[better]; Fp[better] = F[better]
        ib = int(np.argmin(Fp))
        if Fp[ib] < Fg:
            G = P[ib].copy(); Fg = float(Fp[ib])
        # 4. apprentissage elitiste en convergence
        if state == 3:
            sigma = sigma_max - (sigma_max - sigma_min) * it / max(iters - 1, 1)
            k = int(rng.integers(d))
            Gn = G.copy()
            Gn[k] += span[k] * rng.normal(0.0, sigma)
            Gn, _ = reflect(Gn, None, lo, hi)
            fn = float(evaluate(Gn)[0])
            if fn < Fg:
                G, Fg = Gn.copy(), fn
                X[ib] = Gn; F[ib] = fn; P[ib] = Gn; Fp[ib] = fn
            else:
                iw = int(np.argmax(F))
                X[iw] = Gn; F[iw] = fn
                if fn < Fp[iw]:
                    P[iw] = Gn; Fp[iw] = fn
                    ib = int(np.argmin(Fp))
        # journal
        for k, v in zip(keys, (Fg, f_evol, state, w, c1, c2, n_eval[0], spread)):
            hist[k].append(v)
        if verbose and (it % max(iters // 10, 1) == 0 or it == iters - 1):
            print(f"  it {it:5d}  f_best {Fg:12.4e}  f {f_evol:.2f} {STATE_NAMES[state]:12s} "
                  f"w {w:.2f} c1 {c1:.2f} c2 {c2:.2f}  spread {spread:.3f}")
        if callback is not None:
            callback(it, G, Fg, hist)
        # arret anticipe
        if ftol is not None and patience is not None:
            if f_ref - Fg > ftol * max(abs(f_ref), 1e-12):
                f_ref, stall = Fg, 0
            else:
                stall += 1
                if stall >= patience:
                    break

    hist = {k: np.array(v) for k, v in hist.items()}
    hist["x_best"] = G.copy()
    hist["iters_done"] = len(hist["f_best"])
    return G.copy(), Fg, hist


# ---------------------------------------------------------------------------
# Fonctions test (unites physiques via une affine [0,1] -> [lo, hi])
# ---------------------------------------------------------------------------
def _sphere(x):
    return float(np.sum(x ** 2))


def _rosenbrock(x):
    return float(np.sum(100.0 * (x[1:] - x[:-1] ** 2) ** 2 + (1.0 - x[:-1]) ** 2))


def _ackley(x):
    n = len(x)
    return float(-20.0 * np.exp(-0.2 * np.sqrt(np.sum(x ** 2) / n)) - np.exp(np.sum(np.cos(2 * np.pi * x)) / n) + 20.0 + np.e)


def _rastrigin(x):
    return float(10.0 * len(x) + np.sum(x ** 2 - 10.0 * np.cos(2 * np.pi * x)))


def _unit_wrap(f, lo, hi):
    lo = np.asarray(lo, dtype=float); hi = np.asarray(hi, dtype=float)
    return lambda u: f(lo + u * (hi - lo))


def _test():
    import time

    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  [OK]  " if cond else "  [FAIL] ") + msg)
        ok = ok and bool(cond)

    # 0. composantes
    print("0. composantes de l ESE")
    mu = np.array([membership(f) for f in np.linspace(0, 1, 101)])
    check(np.all(mu >= 0) and np.all(mu <= 1), "appartenances dans [0, 1]")
    check(membership(0.0) == (0, 0, 1, 0) and membership(0.35) == (0, 1, 0, 0)
          and membership(0.65) == (1, 0, 0, 0) and membership(0.95) == (0, 0, 0, 1),
          "S3 en f=0, S2 en 0.35, S1 en 0.65, S4 en 0.95")
    check(classify(0.5, 1) == 2 and classify(0.5, 2) == 2 and classify(0.75, 3) == 4 and classify(0.75, 4) == 1,
          "recouvrement : S1->S2 (f=0.5), S2 reste (f=0.5), S3->S4 (f=0.75), S4->S1 (f=0.75)")
    check(abs(inertia(0.0) - 0.4) < 1e-9 and abs(inertia(1.0) - 0.9) < 0.02 and inertia(0.5) > inertia(0.2),
          "inertie : 0.4 en f=0, ~0.9 en f=1, croissante")
    rng = np.random.default_rng(0)
    c1, c2 = 2.0, 2.0
    for _ in range(100):
        c1, c2 = adapt_coefficients(c1, c2, 1, rng)
    check(1.5 <= c2 <= c1 <= 2.5 and c1 + c2 <= 4.0 + 1e-12, f"100 pas en exploration : c1 {c1:.2f} >= c2 {c2:.2f}, somme <= 4")
    for _ in range(100):
        c1, c2 = adapt_coefficients(c1, c2, 3, rng)
    check(c1 + c2 <= 4.0 + 1e-12, "100 pas en convergence : somme bornee a 4")
    Xr, Vr = reflect(np.array([[-0.2, 1.3, 0.5]]), np.array([[-1.0, 1.0, 1.0]]), np.zeros(3), np.ones(3))
    check(np.allclose(Xr, [[0.2, 0.7, 0.5]]) and np.allclose(Vr, [[1.0, -1.0, 1.0]]), "reflexion : position repliee, vitesse inversee")

    # 1-4. fonctions test en d = 6 (optimum ramene en unites [0, 1] par une affine)
    # Sphere decalee : optimum en u = 0.85 (hors centre) pour que l initialisation
    # uniforme ne place pas le meilleur au centre de l essaim (sinon f ~ 0 des le debut,
    # artefact des fonctions test centrees) : la machine a etats est parcourue en entier.
    d = 6
    shift = 0.85 * (5.12 - (-5.12)) + (-5.12)          # 3.584 en unites physiques
    _sphere_shift = lambda x: _sphere(x - shift)         # noqa: E731
    # (nom, f, lo, hi, x*, generations, tolerance, verdict strict)
    cases = [
        ("Sphere", _sphere, -5.12, 5.12, 0.0, 200, 1e-3, True),
        ("Sphere decalee", _sphere_shift, -5.12, 5.12, shift, 200, 1e-3, True),
        ("Rosenbrock", _rosenbrock, -2.048, 2.048, 1.0, 5000, 1e-3, True),
        ("Ackley", _ackley, -32.768, 32.768, 0.0, 400, 1e-2, False),
        ("Rastrigin", _rastrigin, -5.12, 5.12, 0.0, 400, 1.0, False),
    ]
    for name, fn, lo, hi, xstar, iters, tol, strict in cases:
        print(f"\n{name} d={d}, bornes [{lo}, {hi}], {iters} generations x 30 particules")
        fu = _unit_wrap(fn, [lo] * d, [hi] * d)
        t0 = time.time()
        xb, fb, hist = apso(fu, np.column_stack([np.zeros(d), np.ones(d)]), n_particles=30, iters=iters, seed=1, verbose=True)
        xphys = lo + xb * (hi - lo)
        dist = float(np.max(np.abs(xphys - xstar)))
        states = np.bincount(hist["state"].astype(int), minlength=5)[1:]
        print(f"     f_best = {fb:.3e}, max|x - x*| = {dist:.3e}, evaluations {hist['n_eval'][-1]}, {time.time() - t0:.1f} s")
        print(f"     etats visites (S1..S4) : {states.tolist()}  |  f_evol moyen {hist['f_evol'].mean():.2f}")
        check(np.all(hist["w"] >= 0.4 - 1e-12) and np.all(hist["w"] <= 0.9 + 1e-12), "w dans [0.4, 0.9]")
        check(np.all(hist["c1"] + hist["c2"] <= 4.0 + 1e-9), "c1 + c2 <= 4")
        check(np.all(np.diff(hist["f_best"]) <= 0), "f_best monotone decroissante")
        check(np.all(xb >= 0) and np.all(xb <= 1), "x_best dans [0, 1]^d")
        if strict:
            check(fb < tol, f"{name} : f_best < {tol:g}")
            check(dist < 0.1, f"{name} : |x - x*| < 0.1 en unites physiques")
        else:
            print(f"     (multimodale : verdict informatif) f_best < {tol:g} : {fb < tol}")
        if name == "Sphere decalee":
            # en d = 6 l essaim se resserre en 2-3 generations : on exige le passage
            # exploration/saut -> convergence, pas l exploitation intermediaire
            check((states[0] + states[3] > 0) and states[2] > 0 and hist["state"][0] in (1, 4) and hist["state"][-1] == 3,
                  "Sphere decalee : demarre en exploration/saut, finit en convergence")
            check(hist["w"][0] > 0.6 and hist["w"][-1] < 0.45, "w part haut (> 0.6) et finit bas (< 0.45)")

    # 5. reproductibilite, x0, arret anticipe, NaN
    print("\n5. divers")
    fu = _unit_wrap(_sphere, [-5.12] * d, [5.12] * d)
    B = np.column_stack([np.zeros(d), np.ones(d)])
    a = apso(fu, B, iters=50, seed=3); b = apso(fu, B, iters=50, seed=3)
    check(np.allclose(a[0], b[0]) and a[1] == b[1], "reproductible a seed egal")
    x0 = np.full(d, 0.5)
    c = apso(fu, B, iters=5, seed=3, x0=x0)
    check(c[1] <= 1e-12, "x0 = optimum injecte -> f_best = 0 des l initialisation")
    e = apso(fu, B, iters=2000, seed=3, ftol=1e-9, patience=20)
    check(e[2]["iters_done"] < 2000, f"arret anticipe apres {e[2]['iters_done']} generations")

    def f_nan(u):
        return np.nan if u[0] < 0.3 else fu(u)

    g = apso(f_nan, B, iters=100, seed=4)
    check(np.isfinite(g[1]) and g[0][0] >= 0.3, "NaN traite comme +inf, optimum trouve dans la zone finie")
    v = apso(lambda U: np.sum((U - 0.5) ** 2, axis=1), B, iters=100, seed=5, vectorized=True)
    check(v[1] < 1e-6, "mode vectorise")

    print("\nRESULTAT apso.py :", "tous les tests passent" if ok else "ECHEC")
    return ok


if __name__ == "__main__":
    if "--test" in sys.argv:
        sys.exit(0 if _test() else 1)
    print(__doc__)
