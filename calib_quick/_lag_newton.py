#!/usr/bin/env python3
# Replica : poids de Laguerre a AIRES PRESCRITES par Newton amorti sur le dual
# semi-discret (Aurenhammer-Hoffmann-Aronov 1998 : existence et unicite a une
# constante pres ; Kitagawa-Merigot-Thibert 2019 : Newton amorti converge si
# aucune cellule ne se vide ; Bourne-Kok-Roper-Spanjer 2020 : grains Laguerre a
# volumes donnes). Depart w = 0 (Voronoi, toutes cellules non vides).
#   dA_i/dw_i = sum_j l_ij / (2 d_ij) ; dA_i/dw_j = -l_ij / (2 d_ij)   (Laplacien L)
#   Newton : L delta = aT - A ; amortissement : garder min A >= eps0.
import sys, math, numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
W, H, d = 0.02, 0.04, 0.003
s = d * math.sqrt(math.pi / (2 * math.sqrt(3)))


def clip(poly, n, b):
    out = []
    m = len(poly)
    for k in range(m):
        P, Q = poly[k], poly[(k + 1) % m]
        fp, fq = n @ P - b, n @ Q - b
        if fp <= 0: out.append(P)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq); out.append(P + t * (Q - P))
    return out


def area_centroid(poly):
    if len(poly) < 3: return 0.0, None
    a = cx = cy = 0.0
    for k in range(len(poly)):
        P, Q = poly[k], poly[(k + 1) % len(poly)]
        cr = P[0] * Q[1] - Q[0] * P[1]
        a += cr; cx += (P[0] + Q[0]) * cr; cy += (P[1] + Q[1]) * cr
    a *= 0.5
    return a, np.array([cx, cy]) / (6 * a)


def laguerre(seeds, w):
    """cellules + Laplacien (aretes i-j : l_ij / (2 d_ij), voisin par test de puissance au milieu)"""
    N = len(seeds)
    rect = [np.array(v) for v in [(0, 0), (W, 0), (W, H), (0, H)]]
    cells, L = [], np.zeros((N, N))
    for i, p in enumerate(seeds):
        poly = rect
        order = np.argsort(((seeds - p) ** 2).sum(1))
        for j in order[1:]:
            n = seeds[j] - p; d2 = n @ n
            poly = clip(poly, n, n @ p + 0.5 * (d2 + w[i] - w[j]))
            if len(poly) < 3: break
        cells.append(poly)
        if len(poly) < 3: continue
        for k in range(len(poly)):
            P, Q = poly[k], poly[(k + 1) % len(poly)]
            m = 0.5 * (P + Q); ell = np.linalg.norm(Q - P)
            pw = ((seeds - m) ** 2).sum(1) - w
            pw[i] = np.inf
            j = int(np.argmin(pw))
            if abs(pw[j] - (((p - m) ** 2).sum() - w[i])) < 1e-9 * s * s:   # voisin (pas un mur)
                dij = np.linalg.norm(seeds[j] - p)
                L[i, i] += ell / (2 * dij); L[i, j] -= ell / (2 * dij)
    return cells, L


def newton(seeds, aT, iters=30, verbose=False):
    N = len(seeds); w = np.zeros(N)
    cells, L = laguerre(seeds, w)
    A = np.array([area_centroid(c)[0] for c in cells])
    eps0 = 0.5 * min(A.min(), aT.min())
    hist = []
    for k in range(iters):
        g = aT - A
        err = (np.abs(g) / aT).max()
        hist.append(err)
        if err < 1e-3: break
        # L singulier (vecteur constant) : on fixe la moyenne de delta a 0
        Lr = L + np.ones((N, N)) / N
        delta = np.linalg.solve(Lr, g)
        tau = 1.0
        while True:
            w2 = w + tau * delta
            cells2, L2 = laguerre(seeds, w2)
            A2 = np.array([area_centroid(c)[0] for c in cells2])
            if A2.min() >= eps0 or tau < 1e-4: break
            tau *= 0.5
        w, cells, L, A = w2, cells2, L2, A2
        if verbose: print(f"      it {k+1:2d} tau {tau:.3f} err {100*err:.2f}%")
    return w, cells, hist


def seed_sorted(sigma, rng, a_excl=0.35):
    target = round(W * H / (math.sqrt(3) / 2 * s * s * math.exp(sigma ** 2)))
    want = np.sort(s * np.exp(sigma * rng.standard_normal(target) - 0.5 * sigma ** 2))[::-1]
    seeds, sp = [], []
    for spi in want:
        for _ in range(4000):
            p = np.array([rng.uniform(0, W), rng.uniform(0, H)])
            if all(np.linalg.norm(q - p) >= a_excl * (spi + sq) for q, sq in zip(seeds, sp)):
                seeds.append(p); sp.append(spi); break
    return np.array(seeds), np.array(sp), target


for sigma in (0.5, 0.8, 1.0):
    for lloyd in (0, 2):
        rng = np.random.default_rng(1)
        seeds, sp, target = seed_sorted(sigma, rng)
        aT = W * H * sp ** 2 / (sp ** 2).sum()
        nit = []
        for it in range(lloyd + 1):
            w, cells, hist = newton(seeds, aT)
            nit.append(len(hist))
            if it == lloyd: break
            for i, c in enumerate(cells):
                a, cen = area_centroid(c)
                if cen is not None: seeds[i] = np.clip(cen, [1e-8, 1e-8], [W - 1e-8, H - 1e-8])
        A = np.array([area_centroid(c)[0] for c in cells])
        deq = 2 * np.sqrt(A / math.pi); lnd = np.log(deq)
        err = (np.abs(A - aT) / aT).max()
        print(f"sigma {sigma:.1f} lloyd {lloyd} : {len(seeds)}/{target} graines, Newton {nit} it, "
              f"err max {100*err:.2f} %, sd ln d_eq = {lnd.std():.3f}, corr(ln s, ln d) = {np.corrcoef(np.log(sp), lnd)[0,1]:.3f}, "
              f"d_eq = {deq.mean()*1e3:.2f} mm, vides {int((A <= 0).sum())}")
