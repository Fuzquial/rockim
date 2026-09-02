#!/usr/bin/env python3
# Replica numerique de la graine polydisperse de Tessellation::build :
# addition sequentielle triee (exclusion a (s_i + s_j)) + diagramme de Laguerre
# de poids (kappa s_i)^2, avec ou sans Lloyd. Mesure : nombre place / vise,
# ecart-type de ln(d_eq) realise (demande sigma), d_eq moyen.
import sys, math, numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

W, H, d, sigma = 0.02, 0.04, 0.003, float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
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


def laguerre(seeds, w):
    rect = [np.array(v) for v in [(0, 0), (W, 0), (W, H), (0, H)]]
    cells = []
    for i, p in enumerate(seeds):
        poly = rect
        order = np.argsort(((seeds - p) ** 2).sum(1))
        for j in order[1:]:
            n = seeds[j] - p; d2 = n @ n
            b = n @ p + 0.5 * (d2 + w[i] - w[j])
            poly = clip(poly, n, b)
            if len(poly) < 3: break
        cells.append(poly)
    return cells


def area_centroid(poly):
    if len(poly) < 3: return 0.0, None
    a = cx = cy = 0.0
    for k in range(len(poly)):
        P, Q = poly[k], poly[(k + 1) % len(poly)]
        cr = P[0] * Q[1] - Q[0] * P[1]
        a += cr; cx += (P[0] + Q[0]) * cr; cy += (P[1] + Q[1]) * cr
    a *= 0.5
    return a, np.array([cx, cy]) / (6 * a)


def seed_sorted(a_excl, rng):
    target = round(W * H / (math.sqrt(3) / 2 * s * s * math.exp(sigma ** 2)))
    want = np.sort(s * np.exp(sigma * rng.standard_normal(target) - 0.5 * sigma ** 2))[::-1]
    seeds, sp = [], []
    skipped = 0
    for spi in want:
        placed = False
        for _ in range(4000):
            p = np.array([rng.uniform(0, W), rng.uniform(0, H)])
            if all(np.linalg.norm(q - p) >= a_excl * (spi + sq) for q, sq in zip(seeds, sp)):
                seeds.append(p); sp.append(spi); placed = True; break
        if not placed: skipped += 1
    return np.array(seeds), np.array(sp), target, skipped


def run(a_excl, kappa, lloyd, seedno=1):
    rng = np.random.default_rng(seedno)
    seeds, sp, target, skipped = seed_sorted(a_excl, rng)
    w = (kappa * sp) ** 2
    for it in range(lloyd + 1):
        cells = laguerre(seeds, w)
        keep = [i for i, c in enumerate(cells) if len(c) >= 3]
        if len(keep) < len(cells):
            seeds, w, sp = seeds[keep], w[keep], sp[keep]
            cells = laguerre(seeds, w)
        if it == lloyd: break
        for i, c in enumerate(cells):
            a, cen = area_centroid(c)
            seeds[i] = np.clip(cen, [1e-8, 1e-8], [W - 1e-8, H - 1e-8])
    areas = np.array([area_centroid(c)[0] for c in cells])
    deq = 2 * np.sqrt(areas / math.pi)
    lnd = np.log(deq)
    corr = np.corrcoef(np.log(sp), lnd)[0, 1]
    return dict(placed=len(seeds), target=target, skipped=skipped, sd=lnd.std(), dmean=deq.mean() * 1e3,
                corr=corr, cover=areas.sum() / (W * H))


print(f"sigma demande = {sigma} ; s = {s*1e3:.3f} mm ; cible d_eq moyen ~ {d*1e3:.1f} mm x exp(sigma^2/2)")
print(f"{'a_excl':>7s} {'kappa':>6s} {'lloyd':>5s} {'place/vise':>11s} {'sd ln d':>8s} {'corr':>6s} {'d_eq mm':>8s}")
for a_excl, kappa in [(0.35, 0.315), (0.35, 0.35), (0.42, 0.42), (0.5, 0.5), (0.5, 0.45)]:
    for lloyd in (0, 2):
        r = run(a_excl, kappa, lloyd)
        print(f"{a_excl:7.2f} {kappa:6.3f} {lloyd:5d} {r['placed']:5d}/{r['target']:<5d} {r['sd']:8.3f} {r['corr']:6.2f} {r['dmean']:8.2f}")
