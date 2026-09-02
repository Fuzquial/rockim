#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# design.py - plan d experiences (DoE) pour la calibration rockim (calib_quick).
#
#   python calib_quick/calib/design.py --test
#
# Ne lance aucun run : ne fait que produire / lire des tables de parametres.
# ---------------------------------------------------------------------------
"""Plan d experiences : espace de parametres, hypercube latin, augmentation maximin, CSV.

Algorithmes et references
-------------------------
* ParamSpace : chaque parametre est decrit par (lo, hi, echelle), echelle 'lin' ou 'log'.
  La coordonnee normalisee u dans [0, 1] vaut (x - lo) / (hi - lo) en 'lin' et
  log(x / lo) / log(hi / lo) en 'log' (utile pour ft, c, l_cz qui couvrent une decade).
  Tout l outillage aval (emulateur, optimiseur) travaille dans [0, 1]^d ; les runs
  rockim recoivent les unites physiques via from_unit.
* lhs(n, seed) : hypercube latin brouille (McKay, Beckman & Conover 1979) tire par
  scipy.stats.qmc.LatinHypercube(scramble=True, optimization='random-cd') : descente de
  coordonnees aleatoire qui minimise la discrepance L2 centree (Fang, Ma & Winker 2002),
  donc un plan stratifie par dimension ET bien reparti dans le cube.
* augment(existing, n_new, seed) : ajout sequentiel par critere maximin (Johnson, Moore &
  Ylvisaker 1990). Un grand nuage de candidats (Sobol brouille, Owen 1998) est tire dans
  [0, 1]^d, puis on retient un a un le candidat dont la distance euclidienne minimale aux
  points deja retenus (existants + nouveaux) est la plus grande (glouton). Les distances
  sont evaluees en coordonnees normalisees, donc l echelle log est respectee.
* write_csv / read_csv : colonnes = parametres (unites physiques) + colonne `id`.
  Le point decimal est toujours '.', quelle que soit la locale.

Dependances : numpy, scipy uniquement.
"""
import csv
import math
import sys

import numpy as np
from scipy.stats import qmc


class ParamSpace:
    """Espace de parametres : dict ordonne nom -> (lo, hi, 'lin' | 'log')."""

    def __init__(self, params):
        if not params:
            raise ValueError("ParamSpace vide")
        self.names = list(params.keys())
        lo, hi, scale = [], [], []
        for k in self.names:
            spec = params[k]
            if len(spec) == 2:
                a, b, s = spec[0], spec[1], "lin"
            else:
                a, b, s = spec
            a, b = float(a), float(b)
            s = str(s).lower()
            if s not in ("lin", "log"):
                raise ValueError(f"{k}: echelle '{s}' inconnue (lin | log)")
            if not a < b:
                raise ValueError(f"{k}: lo={a} doit etre < hi={b}")
            if s == "log" and a <= 0:
                raise ValueError(f"{k}: echelle log exige lo > 0 (lo={a})")
            lo.append(a); hi.append(b); scale.append(s)
        self.lo = np.array(lo); self.hi = np.array(hi); self.scale = scale
        self.is_log = np.array([s == "log" for s in scale])

    # -- proprietes ---------------------------------------------------------
    @property
    def d(self):
        return len(self.names)

    def __len__(self):
        return self.d

    def __repr__(self):
        rows = [f"  {k:12s} [{a:g}, {b:g}] {s}" for k, a, b, s in zip(self.names, self.lo, self.hi, self.scale)]
        return "ParamSpace(\n" + "\n".join(rows) + "\n)"

    def bounds_unit(self):
        """Bornes (d, 2) en coordonnees normalisees, pour l optimiseur."""
        return np.column_stack([np.zeros(self.d), np.ones(self.d)])

    # -- conversions ----------------------------------------------------------
    def _check_shape(self, A):
        A = np.asarray(A, dtype=float)
        if A.ndim == 1:
            if A.shape[0] != self.d:
                raise ValueError(f"attendu {self.d} parametres, recu {A.shape[0]}")
            return A[None, :], True
        if A.ndim != 2 or A.shape[1] != self.d:
            raise ValueError(f"attendu un tableau (n, {self.d}), recu {A.shape}")
        return A, False

    def to_unit(self, X):
        """Unites physiques -> [0, 1]^d (memes formes : (d,) -> (d,), (n, d) -> (n, d))."""
        X, flat = self._check_shape(X)
        U = np.empty_like(X)
        lin = ~self.is_log
        U[:, lin] = (X[:, lin] - self.lo[lin]) / (self.hi[lin] - self.lo[lin])
        if self.is_log.any():
            lg = self.is_log
            U[:, lg] = np.log(X[:, lg] / self.lo[lg]) / np.log(self.hi[lg] / self.lo[lg])
        return U[0] if flat else U

    def from_unit(self, U):
        """[0, 1]^d -> unites physiques."""
        U, flat = self._check_shape(U)
        X = np.empty_like(U)
        lin = ~self.is_log
        X[:, lin] = self.lo[lin] + U[:, lin] * (self.hi[lin] - self.lo[lin])
        if self.is_log.any():
            lg = self.is_log
            X[:, lg] = self.lo[lg] * np.exp(U[:, lg] * np.log(self.hi[lg] / self.lo[lg]))
        return X[0] if flat else X

    def clip(self, X):
        """Ramene dans les bornes physiques."""
        return np.clip(np.asarray(X, dtype=float), self.lo, self.hi)

    def as_dict(self, x):
        """Vecteur physique (d,) -> dict nom -> valeur (float natif)."""
        x = np.asarray(x, dtype=float).ravel()
        return {k: float(v) for k, v in zip(self.names, x)}

    # -- plans ------------------------------------------------------------------
    def lhs(self, n, seed=None, optimization="random-cd"):
        """Hypercube latin brouille, optimise (discrepance centree), en unites physiques."""
        if n < 1:
            raise ValueError("n >= 1")
        # random-cd n a pas de sens pour un seul point
        opt = optimization if n > 1 else None
        sampler = qmc.LatinHypercube(d=self.d, scramble=True, optimization=opt, seed=seed)
        return self.from_unit(sampler.random(n))

    def augment(self, existing, n_new, seed=None, n_cand=None):
        """Ajoute n_new points (unites physiques) par maximin vis-a-vis de `existing`."""
        if existing is None or len(existing) == 0:
            E = np.zeros((0, self.d))
        else:
            E, _ = self._check_shape(existing)
            E = self.to_unit(E)
        U = augment_unit(E, n_new, seed=seed, n_cand=n_cand)
        return self.from_unit(U)

    # -- CSV -----------------------------------------------------------------------
    def write_csv(self, path, X, ids=None, prefix="d", start=0):
        """Ecrit un design (n, d) physique : colonnes = parametres + id."""
        X, _ = self._check_shape(X)
        n = X.shape[0]
        if ids is None:
            ids = [f"{prefix}{start + i:03d}" for i in range(n)]
        if len(ids) != n:
            raise ValueError("ids et X de tailles differentes")
        with open(path, "w", newline="", encoding="ascii") as fh:
            w = csv.writer(fh)
            w.writerow(self.names + ["id"])
            for i in range(n):
                w.writerow([f"{v:.10g}" for v in X[i]] + [str(ids[i])])

    def read_csv(self, path):
        """Lit un design ecrit par write_csv -> (X (n, d) physique, ids list)."""
        with open(path, newline="", encoding="ascii") as fh:
            rows = list(csv.reader(fh))
        if not rows:
            return np.zeros((0, self.d)), []
        header = [h.strip() for h in rows[0]]
        missing = [k for k in self.names if k not in header]
        if missing:
            raise ValueError(f"colonnes absentes dans {path}: {missing}")
        col = {k: header.index(k) for k in self.names}
        icol_id = header.index("id") if "id" in header else None
        X, ids = [], []
        for r in rows[1:]:
            if not r or all(not c.strip() for c in r):
                continue
            X.append([float(r[col[k]]) for k in self.names])
            ids.append(r[icol_id].strip() if icol_id is not None else f"d{len(ids):03d}")
        return np.array(X, dtype=float).reshape(-1, self.d), ids


# ---------------------------------------------------------------------------
# Augmentation maximin dans le cube unite
# ---------------------------------------------------------------------------
def _min_dist_to(cand, pts):
    """Distance euclidienne minimale de chaque candidat (m, d) aux points (k, d)."""
    if len(pts) == 0:
        return np.full(len(cand), np.inf)
    # (m, k) via identite ||a-b||^2 = |a|^2 + |b|^2 - 2 a.b, borne a 0 pour l arrondi
    a2 = (cand ** 2).sum(1)[:, None]; b2 = (pts ** 2).sum(1)[None, :]
    d2 = np.maximum(a2 + b2 - 2.0 * cand @ pts.T, 0.0)
    return np.sqrt(d2.min(axis=1))


def augment_unit(existing_unit, n_new, seed=None, n_cand=None):
    """Ajout glouton maximin de n_new points dans [0,1]^d, existants en coordonnees unite.

    n_cand : taille du nuage candidat (Sobol brouille) ; defaut max(4096, 512 * n_new),
    arrondi a la puissance de 2 superieure (equilibre de Sobol).
    """
    E = np.asarray(existing_unit, dtype=float)
    if E.ndim != 2:
        raise ValueError("existing_unit doit etre (k, d)")
    d = E.shape[1]
    if n_new < 1:
        return np.zeros((0, d))
    if n_cand is None:
        n_cand = max(4096, 512 * n_new)
    m = 2 ** int(math.ceil(math.log2(n_cand)))
    cand = qmc.Sobol(d=d, scramble=True, seed=seed).random(m)
    dmin = _min_dist_to(cand, E)
    if len(E) == 0:
        # sans point existant, on part du candidat le plus proche du centre
        dmin = -np.linalg.norm(cand - 0.5, axis=1)
    chosen = []
    for _ in range(n_new):
        i = int(np.argmax(dmin))
        p = cand[i]
        chosen.append(p.copy())
        # mise a jour : distance au nouveau point retenu
        dmin = np.minimum(dmin, np.linalg.norm(cand - p, axis=1))
        dmin[i] = -np.inf
    return np.array(chosen)


def min_pairwise_distance(U):
    """Plus petite distance entre deux points distincts d un nuage (n, d)."""
    U = np.asarray(U, dtype=float)
    n = len(U)
    if n < 2:
        return np.inf
    d2 = ((U[:, None, :] - U[None, :, :]) ** 2).sum(-1)
    d2[np.arange(n), np.arange(n)] = np.inf
    return float(np.sqrt(d2.min()))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def _test():
    import os
    import tempfile

    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  [OK]  " if cond else "  [FAIL] ") + msg)
        ok = ok and bool(cond)

    ps = ParamSpace({
        "ft": (10.0, 60.0, "lin"),
        "c": (5.0, 50.0, "log"),
        "phi": (10.0, 40.0, "lin"),
        "l_cz": (0.5, 5.0, "log"),
    })
    print(ps)

    # 1. aller-retour des conversions, y compris sur les bornes
    print("1. conversions to_unit / from_unit")
    rng = np.random.default_rng(0)
    U = rng.random((50, ps.d))
    X = ps.from_unit(U)
    check(np.allclose(ps.to_unit(X), U, atol=1e-12), "aller-retour unite -> physique -> unite (50 pts)")
    check(np.allclose(ps.from_unit(np.zeros(ps.d)), ps.lo) and np.allclose(ps.from_unit(np.ones(ps.d)), ps.hi),
          "0 -> lo, 1 -> hi")
    check(abs(ps.from_unit([0.0, 0.5, 0.0, 0.5])[1] - math.sqrt(5.0 * 50.0)) < 1e-12,
          "echelle log : u = 0.5 -> moyenne geometrique (c = sqrt(5 * 50))")
    check(ps.to_unit(ps.from_unit([0.3, 0.7, 0.1, 0.9])).shape == (4,), "forme (d,) preservee")

    # 2. LHS : bornes, stratification, reproductibilite
    print("2. hypercube latin")
    n = 24
    X1 = ps.lhs(n, seed=7); X2 = ps.lhs(n, seed=7); X3 = ps.lhs(n, seed=8)
    check(X1.shape == (n, ps.d), f"forme ({n}, {ps.d})")
    check(np.all(X1 >= ps.lo) and np.all(X1 <= ps.hi), "dans les bornes physiques")
    U1 = ps.to_unit(X1)
    strat = all(sorted(np.floor(U1[:, j] * n).astype(int).tolist()) == list(range(n)) for j in range(ps.d))
    check(strat, "stratification : une cellule par strate et par dimension")
    check(np.allclose(X1, X2), "reproductible a seed egal")
    check(not np.allclose(X1, X3), "differe a seed different")
    check(ps.lhs(1, seed=0).shape == (1, ps.d), "n = 1 accepte (sans optimisation)")
    print(f"     distance min entre points (unite) : {min_pairwise_distance(U1):.3f}")

    # 3. augmentation maximin vs ajout aleatoire
    print("3. augmentation maximin")
    Xe = ps.lhs(10, seed=1)
    Xn = ps.augment(Xe, 5, seed=2)
    check(Xn.shape == (5, ps.d) and np.all(Xn >= ps.lo) and np.all(Xn <= ps.hi), "5 points dans les bornes")
    Ue = ps.to_unit(Xe); Un = ps.to_unit(Xn)
    d_maximin = min_pairwise_distance(np.vstack([Ue, Un]))
    d_rand = [min_pairwise_distance(np.vstack([Ue, np.random.default_rng(s).random((5, ps.d))])) for s in range(30)]
    print(f"     distance min apres augmentation : maximin {d_maximin:.3f} | aleatoire "
          f"moy {np.mean(d_rand):.3f} max {np.max(d_rand):.3f}")
    check(d_maximin >= np.max(d_rand), "maximin >= meilleur des 30 ajouts aleatoires")
    check(d_maximin > min_pairwise_distance(Ue) * 0.5, "n ecrase pas le nuage existant")
    X0 = ps.augment(None, 4, seed=3)
    check(X0.shape == (4, ps.d), "augmentation depuis un ensemble vide")
    check(np.allclose(ps.augment(Xe, 5, seed=2), Xn), "reproductible a seed egal")

    # 4. CSV
    print("4. ecriture / lecture CSV")
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "design.csv")
    ids = [f"d{i:03d}" for i in range(len(Xe))]
    ps.write_csv(path, Xe, ids)
    head = open(path, encoding="ascii").readline().strip()
    check(head == "ft,c,phi,l_cz,id", f"en-tete = {head}")
    Xr, idr = ps.read_csv(path)
    check(np.allclose(Xr, Xe, rtol=1e-9) and idr == ids, "aller-retour valeurs + ids")
    check("," in open(path).read() and "\t" not in open(path).read(), "separateur virgule, decimal point")
    # colonnes dans un autre ordre : lecture par nom
    with open(os.path.join(tmp, "shuffled.csv"), "w", newline="", encoding="ascii") as fh:
        w = csv.writer(fh); w.writerow(["id", "l_cz", "ft", "phi", "c"])
        for i, x in zip(ids, Xe):
            w.writerow([i, x[3], x[0], x[2], x[1]])
    Xs, ids_s = ps.read_csv(os.path.join(tmp, "shuffled.csv"))
    check(np.allclose(Xs, Xe) and ids_s == ids, "lecture par nom de colonne (ordre quelconque)")
    # design augmente : ids qui continuent
    Xall = np.vstack([Xe, Xn])
    ps.write_csv(path, Xall, ids + [f"d{10 + i:03d}" for i in range(5)])
    Xr, idr = ps.read_csv(path)
    check(len(idr) == 15 and idr[-1] == "d014", "design augmente relu (15 lignes, dernier id d014)")

    print("\nRESULTAT design.py :", "tous les tests passent" if ok else "ECHEC")
    return ok


if __name__ == "__main__":
    if "--test" in sys.argv:
        sys.exit(0 if _test() else 1)
    print(__doc__)
