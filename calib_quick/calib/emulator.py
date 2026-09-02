#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# emulator.py - substitut (surrogate) des observables rockim en fonction des
# parametres NORMALISES dans [0, 1]^d (voir design.ParamSpace.to_unit).
#
#   python calib_quick/calib/emulator.py --test
# ---------------------------------------------------------------------------
"""Emulateur d observables : processus gaussien (GP) ou foret aleatoire (RF).

Pourquoi
--------
Un run rockim coute des heures ; l optimiseur (apso.py) a besoin de milliers d evaluations.
On apprend donc un substitut y(u) sur les runs deja faits, u = parametres normalises dans
[0, 1]^d, et on optimise le substitut. L incertitude predite sert a choisir ou relancer.

Modeles
-------
* kind='gp' : sklearn GaussianProcessRegressor, noyau
      ConstantKernel * Matern(nu=2.5, ARD : une longueur par parametre) + WhiteKernel
  (Rasmussen & Williams 2006, ch. 2, 4, 5). Le WhiteKernel est un nugget APPRIS : le bruit
  de realisation de maillage (tirage des grains) est reel et non reproductible, il ne faut
  pas interpoler exactement les points. normalize_y=True, n_restarts_optimizer=8 (maximum
  de la log-vraisemblance marginale, L-BFGS-B avec redemarrages aleatoires). Si un bruit de
  mesure connu `noise` (ecart-type par point) est fourni, il s ajoute sur la diagonale via
  alpha (converti en unites de y normalisees), en plus du nugget appris.
  predict renvoie par defaut l ecart-type de la fonction LATENTE (nugget retire) ;
  include_noise=True redonne l ecart-type d une nouvelle realisation (utile pour comparer
  a un run futur).
* kind='rf' : sklearn RandomForestRegressor(n_estimators=400, min_samples_leaf=2)
  (Breiman 2001) ; incertitude = ecart-type des predictions entre arbres (proxy grossier,
  pas un intervalle calibre). Si `noise` est fourni, poids d echantillon 1/noise^2.

Validation
----------
loo_rmse() : validation croisee leave-one-out (KFold a 10 plis si n > 60). Pour le GP,
chaque pli est re-ajuste en partant des hyperparametres du modele complet (une descente
L-BFGS-B, sans redemarrage) ; refit_hyper=False fige les hyperparametres (LOO classique
a hyperparametres fixes, Rasmussen & Williams 5.4.2, plus optimiste).
importance() : RF = permutation importance (Breiman 2001, sklearn.inspection) ;
GP = 1 / longueur ARD, normalise a somme 1 (une longueur courte = parametre influent).

Valeurs manquantes : les lignes ou y (ou une coordonnee de X) n est pas finie sont
ECARTEES a l ajustement (n_dropped les compte) - un run non rompu donne des d_CD = NaN,
il ne faut pas les imputer a zero. C est a l objectif (objective.py) de penaliser.

MultiEmulator : une instance par observable, fit(X, {obs: y}) / predict -> {obs: (mu, sd)}.

Dependances : numpy, scipy, scikit-learn.
"""
import sys
import warnings

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold, LeaveOneOut

KINDS = ("gp", "rf")


def _clean(X, y, noise):
    """Ecarte les lignes non finies ; renvoie X, y, noise, nombre de lignes ecartees."""
    X = np.atleast_2d(np.asarray(X, dtype=float))
    y = np.asarray(y, dtype=float).ravel()
    if X.shape[0] != y.shape[0]:
        raise ValueError(f"X ({X.shape[0]} lignes) et y ({y.shape[0]}) incompatibles")
    mask = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if noise is not None:
        noise = np.asarray(noise, dtype=float).ravel()
        if noise.shape[0] != y.shape[0]:
            raise ValueError("noise doit avoir une valeur par point")
        mask &= np.isfinite(noise) & (noise >= 0)
        noise = noise[mask]
    return X[mask], y[mask], noise, int((~mask).sum())


class Emulator:
    """Substitut scalaire y(u), u dans [0, 1]^d."""

    def __init__(self, kind="gp", seed=0, n_restarts=8, n_trees=400, min_samples_leaf=2,
                 names=None, length_scale_bounds=(1e-2, 1e2), noise_bounds=(1e-8, 1e1)):
        if kind not in KINDS:
            raise ValueError(f"kind {kind!r} inconnu : {KINDS}")
        self.kind = kind
        self.seed = seed
        self.n_restarts = n_restarts
        self.n_trees = n_trees
        self.min_samples_leaf = min_samples_leaf
        self.names = list(names) if names is not None else None
        self.length_scale_bounds = length_scale_bounds
        self.noise_bounds = noise_bounds
        self.model = None
        self.X_ = None; self.y_ = None; self.noise_ = None
        self.n_dropped = 0

    # -- ajustement -----------------------------------------------------------
    def _check_unit(self, X):
        X = np.atleast_2d(np.asarray(X, dtype=float))
        if self.X_ is not None and X.shape[1] != self.X_.shape[1]:
            raise ValueError(f"attendu {self.X_.shape[1]} parametres, recu {X.shape[1]}")
        fin = X[np.isfinite(X).all(axis=1)]
        if fin.size and (fin.min() < -0.05 or fin.max() > 1.05):
            raise ValueError("X doit etre en coordonnees normalisees [0, 1] (ParamSpace.to_unit)")
        return X

    def _make_gp(self, d, alpha, kernel=None, n_restarts=None):
        if kernel is None:
            kernel = (ConstantKernel(1.0, (1e-3, 1e3))
                      * Matern(length_scale=np.full(d, 0.5), length_scale_bounds=self.length_scale_bounds, nu=2.5)
                      + WhiteKernel(noise_level=1e-2, noise_level_bounds=self.noise_bounds))
        return GaussianProcessRegressor(kernel=kernel, alpha=alpha, normalize_y=True,
                                        n_restarts_optimizer=self.n_restarts if n_restarts is None else n_restarts,
                                        random_state=self.seed)

    def _make_rf(self, n_trees=None):
        return RandomForestRegressor(n_estimators=self.n_trees if n_trees is None else n_trees,
                                     min_samples_leaf=self.min_samples_leaf, random_state=self.seed, n_jobs=1)

    @staticmethod
    def _alpha(y, noise):
        """Bruit connu (ecart-type physique) -> alpha en unites de y normalisees."""
        if noise is None:
            return 1e-10
        sd = float(np.std(y))
        sd = sd if sd > 0 else 1.0
        return (noise / sd) ** 2 + 1e-10

    def fit(self, X, y, noise=None):
        """Ajuste sur X (n, d) dans [0, 1]^d et y (n,) ; noise = ecart-type connu par point."""
        X = self._check_unit(X)
        X, y, noise, self.n_dropped = _clean(X, y, noise)
        if len(y) < 3:
            raise ValueError(f"trop peu de points finis pour ajuster ({len(y)})")
        self.X_, self.y_, self.noise_ = X, y, noise
        if self.kind == "gp":
            self.model = self._make_gp(X.shape[1], self._alpha(y, noise))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ConvergenceWarning)
                self.model.fit(X, y)
        else:
            self.model = self._make_rf()
            sw = None if noise is None else 1.0 / np.maximum(noise, 1e-12) ** 2
            if sw is not None:
                sw = sw / sw.mean()
            self.model.fit(X, y, sample_weight=sw)
        return self

    # -- prediction --------------------------------------------------------------
    @property
    def y_std_(self):
        """Ecart-type de y d entrainement utilise par normalize_y (sklearn)."""
        sd = getattr(self.model, "_y_train_std", None)
        if sd is None:
            sd = np.std(self.y_)
        return float(sd) if float(sd) > 0 else 1.0

    @property
    def nugget_std(self):
        """Ecart-type du bruit de realisation appris (GP, unites physiques) ; NaN pour RF."""
        if self.kind != "gp" or self.model is None:
            return float("nan")
        return float(np.sqrt(self.model.kernel_.k2.noise_level) * self.y_std_)

    def predict(self, X, include_noise=False):
        """-> (mean (m,), std (m,)). std = fonction latente ; include_noise ajoute le nugget."""
        if self.model is None:
            raise RuntimeError("appeler fit d abord")
        X = self._check_unit(X)
        if self.kind == "gp":
            mean, std = self.model.predict(X, return_std=True)
            # sklearn inclut le WhiteKernel dans la variance predite (diag du noyau)
            nug = self.model.kernel_.k2.noise_level * self.y_std_ ** 2
            var = std ** 2 - nug
            if include_noise:
                var = var + nug
            return mean, np.sqrt(np.maximum(var, 0.0))
        per_tree = np.stack([t.predict(X) for t in self.model.estimators_], axis=0)
        return per_tree.mean(axis=0), per_tree.std(axis=0)

    def __call__(self, X):
        return self.predict(X)[0]

    # -- validation croisee --------------------------------------------------------
    def loo_rmse(self, max_loo=60, n_splits=10, refit_hyper=True, return_pred=False):
        """RMSE leave-one-out (KFold(n_splits) si n > max_loo). return_pred -> (rmse, yhat)."""
        if self.model is None:
            raise RuntimeError("appeler fit d abord")
        X, y, noise = self.X_, self.y_, self.noise_
        n = len(y)
        cv = LeaveOneOut() if n <= max_loo else KFold(n_splits=min(n_splits, n), shuffle=True, random_state=self.seed)
        yhat = np.full(n, np.nan)
        for tr, te in cv.split(X):
            if self.kind == "gp":
                kern = self.model.kernel_
                m = self._make_gp(X.shape[1], self._alpha(y[tr], None if noise is None else noise[tr]),
                                  kernel=kern, n_restarts=0)
                if not refit_hyper:
                    m.set_params(optimizer=None)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", ConvergenceWarning)
                    m.fit(X[tr], y[tr])
                yhat[te] = m.predict(X[te])
            else:
                m = self._make_rf()
                sw = None if noise is None else 1.0 / np.maximum(noise[tr], 1e-12) ** 2
                m.fit(X[tr], y[tr], sample_weight=None if sw is None else sw / sw.mean())
                yhat[te] = m.predict(X[te])
        rmse = float(np.sqrt(np.mean((yhat - y) ** 2)))
        return (rmse, yhat) if return_pred else rmse

    def r2_loo(self, **kw):
        """Coefficient de determination en validation croisee (Q2)."""
        rmse, yhat = self.loo_rmse(return_pred=True, **kw)
        ss = np.sum((self.y_ - self.y_.mean()) ** 2)
        return float(1.0 - np.sum((yhat - self.y_) ** 2) / ss) if ss > 0 else float("nan")

    # -- importance ------------------------------------------------------------------
    def importance(self, n_repeats=10):
        """Importance normalisee (somme 1) par parametre ; dict si names, sinon ndarray."""
        if self.model is None:
            raise RuntimeError("appeler fit d abord")
        if self.kind == "gp":
            ls = np.atleast_1d(self.model.kernel_.k1.k2.length_scale)
            imp = 1.0 / ls
        else:
            pi = permutation_importance(self.model, self.X_, self.y_, n_repeats=n_repeats, random_state=self.seed)
            imp = np.maximum(pi.importances_mean, 0.0)
        s = imp.sum()
        imp = imp / s if s > 0 else np.full_like(imp, 1.0 / len(imp))
        return dict(zip(self.names, imp.tolist())) if self.names else imp

    def describe(self):
        if self.kind == "gp":
            return f"GP n={len(self.y_)} (ecartes {self.n_dropped}) noyau {self.model.kernel_} nugget_sd={self.nugget_std:.4g}"
        return f"RF n={len(self.y_)} (ecartes {self.n_dropped}) arbres {self.n_trees} feuille>={self.min_samples_leaf}"


class MultiEmulator:
    """Un Emulator par observable ; memes X, y sous forme de dict observable -> (n,)."""

    def __init__(self, kind="gp", names=None, **kw):
        self.kind = kind
        self.names = list(names) if names is not None else None
        self.kw = kw
        self.models = {}

    def fit(self, X, ydict, noise=None):
        """ydict : {observable: y (n,)} ; noise : None ou {observable: ecart-type (n,)}."""
        self.models = {}
        for k, y in ydict.items():
            nz = None if noise is None else noise.get(k)
            self.models[k] = Emulator(self.kind, names=self.names, **self.kw).fit(X, y, nz)
        return self

    def keys(self):
        return list(self.models.keys())

    def __getitem__(self, k):
        return self.models[k]

    def __contains__(self, k):
        return k in self.models

    def predict(self, X, include_noise=False):
        """-> {observable: (mean, std)}."""
        return {k: m.predict(X, include_noise=include_noise) for k, m in self.models.items()}

    def predict_point(self, u):
        """Un seul point (d,) -> {observable: mean (float)} - pratique pour l objectif."""
        u = np.atleast_2d(np.asarray(u, dtype=float))
        return {k: float(m.predict(u)[0][0]) for k, m in self.models.items()}

    def loo_rmse(self, **kw):
        return {k: m.loo_rmse(**kw) for k, m in self.models.items()}

    def importance(self, **kw):
        return {k: m.importance(**kw) for k, m in self.models.items()}

    def describe(self):
        return "\n".join(f"{k:12s} {m.describe()}" for k, m in self.models.items())


# ---------------------------------------------------------------------------
# Tests sur fonctions analytiques
# ---------------------------------------------------------------------------
def _test():
    import time
    from scipy.stats import qmc

    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  [OK]  " if cond else "  [FAIL] ") + msg)
        ok = ok and bool(cond)

    rng = np.random.default_rng(0)
    d = 3
    names = ["ft", "c", "inerte"]

    # fonction test : lisse, 2 parametres actifs, 1 inerte, bruit gaussien 0.05 (nugget)
    def f_true(U):
        U = np.atleast_2d(U)
        return np.sin(2 * np.pi * U[:, 0]) + 2.0 * (U[:, 1] - 0.3) ** 2

    sd_noise = 0.05
    n = 40
    X = qmc.LatinHypercube(d=d, scramble=True, seed=1).random(n)
    y = f_true(X) + rng.normal(0, sd_noise, n)
    Xt = qmc.Sobol(d=d, scramble=True, seed=2).random(256)
    yt = f_true(Xt)

    # 1. GP
    print("1. GP (Matern 5/2 ARD + nugget) sur f = sin(2 pi u0) + 2 (u1 - 0.3)^2 + N(0, 0.05^2), n = 40")
    t0 = time.time()
    gp = Emulator("gp", names=names, seed=0).fit(X, y)
    print(f"     {gp.describe()}  ({time.time() - t0:.1f} s)")
    mu, sd = gp.predict(Xt)
    rmse_t = float(np.sqrt(np.mean((mu - yt) ** 2)))
    print(f"     RMSE test (256 pts Sobol, sans bruit) = {rmse_t:.4f}  | nugget appris = {gp.nugget_std:.4f} (vrai 0.05)")
    check(rmse_t < 0.10, "RMSE test < 0.10 (amplitude ~ 2)")
    check(0.02 < gp.nugget_std < 0.12, "nugget appris dans [0.02, 0.12] autour de 0.05")
    mu_tr, sd_tr = gp.predict(X)
    mu_n, sd_n = gp.predict(X, include_noise=True)
    check(np.all(sd_n >= sd_tr - 1e-12) and np.mean(sd_n) > np.mean(sd_tr), "std avec bruit >= std latente")
    check(np.mean(sd_tr) < np.mean(sd), "std latente plus petite aux points d entrainement qu en test")
    t0 = time.time()
    rmse_loo, yhat = gp.loo_rmse(return_pred=True)
    print(f"     LOO RMSE (refit) = {rmse_loo:.4f}, Q2 = {1 - np.sum((yhat - y) ** 2) / np.sum((y - y.mean()) ** 2):.3f}  ({time.time() - t0:.1f} s)")
    check(rmse_loo < 0.15, "LOO RMSE < 0.15")
    rmse_fix = gp.loo_rmse(refit_hyper=False)
    print(f"     LOO RMSE (hyperparametres figes) = {rmse_fix:.4f}")
    check(abs(rmse_fix - rmse_loo) < 0.1, "LOO fige et LOO refit coherents")
    imp = gp.importance()
    print("     importance GP :", {k: round(v, 3) for k, v in imp.items()})
    check(imp["inerte"] == min(imp.values()) and imp["inerte"] < 0.15, "parametre inerte = importance minimale (< 0.15)")
    check(abs(sum(imp.values()) - 1) < 1e-9, "importances normalisees a 1")

    # 2. bruit connu par point
    print("2. GP avec bruit connu (noise=ecart-type par point)")
    gp2 = Emulator("gp", names=names, seed=0).fit(X, y, noise=np.full(n, sd_noise))
    mu2, _ = gp2.predict(Xt)
    rmse2 = float(np.sqrt(np.mean((mu2 - yt) ** 2)))
    print(f"     RMSE test = {rmse2:.4f}, nugget residuel appris = {gp2.nugget_std:.4f}")
    check(rmse2 < 0.12, "RMSE test < 0.12 avec bruit impose")

    # 3. NaN
    print("3. valeurs manquantes")
    y_nan = y.copy(); y_nan[[3, 7, 11]] = np.nan
    gp3 = Emulator("gp", seed=0, n_restarts=2).fit(X, y_nan)
    check(gp3.n_dropped == 3 and len(gp3.y_) == n - 3, "3 lignes NaN ecartees a l ajustement")
    mu3, _ = gp3.predict(Xt)
    check(np.isfinite(mu3).all(), "predictions finies apres NaN")
    try:
        Emulator("gp").fit(X * 3.0, y)
        check(False, "X hors [0, 1] doit lever ValueError")
    except ValueError:
        check(True, "X hors [0, 1] refuse (ValueError)")

    # 4. RF
    print("4. RF (400 arbres, feuille >= 2)")
    t0 = time.time()
    rf = Emulator("rf", names=names, seed=0, n_trees=200).fit(X, y)
    mu_r, sd_r = rf.predict(Xt)
    rmse_r = float(np.sqrt(np.mean((mu_r - yt) ** 2)))
    print(f"     {rf.describe()}  RMSE test = {rmse_r:.4f}  std inter-arbres moy = {sd_r.mean():.3f}  ({time.time() - t0:.1f} s)")
    check(rmse_r < 0.45, "RMSE test RF < 0.45 (plus grossier que le GP sur une fonction lisse)")
    check(rmse_t < rmse_r, "GP meilleur que RF sur fonction lisse")
    t0 = time.time()
    rmse_r_loo = rf.loo_rmse()
    print(f"     LOO RMSE RF = {rmse_r_loo:.4f}  ({time.time() - t0:.1f} s)")
    check(rmse_r_loo < 0.6, "LOO RMSE RF < 0.6")
    imp_r = rf.importance()
    print("     importance RF (permutation) :", {k: round(v, 3) for k, v in imp_r.items()})
    check(imp_r["inerte"] == min(imp_r.values()), "parametre inerte = importance minimale (RF)")
    check(np.all(sd_r >= 0), "std RF >= 0")

    # 5. MultiEmulator
    print("5. MultiEmulator (2 observables, dont une avec NaN)")
    y2 = np.exp(-3 * X[:, 2]) + rng.normal(0, 0.02, n)
    y2[[0, 5]] = np.nan
    me = MultiEmulator("gp", names=names, seed=0, n_restarts=2).fit(X, {"q_peak": y, "d_CD": y2})
    pr = me.predict(Xt)
    check(set(pr) == {"q_peak", "d_CD"} and pr["d_CD"][0].shape == (256,), "predict -> dict observable -> (mean, std)")
    check(me["d_CD"].n_dropped == 2, "NaN ecartees par observable")
    imp2 = me.importance()["d_CD"]
    check(max(imp2, key=imp2.get) == "inerte", "d_CD depend du 3e parametre : importance max sur 'inerte'")
    p = me.predict_point([0.25, 0.3, 0.5])
    check(abs(p["q_peak"] - 1.0) < 0.2, f"predict_point q_peak({{0.25, 0.3, 0.5}}) = {p['q_peak']:.3f} ~ 1.0")
    loo = me.loo_rmse(refit_hyper=False)
    print("     LOO (figes) :", {k: round(v, 4) for k, v in loo.items()})
    check(all(np.isfinite(list(loo.values()))), "LOO multi finies")

    # 6. KFold pour n > 60
    print("6. n = 80 -> KFold 10 plis")
    X6 = qmc.LatinHypercube(d=d, scramble=True, seed=5).random(80)
    y6 = f_true(X6) + rng.normal(0, sd_noise, 80)
    gp6 = Emulator("gp", seed=0, n_restarts=2).fit(X6, y6)
    t0 = time.time()
    r6 = gp6.loo_rmse()
    print(f"     KFold RMSE = {r6:.4f}  ({time.time() - t0:.1f} s)")
    check(r6 < 0.12, "KFold RMSE < 0.12")

    print("\nRESULTAT emulator.py :", "tous les tests passent" if ok else "ECHEC")
    return ok


if __name__ == "__main__":
    if "--test" in sys.argv:
        sys.exit(0 if _test() else 1)
    print(__doc__)
