#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# objective.py - fonction objectif de la calibration triaxiale rockim :
# combine les ecarts d extract.observables() sur les confinements {20, 50} MPa.
#
#   python calib_quick/calib/objective.py --test
#   python calib_quick/calib/objective.py out_a_P020 out_a_P050     (J d un jeu de runs)
# ---------------------------------------------------------------------------
"""Objectif de calibration : somme ponderee de residus normalises au carre.

Pour chaque confinement sigma3 dans `confinements` (defaut 20 et 50 MPa) et chaque terme k :
    J_k = w_k * (d_k / tol_k)^2
avec les observables relatives d extract.py (d_peak = q_pic/q_pic_exp - 1, d_eps_peak,
d_CI, d_CD, d_drop = chute - chute_exp) et rmse_curve (MPa, sur la grille des cibles) :

    terme        tolerance   poids   sens
    d_peak         0.03       1.0    pic a 3 %
    d_eps_peak     0.10       0.5    deformation au pic a 10 %
    d_CI           0.10       1.0    seuil d initiation (1re insertion adaptative) a 10 %
    d_CD           0.10       0.5    seuil de croissance instable (1re rupture) a 10 %
    d_drop         0.10       0.5    chute post-pic a 0.10 (absolu, fraction du pic)
    rmse_curve    15 MPa      1.0    ecart de forme sur la montee

Une tolerance = 'ecart acceptable' : un terme vaut w_k quand l ecart egale sa tolerance.
J = somme sur les confinements et les termes (un chi2 pondere, sans racine) ; J ~ 1-2 par
confinement signifie 'dans les tolerances'.

Penalites :
* run non rompu (nb_end == 0, ou nb_end absent) : + unbroken_penalty (defaut 10) - un run qui
  ne casse pas n a ni CD ni chute et ne doit jamais gagner ;
* terme non fini (NaN : pas d insertion -> d_CI NaN, pas de rupture -> d_CD NaN, courbe trop
  courte -> rmse NaN) : w_k * nan_penalty (defaut 9 = un ecart de 3 tolerances) ;
* confinement absent du dict (run plante) : tous les termes NaN + penalite non rompu.

Interface
---------
objective(obs_by_conf, weights=None, tol=None, confinements=(20, 50), ...) -> (J, detail)
    obs_by_conf : {20: obs, 50: obs} (cles int, float ou str acceptees), obs = dict
    d extract.observables (ou toute source qui fournit les memes cles).
    detail : {'20': {'terms': {k: J_k}, 'resid': {k: d_k / tol_k}, 'unbroken': p, 'missing': b,
              'sum': J_conf}, '50': ..., 'total': J}.
observables_by_confinement(out_dirs) : construit obs_by_conf a partir de dossiers de runs
    (import paresseux d extract.py, qui lit history.csv et les cibles).

Dependances : numpy uniquement (extract.py seulement pour la lecture de runs).
"""
import json
import math
import os
import sys

import numpy as np

TERMS = ("d_peak", "d_eps_peak", "d_CI", "d_CD", "d_drop", "rmse_curve")
TOL_DEFAULT = {"d_peak": 0.03, "d_eps_peak": 0.10, "d_CI": 0.10, "d_CD": 0.10, "d_drop": 0.10, "rmse_curve": 15.0}
W_DEFAULT = {"d_peak": 1.0, "d_eps_peak": 0.5, "d_CI": 1.0, "d_CD": 0.5, "d_drop": 0.5, "rmse_curve": 1.0}
CONF_DEFAULT = (20, 50)


def _conf_key(k):
    """Normalise une cle de confinement (20, 20.0, '20', '020', 'P020') -> '20'."""
    s = str(k).strip()
    if s.upper().startswith("P"):
        s = s[1:]
    try:
        return str(int(round(float(s))))
    except ValueError:
        return s


def _lookup(obs_by_conf, conf):
    want = _conf_key(conf)
    for k, v in obs_by_conf.items():
        if _conf_key(k) == want:
            return v
    return None


def _finite(v):
    try:
        return v is not None and math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


BTS_TOL = 0.10        # tolerance relative sur le BTS (exp 10,27 +- 0,98 MPa)
W_BTS = 1.0


def objective_with_bts(obs_by_conf, bts_obs=None, weights=None, tol=None, confinements=CONF_DEFAULT,
                       bts_tol=BTS_TOL, w_bts=W_BTS, **kw):
    """objective() + terme BTS (confinement-independant) : w_bts (d_bts / bts_tol)^2 ;
    bts_obs = dict de extract.bts_observables (cle d_bts) ou None (terme absent)."""
    J, detail = objective(obs_by_conf, weights=weights, tol=tol, confinements=confinements, **kw)
    if bts_obs is not None:
        d = bts_obs.get("d_bts", float("nan"))
        term = w_bts * (d / bts_tol) ** 2 if (d == d) else w_bts * 9.0
        detail["BTS"] = {"terms": {"d_bts": term}, "resid": {"d_bts": d / bts_tol if d == d else float("nan")}, "sum": term}
        J += term
        detail["total"] = J
    return J, detail


def objective(obs_by_conf, weights=None, tol=None, confinements=CONF_DEFAULT,
              unbroken_penalty=10.0, nan_penalty=9.0):
    """-> (J, detail). Voir la docstring du module."""
    w = dict(W_DEFAULT); w.update(weights or {})
    t = dict(TOL_DEFAULT); t.update(tol or {})
    bad = [k for k in TERMS if not (t[k] > 0)]
    if bad:
        raise ValueError(f"tolerances non positives : {bad}")
    detail = {}
    total = 0.0
    for conf in confinements:
        obs = _lookup(obs_by_conf, conf)
        terms, resid = {}, {}
        missing = obs is None
        for k in TERMS:
            v = None if missing else obs.get(k)
            if _finite(v):
                r = float(v) / t[k]
                resid[k] = r
                terms[k] = w[k] * r * r
            else:
                resid[k] = float("nan")
                terms[k] = w[k] * nan_penalty
        if missing:
            unbroken = unbroken_penalty
        else:
            nb_end = obs.get("nb_end")
            unbroken = unbroken_penalty if (not _finite(nb_end) or float(nb_end) <= 0) else 0.0
        s = float(sum(terms.values()) + unbroken)
        detail[_conf_key(conf)] = {"terms": terms, "resid": resid, "unbroken": unbroken, "missing": missing, "sum": s}
        total += s
    detail["total"] = float(total)
    return float(total), detail


def format_detail(detail):
    """Tableau lisible des termes par confinement."""
    confs = [k for k in detail if k != "total"]
    lines = [f"{'terme':12s} " + " ".join(f"{'P' + c:>16s}" for c in confs)]
    for k in TERMS:
        cells = []
        for c in confs:
            d = detail[c]
            cells.append(f"{d['terms'][k]:8.3f} ({d['resid'][k]:+6.2f})")
        lines.append(f"{k:12s} " + " ".join(f"{x:>16s}" for x in cells))
    lines.append(f"{'non rompu':12s} " + " ".join(f"{detail[c]['unbroken']:>16.1f}" for c in confs))
    lines.append(f"{'somme':12s} " + " ".join(f"{detail[c]['sum']:>16.3f}" for c in confs))
    lines.append(f"{'TOTAL':12s} {detail['total']:>16.3f}      (terme (residu / tolerance))")
    return "\n".join(lines)


def observables_by_confinement(out_dirs):
    """Dossiers de runs -> {s3 (int MPa): obs} via extract.py (lecture seule, aucun run)."""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from extract import TARGETS, load_run, observables  # noqa: E402
    targets = json.load(open(TARGETS)) if os.path.exists(TARGETS) else None
    out = {}
    for d in out_dirs:
        if not os.path.exists(os.path.join(d, "history.csv")):
            continue
        r = load_run(d)
        out[int(round(r["s3"]))] = observables(r, targets)
    return out


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def _test():
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  [OK]  " if cond else "  [FAIL] ") + msg)
        ok = ok and bool(cond)

    perfect = {k: 0.0 for k in TERMS}
    perfect["nb_end"] = 12
    print("1. run parfait sur les deux confinements")
    J, det = objective({20: dict(perfect), 50: dict(perfect)})
    check(J == 0.0, f"J = {J}")
    check(set(det) == {"20", "50", "total"} and det["20"]["sum"] == 0.0, "detail par confinement + total")

    print("2. un ecart d une tolerance par terme -> J = poids")
    for k in TERMS:
        o = dict(perfect); o[k] = TOL_DEFAULT[k]
        J, det = objective({20: o, 50: dict(perfect)})
        check(abs(J - W_DEFAULT[k]) < 1e-12 and abs(det["20"]["resid"][k] - 1.0) < 1e-12,
              f"{k} = tol -> J = {J:.3f} = w ({W_DEFAULT[k]}), residu 1")
    o = dict(perfect); o["d_peak"] = -0.06                    # -2 tolerances
    J, det = objective({20: o, 50: dict(perfect)})
    check(abs(J - 4.0) < 1e-12 and det["20"]["resid"]["d_peak"] == -2.0, "d_peak = -2 tol -> J = 4, residu signe -2")

    print("3. tolerances et poids en arguments")
    o = dict(perfect); o["rmse_curve"] = 30.0
    J1, _ = objective({20: o, 50: dict(perfect)})
    J2, _ = objective({20: o, 50: dict(perfect)}, tol={"rmse_curve": 30.0})
    J3, _ = objective({20: o, 50: dict(perfect)}, weights={"rmse_curve": 0.25})
    check(abs(J1 - 4.0) < 1e-12 and abs(J2 - 1.0) < 1e-12 and abs(J3 - 1.0) < 1e-12,
          f"rmse 30 MPa : J = {J1} (tol 15) ; {J2} (tol 30) ; {J3} (poids 0.25)")
    try:
        objective({20: o, 50: o}, tol={"d_peak": 0.0})
        check(False, "tolerance nulle doit lever ValueError")
    except ValueError:
        check(True, "tolerance nulle refusee")

    print("4. penalites")
    o = dict(perfect); o["nb_end"] = 0
    J, det = objective({20: o, 50: dict(perfect)})
    check(J == 10.0 and det["20"]["unbroken"] == 10.0, f"non rompu (nb_end = 0) -> +10 : J = {J}")
    J, _ = objective({20: o, 50: dict(perfect)}, unbroken_penalty=100.0)
    check(J == 100.0, "penalite non rompu en argument")
    o = dict(perfect); o["d_CI"] = float("nan")
    J, det = objective({20: o, 50: dict(perfect)})
    check(abs(J - 9.0 * W_DEFAULT["d_CI"]) < 1e-12 and math.isnan(det["20"]["resid"]["d_CI"]), f"d_CI NaN -> w * 9 : J = {J}")
    o = dict(perfect); del o["d_CD"]
    J, _ = objective({20: o, 50: dict(perfect)})
    check(abs(J - 9.0 * W_DEFAULT["d_CD"]) < 1e-12, "cle absente traitee comme NaN")
    o = dict(perfect); del o["nb_end"]
    J, _ = objective({20: o, 50: dict(perfect)})
    check(J == 10.0, "nb_end absent -> penalite non rompu (prudence)")
    J, det = objective({20: dict(perfect)})
    n_missing = 9.0 * sum(W_DEFAULT.values()) + 10.0
    check(abs(J - n_missing) < 1e-12 and det["50"]["missing"], f"confinement 50 absent -> {n_missing:.1f} (tous NaN + non rompu)")

    print("5. cles de confinement et liste de confinements")
    J1, _ = objective({"20": dict(perfect), "P050": dict(perfect)})
    J2, _ = objective({20.0: dict(perfect), 50: dict(perfect)})
    check(J1 == 0.0 and J2 == 0.0, "cles '20', 'P050', 20.0 reconnues")
    o = dict(perfect); o["d_peak"] = 0.03
    J, det = objective({20: dict(perfect), 50: dict(perfect), 75: o}, confinements=(20, 50, 75))
    check(abs(J - 1.0) < 1e-12 and "75" in det, "confinements=(20, 50, 75)")
    J, _ = objective({20: dict(perfect), 50: dict(perfect), 75: o})
    check(J == 0.0, "75 ignore par defaut")

    print("6. cas realiste (chiffres du type extract.py) et tableau")
    real20 = dict(d_peak=-0.045, d_eps_peak=0.12, d_CI=0.08, d_CD=-0.15, d_drop=0.05, rmse_curve=22.0, nb_end=31)
    real50 = dict(d_peak=0.02, d_eps_peak=-0.04, d_CI=float("nan"), d_CD=0.10, d_drop=-0.12, rmse_curve=11.0, nb_end=18)
    J, det = objective({20: real20, 50: real50})
    expected = (W_DEFAULT["d_peak"] * (0.045 / 0.03) ** 2 + 0.5 * 1.2 ** 2 + 1.0 * 0.8 ** 2 + 0.5 * 1.5 ** 2
                + 0.5 * 0.5 ** 2 + (22 / 15) ** 2
                + 1.0 * (0.02 / 0.03) ** 2 + 0.5 * 0.4 ** 2 + 9.0 + 0.5 * 1.0 ** 2 + 0.5 * 1.2 ** 2 + (11 / 15) ** 2)
    check(abs(J - expected) < 1e-9, f"J = {J:.4f} = somme calculee a la main {expected:.4f}")
    print(format_detail(det))

    print("\nRESULTAT objective.py :", "tous les tests passent" if ok else "ECHEC")
    return ok


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if "--test" in sys.argv:
        sys.exit(0 if _test() else 1)
    dirs = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not dirs:
        print(__doc__)
        return
    obs = observables_by_confinement(dirs)
    if not obs:
        raise SystemExit("aucun history.csv trouve")
    J, det = objective(obs)
    print(format_detail(det))
    if "--json" in sys.argv:
        print(json.dumps(det, indent=1))


if __name__ == "__main__":
    main()
