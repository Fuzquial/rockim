#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# campaign.py - orchestration de la calibration par substitut (calib_quick) :
# collecte des runs -> table de base -> GP par observable -> MAP par APSO ->
# enrichissement LCB / Kriging believer -> rapport.
#
#   python calib_quick/calib/campaign.py init     --tag h1 --template calib_quick/q1v070_P050.cfg
#                                                 --bts calib_quick/bts_v070b.cfg --extra k=v ...
#   python calib_quick/calib/campaign.py collect  --tag h1 [--runs calib_quick/runs_h1] [--weights w.json]
#   python calib_quick/calib/campaign.py fit      --tag h1 [--rf] [--fast]
#   python calib_quick/calib/campaign.py optimize --tag h1 [--n-seeds 20 --particles 40 --iters 300]
#   python calib_quick/calib/campaign.py enrich   --tag h1 --n 8 --kappa 2
#   python calib_quick/calib/campaign.py report   --tag h1
#   python calib_quick/calib/campaign.py --test
#
# NE LANCE AUCUN run rockim : lecture des sorties, decks via make_decks.py.
# ---------------------------------------------------------------------------
"""Campagne de calibration assistee par substitut (PLAN_calibration.md 4a, A4 4-5-6-8).

Fichiers d une campagne <tag> (tous dans le dossier des runs, defaut calib_quick/runs_<tag>) :
  campaign_<tag>.json   configuration (template, bts, extra, conf, poids, nu) ecrite par `init`
  space_<tag>.json      espace {nom: [lo, hi, 'lin'|'log']} (cree depuis PLAN 4a si absent)
  design_<tag>.csv      plan initial (ParamSpace.write_csv) ; design_<tag>_r<k>.csv = lots d enrichissement
  jobs_<tag>.json       liste des runs (make_decks.py) ; jobs_<tag>_r<k>.json pour les lots
  base_<tag>.csv        table de base : une ligne par point (parametres, observables par
                        confinement P020_/P050_/BTS_, J et detail par terme, complete)
  models_<tag>.pkl      emulateurs GP (joblib) + classifieurs + LOO + importances
  map_<tag>.json        optimum du substitut (APSO x n graines), top 5 distincts
  enrich_<tag>_r<k>.json  trace du lot k (u, LCB, moyenne, ecart-type, distance)
  report_<tag>.png/.md  figure et tableaux
Cles optionnelles de campaign_<tag>.json : weights (JSON poids/tolerances, voir
objective_settings), nu (defaut 0.29), nb_min (defaut 1 : seuil 'rompu' en joints rompus),
targets / seuils (chemins des cibles si differents d extract.py).

Observables emulees (une GP par observable et par confinement, en coordonnees [0, 1]^d) :
  lq_peak = log(q_peak)        eps_peak (%)      CI_frac      CD_frac      drop
  lrmse   = log1p(rmse_curve)  lnb_end = log1p(nb_end)       BTS : lbts = log(bts)
Les transformations log rendent les reponses plus proches d une gaussienne stationnaire
(pic et BTS varient d un facteur 3 sur l espace, rmse de 5 a 200 MPa) ; les fractions et la
chute restent lineaires. Un run NON ROMPU (nb_end < nb_min, defaut 1 : le modele n a pas casse
dans la fenetre imposee) donne q_peak = q_end, CI/CD NaN, chute 0 : c est une discontinuite
qu un GP stationnaire ne sait pas suivre (A4 4.2). Traitement : (i) les GP de regression sont
ajustes sur les seuls runs rompus (les autres sont mis a NaN et ecartes par emulator.py, pas
imputes) ; (ii) si au moins 2 runs sont bloques a un confinement, un classifieur GP
rompu/bloque (GaussianProcessClassifier, Matern ARD) est ajuste et le substitut impose
nb_end = 0 (penalite 'non rompu' de objective.py) la ou P(rompu) < 0,5. lnb_end est emule sur
les runs rompus pour rendre le nombre de joints rompus.

Objectif emule J_hat(u) : les observables relatives (d_peak, d_eps_peak, d_CI, d_CD, d_drop)
sont RECONSTRUITES a partir des predictions inverses-transformees et des cibles
experimentales (targets_triax_bohus.json, seuils_sbm_bohus.json, BTS_EXP d extract.py),
puis passees a objective.objective_with_bts (poids et tolerances de la campagne). Une version
numpy vectorisee (j_numpy) sert a l optimiseur et a la propagation Monte-Carlo ; son egalite
avec objective_with_bts est verifiee a chaque `optimize` et dans --test.

Enrichissement (A4 5.4, 8.7) : LCB(u) = mean_J - kappa std_J, moments par 200 tirages des GP
(normales independantes par observable, tirage commun a tous les u pour une acquisition
deterministe), lots par Kriging believer (la prediction moyenne devient pseudo-observation,
re-ajustement a hyperparametres figes), distance minimale 0,15 en coordonnees normalisees
entre nouveaux points et points existants (penalite dans l acquisition + verification).

Dependances : numpy, scipy, scikit-learn, joblib, matplotlib (report) + modules du dossier.
"""
import argparse
import copy
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys
import time
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from design import ParamSpace, min_pairwise_distance  # noqa: E402
from emulator import Emulator, MultiEmulator  # noqa: E402
from apso import apso  # noqa: E402
from objective import (objective_with_bts, format_detail, TERMS, TOL_DEFAULT, W_DEFAULT,  # noqa: E402
                       BTS_TOL, W_BTS)

# espace par defaut : PLAN_calibration.md 4a (homogene, 5 parametres)
SPACE_DEFAULT = {
    "ft": [4e6, 20e6, "log"],
    "c": [8e6, 40e6, "log"],
    "phi": [35.0, 55.0, "lin"],
    "Gf": [10.0, 60.0, "log"],
    "gII": [2.0, 30.0, "log"],
}
CONF_DEFAULT = [20, 50]
NU_DEFAULT = 0.29          # nu des decks (E physique) : cibles eps x (1 - nu^2)

# colonnes conservees dans base_<tag>.csv par confinement / pour le bresilien
TRIAX_COLS = ["q_peak", "eps_peak", "E_GPa", "drop", "q_CI", "q_CD", "CI_frac", "CD_frac", "ni_peak",
              "nb_peak", "nb_end", "rmse_curve", "rmse_curve_norm", "d_peak", "d_eps_peak", "d_drop",
              "d_CI", "d_CD"]
BTS_COLS = ["bts", "d_bts", "bts_nominal", "k_band", "nb_end"]


# ---------------------------------------------------------------------------
# Transformations des observables
# ---------------------------------------------------------------------------
def _log_pos(x):
    x = np.asarray(x, dtype=float)
    out = np.full(x.shape, np.nan)
    m = np.isfinite(x) & (x > 0)
    out[m] = np.log(x[m])
    return out


def _log1p_pos(x):
    x = np.asarray(x, dtype=float)
    out = np.full(x.shape, np.nan)
    m = np.isfinite(x) & (x > -1)
    out[m] = np.log1p(x[m])
    return out


def _ident(x):
    return np.asarray(x, dtype=float)


TRANSFORMS = {"id": (_ident, _ident), "log": (_log_pos, np.exp), "log1p": (_log1p_pos, np.expm1)}
# (cle extract, nom emule, transformation)
OBS_TRIAX = [("q_peak", "lq_peak", "log"), ("eps_peak", "eps_peak", "id"), ("CI_frac", "CI_frac", "id"),
             ("CD_frac", "CD_frac", "id"), ("drop", "drop", "id"), ("rmse_curve", "lrmse", "log1p"),
             ("nb_end", "lnb_end", "log1p")]
OBS_BTS = [("bts", "lbts", "log")]


# ---------------------------------------------------------------------------
# Chemins, configuration, espace
# ---------------------------------------------------------------------------
def conf_prefix(conf):
    """20 / 20.0 / '20' -> 'P020' ; 'BTS' -> 'BTS'."""
    if isinstance(conf, str) and conf.strip().upper() == "BTS":
        return "BTS"
    return f"P{int(round(float(conf))):03d}"


def prefix_conf(prefix):
    """'P020' -> 20 ; 'BTS' -> 'BTS'."""
    return "BTS" if prefix == "BTS" else int(prefix[1:])


def resolve_path(p):
    """Chemin absolu, relatif au cwd s il existe, sinon relatif a ROOT (rockim_f2)."""
    if p is None:
        return None
    if os.path.isabs(p):
        return p
    if os.path.exists(p):
        return os.path.abspath(p)
    return os.path.join(ROOT, p)


def resolve_runs(tag, runs):
    if runs is None:
        return os.path.join(ROOT, "calib_quick", f"runs_{tag}")
    return resolve_path(runs)


class Paths:
    """Fichiers d une campagne dans son dossier de runs."""

    def __init__(self, tag, runs):
        self.tag = tag
        self.runs = runs
        j = lambda name: os.path.join(runs, name)  # noqa: E731
        self.cfg = j(f"campaign_{tag}.json")
        self.space = j(f"space_{tag}.json")
        self.base = j(f"base_{tag}.csv")
        self.models = j(f"models_{tag}.pkl")
        self.map = j(f"map_{tag}.json")
        self.report_png = j(f"report_{tag}.png")
        self.report_md = j(f"report_{tag}.md")
        self.design0 = j(f"design_{tag}.csv")

    def designs(self):
        """[(round, chemin)] : design_<tag>.csv (0) puis design_<tag>_r<k>.csv tries par k."""
        out = [(0, self.design0)] if os.path.exists(self.design0) else []
        pat = re.compile(re.escape(f"design_{self.tag}_r") + r"(\d+)\.csv$")
        for p in glob.glob(os.path.join(self.runs, f"design_{self.tag}_r*.csv")):
            m = pat.search(os.path.basename(p))
            if m:
                out.append((int(m.group(1)), p))
        return sorted(out)

    def jobs_files(self):
        out = []
        p0 = os.path.join(self.runs, f"jobs_{self.tag}.json")
        if os.path.exists(p0):
            out.append((0, p0))
        pat = re.compile(re.escape(f"jobs_{self.tag}_r") + r"(\d+)\.json$")
        for p in glob.glob(os.path.join(self.runs, f"jobs_{self.tag}_r*.json")):
            m = pat.search(os.path.basename(p))
            if m:
                out.append((int(m.group(1)), p))
        return [p for _, p in sorted(out)]

    def next_round(self):
        rounds = [k for k, _ in self.designs()]
        return (max(rounds) + 1) if rounds else 1

    def design_round(self, k):
        return os.path.join(self.runs, f"design_{self.tag}_r{k}.csv")


def load_config(P):
    if os.path.exists(P.cfg):
        return json.load(open(P.cfg))
    return {}


def save_config(P, cfg):
    os.makedirs(P.runs, exist_ok=True)
    json.dump(cfg, open(P.cfg, "w"), indent=1)


def load_space(P, create=True):
    """ParamSpace depuis space_<tag>.json ; cree le fichier depuis PLAN 4a si absent."""
    if os.path.exists(P.space):
        spec = json.load(open(P.space))
    else:
        if not create:
            raise FileNotFoundError(P.space)
        spec = dict(SPACE_DEFAULT)
        os.makedirs(P.runs, exist_ok=True)
        json.dump(spec, open(P.space, "w"), indent=1)
        print(f"[space] {P.space} cree depuis PLAN 4a : " + ", ".join(f"{k} [{v[0]:g}, {v[1]:g}] {v[2]}" for k, v in spec.items()))
    return ParamSpace({k: tuple(v) for k, v in spec.items()}), spec


def objective_settings(cfg, weights_path=None):
    """Poids / tolerances / penalites de la campagne -> dict pour objective_with_bts.

    JSON optionnel : {"weights": {terme: w}, "tol": {terme: tol}, "bts_tol": 0.1, "w_bts": 1.0,
    "unbroken_penalty": 10, "nan_penalty": 9, "confinements": [20, 50]}. Un JSON plat dont les
    cles sont des termes est lu comme "weights"."""
    S = dict(weights={}, tol={}, confinements=[int(c) for c in cfg.get("conf", CONF_DEFAULT)],
             bts_tol=BTS_TOL, w_bts=W_BTS, unbroken_penalty=10.0, nan_penalty=9.0)
    src = weights_path or cfg.get("weights")
    if src:
        raw = json.load(open(resolve_path(src))) if isinstance(src, str) else dict(src)
        if any(k in TERMS for k in raw) and "weights" not in raw:
            raw = {"weights": {k: v for k, v in raw.items() if k in TERMS}}
        for k in S:
            if k in raw:
                S[k] = raw[k]
        S["confinements"] = [int(c) for c in S["confinements"]]
    return S


def eval_objective(obs_by_conf, bts_obs, S):
    return objective_with_bts(obs_by_conf, bts_obs, weights=S["weights"], tol=S["tol"],
                              confinements=tuple(S["confinements"]), bts_tol=S["bts_tol"], w_bts=S["w_bts"],
                              unbroken_penalty=S["unbroken_penalty"], nan_penalty=S["nan_penalty"])


def format_detail_bts(detail):
    """format_detail d objective.py + ligne BTS."""
    core = {k: v for k, v in detail.items() if k != "BTS"}
    txt = format_detail(core)
    if "BTS" in detail:
        b = detail["BTS"]
        txt += f"\n{'BTS d_bts':12s} {b['terms']['d_bts']:8.3f} ({b['resid']['d_bts']:+6.2f})"
    return txt


# ---------------------------------------------------------------------------
# Cibles experimentales (memes sources qu extract.py)
# ---------------------------------------------------------------------------
def _seuils_from(path):
    raw = json.load(open(path))
    acc = {}
    for t in raw.values():
        k = str(int(round(t["sigma3"])))
        acc.setdefault(k, []).append((t["CI_MPa"], t["CD_MPa"], t["q_peak"]))
    return {k: dict(q_CI=float(np.mean([v[0] for v in a])), q_CD=float(np.mean([v[1] for v in a])),
                    CI_frac=float(np.mean([v[0] / v[2] for v in a])), CD_frac=float(np.mean([v[1] / v[2] for v in a])))
            for k, a in acc.items()}


def exp_targets(nu=NU_DEFAULT, targets_path=None, seuils_path=None):
    """{'20': {q_peak, eps_peak (fraction, x (1-nu^2)), chute, q_CI, q_CD, CI_frac, CD_frac,
    eps_grid, q_mean, q_std}, '50': ..., 'BTS': BTS_EXP}."""
    import extract
    tp = targets_path or extract.TARGETS
    sp = seuils_path or extract.SEUILS
    if not os.path.exists(tp):
        raise FileNotFoundError(f"cibles introuvables : {tp}")
    T = json.load(open(tp))
    f = 1.0 - nu * nu
    seuils = _seuils_from(sp) if os.path.exists(sp) else {}
    out = {}
    for key, tg in T["confinements"].items():
        d = dict(q_peak=float(tg["q_peak_mean_MPa"]), eps_peak=float(tg["eps_peak_microstrain"]) * 1e-6 * f,
                 chute=float(tg["chute_fraction_moyenne"]),
                 eps_grid=np.array(tg["eps_grid_microstrain"], dtype=float) * 1e-6 * f,
                 q_mean=np.array(tg["q_mean_MPa"], dtype=float), q_std=np.array(tg["q_std_MPa"], dtype=float))
        if key in seuils:
            d.update(seuils[key])
        out[key] = d
    out["BTS"] = float(extract.BTS_EXP)
    out["_nu"] = float(nu)
    return out


# ---------------------------------------------------------------------------
# Tables CSV (base_<tag>.csv)
# ---------------------------------------------------------------------------
def _fmt(v):
    if isinstance(v, str):
        return v
    if v is None:
        return "nan"
    if isinstance(v, (bool, np.bool_)):
        return str(int(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    return f"{f:.10g}" if math.isfinite(f) else "nan"


def write_table(path, header, rows):
    with open(path, "w", newline="", encoding="ascii") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow([_fmt(r.get(h)) for h in header])


def read_table(path):
    """-> (header, [dict]) ; 'id' reste une chaine, le reste est converti en float si possible."""
    with open(path, newline="", encoding="ascii") as fh:
        rows = list(csv.reader(fh))
    header = [h.strip() for h in rows[0]]
    out = []
    for r in rows[1:]:
        if not r or all(not c.strip() for c in r):
            continue
        d = {}
        for h, v in zip(header, r):
            if h == "id":
                d[h] = v.strip()
            else:
                try:
                    d[h] = float(v)
                except ValueError:
                    d[h] = v
        out.append(d)
    return header, out


def _col(rows, key):
    return np.array([float(r.get(key, np.nan)) if r.get(key) not in (None, "") else np.nan for r in rows], dtype=float)


# ---------------------------------------------------------------------------
# Lecture des runs
# ---------------------------------------------------------------------------
def run_done(out):
    """Run termine : _run.json avec rc == 0 et history.csv (ou _metrics.json deja calcule)."""
    st = os.path.join(out, "_run.json")
    if not os.path.exists(st):
        return False
    try:
        rc = json.load(open(st)).get("rc")
    except Exception:
        return False
    return rc == 0 and (os.path.exists(os.path.join(out, "history.csv")) or os.path.exists(os.path.join(out, "_metrics.json")))


def run_metrics(out, conf, targets, cache=True):
    """Observables d un run (extract.py), lues dans out/_metrics.json si present, sinon calculees
    (et ecrites si cache). None si le run n est pas termine."""
    if not run_done(out):
        return None
    mp = os.path.join(out, "_metrics.json")
    if os.path.exists(mp):
        try:
            return json.load(open(mp))
        except Exception:
            pass
    import extract
    try:
        if conf_prefix(conf) == "BTS":
            m = extract.bts_observables(out)
        else:
            m = extract.observables(extract.load_run(out), targets)
    except Exception as e:  # history.csv corrompu, cfg absent...
        print(f"  [collect] {out}: extraction impossible ({e})")
        return None
    if cache:
        try:
            json.dump(m, open(mp, "w"), indent=1)
        except OSError:
            pass
    return m


def merge_seeds(lst):
    """Plusieurs graines pour un meme (point, confinement) : moyenne des valeurs finies."""
    if len(lst) == 1:
        return dict(lst[0])
    out = dict(lst[0])
    keys = set().union(*[set(m.keys()) for m in lst])
    for k in keys:
        vals = []
        for m in lst:
            v = m.get(k)
            try:
                v = float(v)
            except (TypeError, ValueError):
                continue
            if math.isfinite(v):
                vals.append(v)
        if vals:
            out[k] = float(np.mean(vals))
        elif k not in out:
            out[k] = float("nan")
    out["n_seeds"] = len(lst)
    return out


# ---------------------------------------------------------------------------
# init
# ---------------------------------------------------------------------------
def cmd_init(a):
    runs = resolve_runs(a.tag, a.runs)
    P = Paths(a.tag, runs)
    cfg = load_config(P)
    cfg.update(dict(tag=a.tag, runs=os.path.relpath(runs, ROOT) if runs.startswith(ROOT) else runs,
                    template=a.template, bts=a.bts, extra=list(a.extra or []), conf=[float(c) for c in a.conf],
                    seeds=list(a.seeds or []), nu=float(a.nu), created=time.strftime("%Y-%m-%d %H:%M:%S")))
    if a.weights:
        cfg["weights"] = a.weights
    if a.targets:
        cfg["targets"] = a.targets
    if a.seuils:
        cfg["seuils"] = a.seuils
    save_config(P, cfg)
    load_space(P)
    print(f"[init] {P.cfg}")
    print(json.dumps(cfg, indent=1))
    return P


# ---------------------------------------------------------------------------
# collect
# ---------------------------------------------------------------------------
def cmd_collect(a):
    P = Paths(a.tag, resolve_runs(a.tag, a.runs))
    cfg = load_config(P)
    if a.weights:
        cfg["weights"] = os.path.abspath(a.weights)
        save_config(P, cfg)
    space, _ = load_space(P)
    S = objective_settings(cfg)
    prefixes = [conf_prefix(c) for c in S["confinements"]]

    designs = P.designs()
    if not designs:
        raise SystemExit(f"aucun design_{a.tag}.csv dans {P.runs}")
    points = []
    for k, path in designs:
        X, ids = space.read_csv(path)
        for x, i in zip(X, ids):
            points.append(dict(id=i, round=k, x=x))
    jobs = []
    for jf in P.jobs_files():
        jobs += json.load(open(jf))
    by_id = {}
    for j in jobs:
        by_id.setdefault(str(j["id"]), []).append(j)

    import extract
    tp = resolve_path(cfg.get("targets")) if cfg.get("targets") else extract.TARGETS
    targets = json.load(open(tp)) if os.path.exists(tp) else None

    rows = []
    n_done_tot = 0
    for p in points:
        js = by_id.get(p["id"], [])
        per_conf = {}
        n_done = 0
        for j in js:
            out = resolve_path(j["out"])
            m = run_metrics(out, j["conf"], targets, cache=not a.no_cache)
            if m is None:
                continue
            n_done += 1
            per_conf.setdefault(conf_prefix(j["conf"]), []).append(m)
        merged = {pref: merge_seeds(lst) for pref, lst in per_conf.items()}
        complete = len(js) > 0 and n_done == len(js)
        n_done_tot += n_done
        row = dict(id=p["id"], round=p["round"])
        for name, v in zip(space.names, p["x"]):
            row[name] = float(v)
        row.update(complete=int(complete), n_jobs=len(js), n_done=n_done)
        for pref in prefixes:
            o = merged.get(pref, {})
            for c in TRIAX_COLS:
                row[f"{pref}_{c}"] = o.get(c, float("nan"))
        ob = merged.get("BTS", {})
        for c in BTS_COLS:
            row[f"BTS_{c}"] = ob.get(c, float("nan"))
        # objectif (points complets seulement : un run manquant fausserait J par la penalite 'absent')
        row["J"] = float("nan")
        for pref in prefixes + ["BTS"]:
            row[f"J_{pref}"] = float("nan")
            terms = TERMS if pref != "BTS" else ("d_bts",)
            for t in terms:
                row[f"J_{pref}_{t}"] = float("nan")
            if pref != "BTS":
                row[f"unbroken_{pref}"] = float("nan")
        if complete:
            obs_by_conf = {prefix_conf(pref): merged[pref] for pref in merged if pref != "BTS"}
            J, det = eval_objective(obs_by_conf, merged.get("BTS"), S)
            row["J"] = J
            for key, d in det.items():
                if key == "total":
                    continue
                pref = "BTS" if key == "BTS" else conf_prefix(key)
                row[f"J_{pref}"] = d["sum"]
                for t, v in d["terms"].items():
                    row[f"J_{pref}_{t}"] = v
                if "unbroken" in d:
                    row[f"unbroken_{pref}"] = d["unbroken"]
        rows.append(row)

    header = list(rows[0].keys())
    write_table(P.base, header, rows)
    n_comp = sum(int(r["complete"]) for r in rows)
    print(f"[collect] {len(points)} points ({len(designs)} designs), {len(jobs)} jobs, {n_done_tot} runs termines, "
          f"{n_comp} points complets -> {P.base}")
    comp = sorted([r for r in rows if r["complete"] and math.isfinite(r["J"])], key=lambda r: r["J"])
    if comp:
        print(f"  J : min {comp[0]['J']:.3f} ({comp[0]['id']}), mediane {np.median([r['J'] for r in comp]):.3f}")
        for r in comp[:5]:
            print("   " + f"{r['id']:14s} J {r['J']:8.3f}  " + " ".join(f"{k}={r[k]:.4g}" for k in space.names))
    return rows


# ---------------------------------------------------------------------------
# fit
# ---------------------------------------------------------------------------
def build_ydict(rows, confinements, with_bts=True, nb_min=1.0):
    """Lignes de base -> ({cle emulee: y transforme (n,)}, {prefixe: rompu (n,) bool}).
    Les runs non rompus (nb_end < nb_min) sont mis a NaN pour toutes les observables du
    confinement : la regression n apprend que la branche 'rompu' (A4 4.2)."""
    y, broken = {}, {}
    for conf in confinements:
        pref = conf_prefix(conf)
        nb = _col(rows, f"{pref}_nb_end")
        brk = np.isfinite(nb) & (nb >= nb_min)
        broken[pref] = brk
        for key, name, tr in OBS_TRIAX:
            fwd = TRANSFORMS[tr][0]
            v = fwd(_col(rows, f"{pref}_{key}"))
            v[~brk] = np.nan
            y[f"{pref}_{name}"] = v
    if with_bts:
        for key, name, tr in OBS_BTS:
            v = TRANSFORMS[tr][0](_col(rows, f"BTS_{key}"))
            if np.isfinite(v).any():
                y[f"BTS_{name}"] = v
    return y, broken


def fit_classifiers(U, broken, names, min_locked=2, seed=0):
    """Classifieur GP rompu / bloque par confinement (seulement si >= min_locked runs bloques
    ET >= min_locked rompus) -> {prefixe: GaussianProcessClassifier}."""
    from sklearn.gaussian_process import GaussianProcessClassifier
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern
    clf = {}
    for pref, brk in broken.items():
        m = np.isfinite(U).all(axis=1)
        n_lock = int((~brk[m]).sum()); n_brk = int(brk[m].sum())
        if n_lock < min_locked or n_brk < min_locked:
            print(f"  [fit] {pref} : {n_brk} rompus, {n_lock} bloques -> pas de classifieur"
                  + (" (tous rompus)" if n_lock == 0 else ""))
            continue
        kern = ConstantKernel(1.0, (1e-2, 1e2)) * Matern(length_scale=np.full(U.shape[1], 0.5), length_scale_bounds=(1e-2, 1e2), nu=2.5)
        g = GaussianProcessClassifier(kernel=kern, n_restarts_optimizer=2, random_state=seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            g.fit(U[m], brk[m].astype(int))
        p = g.predict_proba(U[m])[:, 1]
        acc = float(np.mean((p >= 0.5) == brk[m]))
        print(f"  [fit] {pref} : {n_brk} rompus, {n_lock} bloques -> classifieur GP (precision apprentissage {acc:.2f}, "
              f"P(rompu) moyen sur les bloques {p[~brk[m]].mean():.2f})")
        clf[pref] = g
    return clf


def et_loo(X, y, n_trees=100, seed=0):
    """LOO d un ExtraTrees (controle sans hypothese de lissite) -> (rmse, q2)."""
    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.model_selection import KFold, LeaveOneOut
    m = np.isfinite(y) & np.isfinite(X).all(axis=1)
    X, y = X[m], y[m]
    n = len(y)
    if n < 5:
        return float("nan"), float("nan")
    cv = LeaveOneOut() if n <= 60 else KFold(10, shuffle=True, random_state=seed)
    yhat = np.full(n, np.nan)
    for tr, te in cv.split(X):
        et = ExtraTreesRegressor(n_estimators=n_trees, min_samples_leaf=2, random_state=seed, n_jobs=1)
        yhat[te] = et.fit(X[tr], y[tr]).predict(X[te])
    rmse = float(np.sqrt(np.mean((yhat - y) ** 2)))
    ss = float(np.sum((y - y.mean()) ** 2))
    return rmse, (1.0 - float(np.sum((yhat - y) ** 2)) / ss if ss > 0 else float("nan"))


def cmd_fit(a):
    import joblib
    P = Paths(a.tag, resolve_runs(a.tag, a.runs))
    cfg = load_config(P)
    space, spec = load_space(P)
    S = objective_settings(cfg)
    if not os.path.exists(P.base):
        raise SystemExit(f"{P.base} absent : lancer collect d abord")
    header, rows = read_table(P.base)
    comp = [r for r in rows if int(r.get("complete", 0)) == 1]
    if len(comp) < 5:
        raise SystemExit(f"{len(comp)} points complets : trop peu pour ajuster")
    Xp = np.array([[r[k] for k in space.names] for r in comp], dtype=float)
    U = space.to_unit(Xp)
    bad = [r["id"] for r, u in zip(comp, U) if (u < -1e-6).any() or (u > 1 + 1e-6).any()]
    if bad:
        raise SystemExit(f"points hors de l espace {P.space} : {bad[:5]}... (design d un autre espace ?)")
    nb_min = float(cfg.get("nb_min", 1.0))
    ydict, broken = build_ydict(comp, S["confinements"], nb_min=nb_min)
    # observables inexploitables (< 3 valeurs finies) ecartees avec avertissement
    for k in list(ydict):
        if np.isfinite(ydict[k]).sum() < 3:
            print(f"  [fit] {k} : {np.isfinite(ydict[k]).sum()} valeurs finies, observable ecartee")
            del ydict[k]
    n_restarts = 2 if a.fast else 8
    t0 = time.time()
    emu = MultiEmulator("gp", names=space.names, seed=0, n_restarts=n_restarts).fit(U, ydict)
    print(f"[fit] {len(comp)} points complets, {len(ydict)} observables, GP Matern 5/2 ARD + nugget "
          f"({n_restarts} redemarrages) en {time.time() - t0:.1f} s ; regression sur les runs rompus (nb_end >= {nb_min:g})")
    clf = fit_classifiers(U, broken, space.names)

    # controle de coherence cibles / nu : d_peak reconstruit vs d_peak stocke
    try:
        T = exp_targets(cfg.get("nu", NU_DEFAULT), resolve_path(cfg.get("targets")), resolve_path(cfg.get("seuils")))
        for conf in S["confinements"]:
            pref = conf_prefix(conf)
            q = _col(comp, f"{pref}_q_peak"); dp = _col(comp, f"{pref}_d_peak")
            m = np.isfinite(q) & np.isfinite(dp)
            if m.any():
                err = np.max(np.abs(q[m] / T[str(conf)]["q_peak"] - 1 - dp[m]))
                if err > 1e-3:
                    print(f"  [fit] ATTENTION {pref} : d_peak reconstruit s ecarte de {err:.2e} du d_peak stocke (cibles differentes ?)")
            e = _col(comp, f"{pref}_eps_peak"); de = _col(comp, f"{pref}_d_eps_peak")
            m = np.isfinite(e) & np.isfinite(de)
            if m.any():
                err = np.max(np.abs((e[m] / 100.0) / T[str(conf)]["eps_peak"] - 1 - de[m]))
                if err > 1e-3:
                    print(f"  [fit] ATTENTION {pref} : d_eps_peak reconstruit s ecarte de {err:.2e} (nu de la campagne = {cfg.get('nu', NU_DEFAULT)} ?)")
    except FileNotFoundError as e:
        print(f"  [fit] cibles non verifiees ({e})")

    # LOO, Q2, nugget, ExtraTrees en controle
    loo, q2, et = {}, {}, {}
    t0 = time.time()
    print(f"\n{'observable':16s} {'n':>4s} {'ecart':>5s} {'sd(y)':>8s} {'nugget':>8s} {'LOO rmse':>9s} {'Q2':>7s}"
          + (f" {'ET rmse':>8s} {'ET Q2':>7s}" if a.rf else ""))
    for k, m in emu.models.items():
        rmse, yhat = m.loo_rmse(return_pred=True, refit_hyper=not a.fast)
        ss = float(np.sum((m.y_ - m.y_.mean()) ** 2))
        loo[k] = rmse
        q2[k] = float(1.0 - np.sum((yhat - m.y_) ** 2) / ss) if ss > 0 else float("nan")
        line = f"{k:16s} {len(m.y_):4d} {m.n_dropped:5d} {np.std(m.y_):8.4f} {m.nugget_std:8.4f} {rmse:9.4f} {q2[k]:7.3f}"
        if a.rf:
            et[k] = et_loo(U, ydict[k], n_trees=a.rf_trees)
            line += f" {et[k][0]:8.4f} {et[k][1]:7.3f}"
        print(line)
    print(f"(validation croisee en {time.time() - t0:.1f} s ; LOO {'a hyperparametres figes' if a.fast else 'avec re-ajustement'})")
    if a.rf:
        worse = [k for k in et if np.isfinite(et[k][0]) and et[k][0] < 0.8 * loo[k]]
        if worse:
            print(f"  ExtraTrees nettement meilleur que le GP sur {worse} : discontinuite probable (seuil de rupture), voir A4 4.2")

    imp = emu.importance()
    print(f"\nimportances ARD (1/longueur, somme 1)\n{'observable':16s} " + " ".join(f"{n:>8s}" for n in space.names))
    for k, d in imp.items():
        print(f"{k:16s} " + " ".join(f"{d[n]:8.3f}" for n in space.names))
    mean_imp = {n: float(np.mean([d[n] for d in imp.values()])) for n in space.names}
    print(f"{'moyenne':16s} " + " ".join(f"{mean_imp[n]:8.3f}" for n in space.names))

    bundle = dict(tag=a.tag, names=space.names, space=spec, emu=emu, clf=clf, X=U, ids=[r["id"] for r in comp],
                  broken={k: v.tolist() for k, v in broken.items()}, nb_min=nb_min,
                  obs_keys=list(ydict.keys()), loo=loo, q2=q2, et=et, importance=imp, objective=S,
                  nu=float(cfg.get("nu", NU_DEFAULT)), targets=cfg.get("targets"), seuils=cfg.get("seuils"),
                  fitted=time.strftime("%Y-%m-%d %H:%M:%S"), fast=bool(a.fast))
    joblib.dump(bundle, P.models)
    print(f"\n[fit] modeles -> {P.models}")
    return bundle


# ---------------------------------------------------------------------------
# Substitut : predictions -> observables -> objectif
# ---------------------------------------------------------------------------
def j_numpy(arrs, S):
    """Version vectorisee de objective_with_bts sur des tableaux (N,) par observable.
    arrs = {conf: {d_peak, d_eps_peak, d_CI, d_CD, d_drop, rmse_curve, nb_end}, 'BTS': {d_bts}}."""
    w = dict(W_DEFAULT); w.update(S["weights"])
    t = dict(TOL_DEFAULT); t.update(S["tol"])
    total = None
    for conf in S["confinements"]:
        o = arrs[conf]
        for k in TERMS:
            v = np.asarray(o[k], dtype=float)
            fin = np.isfinite(v)
            r = np.where(fin, v / t[k], 0.0)
            term = np.where(fin, w[k] * r * r, w[k] * S["nan_penalty"])
            total = term if total is None else total + term
        nb = np.asarray(o["nb_end"], dtype=float)
        total = total + np.where(np.isfinite(nb) & (nb > 0), 0.0, S["unbroken_penalty"])
    if "BTS" in arrs:
        d = np.asarray(arrs["BTS"]["d_bts"], dtype=float)
        fin = np.isfinite(d)
        total = total + np.where(fin, S["w_bts"] * (np.where(fin, d, 0.0) / S["bts_tol"]) ** 2, S["w_bts"] * 9.0)
    return total


class Surrogate:
    """Emulateurs + cibles -> predict_obs(u), J_hat(u), LCB(u)."""

    def __init__(self, bundle, S=None):
        self.emu = bundle["emu"]
        self.space = ParamSpace({k: tuple(v) for k, v in bundle["space"].items()})
        self.S = S or bundle["objective"]
        self.confs = [int(c) for c in self.S["confinements"]]
        self.T = exp_targets(bundle.get("nu", NU_DEFAULT), resolve_path(bundle.get("targets")), resolve_path(bundle.get("seuils")))
        self.keys = list(self.emu.keys())
        self.has_bts = "BTS_lbts" in self.keys
        self.clf = bundle.get("clf", {}) or {}

    # -- reconstruction des observables relatives ------------------------------
    def obs_from_raw(self, vals, pbrk=None):
        """vals = {cle emulee: tableau (N,) transforme}, pbrk = {prefixe: P(rompu) (N,)} optionnel
        -> {conf: {obs: (N,)}, 'BTS': {...}} ; nb_end = 0 la ou P(rompu) < 0,5."""
        n = len(next(iter(vals.values())))
        nan = np.full(n, np.nan)

        def get(pref, name, tr):
            k = f"{pref}_{name}"
            return TRANSFORMS[tr][1](vals[k]) if k in vals else nan

        out = {}
        for conf in self.confs:
            pref = conf_prefix(conf)
            tg = self.T.get(str(conf), {})
            q = get(pref, "lq_peak", "log")
            eps = get(pref, "eps_peak", "id")
            ci = get(pref, "CI_frac", "id")
            cd = get(pref, "CD_frac", "id")
            drop = get(pref, "drop", "id")
            rmse = get(pref, "lrmse", "log1p")
            nb = np.rint(get(pref, "lnb_end", "log1p")) if f"{pref}_lnb_end" in vals else np.ones(n)
            if pbrk is not None and pref in pbrk:
                nb = np.where(pbrk[pref] >= 0.5, nb, 0.0)
            o = dict(q_peak=q, eps_peak=eps, CI_frac=ci, CD_frac=cd, drop=drop, rmse_curve=rmse, nb_end=nb,
                     q_CI=ci * q, q_CD=cd * q)
            o["d_peak"] = q / tg["q_peak"] - 1 if "q_peak" in tg else nan
            o["d_eps_peak"] = (eps / 100.0) / tg["eps_peak"] - 1 if "eps_peak" in tg else nan
            o["d_CI"] = ci * q / tg["q_CI"] - 1 if "q_CI" in tg else nan
            o["d_CD"] = cd * q / tg["q_CD"] - 1 if "q_CD" in tg else nan
            o["d_drop"] = drop - tg["chute"] if "chute" in tg else nan
            out[conf] = o
        if self.has_bts:
            b = np.exp(vals["BTS_lbts"])
            out["BTS"] = dict(bts=b, d_bts=b / self.T["BTS"] - 1)
        return out

    def _means(self, U):
        """-> (moyennes, ecarts-types latents, P(rompu)) par cle / prefixe sur (N, d)."""
        U = np.atleast_2d(U)
        pr = self.emu.predict(U)
        pbrk = {pref: g.predict_proba(U)[:, 1] for pref, g in self.clf.items()}
        return {k: v[0] for k, v in pr.items()}, {k: v[1] for k, v in pr.items()}, pbrk

    def predict_obs(self, u):
        """u (d,) -> {20: obs, 50: obs, 'BTS': {bts, d_bts}} (floats) pour objective_with_bts."""
        means, _, pbrk = self._means(np.atleast_2d(np.asarray(u, dtype=float)))
        arrs = self.obs_from_raw(means, pbrk)
        return {c: {k: float(v[0]) for k, v in o.items()} for c, o in arrs.items()}

    def predict_obs_batch(self, U):
        means, _, pbrk = self._means(U)
        return self.obs_from_raw(means, pbrk)

    def p_broken(self, U):
        """{prefixe: P(rompu) (N,)} (vide sans classifieur)."""
        return self._means(U)[2]

    def j_batch(self, U):
        """(N, d) -> J_hat (N,) sur les predictions moyennes."""
        return j_numpy(self.predict_obs_batch(U), self.S)

    def jhat(self, u):
        return float(self.j_batch(np.atleast_2d(np.asarray(u, dtype=float)))[0])

    def j_objective(self, u):
        """J_hat par objective_with_bts (reference, un point) -> (J, detail)."""
        obs = self.predict_obs(u)
        return eval_objective({c: o for c, o in obs.items() if c != "BTS"}, obs.get("BTS"), self.S)

    def draw_matrix(self, n_draws=200, seed=1):
        """Tirages normaux communs (n_obs, n_draws) pour la propagation Monte-Carlo."""
        return np.random.default_rng(seed).standard_normal((len(self.keys), n_draws))

    def lcb_batch(self, U, kappa, Z):
        """-> (LCB, mean_J, std_J) sur (N, d) par propagation de Z (n_obs, ns) a travers les GP."""
        U = np.atleast_2d(U)
        means, stds, pbrk = self._means(U)
        N, ns = len(U), Z.shape[1]
        vals = {}
        for i, k in enumerate(self.keys):
            vals[k] = (means[k][:, None] + stds[k][:, None] * Z[i][None, :]).ravel()
        pb = {pref: np.repeat(p, ns) for pref, p in pbrk.items()}
        J = j_numpy(self.obs_from_raw(vals, pb), self.S).reshape(N, ns)
        mj, sj = J.mean(axis=1), J.std(axis=1)
        return mj - kappa * sj, mj, sj

    def believer_update(self, u):
        """Kriging believer : la prediction moyenne en u devient une pseudo-observation de
        chaque GP, re-ajuste a hyperparametres FIGES (Ginsbourger 2010)."""
        U = np.atleast_2d(np.asarray(u, dtype=float))
        means, _, _ = self._means(U)
        for k, m in self.emu.models.items():
            Xa = np.vstack([m.X_, U])
            ya = np.append(m.y_, means[k][0])
            gp = m._make_gp(Xa.shape[1], m._alpha(ya, None), kernel=m.model.kernel_, n_restarts=0)
            gp.set_params(optimizer=None)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gp.fit(Xa, ya)
            m.model, m.X_, m.y_ = gp, Xa, ya

    def check_numpy_objective(self, n=5, seed=0):
        """Egalite j_numpy / objective_with_bts sur n points aleatoires -> ecart relatif max."""
        U = np.random.default_rng(seed).random((n, self.space.d))
        jb = self.j_batch(U)
        jo = np.array([self.j_objective(u)[0] for u in U])
        return float(np.max(np.abs(jb - jo) / np.maximum(1.0, np.abs(jo))))


def load_surrogate(tag, runs=None, weights=None):
    import joblib
    P = Paths(tag, resolve_runs(tag, runs))
    if not os.path.exists(P.models):
        raise SystemExit(f"{P.models} absent : lancer fit d abord")
    bundle = joblib.load(P.models)
    cfg = load_config(P)
    S = objective_settings(cfg, weights)
    return P, cfg, bundle, Surrogate(bundle, S)


# ---------------------------------------------------------------------------
# optimize
# ---------------------------------------------------------------------------
def distinct_best(results, dmin=0.1, n=5):
    """results = [(J, u)] -> les n meilleurs a distance mutuelle > dmin."""
    out = []
    for J, u in sorted(results, key=lambda r: r[0]):
        if all(np.linalg.norm(u - v) > dmin for _, v in out):
            out.append((J, u))
        if len(out) >= n:
            break
    return out


def cmd_optimize(a):
    P, cfg, bundle, surr = load_surrogate(a.tag, a.runs, a.weights)
    err = surr.check_numpy_objective()
    if err > 1e-9:
        raise SystemExit(f"j_numpy et objective_with_bts divergent (ecart relatif {err:.2e})")
    print(f"[optimize] {len(surr.keys)} emulateurs, {len(surr.clf)} classifieurs rompu/bloque, objectif numpy verifie (ecart {err:.1e}) ; "
          f"APSO x {a.n_seeds} graines, {a.particles} particules, {a.iters} generations")
    B = surr.space.bounds_unit()
    results = []
    t0 = time.time()
    for s in range(a.n_seeds):
        xb, fb, hist = apso(surr.j_batch, B, n_particles=a.particles, iters=a.iters, seed=a.seed0 + s, vectorized=True)
        results.append((float(fb), np.asarray(xb)))
        print(f"  graine {a.seed0 + s:4d}  J* {fb:10.4f}  u* " + " ".join(f"{v:.3f}" for v in xb)
              + f"  ({hist['n_eval'][-1]} evaluations)")
    print(f"  ({time.time() - t0:.1f} s)")
    Jb, ub = min(results, key=lambda r: r[0])
    xb = surr.space.from_unit(ub)
    J, det = surr.j_objective(ub)
    Z = surr.draw_matrix(200, seed=1)
    _, mj, sj = surr.lcb_batch(ub[None, :], 0.0, Z)
    print(f"\n[optimize] MAP du substitut : J* = {Jb:.4f} (Monte-Carlo : {mj[0]:.3f} +- {sj[0]:.3f})")
    print("  u* : " + " ".join(f"{n}={v:.3f}" for n, v in zip(surr.space.names, ub)))
    print("  x* : " + " ".join(f"{n}={v:.5g}" for n, v in zip(surr.space.names, xb)))
    print(format_detail_bts(det))
    obs = surr.predict_obs(ub)
    for c in surr.confs:
        o = obs[c]
        print(f"  P{c:03d} predit : q_peak {o['q_peak']:.1f} MPa, eps_peak {o['eps_peak']:.3f} %, CI {o['CI_frac']:.3f}, "
              f"CD {o['CD_frac']:.3f}, chute {o['drop']:.3f}, rmse {o['rmse_curve']:.1f} MPa, rompus {o['nb_end']:.0f}")
    if "BTS" in obs:
        print(f"  BTS predit : {obs['BTS']['bts']:.2f} MPa ({100 * obs['BTS']['d_bts']:+.1f} %)")
    top = distinct_best(results, dmin=0.1, n=5)
    print(f"\n  {len(top)} meilleurs distincts (distance > 0.1 en unites) :")
    for J_i, u_i in top:
        x_i = surr.space.from_unit(u_i)
        print(f"    J {J_i:9.4f}  " + " ".join(f"{n}={v:.4g}" for n, v in zip(surr.space.names, x_i))
              + f"   |u - u*| = {np.linalg.norm(u_i - ub):.3f}")
    spread = float(np.max([np.linalg.norm(u - ub) for _, u in results]))
    n_close = int(sum(np.linalg.norm(u - ub) < 0.05 for _, u in results))
    print(f"  {n_close}/{len(results)} graines a moins de 0.05 de u* ; dispersion max {spread:.3f}")

    def _json_obs(o):
        return {c: {k: float(v) for k, v in d.items()} for c, d in o.items()}

    out = dict(tag=a.tag, J=Jb, J_mc_mean=float(mj[0]), J_mc_std=float(sj[0]), u=ub.tolist(),
               x=surr.space.as_dict(xb), names=surr.space.names, detail=det, obs_pred=_json_obs(obs),
               top=[dict(J=J_i, u=u_i.tolist(), x=surr.space.as_dict(surr.space.from_unit(u_i))) for J_i, u_i in top],
               runs=[dict(J=J_i, u=u_i.tolist()) for J_i, u_i in results], objective=surr.S,
               apso=dict(n_seeds=a.n_seeds, particles=a.particles, iters=a.iters, seed0=a.seed0),
               models=os.path.basename(P.models), fitted=bundle.get("fitted"), date=time.strftime("%Y-%m-%d %H:%M:%S"))
    json.dump(out, open(P.map, "w"), indent=1)
    print(f"[optimize] -> {P.map}")
    return out


# ---------------------------------------------------------------------------
# enrich
# ---------------------------------------------------------------------------
def _min_dist(u, pts):
    if len(pts) == 0:
        return float("inf")
    return float(np.min(np.linalg.norm(np.asarray(pts) - u, axis=1)))


def cmd_enrich(a):
    P, cfg, bundle, surr = load_surrogate(a.tag, a.runs, a.weights)
    space = surr.space
    # points existants = tous les points du design (complets ou en attente), en unites
    existing = []
    if os.path.exists(P.base):
        _, rows = read_table(P.base)
        existing = [space.to_unit(np.array([r[k] for k in space.names], dtype=float)) for r in rows]
    else:
        existing = [u for u in bundle["X"]]
    existing = [np.clip(u, 0.0, 1.0) for u in existing]
    Z = surr.draw_matrix(a.n_draws, seed=a.seed0)
    believer = copy.deepcopy(surr)
    B = space.bounds_unit()
    batch, info = [], []
    print(f"[enrich] {len(existing)} points existants, lot de {a.n} par LCB (kappa {a.kappa}, {a.n_draws} tirages), "
          f"Kriging believer, distance min {a.dmin}")
    t0 = time.time()
    for i in range(a.n):
        pts = existing + batch

        def acq(U, pts=pts):
            lcb = believer.lcb_batch(U, a.kappa, Z)[0]
            # penalite continue dans la boule d exclusion (dmin) autour de chaque point
            D = np.sqrt(((U[:, None, :] - np.asarray(pts)[None, :, :]) ** 2).sum(-1)).min(axis=1)
            return lcb + 1e3 * np.clip(1.0 - D / a.dmin, 0.0, None)

        chosen = None
        for attempt in range(3):
            xb, fb, _ = apso(acq, B, n_particles=a.particles, iters=a.iters, seed=a.seed0 + 10 * i + attempt, vectorized=True)
            if _min_dist(xb, pts) >= a.dmin:
                chosen = xb
                break
        if chosen is None:
            print(f"  point {i}: aucun candidat a plus de {a.dmin} des points existants apres 3 essais, lot arrete a {len(batch)}")
            break
        lcb, mj, sj = surr.lcb_batch(chosen[None, :], a.kappa, Z)          # moments du VRAI substitut
        lcb_b, mj_b, sj_b = believer.lcb_batch(chosen[None, :], a.kappa, Z)
        batch.append(chosen)
        info.append(dict(lcb=float(lcb[0]), mean=float(mj[0]), std=float(sj[0]), dist=_min_dist(chosen, pts)))
        print(f"  {i + 1:2d}. u = " + " ".join(f"{v:.3f}" for v in chosen)
              + f"  LCB {lcb[0]:8.3f} (mean {mj[0]:7.3f} sd {sj[0]:6.3f} ; believer sd {sj_b[0]:6.3f})  dist {info[-1]['dist']:.3f}")
        believer.believer_update(chosen)
    print(f"  ({time.time() - t0:.1f} s)")
    if not batch:
        raise SystemExit("lot vide")
    k = P.next_round()
    ids = [f"{a.tag}_r{k}_{i:03d}" for i in range(len(batch))]
    Xn = space.from_unit(np.array(batch))
    path = P.design_round(k)
    space.write_csv(path, Xn, ids)
    dmin_all = min(_min_dist(u, existing) for u in batch)
    print(f"[enrich] lot r{k} : {len(batch)} points -> {path} (distance min aux existants {dmin_all:.3f}, "
          f"entre nouveaux {min_pairwise_distance(np.array(batch)):.3f})")
    for i, x in zip(ids, Xn):
        print("   " + f"{i:14s} " + " ".join(f"{n}={v:.4g}" for n, v in zip(space.names, x)))
    json.dump(dict(round=k, ids=ids, u=[u.tolist() for u in batch], info=info, kappa=a.kappa, n_draws=a.n_draws),
              open(os.path.join(P.runs, f"enrich_{a.tag}_r{k}.json"), "w"), indent=1)

    # decks par make_decks.py (memes template / bts / extra que le lot initial)
    tpl = resolve_path(cfg.get("template")) if cfg.get("template") else None
    if a.no_decks or tpl is None or not os.path.exists(tpl):
        print("[enrich] template absent ou --no-decks : CSV seul, pas de decks")
        return path
    cmd = [sys.executable, os.path.join(HERE, "make_decks.py"), path, "--template", tpl, "--tag", f"{a.tag}_r{k}",
           "--outdir", P.runs, "--conf"] + [f"{float(c):g}" for c in cfg.get("conf", CONF_DEFAULT)]
    if cfg.get("bts"):
        cmd += ["--bts", resolve_path(cfg["bts"])]
    if cfg.get("seeds"):
        cmd += ["--seeds"] + [str(s) for s in cfg["seeds"]]
    if cfg.get("extra"):
        cmd += ["--extra"] + list(cfg["extra"])
    print("[enrich] " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        print(r.stderr)
        raise SystemExit(f"make_decks.py a echoue (rc {r.returncode})")
    print(f"[enrich] decks du lot r{k} prets : jobs_{a.tag}_r{k}.json dans {P.runs} (RIEN n est lance)")
    return path


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def _md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def cmd_report(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import joblib

    P = Paths(a.tag, resolve_runs(a.tag, a.runs))
    cfg = load_config(P)
    space, _ = load_space(P)
    S = objective_settings(cfg)
    if not os.path.exists(P.base):
        raise SystemExit(f"{P.base} absent : lancer collect d abord")
    header, rows = read_table(P.base)
    comp = sorted([r for r in rows if int(r.get("complete", 0)) == 1 and math.isfinite(r["J"])], key=lambda r: r["J"])
    if not comp:
        raise SystemExit("aucun point complet avec J fini")
    bundle = joblib.load(P.models) if os.path.exists(P.models) else None
    mapj = json.load(open(P.map)) if os.path.exists(P.map) else None
    try:
        T = exp_targets(cfg.get("nu", NU_DEFAULT), resolve_path(cfg.get("targets")), resolve_path(cfg.get("seuils")))
    except FileNotFoundError:
        T = None
    jobs = []
    for jf in P.jobs_files():
        jobs += json.load(open(jf))
    out_of = {}
    for j in jobs:
        out_of.setdefault((str(j["id"]), conf_prefix(j["conf"])), resolve_path(j["out"]))

    J = np.array([r["J"] for r in comp])
    n = len(comp)
    n_top = max(5, int(math.ceil(0.10 * n)))
    top = comp[:min(n_top, n)]
    Utop = space.to_unit(np.array([[r[k] for k in space.names] for r in top], dtype=float))
    corr = np.corrcoef(Utop.T) if len(top) >= 3 else np.full((space.d, space.d), np.nan)
    with np.errstate(all="ignore"):
        cov = np.cov(Utop.T) if len(top) >= 3 else None
        eig = np.sort(np.linalg.eigvalsh(cov))[::-1] if cov is not None else None

    plt.rcParams.update({"font.family": "serif", "font.size": 9})
    fig = plt.figure(figsize=(15, 9.5))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    # (a) J par rang
    ax = fig.add_subplot(gs[0, 0])
    ax.semilogy(np.arange(1, n + 1), J, "o", ms=3, color="0.3")
    for i in range(min(3, n)):
        ax.annotate(comp[i]["id"], (i + 1, J[i]), fontsize=7, xytext=(6, 2 + 9 * i), textcoords="offset points")
    if mapj is not None:
        ax.axhline(mapj["J"], color="C3", ls="--", lw=1, label=f"J_hat(MAP) = {mapj['J']:.2f}")
        ax.legend(fontsize=8)
    ax.set_xlabel("rang"); ax.set_ylabel("J"); ax.set_title(f"(a) J des {n} points complets")
    ax.grid(alpha=0.3, which="both")

    # (b) courbes des 3 meilleurs vs bande experimentale
    best3 = comp[:3]
    curve_note = []
    for ic, conf in enumerate(S["confinements"][:2]):
        ax = fig.add_subplot(gs[0, 1 + ic])
        key = str(int(conf))
        if T is not None and key in T:
            tg = T[key]
            e = 100 * tg["eps_grid"]
            ax.fill_between(e, tg["q_mean"] - tg["q_std"], tg["q_mean"] + tg["q_std"], color="0.8", label="exp +- 1 sd")
            ax.plot(e, tg["q_mean"], "k-", lw=1, label="exp moyenne")
        import extract
        for i, r in enumerate(best3):
            od = out_of.get((r["id"], conf_prefix(conf)))
            if od and os.path.exists(os.path.join(od, "history.csv")):
                try:
                    run = extract.load_run(od)
                    ax.plot(100 * run["eps"], run["q"], color=f"C{i}", lw=1.2, label=f"{r['id']} (J {r['J']:.2f})")
                except Exception as e:
                    curve_note.append(f"{r['id']} P{int(conf):03d}: {e}")
            else:
                curve_note.append(f"{r['id']} P{int(conf):03d}: pas de history.csv")
        ax.set_xlabel("eps axial [%]"); ax.set_ylabel("q [MPa]")
        ax.set_title(f"(b) q(eps) a sigma3 = {int(conf)} MPa, 3 meilleurs")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)
        if curve_note and ic == 0:
            ax.text(0.02, 0.02, "\n".join(curve_note[:6]), transform=ax.transAxes, fontsize=6, va="bottom")

    # (c) correlations des 10 % meilleurs
    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(space.d)); ax.set_xticklabels(space.names, rotation=45)
    ax.set_yticks(range(space.d)); ax.set_yticklabels(space.names)
    for i in range(space.d):
        for j in range(space.d):
            if np.isfinite(corr[i, j]):
                ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", fontsize=7)
    ttl = f"(c) correlations des {len(top)} meilleurs (unites)"
    if eig is not None and eig[-1] > 0:
        ttl += f"\nlambda_max/lambda_min = {eig[0] / eig[-1]:.0f}"
    ax.set_title(ttl); fig.colorbar(im, ax=ax, shrink=0.8)

    # (d) importances ARD
    ax = fig.add_subplot(gs[1, 1])
    if bundle is not None:
        keys = list(bundle["importance"].keys())
        M = np.array([[bundle["importance"][k][nm] for nm in space.names] for k in keys])
        im = ax.imshow(M, vmin=0, vmax=max(1e-9, M.max()), cmap="viridis", aspect="auto")
        ax.set_xticks(range(space.d)); ax.set_xticklabels(space.names, rotation=45)
        ax.set_yticks(range(len(keys))); ax.set_yticklabels(keys, fontsize=7)
        for i in range(len(keys)):
            for j in range(space.d):
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=6, color="w" if M[i, j] < 0.5 * M.max() else "k")
        fig.colorbar(im, ax=ax, shrink=0.8)
        ax.set_title("(d) importances ARD (1/longueur, somme 1)")
    else:
        ax.text(0.5, 0.5, "models absents (fit)", ha="center"); ax.set_title("(d) importances")

    # (e) LOO / Q2
    ax = fig.add_subplot(gs[1, 2])
    if bundle is not None:
        keys = list(bundle["q2"].keys())
        q2 = np.array([bundle["q2"][k] for k in keys])
        ax.barh(range(len(keys)), np.clip(q2, -0.2, 1), color="C0")
        ax.set_yticks(range(len(keys))); ax.set_yticklabels(keys, fontsize=7)
        ax.set_xlim(-0.2, 1); ax.axvline(0.7, color="0.5", ls="--", lw=0.8)
        ax.set_xlabel("Q2 (validation croisee)"); ax.set_title("(e) qualite des emulateurs")
        ax.grid(alpha=0.3, axis="x")
    else:
        ax.axis("off")
    fig.suptitle(f"campagne {a.tag} - {n} points complets / {len(rows)} ; J min {J[0]:.3f} ({comp[0]['id']})", fontsize=11)
    fig.savefig(P.report_png, dpi=130, bbox_inches="tight")
    plt.close(fig)

    # markdown
    md = [f"# Campagne {a.tag} - rapport ({time.strftime('%Y-%m-%d %H:%M')})", "",
          f"Dossier : `{P.runs}` ; {len(rows)} points, {n} complets ({len(P.designs())} designs), "
          f"J min {J[0]:.3f} ({comp[0]['id']}), mediane {np.median(J):.3f}.", "",
          f"![figure]({os.path.basename(P.report_png)})", "", "## Meilleurs points", ""]
    jcols = [c for c in header if c.startswith("J_") and c.count("_") == 1]
    hdr = ["rang", "id", "round"] + space.names + ["J"] + [c for c in jcols]
    body = []
    for i, r in enumerate(comp[:10]):
        body.append([i + 1, r["id"], int(r["round"])] + [f"{r[k]:.4g}" for k in space.names] + [f"{r['J']:.3f}"]
                    + [f"{r[c]:.3f}" if math.isfinite(r[c]) else "-" for c in jcols])
    md += [_md_table(hdr, body), ""]
    if bundle is not None:
        md += ["## Emulateurs (GP Matern 5/2 ARD + nugget)", ""]
        hdr = ["observable", "n", "ecartes", "nugget sd", "LOO rmse", "Q2"] + (["ET rmse", "ET Q2"] if bundle.get("et") else [])
        body = []
        for k, m in bundle["emu"].models.items():
            row = [k, len(m.y_), m.n_dropped, f"{m.nugget_std:.4f}", f"{bundle['loo'][k]:.4f}", f"{bundle['q2'][k]:.3f}"]
            if bundle.get("et"):
                e = bundle["et"].get(k, (float("nan"), float("nan")))
                row += [f"{e[0]:.4f}", f"{e[1]:.3f}"]
            body.append(row)
        md += [_md_table(hdr, body), "", "## Importances ARD (1/longueur, somme 1)", ""]
        body = [[k] + [f"{bundle['importance'][k][nm]:.3f}" for nm in space.names] for k in bundle["importance"]]
        md += [_md_table(["observable"] + space.names, body), ""]
    md += [f"## Correlations des {len(top)} meilleurs points (coordonnees normalisees)", ""]
    body = [[nm] + [f"{corr[i, j]:.2f}" if np.isfinite(corr[i, j]) else "-" for j in range(space.d)] for i, nm in enumerate(space.names)]
    md += [_md_table([""] + space.names, body), ""]
    if eig is not None:
        md += ["valeurs propres de la covariance : " + ", ".join(f"{v:.3g}" for v in eig)
               + (f" ; lambda_max/lambda_min = {eig[0] / eig[-1]:.0f} (> 100 = direction molle, A4 6.3)" if eig[-1] > 0 else ""), ""]
    if mapj is not None:
        md += ["## MAP du substitut (APSO)", "",
               f"J_hat* = {mapj['J']:.4f} (Monte-Carlo {mapj['J_mc_mean']:.3f} +- {mapj['J_mc_std']:.3f}), "
               f"{mapj['apso']['n_seeds']} graines x {mapj['apso']['particles']} particules x {mapj['apso']['iters']} generations.", "",
               _md_table(["parametre", "x*", "u*"], [[nm, f"{mapj['x'][nm]:.5g}", f"{u:.3f}"] for nm, u in zip(mapj["names"], mapj["u"])]), "",
               "```", format_detail_bts(mapj["detail"]), "```", "",
               "Meilleurs distincts (distance > 0.1) :", "",
               _md_table(["J", *mapj["names"]], [[f"{t['J']:.4f}"] + [f"{t['x'][nm]:.4g}" for nm in mapj["names"]] for t in mapj["top"]]), ""]
    if curve_note:
        md += ["Courbes non tracees : " + " ; ".join(curve_note), ""]
    open(P.report_md, "w", encoding="ascii", errors="replace").write("\n".join(md))
    print(f"[report] -> {P.report_png}\n[report] -> {P.report_md}")
    return P.report_png


# ---------------------------------------------------------------------------
# Test synthetique
# ---------------------------------------------------------------------------
U_STAR_TEST = np.array([0.55, 0.40, 0.60, 0.35, 0.50])


def synth_triax(u, conf, T, rng, u_star=U_STAR_TEST, noise=1.0, unbroken=False):
    """Fonction analytique lisse u -> observables d un triaxial (cle du _metrics.json), egales
    aux cibles en u = u_star (J = 2 (5/15)^2 au minimum : rmse plancher 5 MPa)."""
    du = np.asarray(u, dtype=float) - u_star
    tg = T[str(int(conf))]
    s = (float(conf) - 20.0) / 30.0
    N = lambda sd: rng.normal(0.0, sd * noise)  # noqa: E731
    # zone 'bloquee' : resistances trop hautes, l eprouvette ne casse pas dans la fenetre
    # (q depasserait 1,5 x le pic experimental) -> run non rompu, comme A4 6.1
    if 0.5 * du[0] + 0.7 * du[1] + (0.4 + 0.5 * s) * du[2] + 0.1 * du[3] > math.log(1.5):
        unbroken = True
    if unbroken:
        q = 8.0 * (1 + N(0.01)); eps = 0.01 * tg["eps_peak"]
        return dict(q_peak=q, eps_peak=100 * eps, E_GPa=85.0, drop=0.0, q_CI=float("nan"), q_CD=float("nan"),
                    CI_frac=float("nan"), CD_frac=float("nan"), ni_peak=0, nb_peak=0, nb_end=0, rmse_curve=250.0,
                    rmse_curve_norm=20.0, d_peak=q / tg["q_peak"] - 1, d_eps_peak=eps / tg["eps_peak"] - 1,
                    d_drop=-tg["chute"], d_CI=float("nan"), d_CD=float("nan"), s3=float(conf))
    q = tg["q_peak"] * math.exp(0.5 * du[0] + 0.7 * du[1] + (0.4 + 0.5 * s) * du[2] + 0.1 * du[3]) * (1 + N(0.01))
    eps = tg["eps_peak"] * math.exp(0.3 * du[1] + 0.6 * du[3] - 0.2 * du[4]) * (1 + N(0.02))
    ci = tg["CI_frac"] + 0.25 * du[2] - 0.15 * du[0] + 0.10 * du[4] + N(0.01)
    cd = tg["CD_frac"] + 0.20 * du[1] + 0.15 * du[4] - 0.10 * du[3] + N(0.01)
    drop = tg["chute"] + 0.30 * du[4] - 0.20 * du[3] + 0.10 * du[0] + N(0.02)
    rmse = (5.0 + 60.0 * float(np.sum(du ** 2))) * (1 + N(0.05))
    nb = max(1, int(round(300 + 200 * du[0] + 100 * du[2] + N(20))))
    return dict(q_peak=q, eps_peak=100 * eps, E_GPa=84.8 + N(1.0), drop=drop, q_CI=ci * q, q_CD=cd * q, CI_frac=ci,
                CD_frac=cd, ni_peak=int(0.5 * nb), nb_peak=int(0.3 * nb), nb_end=nb, rmse_curve=rmse,
                rmse_curve_norm=rmse / 15.0, d_peak=q / tg["q_peak"] - 1, d_eps_peak=eps / tg["eps_peak"] - 1,
                d_drop=drop - tg["chute"], d_CI=ci * q / tg["q_CI"] - 1, d_CD=cd * q / tg["q_CD"] - 1, s3=float(conf))


def synth_bts(u, T, rng, u_star=U_STAR_TEST, noise=1.0):
    du = np.asarray(u, dtype=float) - u_star
    bts = T["BTS"] * math.exp(0.8 * du[0] + 0.2 * du[1]) * (1 + rng.normal(0.0, 0.015 * noise))
    return dict(bts_nominal=bts / 0.89, k_band=0.89, bts=bts, d_bts=bts / T["BTS"] - 1, nb_end=150)


def _mock_runs(runs, tag, ids, U, T, rng, confs=(20, 50), skip=(), fail=(), unbroken=()):
    """Faux dossiers de runs (_run.json + _metrics.json, sans history.csv) et jobs json."""
    jobs = []
    for i, (pid, u) in enumerate(zip(ids, U)):
        for conf in list(confs) + ["BTS"]:
            pref = conf_prefix(conf)
            out = os.path.join(runs, f"out_{tag}_{pid}_{pref}")
            jobs.append(dict(cfg=os.path.join(runs, f"{tag}_{pid}_{pref}.cfg"), out=out, id=pid, conf=conf, seed=None))
            if (pid, pref) in skip:
                continue
            os.makedirs(out, exist_ok=True)
            rc = 1 if (pid, pref) in fail else 0
            json.dump(dict(cfg=jobs[-1]["cfg"], out=out, rc=rc, wall=1.0), open(os.path.join(out, "_run.json"), "w"))
            if rc != 0:
                continue
            m = synth_bts(u, T, rng) if pref == "BTS" else synth_triax(u, conf, T, rng, unbroken=(pid, pref) in unbroken)
            json.dump(m, open(os.path.join(out, "_metrics.json"), "w"), indent=1)
    return jobs


def _test():
    import tempfile
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  [OK]  " if cond else "  [FAIL] ") + msg)
        ok = ok and bool(cond)

    t_all = time.time()
    try:
        T = exp_targets(NU_DEFAULT)
    except FileNotFoundError as e:
        print(f"test impossible : {e}")
        return False
    tmp = tempfile.mkdtemp(prefix="campaign_test_")
    runs = os.path.join(tmp, "runs_t1")
    os.makedirs(runs)
    tag = "t1"
    rng = np.random.default_rng(3)
    space = ParamSpace({k: tuple(v) for k, v in SPACE_DEFAULT.items()})
    print(f"0. jeu synthetique dans {runs}\n   optimum vrai u* = {U_STAR_TEST.tolist()}, "
          f"x* = " + " ".join(f"{n}={v:.4g}" for n, v in zip(space.names, space.from_unit(U_STAR_TEST))))
    n0 = 50
    X0 = space.lhs(n0, seed=11)
    U0 = space.to_unit(X0)
    ids0 = [f"{tag}_{i:03d}" for i in range(n0)]
    P = Paths(tag, runs)
    space.write_csv(P.design0, X0, ids0)
    skip = {(ids0[3], "P050")}
    fail = {(ids0[7], "P020")}
    unbroken = {(ids0[12], "P020")}
    jobs = _mock_runs(runs, tag, ids0, U0, T, rng, skip=skip, fail=fail, unbroken=unbroken)
    json.dump(jobs, open(os.path.join(runs, f"jobs_{tag}.json"), "w"), indent=1)

    # objectif numpy = objective_with_bts sur les metriques synthetiques brutes
    S = objective_settings({})
    o20 = synth_triax(U0[0], 20, T, rng); o50 = synth_triax(U0[0], 50, T, rng); ob = synth_bts(U0[0], T, rng)
    J_ref, _ = eval_objective({20: o20, 50: o50}, ob, S)
    arrs = {20: {k: np.array([o20[k]]) for k in list(TERMS) + ["nb_end"]}, 50: {k: np.array([o50[k]]) for k in list(TERMS) + ["nb_end"]},
            "BTS": {"d_bts": np.array([ob["d_bts"]])}}
    check(abs(float(j_numpy(arrs, S)[0]) - J_ref) < 1e-10, f"j_numpy = objective_with_bts sur un point brut ({J_ref:.4f})")
    o20u = synth_triax(U0[0], 20, T, rng, unbroken=True)
    J_u, _ = eval_objective({20: o20u, 50: o50}, ob, S)
    arrs[20] = {k: np.array([o20u[k]]) for k in list(TERMS) + ["nb_end"]}
    check(abs(float(j_numpy(arrs, S)[0]) - J_u) < 1e-10, f"j_numpy = objective_with_bts avec NaN + non rompu ({J_u:.3f})")

    print("\n1. init + collect")
    main(["init", "--tag", tag, "--runs", runs, "--template", os.path.join(tmp, "absent.cfg")])
    check(os.path.exists(P.cfg) and os.path.exists(P.space), "campaign_t1.json + space_t1.json (PLAN 4a) crees")
    rows = main(["collect", "--tag", tag, "--runs", runs])
    header, rows = read_table(P.base)
    by = {r["id"]: r for r in rows}
    check(len(rows) == n0, f"base_t1.csv : {len(rows)} lignes")
    check(int(by[ids0[3]]["complete"]) == 0 and math.isnan(by[ids0[3]]["P050_q_peak"]) and math.isnan(by[ids0[3]]["J"]),
          "run manquant -> complete 0, NaN, J NaN")
    check(int(by[ids0[7]]["complete"]) == 0 and by[ids0[7]]["n_done"] == 2, "run rc != 0 -> incomplet (2/3 termines)")
    check(int(by[ids0[12]]["complete"]) == 1 and by[ids0[12]]["unbroken_P020"] == 10.0 and math.isnan(by[ids0[12]]["P020_CI_frac"]),
          f"run non rompu -> complet, penalite 10, CI NaN (J = {by[ids0[12]]['J']:.1f})")
    check(sum(int(r["complete"]) for r in rows) == n0 - 2, f"{n0 - 2} points complets")
    n_lock = {pref: int(sum(1 for r in rows if int(r["complete"]) and r[f"{pref}_nb_end"] < 1)) for pref in ("P020", "P050")}
    print(f"   runs bloques (zone de resistances trop hautes + point 12 force) : {n_lock}")
    r0 = by[ids0[0]]
    Jc, _ = eval_objective({20: {k: r0[f"P020_{k}"] for k in list(TERMS) + ["nb_end"]}, 50: {k: r0[f"P050_{k}"] for k in list(TERMS) + ["nb_end"]}},
                           {"d_bts": r0["BTS_d_bts"]}, S)
    check(abs(Jc - r0["J"]) < 1e-6 and abs(r0["J_P020"] + r0["J_P050"] + r0["J_BTS"] - r0["J"]) < 1e-6,
          f"J recompose depuis la table = J stocke ({r0['J']:.3f}), somme des confinements + BTS")
    check(all(os.path.exists(os.path.join(resolve_path(j["out"]), "_metrics.json")) for j in jobs if run_done(resolve_path(j["out"]))),
          "_metrics.json presents (lus, pas recalcules : aucun history.csv)")

    print("\n2. fit (--fast, --rf 20 arbres)")
    t0 = time.time()
    bundle = main(["fit", "--tag", tag, "--runs", runs, "--fast", "--rf", "--rf-trees", "20"])
    print(f"   ({time.time() - t0:.1f} s)")
    check(os.path.exists(P.models), "models_t1.pkl ecrit")
    check(len(bundle["obs_keys"]) == 15, f"{len(bundle['obs_keys'])} observables emulees (7 x 2 + BTS)")
    check(bundle["emu"]["P020_CI_frac"].n_dropped == n_lock["P020"] and bundle["emu"]["P020_lq_peak"].n_dropped == n_lock["P020"],
          f"P020 : {n_lock['P020']} runs bloques ecartes de la regression (CI NaN et pic)")
    check(all(pref in bundle["clf"] for pref in n_lock if n_lock[pref] >= 2), f"classifieurs rompu/bloque ajustes : {sorted(bundle['clf'])}")
    q2 = bundle["q2"]
    key_obs = ["P020_lq_peak", "P050_lq_peak", "P020_eps_peak", "BTS_lbts", "P020_drop"]
    check(all(q2[k] > 0.9 for k in key_obs), "Q2 > 0.9 sur pic, eps, BTS, chute : " + ", ".join(f"{k} {q2[k]:.3f}" for k in key_obs))
    imp = bundle["importance"]
    check(max(imp["BTS_lbts"], key=imp["BTS_lbts"].get) == "ft", "BTS : ft parametre dominant (ARD)")
    check(imp["P050_lq_peak"]["phi"] > imp["P020_lq_peak"]["phi"], "pic : phi plus influent a 50 qu a 20 MPa (ARD)")

    print("\n3. optimize (4 graines x 30 particules x 150 generations)")
    t0 = time.time()
    mapj = main(["optimize", "--tag", tag, "--runs", runs, "--n-seeds", "4", "--particles", "30", "--iters", "150"])
    print(f"   ({time.time() - t0:.1f} s)")
    u_hat = np.array(mapj["u"])
    err = np.max(np.abs(u_hat - U_STAR_TEST))
    check(os.path.exists(P.map), "map_t1.json ecrit")
    check(err <= 0.05, f"optimum retrouve a {err:.3f} (max |du|) du vrai, tolerance 0.05 ; |du| = {np.linalg.norm(u_hat - U_STAR_TEST):.3f}")
    check(mapj["J"] < 1.0, f"J_hat* = {mapj['J']:.3f} < 1")
    P_, cfg_, b_, surr = load_surrogate(tag, runs)
    po = surr.predict_obs(U_STAR_TEST)
    check(set(po) == {20, 50, "BTS"} and abs(po[20]["d_peak"]) < 0.03 and abs(po["BTS"]["d_bts"]) < 0.03,
          f"predict_obs(u*) : d_peak20 {po[20]['d_peak']:+.3f}, d_bts {po['BTS']['d_bts']:+.3f}")
    Jo, _ = surr.j_objective(U_STAR_TEST)
    check(abs(surr.jhat(U_STAR_TEST) - Jo) < 1e-9, f"jhat(u*) = objective_with_bts(predict_obs(u*)) = {Jo:.3f}")
    u_lock = np.array([1.0, 1.0, 1.0, 0.5, 0.5])
    pl = surr.predict_obs(u_lock); pb = surr.p_broken(u_lock[None, :])
    check(any(pl[c]["nb_end"] == 0 for c in (20, 50)) and surr.jhat(u_lock) > surr.jhat(U_STAR_TEST) + 10,
          "coin des resistances hautes : classifieur -> nb_end 0, penalite non rompu dans J_hat "
          + "(P(rompu) " + ", ".join(f"{k} {v[0]:.2f}" for k, v in pb.items()) + ")")

    print("\n4. enrich (8 points, kappa 2, sans template -> CSV seul)")
    t0 = time.time()
    path_r1 = main(["enrich", "--tag", tag, "--runs", runs, "--n", "8", "--kappa", "2", "--particles", "30", "--iters", "100"])
    print(f"   ({time.time() - t0:.1f} s)")
    X1, ids1 = space.read_csv(path_r1)
    U1 = space.to_unit(X1)
    check(os.path.exists(P.design_round(1)) and len(ids1) == 8 and ids1[0] == "t1_r1_000", f"design_t1_r1.csv : {len(ids1)} points, ids t1_r1_NNN")
    dmin_ex = min(_min_dist(u, list(U0)) for u in U1)
    check(dmin_ex >= 0.15 and min_pairwise_distance(U1) >= 0.15, f"distance min 0.15 : existants {dmin_ex:.3f}, entre nouveaux {min_pairwise_distance(U1):.3f}")
    check(not glob.glob(os.path.join(runs, "jobs_t1_r1.json")), "pas de decks (template absent)")
    near = min(np.linalg.norm(u - U_STAR_TEST) for u in U1)
    print(f"   point du lot le plus proche de u* : {near:.3f}")

    print("\n5. deuxieme tour : faux runs du lot r1 -> collect -> fit -> optimize")
    jobs1 = _mock_runs(runs, f"{tag}_r1", ids1, U1, T, rng)
    json.dump(jobs1, open(os.path.join(runs, f"jobs_{tag}_r1.json"), "w"), indent=1)
    main(["collect", "--tag", tag, "--runs", runs])
    _, rows = read_table(P.base)
    n1 = len(ids1)
    check(len(rows) == n0 + n1 and sum(int(r["complete"]) for r in rows) == n0 - 2 + n1 and any(int(r["round"]) == 1 for r in rows),
          f"base : {len(rows)} lignes, lot r1 integre (round = 1)")
    main(["fit", "--tag", tag, "--runs", runs, "--fast"])
    mapj2 = main(["optimize", "--tag", tag, "--runs", runs, "--n-seeds", "4", "--particles", "30", "--iters", "150"])
    u2 = np.array(mapj2["u"])
    err2 = np.max(np.abs(u2 - U_STAR_TEST))
    check(err2 <= 0.05, f"apres enrichissement : optimum a {err2:.3f} du vrai (avant {err:.3f})")
    check(P.next_round() == 2, "prochain lot = r2")

    print("\n6. report")
    main(["report", "--tag", tag, "--runs", runs])
    check(os.path.exists(P.report_png) and os.path.getsize(P.report_png) > 10000, "report_t1.png ecrit")
    md = open(P.report_md).read()
    check("## Meilleurs points" in md and "## Emulateurs" in md and "## MAP du substitut" in md and "Correlations" in md,
          "report_t1.md : tableaux meilleurs points, emulateurs, correlations, MAP")

    # 7. decks par make_decks.py si les templates reels sont presents (ecrit des .cfg, ne lance rien)
    tpl = os.path.join(ROOT, "calib_quick", "q1v070_P050.cfg")
    bts = os.path.join(ROOT, "calib_quick", "bts_v070b.cfg")
    if os.path.exists(tpl):
        print("\n7. enrich avec decks (make_decks.py, template reel, 3 points)")
        main(["init", "--tag", tag, "--runs", runs, "--template", tpl] + (["--bts", bts] if os.path.exists(bts) else [])
             + ["--extra", "seed=777"])
        main(["enrich", "--tag", tag, "--runs", runs, "--n", "3", "--kappa", "1", "--particles", "20", "--iters", "60"])
        jf = os.path.join(runs, f"jobs_{tag}_r2.json")
        check(os.path.exists(jf), "jobs_t1_r2.json ecrit par make_decks.py")
        if os.path.exists(jf):
            j2 = json.load(open(jf))
            n_bts = sum(1 for j in j2 if conf_prefix(j["conf"]) == "BTS")
            check(len(j2) == 3 * (2 + int(os.path.exists(bts))) and all(j["id"].startswith("t1_r2_") for j in j2),
                  f"{len(j2)} decks (3 points x (20, 50{', BTS' if n_bts else ''})), ids t1_r2_NNN")
            cfg0 = open(resolve_path(j2[-1]["cfg"])).read()
            check("seed = 777" in cfg0 and "confiningPressure = 50e6" in cfg0, "deck : --extra seed=777 applique, confinement 50e6")
            check(not any(os.path.exists(os.path.join(resolve_path(j["out"]), "_run.json")) for j in j2), "aucun run lance")
    else:
        print("\n7. (template reel absent : integration make_decks.py non testee)")

    print(f"\nRESULTAT campaign.py : {'tous les tests passent' if ok else 'ECHEC'}  ({time.time() - t_all:.0f} s, {tmp})")
    return ok


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    ap = argparse.ArgumentParser(description="campagne de calibration par substitut (calib_quick)")
    ap.add_argument("--test", action="store_true", help="test synthetique complet")
    sub = ap.add_subparsers(dest="cmd")

    def common(p):
        p.add_argument("--tag", required=True)
        p.add_argument("--runs", default=None, help="dossier des runs (defaut calib_quick/runs_<tag>)")

    p = sub.add_parser("init", help="configuration de la campagne (template, bts, extra, conf, poids)")
    common(p)
    p.add_argument("--template", required=True)
    p.add_argument("--bts", default=None)
    p.add_argument("--extra", nargs="*", default=[])
    p.add_argument("--conf", type=float, nargs="+", default=CONF_DEFAULT)
    p.add_argument("--seeds", type=int, nargs="*", default=[])
    p.add_argument("--nu", type=float, default=NU_DEFAULT, help="nu des decks (cibles eps x (1 - nu^2))")
    p.add_argument("--weights", default=None, help="JSON poids/tolerances (voir objective_settings)")
    p.add_argument("--targets", default=None, help="cibles (defaut extract.TARGETS)")
    p.add_argument("--seuils", default=None, help="seuils CI/CD (defaut extract.SEUILS)")

    p = sub.add_parser("collect", help="runs termines -> base_<tag>.csv")
    common(p)
    p.add_argument("--weights", default=None)
    p.add_argument("--no-cache", action="store_true", help="ne pas ecrire out/_metrics.json")

    p = sub.add_parser("fit", help="GP par observable -> models_<tag>.pkl")
    common(p)
    p.add_argument("--rf", action="store_true", help="ExtraTrees en controle (LOO compare)")
    p.add_argument("--rf-trees", type=int, default=100)
    p.add_argument("--fast", action="store_true", help="2 redemarrages, LOO a hyperparametres figes")

    p = sub.add_parser("optimize", help="MAP du substitut par APSO -> map_<tag>.json")
    common(p)
    p.add_argument("--weights", default=None)
    p.add_argument("--n-seeds", type=int, default=20)
    p.add_argument("--particles", type=int, default=40)
    p.add_argument("--iters", type=int, default=300)
    p.add_argument("--seed0", type=int, default=1000)

    p = sub.add_parser("enrich", help="lot LCB / Kriging believer -> design_<tag>_r<k>.csv + decks")
    common(p)
    p.add_argument("--weights", default=None)
    p.add_argument("--n", type=int, default=8)
    p.add_argument("--kappa", type=float, default=2.0)
    p.add_argument("--dmin", type=float, default=0.15)
    p.add_argument("--n-draws", type=int, default=200)
    p.add_argument("--particles", type=int, default=30)
    p.add_argument("--iters", type=int, default=150)
    p.add_argument("--seed0", type=int, default=2000)
    p.add_argument("--no-decks", action="store_true", help="CSV seul, sans make_decks.py")

    p = sub.add_parser("report", help="figure + tableaux -> report_<tag>.png/.md")
    common(p)
    return ap


def main(argv=None):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = build_parser()
    a = ap.parse_args(argv)
    if a.test:
        return _test()
    if a.cmd is None:
        ap.print_help()
        return None
    return {"init": cmd_init, "collect": cmd_collect, "fit": cmd_fit, "optimize": cmd_optimize,
            "enrich": cmd_enrich, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    if "--test" in sys.argv:
        sys.exit(0 if _test() else 1)
    main()
