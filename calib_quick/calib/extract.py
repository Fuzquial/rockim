#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# extract.py — observables d un run triaxial rockim (scenario tension confine)
# pour la calibration : pic, eps_pic, module, chute, joints rompus, debut de
# rupture (proxy de sigma_ci), courbe q(eps) sur la grille des cibles, RMSE.
#
#   python calib_quick/calib/extract.py out_dir [out_dir ...] [--json]
#
# q, eps : memes definitions que plot_quick.py (offset de consolidation,
# deplacement des mors analytique / H). Cibles : targets_triax_bohus.json
# (grille eps commune, q_mean, q_std) au confinement du deck.
# ---------------------------------------------------------------------------
import argparse
import json
import math
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(HERE))
from plot_quick import parse_cfg, grip_displacement, fit_E  # noqa: E402

TARGETS = os.path.join(ROOT, "..", "rockim_f1", "calib_triax3d", "targets_triax_bohus.json")
SEUILS = os.path.join(ROOT, "..", "..", "CONTINUUM", "calib_bohus_triax", "exp_qc", "seuils_sbm_bohus.json")


def seuils_exp():
    """CI / CD experimentaux (SBM, Martin & Chandler) moyennes par confinement, en q [MPa]"""
    if not os.path.exists(SEUILS):
        return {}
    raw = json.load(open(SEUILS))
    acc = {}
    for t in raw.values():
        k = str(int(round(t["sigma3"])))
        acc.setdefault(k, []).append((t["CI_MPa"], t["CD_MPa"], t["q_peak"]))
    return {k: dict(q_CI=float(np.mean([v[0] for v in a])), q_CD=float(np.mean([v[1] for v in a])),
                    CI_frac=float(np.mean([v[0] / v[2] for v in a])), CD_frac=float(np.mean([v[1] / v[2] for v in a])))
            for k, a in acc.items()}


def load_run(out_dir, cfg_path=None):
    if cfg_path is None:
        run = json.load(open(os.path.join(out_dir, "_run.json"))) if os.path.exists(os.path.join(out_dir, "_run.json")) else {}
        cfg_path = run.get("cfg") or os.path.join(ROOT, "calib_quick", os.path.basename(out_dir).replace("out_", "") + ".cfg")
        if not os.path.isabs(cfg_path):
            cfg_path = os.path.join(ROOT, cfg_path)
    cfg = parse_cfg(cfg_path)
    H = float(cfg.get("H", 0.04))
    hist = np.genfromtxt(os.path.join(out_dir, "history.csv"), delimiter=",", names=True, invalid_raise=False)
    t = hist["t"]; sig = np.abs(hist["sigma"]); nb = hist["nBroken"]
    ni = hist["nInserted"] if "nInserted" in hist.dtype.names else np.full_like(nb, np.nan)
    ok = np.isfinite(t) & np.isfinite(sig)
    t, sig, nb, ni = t[ok], sig[ok], nb[ok], ni[ok]
    pullV = float(cfg.get("pullV", -0.25)); ramp = float(cfg.get("pullRamp", 0.0)); delay = float(cfg.get("pullDelay", 0.0))
    pre = (t > 0.7 * delay) & (t <= delay)
    nb0 = int(nb[pre][-1]) if pre.sum() else 0
    # CORRECTION C0 (critique A5, 2026-09-02) : pendant pullDelay les mors sont
    # bloques (eps_yy = 0) donc sigma_yy = nu/(1-nu) sigma3 (16,7 MPa a 50) et non
    # sigma3 : l etat de depart n est PAS isotrope. Le deviateur comparable a
    # l essai est q = sigma - sigma3, et eps = 0 a l instant ou sigma atteint
    # sigma3 (exact pour un bulk elastique).
    s3 = float(cfg.get("confiningPressure", 0))
    eps_raw = grip_displacement(t, pullV, ramp, delay) / H
    q = (sig - s3) / 1e6
    load = t > delay
    reach = np.where(load & (sig >= s3))[0]
    i0 = int(reach[0]) if len(reach) else int(np.argmax(load))
    eps = eps_raw - eps_raw[i0]
    load = load & (np.arange(len(t)) >= i0)
    ni0 = float(ni[pre][-1]) if pre.sum() and np.isfinite(ni[pre][-1]) else 0.0
    nu = float(cfg.get("nu", 0.25))
    return dict(cfg=cfg, eps=eps[load], q=q[load], nb=nb[load] - nb0, ni=ni[load] - ni0,
                s3=s3 / 1e6, T=float(cfg.get("T", 0)), nj0=nb0, nu=nu,
                eps_factor=(1.0 - nu * nu))   # deformation plane : eps_2D = eps_exp (1 - nu^2) a E physique


def observables(r, targets=None):
    eps, q, nb = r["eps"], r["q"], r["nb"]
    ipk = int(np.nanargmax(q)); qpk = float(q[ipk]); epk = float(eps[ipk])
    E = float(fit_E(eps, q))
    drop = float((qpk - q[-1]) / qpk) if qpk > 0 else float("nan")
    # CI (initiation de la microfissuration) = premiere INSERTION adaptative (le joint
    # nait sur l enveloppe de Mohr-Coulomb : c est l amorcage) ; CD (croissance
    # instable) = premiere rupture complete D >= 1. Experience (SBM) : CI 0,55, CD 0,72 du pic.
    ni = r["ni"]
    i_ci = int(np.argmax(ni >= 1)) if np.isfinite(ni).any() and (ni >= 1).any() else None
    q_CI = float(q[i_ci]) if i_ci is not None else float("nan")
    i1 = int(np.argmax(nb >= 1)) if (nb >= 1).any() else None
    q_onset = float(q[i1]) if i1 is not None and nb[i1] >= 1 else float("nan")
    q_CD = q_onset
    ni_pk = int(ni[ipk]) if np.isfinite(ni[ipk]) else -1
    nb_pk = int(nb[ipk]); nb_end = int(nb[-1])
    # nombre de joints rompus quand q atteint 50 / 80 % du pic (montee)
    def nb_at(frac):
        m = np.where(q[:ipk + 1] >= frac * qpk)[0]
        return int(nb[m[0]]) if len(m) else 0
    obs = dict(q_peak=qpk, eps_peak=100 * epk, E_GPa=E, drop=drop, q_CI=q_CI, q_CD=q_CD, q_onset=q_onset,
               CI_frac=q_CI / qpk if qpk > 0 else float("nan"), CD_frac=q_CD / qpk if qpk > 0 else float("nan"),
               ni_peak=ni_pk, nb_peak=nb_pk, nb_end=nb_end, nb_50=nb_at(0.5), nb_80=nb_at(0.8), s3=r["s3"])
    if targets is not None:
        key = str(int(round(r["s3"])))
        if key in targets["confinements"]:
            tg = targets["confinements"][key]
            # cibles en deformation converties pour la deformation plane a E physique (critique A5, R1)
            f = r.get("eps_factor", 1.0)
            eT = np.array(tg["eps_grid_microstrain"]) * 1e-6 * f; qT = np.array(tg["q_mean_MPa"]); sT = np.array(tg["q_std_MPa"])
            epk_exp = tg["eps_peak_microstrain"] * 1e-6 * f
            m = (eT <= 1.2 * epk_exp) & np.isfinite(qT)
            qi = np.interp(eT[m], eps, q, left=0.0, right=np.nan)
            good = np.isfinite(qi)
            obs["rmse_curve"] = float(np.sqrt(np.nanmean((qi[good] - qT[m][good]) ** 2))) if good.any() else float("nan")
            obs["rmse_curve_norm"] = float(np.sqrt(np.nanmean(((qi[good] - qT[m][good]) / np.maximum(sT[m][good], 5.0)) ** 2))) if good.any() else float("nan")
            obs["d_peak"] = qpk / tg["q_peak_mean_MPa"] - 1
            obs["d_eps_peak"] = epk / epk_exp - 1
            obs["d_drop"] = drop - tg["chute_fraction_moyenne"]
        se = seuils_exp().get(key)
        if se:
            obs["q_CI_exp"] = se["q_CI"]; obs["q_CD_exp"] = se["q_CD"]
            obs["CI_frac_exp"] = se["CI_frac"]; obs["CD_frac_exp"] = se["CD_frac"]
            obs["d_CI"] = q_CI / se["q_CI"] - 1
            obs["d_CD"] = q_CD / se["q_CD"] - 1
    return obs


BTS_EXP = 10.27      # MPa, 4 essais (2P/(pi D t), D 49,4 mm, disque plein), sd 0,98


def bts_observables(out_dir):
    """Bresilien rockim (disque Gmsh a meplats 2 x 20 deg, plateaux) :
    - bts_nominal = max de sigmaT = 2 P / (pi D t) ;
    - k_band = rapport (contrainte au centre / sigma_t nominal) mesure par la jauge
      elastique du solveur dans la bande [0,3 ; 0,8] ft (ligne 'mean sigma_t ... ratio'
      du journal) : la correction des meplats (Wang 2004) mesuree sur le maillage ;
    - bts = k_band x bts_nominal = traction VRAIE au centre, comparable au BTS
      experimental d un disque plein (ou la formule nominale est exacte)."""
    hist = np.genfromtxt(os.path.join(out_dir, "history.csv"), delimiter=",", names=True, invalid_raise=False)
    nom = float(np.nanmax(hist["sigmaT"])) / 1e6 if "sigmaT" in hist.dtype.names else float("nan")
    k = float("nan")
    for cand in ("_log.txt", "run.log"):
        p = os.path.join(out_dir, cand)
        if os.path.exists(p):
            for line in open(p, encoding="utf-8", errors="replace"):
                if "mean sigma_t" in line and "ratio" in line:
                    try:
                        k = float(line.split("ratio")[1].split()[0])
                    except Exception:
                        pass
    bts = k * nom if np.isfinite(k) else nom
    return dict(bts_nominal=nom, k_band=k, bts=bts, d_bts=bts / BTS_EXP - 1,
                nb_end=int(np.nanmax(hist["nBroken"])) if "nBroken" in hist.dtype.names else -1)


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--json", action="store_true", help="ecrit out/_metrics.json")
    a = ap.parse_args()
    targets = json.load(open(TARGETS)) if os.path.exists(TARGETS) else None
    cols = ["s3", "q_peak", "eps_peak", "E_GPa", "drop", "q_CI", "CI_frac", "q_CD", "CD_frac", "ni_peak", "nb_end", "rmse_curve", "d_peak", "d_eps_peak", "d_drop", "d_CI", "d_CD"]
    print(f"{'run':34s} " + " ".join(f"{c:>9s}" for c in cols))
    for d in a.dirs:
        if not os.path.exists(os.path.join(d, "history.csv")):
            print(f"{d:34s} (pas de history.csv)"); continue
        hdr = open(os.path.join(d, "history.csv")).readline()
        if "sigmaT" in hdr:                       # bresilien
            b = bts_observables(d)
            if a.json:
                json.dump(b, open(os.path.join(d, "_metrics.json"), "w"), indent=1)
            print(f"{os.path.basename(d):34s} BTS nominal {b['bts_nominal']:6.2f}  k_band {b['k_band']:5.3f}  "
                  f"BTS corrige {b['bts']:6.2f} MPa (exp {BTS_EXP}, {100*b['d_bts']:+.0f} %)  rompus {b['nb_end']}")
            continue
        r = load_run(d)
        obs = observables(r, targets)
        if a.json:
            json.dump(obs, open(os.path.join(d, "_metrics.json"), "w"), indent=1)
        print(f"{os.path.basename(d):34s} " + " ".join(f"{obs.get(c, float('nan')):9.3f}" if isinstance(obs.get(c), float) else f"{obs.get(c, ''):>9}" for c in cols))


if __name__ == "__main__":
    main()
