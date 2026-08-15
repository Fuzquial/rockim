# -*- coding: utf-8 -*-
"""Phase 0 — extraction des CIBLES experimentales Red Bohus.

Source : experimental_data_red_bohus.json (Dumoulin et al. 2024,
Geomechanics for Energy and the Environment 40, 100592 ; dataset Zenodo
10.5281/zenodo.10617548) + les .ASC bruts des essais bresiliens.

Produit targets/targets_redbohus.json :
  - scalaires (pic, module, dispersion) par essai et en moyenne,
  - COURBES moyennes reechantillonnees (l'objectif de Ye et al. 2025 porte
    sur m points de la courbe sigma-eps, pas sur le seul pic),
  - l'enveloppe q(sigma3) avec ecarts-types (ponderation de la vraisemblance).

usage : python extract_targets.py [chemin_json] [dossier_brazilian_asc]
"""
import json, os, sys, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "targets", "targets_redbohus.json")
ARCH = (r"C:\Users\fuzquianoalricabi\OneDrive - Université Paris Sciences "
        r"et Lettres\Documents\phd_geothermie")
JSON = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    ARCH, "FDEM", "calib_cdp_bohus", "experimental_data_red_bohus.json")
BZDIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    ARCH, "DEM", "triax_polycristal_3d", "Mechanical tests granites", "Brazilian")

NPT = 40                      # points de la courbe moyenne (objectif Ye)


def modulus(sig, eps, lo=0.4, hi=0.6):
    """Module secant sur la bande [lo, hi] x pic — bande lineaire franche,
    au-dessus de la fermeture des microfissures et sous l'endommagement."""
    ip = int(np.argmax(sig))
    m = (sig > lo * sig[ip]) & (sig < hi * sig[ip]) & (np.arange(len(sig)) < ip)
    if m.sum() < 3:
        return float("nan")
    return float(np.polyfit(eps[m], sig[m], 1)[0])


def _colstat(curves, fn):
    """Statistique par colonne en ignorant les NaN — les essais ne vont pas
    tous aussi loin en post-pic, une colonne peut n'avoir aucune valeur."""
    a = np.array(curves, dtype=float)
    out = []
    for j in range(a.shape[1]):
        col = a[:, j]
        out.append(float(fn(col, axis=0)) if np.any(np.isfinite(col))
                   else None)
    return out


def resample(sig, eps, npt=NPT):
    """Courbe normalisee : eps/eps_pic en abscisse de 0 a 1.5, sigma en MPa.
    Permet de moyenner des essais de raideurs differentes sans les deformer."""
    ip = int(np.argmax(sig))
    if eps[ip] <= 0:
        return None
    x = eps / eps[ip]
    xs = np.linspace(0.0, 1.5, npt)
    keep = np.concatenate(([True], np.diff(x) > 0))       # x croissant
    return np.interp(xs, x[keep], sig[keep], left=0.0, right=np.nan).tolist()


def main():
    d = json.load(open(JSON, encoding="utf-8"))
    T = {"source": {
        "article": "Dumoulin et al. (2024), Geomechanics for Energy and the "
                   "Environment 40, 100592",
        "dataset": "Zenodo 10.5281/zenodo.10617548",
        "json": os.path.basename(JSON),
        "note": "granite Red Bohus, sud-ouest Suede ; 60 % feldspath, "
                "35 % quartz, 5 % biotite (poids)"}}

    # ---- UCS ------------------------------------------------------------
    peaks, mods, curves = [], [], []
    for k, s in d["UCS"].items():
        sig = np.array(s["stress_MPa"], dtype=float)
        eps = np.array(s["eps_axial_microstrain"], dtype=float) * 1e-6
        peaks.append(float(np.max(sig)))
        E = modulus(sig, eps)
        if np.isfinite(E):
            mods.append(E / 1e3)                           # MPa/1 -> GPa
        c = resample(sig, eps)
        if c is not None:
            curves.append(c)
    T["UCS"] = {
        "n": len(peaks), "peaks_MPa": peaks,
        "mean_MPa": float(np.mean(peaks)), "std_MPa": float(np.std(peaks, ddof=1)),
        "median_MPa": float(np.median(peaks)),
        "E_per_test_GPa": mods,
        "curve_x_eps_over_epspeak": np.linspace(0, 1.5, NPT).tolist(),
        "curve_mean_MPa": _colstat(curves, np.nanmean),
        "curve_std_MPa": _colstat(curves, lambda a, axis: np.nanstd(a, axis=axis)),
        "curve_n_valid": [int(np.sum(np.isfinite(np.array(curves)[:, j])))
                          for j in range(NPT)]}

    # ---- TRIAXIAL --------------------------------------------------------
    by_s3 = {}
    for k, s in d["triaxial"].items():
        s3 = float(s["sigma3_MPa"])
        by_s3.setdefault(s3, []).append(float(s["sigma_peak_deviatoric"]))
    T["triaxial"] = {
        "sigma3_MPa": sorted(by_s3),
        "q_peak_mean_MPa": [float(np.mean(by_s3[s])) for s in sorted(by_s3)],
        "q_peak_std_MPa": [float(np.std(by_s3[s], ddof=1)) for s in sorted(by_s3)],
        "n_per_level": [len(by_s3[s]) for s in sorted(by_s3)],
        "calibration_levels_MPa": [20.0, 50.0],
        "validation_levels_MPa": [75.0, 100.0],
        "note": "enveloppe CONCAVE (pente locale 13,9 -> 3,8) : Mohr-Coulomb "
                "lineaire ne peut pas la suivre sur toute la plage — 75 et 100 "
                "sont gardes en PREDICTION pure"}

    # ---- BTS (recalcule depuis les .ASC : le json ne stocke pas les pics) -
    # Serie 2 = Red Bohus (serie 1 = Sidobre, autre granite). Dimensions et
    # conventions de lecture reprises telles quelles de plot_curves_BD.py :
    # separateur TAB, decimale VIRGULE, et pour le fichier contenant « 4 »
    # l'ordre des colonnes est Time / Axial / Radial / Load (Load en 4e).
    dims = {"2-1": (49.4, 24.8), "2-2": (49.4, 24.8),
            "2-3": (49.4, 24.6), "2-4": (50.2, 24.8)}      # D, t en mm
    bts, pmax, names = [], [], []
    for key in sorted(dims):
        f = os.path.join(BZDIR, key + ".ASC")
        if not os.path.exists(f):
            continue
        col = 3 if "4" in key else 1                        # cf. read_data()
        vals = []
        with open(f, "r", encoding="latin-1") as fh:
            fh.readline()                                   # en-tete
            for line in fh:
                p = line.replace(",", ".").split("\t")
                if len(p) <= col:
                    continue
                try:
                    vals.append(float(p[col]))
                except ValueError:
                    continue
        if not vals:
            continue
        P = float(np.max(np.abs(np.array(vals)))) * 1e3      # kN -> N
        D, t = dims[key]
        names.append(key)
        pmax.append(P / 1e3)
        bts.append(2.0 * P / (np.pi * (D * 1e-3) * (t * 1e-3)) / 1e6)
    if bts:
        T["BTS"] = {"n": len(bts), "tests": names,
                    "P_max_kN": pmax, "sigma_t_MPa": bts,
                    "mean_MPa": float(np.mean(bts)),
                    "std_MPa": float(np.std(bts, ddof=1)) if len(bts) > 1 else None,
                    "note": "RECALCULE par 2P/(pi D t) depuis les .ASC — aucun "
                            "fichier de synthese ne donnait ces pics"}
    else:
        T["BTS"] = {"n": 0, "note": "ASC illisibles — a reprendre"}

    # ---- elasticite retenue ---------------------------------------------
    T["elasticity"] = {
        "E_GPa": 77.66, "nu": 0.29,
        "source": "fit des branches de charge des 12 essais triaxiaux "
                  "(CONTINUUM/calib_bohus_triax/README.md) — le plus stable",
        "alternatives": {"litterature_DP_DFH": [52.0, 0.25],
                         "PSO_local_moyen": [57.3, 0.17]},
        "note": "E et nu sont des SORTIES a reproduire par le modele (methode "
                "Bu 2026 / Jiang 2025), pas des entrees figees"}
    T["density_kg_m3"] = 2620.0

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(T, open(OUT, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("ecrit", os.path.normpath(OUT))
    print("  UCS    : %.1f +- %.1f MPa (n=%d)"
          % (T["UCS"]["mean_MPa"], T["UCS"]["std_MPa"], T["UCS"]["n"]))
    print("  BTS    : %s" % ("%.2f +- %.2f MPa (n=%d)"
          % (T["BTS"]["mean_MPa"], T["BTS"]["std_MPa"] or 0, T["BTS"]["n"])
          if T["BTS"]["n"] else "NON EXTRAIT"))
    for s3, q, sd in zip(T["triaxial"]["sigma3_MPa"],
                         T["triaxial"]["q_peak_mean_MPa"],
                         T["triaxial"]["q_peak_std_MPa"]):
        print("  triax  : s3 = %5.1f -> q = %6.1f +- %.1f MPa" % (s3, q, sd))


if __name__ == "__main__":
    main()
