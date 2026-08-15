# -*- coding: utf-8 -*-
"""Extraction PROPRE des courbes experimentales Red Bohus, fidele aux scripts
de depouillement de la these (Mechanical tests granites/*/plot_curves_*.py).

Pourquoi ce script : le JSON experimental_data_red_bohus.json contient les
jauges LOCALES NON FILTREES — elles decrochent en fin d'essai et produisent
des courbes qui reviennent en arriere (non physiques). Les scripts d'origine
tronquent explicitement la courbe au premier decrochage :
  - UC       : |diff2(Local long. strain)| >= 400  ou  |diff2(Stress)| >= 40
  - triaxial : diff2(Ea) <= -25  ou  diff2(Er) <= -100
et tracent la deformation GLOBALE pour les courbes d'ensemble (fig. 5) et les
jauges locales, tronquees, pour les courbes locales (fig. 6).

Sortie : targets/curves_redbohus.json (courbes propres + pics + modules).
"""
import json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "targets",
                                    "curves_redbohus.json"))
ROOT = (r"C:\Users\fuzquianoalricabi\OneDrive - Université Paris Sciences "
        r"et Lettres\Documents\phd_geothermie\DEM\triax_polycristal_3d"
        r"\Mechanical tests granites")
UCDIR = os.path.join(ROOT, "UC")
TXFILE = os.path.join(ROOT, "Triaxial", "Raw data triax tests.xlsx")
UC_TESTS = ["RED_1", "RED_2", "RED_3", "RED_4"]
# essais triaxiaux Red Bohus et leur confinement (plot_curves_triax.py)
TX = {"2_2": 20, "2_3": 20, "2_4": 20, "2_5": 50, "2_6": 50, "2_7": 50,
      "2_8": 75, "2_9": 75, "2_10": 75, "2_11": 100, "2_12": 100, "2_13": 100}
COLS = ["Global strain", "Stress", "Local long. strain", "Local trans. strain"]


def read_uc(test):
    """Retourne (global, local). ATTENTION a la distinction du script
    d'origine : la courbe d'ensemble (fig. 5) est tracee avec la deformation
    GLOBALE et SANS filtre — c'est elle qui porte le pic. Les filtres de
    troncature ne s'appliquent qu'aux JAUGES LOCALES (fig. 6), qui decrochent
    en fin d'essai ; les appliquer au global couperait avant le pic (mesure :
    RED_1 tomberait de 112 a 80 MPa)."""
    f = os.path.join(UCDIR, "Resultats Compression Simple (UCS_%s).xlsm" % test)
    df = pd.read_excel(f, sheet_name="Traitement", usecols=[3, 4, 6, 8],
                       skiprows=2, header=0, names=COLS)
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    loc = df.copy()
    for col, thr in (("Local long. strain", 400.0), ("Stress", 40.0)):
        bad = loc[loc[col].diff(2).abs() >= thr].index
        if len(bad):
            loc = loc.loc[:bad[0]]
    return df, loc


def read_tx(sheet):
    names = ["Sd - Deviator Stress", "Ea - Axial Strain",
             "Er - Radial Strain", "Ev - Volumetric Strain"]
    df = pd.read_excel(TXFILE, sheet_name=sheet, usecols=[0, 1, 3, 4, 5, 6],
                       header=0, skiprows=100)
    # l'en-tete reel commence a la ligne ou la 1re colonne vaut 'Time'
    k = df[df.iloc[:, 0].astype(str).eq("Time")].index
    if len(k) == 0:
        return None
    df = df.iloc[k[0]:, :].reset_index(drop=True)
    df.columns = df.iloc[0]
    df = df.drop([0, 1]).reset_index(drop=True)
    if not all(n in df.columns for n in names):
        return None
    df = df[names].apply(pd.to_numeric, errors="coerce").dropna()
    for col, thr in ((names[1], -25.0), (names[2], -100.0)):
        bad = df[df[col].diff(2) <= thr].index
        if len(bad):
            df = df.loc[:bad[0]]
    return df


def main():
    out = {"note": "courbes filtrees selon plot_curves_UC.py / "
                   "plot_curves_triax.py (troncature au decrochage des jauges)",
           "UC": {}, "triaxial": {}}
    print("=== compression simple (UC) ===")
    for t in UC_TESTS:
        try:
            df, loc = read_uc(t)
        except Exception as e:
            print("  %s : ECHEC (%s)" % (t, e)); continue
        s = df["Stress"].to_numpy()
        eg = np.abs(df["Global strain"].to_numpy()) * 1e-4   # microstrain -> %
        ip = int(np.argmax(s))
        # module et Poisson sur les JAUGES LOCALES tronquees, bande 40-60 %
        sl = loc["Stress"].to_numpy()
        el = np.abs(loc["Local long. strain"].to_numpy()) * 1e-4
        et = np.abs(loc["Local trans. strain"].to_numpy()) * 1e-4
        ipl = int(np.argmax(sl))
        m = (sl > 0.4 * s[ip]) & (sl < 0.6 * s[ip]) & (np.arange(len(sl)) <= ipl)
        E = nu = float("nan")
        if m.sum() > 3:
            E = np.polyfit(el[m] / 100.0, sl[m], 1)[0] / 1e3      # GPa
            nu = abs(np.polyfit(el[m], et[m], 1)[0])
        out["UC"][t] = {"n": len(s), "peak_MPa": float(s[ip]),
                        "eps_global_pct": eg.tolist(),
                        "stress_MPa": s.tolist(),
                        "eps_local_pct": el.tolist(),
                        "stress_local_MPa": sl.tolist(),
                        "E_GPa": None if np.isnan(E) else float(E),
                        "nu": None if np.isnan(nu) else float(nu)}
        print("  %-6s %5d pts  pic %6.1f MPa a eps_g %.3f %%  E %5.1f GPa  nu %.2f"
              % (t, len(s), s[ip], eg[ip], E, nu))

    print("=== triaxiaux ===")
    for sheet, s3 in TX.items():
        try:
            df = read_tx(sheet)
        except Exception as e:
            print("  %-5s : ECHEC (%s)" % (sheet, e)); continue
        if df is None or df.empty:
            print("  %-5s : vide" % sheet); continue
        q = df["Sd - Deviator Stress"].to_numpy()
        ea = df["Ea - Axial Strain"].to_numpy()
        ip = int(np.argmax(q))
        out["triaxial"][sheet] = {
            "sigma3_MPa": s3, "n": len(q), "q_peak_MPa": float(q[ip]),
            "eps_axial_pct": (ea * 1e-4).tolist(), "q_MPa": q.tolist()}
        print("  %-5s s3 = %3d  %5d pts  q_pic %6.1f MPa a eps %.3f %%"
              % (sheet, s3, len(q), q[ip], ea[ip] * 1e-4))

    json.dump(out, open(OUT, "w"), indent=0)
    print("\necrit", OUT)


if __name__ == "__main__":
    main()
