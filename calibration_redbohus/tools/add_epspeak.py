# -*- coding: utf-8 -*-
"""Ajoute la DEFORMATION AU PIC aux observables de la base — sans relancer un
seul calcul : l'information est deja dans les history.csv.

Motivation (observation de Fernando, 2026-08-16) : a raideur identique
(module secant simule 79-85 GPa contre 73-76 GPa mesure), les eprouvettes
simulees cassent a 0,23-0,39 % de deformation axiale contre 0,65-0,67 % dans
les essais. Deux jeux peuvent donner le meme PIC avec des deformations au pic
tres differentes : en ne calibrant que sur les pics, l'emulateur ne les
distingue pas. C'est l'argument de Ye et al. (2025) — calibrer sur la courbe,
pas sur un scalaire.
"""
import csv, glob, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUNS = os.path.join(BASE, "runs")


def eps_peak(run):
    """Deformation de l'eprouvette au pic de contrainte (%)."""
    p = os.path.join(RUNS, run, "history.csv")
    if not os.path.exists(p):
        return ""
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if len(rows) < 10 or "epsSpec" not in rows[0]:
        return ""
    s = np.abs(np.array([float(r["sigma"]) for r in rows]))
    e = np.abs(np.array([float(r["epsSpec"]) for r in rows]))
    return round(float(e[int(np.argmax(s))]) * 100.0, 4)


def main():
    src = os.path.join(BASE, "lhs_results.csv")
    rows = list(csv.DictReader(open(src)))
    for r in rows:
        r["ucs_eps_pk"] = eps_peak("%s_ucs_s4211" % r["tag"])
        r["tx20_eps_pk"] = eps_peak("%s_tx20_s4211" % r["tag"])
    with open(src, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    ok = [r for r in rows if r["tx20_eps_pk"] != ""]
    print("enrichi %s : %d jeux avec eps_pic" % (os.path.basename(src), len(ok)))
    v = np.array([float(r["tx20_eps_pk"]) for r in ok])
    u = np.array([float(r["ucs_eps_pk"]) for r in rows if r["ucs_eps_pk"] != ""])
    print("  eps_pic tx20 simule : %.3f a %.3f %% (mediane %.3f)"
          % (v.min(), v.max(), np.median(v)))
    print("  eps_pic ucs  simule : %.3f a %.3f %% (mediane %.3f)"
          % (u.min(), u.max(), np.median(u)))

    # --- cibles experimentales -------------------------------------------
    d = json.load(open(os.path.join(BASE, "targets", "curves_redbohus.json")))
    ep = [np.array(s["eps_axial_pct"])[int(np.argmax(s["q_MPa"]))]
          for s in d["triaxial"].values() if s["sigma3_MPa"] == 20]
    print("\ncible experimentale tx20 : eps_pic = %.3f +- %.3f %% (n=%d)"
          % (np.mean(ep), np.std(ep, ddof=1), len(ep)))
    # UCS : la deformation GLOBALE inclut la complaisance machine et n'est pas
    # comparable a la simulation ; les jauges LOCALES sont tronquees avant le
    # pic sur certains essais — on ne retient donc PAS de cible eps_pic en UCS.
    for k, s in d["UC"].items():
        sl = np.array(s["stress_local_MPa"]); el = np.array(s["eps_local_pct"])
        ip = int(np.argmax(sl))
        flag = "" if sl[ip] > 0.97 * s["peak_MPa"] else "  (tronquee avant le pic)"
        print("  UC %-6s local : pic %.1f MPa a %.3f %%%s"
              % (k, sl[ip], el[ip], flag))


if __name__ == "__main__":
    main()
