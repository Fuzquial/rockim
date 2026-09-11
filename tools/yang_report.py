#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# yang_report.py — dépouillement d'un run d'impact yang2026 (rockim fdem3d)
#
#   python tools/yang_report.py out_dir [journal.log]
#
# Lit history.csv (colonnes trackGroups / gauge / bulkDamage) et, s'il est
# donné, le journal du solveur (dt, résidu d'énergie, mu résiduel, DIF).
# Imprime les sept critères de Yang et al. (fig. 8) tels que rockim les
# mesure, avec le SENS de chaque chiffre, puis les contrôles de santé
# (énergie créée, rupture dans l'acier, joints gelés). Aucun tracé.
# ---------------------------------------------------------------------------
import csv
import re
import sys


def load_history(path):
    with open(path, newline="") as f:
        rd = csv.reader(f)
        head = next(rd)
        rows = [[float(x) for x in r] for r in rd if r]
    col = {h: i for i, h in enumerate(head)}
    return head, col, rows


def main():
    if len(sys.argv) < 2:
        print(__doc__ or "usage: yang_report.py out_dir [journal.log]")
        return 1
    out = sys.argv[1]
    log = sys.argv[2] if len(sys.argv) > 2 else None
    head, col, rows = load_history(out + "/history.csv")
    if not rows:
        print("history.csv vide")
        return 1
    t = [r[col["t"]] for r in rows]
    T = t[-1]
    us = lambda x: x * 1e6

    def series(name):
        return [r[col[name]] for r in rows] if name in col else None

    print("=== run %s : %d lignes, t final %.1f us ===" % (out, len(rows), us(T)))

    # ---- critere 1 : contrainte max au bit (jauge) ------------------------
    szz = series("szz_bit")
    if szz:
        # compression negative dans rockim ; Yang trace la valeur absolue
        imin = min(range(len(szz)), key=lambda i: szz[i])
        print("1. contrainte max au bit (jauge)      : %.1f MPa a %.1f us"
              "   [Yang 9 m/s : ~160 MPa ; onde 1D rho c v/2 = 178 MPa]"
              % (-szz[imin] / 1e6, us(t[imin])))
    # ---- criteres 2-4 : cinematique du bit -------------------------------
    vzb, zb = series("vz_bit"), series("z_bit")
    if vzb:
        imin = min(range(len(vzb)), key=lambda i: vzb[i])
        vind = -vzb[imin]
        # rebond : max de vz APRES la vitesse d indentation max
        after = list(range(imin, len(vzb)))
        imax = max(after, key=lambda i: vzb[i]) if after else imin
        vreb = vzb[imax]
        print("2. vitesse d indentation du bit       : %.2f m/s a %.1f us"
              "   [Yang 9 m/s, texte p. 11 : 5,62 m/s]" % (vind, us(t[imin])))
        print("3. vitesse de rebond du bit           : %.2f m/s a %.1f us"
              "   [Yang 9 m/s : 4,65 m/s, rapport 0,83 ; KE bit 23,84 -> 16,31 J, "
              "bit 1,509 kg ; retournement ~255 us]" % (vreb, us(t[imax])))
        if vind > 0:
            print("   rapport rebond / indentation       : %.3f" % (vreb / vind))
    zi = series("z_insert")
    if zi:
        depth = zi[0] - min(zi)
        print("4. profondeur d indentation max       : %.3f mm"
              "   [Yang 9 m/s, fig. 9b : ~1,0 mm ; l insert part POSE, jeu 0,02 mm]"
              % (depth * 1e3))
    # ---- criteres 5-7 : fissuration (proxies de history.csv) --------------
    nb = series("nBroken")
    npv = series("nPulv")
    dv = series("detachedVol")
    nf = series("nFrag")
    if nb:
        print("5. joints rompus (roche)              : %d  (fragments %d, volume "
              "detache %.1f mm3)" % (nb[-1], nf[-1] if nf else -1,
                                     (dv[-1] if dv else 0) * 1e9))
    if npv:
        print("6. elements pulverises (D = Dmax)     : %d   [Yang fig. 18 : ~360 a "
              "9 m/s, maille 1 mm]" % npv[-1])
    print("7. radiales / cratere                 : voir crater_metrics.py sur le "
          "dernier VTU (geometrie, pas dans history.csv)   [Yang 9 m/s, fig. 10 : "
          "fragments ~0,3 g, radiale ~9-10 mm, rayon de cratere ~7 mm]")

    # ---- bilan d energie -------------------------------------------------
    ke0 = None
    piston = series("vz_piston")
    print("--- energie (J) ---")
    for k in ("eEl", "eJnt", "eGc", "eFric", "eCund", "eLys", "bdWork", "toolKE"):
        s = series(k)
        if s is not None:
            print("   %-8s final %12.4f   min %12.4f   max %12.4f"
                  % (k, s[-1], min(s), max(s)))
    if piston:
        print("   piston vz : %.2f -> %.2f m/s ; bit vz final %.2f m/s"
              % (piston[0], piston[-1], vzb[-1] if vzb else float("nan")))

    # ---- journal du solveur ---------------------------------------------
    if log:
        txt = open(log, errors="replace").read()
        print("--- journal %s ---" % log)
        for pat in (r"dt = [0-9.e+-]+ s, steps = \d+",
                    r"groupContinuum : .*",
                    r"residu\s*:.*",
                    r"joints\s*: .*cohesif.*",
                    r"contact\s*: .*",
                    r"contact residuel\s*:.*",
                    r"DIF intrinseque.*",
                    r"strainRateDIF \(sur.*",
                    r"broken joints\s*:.*",
                    r"joints intra/homo/hetero:.*",
                    r"\[FDEM3D\] corps '.*",
                    r"\*\*\* AVERTISSEMENT \*\*\*.*",
                    r"ABORT.*|abort.*|budget.*depasse.*",
                    r"wall time:.*",
                    r"prof3d.*"):
            for m in re.finditer(pat, txt):
                print("   " + m.group(0).strip()[:200])
    return 0


if __name__ == "__main__":
    sys.exit(main())
