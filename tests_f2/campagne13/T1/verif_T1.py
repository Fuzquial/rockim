#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# verif_T1.py — controles FALSIFIANTS de la tache T1 (campagne du 13/09/2026)
# au-dela du --selftest de tools/crack_paths.py :
#   1. imp_lib.broken_mask() est un pur refactor : le masque qu'il rend est
#      IDENTIQUE au code en ligne de broken() d'avant (egalite exacte sur les
#      booleens), et broken() rend les memes sommets a 0 (identite binaire) ;
#   2. variante qui DOIT echouer : l'ANCIEN filtre damage >= 0,999 compte plus
#      de facettes que tBreak >= 0 (faux positifs de jointFailRule = majority)
#      et le masque n'est PAS le meme -> le refactor n'a pas ramene l'ancien
#      comportement par accident ;
#   3. le masque tBreak est inclus dans l'ancien (aucune facette rompue perdue) ;
#   4. verification de l'affirmation « le r max du cratere du s = 2,5 est un
#      point aberrant du maillage grossier » : nombre de facettes de surface du
#      noyau a moins de 0,5 mm du r max, aire de la ou des facettes porteuses,
#      et rapport moyenne par 12 secteurs / r max, COMPARE au s = 1 (variante
#      de controle : sur le maillage fin le r max ne doit PAS etre aberrant) ;
#   5. la police effectivement choisie pour la figure (Computer Modern attendu).
# Usage : python tests_f2/campagne13/T1/verif_T1.py
# ---------------------------------------------------------------------------
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bench_impact", "tools"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from imp_lib import broken, broken_mask, read_vtu  # noqa: E402
import crack_paths as cp                           # noqa: E402

OK = [True]


def check(name, cond, detail=""):
    print("   [%s] %s %s" % ("OK " if cond else "ECHEC", name, detail))
    OK[0] = OK[0] and bool(cond)


def ancien_inline(f):
    """Le corps EXACT de broken() avant l'extraction de broken_mask()."""
    if "tBreak" in f:
        return (f["tBreak"] >= 0.0) & (f["bonded"] < 0.5)
    return (f["damage"] >= 0.999) & (f["bonded"] < 0.5)


def filtre_du_12_09(f):
    """Le filtre d'AVANT le correctif du 13/09 (damage, diagnostic §5)."""
    return (f["damage"] >= 0.999) & (f["bonded"] < 0.5)


def main():
    v3 = os.path.join(ROOT, "out_yang2026_v3", "fdem3d_joints_0018.vtu")
    s25 = os.path.join(ROOT, "out_yang_bench_s25_v3P", "fdem3d_joints_0009.vtu")
    print("== 1-3. imp_lib.broken_mask : refactor pur, et variante falsifiante")
    for path, tag in ((v3, "s = 1 trame 18"), (s25, "s = 2,5 trame 9")):
        pts, con, f = read_vtu(path)
        m_new = broken_mask(f)
        m_old = ancien_inline(f)
        m_dmg = filtre_du_12_09(f)
        c1, n1, mo1, P1 = broken(pts, con, f)
        check("%s : masque = code en ligne d'avant (exact)" % tag,
              bool(np.array_equal(m_new, m_old)),
              "(%d facettes)" % int(m_new.sum()))
        check("%s : broken() rend pts[con[masque]] a 0 exactement" % tag,
              bool(np.array_equal(P1, pts[con[m_new]])))
        check("%s : le filtre damage du 12/09 compte PLUS et differe "
              "(variante qui DOIT echouer)" % tag,
              int(m_dmg.sum()) > int(m_new.sum())
              and not np.array_equal(m_dmg, m_new),
              "(damage %d contre tBreak %d : %d faux positifs, %.1f %%)"
              % (int(m_dmg.sum()), int(m_new.sum()),
                 int(m_dmg.sum()) - int(m_new.sum()),
                 100.0 * (int(m_dmg.sum()) - int(m_new.sum())) / int(m_new.sum())))
        check("%s : tBreak inclus dans damage (aucune rompue perdue)" % tag,
              bool((m_new & ~m_dmg).sum() == 0),
              "(%d hors)" % int((m_new & ~m_dmg).sum()))
    print("== 4. le r max du cratere : point aberrant (s = 2,5) ou non (s = 1) ?")
    stat = {}
    for run, frame, tag in ((os.path.join(ROOT, "out_yang_bench_s25_v3P"), -1,
                             "s = 2,5"),
                            (os.path.join(ROOT, "out_yang2026_v3"), 18, "s = 1")):
        d = cp.load_run(run, frame, False)
        res = cp.analyse(d["P"], d["tri"], rcore=6e-3, skin=1e-3)
        g, lab = res["g"], res["lab"]
        surf = np.where((lab == res["noyau"]) & (g["depth"] < 1e-3))[0]
        rv = g["rvmax"][surf]
        rmax = float(rv.max())
        n05 = int((rv > rmax - 0.5e-3).sum())
        port = surf[rv > rmax - 0.5e-3]
        ar = [float(x) * 1e6 for x in g["area"][port]]
        med = float(np.median(g["area"][surf])) * 1e6
        ratio = res["crater"]["rmean_sect"] / rmax if rmax > 0 else 0.0
        # les facettes porteuses partagent-elles le sommet exterieur ?
        P = d["P"][port]
        rp = np.hypot(P[:, :, 0] - cp.CX, P[:, :, 1] - cp.CY)
        tips = [tuple(np.round(P[k, int(np.argmax(rp[k]))] * 1e9).astype(np.int64))
                for k in range(len(port))]
        print("   %s : r max %.3f mm ; facettes de surface du noyau %d ; a "
              "moins de 0,5 mm du max %d ; aires porteuses %s mm2 (mediane "
              "%.3f, rapport max %.2f) ; sommet exterieur commun : %s ; "
              "moyenne 12 secteurs / r max = %.3f"
              % (tag, rmax * 1e3, len(surf), n05,
                 " ".join("%.3f" % x for x in ar), med, max(ar) / med,
                 "oui" if len(set(tips)) == 1 else "non",
                 ratio))
        stat[tag] = dict(n05=n05, ratio=ratio, ar=max(ar), med=med,
                         one_tip=len(set(tips)) == 1, rmax=rmax)
    a, b = stat["s = 2,5"], stat["s = 1"]
    check("s = 2,5 : le r max est un point aberrant (moins de 3 facettes sur "
          "%d l'atteignent, porteuses plus grosses que la mediane, un seul "
          "sommet exterieur)" % 105,
          a["n05"] < 3 and a["ar"] > a["med"] and a["one_tip"],
          "(%d facettes, %.2f x la mediane)" % (a["n05"], a["ar"] / a["med"]))
    check("s = 2,5 : moyenne par secteur tres en dessous du r max (< 0,65), "
          "s = 1 non (> 0,75) : le r max seul n'est pas une mesure de cratere",
          a["ratio"] < 0.65 < 0.75 < b["ratio"],
          "(s = 2,5 %.3f contre s = 1 %.3f)" % (a["ratio"], b["ratio"]))
    print("== 5. police de la figure")
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.font_manager import findfont, FontProperties
    fp = FontProperties(family="serif")
    matplotlib.rcParams["font.serif"] = ["CMU Serif", "STIXGeneral",
                                         "DejaVu Serif"]
    print("   police serif resolue : %s" % os.path.basename(findfont(fp)))
    print("== verif_T1 : %s" % ("TOUT OK" if OK[0] else "ECHEC"))
    return 0 if OK[0] else 1


if __name__ == "__main__":
    sys.exit(main())
