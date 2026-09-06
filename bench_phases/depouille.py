# -*- coding: utf-8 -*-
"""depouille.py — verdict du banc « materiau par phase en fem3d » (2026-09-06).

Cinq controles, chacun avec son critere FALSIFIANT et, quand c'est possible,
la variante qui DOIT echouer :

  F2  NEUTRALITE  — deux phases strictement identiques a la fiche globale :
      history.csv doit etre OCTET-IDENTIQUE au deck sans `phases`. Falsifie le
      cablage : l'indice de phase traverse masse, CFL, penalite, Lysmer,
      viscosite de volume et loi ; un seul produit reassocie et l'egalite
      tombe. Variante qui doit ECHOUER : f3c (contraste) contre le meme f3a.

  F3  BORNE DE REUSS — barre a deux couches en serie, loi elastique. La
      solution est fermee : 1/E_app = f1/E1 + f2/E2. C'est le SEUL controle
      qui distingue une capacite qui MARCHE d'une capacite qui AFFICHE : si
      laws_[e.phase]->stress n'etait pas branche, les champs VTU seraient
      parfaits et E_app vaudrait exactement E1.
      Les mors (bas encastre, haut a x,y bloques) raidissent les extremites ;
      on ne compare donc pas a la theorie nue mais aux DEUX runs homogenes sur
      le MEME maillage, ce qui elimine l'effet de mors au premier ordre :
          1/sigma_serie == f1/sigma_dur + f2/sigma_mou
      (a deplacement impose identique, sigma est proportionnel a E_app).

  F3d COMMUTATIVITE — couches permutees par groupPhase : E_app inchange.

  F4  MASSE ET PAS DE TEMPS — audit sum(m) == sum_p rho_p V_p (dans le code,
      seuil 1e-9) et dt pris sur c_P MAXIMALE des phases.

  F5  VORONOI — pavage, soudure des grains, fractions realisees.

usage : python bench_phases/depouille.py bench_phases
"""
import csv, os, sys


def last_row(out):
    p = os.path.join(out, "history.csv")
    with open(p, newline="") as f:
        rows = list(csv.DictReader(f))
    return rows[-1]


def sig_eps(out):
    r = last_row(out)
    return float(r["sigZZmid"]), float(r["epsAxGrip"])


def main(root):
    ok = True

    def verdict(name, cond, detail):
        nonlocal ok
        if not cond:
            ok = False
        print("[%s] %-34s %s" % ("OK  " if cond else "ECHEC", name, detail))

    # ---- F2 : neutralite octet a octet -----------------------------------
    a = open(os.path.join(root, "out_f3a_homogene_dur", "history.csv"), "rb").read()
    b = open(os.path.join(root, "out_f2_deux_phases_identiques", "history.csv"), "rb").read()
    c = open(os.path.join(root, "out_f3c_serie_reuss", "history.csv"), "rb").read()
    verdict("F2 neutralite (2 phases = globale)", a == b,
            "history.csv identique octet a octet : %s" % (a == b))
    verdict("F2 variante qui DOIT echouer", a != c,
            "le contraste change bien history.csv : %s" % (a != c))

    # ---- F3 : borne de Reuss ---------------------------------------------
    sa, ea = sig_eps(os.path.join(root, "out_f3a_homogene_dur"))
    sb, eb = sig_eps(os.path.join(root, "out_f3b_homogene_mou"))
    sc, ec = sig_eps(os.path.join(root, "out_f3c_serie_reuss"))
    sd, ed = sig_eps(os.path.join(root, "out_f3d_serie_permutee"))
    # Les quatre runs imposent le MEME deplacement de mors, mais leur pas de
    # temps differe (dt = hmin / c_P max des phases) : la derniere ligne ne
    # tombe pas exactement au meme instant, d'ou un ecart de l'ordre de dt/T.
    # C'est sans effet sur E_app, calcule ligne a ligne comme sigma / eps.
    verdict("F3 meme deformation imposee",
            max(abs(ea - ec), abs(ea - eb), abs(ea - ed)) / abs(ea) < 1e-3,
            "epsAxGrip = %.6e (dur) %.6e (mou) %.6e (serie) %.6e (permutee)"
            % (ea, eb, ec, ed))
    Ea, Eb, Ec, Ed = sa / ea, sb / eb, sc / ec, sd / ed
    pred = 1.0 / (0.5 / Ea + 0.5 / Eb)          # Reuss, calibre sur CE maillage
    err = 100.0 * (Ec - pred) / pred
    verdict("F3 borne de Reuss (serie)", abs(err) < 3.0,
            "E_app mesure %.3f GPa contre %.3f GPa attendu (%.2f %%) ; "
            "homogenes %.3f / %.3f GPa"
            % (Ec / 1e9, pred / 1e9, err, Ea / 1e9, Eb / 1e9))
    verdict("F3 le contraste EST vu par la loi", abs(Ec - Ea) / Ea > 0.2,
            "E_app(serie) / E_app(dur) = %.4f — vaudrait 1 exactement si "
            "laws_[e.phase] n'etait pas branche" % (Ec / Ea))
    verdict("F3d commutativite (groupPhase)", abs(Ed - Ec) / Ec < 0.05,
            "permutee %.3f GPa contre %.3f GPa (%.2f %%)"
            % (Ed / 1e9, Ec / 1e9, 100.0 * (Ed - Ec) / Ec))

    # ---- F5 : chemin mesh = voronoi ---------------------------------------
    va = os.path.join(root, "out_f5a_voronoi_1phase", "history.csv")
    vb = os.path.join(root, "out_f5b_voronoi_3phases_egales", "history.csv")
    vc = os.path.join(root, "out_f5c_voronoi_3phases_contrastees", "history.csv")
    if os.path.exists(va):
        ha, hb, hc = (open(p, "rb").read() for p in (va, vb, vc))
        verdict("F5 voronoi : 3 phases egales", ha == hb,
                "history.csv identique octet a octet au deck sans `phases`")
        verdict("F5 voronoi : variante qui DOIT echouer", ha != hc,
                "le contraste quartz/feldspath/biotite change history.csv")
        # F4 : le pas de temps doit suivre la celerite MAXIMALE des phases.
        # E passe de 60 GPa (uniforme) a 83,1 GPa (quartz) : dt doit CHUTER du
        # rapport des c_P, sinon la CFL a ete laissee sur la fiche globale et
        # le schema tourne 18 % trop vite — il ne plante pas, il devient
        # bruyant, et ce bruit se lit comme de la fissuration.
        def dt_of(log):
            for l in open(os.path.join(root, log), encoding="utf-8",
                          errors="replace"):
                if "dt = " in l and "steps" in l:
                    return float(l.split("dt = ")[1].split(" s")[0])
            return 0.0
        d1 = dt_of("log_f5a_voronoi_1phase.txt")
        d3 = dt_of("log_f5c_voronoi_3phases_contrastees.txt")
        att = (60.0 / 83.1) ** 0.5          # rapport des c_P a rho egal
        verdict("F4 CFL sur c_P MAX des phases",
                abs(d3 / d1 - att) < 0.02,
                "dt %.4e -> %.4e (rapport %.4f, attendu %.4f = "
                "sqrt(60/83,1))" % (d1, d3, d3 / d1, att))

    print("\nVERDICT GLOBAL : %s" % ("TOUS LES CONTROLES PASSENT" if ok
                                     else "AU MOINS UN CONTROLE ECHOUE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "bench_phases"))
