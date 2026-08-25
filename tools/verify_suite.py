#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# verify_suite.py — runner de non-régression rockim (patron Yade : références
# en dur, tolérances PAR NATURE de test, un point d'entrée, sortie table +
# code retour). Écrit le 2026-08-11 pour verrouiller la campagne Yan avant
# toute modification du code.
#
#   python3 tools/verify_suite.py --exe build/rockim               # tier fast
#   python3 tools/verify_suite.py --exe build/rockim --tier full
#   python3 tools/verify_suite.py --exe build/rockim --only fdem3d
#   python3 tools/verify_suite.py ... --update-refs refs_linux.json
#
# Règles maison encodées ici :
#   * OMP_NUM_THREADS = 1 par défaut (certification à 1 thread) ;
#   * les contrôles à CHARGE NULLE (pullV = 1e-12) exigent 0 joint cassé —
#     le test le moins cher et le plus discriminant du dépôt ;
#   * dampWork <= 0 partout où le solveur l'imprime (un dashpot ne peut que
#     dissiper) ;
#   * les références Voronoï sont PROPRES À LA PLATEFORME (les distributions
#     std divergent entre libstdc++ et MSVC à graine égale) : le jeu embarqué
#     ici est la baseline Linux/libstdc++ du 2026-08-11.
# ---------------------------------------------------------------------------
import argparse, json, os, re, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CFG = os.path.join(ROOT, "configs")

# ---- extracteurs (regex sur stdout) ---------------------------------------
RX = {
    "err_pct":   r"error = (-?[\d.eE+]+) %",
    "peak_mpa":  r"peak macro stress = (-?[\d.eE+]+) MPa",
    "ucs_mpa":   r"peak axial stress = (-?[\d.eE+]+) MPa",
    "shear_pct": r"([\d.eE+-]+) % shear\)",
    "broken":    r"broken joints\s*[:=]\s*(\d+)",
    "dampwork":  r"joint dashpot work: (-?[\d.eE+-]+) J",
    "yan_int":   r"int f\(D\) dD = ([\d.eE+-]+)",
    "inserted":  r"adaptive insertion: (\d+) / \d+ joints inserted",
    "ratio":     r"measured/expected ratio = ([\d.eE+-]+)",
    "gcact":     r"adaptive contact activation: (\d+) / \d+",
    "gcwork":    r"net work (?:injected )?by general contact: (-?[\d.eE+-]+)",
    "budget":    r"residu\s+: (-?[\d.eE+-]+) J",
    "pot_ke":    r"pot_ke_rel = ([\d.eE+-]+)",
    "pot_mom":   r"pot_mom_rel = ([\d.eE+-]+)",
    "pot3_ke":   r"pot3_ke_rel = ([\d.eE+-]+)",
    "pot3_mom":  r"pot3_mom_rel = ([\d.eE+-]+)",
    "szfac":     r"facteur mean/min/max = ([\d.eE+-]+)",
    # --- viscosite de Yan (eq. 6) et DIF de Yang (eq. 2-3), ajoutes 2026-08-18
    "viscwork":  r"dont visqueux \(2 mu D\) : (-?[\d.eE+-]+) J",
    "edotmed":   r"insertion, mediane ([\d.eE+-]+) /s",
    "difmed":    r"DIF_traction median ([\d.eE+-]+)",
    # --- DIF en schema INTRINSEQUE (strainRateDIFArm = envelope), 2026-08-25
    # difsansdif est LE controle falsifiable de l armement : un joint qui
    # s endommage sans avoir recu son DIF est un joint que l armement a rate.
    # Il doit valoir 0 par construction (le gel a lieu a l instant meme ou le
    # joint quitte sa branche elastique). Une valeur non nulle signalerait que
    # le critere d armement a derive par rapport a la loi de joint.
    "difarmed":   r"armement a l enveloppe\): (\d+) /",
    "difsansdif": r"; (\d+) joints endommages SANS DIF",
    # --- jointDeath : le relais joint -> contact, 2026-08-25 ---------------
    # deadcomp compte les joints morts EN COMPRESSION, c est-a-dire ceux dont
    # la charge normale doit etre reprise par le contact. deadload est cette
    # charge cumulee. Les deux sont des MESURES du relais, pas des verdicts.
    "dead":       r"relais joint->contact: (\d+) joints morts",
    "deadcomp":   r"dont (\d+) EN COMPRESSION",
    "deadload":   r"relais ([\d.eE+-]+) kN",
}

# ---- définition des tests --------------------------------------------------
# check = (extracteur, référence, tolérance ABSOLUE, obligatoire?)
# ref None => seулement présence/PASS ; "PASS" => chercher [PASS] dans stdout.
# Baseline Linux g++13/libstdc++/Eigen 3.4, OMP=1, source rockim_partage_2026-08-11.
TESTS = [
    # --- tier fast : lois, 2D, charge nulle 2D (≈ 1 min) --------------------
    dict(name="selftest_saksala2011", tier="fast", selftest="selftest-saksala2011",
         checks=[]),                                    # exit 0 == PASS (throw sinon)
    dict(name="selftest_dpdfh", tier="fast", selftest="selftest-dpdfh", checks=[]),
    # loi mc (Mohr-Coulomb de Ye et al. 2025) : plateaux plastiques compares
    # aux formules exactes en traction, compression et 3 confinements —
    # exit 0 si l'ecart max est < 1 % (mesure : 0,13 %, et 1e-9 % en compression)
    dict(name="selftest_mc", tier="fast", selftest="selftest-mc", checks=[]),
    # A3 : LE test decisif du contact par potentiel — collision elastique sans
    # frottement, conservation jugee sur dKE (3.7e-12 frontale, 5.1e-7 oblique)
    # et quantite de mouvement machine (3e loi exacte par construction)
    dict(name="selftest_potential2d", tier="fast", selftest="selftest-potential2d",
         checks=[("pot_ke", 0.0, 1e-5, True), ("pot_mom", 0.0, 1e-12, True),
                 ("pass_tag", None, 0, True)]),
    # ... et son miroir 3D (tet-tet, A3 phase 2) : pointe-contre-face puis
    # oblique, dKE 2e-8, quantite de mouvement machine
    dict(name="selftest_potential3d", tier="fast", selftest="selftest-potential3d",
         checks=[("pot3_ke", 0.0, 1e-5, True), ("pot3_mom", 0.0, 1e-12, True),
                 ("pass_tag", None, 0, True)]),
    dict(name="yan_integral", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["jointSoftening = yan", "T = 2e-6"],
         # stdout n'imprime que 6 décimales : la vérification à 1e-12 vit dans
         # tools/yan_point.cpp ; ici on verrouille la valeur imprimée
         checks=[("yan_int", 0.386307294744, 5e-7, True)]),
    dict(name="fem_bar", tier="fast", cfg="verify_fem_bar.cfg",
         checks=[("err_pct", 0.889143, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="dem_tension", tier="fast", cfg="verify_dem_tension.cfg",
         checks=[("err_pct", -0.109232, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="fdem_voronoi_tension", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         checks=[("err_pct", 11.5981, 0.05, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),      # <= 0 (tol vers le +)
    dict(name="zeroload_2d_fan", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # même contrôle sous la loi de décharge en cisaillement de l'eq. 18 : la
    # sécante à l'origine ne doit rien créer au repos non plus (A2, 2026-08-13)
    dict(name="zeroload_2d_origin", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false", "jointShearUnload = origin"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # ... ni l'activation adaptative du contact (A1, 2026-08-13)
    dict(name="zeroload_2d_gcadaptive", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false", "gcActivation = adaptive"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # ... ni le contact par potentiel (A3, 2026-08-13)
    dict(name="zeroload_2d_potential", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false", "contact = potential"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # le cas delaunay est LE cas historique du cliquet (158 joints cassés a
    # charge nulle avant correctif) mais coute ~2.5 min : tier full
    dict(name="zeroload_2d_delaunay", tier="full", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "grainMesh = delaunay", "verifyFt = false"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # effet d'echelle statistique (eq. 42 du rapport DP-DFH, convention des
    # VUMAT sig_k = sigw*(Zeff/V_el)^(1/m)) : maillage UNIFORME, donc le
    # facteur est analytique et unique. Deux Zeff a trois decades d'ecart
    # verifient l'exposant, pas seulement le niveau.
    dict(name="sizeeffect_2d_zeff1mm3", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["jointSizeEffect = 1", "jointSizeEffectM = 24",
               "jointZeff = 1e-9", "T = 1e-9", "verifyFt = false"],
         checks=[("szfac", 0.736079, 1e-5, True)]),
    dict(name="sizeeffect_2d_zeff1cm3", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["jointSizeEffect = 1", "jointSizeEffectM = 24",
               "jointZeff = 1e-6", "T = 1e-9", "verifyFt = false"],
         checks=[("szfac", 0.981577, 1e-5, True)]),
    # ---- viscosite newtonienne de Yan et al. 2023 (leur eq. 6, terme 2 mu D)
    # bulkViscosityXi = 2 EST l amortissement critique de Munjiza
    # mu = 2 h sqrt(E rho) : applique aux chiffres de la Table 1 de Yan
    # (h = 0,75 mm, E = 15 GPa, rho = 1704) il redonne 7583 Pa.s contre les
    # 7600 publies. Le test verrouille DEUX choses : la valeur dissipee, et le
    # SIGNE de la dissipation. La reference est POSITIVE (le solveur imprime -viscWork_) :
    # un terme dissipatif qui INJECTERAIT ressortirait negatif et ferait
    # meme patron que dampWork_.
    dict(name="visc_yan_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "bulkViscosityXi = 2.0",
               "pullV = 0.5", "T = 3e-5"],
         checks=[("viscwork", 0.793632, 1e-4, True), ("broken", 0, 0, True)]),
    # ---- DIF de Yang et al. 2025 (leurs eq. 2-3) --------------------------
    # LE TRIPLET CI-DESSOUS TESTE LA DISCONTINUITE, PAS SEULEMENT UNE VALEUR.
    # Leur eq. 3 imprimee (exposant 0,07) ne se raccorde a aucune de ses deux
    # bornes : elle saute de 1,516 a 1,85 en edot = 1e2 /s. Dans un schema
    # d insertion EXTRINSEQUE ce saut est un attracteur — un joint qui
    # franchit 1e2 voit son seuil bondir de 22 % et cesse de s inserer, si
    # bien que la population inseree s empile JUSTE SOUS 1e2. Mesure du
    # 2026-08-18, meme config a l exposant pres : mediane 99,36 /s (max
    # 99,9988) avec l exposant litteral, contre 40,22 /s avec l exposant
    # 0,1707 deduit de LEUR figure 2b, qui rend la loi continue.
    # Si un jour la premiere mediane s eloigne de 99,4 le collage a disparu :
    # ce test est la pour que ce ne soit pas silencieux.
    dict(name="dif_yang_litteral_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "insertion = adaptive",
               "strainRateDIF = yang", "pullV = 0.5", "T = 3e-5"],
         checks=[("edotmed", 99.3556, 1e-3, True),
                 ("difmed", 1.5157, 1e-4, True)]),
    dict(name="dif_yang_fig2_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "insertion = adaptive",
               "strainRateDIF = yang-fig2", "pullV = 0.5", "T = 3e-5"],
         checks=[("edotmed", 40.2157, 1e-3, True),
                 ("difmed", 1.72029, 1e-4, True)]),
    # Meme loi, chargement x3 : le taux mesure passe de 40 a 316 /s (facteur
    # 7,9) et le DIF sature a son plateau 1,85. C est le taux QUI SUIT LE
    # CHARGEMENT qui est teste ici, pas un niveau isole.
    dict(name="dif_yang_fig2_plateau_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "insertion = adaptive",
               "strainRateDIF = yang-fig2", "pullV = 1.5", "T = 1.2e-5"],
         checks=[("edotmed", 315.75, 1e-2, True),
                 ("difmed", 1.85, 1e-9, True)]),
    # ---- DIF en schema INTRINSEQUE (strainRateDIFArm = envelope) ----------
    # Le TEMOIN de l etude « insertion adaptative vs intrinseque » : meme deck
    # que dif_yang_fig2_2d ci-dessus, SEUL le schema change. Il faut donc que
    # le DIF sache s armer sans instant d insertion — c est ce que ce repere
    # verrouille.
    # difsansdif = 0 est le controle FALSIFIABLE : le gel a lieu a l instant
    # meme ou le joint quitte sa branche elastique, donc aucun joint ne peut
    # s endommager sans DIF. Mesure du 2026-08-25 : 547 joints armes sur 6840,
    # 0 sans DIF. Le premier essai, qui armait sur la contrainte d ELEMENT (le
    # critere de insertionSweep), donnait 0 arme et 24 sans DIF — le joint
    # ecrete la contrainte que ce critere surveille. Ce test existe pour que
    # ce mode d echec ne revienne pas en silence.
    dict(name="dif_intrinseque_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "insertion = intrinsic",
               "strainRateDIF = yang-fig2", "strainRateDIFArm = envelope",
               "pullV = 0.5", "T = 3e-5"],
         checks=[("difarmed", 547, 0, True),
                 ("difsansdif", 0, 0, True),
                 ("edotmed", 57.1995, 1e-3, True),
                 ("difmed", 1.76803, 1e-4, True)]),
    # Charge nulle AVEC le DIF intrinseque arme : personne ne doit franchir
    # l enveloppe, donc AUCUN joint arme, aucun joint casse, dampWork <= 0.
    # C est le patron zeroload applique a la capacite nouvelle (principe II).
    dict(name="zeroload_dif_intrinseque_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12", "insertion = intrinsic",
               "strainRateDIF = yang-fig2", "strainRateDIFArm = envelope"],
         checks=[("difarmed", 0, 0, True),
                 ("difsansdif", 0, 0, True),
                 ("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # ---- jointDeath = damage : le relais joint -> contact de Guo §2.3.3 ----
    # INVARIANCE EN TRACTION PURE. Un joint qui atteint D = 1 en traction ne
    # transmet deja plus rien (f(1) = 0) et ses levres s ecartent : le tuer
    # tout de suite ou attendre dnMax > 3 dnF revient EXACTEMENT au meme. Ce
    # repere verrouille cette equivalence — mesure du 2026-08-25 : err_pct et
    # nombre de casses identiques a la reference fdem_tension, 24 joints morts
    # dont ZERO en compression. Si un jour ce test devie, c est que la mort du
    # joint a cesse d etre neutre la ou elle doit l etre.
    dict(name="jointdeath_tension_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["jointDeath = damage"],
         checks=[("err_pct", -1.38789, 0.01, True),
                 ("dead", 24, 0, True),
                 ("deadcomp", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # Charge nulle AVEC le relais arme : personne ne casse, donc personne ne
    # meurt, donc aucune charge lachee. Patron zeroload (principe II).
    dict(name="zeroload_jointdeath_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12", "jointDeath = damage"],
         checks=[("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # --- tier full : bit-repères 2D longs, adaptatif, 3D grille -------------
    # charge nulle AVEC viscosite : le terme dissipatif ne doit rien casser ni
    # rien injecter quand il n y a pas de chargement (patron zeroload).
    dict(name="visc_zeroload_2d", tier="full", cfg="verify_fdem_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false", "bulkViscosityXi = 2.0"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True),
                 ]),
    dict(name="fdem_tension", tier="full", cfg="verify_fdem_tension.cfg",
         checks=[("err_pct", -1.38789, 0.01, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem_tension_adaptive", tier="full", cfg="verify_fdem_tension.cfg",
         over=["insertion = adaptive"],
         checks=[("err_pct", -4.08193, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="zeroload_2d_grid", tier="full", cfg="verify_fdem_tension.cfg",
         over=["pullV = 1e-12"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # ---- portage 3D des deux capacites (2026-08-19) -----------------------
    # Tier full : 316 s et 121 s en mono-thread, trop lourds pour le fast.
    # visc_yan_3d verrouille la valeur ET le signe : la reference est POSITIVE
    # (le solveur imprime -viscWork_), un terme dissipatif qui injecterait
    # ressortirait negatif et ferait echouer le test.
    dict(name="visc_yan_3d", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["bulkViscosityXi = 2.0"],
         checks=[("viscwork", 0.00316169, 1e-8, True),
                 ("broken", 200, 0, True)]),
    # INVARIANT : sur un maillage UNIFORME, mu gradue par element et mu global
    # (mediane) doivent donner le meme resultat, puisque tous les hEl_ sont
    # egaux. Si ce test diverge de visc_yan_3d, le mu par element est faux —
    # c est le seul controle qui distingue les deux chemins sans avoir besoin
    # d un maillage gradue de reference.
    dict(name="visc_yan_graded_3d", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["bulkViscosityXi = 2.0", "bulkViscosityGraded = 1"],
         checks=[("viscwork", 0.00316169, 1e-8, True),
                 ("broken", 200, 0, True)]),
    # DIF 3D : la mesure du taux passe par maxAbsEigSym3 (forme fermee de
    # Smith 1961) la ou le 2D ecrit un cercle de Mohr a la main. Ce test est
    # le seul qui exerce ce chemin ; sans lui une erreur de spectre 3x3 serait
    # muette (elle ne deplacerait que le seuil d insertion).
    dict(name="dif_yang_fig2_3d", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["insertion = adaptive", "strainRateDIF = yang-fig2"],
         checks=[("edotmed", 1.57523, 1e-4, True),
                 ("difmed", 1.39307, 1e-5, True),
                 ("broken", 198, 0, True)]),
    # Le meme, en schema INTRINSEQUE : le DIF s arme a l enveloppe DU JOINT
    # faute d instant d insertion (principe III — la capacite nait dans les
    # deux solveurs au meme chantier). Mesure du 2026-08-25 : 1800 joints
    # armes sur 11400, 0 endommage sans DIF. A comparer au repere adaptatif
    # ci-dessus, dont il ne differe QUE par le schema : taux median 1,197 /s
    # contre 1,575, DIF 1,373 contre 1,393, casse 200 contre 198.
    # jointDeath = damage en 3D : meme invariance qu en 2D (principe III).
    # Mesure du 2026-08-25 : la reference en `separation` ne tue AUCUN joint
    # (aucun n atteint dnMax > 3 dnF dans la duree du test) tandis que `damage`
    # en tue 200 — et pourtant err_pct et le nombre de casses sont IDENTIQUES
    # au chiffre pres (-4,76678 %, 200). Le relais est donc bien neutre en
    # traction pure, ou les levres s ecartent et ne portent plus rien.
    dict(name="jointdeath_tension_3d", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointDeath = damage"],
         checks=[("err_pct", -4.76678, 0.02, True),
                 ("dead", 200, 0, True),
                 ("deadcomp", 0, 0, True)]),
    dict(name="dif_intrinseque_3d", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["insertion = intrinsic", "strainRateDIF = yang-fig2",
               "strainRateDIFArm = envelope"],
         checks=[("difarmed", 1800, 0, True),
                 ("difsansdif", 0, 0, True),
                 ("edotmed", 1.19703, 1e-4, True),
                 ("difmed", 1.37278, 1e-5, True),
                 ("broken", 200, 0, True)]),
    dict(name="fdem3d_tension", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointXi = 0"],                          # doctrine 2026-08-05 :
         # une vérification de la LOI ne mesure pas la dissipation visqueuse
         checks=[("err_pct", -3.04802, 0.02, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem3d_tension_adaptive", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointXi = 0", "insertion = adaptive"],
         checks=[("err_pct", -0.034309, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="dem3d_tension", tier="full", cfg="verify_dem3d_tension.cfg",
         checks=[("err_pct", -0.0588759, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="zeroload_3d_grid", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["pullV = 1e-12"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem3d_tension_yan", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointXi = 0", "jointSoftening = yan"],
         checks=[("err_pct", -1.10478, 0.05, True), ("pass_tag", None, 0, True),
                 ("yan_int", 0.386307294744, 1e-6, True)]),
    # --- A2 : décharge en cisaillement sur sécante à l'origine (eq. 18) -----
    # forme LITTÉRALE de l'article = origin + jointFrictionScaled = 1.
    # En traction (mode I dominant) l'écart au retour radial doit rester
    # marginal : c'est ce que ces deux repères verrouillent.
    dict(name="fdem_tension_origin", tier="full", cfg="verify_fdem_tension.cfg",
         over=["jointSoftening = yan", "jointFrictionScaled = 1",
               "jointShearUnload = origin"],
         checks=[("err_pct", -1.721, 0.02, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem3d_tension_origin", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointXi = 0", "jointSoftening = yan", "jointFrictionScaled = 1",
               "jointShearUnload = origin"],
         checks=[("err_pct", -1.08522, 0.05, True), ("pass_tag", None, 0, True)]),
    dict(name="zeroload_3d_origin", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["pullV = 1e-12", "jointShearUnload = origin"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # UCS de la fig. 17 de l'article : le SEUL repère du parc dont la rupture
    # soit pilotée par le CISAILLEMENT (~47 % des joints), et donc le seul qui
    # exerce vraiment la loi de mode II — insertion adaptative + f(D) comprises.
    dict(name="ucs_yan_adaptive", tier="full", cfg="../configs_yan/ucs_adap.cfg",
         checks=[("ucs_mpa", 51.0735, 0.15, True), ("broken", 327, 0, True),
                 ("inserted", 1288, 0, True)]),
    dict(name="ucs_yan_origin", tier="full", cfg="../configs_yan/ucs_adap.cfg",
         over=["jointShearUnload = origin", "jointFrictionScaled = 1"],
         checks=[("ucs_mpa", 50.3671, 0.15, True), ("broken", 318, 0, True),
                 ("inserted", 1131, 0, True)]),
    # --- A1 : activation adaptative du contact (gcActivation) ---------------
    # La paire percussion 2D est BIT-IDENTIQUE mode a mode (174 joints, memes
    # energies, fichiers identiques — verifie le 2026-08-13, phase debris
    # comprise, apres le correctif du cache des faces mortes). La paire SHPB
    # (multi-corps, contact necessaire des t = 0) est bit-identique a
    # T = 9e-5 ; au-dela la tempete de fragments entre dans l'enveloppe
    # chaotique (le mode full lui-meme diverge 8x plus tot sous OMP = 2).
    dict(name="percussion_2d", tier="full", cfg="fdem_percussion.cfg",
         checks=[("broken", 174, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    dict(name="percussion_2d_gcadaptive", tier="full", cfg="fdem_percussion.cfg",
         over=["gcActivation = adaptive"],
         checks=[("broken", 174, 0, True), ("gcact", 21, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="shpb_mini", tier="full", cfg="../configs_yan/shpb_mini.cfg",
         over=["T = 9e-5", "frames = 4"],
         checks=[("broken", 812, 0, True)]),
    dict(name="shpb_mini_gcadaptive", tier="full", cfg="../configs_yan/shpb_mini.cfg",
         over=["T = 9e-5", "frames = 4", "gcActivation = adaptive"],
         checks=[("broken", 812, 0, True), ("gcact", 122, 0, True)]),
    # --- A3 : contact par potentiel de Munjiza (contact = potential) --------
    # Conservation au niveau SOLVEUR : assemblage SHPB incassable sans
    # frottement — |gcWork| doit rester au niveau du biais du compteur
    # (mesure -2.6 J/m pour 766 J/m en jeu, 0.3 %), la ou le penalty
    # quasi-plastique par defaut dissiperait massivement.
    dict(name="shpb_elastic_potential", tier="full", cfg="../configs_yan/shpb_mini.cfg",
         over=["T = 9e-5", "frames = 4", "phase.rock.ft = 1e12",
               "phase.rock.cohesion = 1e12", "phase.rock.Gf = 1e6",
               "contactMu = 0", "contact = potential"],
         checks=[("broken", 0, 0, True), ("gcwork", 0.0, 8.0, True)]),
    # SHPB reel en potentiel : reference de plateforme (chaos de broyage
    # verrouille, comme les autres). L'onde incidente colle au penalty a
    # 3e-6 pres ; le disque casse moins (595 vs 812) — loi de contact
    # differente a l'interface, ecart PHYSIQUE assume et documente.
    # ref 595 -> 578 le 13/08 au soir : l'ordre CANONIQUE des paires (N1) a
    # change une fois pour toutes l'ordre des sommes de forces en tempete de
    # fragments (enveloppe chaotique, cf. controle OMP=2). C'est le DERNIER
    # changement d'ordre possible : l'ordre est desormais independant de
    # l'implementation de la detection.
    dict(name="shpb_mini_potential", tier="full", cfg="../configs_yan/shpb_mini.cfg",
         over=["T = 9e-5", "frames = 4", "contact = potential"],
         checks=[("broken", 578, 0, True)]),
    # Percussion 2D en potentiel : verrouille (a) la borne du residu
    # d'injection de la releve de naissance (+27 J/m mesure, garde a 60 —
    # la regression a +936 sans releve doit rester impossible) et (b) le
    # comptage debris de plateforme.
    dict(name="percussion_2d_potential", tier="full", cfg="fdem_percussion.cfg",
         over=["contact = potential"],
         checks=[("broken", 5, 0, True), ("gcwork", 0.0, 60.0, True)]),
    # A3 phase 2 (3D) : charge nulle sous le potentiel tet-tet — c'est CE
    # controle qui a attrape les slivers de tets tangents (5 joints casses au
    # repos avant les gardes plancher-de-volume + fermeture du polyedre)
    # (T reduit a 1e-4 : les slivers existent des le pas 0, la fenetre courte
    # suffit a les attraper et le controle coute ~5 min au lieu de ~11)
    dict(name="zeroload_3d_potential", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["pullV = 1e-12", "T = 1e-4", "contact = potential"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # charge nulle sur le maillage IMPORTE (mesh = file, tets Gmsh) : certifie
    # la machinerie d'import comme les autres controles certifient les leurs
    dict(name="zeroload_3d_filemesh", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["mesh = file",
               "meshFile = " + os.path.join(ROOT, "meshes", "box3d_h45.msh"),
               "pullV = 1e-12", "T = 1e-4"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # --- tier all : Voronoï 3D + charge nulle 3D Voronoï (~40 min) ----------
    dict(name="fdem3d_voronoi_tension", tier="all", cfg="verify_fdem3d_voronoi_tension.cfg",
         over=["jointXi = 0"],
         checks=[("err_pct", 8.45551, 0.10, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="zeroload_3d_voronoi", tier="all", cfg="verify_fdem3d_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, False)]),
    # A1 en 3D sur le maillage Gmsh importe : paire bit-identique a T = 5e-5
    # (6 joints, memes energies), et le x2.4 de mur est deja la (265 -> 111 s)
    dict(name="percussion_3d_gmsh", tier="all", cfg="fdem3d_percussion_base.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "box3d_h45.msh"),
               "T = 5e-5", "frames = 4"],
         checks=[("broken", 6, 0, True)]),
    dict(name="percussion_3d_gmsh_gcadaptive", tier="all", cfg="fdem3d_percussion_base.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "box3d_h45.msh"),
               "T = 5e-5", "frames = 4", "gcActivation = adaptive"],
         checks=[("broken", 6, 0, True), ("gcact", 0, 0, True)]),
    # A3 phase 2 : percussion 3D Gmsh sous le potentiel tet-tet (memes 6
    # joints que la penalite a ce T — le broyage debute a peine et les
    # debris n'ont pas encore d'influence), gcWork garde pres de zero
    dict(name="percussion_3d_gmsh_potential", tier="all", cfg="fdem3d_percussion_base.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "box3d_h45.msh"),
               "T = 5e-5", "frames = 4", "contact = potential"],
         checks=[("broken", 6, 0, True), ("gcwork", 0.0, 1.0, True)]),
    # V1 : DEUX CORPS (groupes physiques Gmsh — insert sphere + bloc, aucun
    # joint inter-corps) au repos : l'insert immobile a 0.5 mm au-dessus ne
    # doit RIEN ressentir — 0 casse, travail de contact exactement nul.
    dict(name="zeroload_bench1_3d", tier="full", cfg="fdem3d_bench1_insert.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "bench1_insert.msh"),
               "T = 2e-5", "frames = 2", "groupVel.insert = 0 0 0"],
         checks=[("broken", 0, 0, True), ("gcwork", 0.0, 0.0, True),
                 ("budget", 0.0, 1e-12, True)]),   # V2/B4 : bilan clos au repos
    # V1 : l'impact complet insert -> roche (2.53 J a -8 m/s). Reference du
    # 2026-08-14 : 12 joints (2 traction / 10 cisaillement), rebond a
    # +5.67 m/s (e = 0.71 — la restitution du potentiel, cf. percussion 2D),
    # gcWork absorbant. ~45 min : le prix du jalon multi-corps.
    dict(name="bench1_insert_impact", tier="all", cfg="fdem3d_bench1_insert.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "bench1_insert.msh")],
         checks=[("broken", 12, 0, True), ("gcwork", 0.0, 1.0, True),
                 ("budget", 0.0, 0.026, True)]),   # V2/B4 : residu <= 1 % de KE0
]

TIERS = {"fast": ["fast"], "full": ["fast", "full"], "all": ["fast", "full", "all"]}


def run_one(exe, t, outroot, env, timeout):
    out = os.path.join(outroot, t["name"])
    if "selftest" in t:
        cmd = [exe, t["selftest"], os.path.join(out + ".csv")]
    else:
        cfg = os.path.join(CFG, t["cfg"])
        if t.get("over"):
            fd, cfg2 = tempfile.mkstemp(suffix=".cfg", prefix=t["name"] + "_")
            with os.fdopen(fd, "w") as f:
                f.write(open(cfg).read())
                f.write("\n# --- verify_suite overrides ---\n")
                for line in t["over"]:
                    f.write(line + "\n")
            cfg = cfg2
        cmd = [exe, cfg, out]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                           env=env)
    except subprocess.TimeoutExpired:
        return dict(name=t["name"], ok=False, dt=time.time() - t0,
                    detail=[f"TIMEOUT après {timeout}s"])
    dt = time.time() - t0
    text = p.stdout + p.stderr
    detail, ok = [], True
    if p.returncode != 0:
        return dict(name=t["name"], ok=False, dt=dt,
                    detail=[f"exit {p.returncode}: {text.strip().splitlines()[-1] if text.strip() else '?'}"])
    for (kind, ref, tol, required) in t["checks"]:
        if kind == "pass_tag":
            n_fail = text.count("[FAIL]")
            if n_fail or "[PASS]" not in text:
                ok = False
                detail.append(f"attendu [PASS], trouvé {n_fail} [FAIL]")
            continue
        m = re.findall(RX[kind], text)
        if not m:
            if required:
                ok = False
                detail.append(f"{kind}: non trouvé dans la sortie")
            else:
                detail.append(f"{kind}: absent (toléré)")
            continue
        val = float(m[-1])
        if kind == "dampwork":                          # <= 0 exigé
            if val > tol:
                ok = False
                detail.append(f"dampWork = {val:g} > 0 : INJECTION d'énergie")
            else:
                detail.append(f"dampWork = {val:g} <= 0 ok")
            continue
        if ref is None:
            detail.append(f"{kind} = {val:g} (référence à fixer)")
            continue
        if abs(val - ref) > tol:
            ok = False
            detail.append(f"{kind} = {val:g}, attendu {ref:g} ± {tol:g}")
        else:
            detail.append(f"{kind} = {val:g} ok (ref {ref:g})")
    return dict(name=t["name"], ok=ok, dt=dt, detail=detail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=os.path.join(ROOT, "build", "rockim"))
    ap.add_argument("--tier", default="fast", choices=list(TIERS))
    ap.add_argument("--only", default=None, help="filtre sur le nom")
    ap.add_argument("--threads", default="1")
    ap.add_argument("--timeout", type=float, default=3600.0)
    ap.add_argument("--json", default=None, help="rapport JSON")
    args = ap.parse_args()

    env = dict(os.environ, OMP_NUM_THREADS=args.threads)
    outroot = tempfile.mkdtemp(prefix="rockim_suite_")
    sel = [t for t in TESTS if t["tier"] in TIERS[args.tier]
           and (args.only is None or args.only in t["name"])]
    print(f"[suite] {len(sel)} tests, exe = {args.exe}, "
          f"OMP_NUM_THREADS = {args.threads}, sorties {outroot}")
    results, allok = [], True
    for t in sel:
        r = run_one(args.exe, t, outroot, env, args.timeout)
        results.append(r)
        allok &= r["ok"]
        print(f"  {'PASS' if r['ok'] else 'FAIL'}  {r['name']:<28s} "
              f"{r['dt']:7.1f}s  " + "; ".join(r["detail"]))
    print(f"[suite] {'TOUT PASSE' if allok else 'ÉCHECS'} "
          f"({sum(r['ok'] for r in results)}/{len(results)})")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
