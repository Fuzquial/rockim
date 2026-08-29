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
    # --- pas de temps stable, ajoute 2026-08-29 (chantier A11) -------------
    # Le budget de pas de temps du 3D ignorait la raideur TANGENTIELLE du
    # contact par potentiel, alors que le 2D la prend depuis longtemps. Xiang,
    # Munjiza, Latham & Guises, Eng. Comput. 26(6) (2009) 673-687, p. 677 :
    # « in order to reduce the numerical error for calculation of TANGENTIAL
    # FORCES, the smaller time step is required ». La cle dtBudgetTangential
    # l y fait entrer ; ces deux controles verrouillent SON EFFET et son
    # innocuite au defaut.
    "dt":        r"dt = ([\d.eE+-]+) s",
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
    # --- WP6 : mu de contact residuel post-pulverisation, 2026-08-28 ------
    "ctcpulv":    r"contact residuel\s*: (\d+) evaluations",
    "npulvel":    r"(\d+) elements a D = Dmax",
    "dead":       r"relais joint->contact: (\d+) joints morts",
    "deadcomp":   r"dont (\d+) EN COMPRESSION",
    "deadload":   r"relais ([\d.eE+-]+) kN",
    # gcfric = le travail de FROTTEMENT du contact general. C est l observable
    # centrale du point 2 : chez Yang et al. c est ce poste qui porte 65 % de
    # l energie d impact (32,0 J sur 49,3, ARMA 2024), et il ne peut se remplir
    # que si les joints rompus passent la main au contact.
    "gcfric":     r"dont frottement (-?[\d.eE+-]+) J",
    # --- gcBirth = penalty : le facteur de naissance, 2026-08-26 -----------
    # C est LA mesure qui dit si le re-echelonnement a seulement pu s exercer :
    # il ne vaut 1 que si AUCUN joint mort ne portait de charge (rupture en
    # traction pure). Un facteur different de 1 prouve que le relais a repris
    # une charge reelle. `birthn` compte les paires concernees.
    "birthfac":   r"paires calees a la naissance, facteur moyen ([\d.eE+-]+)",
    "birthn":     r"gcBirth = penalty : (\d+) paires calees",
}

# ---- références RE-MESUREES SUR LES DEUX PLATEFORMES (2026-08-28) ---------
# Six repères de a8732cc/9462177 échouaient sous Linux g++13. Première
# hypothèse (écart de plateforme pur) REFUTEE par la mesure : la suite
# rejouée le 28/08 sur la machine MSVC de Fernando à HEAD donne, pour
# QUATRE d'entre eux, exactement la valeur Linux —
#   gcbirth_ramp -1.24549, gcbirth_penalty -1.18459,
#   dif_continuous edotmed 7.65306, jointfailrule -2.79574
# alors que leurs références d'origine disaient -1.67766, -1.85267,
# 7.65775, -2.78505. Ces quatre références étaient donc PERIMEES : la
# valeur MSVC a changé entre leur enregistrement et aujourd'hui, et a
# convergé vers la valeur Linux (cause non élucidée — build MSVC
# différemment configuré à l'enregistrement, ou référence écrite sans
# rejouer la suite au commit final). Elles sont désormais SERREES sur la
# valeur commune : une dérive de 1 % sera détectée, ce que le médian large
# de la première correction aurait laissé passer.
# DEUX repères montrent un vrai écart de plateforme, conservé en médian :
#   srfilter_none edotmed  Linux 4.01048 / MSVC 4.01625  (0,14 %)
#                 difmed   Linux 1.4697  / MSVC 1.46982
#   bulkmodel_neohooke     Linux -2.04987 / MSVC -2.20202 (7 %)
# WP6 : les deux plateformes concordent à 0,02 % (ctcpulv 20774 / 20778,
# npulvel 481 des deux côtés) — tolérances resserrées en conséquence.
# Les MÉCANISMES restent verrouillés par les invariants entiers
# (broken/dead/birthfac/difarmed), identiques partout.
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
    # ---- WP6 : contactResidualMu (spec 005) -------------------------------
    # Jumeau SANS contact (jeu 50 mm >> vT) : la cle posee ne doit RIEN
    # engager — ctcpulv = 0 exact, npulvel = 0 exact, aucun joint rompu.
    dict(name="wp6_zeroload_2d", tier="fast", cfg="fdem_percussion.cfg",
         over=["W = 0.06", "H = 0.05", "nx = 30", "ny = 25", "T = 3.0e-5",
               "frames = 1", "toolGap = 0.05", "bulkDamage = yang",
               "bulkDamageDelta0 = 1.0e-7", "bulkDamageDeltaF = 2.0e-7",
               "contactResidualMu = 0.2"],
         checks=[("ctcpulv", 0, 0, True), ("npulvel", 0, 0, True),
                 ("broken", 0, 0, True)]),
    # Micro-fonctionnel : jeu nul + seuils reduits (delta0 = 0,1 um) pour
    # pulveriser d emblee sous le disque -> la bascule DOIT s engager.
    # References du 2026-08-28 (conteneur, gcc) : 20 774 evaluations,
    # 481 elements pulverises ; tolerances +-50 % (comptages dependants de
    # la plateforme via l arithmetique flottante, l ordre de grandeur est
    # le verdict — un zero est le seul vrai FAIL).
    dict(name="wp6_pulv_2d", tier="fast", cfg="fdem_percussion.cfg",
         over=["W = 0.06", "H = 0.05", "nx = 30", "ny = 25", "T = 3.0e-5",
               "frames = 1", "toolGap = 0", "bulkDamage = yang",
               "bulkDamageDelta0 = 1.0e-7", "bulkDamageDeltaF = 2.0e-7",
               "contactResidualMu = 0.2"],
         checks=[("ctcpulv", 20776, 100, True),
                 ("npulvel", 481, 5, True), ("broken", 0, 0, True)]),
    # ---- bulkModel = neohookean : la loi de volume de Guo (eq. 2.6) -------
    # T = (mu/J)(B - I) + (lambda/J) ln(J) I, avec l assemblage EXACT
    # P = J T F^-T = R sigma cof(U). La loi est hyperelastique et redonne
    # l elasticite lineaire AU PREMIER ORDRE avec les memes lambda et mu :
    # verifie analytiquement hors solveur, l ecart vaut -0,008 % a eps = 1e-4,
    # -0,82 % a 1 %, puis +9,4 % a -10 % et +59,8 % a -40 % de deformation.
    # C est donc un remplacement CONTINU de la branche co-rotationnelle, qui
    # n en diverge qu aux grandes deformations — celles de la zone broyee.
    # Sur cette traction (deformations de l ordre du %), l ecart mesure est
    # de 0,81 point : coherent avec la table analytique, et le nombre de
    # casses est INCHANGE (24).
    dict(name="bulkmodel_neohooke_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["bulkModel = neohookean"],
         checks=[("err_pct", -2.12595, 0.0876, True),
                 ("broken", 24, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # Charge nulle sous la loi neuve : aucune casse, aucune injection. La
    # configuration de reference est NATURELLE pour cette loi (W(I) = 0 et
    # dW/dF(I) = 0), donc un maillage au repos ne doit rien produire.
    dict(name="zeroload_neohooke_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12", "bulkModel = neohookean"],
         checks=[("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # ---- loi de joint : les deux dernieres conventions de Guo -------------
    # jointElastic = parabolic (son eq. 2.31) et jointDeltaC = guo (son
    # eq. 2.30). Mesures du 2026-08-25, toutes sous jointSoftening = yan
    # puisque c est le chemin ou la parabole vit :
    #   yan seul                    -1,70281 %
    #   + parabolic                 -2,49639 %   (0,79 point)
    #   + parabolic + deltaC guo    -2,33108 %
    # Le nombre de casses est INCHANGE (24) dans les trois cas : ces
    # conventions deplacent la COMPLAISANCE et l energie dissipee par
    # fissure, pas le nombre de fissures.
    dict(name="jointlaw_guo_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["jointSoftening = yan", "jointElastic = parabolic",
               "jointDeltaC = guo"],
         checks=[("err_pct", -2.33108, 0.01, True),
                 ("broken", 24, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # Charge nulle sous la loi de joint complete : rien ne casse, rien
    # n injecte (principe II applique a la capacite neuve).
    dict(name="zeroload_jointlaw_guo_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12", "jointSoftening = yan",
               "jointElastic = parabolic", "jointDeltaC = guo"],
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
    # ---- chantier A11 : dtBudgetTangential, 2026-08-29 --------------------
    # Le 3D ignorait potKt_ dans son budget de pas de temps ; le 2D le prend.
    # SOURCE : Xiang, Munjiza, Latham & Guises, « On the validation of DEM and
    # FEM/DEM models in 2D and 3D », Eng. Comput. 26(6) (2009) 673-687, p. 677
    # — publication d ORIGINE de leur loi de frottement (eq. 8-9) : « in order
    # to reduce the numerical error for calculation of TANGENTIAL FORCES, the
    # smaller time step is required », alors que le meme calcul SANS frottement
    # est stable au pas plus grand. Les auteurs qualifient eux-memes le point
    # d alarmant.
    # PIEGE D UNITES verrouille ici : en 3D potP_ est en Pa et potKt_ en N/m
    # (Fdem3dSolver.cpp, l. 696 et 701). Recopier le max(potP_, potKt_) du 2D
    # aurait choisi potP_, mille fois trop grand a hmin = 1 mm, et divise le
    # pas par ~32 sans erreur visible. Seul potKt_ entre.
    # Les deux controles ci-dessous sont un couple : le premier prouve que le
    # DEFAUT est bit-identique, le second que la cle FAIT quelque chose. L un
    # sans l autre ne prouverait rien.
    dict(name="dtbudget_tangential_defaut_3d", tier="fast",
         cfg="verify_fdem3d_tension.cfg",
         over=["T = 1e-9", "contact = potential", "potTangentFactor = 1.4286",
               "dtBudgetTangential = off"],
         checks=[("dt", 1.30191e-08, 1e-5, True)]),
    dict(name="dtbudget_tangential_on_3d", tier="fast",
         cfg="verify_fdem3d_tension.cfg",
         over=["T = 1e-9", "contact = potential", "potTangentFactor = 1.4286",
               "dtBudgetTangential = on"],
         checks=[("dt", 1.27315e-08, 1e-5, True)]),
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
    # ---- jointQuadrature = midedge : la regle de Guo (Table 2.2) ---------
    # Les deux regles a 3 points coincident EXACTEMENT en chargement uniforme
    # — ce qui les rend invisibles a tout essai de traction directe classique.
    # Elles ne different que la ou l ouverture a un gradient a travers la
    # facette. Ces reperes existent donc pour verrouiller la MECANIQUE de la
    # regle (repartition des forces sur les deux paires, ponderation), pas un
    # effet physique que ces essais ne montrent pas.
    # Mesures du 2026-08-26 : 2D -1,38789 -> -2,72922 %, 3D -4,76678 ->
    # -3,65374 %, nombre de casses INCHANGE (24 et 200).
    dict(name="jointquad_midedge_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["jointQuadrature = midedge"],
         checks=[("err_pct", -2.72922, 0.01, True),
                 ("broken", 24, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="zeroload_jointquad_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12",
               "jointQuadrature = midedge"],
         checks=[("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # ---- les trois conventions relevees dans le CODE de Solidity ----------
    # (ImperialCollegeLondon/solidity-solver-open, LGPL-3.0, lu le 2026-08-26).
    # jointDeltaC = solidity : leur ot = MAXIM(2 op, 3 Gf/ft) et la rupture a
    # op + ot (Y3Dfd.c l. 1099). Contre `guo` (3 Gf/ft depuis zero), la plage
    # d adoucissement s allonge : l essai doit donc casser AUSSI PEU mais plus
    # tard, et l ecart au pic theorique se creuser. Le plancher 2 op mord ici.
    dict(name="jointdeltac_solidity_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["jointSoftening = yan", "jointElastic = parabolic",
               "jointDeltaC = solidity"],
         checks=[("err_pct", -2.33345, 0.01, True),
                 ("broken", 24, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="zeroload_jointdeltac_solidity_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12", "jointSoftening = yan",
               "jointElastic = parabolic", "jointDeltaC = solidity"],
         checks=[("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # jointFailRule = majority : leur `nfail>1` (l. 1175). Le joint 2D n a que
    # deux points, la regle exige donc les DEUX. Elle arme l endommagement PAR
    # POINT : la facette ne meurt plus au premier point casse, ce qui doit
    # RETARDER la rupture et donc casser moins ou plus tard.
    dict(name="jointfailrule_majority_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["jointQuadrature = midedge", "jointFailRule = majority"],
         checks=[("err_pct", -2.79574, 0.01, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="zeroload_jointfailrule_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12",
               "jointQuadrature = midedge", "jointFailRule = majority"],
         checks=[("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # strainRateDIFArm = continuous : leur dpeftdif, recalcule a chaque pas
    # (l. 1448-1456). CONTROLE FALSIFIANT : sans composition, le facteur median
    # doit rester du meme ordre que celui des deux armements a gel — s il
    # derivait vers des valeurs enormes, ce serait la signature du bug que
    # snapBase() previent (le facteur applique en place a chaque pas).
    # Meme essai que dif_intrinseque_2d, au SEUL armement pres : le repere se
    # lit donc en regard du sien (edot 57,1995 / dif 1,76803, gel a
    # l enveloppe). L ecart entre les deux EST la difference entre geler le
    # facteur une fois et le suivre a chaque pas.
    dict(name="dif_continuous_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "insertion = intrinsic",
               "strainRateDIF = yang-fig2", "strainRateDIFArm = continuous",
               "pullV = 0.5", "T = 3e-5"],
         checks=[("edotmed", 7.65306, 1e-3, True),
                 ("difmed", 1.53036, 1e-4, True)]),
    dict(name="zeroload_dif_continuous_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12", "strainRateDIF = yang",
               "strainRateDIFArm = continuous"],
         checks=[("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # ---- gcBirth = penalty : la naissance du contact sur un joint mort -----
    # Y3Did.c l. 915-964. Il FAUT contact = potential (la cle y vit) et
    # jointDeath = damage (sinon rien ne meurt et le relais ne se produit
    # jamais). Le temoin a comparer est gcbirth_ramp_2d, meme deck a la seule
    # cle pres : l ecart entre les deux EST le prix de la rampe a force nulle.
    dict(name="gcbirth_ramp_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["contact = potential", "jointDeath = damage"],
         checks=[("err_pct", -1.24549, 0.01, True),
                 ("dead", 24, 0, True)]),
    dict(name="gcbirth_penalty_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["contact = potential", "jointDeath = damage",
               "gcBirth = penalty"],
         # ATTENTION a la lecture : ici les 24 joints meurent en TRACTION PURE
         # (deadcomp = 0, charge relayee 0 kN/m), donc fDeath = 0 et le facteur
         # de naissance retombe a 1 pour toutes les paires. Ce repere mesure
         # donc la SUPPRESSION DE LA RAMPE (-1,67766 -> -1,85267 %), pas le
         # re-echelonnement. Celui-ci est verrouille par
         # gcbirth_penalty_percussion_2d, ou l indenteur fait mourir des joints
         # COMPRIMES.
         checks=[("err_pct", -1.18459, 0.01, True),
                 ("dead", 24, 0, True),
                 ("birthfac", 1.0, 1e-9, True)]),
    dict(name="zeroload_gcbirth_penalty_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12", "contact = potential",
               "jointDeath = damage", "gcBirth = penalty"],
         checks=[("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    # ---- strainRateFilter = none : le taux BRUT, ce que fait leur code -----
    # Meme essai que dif_intrinseque_2d/dif_continuous_2d, au seul filtre
    # pres. Se lit en regard de dif_continuous_2d (edot 7,65775 filtre) :
    # l ecart est exactement ce que le passe-bas retirait.
    dict(name="srfilter_none_2d", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "insertion = intrinsic",
               "strainRateDIF = yang-fig2", "strainRateDIFArm = continuous",
               "strainRateFilter = none", "pullV = 0.5", "T = 3e-5"],
         checks=[("edotmed", 4.01337, 3.4e-3, True),
                 ("difmed", 1.46976, 1e-4, True)]),
    # LE repere qui exerce vraiment le re-echelonnement. En TRACTION pure les
    # joints meurent sans charge a relayer (fDeath = 0) et le facteur retombe
    # a 1 : gcbirth_penalty_2d ne teste alors que la SUPPRESSION de la rampe.
    # Il faut un indenteur pour que des joints meurent EN COMPRESSION. Mesure
    # du 2026-08-26 : 4 joints morts, 100 % en compression, 651 kN/m relayes,
    # facteur moyen 1,03077 — et le travail de contact passe de 0,213 a
    # 0,274 J/m, dont frottement 0,111 -> 0,116.
    dict(name="gcbirth_penalty_percussion_2d", tier="full",
         cfg="fdem_percussion.cfg",
         over=["contact = potential", "jointDeath = damage",
               "gcBirth = penalty"],
         # LE garde-fou. Le releve de naissance par AIRE n existait pas que
         # pour la douceur : il empeche une INJECTION d energie sur les paires
         # nees en recouvrement (l en-tete de PotHist::aRef mesure +936 J/m
         # SANS releve). gcBirth = penalty le supprime — il faut donc verifier
         # que le bilan reste dissipatif. Le solveur l ecrit lui-meme :
         # « en mode penalty tout positif est une injection ».
         # MESURE 2026-08-26 : residu -0,9174 J/m en `ramp`, -0,994861 en
         # `penalty` — NEGATIF donc dissipatif dans les deux cas, et meme
         # legerement PLUS dissipatif. Aucune injection sur cet essai.
         checks=[("birthfac", 1.03077, 1e-4, True),
                 ("birthn", 260, 0, True),
                 ("broken", 4, 0, True),
                 ("budget", -0.994861, 1e-3, True)]),
    dict(name="zeroload_srfilter_none_2d", tier="fast",
         cfg="verify_fdem_tension.cfg",
         over=["verifyFt = false", "pullV = 1e-12", "strainRateDIF = yang",
               "strainRateDIFArm = continuous", "strainRateFilter = none"],
         checks=[("broken", 0, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="jointquad_midedge_3d", tier="full",
         cfg="verify_fdem3d_tension.cfg",
         over=["jointQuadrature = midedge"],
         checks=[("err_pct", -3.65374, 0.02, True),
                 ("broken", 200, 0, True)]),
    # ---- miroirs 3D des trois conventions de Solidity (2026-08-26) --------
    # La constitution impose 2D ET 3D : les memes lignes ont ete ecrites dans
    # les deux solveurs, les deux doivent etre verrouillees. Tier full parce
    # que verify_fdem3d_tension est long, pas parce qu ils seraient secondaires.
    dict(name="jointdeltac_solidity_3d", tier="full",
         cfg="verify_fdem3d_tension.cfg",
         over=["jointSoftening = yan", "jointElastic = parabolic",
               "jointDeltaC = solidity"],
         # broken = 0 et NON 200 : sous `guo` comme sous `solidity` le joint
         # 3D est assez ductile (kI = 3 contre 2) pour survivre au deplacement
         # impose de cet essai. VERIFIE le 2026-08-26 : `guo` seul donne deja
         # 0 casse et -1,11074 %, `exact` en donne 200 et -1,27035 %. L ecart
         # solidity/guo (-1,10893 contre -1,11074) est celui des 0,2 % de dnF
         # en plus, le plancher 2 dnE ne mordant pas sur ce maillage.
         checks=[("err_pct", -1.10893, 0.02, True),
                 ("broken", 0, 0, True)]),
    dict(name="jointfailrule_majority_3d", tier="full",
         cfg="verify_fdem3d_tension.cfg",
         over=["jointQuadrature = midedge", "jointFailRule = majority"],
         checks=[("err_pct", -2.25444, 0.02, True)]),
    dict(name="zeroload_jointfailrule_3d", tier="full",
         cfg="verify_fdem3d_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false",
               "jointQuadrature = midedge", "jointFailRule = majority"],
         checks=[("broken", 0, 0, True)]),
    # gcBirth = penalty en 3D. MESURE HONNETE du 2026-08-26 : sur cet essai de
    # TRACTION les deux modes donnent le MEME err_pct (-4,75889 %) — les 200
    # joints meurent en traction pure, il n y a aucune charge a relayer et les
    # contacts naissants ne pesent pas sur le pic. Le repere ne verrouille donc
    # pas un ECART mais le FAIT que le mecanisme se declenche : 118 paires
    # calees. Si le relais cessait de s armer, birthn changerait.
    # La discrimination reelle est en 2D, sous indenteur :
    # gcbirth_penalty_percussion_2d.
    dict(name="gcbirth_penalty_3d", tier="full",
         cfg="verify_fdem3d_tension.cfg",
         over=["contact = potential", "jointDeath = damage",
               "gcBirth = penalty"],
         checks=[("err_pct", -4.75889, 0.02, True),
                 ("birthn", 118, 0, True),
                 ("birthfac", 1.0, 1e-9, True)]),
    dict(name="zeroload_gcbirth_penalty_3d", tier="full",
         cfg="verify_fdem3d_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false", "contact = potential",
               "jointDeath = damage", "gcBirth = penalty"],
         checks=[("broken", 0, 0, True)]),
    # strainRateFilter = none en 3D : a lire en regard de dif_continuous_3d
    # (edot 0,094298 filtre contre 0,0746516 brut) — l ecart est exactement ce
    # que le passe-bas retirait.
    dict(name="srfilter_none_3d", tier="full",
         cfg="verify_fdem3d_tension.cfg",
         over=["verifyFt = false", "strainRateDIF = yang-fig2",
               "strainRateDIFArm = continuous", "strainRateFilter = none"],
         checks=[("edotmed", 0.0746516, 1e-5, True),
                 ("difmed", 1.21328, 1e-4, True)]),
    dict(name="dif_continuous_3d", tier="full",
         cfg="verify_fdem3d_tension.cfg",
         over=["verifyFt = false", "strainRateDIF = yang-fig2",
               "strainRateDIFArm = continuous"],
         checks=[("edotmed", 0.094298, 1e-3, True),
                 ("difmed", 1.22399, 1e-4, True)]),
    # bulkModel = neohookean en 3D (principe III). L exposant de l ecart a la
    # forme co-rotationnelle DIFFERE entre dimensions — J^(-2/3) en 3D contre
    # J^(-1/2) en deformation plane — d ou l importance d avoir ce repere dans
    # les deux : un exposant ecrit en dur au lieu de la forme generique
    # cof(U) = J U^-1 passerait le test 2D et casserait celui-ci.
    # Mesure du 2026-08-25 : -4,55926 % contre -4,76678 en co-rotationnel,
    # 200 casses des deux cotes.
    dict(name="bulkmodel_neohooke_3d", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["bulkModel = neohookean"],
         checks=[("err_pct", -4.55926, 0.02, True),
                 ("broken", 200, 0, True)]),
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
    # ---- LE repere du point 2 : le relais ACHEMINE-T-IL la dissipation ? ---
    # Chez Yang et al. (ARMA 2024) le frottement entre fragments porte 65 % de
    # l energie d impact, et il ne peut se remplir que si les joints rompus
    # passent la main au contact. Mesure du 2026-08-25 sur cet UCS, separation
    # -> damage : travail de contact 0,765 -> 24,68 J/m (x32) et sa part de
    # FROTTEMENT 0,0725 -> 2,160 J/m (x30), a bilan d energie inchange
    # (residu 4,7e-13 % de l echelle). C est la demonstration que le relais
    # fonctionne — et cet UCS tourne a contactMu = 0,1 seulement.
    # ⚠️ TOLERANCE LARGE ASSUMEE. Cet essai est post-pic et chaotique : ses
    # comptages entiers divergent deja entre MSVC et la baseline Linux (cf.
    # SUITE_full_MSVC). On ne verrouille donc pas une valeur mais un ORDRE DE
    # GRANDEUR : le frottement doit rester au voisinage de 2 J/m et surtout ne
    # pas retomber vers les 0,07 du mode separation. C est un test de
    # MECANISME, pas de chiffre.
    dict(name="jointdeath_friction_2d", tier="full",
         cfg="../configs_yan/ucs_adap.cfg",
         over=["jointDeath = damage"],
         checks=[("gcfric", 2.16, 1.0, True),
                 ("deadcomp", 162, 60, True)]),
    # ---- jointResidualMu : les DEUX equivalences qui prouvent qu il --------
    # ---- GENERALISE l existant au lieu de le remplacer (principe VIII) -----
    # La cle interpole le coefficient de frottement du PIC tan(frictionDeg)
    # vers un RESIDUEL, par la meme f(D) que la cohesion. Elle doit donc
    # redonner EXACTEMENT les deux comportements deja en place aux deux bouts
    # de son intervalle. Mesure du 2026-08-25 sur cet UCS (frictionDeg = 23,
    # donc tan = 0,4244748162096047), a OMP = 1 :
    #   jointResidualMu = tan(frictionDeg) == defaut
    #   jointResidualMu = 0               == jointFrictionScaled = 1
    #                                        (0,13478 J/m, 50,944 MPa, 361)
    # Les deux egalites sont EXACTES, pas approchees.
    # L observable discriminante est gcfric et non ucs_mpa : entre le defaut
    # (0,0725 J/m) et jointFrictionScaled = 1 (0,13478) le pic ne bouge que de
    # 0,1 MPa, sous la tolerance de plateforme, tandis que le frottement fait
    # un facteur 1,9. Tolerance 0,03 : assez serree pour attraper un retour au
    # defaut, assez lache pour la derive MSVC/Linux de cet essai post-pic.
    dict(name="residualmu_equiv_defaut_2d", tier="full",
         cfg="../configs_yan/ucs_adap.cfg",
         over=["jointResidualMu = 0.4244748162096047"],
         checks=[("gcfric", 0.0724791, 0.03, True),
                 ("ucs_mpa", 51.0735, 0.15, True)]),
    dict(name="residualmu_equiv_scaled_2d", tier="full",
         cfg="../configs_yan/ucs_adap.cfg",
         over=["jointResidualMu = 0"],
         checks=[("gcfric", 0.13478, 0.03, True),
                 ("ucs_mpa", 50.944, 0.15, True)]),
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
    # ---- bench_polyaxial : le repere COMPRESSED-SHEAR-TO-FAILURE ----------
    # (Guo 2014 §3.5 en triaxial equivalent, cf. bench_polyaxial/LISEZMOI.md).
    # LE regime que la suite n'a jamais charge : un joint COMPRIME mene a la
    # rupture en cisaillement — son absence a laisse la plage de mode II
    # cohesion-seule passer 40 tests verts pendant des mois. Mesures :
    #   2026-08-27 (48fd6ba, OMP=2) : t1 coulomb 11 893 / cohesion 0 ;
    #                                 t3 coulomb 7 824 / cohesion 5 104
    #   2026-08-28 (HEAD,    OMP=4) : t1 coulomb 11 773 / cohesion 0 ;
    #                                 t3 coulomb 7 821 / cohesion 5 154
    # Les comptages bougent de ~1 % avec le nombre de threads (chaos FP d'une
    # rupture massive) : fenetres LARGES, le verdict est porte par les
    # invariants — t1/cohesion casse ZERO joint a 33 MPa = 2,1 x son seuil de
    # Mohr-Coulomb (« plus la cohesion est faible, plus la roche est
    # incassable »), t1/coulomb en casse > 5 000, quasi tous en cisaillement.
    # t1 = paire verdict, tier full (~2 x 5-25 min selon threads) ; t3 =
    # confirmation (pics des DEUX lois dans [15 ; 18] MPa), tier all.
    dict(name="polyaxial_t1_coulomb", tier="full",
         cfg=os.path.join(ROOT, "bench_polyaxial", "polyaxial_guo_t1.cfg"),
         checks=[("broken", 11800, 6800, True),     # > 5 000 : la loi agit
                 ("shear_pct", 99.9, 5.0, True),    # cisaillement dominant
                 ("peak_mpa", 23.3, 2.5, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="polyaxial_t1_cohesion", tier="full",
         cfg=os.path.join(ROOT, "bench_polyaxial", "polyaxial_guo_t1.cfg"),
         over=["jointShearRange = cohesion"],       # la derniere cle gagne
         checks=[("broken", 0, 0, True),            # l'invariant du repere
                 ("peak_mpa", 33.15, 3.0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="polyaxial_t3_coulomb", tier="all",
         cfg=os.path.join(ROOT, "bench_polyaxial", "polyaxial_guo_t3.cfg"),
         checks=[("broken", 7820, 2800, True),
                 ("shear_pct", 100.0, 5.0, True),
                 ("peak_mpa", 16.7, 1.5, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="polyaxial_t3_cohesion", tier="all",
         cfg=os.path.join(ROOT, "bench_polyaxial", "polyaxial_guo_t3.cfg"),
         over=["jointShearRange = cohesion"],
         checks=[("broken", 5130, 2600, True),
                 ("peak_mpa", 17.3, 1.5, True),
                 ("dampwork", 0.0, 1e-12, True)]),
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
