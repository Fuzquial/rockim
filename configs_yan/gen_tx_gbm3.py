#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# gen_tx_gbm3.py — triaxiaux 2D FEMDEM sur un GBM a TROIS phases polydisperses
#
#   python configs_yan/gen_tx_gbm3.py            # ecrit tx_gbm3_{0,10,20,40}.cfg
#
# CE QUI EST REPRIS, ET D'OU (regle « chercher l'existant avant de construire ») :
#   * le BANC : configs_yan/tx_adap_20.cfg — eprouvette 36 x 72 mm, plateaux
#     rigides frottants fermant a 0,1 m/s, confinement lateral rampe avant le
#     chargement, arret automatique apres le pic. C'est le banc de Yan, Zheng &
#     Wang (IJRMMS 169, 2023) deja valide dans ce depot (ancre bitid
#     fdem_ucs_yan_adaptive_court).
#   * la MINERALOGIE : configs/fdem_percussion_gbm_fin.cfg — quartz / feldspar /
#     biotite, fractions 0,33 / 0,59 / 0,08, avec E, nu, rho, ft, cohesion,
#     frictionDeg, Gf par mineral. Carte deja calibree, non retouchee.
#   * le SCHEMA : insertion = adaptive + jointSoftening = yan (loi de Yan).
#
# CE QUI EST AJOUTE ICI (les capacites que rockim avait deja mais que ces decks
# n'armaient pas) :
#   * phase.<nom>.grainSize — TAILLE DE GRAIN PROPRE A CHAQUE PHASE, fractions
#     surfaciques conservees (FdemSolver.cpp:2202, diagramme de Laguerre a
#     aires prescrites, Kitagawa-Merigot-Thibert 2019). Contraste choisi d'apres
#     la petrographie d'un granite : feldspath le plus grossier, quartz
#     intermediaire, biotite en paillettes fines.
#   * grainSizeSpread — polydispersite log-normale AU SEIN de chaque phase
#     (exige grainSeeding = random). 0,35 = modere.
#   * gbAlpha* = 0,7 — les joints INTER-granulaires portent 0,7 fois la
#     resistance du plus faible des deux grains : c'est chi_t = chi_G de
#     l'eq. 21 de la note de septembre 2026.
#   * gb.<a>.<b>.* — surcharges PAR PAIRE. Les interfaces impliquant la biotite
#     sont plus faibles encore (clivage basal du mica) ; quartz-quartz reste le
#     contact le plus resistant. Ces cles ECRASENT gbAlpha* sur les paires
#     concernees (piege documente du depot : gbAlpha* est inerte la ou une
#     paire est posee).
#   * grainElemSize 0,0015 -> 0,0011 m : 38 elements par grain de feldspath et
#     62 par grain de quartz, dans la fourchette 35-90 de la regle du
#     2026-09-07 (a 0,0015 on n'avait que 26 el./grain, sous la fourchette).
#
#   * grainMeshRandom = true — NON NEGOCIABLE, et ce n'est pas un choix de ce
#     generateur : DOCUMENTATION_rockim.md §8 regle 4 pose « maillage structure
#     en FDEM = CONDITION D'INVALIDITE (trajets biaises, divergence en phase
#     debris) », et DOC 5.16 bis conclut « tout deck GBM de calibration doit le
#     poser ». Mesure a l'origine de la regle (remarque de F. Uzquiano du
#     2026-09-02, CHANTIER_f2.md) : le Delaunay intra-grain place ses points
#     interieurs sur un RESEAU TRIANGULAIRE — orientations d'aretes R6 = 0,548,
#     pic/creux 18,8, soit PIRE que le mailleur frontal banni (0,34) : trois
#     directions de fissure imposees a l'interieur de chaque grain. Avec
#     grainMeshRandom (Poisson-disc dans le polygone) : R6 = 0,007, pic/creux
#     1,37. Cout : +15 % d'elements, dt -15 %, soit ~+30 %.
#
#
# NOTE MESUREE LE 2026-09-11 — LA CIBLE N EST PAS LA TAILLE REALISEE.
# `phase.<nom>.grainSize` ne cree pas de grains : la tessellation tire UNE
# log-normale de reference, et la cle ne fait que CHOISIR quels grains vont a
# quelle phase (Tessellation.cpp:738-778, affinite log-normale x deficit
# RELATIF d aire, grains pris du plus grand au plus petit). Consequence
# mesuree : une phase MINORITAIRE gagne l arbitrage sur son deficit relatif
# persistant AVANT que les petits grains n arrivent, et recolte donc des gros.
# Il faut SUR-SPECIFIER : cible 0,8 mm pour realiser 3,1 mm sur la biotite.
# Balayage fait (reference 4,5 mm, 144 grains) :
#     cible biotite  3,5 / 2,5 / 1,5 / 0,8 mm
#     realise        5,1 / 4,5 / 3,7 / 3,1 mm   (11 / 14 / 20 / 28 grains)
# Et une phase MAJORITAIRE ne peut pas etre plus grossiere que les autres :
# le feldspath, a 59 % d aire, doit prendre 92 grains sur 144, donc sa taille
# moyenne est bornee par aire/nombre — cible 8, 12 ou 16 mm donnent toutes
# ~5,4 mm. Contraste REALISABLE ici : quartz 6,9 / feldspath 5,4 / biotite 3,1
# (facteur 2,2 sur la phase fine), ce qui est le bon ordre petrographique.
# Un contraste plus fort demanderait un changement de code (appariement par
# rang plutot qu arbitrage glouton), pas une cle de deck.
#
# RESOLUTION : feldspath 38 el./grain, quartz 62 — dans la fourchette 35-90 de
# la regle du 2026-09-07. La biotite n en a que 13 : c est intrinseque (a
# 35 el./grain sur une phase 2,2 fois plus fine il faudrait 49 000 elements et
# dt/3,6, soit x47 sur le cout). Elle ne porte que 8 % du volume.
#
# LECTURE DU JOURNAL : le solveur imprime la distribution de taille REALISEE
# par phase (apres Lloyd et contraction), pas celle demandee. C'est ce
# chiffre-la qu'il faut regarder.
# ---------------------------------------------------------------------------
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- mineralogie : configs/fdem_percussion_gbm_fin.cfg, non retouchee -------
# (fraction, E [Pa], nu, rho, ft [Pa], cohesion [Pa], frictionDeg, Gf [J/m2],
#  grainSize [m] — AJOUTE ICI)
PHASES = {
    "quartz":   dict(fraction=0.33, E=94e9, nu=0.08, rho=2650,
                     ft=13e6, cohesion=30e6, frictionDeg=45, Gf=90,
                     grainSize=0.0060),   # realise 6,9 mm
    "feldspar": dict(fraction=0.59, E=70e9, nu=0.29, rho=2560,
                     ft=10e6, cohesion=25e6, frictionDeg=40, Gf=70,
                     grainSize=0.0120),   # realise 5,4 mm (voir NOTE)
    "biotite":  dict(fraction=0.08, E=34e9, nu=0.25, rho=3050,
                     ft=7e6,  cohesion=15e6, frictionDeg=30, Gf=40,
                     grainSize=0.0008),   # realise 3,1 mm (voir NOTE)
}

# ---- interfaces par paire : le mica clive, le quartz tient ------------------
# Valeurs ABSOLUES (Pa, deg, J/m2) ; elles ecrasent la regle gbAlpha* sur la
# paire concernee. gfShearFactor exige que Gf soit pose sur la meme paire.
PAIRS = [
    # (a, b, ft, cohesion, Gf, frictionDeg)
    ("biotite", "biotite",  2.0e6,  5.0e6, 12.0, 18.0),   # clivage basal
    ("biotite", "feldspar", 3.5e6,  8.0e6, 20.0, 24.0),
    ("biotite", "quartz",   3.5e6,  8.0e6, 20.0, 24.0),
    ("quartz",  "quartz",   9.0e6, 21.0e6, 60.0, 42.0),   # le plus resistant
]

CONFINEMENTS = [0, 10, 20, 40]          # MPa

HEADER = """# ---------------------------------------------------------------------------
# TRIAXIAL 2D FEMDEM — GBM A TROIS PHASES POLYDISPERSES, confinement {conf} MPa
# Genere par configs_yan/gen_tx_gbm3.py (y lire le detail des choix).
#
# Banc      : Yan, Zheng & Wang (IJRMMS 169, 2023), eprouvette 36 x 72 mm,
#             plateaux rigides frottants, confinement lateral rampe puis
#             chargement axial a 0,1 m/s, arret automatique apres le pic.
# Schema    : insertion ADAPTATIVE des cohesifs + loi de Yan (jointSoftening).
# Micro     : Voronoi/Laguerre a TROIS phases, TAILLES DE GRAIN DISTINCTES
#             REALISEES quartz 6,9 > feldspar 5,4 > biotite 3,1 mm (facteur 2,2
#             sur la phase fine), polydispersite log-normale 0,35 au sein de
#             chaque phase, 38 elements par grain de feldspath et 62 de quartz.
#             ATTENTION : les cles `phase.<nom>.grainSize` du deck sont des
#             cibles SUR-SPECIFIEES, pas les tailles obtenues — voir la NOTE du
#             generateur.
# Maillage  : intra-grain Delaunay NON STRUCTURE (grainMeshRandom) — regle 8.4
#             de DOCUMENTATION_rockim.md : un maillage structure en FDEM est une
#             condition d'invalidite.
# Interfaces: gbAlpha* = 0,7 (chi_t = chi_G de l'eq. 21 de la note 2026) et
#             surcharges par paire — les joints au contact de la biotite sont
#             les plus faibles, quartz-quartz le plus resistant.
#
#   rockim_g1fix.exe configs_yan/tx_gbm3_{conf}.cfg out_tx_gbm3_{conf}
#
# A LIRE DANS LE JOURNAL : la distribution de taille REALISEE par phase, le
# nombre de joints intra / homophase / heterophase, et le mode de rupture
# (tensile / shear) — c'est la que se voit l'effet de la microstructure.
# ---------------------------------------------------------------------------
"""


def deck(conf):
    L = [HEADER.format(conf=conf)]
    a = L.append
    a("mode = fdem")
    a("scenario = tension")
    a("loading = platens")
    a("T = 1.0e-2")
    a("frames = 24")
    a("")
    a("# --- eprouvette (banc de Yan) ---------------------------------------")
    a("W = 0.036")
    a("H = 0.072")
    a("thickness = 1.0")
    a("")
    a("# --- microstructure : Voronoi/Laguerre 3 phases polydisperses -------")
    a("mesh = voronoi")
    a("grainSize = 0.0045                # taille de reference (phases la surchargent)")
    a("grainSeeding = random             # exige par grainSizeSpread")
    a("grainJitter = 0.5")
    a("lloydIters = 2")
    a("grainMesh = delaunay")
    a("grainMeshRandom = true          # OBLIGATOIRE en GBM (DOC 5.16 bis, regle 8.4)")
    a("grainElemSize = 0.0011            # 38 el./grain feldspath, 62 quartz (regle 35-90)")
    a("grainSizeSpread = 0.35            # ecart-type de ln(taille) DANS chaque phase")
    a("seed = 4211")
    a("")
    a("# --- materiau de reference (herite par les phases non surchargees) --")
    a("rho = 2600")
    a("E = 70e9")
    a("nu = 0.25")
    a("ft = 10e6")
    a("cohesion = 25e6")
    a("frictionDeg = 40")
    a("Gf = 70")
    a("gfShearFactor = 22.105            # banc de Yan, inchange")
    a("")
    a("# --- les TROIS phases -----------------------------------------------")
    a("phases = " + " ".join(PHASES))
    for nm, p in PHASES.items():
        a("")
        for k in ("fraction", "grainSize", "rho", "E", "nu",
                  "ft", "cohesion", "frictionDeg", "Gf"):
            a("phase.%s.%s = %g" % (nm, k, p[k]))
    a("")
    a("# --- joints de grain : regle generale puis surcharges par paire -----")
    a("gbAlphaTen = 0.7                  # chi_t, eq. 21 de la note 2026")
    a("gbAlphaCoh = 0.7")
    a("gbAlphaGf = 0.7                   # chi_G")
    a("gbAlphaE = 1.0                    # la raideur n'est pas affaiblie")
    a("gbAlphaFric = 1.0")
    for A, B, ft, coh, gf, phi in PAIRS:
        a("")
        a("gb.%s.%s.ft = %g" % (A, B, ft))
        a("gb.%s.%s.cohesion = %g" % (A, B, coh))
        a("gb.%s.%s.Gf = %g" % (A, B, gf))
        a("gb.%s.%s.frictionDeg = %g" % (A, B, phi))
    a("")
    a("# --- schema : insertion adaptative + loi de Yan ----------------------")
    a("insertion = adaptive")
    a("jointSoftening = yan")
    a("insertionPenaltyFactor = 4")
    a("contactMu = 0.1")
    a("dampingLocal = 0.7                # quasi-statique : amortit les ondes")
    a("jointXi = 0.0")
    a("")
    a("# --- chargement et confinement (protocole du banc) -------------------")
    a("pullV = -0.1")
    a("pullRamp = 2e-4")
    a("pullDelay = 6.0e-4                # laisse le confinement s'etablir")
    a("confiningPressure = %ge6" % conf)
    a("confiningRamp = 2e-4")
    a("confineFaces = sides")
    a("confineGaugeTime = 5.0e-4")
    a("ucsStopAfterPeak = true")
    a("ucsStopDelay = 2.0e-4")
    a("gaugeLoFrac = 0.25")
    a("gaugeHiFrac = 0.75")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    for c in CONFINEMENTS:
        p = os.path.join(HERE, "tx_gbm3_%d.cfg" % c)
        io.open(p, "w", encoding="utf-8", newline="\n").write(deck(c))
        print("ecrit", p)
