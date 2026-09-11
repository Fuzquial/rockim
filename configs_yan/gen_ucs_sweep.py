#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# gen_ucs_sweep.py — compression UNIAXIALE 2D FEMDEM, eprouvette HOMOGENE,
# balayage des parametres de la loi cohesive.
#
#   python configs_yan/gen_ucs_sweep.py      # ecrit configs_yan/ucs_sw_*.cfg
#
# QUESTION POSEE : quel parametre gouverne la REPONSE de l'eprouvette —
# raideur, pic, et surtout FRAGILITE (chute post-pic) ?
#
# CE QUI EST REPRIS : le banc de Yan, Zheng & Wang (IJRMMS 169, 2023) tel qu'il
# est deja dans configs_yan/ucs_adap.cfg et tx_adap_*.cfg — eprouvette
# 36 x 72 mm, plateaux rigides frottants a 0,1 m/s, insertion ADAPTATIVE des
# cohesifs + loi de Yan (jointSoftening). Confinement NUL : uniaxial.
#
# VRAIMENT HOMOGENE : `mesh = file`, un maillage Delaunay/Netgen importe.
# AUCUNE microstructure — ni phases, ni grains, ni frontieres. Correction du
# 2026-09-11 : une premiere version utilisait `mesh = voronoi` a une seule
# phase en supposant que le Voronoi n'etait qu'un support de mailleur. C'est
# faux : meme a materiau identique, la tessellation cree un RESEAU DE JOINTS
# ALIGNE sur les frontieres de cellules, donc des chemins de fissuration
# privilegies — et le facies conjugue observe pouvait en partie en venir.
# La doc est explicite : « la regle de la these impose mesh = file pour tout
# essai lu en FACIES ou en energie de bande ».
#
# MAILLAGE NON STRUCTURE : c'est aussi la regle 8.4 — « maillage structure en
# FDEM = CONDITION D'INVALIDITE (trajets biaises, divergence en phase debris) ».
# `mesh = grid` (le defaut) est structure et donc exclu ; le maillage importe
# est un Delaunay optimise Netgen, sans direction privilegiee.
#
# POURQUOI CE BALAYAGE-LA. Mesure du 2026-09-11 sur la serie triaxiale GBM :
# au pic, 96-99 % des joints inseres sont a D < 0,2 — ils portent de la charge
# et ne cassent pas. Cause chiffree : la LONGUEUR COHESIVE de mode II,
#     l_ch,II = E G_fII / c^2 ,
# valait 173 mm pour une eprouvette de 36 x 72 mm. La zone de process ne tient
# pas dans l'eprouvette, donc aucune bande de cisaillement ne peut se former et
# la courbe plafonne au lieu de chuter. Origine de l'erreur : gfShearFactor =
# 22,105 a ete herite du deck de Yan — ou G_fI = 3,8 J/m2, ce qui donne
# G_fII = 84 J/m2 et l_ch,II = 4,7 mm, PLUS PETIT QUE SON GRAIN — puis applique
# tel quel a un G_fI de granite de 70 J/m2, soit G_fII = 1547 J/m2.
#
# En compression la rupture est en CISAILLEMENT (Yan : bande conjuguee a 62
# deg, « cisaillement dominant »), donc c'est l_ch,II qui gouverne, pas l_ch,I.
# Le balayage teste cette lecture : si elle est juste, la fragilite doit suivre
# l_ch,II et NON G_fI seul.
#
# LES HUIT CAS (un seul parametre change a la fois, sauf le dernier)
#   ref      : la carte granite telle qu'utilisee jusqu'ici          l_ch,II = 173 mm
#   gfs10    : gfShearFactor 22,105 -> 10  (borne haute de la note)   78 mm
#   gfs5     : gfShearFactor -> 5          (borne basse de la note)   39 mm
#   gfs1.3   : gfShearFactor -> 1,3        (l_ch,II ~ taille de grain) 10 mm
#   gfI10    : G_fI 70 -> 10 J/m2, gfShearFactor inchange             25 mm
#   ft20     : f_t 10 -> 20 MPa  (l_ch,I /4 ; l_ch,II inchange)      173 mm
#   coh40    : cohesion 25 -> 40 MPa                                  68 mm
#   brittle  : G_fI = 10 ET gfShearFactor = 5                        5,6 mm
#
# LECTURE : la raideur initiale ne doit dependre d'AUCUN de ces parametres (ils
# ne touchent ni E ni nu) ; le pic doit suivre f_t et la cohesion ; la chute
# post-pic doit suivre l_ch,II. Tout ecart a cette lecture est un resultat.
# ---------------------------------------------------------------------------
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# carte de reference : granite, moyenne des trois mineraux de
# configs/fdem_percussion_gbm_fin.cfg (l'eprouvette est HOMOGENE ici)
REF = dict(E=70e9, nu=0.25, rho=2600, ft=10e6, cohesion=25e6,
           frictionDeg=40, Gf=70.0, gfShearFactor=22.105)

CASES = [
    ("ref",     {}),
    ("gfs10",   dict(gfShearFactor=10.0)),
    ("gfs5",    dict(gfShearFactor=5.0)),
    ("gfs1p3",  dict(gfShearFactor=1.3)),
    ("gfI10",   dict(Gf=10.0)),
    ("ft20",    dict(ft=20e6)),
    ("coh40",   dict(cohesion=40e6)),
    ("brittle", dict(Gf=10.0, gfShearFactor=5.0)),
]

HEADER = """# ---------------------------------------------------------------------------
# UCS 2D FEMDEM, eprouvette HOMOGENE, maillage NON STRUCTURE — cas « {tag} »
# Genere par configs_yan/gen_ucs_sweep.py (y lire le detail du balayage).
#
# Banc    : Yan, Zheng & Wang (IJRMMS 169, 2023) — 36 x 72 mm, plateaux
#           rigides frottants a 0,1 m/s, confinement NUL.
# Schema  : insertion ADAPTATIVE + loi de Yan (jointSoftening).
# Maillage: NON STRUCTURE importe (Gmsh Delaunay + Netgen), AUCUNE
#           microstructure — ni phases, ni grains, ni frontieres (regle 8.4 :
#           un maillage structure est une condition d'invalidite ; et la regle
#           de la these impose mesh = file pour tout essai lu en facies).
#
# Ce cas : {desc}
#   l_ch,I  = E G_fI  / f_t^2 = {lchI:6.1f} mm
#   l_ch,II = E G_fII / c^2   = {lchII:6.1f} mm      (eprouvette : 36 x 72 mm)
#
#   rockim_g1fix.exe configs_yan/ucs_sw_{tag}.cfg out_ucs_sw_{tag}
# ---------------------------------------------------------------------------
"""


def deck(tag, over):
    m = dict(REF)
    m.update(over)
    GfII = m["Gf"] * m["gfShearFactor"]
    lchI = m["E"] * m["Gf"] / m["ft"] ** 2 * 1e3
    lchII = m["E"] * GfII / m["cohesion"] ** 2 * 1e3
    desc = "reference" if not over else ", ".join(
        "%s = %g" % (k, v) for k, v in over.items())
    L = [HEADER.format(tag=tag, desc=desc, lchI=lchI, lchII=lchII)]
    a = L.append
    a("mode = fdem")
    a("scenario = tension")
    a("loading = platens")
    # T couvre le pic (~1,7 ms a 0,23 % et 1,39 /s) plus le post-pic ; l'arret
    # automatique coupe avant si la chute est franche (voir plus bas).
    a("T = 3.5e-3")
    a("frames = 20")
    a("historyFlush = true")
    a("")
    a("# --- eprouvette -----------------------------------------------------")
    a("W = 0.036")
    a("H = 0.072")
    a("thickness = 1.0")
    a("")
    a("# --- maillage : NON STRUCTURE IMPORTE, aucune microstructure --------")
    a("# `mesh = voronoi` a UNE phase n'est PAS un milieu homogene : la")
    a("# tessellation cree un reseau de joints ALIGNE sur les frontieres de")
    a("# cellules, chemins de fissuration privilegies meme a materiau")
    a("# identique. DOCUMENTATION_rockim.md : « la regle de la these impose")
    a("# mesh = file pour tout essai lu en FACIES ou en energie de bande ».")
    a("# Maillage genere par :")
    a("#   python tools/make_unstructured_mesh.py box2d 0.036 0.072 0.0015 \\")
    a("#          meshes/ucs2d_h15.msh 4211")
    a("# Delaunay + optimisation Netgen (Gmsh), 3192 triangles, 4796 noeuds.")
    a("mesh = file")
    a("meshFile = meshes/ucs2d_h15.msh")
    a("")
    a("# --- materiau HOMOGENE (une seule phase) ----------------------------")
    a("rho = %g" % m["rho"])
    a("E = %g" % m["E"])
    a("nu = %g" % m["nu"])
    a("ft = %g" % m["ft"])
    a("cohesion = %g" % m["cohesion"])
    a("frictionDeg = %g" % m["frictionDeg"])
    a("Gf = %g" % m["Gf"])
    a("gfShearFactor = %g               # G_fII = %.0f J/m2" % (m["gfShearFactor"], GfII))
    a("")
    a("# --- schema ---------------------------------------------------------")
    a("insertion = adaptive")
    a("jointSoftening = yan")
    a("insertionPenaltyFactor = 4")
    a("contactMu = 0.1")
    a("dampingLocal = 0.7               # quasi-statique")
    a("jointXi = 0.0")
    a("")
    a("# --- chargement (UNIAXIAL : aucun confinement) ----------------------")
    a("pullV = -0.1")
    a("pullRamp = 2e-4")
    a("gaugeLoFrac = 0.25")
    a("gaugeHiFrac = 0.75")
    a("ucsStopAfterPeak = true")
    a("ucsStopDelay = 2.0e-4")
    # Le defaut attend une chute de 70 % : sur les cas ductiles il n'arrive
    # jamais et le run va jusqu'a T. 35 % suffit a caracteriser la chute et
    # libere la machine (mesure du 11/09 : les triaxiaux n'avaient pas verrouille).
    a("stopPeakDrop = 0.35")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    print(" %-9s %10s %10s %10s %10s" % ("cas", "G_fI", "G_fII", "l_ch,I", "l_ch,II"))
    for tag, over in CASES:
        p = os.path.join(HERE, "ucs_sw_%s.cfg" % tag)
        io.open(p, "w", encoding="utf-8", newline="\n").write(deck(tag, over))
        m = dict(REF)
        m.update(over)
        GfII = m["Gf"] * m["gfShearFactor"]
        print(" %-9s %10.1f %10.0f %8.1f mm %8.1f mm"
              % (tag, m["Gf"], GfII,
                 m["E"] * m["Gf"] / m["ft"] ** 2 * 1e3,
                 m["E"] * GfII / m["cohesion"] ** 2 * 1e3))
