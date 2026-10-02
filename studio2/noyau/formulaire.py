"""Choix de l'écran Essai -> Essai complet.

L'utilisateur choisit en termes physiques (« trois minéraux », « tailles par phase »,
« joints inter/intra »...) ; ce module en déduit l'essai, avec les préréglages du noyau.
C'est la seule traduction : le serveur et les tests passent par elle. Le test
tests/test_formulaire.py reconstruit les 81 runs de etude_triax_hetero à partir de
leurs seuls choix et retrouve les decks de la campagne.

Un choix est un dictionnaire JSON ; toute clé absente prend la valeur de CHOIX_DEFAUT.
"""
import copy
import math

from . import materiaux
from .essai import (Chargement, Discontinuites, Eprouvette, Essai, FamillePlans, JointsGrain,
                    Maillage, Schema, Sorties)

CHOIX_DEFAUT = {
    "nom": "essai",
    "description": "",
    "type_essai": "triaxial",            # triaxial | traction
    "sigma3_MPa": 20.0,
    "vitesse": None,                     # m/s ; None : 0,1 en triaxial, 0,08 en traction
    "delai_axial": 6e-4,                 # s, fin de consolidation avant la charge axiale
    "chute_arret": 0.5,                  # arrêt quand sigma < (1 - x) pic
    "W_mm": 36.0, "H_mm": 72.0,
    "materiau": "fragile",               # bohus | fragile
    "surcharges": {},                    # propriétés de volume/joint modifiées (E, ft, Gf...)
    "maillage": "voronoi",               # voronoi | gmsh
    "taille_grain_mm": 3.0,
    "taille_element_mm": 0.7,
    "graine": 4211,
    "fichier_msh": None,
    "phases": "une",                     # une | trois
    "contraste": "complet",              # elastique | complet
    "tailles": "uniformes",              # uniformes | dispersees | par_phase
    "dispersion": 0.35,
    "joints_grain": "niveau",            # niveau | identiques | paires | inter_intra
    "alphas": {"Ten": 0.67, "Coh": 0.73, "Gf": 0.55, "E": 1.30, "Fric": 1.0},
    "diffus": False, "fraction_diffuse": 0.05,
    "plans": False, "pendage": 60.0, "espacement_mm": 9.0, "facteur_plans": 0.3, "fraction_rompue": 0.0,
    "segments": [],                      # [[x1, y1, x2, y2], ...] en m, dessinés sur l'aperçu
    "loi": {},                           # montee, adoucissement, yan_c, penalite (éditeur de loi)
    "T": None,                           # None : durée par défaut du type d'essai
    "frames": 24,
    "deformations_historique": True,
    "champs_deformation": True,
}

T_DEFAUT = {"triaxial": 1.0e-2, "traction": 1.2e-3}


def complet(choix):
    c = copy.deepcopy(CHOIX_DEFAUT)
    c.update(copy.deepcopy(choix or {}))
    return c


def essai_depuis_choix(choix):
    c = complet(choix)
    niveau = c["materiau"]
    traction = c["type_essai"] == "traction"
    mat = materiaux.materiau(niveau)
    for k, v in c["surcharges"].items():
        setattr(mat, k, v)

    phases = []
    if c["phases"] == "trois" and c["maillage"] == "voronoi":
        phases = materiaux.phases(niveau, contraste_resistance=c["contraste"] == "complet",
                                  tailles=c["tailles"] == "par_phase",
                                  taille_grain=c["taille_grain_mm"] / 1000)

    jg = materiaux.joints_grain(niveau, paires=c["joints_grain"] == "paires" and bool(phases))
    if c["joints_grain"] == "identiques":
        jg = JointsGrain(1.0, 1.0, 1.0, 1.0, 1.0)
    elif c["joints_grain"] == "inter_intra":
        a = c["alphas"]
        jg = JointsGrain(a["Ten"], a["Coh"], a["Gf"], a["E"], a["Fric"])

    disc = Discontinuites(fraction_diffuse=c["fraction_diffuse"] if c["diffus"] else 0.0,
                          segments=[tuple(s) for s in c["segments"]])
    if c["plans"]:
        disc.plans = FamillePlans(pendage=c["pendage"], espacement=c["espacement_mm"] / 1000,
                                  facteur=c["facteur_plans"], fraction_rompue=c["fraction_rompue"])

    loi = c["loi"]
    schema = Schema(montee=loi.get("montee", "linear"), adoucissement=loi.get("adoucissement", "yan"),
                    yan_c=loi.get("yan_c", 6.0), penalite_insertion=loi.get("penalite", 4.0))

    return Essai(
        nom=c["nom"], description=c["description"],
        eprouvette=Eprouvette(W=c["W_mm"] / 1000, H=c["H_mm"] / 1000),
        chargement=Chargement(type_essai=c["type_essai"], sigma3_MPa=0.0 if traction else c["sigma3_MPa"],
                              vitesse=c["vitesse"] or (0.08 if traction else 0.1),
                              chute_arret=c["chute_arret"], delai_axial=c["delai_axial"]),
        maillage=Maillage(type=c["maillage"], taille_grain=c["taille_grain_mm"] / 1000,
                          taille_element=c["taille_element_mm"] / 1000, graine=int(c["graine"]),
                          dispersion_tailles=c["dispersion"] if c["tailles"] in ("dispersees", "par_phase") else None,
                          fichier_msh=c["fichier_msh"]),
        materiau=mat, phases=phases, joints_grain=jg, discontinuites=disc, schema=schema,
        sorties=Sorties(T=c["T"] or T_DEFAUT[c["type_essai"]], frames=int(c["frames"]),
                        deformations_historique=c["deformations_historique"],
                        champs_deformation=c["champs_deformation"]),
    )


# ---------------------------------------------------------------- estimation du coût
# Mesures du 2026-10-02 sous rockim_j3, un fil : F7 (12 995 éléments, h = 0,7 mm) 4,4 ms par
# pas, dt = 6,0e-9 s ; banc éclair (367 éléments, h = 2,3 mm) dt = 2,4e-8 s, 80 µs par pas en
# triaxial, 43 en compression simple, 29 en traction. D'où dt ~ 0,028 h / c (c = sqrt(E/rho))
# et un coût par pas de 0,2 à 0,35 µs par élément. À 4 fils le 2D ne gagne que ~30 %
# (campagne : 3,2 ms par pas pour F7). C'est un ordre de grandeur, pas une promesse.
def estimation_cout(essai, fils=1):
    m, e, mat = essai.maillage, essai.eprouvette, essai.materiau
    from .geometrie import nombre_elements
    n_el = nombre_elements(e.W, e.H, m.taille_element)
    dt = 0.028 * m.taille_element / math.sqrt(mat.E / mat.rho)
    pas = essai.sorties.T / dt
    par_el = 0.1e-6 if essai.chargement.type_essai == "traction" else 0.3e-6
    s = pas * n_el * par_el * (1.0 if fils == 1 else 0.73)
    return {"elements": int(round(n_el)), "dt": dt, "pas_max": int(pas), "duree_max_s": s}
