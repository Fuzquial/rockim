"""Pont entre la campagne etude_triax_hetero et le noyau.

Traduit chaque cas de gen_decks.py (un dictionnaire de drapeaux) en objet Essai
construit avec les préréglages du noyau. Le test d'identité compare ensuite le
deck écrit par le noyau à celui de la campagne. Si ce pont doit tricher (lire
une valeur dans gen_decks au lieu de la recalculer), c'est que le noyau ne sait
pas exprimer le cas : seuls les choix propres AU CAS viennent d'ici.
"""
import importlib.util
import os

from noyau import materiaux
from noyau.essai import (Chargement, Discontinuites, Essai, FamillePlans, JointsGrain,
                         Maillage, Sorties)

RACINE_G1 = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CAMPAGNE = os.path.join(RACINE_G1, "etude_triax_hetero")


def gen_decks():
    """Le générateur de la campagne, importé SANS l'exécuter (il n'écrit rien à l'import)."""
    spec = importlib.util.spec_from_file_location("gen_decks", os.path.join(CAMPAGNE, "gen_decks.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def essai_du_cas(G, nom_run, cas, sigma3, graine):
    c = G.CASES[cas]
    niveau = c.get("materiau", "bohus")
    traction = bool(c.get("traction"))

    jg = materiaux.joints_grain(niveau, paires=c["gbWeak"] and not c.get("gbAlpha"))
    if c.get("gbAlpha"):
        a = c["gbAlpha"]
        jg = JointsGrain(alpha_ten=a["Ten"], alpha_coh=a["Coh"], alpha_gf=a["Gf"],
                         alpha_e=a["E"], alpha_fric=a["Fric"])

    disc = Discontinuites(fraction_diffuse=c.get("prebrokenFrac") or 0.0)
    if c.get("beta") is not None:
        disc.plans = FamillePlans(pendage=c["beta"], espacement=G.SPACING, facteur=G.WP_FACTOR,
                                  fraction_rompue=c.get("brokenFrac", 0.0))

    return Essai(
        nom=nom_run, description=c["titre"],
        chargement=Chargement(type_essai="traction" if traction else "triaxial",
                              sigma3_MPa=sigma3,
                              vitesse=0.08 if traction else abs(c.get("pullV", -0.1))),
        maillage=Maillage(type="voronoi", taille_grain=G.GRAIN, taille_element=G.ELEM,
                          dispersion_tailles=c.get("spread"), graine=graine),
        materiau=materiaux.materiau(niveau),
        phases=materiaux.phases(niveau, contraste_resistance=c["phaseStrength"],
                                tailles=bool(c.get("phaseSize")), taille_grain=G.GRAIN)
        if c["phases"] else [],
        joints_grain=jg,
        discontinuites=disc,
        sorties=Sorties(T=float(c.get("T", G.T_RUN)), frames=G.FRAMES,
                        deformations_historique=False),
    )
