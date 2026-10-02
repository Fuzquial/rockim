"""La validation laisse passer la campagne, et chaque règle se déclenche sur son cas."""
import copy

import pytest

from campagne_triax_hetero import essai_du_cas, gen_decks
from noyau import materiaux
from noyau.essai import Essai, FamillePlans, Paire
from noyau.validation import erreurs, verifier

G = gen_decks()
RUNS = G.runs()


def codes(essai, niveau=None):
    return {a.code for a in verifier(essai) if niveau is None or a.niveau == niveau}


@pytest.mark.parametrize("lot,nom,cas,sigma3,graine", RUNS, ids=[r[1] for r in RUNS])
def test_la_campagne_passe_sans_erreur_ni_alerte(lot, nom, cas, sigma3, graine):
    e = essai_du_cas(G, nom, cas, sigma3, graine)
    assert codes(e, "erreur") == set() and codes(e, "alerte") == set()


def gbm():
    e = Essai(phases=materiaux.phases("bohus"), joints_grain=materiaux.joints_grain("bohus", paires=True))
    assert erreurs(verifier(e)) == []
    return e


CAS_FAUTIFS = {
    "gmsh_phases": lambda e: setattr(e.maillage, "type", "gmsh"),
    "fractions": lambda e: setattr(e.phases[0], "fraction", 0.5),
    "phases_doublon": lambda e: setattr(e.phases[1], "nom", e.phases[0].nom),
    "paire_inconnue": lambda e: e.joints_grain.paires.append(Paire("granite", "quartz", 1, 1, 1, 1)),
    "taille_element": lambda e: setattr(e.maillage, "taille_element", 0.004),
    "elements_par_grain": lambda e: setattr(e.maillage, "taille_element", 0.0002),
    "dispersion": lambda e: setattr(e.maillage, "dispersion_tailles", 2.0),
    "dispersion_semis": lambda e: (setattr(e.maillage, "dispersion_tailles", 0.3),
                                   setattr(e.maillage, "semis", "hex")),
    "mu_residuel": lambda e: (setattr(e.discontinuites, "fraction_diffuse", 0.05),
                              setattr(e.discontinuites, "mu_residuel", 0.58)),
    "fraction_diffuse": lambda e: setattr(e.discontinuites, "fraction_diffuse", 1.0),
    "plans": lambda e: setattr(e.discontinuites, "plans", FamillePlans(60, 0.0)),
    "segment_hors": lambda e: e.discontinuites.segments.append((0.0, 0.0, 0.05, 0.01)),
    "sigma3": lambda e: setattr(e.chargement, "sigma3_MPa", -5),
    "consolidation": lambda e: setattr(e.chargement, "delai_axial", 1e-4),
    "jauge_confinement": lambda e: setattr(e.chargement, "temps_jauge_confinement", 1e-4),
    "vitesse": lambda e: setattr(e.chargement, "vitesse", 0.0),
    "sorties": lambda e: setattr(e.sorties, "frames", 0),
}


@pytest.mark.parametrize("code", sorted(CAS_FAUTIFS))
def test_chaque_regle_se_declenche(code):
    e = gbm()
    CAS_FAUTIFS[code](e)
    assert code in codes(e)


def test_gmsh_homogene_signale_gb_inerte_et_fichier_manquant():
    e = Essai()
    e.maillage.type = "gmsh"
    c = codes(e)
    assert {"gmsh_fichier", "gmsh_gb_inerte"} <= c
    e.maillage.fichier_msh = "x.msh"
    e.joints_grain.alpha_ten = e.joints_grain.alpha_coh = e.joints_grain.alpha_gf = 1.0
    assert codes(e, "erreur") == set() and "gmsh_gb_inerte" not in codes(e)


def test_mu_residuel_par_defaut_vaut_le_frottement_de_pic():
    e = gbm()
    e.discontinuites.fraction_diffuse = 0.05
    assert "mu_residuel" not in codes(e)
    assert e.mu_residuel() == round(__import__("math").tan(__import__("math").radians(13.4)), 4)
