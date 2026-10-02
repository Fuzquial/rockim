"""Les choix de l'écran Essai suffisent à reconstruire toute la campagne etude_triax_hetero."""
import os

import pytest

from campagne_triax_hetero import CAMPAGNE, gen_decks
from noyau import cfg
from noyau.formulaire import essai_depuis_choix, estimation_cout
from noyau.validation import erreurs, verifier

G = gen_decks()
RUNS = G.runs()


def choix_du_cas(nom, cas, sigma3, graine):
    """Ce qu'un utilisateur cocherait à l'écran pour obtenir ce run."""
    c = G.CASES[cas]
    if c.get("gbAlpha"):
        jg = "inter_intra"
    elif c["gbWeak"]:
        jg = "paires"
    else:
        jg = "niveau"
    tailles = "par_phase" if c.get("phaseSize") else "dispersees" if c.get("spread") else "uniformes"
    ch = dict(nom=nom, description=c["titre"], type_essai="traction" if c.get("traction") else "triaxial",
              sigma3_MPa=sigma3, materiau=c.get("materiau", "bohus"), graine=graine,
              phases="trois" if c["phases"] else "une", contraste="complet" if c["phaseStrength"] else "elastique",
              tailles=tailles, joints_grain=jg, diffus=bool(c.get("prebrokenFrac")),
              fraction_diffuse=c.get("prebrokenFrac") or 0.05, T=float(c.get("T", G.T_RUN)),
              deformations_historique=False, champs_deformation=False)
    if c.get("gbAlpha"):
        ch["alphas"] = c["gbAlpha"]
    if c.get("pullV"):
        ch["vitesse"] = abs(c["pullV"])
    if c.get("beta") is not None:
        ch.update(plans=True, pendage=c["beta"], fraction_rompue=c.get("brokenFrac", 0.0))
    return ch


@pytest.mark.parametrize("lot,nom,cas,sigma3,graine", RUNS, ids=[r[1] for r in RUNS])
def test_le_formulaire_reconstruit_la_campagne(lot, nom, cas, sigma3, graine):
    e = essai_depuis_choix(choix_du_cas(nom, cas, sigma3, graine))
    assert erreurs(verifier(e)) == []
    disque = cfg.lire(os.path.join(CAMPAGNE, "decks", nom + ".cfg"))
    assert cfg.differences(cfg.lire_texte(e.vers_cfg()), disque) == []


def test_choix_par_defaut_valide_et_complet():
    e = essai_depuis_choix({})
    assert erreurs(verifier(e)) == []
    d = cfg.lire_texte(e.vers_cfg())
    assert d["writeStrainFields"] == "true" and d["historyStrains"] == "true"
    assert "phases" not in d and d["insertion"] == "adaptive"


def test_loi_editee_ecrit_les_cles_hors_defaut_seulement():
    d0 = cfg.lire_texte(essai_depuis_choix({}).vers_cfg())
    assert not {"yanC", "jointElastic"} & set(d0)
    d = cfg.lire_texte(essai_depuis_choix({"loi": {"yan_c": 3.0, "montee": "parabolic", "penalite": 8},
                                           "surcharges": {"ft": 1.7e6}}).vers_cfg())
    assert (d["yanC"], d["jointElastic"], d["insertionPenaltyFactor"], d["ft"]) == ("3", "parabolic", "8", "1700000")


def test_estimation_du_cout_dans_le_bon_ordre_de_grandeur():
    """Mesures du 2026-10-02 : banc éclair triaxial 16,6 s, F7 tronqué à 0,8 ms 585 s (1 fil)."""
    ecl = estimation_cout(essai_depuis_choix(dict(W_mm=20, H_mm=40, taille_grain_mm=10, taille_element_mm=2.3, T=5e-3)))
    assert 8 < ecl["duree_max_s"] < 40 and ecl["dt"] == pytest.approx(2.41e-8, rel=0.15)
    f7 = estimation_cout(essai_depuis_choix(dict(T=8e-4)))
    assert 300 < f7["duree_max_s"] < 1200 and f7["dt"] == pytest.approx(6.0e-9, rel=0.15)
