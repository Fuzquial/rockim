"""Critère J2 n°1 : le noyau réécrit les 81 decks de la campagne à l'identique (sémantique).

Deux oracles :
  1. les decks présents sur disque (etude_triax_hetero/decks/*.cfg), ceux qui ont tourné ;
  2. la fonction deck() de gen_decks.py, appelée en mémoire (rien n'est écrit).
"""
import os

import pytest

from campagne_triax_hetero import CAMPAGNE, essai_du_cas, gen_decks
from noyau import cfg
from noyau.essai import Essai

G = gen_decks()
RUNS = G.runs()


def test_la_matrice_compte_81_runs():
    # decks/ en contient 86 : les 5 de plus (DIAG_gf7, DIAG_phi23, GF15, GF25, GF40)
    # sont des diagnostics écrits à la main, hors du générateur et de MATRICE.csv.
    assert len(RUNS) == 81


@pytest.mark.parametrize("lot,nom,cas,sigma3,graine", RUNS, ids=[r[1] for r in RUNS])
def test_deck_identique(lot, nom, cas, sigma3, graine):
    essai = essai_du_cas(G, nom, cas, sigma3, graine)
    noyau = cfg.lire_texte(essai.vers_cfg())
    disque = cfg.lire(os.path.join(CAMPAGNE, "decks", nom + ".cfg"))
    memoire = cfg.lire_texte(G.deck(G.CASES[cas], sigma3, seed=graine))
    assert cfg.differences(noyau, disque) == []
    assert cfg.differences(noyau, memoire) == []


@pytest.mark.parametrize("lot,nom,cas,sigma3,graine", RUNS[::9], ids=[r[1] for r in RUNS[::9]])
def test_aller_retour_json(lot, nom, cas, sigma3, graine):
    """Un essai sérialisé puis relu redonne exactement le même deck."""
    e = essai_du_cas(G, nom, cas, sigma3, graine)
    assert Essai.depuis_dict(e.vers_dict()).vers_cfg() == e.vers_cfg()


def test_le_comparateur_detecte_un_ecart():
    """Contrôle falsifiant : une valeur modifiée doit faire échouer la comparaison."""
    e = essai_du_cas(G, *RUNS[0][1:])
    d = cfg.lire_texte(e.vers_cfg())
    e.materiau.E *= 1.0001
    assert cfg.differences(cfg.lire_texte(e.vers_cfg()), d) == ["E : %s / %s" % (cfg.nombre(e.materiau.E), d["E"])]
