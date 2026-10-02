"""Essai éclair et calculs courts."""
import os
import sys

import pytest

from noyau import cfg, courts
from noyau.formulaire import choix_eclair, essai_depuis_choix

ECLAIR = os.path.join(os.path.dirname(__file__), "..", "_travail", "eclair")
FAUX = [sys.executable, os.path.join(os.path.dirname(__file__), "faux_rockim.py")]


@pytest.mark.skipif(not os.path.isdir(ECLAIR), reason="decks de l'essai éclair absents")
@pytest.mark.parametrize("chargement,nom", [("traction", "eclair_traction"), ("ucs", "eclair_ucs"), ("tx20", "eclair_tx20")])
def test_eclair_identique_aux_decks_valides_le_2026_10_02(chargement, nom):
    e = essai_depuis_choix(choix_eclair({"materiau": "fragile"}, chargement))
    ref = cfg.lire(os.path.join(ECLAIR, nom, "decks", nom + ".cfg"))
    assert cfg.differences(cfg.lire_texte(e.vers_cfg()), ref) == []


def test_eclair_garde_la_loi_editee_et_rien_d_autre():
    ch = {"materiau": "fragile", "phases": "trois", "diffus": True, "W_mm": 50,
          "surcharges": {"ft": 1.7e6}, "loi": {"yan_c": 3.0}}
    d = cfg.lire_texte(essai_depuis_choix(choix_eclair(ch, "tx20")).vers_cfg())
    assert (d["ft"], d["yanC"], d["W"], d["grainSize"]) == ("1700000", "3", "0.02", "0.01")
    assert "phases" not in d and "jointPrebrokenFrac" not in d


def test_cle_ignore_nom_et_commentaires_mais_pas_les_valeurs():
    a = essai_depuis_choix({"nom": "a"})
    b = essai_depuis_choix({"nom": "b", "description": "autre"})
    c = essai_depuis_choix({"nom": "a", "graine": 7})
    assert courts.cle(a) == courts.cle(b) != courts.cle(c)


def test_un_calcul_identique_n_est_pas_relance(tmp_path):
    e = essai_depuis_choix({"T": 1e-3})
    c1, code = courts.lancer(e, str(tmp_path), FAUX)
    assert code == 0 and courts.termine(str(tmp_path), c1)
    mtime = os.path.getmtime(os.path.join(tmp_path, c1, "run.log"))
    c2, _ = courts.lancer(e, str(tmp_path), FAUX)
    assert c2 == c1 and os.path.getmtime(os.path.join(tmp_path, c1, "run.log")) == mtime
