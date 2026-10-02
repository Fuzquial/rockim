"""Maillage Gmsh (régime homogène), aperçu Voronoï, et cache binaire des résultats."""
import json
import os

import numpy as np
import pytest

from campagne_triax_hetero import CAMPAGNE
from noyau import cfg, maillage, resultats
from noyau.essai import Essai
from noyau.validation import erreurs, verifier

RUN = os.path.join(CAMPAGNE, "out", "F7_disc_gbm_P020")
CACHE_J1 = os.path.join(os.path.dirname(__file__), "..", "proto", "_cache", "F7_disc_gbm_P020")


def test_gmsh_box2d(tmp_path):
    msh = str(tmp_path / "box.msh")
    s = maillage.generer_gmsh(0.036, 0.072, 0.0014, 7, msh)
    assert (s["W"], s["H"]) == (pytest.approx(0.036), pytest.approx(0.072))
    assert s["aire"] == pytest.approx(0.036 * 0.072, rel=1e-9)
    assert 2500 < s["triangles"] < 5000 and s["hauteur_min"] > 0
    # même graine, même maillage : la reproductibilité qu'on attend d'un essai
    s2 = maillage.generer_gmsh(0.036, 0.072, 0.0014, 7, str(tmp_path / "bis.msh"))
    assert open(msh).read() == open(str(tmp_path / "bis.msh")).read()


def test_essai_gmsh_ecrit_un_deck_valide(tmp_path):
    e = Essai(nom="homog_gmsh")
    e.maillage.type = "gmsh"
    e.maillage.fichier_msh = str(tmp_path / "m.msh")
    e.joints_grain.alpha_ten = e.joints_grain.alpha_coh = e.joints_grain.alpha_gf = 1.0
    assert erreurs(verifier(e)) == []
    d = cfg.lire_texte(e.vers_cfg())
    assert d["mesh"] == "file" and d["meshFile"] == e.maillage.fichier_msh
    assert "W" not in d and "grainSize" not in d       # W, H : le solveur les prend du maillage
    assert d["historyStrains"] == "true"                # défaut des nouveaux essais (spec P6)


def test_apercu_tronque_sans_toucher_l_essai():
    e = Essai(nom="x")
    a = maillage.essai_apercu(e)
    assert (a.sorties.T, a.sorties.frames, a.nom) == (2e-6, 1, "x_apercu")
    assert (e.sorties.T, e.sorties.frames) == (1e-2, 24)


def test_estimation_voronoi_sur_la_campagne():
    est = maillage.estimation_voronoi(Essai())
    assert est["elements"] == pytest.approx(12995, rel=0.02)      # mesuré : 12 995
    assert est["elements_par_grain"] == pytest.approx(35.4, rel=0.01)
    assert est["grains"] == pytest.approx(367, rel=0.05)           # mesuré : 367


@pytest.mark.skipif(not os.path.exists(os.path.join(CACHE_J1, "meta.json")), reason="cache J1 absent")
def test_cache_identique_a_celui_du_jalon_j1(tmp_path):
    m = resultats.convertir(RUN, str(tmp_path))
    j1 = json.load(open(os.path.join(CACHE_J1, "meta.json")))
    for k in ("nFrames", "nVert", "nTri", "nJoint", "nHist", "frameVersHist", "bornes", "bbox", "temps"):
        assert m[k] == j1[k], k
    for nom in ["pts", "jseg", "jmode", "phase", "hist"] + m["champs"]:
        a = open(os.path.join(tmp_path, nom + ".bin"), "rb").read()
        b = open(os.path.join(CACHE_J1, nom + ".bin"), "rb").read()
        assert a == b, nom
