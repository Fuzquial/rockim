"""Critère J2 n°2 : le dépouillement du noyau rend les chiffres de depouille.py.

Oracle : les fonctions mesures / joints / diagnostics de la campagne, appelées en
mémoire (rien n'est écrit), sur chaque run de etude_triax_hetero/out qui a un
history.csv. Une seule différence voulue : le noyau ajoute `grains` aux
diagnostics ; elle est retirée avant comparaison.
"""
import importlib.util
import math
import os

import pytest

from campagne_triax_hetero import CAMPAGNE
from noyau import depouillement as D
from noyau import resultats as R

spec = importlib.util.spec_from_file_location("depouille", os.path.join(CAMPAGNE, "depouille.py"))
ORACLE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ORACLE)

OUT = os.path.join(CAMPAGNE, "out")
RUNS = sorted(d for d in os.listdir(OUT) if os.path.exists(os.path.join(OUT, d, "history.csv")))


def egaux(a, b):
    if isinstance(a, float) or isinstance(b, float):
        if math.isnan(a) and math.isnan(b):
            return True
        return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)
    return a == b


def comparer(noyau, oracle):
    if oracle is None:                   # run arrêté avant sa première ligne d'historique
        assert noyau is None
        return
    # Le noyau peut rendre PLUS que l'oracle (pic filtré), jamais moins ni autre chose.
    assert set(oracle) <= set(noyau)
    faux = {k: (noyau[k], oracle[k]) for k in oracle if not egaux(noyau[k], oracle[k])}
    assert faux == {}


def test_il_y_a_des_runs_a_comparer():
    assert len(RUNS) >= 15


@pytest.mark.parametrize("run", RUNS)
def test_mesures(run):
    d = os.path.join(OUT, run)
    comparer(D.mesures(R.lire_historique(d), R.delai_consolidation(d)),
             ORACLE.mesures(ORACLE.lire_history(os.path.join(d, "history.csv")),
                            ORACLE.delai_consolidation(d)))


@pytest.mark.parametrize("run", RUNS)
def test_joints(run):
    p = os.path.join(OUT, run, "fdem_final_joints.csv")
    comparer(D.joints(p), ORACLE.joints(p))


@pytest.mark.parametrize("run", RUNS)
def test_diagnostics(run):
    p = os.path.join(CAMPAGNE, "logs", run + ".log")
    n = D.diagnostics(p)
    n.pop("grains", None)
    comparer(n, ORACLE.diagnostics(p))


def test_historique_en_cours_d_ecriture(tmp_path):
    """Une dernière ligne coupée (run en cours) est écartée, pas fatale."""
    src = open(os.path.join(OUT, RUNS[0], "history.csv"), encoding="utf-8").read().splitlines()
    (tmp_path / "history.csv").write_text("\n".join(src[:50]) + "\n" + src[50][:7], encoding="utf-8")
    h = R.lire_historique(str(tmp_path))
    assert len(h["t"]) == 49
    assert R.dernier_temps(str(tmp_path)) == h["t"][-1]


ECLAIR = os.path.join(os.path.dirname(__file__), "..", "_travail", "eclair")


@pytest.mark.skipif(not os.path.isdir(ECLAIR), reason="runs de l'essai éclair absents")
def test_pic_filtre_ecarte_les_pointes_du_banc_eclair():
    """Mesure du 2026-10-02 : maximum brut 74,8 MPa, enveloppe 48,3 MPa en UCS."""
    d = os.path.join(ECLAIR, "eclair_ucs", "out", "eclair_ucs")
    m = D.mesures(R.lire_historique(d), R.delai_consolidation(d))
    assert m["q_pic_MPa"] == pytest.approx(74.8, abs=0.1)
    assert m["q_pic_filtre_MPa"] == pytest.approx(48.3, abs=0.2)


@pytest.mark.skipif(not os.path.isdir(ECLAIR), reason="runs de l'essai éclair absents")
def test_traction_par_mors_lit_epsAx():
    d = os.path.join(ECLAIR, "eclair_traction", "out", "eclair_traction")
    m = D.mesures(R.lire_historique(d), R.delai_consolidation(d))
    assert m["q_pic_MPa"] == pytest.approx(1.379, abs=1e-3)
    assert 0 < m["eps_pic_pct"] < 0.02                  # epsAx lu, plus de NaN


def test_pic_filtre_sur_courbe_lisse_egale_le_brut():
    import numpy as np
    t = np.linspace(0, 1, 500)
    h = {"t": t, "sigma": -np.sin(np.pi * t) * 1e8, "epsGauge": t * 1e-2}
    m = D.mesures(h)
    # Sur un pic lisse la médiane de ±7 points rend la 8e valeur sur 15 : sin(pi (0,5 ± 0,007))
    # = 1 - 2,4e-4. Le biais est de cet ordre, négligeable devant les pointes qu'on écarte.
    assert m["q_pic_filtre_MPa"] == pytest.approx(m["q_pic_MPa"], rel=1e-3)
    assert m["q_pic_filtre_MPa"] <= m["q_pic_MPa"]
