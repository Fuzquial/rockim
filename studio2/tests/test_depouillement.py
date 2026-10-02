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
    assert set(noyau) == set(oracle)
    faux = {k: (noyau[k], oracle[k]) for k in noyau if not egaux(noyau[k], oracle[k])}
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
