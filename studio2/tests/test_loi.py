"""La loi de l'éditeur (noyau/loi.py) est celle du solveur : comparaison à yan_point.

Référence : tests/donnees/yan_point_sn*.tsv, sorties de tools/yan_point.cpp (le programme qui
appelle YanSoftening.hpp, la fonction même du solveur), dans une copie qui accepte la
contrainte normale du chemin de cisaillement (tests/donnees/yan_point_sn.cpp, un paramètre
ajouté). N = 20 000 pas, jeu de la Table 1 de Yan 2023, pj = 4 x 15 GPa / 2 mm.
"""
import os

import numpy as np
import pytest

from noyau.loi import Loi, integrale

ICI = os.path.join(os.path.dirname(__file__), "donnees")
PJ, N = 4.0 * 15e9 / 2e-3, 20000


def sections(nom):
    out, cur, val = {}, None, {}
    for l in open(os.path.join(ICI, nom), encoding="utf-8"):
        if l.startswith("# SECTION"):
            cur = l.split()[2]
            out[cur] = []
        elif l.startswith("#"):
            p = l[1:].split()
            if len(p) == 2:
                try:
                    val[p[0]] = float(p[1])
                except ValueError:
                    pass
        elif cur:
            out[cur].append([float(x) for x in l.split()])
    return {k: np.array(v) for k, v in out.items()}, val


def loi():
    return Loi(PJ, 1.3e6, 16.4e6, 23.0, 3.8, 84.0)


def test_integrale_de_f():
    _, v = sections("yan_point_sn0.tsv")
    assert integrale() == pytest.approx(v["integral_fD"], rel=1e-12)


def test_mode_I_point_par_point():
    s, v = sections("yan_point_sn0.tsv")
    L = loi()
    assert (L.dnE, L.ot) == (pytest.approx(v["dnE"], rel=1e-12), pytest.approx(v["ot"], rel=1e-12))
    ref = s["modeI"]
    calc = np.array([L.sigma(x) for x in ref[:, 0]])
    assert np.max(np.abs(calc - ref[:, 1])) / 1.3e6 < 1e-9


@pytest.mark.parametrize("fichier,sn", [("yan_point_sn0.tsv", 0.0), ("yan_point_sn-20e6.tsv", -20e6)])
def test_mode_II_point_par_point(fichier, sn):
    s, v = sections(fichier)
    L = loi()
    sE = (16.4e6 + L.tan_phi * max(0.0, -sn)) / PJ
    smax = 1.6 * (sE + L.st)
    dtg = smax * np.arange(N + 1) / N                  # le même chemin que yan_point
    tau = L.mode2(dtg, sn)
    ref = s["modeII"]
    idx = np.rint(ref[:, 0] / smax * N).astype(int)
    assert np.max(np.abs(tau[idx] - ref[:, 1])) / 16.4e6 < 1e-9


def test_le_frottement_reste_apres_la_cohesion():
    """Sous 20 MPa de compression, tau tend vers tan(23°) x 20 MPa = 8,49 MPa, pas vers 0."""
    L = loi()
    tau = L.mode2(np.linspace(0, 3 * L.st, 4000), -20e6)
    assert tau[-1] == pytest.approx(np.tan(np.radians(23)) * 20e6, rel=1e-9)
    assert loi().mode2(np.linspace(0, 3 * L.st, 4000), 0.0)[-1] == 0.0
