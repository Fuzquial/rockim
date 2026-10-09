"""Dépouillement du brésilien : historique et verdicts du solveur.

Le journal synthétique reprend mot pour mot les formats de FdemSolver.cpp (résumé
« brazilian », l. 10239-10316) ; l'historique a les colonnes de l'en-tête l. 9634.
"""
import numpy as np
import pytest

from noyau import depouillement as D
from noyau import resultats as R

JOURNAL = """[FDEM] brazilian: two rigid platens closing at 0.1 m/s total
[rockim] solver requested an early stop at t = 0.000583 s (96589 / 180471 steps)
[FDEM] ---- summary ----
[FDEM] brazilian ELASTIC-BAND gauge (sigma_t in [0.3, 0.8] x ft, disc intact, 214 history steps, 0.4 % of joints damaged at the top of the band):
[FDEM]   mean sigma_t = 0.71 MPa, mean centre sigma_xx = 0.7107 MPa -> ratio 1.001  [PASS]  (band 0.85-1.25; mean sigma_yy = -1.62 MPa, sigma_yy/sigma_xx = -2.28)
[FDEM] brazilian elastic gauge (last step before first breakage, sigma_t = 2.39 MPa):
[FDEM]   centre sigma_xx = 0.54 MPa, expected +2.39 MPa -> ratio 0.227  [FAIL]  (this IS what sigma_t reports, band 0.85-1.25)
[FDEM] brazilian (D = 0.0508 m, t = 1 m, bearing width 0.005 m):
[FDEM]   peak force P = 115600 N at t = 0.00051 s (431 joints already broken; peak LOCKED at the post-failure load drop)
[FDEM]   indirect tensile strength sigma_t = 2P/(pi D t) = 1.449 MPa
[FDEM]   ratio to the bulk ft (1.3 MPa) = 1.115
[FDEM]   crack location: 812 / 1004 broken joints within 15 % of R from the load axis (80.9 % diametral), mean |x-xc|/R = 0.07
"""


def test_verdicts_du_solveur(tmp_path):
    p = tmp_path / "run.log"
    p.write_text(JOURNAL, encoding="utf-8")
    d = D.diagnostics_bresilien(str(p))
    # la jauge de BANDE (1,001 PASS), pas celle du dernier pas intact (0,227 FAIL)
    assert (d["jauge_elastique"], d["jauge_elastique_verdict"], d["syy_sur_sxx"]) == (1.001, "PASS", -2.28)
    assert (d["sigma_t_solveur_MPa"], d["sigma_t_sur_ft"], d["diametral_pct"]) == (1.449, 1.115, 80.9)
    assert d["P_pic_solveur_N"] == 115600 and d["pic_verrouille_solveur"]
    assert (d["pas_effectues"], d["pas_plafond"]) == (96589, 180471)


def test_jauge_non_mesuree(tmp_path):
    p = tmp_path / "run.log"
    p.write_text("[FDEM] brazilian ELASTIC-BAND gauge: NOT measured — the disc broke before", encoding="utf-8")
    assert D.diagnostics_bresilien(str(p))["jauge_elastique_verdict"] == "non mesurée"


def historique():
    t = np.linspace(0, 6e-4, 400)
    s = 1.4e6 * np.sin(np.pi * t / 1.2e-3)
    s[150] = 3.0e6                                    # une pointe parasite
    return {"t": t, "P": s * np.pi * 0.0508 / 2, "Pbot": s, "drive": 0.0508 - 0.1 * t / 2,
            "sigmaT": s, "sigmaTpeak": np.maximum.accumulate(s), "nBroken": np.arange(400.0),
            "nFrag": np.ones(400), "sxxC": s, "syyC": -3 * s, "peakLocked": (t > 5.9e-4).astype(float)}


def test_mesures_et_courbe_du_bresilien():
    h = historique()
    m = D.mesures_bresilien(h)
    assert m["sigma_t_pic_MPa"] == pytest.approx(3.0)            # le maximum brut prend la pointe
    assert m["sigma_t_pic_filtre_MPa"] == pytest.approx(1.4, rel=1e-3)
    assert m["pic_verrouille"] == 1 and m["essai"] == "bresilien"
    x, y = R.courbe(h)
    assert x[0] == 0 and x[-1] == pytest.approx(1e3 * 0.1 * 6e-4 / 2) and y.max() == pytest.approx(3.0)
    assert R.etiquettes(h)[1].startswith("σt")
