"""API du serveur (J4.1) : espace jetable, faux solveur, runs réels de la campagne en lecture."""
import json
import os
import sys
import threading
import time
import urllib.error
import urllib.request

import pytest

ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ICI, "..", "app"))
import serveur as S                                   # noqa: E402
from noyau.essai import Essai                          # noqa: E402

FAUX = [sys.executable, os.path.join(ICI, "faux_rockim.py")]


@pytest.fixture(scope="module")
def api(tmp_path_factory):
    studio, srv = S.demarrer(str(tmp_path_factory.mktemp("espace")), port=0, commande=FAUX, periode=0.1)
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    base = "http://127.0.0.1:%d" % srv.server_address[1]

    def appel(chemin, corps=None, code=200):
        req = urllib.request.Request(base + chemin, method="POST" if corps is not None else "GET",
                                     data=None if corps is None else json.dumps(corps).encode())
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                assert r.status == code
                b = r.read()
                return json.loads(b) if r.headers["Content-Type"].startswith("application/json") else b
        except urllib.error.HTTPError as e:
            assert e.code == code, e.read()
            return json.loads(e.read() or b"{}")
    yield appel
    srv.shutdown()
    studio.fermer()


def attendre(appel, cond, delai=60):
    t0 = time.time()
    while time.time() - t0 < delai:
        v = appel("/api/file")
        if cond(v):
            return v
        time.sleep(0.2)
    raise TimeoutError


def test_etat_et_exes(api):
    e = api("/api/etat")
    assert "rockim_j3.exe" in e["exes"] and e["reglages"]["jobs"] == 4


def test_runs_de_la_campagne_et_synthese(api):
    runs = {r["nom"]: r for r in api("/api/runs")}
    f7 = runs["F7_disc_gbm_P020"]
    assert f7["source"] == "etude_triax_hetero"
    s = api("/api/runs/%s/synthese" % f7["id"])
    assert s["q_pic_MPa"] == pytest.approx(90.3, abs=0.05)
    assert s["prerompus_libres"] == 424
    h = api("/api/runs/%s/historique" % f7["id"])
    assert len(h["eps"]) == len(h["q"]) == 2003


def test_run_inconnu_404(api):
    assert "erreur" in api("/api/runs/b0~nexistepas/synthese", code=404)


def test_cache_a_la_demande(api):
    ident = "b0~F1_homog_P020"
    e = api("/api/runs/%s/cache" % ident, corps={})
    assert e["etat"] in ("en_cours", "pret")
    t0 = time.time()
    while api("/api/runs/%s/cache" % ident)["etat"] != "pret":
        assert time.time() - t0 < 120
        time.sleep(0.5)
    meta = api("/cache/%s/meta.json" % ident)
    assert meta["nFrames"] == 26 and meta["nTri"] == 12995
    assert len(api("/cache/%s/sigmaYY.bin" % ident)) == 26 * 12995 * 4


def test_verifier_puis_ajouter_une_variation(api):
    e = Essai(nom="var").vers_dict()
    v = api("/api/essai/verifier", corps=e)
    assert "mode = fdem" in v["deck"] and v["estimation"]["grains"] > 100
    r = api("/api/essai/ajouter", corps={"essai": e, "sigma3_liste": [0, 20], "lot": "variation σ3"})
    assert len(r["ajoutes"]) == 2
    f = attendre(api, lambda v: all(t["etat"] == "fini" for t in v["travaux"] if t["nom"].startswith("var_")))
    noms = sorted(t["nom"] for t in f["travaux"] if t["lot"] == "variation σ3")
    assert noms == ["var_P000", "var_P020"]
    assert all(t["avancement"] == 1.0 for t in f["travaux"] if t["nom"].startswith("var_"))
    tid = next(t["id"] for t in f["travaux"] if t["nom"] == "var_P000")
    assert any("12995 elements" in l for l in api("/api/file/%s/journal" % tid)["lignes"])
    assert "espace~var_P020" in {r["id"] for r in api("/api/runs")}


def test_essai_fautif_refuse(api):
    e = Essai(nom="faux")
    e.chargement.sigma3_MPa = -5
    assert "erreurs" in api("/api/essai/ajouter", corps={"essai": e.vers_dict()}, code=400)["erreur"]


def test_ordre_et_retrait(api, monkeypatch):
    for n in ("a", "b"):
        api("/api/essai/ajouter", corps={"essai": Essai(nom="ordre_" + n).vers_dict()})
    f = attendre(api, lambda v: all(t["etat"] in ("fini",) for t in v["travaux"]))
    ids = [t["id"] for t in f["travaux"]]
    f2 = api("/api/file/%s/monter" % ids[-1], corps={})
    assert f2["travaux"][-2]["id"] == ids[-1]
    f3 = api("/api/file/%s/retirer" % ids[-1], corps={})
    assert ids[-1] not in [t["id"] for t in f3["travaux"]]


def test_page_d_accueil(api):
    assert b"<!doctype html>" in api("/").lower()
