"""File de calculs, testée avec un faux solveur (aucun lancement de g1)."""
import json
import os
import subprocess
import sys
import time

import pytest

from noyau.essai import Essai
from noyau.file import File, pid_vivant

FAUX = [sys.executable, os.path.join(os.path.dirname(__file__), "faux_rockim.py")]


def essai(nom):
    return Essai(nom=nom)


def intervalles(f):
    return [json.load(open(os.path.join(t["out"], "faux.json"))) for t in f.travaux
            if os.path.exists(os.path.join(t["out"], "faux.json"))]


def test_parallelisme_borne_et_ordre(tmp_path):
    f = File(tmp_path, FAUX, jobs=2, fils=1)
    for i in range(5):
        f.ajouter(essai("r%d" % i))
    pic = 0
    t0 = time.time()
    while any(t["etat"] in ("en_attente", "en_cours") for t in f.travaux):
        f.pas()
        pic = max(pic, len(f.en_cours()))
        assert time.time() - t0 < 60
        time.sleep(0.05)
    assert pic == 2
    assert [t["etat"] for t in f.travaux] == ["fini"] * 5
    debuts = [json.load(open(os.path.join(t["out"], "faux.json")))["debut"] for t in f.travaux]
    assert debuts == sorted(debuts)                     # lancés dans l'ordre de la file
    # jamais plus de deux runs en même temps, d'après leurs propres horloges
    iv = intervalles(f)
    for x in iv:
        chevauchent = sum(1 for y in iv if y["debut"] < x["fin"] and x["debut"] < y["fin"])
        assert chevauchent <= 2


def test_echec_et_journal(tmp_path, monkeypatch):
    monkeypatch.setenv("FAUX_CODE", "3")
    f = File(tmp_path, FAUX, jobs=1)
    t = f.ajouter(essai("casse"))
    f.attendre(delai=30)
    assert t["etat"] == "echec" and t["code"] == 3
    assert "12995 elements" in open(t["log"], encoding="utf-8").read()


def test_deja_fait_n_est_pas_relance_et_relancer_renomme(tmp_path):
    f = File(tmp_path, FAUX, jobs=1)
    t = f.ajouter(essai("a"))
    f.attendre(delai=30)
    assert t["etat"] == "fini"
    t2 = f.ajouter(essai("a"))                          # même nom, dossier déjà rempli
    f.attendre(delai=30)
    assert t2["etat"] == "deja_fait"
    f.relancer(t2["id"])
    f.attendre(delai=30)
    assert t2["etat"] == "fini"
    anciens = [d for d in os.listdir(os.path.join(tmp_path, "out")) if d.startswith("a__ancien_")]
    assert len(anciens) == 1                            # l'ancien résultat est conservé


def test_arret_d_un_run(tmp_path, monkeypatch):
    monkeypatch.setenv("FAUX_DUREE", "20")
    f = File(tmp_path, FAUX, jobs=1)
    t = f.ajouter(essai("long"))
    f.pas()
    pid = t["pid"]
    assert pid_vivant(pid)
    f.arreter(t["id"])
    assert t["etat"] == "arrete" and not pid_vivant(pid)


def test_avancement(tmp_path, monkeypatch):
    monkeypatch.setenv("FAUX_DUREE", "3")
    f = File(tmp_path, FAUX, jobs=1)
    t = f.ajouter(essai("suivi"))
    f.pas()
    vus = []
    while t["etat"] == "en_cours":
        vus.append(f.avancement(t["id"]))
        time.sleep(0.25)
        f.pas()
    assert vus == sorted(vus) and 0 < max(vus) < 1 and f.avancement(t["id"]) == 1.0


def test_reouverture_rattache_le_run_vivant(tmp_path, monkeypatch):
    """L'interface se ferme pendant un calcul : à la réouverture, le run est
    reconnu vivant, puis son état final est déduit de son journal."""
    monkeypatch.setenv("FAUX_DUREE", "2")
    f = File(tmp_path, FAUX, jobs=1)
    t = f.ajouter(essai("persistant"))
    f.pas()
    f._fermer_journal(t["id"])                          # simule la fermeture : on lâche le Popen
    del f
    g = File(tmp_path, FAUX, jobs=1)
    t2 = g.trouver(t["id"])
    assert t2["etat"] == "en_cours" and pid_vivant(t2["pid"])
    g.attendre(delai=30)
    assert t2["etat"] == "fini"


def test_reouverture_apres_mort_du_run(tmp_path):
    f = File(tmp_path, FAUX, jobs=1)
    t = f.ajouter(essai("mort"))
    t.update(etat="en_cours", pid=999999)               # PID qui n'existe pas, aucun journal
    f._sauver()
    g = File(tmp_path, FAUX, jobs=1)
    assert g.trouver(t["id"])["etat"] == "interrompu"


def test_retirer_et_deplacer(tmp_path):
    f = File(tmp_path, FAUX, jobs=1)
    a, b, c = (f.ajouter(essai(n)) for n in "abc")
    f.deplacer(c["id"], -2)
    assert [t["nom"] for t in f.travaux] == ["c", "a", "b"]
    f.retirer(a["id"])
    assert [t["nom"] for t in File(tmp_path, FAUX).travaux] == ["c", "b"]
    with pytest.raises(ValueError):
        f.ajouter(essai("b"))                            # doublon en attente refusé
