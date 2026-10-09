"""Calculs courts lancés par un clic : aperçu de maillage, essai éclair.

Ils ne passent pas par la file (ils doivent répondre en secondes) mais suivent la même
règle : un dossier par calcul, journal conservé, rien de supprimé. Le résultat est mis
en cache sous une CLÉ tirée du contenu du deck (commentaires et nom exclus) : relancer un
calcul identique rend le résultat déjà calculé, sans relancer le solveur.
"""
import hashlib
import os
import subprocess

from . import cfg


def cle(essai):
    d = cfg.lire_texte(essai.vers_cfg())
    texte = "\n".join("%s=%s" % kv for kv in sorted(d.items()))
    return hashlib.sha1(texte.encode("utf-8")).hexdigest()[:16]


def dossier(racine, c):
    return os.path.join(racine, c)


def termine(racine, c):
    """Fini avec succès : le solveur a écrit son résumé."""
    try:
        return "---- summary ----" in open(os.path.join(dossier(racine, c), "run.log"), encoding="utf-8", errors="replace").read()
    except OSError:
        return False


def lancer(essai, racine, commande, cwd=None, fils=1, delai=900):
    """Exécute le calcul (bloquant) et rend (clé, code). Ne relance pas un calcul déjà fait."""
    c = cle(essai)
    d = dossier(racine, c)
    if termine(racine, c):
        return c, 0
    os.makedirs(d, exist_ok=True)
    deck = os.path.join(d, "deck.cfg")
    with open(deck, "w", encoding="utf-8", newline="\n") as f:
        f.write(essai.vers_cfg())
    env = dict(os.environ, OMP_NUM_THREADS=str(fils))
    with open(os.path.join(d, "run.log"), "w", encoding="utf-8", errors="replace") as log:
        p = subprocess.run(list(commande) + [deck, os.path.join(d, "out")], cwd=cwd, env=env,
                           stdout=log, stderr=subprocess.STDOUT, timeout=delai)
    return c, p.returncode
