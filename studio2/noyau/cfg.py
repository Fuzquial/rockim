"""Lecture et écriture des fichiers .cfg, avec la sémantique exacte du solveur.

Config::load (src/Config.cpp:19-37) : tout ce qui suit `#` est ignoré, la ligne
est coupée au premier `=`, clé et valeur sont rognées, et la DERNIÈRE occurrence
d'une clé l'emporte. Ce module reproduit ces règles et rien de plus.

Les nombres sont écrits avec un point décimal quelle que soit la locale Windows
(repr de Python ne dépend pas de la locale), au plus court sans perte.
"""
import math


def lire(chemin):
    """Dictionnaire clé -> valeur (chaînes), comme le solveur le voit."""
    with open(chemin, encoding="utf-8", errors="replace") as f:
        return lire_texte(f.read())


def lire_texte(texte):
    d = {}
    for ligne in texte.splitlines():
        ligne = ligne.split("#", 1)[0].strip()
        if "=" not in ligne:
            continue
        cle, val = ligne.split("=", 1)
        if cle.strip():
            d[cle.strip()] = val.strip()
    return d


def nombre(x):
    """Écriture d'une valeur numérique ou booléenne, indépendante de la locale."""
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        if not math.isfinite(x):
            raise ValueError("valeur non finie : %r" % x)
        if x == int(x) and abs(x) < 1e15:
            return str(int(x))
        return repr(x)
    return str(x)


def _jetons(v):
    """Une valeur vue comme suite de nombres (None si l'un des jetons n'en est pas)."""
    jetons = v.replace(";", " ").split()
    try:
        return [float(j) for j in jetons]
    except ValueError:
        return None


def valeurs_egales(a, b):
    """Deux valeurs sont égales si le solveur les lit pareil : numériquement
    quand ce sont des nombres (« 20e6 » = « 20000000 »), en texte sinon."""
    if a == b:
        return True
    na, nb = _jetons(a), _jetons(b)
    if na is not None and nb is not None:
        return len(na) == len(nb) and all(x == y for x, y in zip(na, nb))
    return a.lower() == b.lower() and a.lower() in ("true", "false")


def differences(d1, d2):
    """Liste lisible des écarts sémantiques entre deux decks lus."""
    out = []
    for k in sorted(set(d1) | set(d2)):
        if k not in d1:
            out.append("%s : absent / %s" % (k, d2[k]))
        elif k not in d2:
            out.append("%s : %s / absent" % (k, d1[k]))
        elif not valeurs_egales(d1[k], d2[k]):
            out.append("%s : %s / %s" % (k, d1[k], d2[k]))
    return out


class Ecrivain:
    """Construit un deck lisible : sections commentées, une clé par ligne.

    Une clé écrite deux fois est une erreur : chez le solveur la seconde
    écraserait la première en silence (piège de preBrokenJoints, DOC §5.11).
    """

    def __init__(self, titre_lignes=()):
        self.lignes = ["# " + "-" * 73]
        self.lignes += ["# " + t for t in titre_lignes]
        self.lignes.append("# Genere par studio2/noyau — NE PAS EDITER A LA MAIN.")
        self.lignes.append("# " + "-" * 73)
        self.cles = {}

    def section(self, titre):
        self.lignes += ["", "# --- %s %s" % (titre, "-" * max(3, 66 - len(titre)))]

    def commentaire(self, texte):
        self.lignes.append("# " + texte)

    def cle(self, cle, valeur, note=None):
        if cle in self.cles:
            raise ValueError("clé écrite deux fois : %s" % cle)
        v = nombre(valeur)
        self.cles[cle] = v
        self.lignes.append("%s = %s" % (cle, v) + ("    # " + note if note else ""))

    def texte(self):
        return "\n".join(self.lignes) + "\n"
