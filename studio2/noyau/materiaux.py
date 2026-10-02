"""Préréglages de matériau : granite Red Bohus en une ou trois phases minérales.

Repris de etude_triax_hetero/gen_decks.py (§2), qui documente les sources :
fractions massiques de Dumoulin 2024 converties en fractions SURFACIQUES,
constantes minérales du dépôt, puis mise à l'échelle de chaque propriété par un
facteur unique pour que la moyenne pondérée par les aires tombe sur le niveau
choisi (« bohus » calibré, ou « fragile » de Yan 2023 Table 1). Le contraste
entre minéraux est donc celui de la littérature, le niveau celui du matériau.

Toutes les valeurs rendues sont arrondies à 6 chiffres significatifs, la
précision des decks de la campagne : un deck écrit par le noyau relit ainsi les
mêmes nombres.
"""
from .essai import JointsGrain, Materiau, Paire, Phase

NIVEAUX = {
    "bohus":   dict(E=52e9, rho=2624.0, ft=34e6, cohesion=13.6e6, frictionDeg=13.4,
                    Gf=70.0, gfs=10.0, crushCap=400e6, gbAlpha=0.5, nu=None),
    "fragile": dict(E=15e9, rho=1704.0, ft=1.3e6, cohesion=16.4e6, frictionDeg=23.0,
                    Gf=3.8, gfs=22.105, crushCap=None, gbAlpha=1.0, nu=0.24),
}

_MASSE = dict(feldspar=0.60, quartz=0.35, biotite=0.05)
MINERAUX = {
    "quartz":   dict(E=94e9, nu=0.08, rho=2650, ft=13e6, cohesion=30e6, frictionDeg=45, Gf=90, taille=1.0),
    "feldspar": dict(E=70e9, nu=0.29, rho=2560, ft=10e6, cohesion=25e6, frictionDeg=40, Gf=70, taille=2.0),
    "biotite":  dict(E=34e9, nu=0.25, rho=3050, ft=7e6, cohesion=15e6, frictionDeg=30, Gf=40, taille=0.13),
}
# Interfaces par paire, valeurs de référence avant mise à l'échelle (le mica clive,
# quartz-quartz tient) : (a, b, ft, cohesion, Gf, frictionDeg).
_PAIRES = (("biotite", "biotite", 2.0e6, 5.0e6, 12.0, 18.0),
           ("biotite", "feldspar", 3.5e6, 8.0e6, 20.0, 24.0),
           ("biotite", "quartz", 3.5e6, 8.0e6, 20.0, 24.0),
           ("quartz", "quartz", 9.0e6, 21.0e6, 60.0, 42.0))

_v = {k: _MASSE[k] / MINERAUX[k]["rho"] for k in _MASSE}
FRACTIONS = {k: _v[k] / sum(_v.values()) for k in _v}


def arrondi(x):
    return float("%g" % x)


def _echelle(prop, cible):
    return cible / sum(FRACTIONS[k] * MINERAUX[k][prop] for k in MINERAUX)


def _phases_brutes(niveau):
    c = NIVEAUX[niveau]
    ph = {k: dict(fraction=FRACTIONS[k], nu=v["nu"]) for k, v in MINERAUX.items()}
    for prop in ("E", "rho", "ft", "cohesion", "frictionDeg", "Gf"):
        e = _echelle(prop, c[prop])
        for k in ph:
            ph[k][prop] = MINERAUX[k][prop] * e
    return ph


def materiau(niveau="bohus"):
    """Matériau de volume homogène : moyenne pondérée par les aires (Voigt)."""
    c, ph = NIVEAUX[niveau], _phases_brutes(niveau)
    moy = {k: sum(p["fraction"] * p[k] for p in ph.values())
           for k in ("E", "nu", "rho", "ft", "cohesion", "frictionDeg", "Gf")}
    if c["nu"]:
        moy["nu"] = c["nu"]
    return Materiau(nom=niveau, gfShearFactor=c["gfs"], crushCap=c["crushCap"],
                    **{k: arrondi(v) for k, v in moy.items()})


def phases(niveau="bohus", contraste_resistance=True, tailles=False, taille_grain=0.003):
    """Les trois minéraux de Red Bohus, dans l'ordre quartz, feldspath, biotite.

    contraste_resistance=False : seul le contraste ÉLASTIQUE (E, nu, rho) est posé,
    les résistances restent celles du volume. tailles=True : affinité de taille
    par phase, sur-spécifiée (feldspath > quartz > biotite en paillettes).
    """
    out = []
    for nom, p in _phases_brutes(niveau).items():
        res = {k: arrondi(p[k]) for k in ("ft", "cohesion", "frictionDeg", "Gf")} \
            if contraste_resistance else {}
        out.append(Phase(nom=nom, fraction=arrondi(p["fraction"]), rho=arrondi(p["rho"]),
                         E=arrondi(p["E"]), nu=p["nu"],
                         taille_grain=arrondi(MINERAUX[nom]["taille"] * taille_grain) if tailles else None,
                         **res))
    return out


def joints_grain(niveau="bohus", paires=False):
    """Atténuation des joints de grain du niveau, avec ou sans surcharges par paire.

    Les paires sont ramenées au niveau BOHUS quel que soit le niveau demandé :
    c'est le choix de gen_decks.py, où elles ne servent qu'aux séries Bohus.
    """
    a = NIVEAUX[niveau]["gbAlpha"]
    jg = JointsGrain(alpha_ten=a, alpha_coh=a, alpha_gf=a, alpha_e=1.0, alpha_fric=1.0)
    if paires:
        b = NIVEAUX["bohus"]
        s = {k: _echelle(k, b[k]) for k in ("ft", "cohesion", "Gf", "frictionDeg")}
        jg.paires = [Paire(a=x, b=y, ft=arrondi(ft * s["ft"]), cohesion=arrondi(coh * s["cohesion"]),
                           Gf=arrondi(gf * s["Gf"]), frictionDeg=arrondi(phi * s["frictionDeg"]))
                     for x, y, ft, coh, gf, phi in _PAIRES]
    return jg
