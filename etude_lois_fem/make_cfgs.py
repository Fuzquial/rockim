# -*- coding: utf-8 -*-
"""Générateur de decks fem3d de l'étude « briques » avec LISTE BLANCHE de clés.

Le lecteur de configuration de rockim ignore les clés inconnues EN SILENCE ; ce
générateur refuse toute clé absente de la liste blanche de la loi choisie, et
impose les invariants (saksalaEta > 0, erodeDc ⇒ compDamage, blade ⇒ shear…).

usage (module) :
    from make_cfgs import deck, COMMON_T1, COMMON_T2
    deck("C_ref_P100", COMMON_T1, {"law": "dpr", "meridian": "power", ...}, out_dir)
"""
import io, os

COMMON_KEYS = {
    "mode", "scenario", "T", "frames", "outputDir", "historyFlush",
    "W", "D", "H", "nx", "ny", "nz", "meshJitter", "meshMirror", "seed", "geometry", "mesh", "meshFile",
    "rho", "E", "nu", "ft", "cohesion", "frictionDeg", "Gf", "gfShearFactor",
    "toolShape", "toolRadius", "toolMass", "impactSpeed", "toolGap", "toolX", "toolY",
    "cutDepth", "cutSpeed", "contactMu", "contactXi", "contactVreg",
    "absorbing", "absorbSpringFactor", "absorbSpringR", "dampingLocal", "dtFactor",
    "pullV", "pullRamp", "gripLateralFree",
    "matWeibullM", "fieldSeed", "strengthCorrLength", "strengthCorrLengthB", "strengthCorrAngleDeg",
    "weibullScope",
    # étude briques (solveur)
    "confiningPressure", "confiningRamp", "confineFaces", "topPressure", "bottomPressure", "confineGaugeTime",
    "toolDelay", "symmetryY", "activeNodes", "erodeDetMin", "erodeStrainMax", "fieldStats",
    "backRakeDeg", "clearanceDeg", "bladeHeight",
}
LAW_KEYS = {
    "elastic": set(),
    "dpr": {"erodeD", "erodeEpv", "meridian", "merB", "merN", "merFc0", "compDamage", "compAc",
            "compGIIc", "erodeDc", "erodeWfrac", "dpApex", "dpTension", "rankineDrive", "capP0", "capH", "dprCap"},
    "saksala": {"erodeD", "erodeEpv", "saksalaEta", "capP0", "capH", "meridian", "merB", "merN",
                "merFc0", "compDamage", "compAc", "compGIIc", "erodeDc", "erodeWfrac", "dpApex", "dpTension", "rankineDrive"},
    "mc": {"mcCohesion", "mcFrictionDeg", "mcDilationDeg"},
    "saksala2011": {"skBetaDP", "skCres", "skHdp", "skSdp", "skSmr", "skAt", "skBetaT", "skPp0",
                    "skPtr0", "skDcap", "skWcap", "skNd"},
    "dpdfh": {"dfhBetaDeg", "dfhDCoh", "dfhPsiDeg", "dfhWeibullM", "dfhSigW", "dfhZeff", "dfhK",
              "dfhS", "dfhDeld"},
}
FORBIDDEN_HINTS = {"psi": "dpr/saksala n'ont pas de dilatance (utiliser mcDilationDeg ou dfhPsiDeg)",
                   "jointFt": "clé fdem inexistante", "jointCoh": "clé fdem inexistante",
                   "mhN": "la loi mh n'est pas portée : meridian = power", "mhB": "idem"}


def check(keys):
    law = keys.get("law", "dpr")
    if law not in LAW_KEYS:
        raise ValueError("loi inconnue : %s" % law)
    allowed = COMMON_KEYS | LAW_KEYS[law] | {"law"}
    for k in keys:
        if k not in allowed:
            hint = FORBIDDEN_HINTS.get(k, "")
            raise ValueError("clé '%s' refusée pour law = %s (%s)" % (k, law, hint or "absente de la liste blanche"))
    if law == "saksala" and float(keys.get("saksalaEta", 0.05e6)) <= 0:
        raise ValueError("saksalaEta doit etre > 0 (eta = 1 Pa.s ~ dpr)")
    if float(keys.get("erodeDc", 0)) > 0 and keys.get("compDamage", "none") != "crackband":
        raise ValueError("erodeDc exige compDamage = crackband")
    if keys.get("toolShape") == "blade" and keys.get("scenario") != "shear":
        raise ValueError("toolShape = blade exige scenario = shear")
    if law == "dpr" and float(keys.get("capP0", 0)) > 0 and keys.get("dprCap") != "true":
        raise ValueError("capP0 sur dpr exige dprCap = true")
    if keys.get("meridian") == "power" and law not in ("dpr", "saksala"):
        raise ValueError("meridian = power n'existe que pour dpr/saksala")
    if float(keys.get("confiningPressure", 0)) > 0 and float(keys.get("absorbSpringFactor", 1)) > 0:
        raise ValueError("confinement : poser absorbSpringFactor = 0")
    return True


def deck(name, common, overrides, out_dir, comment=""):
    keys = dict(common)
    keys.update(overrides)
    keys["outputDir"] = "out_" + name
    check(keys)
    lines = ["# %s" % name]
    if comment:
        lines.append("# " + comment)
    for k, v in keys.items():
        lines.append("%s = %s" % (k, v))
    path = os.path.join(out_dir, name + ".cfg")
    io.open(path, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return path


# ---- gabarits de l'étude (SI) -------------------------------------------------
BOHUS = {"rho": 2620, "E": "77.66e9", "nu": 0.29, "ft": "9e6", "cohesion": "22.77e6",
         "frictionDeg": 50.4, "Gf": 100, "gfShearFactor": 10}

# T1 percussion 48 x 48 x 32 mm, Gmsh Delaunay gradue (0,5 mm dans r < 8 mm -> 3 mm a 30 mm ; convergence W_abs 6,7/9,5/12,0/12,9 J a 1,5/1,0/0,75/0,5 mm), R 7,94 mm, 16 J, Lysmer partout, ressorts nuls
COMMON_T1 = dict(mode="fem3d", scenario="percussion", T="3.0e-4", frames=6,
                 mesh="file", meshFile="../meshes/T1_c05.msh",
                 **BOHUS,
                 toolShape="sphere", toolRadius="7.94e-3", toolMass=0.5, impactSpeed=8.0,
                 toolGap="1e-4", contactMu=0.5, contactXi=0.05,
                 absorbing="all", absorbSpringFactor=0, dampingLocal=0.05, dtFactor=0.7,
                 fieldStats="true", activeNodes="true", erodeDetMin=0.3, erodeStrainMax=0.3, erodeEpv=0, erodeD=2, erodeWfrac=0.98)

# T2 coupe en tranche 22 x 2 x 10 mm, Gmsh Delaunay uniforme 0,5 mm, lame rake 20, 1,5 mm, 4 m/s
COMMON_T2 = dict(mode="fem3d", scenario="shear", T="3.5e-3", frames=6,
                 mesh="file", meshFile="../meshes/T2_h05.msh", symmetryY="true",
                 **BOHUS,
                 toolShape="blade", backRakeDeg=20, clearanceDeg=10, bladeHeight=0.02,
                 cutDepth="1.5e-3", cutSpeed=4.0, toolX="-1e-4", contactMu=0.4, contactXi=0.05,
                 absorbing="sides", absorbSpringFactor=0, dampingLocal=0.05, dtFactor=0.7,
                 fieldStats="true", activeNodes="true", erodeDetMin=0.3, erodeStrainMax=0.3, erodeEpv=0, erodeD=2, erodeWfrac=0.98)

# mesh = file : hmin = plus petit diametre inscrit (0,39 h_in ~ h/2,6 pour un tet de Kuhn) -> dtFactor 0,7
# sur cette base equivaut au 0,3 sur l'ARETE que la grille utilisait (marge identique, dt x2,3 plus grand).
# la référence R des briques
REF_R = {"law": "dpr", "dpApex": "true", "dpTension": "off", "rankineDrive": "stress", "meridian": "power", "merB": 56.59, "merN": 0.538,
         "compDamage": "crackband", "compAc": 0.98, "compGIIc": "1e4"}


def confined(P_MPa, T=None):
    """clés de confinement pour P (MPa) : rampe 30 µs, jauge et outil à 60 µs ;
    T = durée totale (60 µs de mise en pression + la durée d'impact du deck)."""
    if P_MPa <= 0:
        return {}
    d = {"confiningPressure": "%ge6" % P_MPa, "confiningRamp": "30e-6",
         "confineGaugeTime": "60e-6", "toolDelay": "60e-6"}
    if T is not None:
        d["T"] = T
    return d


if __name__ == "__main__":
    # auto-test de la liste blanche : ces trois decks DOIVENT être refusés
    for bad in ({"law": "dpr", "psi": 15}, {"law": "saksala", "saksalaEta": 0},
                {"law": "dpr", "erodeDc": 0.99}):
        try:
            check(dict(COMMON_T1, **bad))
            print("ERREUR : accepté", bad)
        except ValueError as e:
            print("refuse comme prevu :", str(e).encode("ascii", "replace").decode())
    check(dict(COMMON_T1, **REF_R, **confined(100)))
    print("reference R + P100 : OK")
