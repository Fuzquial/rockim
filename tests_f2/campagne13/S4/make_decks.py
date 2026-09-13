# -*- coding: utf-8 -*-
"""make_decks.py — decks des essais elementaires S4 (campagne du 13/09).

Boite fdem3d 4 x 4 x 4 mm, mesh = grid 2 x 2 x 2 (48 tetraedres de Kuhn,
arete a = 2 mm : h inscrit = 0,41421 a = 0,8284 mm, arete moyenne =
(3 + 2 sqrt 2 + sqrt 3)/6 a = 1,26008 a = 2,5202 mm), joints incassables
(ft = cohesion = 1e12), crushCap et meanTensionCap neutralises,
bulkDamage = yang avec la calibration Kuru (delta0 = 14 um, deltaF = 400 um).

Trois chargements :
  iso  : compression ISOTROPE — scenario = percussion sans outil, base
         bloquee en z seulement (gripLateralFree), pression suiveuse sur les
         6 faces (confineFaces = all). dev(eps) = 0 : D doit rester 0 sous
         `deviatoric`, s armer sous `total` et `principal`.
  biax : confinement LATERAL seul (confineFaces = lateral, sommet libre) :
         etat deviatorique a trace (sigma_xx = sigma_yy = -p, sigma_zz = 0).
  uni  : compression UNIAXIALE — scenario = tension, pullV < 0, mors libres
         lateralement.
Variantes de mesure : ref (aucune cle S4, temoin bit-identique), dev (probe
seule : physique identique a ref), tot, pri, edge, edge_tot ; refus attendus :
bad_nokey (bulkDamageStrain sans bulkDamage), bad_value (valeur inconnue).

usage : python make_decks.py <dossier de sortie>
"""
import io, os, sys

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

COMMON = """mode = fdem3d
mesh = grid
W = 0.004
D = 0.004
H = 0.004
nx = 2
ny = 2
nz = 2
frames = 25
rho = 2620
E = 52e9
nu = 0.25
ft = 1e12
cohesion = 1e12
frictionDeg = 30
Gf = 1e6
gfShearFactor = 10
crushCap = 1e12
meanTensionCapFactor = 0
jointPenaltyFactor = 20
dampingLocal = 0.7
gripLateralFree = true
verifyFt = false
bulkDamage = yang
bulkDamageDelta0 = 1.4e-5
bulkDamageDeltaF = 4.0e-4
bulkDamageDmax = 0.9
"""

LOAD = {
    "iso": """scenario = percussion
toolShape = none
absorbing = none
confiningPressure = 2.5e9
confineFaces = all
confiningRamp = 1.0e-4
T = 2.5e-4
""",
    "biax": """scenario = percussion
toolShape = none
absorbing = none
confiningPressure = 2.5e9
confineFaces = lateral
confiningRamp = 1.0e-4
T = 2.5e-4
""",
    "uni": """scenario = tension
pullV = -0.6
pullRamp = 1.0e-4
T = 2.5e-4
""",
}

MEAS = {
    "ref": "",
    "dev": "bulkDamageProbe = true\n",
    "tot": "bulkDamageProbe = true\nbulkDamageStrain = total\n",
    "pri": "bulkDamageProbe = true\nbulkDamageStrain = principal\n",
    "edge": "bulkDamageProbe = true\nbulkDamageLength = edge\n",
    "edge_tot": "bulkDamageProbe = true\nbulkDamageLength = edge\nbulkDamageStrain = total\n",
}

def write(name, text):
    with io.open(os.path.join(OUT, name + ".cfg"), "w", encoding="ascii", newline="\n") as f:
        f.write("# S4 (campagne du 13/09) - essai elementaire " + name + "\n")
        f.write(text)

for ld, lt in LOAD.items():
    for ms, mt in MEAS.items():
        if ld != "uni" and ms in ("edge",):
            pass
        write(ld + "_" + ms, COMMON + lt + mt)

# refus attendus
write("bad_nokey", COMMON.replace("bulkDamage = yang\n", "").replace(
    "bulkDamageDelta0 = 1.4e-5\nbulkDamageDeltaF = 4.0e-4\nbulkDamageDmax = 0.9\n", "")
      + LOAD["uni"] + "bulkDamageStrain = total\n")
write("bad_value", COMMON + LOAD["uni"] + "bulkDamageStrain = vonmises\n")
write("bad_len", COMMON + LOAD["uni"] + "bulkDamageLength = min\n")

# miroir 2D : plaque fdem 4 x 4 mm, grid, compression uniaxiale par mors
COMMON2D = COMMON.replace("mode = fdem3d\n", "mode = fdem\n").replace(
    "D = 0.004\n", "").replace("nz = 2\n", "")
for ms in ("ref", "dev", "tot", "pri", "edge"):
    write("uni2d_" + ms, COMMON2D + LOAD["uni"] + MEAS[ms])
print("decks ecrits dans", OUT)
