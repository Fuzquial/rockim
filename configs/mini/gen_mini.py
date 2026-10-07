# Genere les decks du banc "petit assemblage" (2026-09-11) : meme carte matrice
# que configs/loi_note_2026.cfg, bloc 6 mm en grille 3x3x3, bouton R = 2 mm.
import itertools, os
BASE = """mode = fdem3d
scenario = percussion
T = {T}
frames = 20
historyFlush = true
mesh = grid
W = {W}
D = {W}
H = {W}
nx = {n}
ny = {n}
nz = {n}
rho = 2630
E = 60e9
nu = 0.2
{matrix}
insertion = adaptive
ft = 11.4e6
Gf = 100
gfShearFactor = 10
cohesion = 30e6
frictionDeg = 35
facetAverage = volume
facetRate = tensor
insertionCriterion = elliptic
insertionHoldSteps = 3
strainRateDIF = yang-fig2
difExpT = 0.333
{joints}
insertionPenaltyFactor = 10
contactMu = 0.55
jointFrictionMobilised = damage
jointViscousInCriterion = off
gbCombine = min
gbAlphaTen = 0.7
gbAlphaGf = 0.7
gbAlphaCoh = 0.7
jointWeibullM = 6
jointWeibullXu = 0.3
jointWeibullScale = area
jointSizeEffect = true
seed = 12345
toolShape = sphere
toolRadius = {R}
toolMass = {M}
impactSpeed = 8.0
toolGap = 1e-5
absorbing = all
dampingLocal = 0
jointXi = 0
dtFactor = 0.15
energyBreakdown = on
"""
MAT_NOTE = """law = saksala
meridian = power
merB = 30.1
merN = 0.676
merFc0 = 235e6
saksalaEta = 1.5e6
capP0 = 0
mhForm = principal
dpFlowForm = mc
dpDilationDeg = 5
compDamage = crackband
compAc = 0.98
compGIIc = 1.0e4
compBandLength = tetEdge
bulkTensionDamage = off"""
MAT_ELAS = "crushCap = 1e12"
JOINTS = {
 "pen":   "jointTSL = penalty",
 "cam":   "jointTSL = camacho\njointTSLRise = 1e-3\njointMixLaw = bk\njointBKEta = 2.0",
 "camnobk": "jointTSL = camacho\njointTSLRise = 1e-3\njointMixLaw = none",
 "cam0":  "jointTSL = camacho\njointTSLRise = 0\njointMixLaw = bk\njointBKEta = 2.0",
 "cam2":  "jointTSL = camacho\njointTSLRise = 1e-2\njointMixLaw = bk\njointBKEta = 2.0",
}
import sys
n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
W = float(sys.argv[2]) if len(sys.argv) > 2 else 0.006
T = sys.argv[3] if len(sys.argv) > 3 else "3e-5"
R = sys.argv[4] if len(sys.argv) > 4 else "0.002"
M = sys.argv[5] if len(sys.argv) > 5 else "0.01"
tag = f"n{n}"
for mat, mname in (("note", MAT_NOTE), ("elas", MAT_ELAS)):
    for jk, jv in JOINTS.items():
        fn = f"configs/mini/mini_{tag}_{mat}_{jk}.cfg"
        with open(fn, "w", newline="\n") as f:
            f.write(BASE.format(T=T, W=W, n=n, matrix=mname, joints=jv, R=R, M=M))
        print(fn)
