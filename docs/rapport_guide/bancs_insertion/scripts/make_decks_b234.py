# -*- coding: utf-8 -*-
"""make_decks_b234.py — decks des bancs 2 (barre), 3 (plaque entaillee SENT), 4 (sensibilite)."""
import os
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M = os.path.join(R, "meshes")

BASE = """# bancs_insertion, banc {banc} — genere par make_decks_b234.py
mode = fdem
scenario = tension
mesh = file
meshFile = {mesh}
thickness = 1.0
T = {T}
frames = 5
rho = 2650
E = 50e9
nu = 0.25
ft = 10e6
cohesion = 25e6
frictionDeg = 40
Gf = 40
gfShearFactor = 10
jointSoftening = yan
pullV = 0.02
pullRamp = 2e-5
gripLateralFree = true
dampingLocal = 0.7
jointXi = 0.0
verifyFt = false
contact = potential
contactMu = 0.3
insertion = {ins}
{pen}
{contact}{extra}"""
CORR = "contactCandidates = vertex\ngcBirth = offset\npotForceExact = true\n"

def deck(path, banc, mesh, T, ins, corr, extra=""):
    pen = ("insertionPenaltyFactor = 100" if ins == "adaptive" else "jointPenaltyFactor = 100")
    open(path, "w").write(BASE.format(banc=banc, mesh=os.path.join(M, mesh), T=T, ins=ins,
                                      pen=pen, contact=CORR if corr else "", extra=extra))
    return path

def main():
    L = {2: [], 3: [], 4: []}
    hs = {"h10": "0.01", "h5": "0.005", "h2.5": "0.0025", "h1.25": "0.00125"}
    for hn, hv in hs.items():
        for ins in ("adaptive", "intrinsic"):
            for corr in (0, 1):
                n = "b2_%s_%s_%s" % (hn, ins[:5], "corr" if corr else "brut")
                L[2].append(deck(os.path.join(R, "b2", n + ".cfg"), 2, "bar_h%s.msh" % hv, 1.25e-3, ins, corr))
    for hn, hv in (("h2.5", "0.0025"), ("h1.25", "0.00125")):
        for ins in ("adaptive", "intrinsic"):
            n = "b3_%s_%s" % (hn, ins[:5])
            L[3].append(deck(os.path.join(R, "b3", n + ".cfg"), 3, "sent_h%s.msh" % hv, 2.0e-3, ins, 1))
    VAR = {"base": "", "elliptic": "insertionCriterion = elliptic\n",
           "volume": "facetAverage = volume\n", "max": "facetAverage = max\n",
           "hold5": "insertionHoldSteps = 5\n", "tip1.6": "insertionTipFactor = 1.6\n"}
    for hn, hv in (("h5", "0.005"), ("h2.5", "0.0025")):
        for v, x in VAR.items():
            n = "b4_%s_%s" % (hn, v)
            L[4].append(deck(os.path.join(R, "b4", n + ".cfg"), 4, "bar_h%s.msh" % hv, 1.25e-3, "adaptive", 1, x))
    for k, v in L.items():
        open(os.path.join(R, "b%d" % k, "liste.txt"), "w").write("\n".join(v) + "\n")
        print("banc", k, len(v), "decks")

if __name__ == "__main__":
    main()
