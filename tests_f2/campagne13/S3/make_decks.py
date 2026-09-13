# -*- coding: utf-8 -*-
"""make_decks.py — decks du mini-test S3 (banc de joint cinematique, scenario = jointbench).

Base commune : joint "Kuru" du deck yang2026_bench_s25_v3P (ft 10,98 MPa, c 29,84 MPa,
Gf 50, GfII 1000, tan phi 1,85), penalite 25 E/edge (edge = arete du tetra regulier = jbEdge),
conventions de Solidity (parabole, z-curve de Munjiza, milieux d aretes, regle 2/3), sans
dashpot. Neuf variantes : trois lois x deux trajets (traction avec decharge a 0,5 ; cisaillement
sous 60 MPa de compression avec decharge a 0,4) + trois bancs de basculement (regle 2/3).

Ce que chaque deck DOIT donner (critere falsifiant, imprime par le solveur [JOINTBENCH]) :
  ten_solidity : (i) PASS (retrace), (iii) W_I = 1,159 Gf (+/- 2 %)
  ten_plastic  : (i) FAIL (secante de decharge), (iii) W_I = Gf (+/- 2 %)
  ten_origin   : idem plastic en mode I (origin ne touche que le mode II)
  sh_plastic   : (ii) glissement residuel > 0 ; (i) FAIL
  sh_solidity  : (ii) residuel = 0 ; (i) PASS
  sh_origin    : (ii) residuel = 0 (secante a l origine) ; (i) FAIL contre la charge
  tilt_sol     : (iv) trois instants distincts, joint mort au 2e point (majority)  -> PASS
  tilt_pl_maj  : idem sous plastic + majority                                      -> PASS
  tilt_pl_any  : jointFailRule = any : mort au PREMIER point (variante qui doit differer)
"""
import io, os

HERE = os.path.dirname(os.path.abspath(__file__))

BASE = """# S3 (campagne du 13/09) : banc de joint cinematique - genere par make_decks.py
mode = fdem3d
scenario = jointbench
frames = 2
# joint "Kuru" (Yang et al. 2026, Table 1), le materiau du deck s = 2,5
rho = 2626
E = 60e9
nu = 0.24
ft = 10.98e6
cohesion = 29.84e6
frictionDeg = 61.61
Gf = 50
gfShearFactor = 20
# penalite 25 E/h avec h = arete de la facette (p0 = 2 pf E = 3000 GPa, cadrage S3)
jointPenaltyFactor = 25
jointPenaltyLength = edge
# conventions de Solidity communes aux trois lois comparees
jointElastic = parabolic
jointSoftening = munjiza
jointQuadrature = midedge
jointXi = 0
dampingLocal = 0
# banc
jbEdge = 1e-3
jbRate = 1.0
jbSteps = 4000
"""

LAW = {
    "solidity": """jointShearUnload = solidity
jointDeltaC = solidity
jointFailRule = majority
""",
    "plastic": """jointShearUnload = plastic
jointShearRange = coulomb
jointSecantRatchet = on
jointNormalProxy = law
jointResidualMu = 0
jointDeath = damage
jointFailRule = majority
""",
    "origin": """jointShearUnload = origin
jointShearRange = coulomb
jointSecantRatchet = on
jointNormalProxy = law
jointResidualMu = 0
jointDeath = damage
jointFailRule = majority
""",
}

PATH = {
    # traction : rupture a dnF (13,7 um solidity ; 11,8 um exact), decharge a 10 um
    "ten": """jbMode = tension
jbAmp = 2.0e-5
jbUnloadAt = 0.5
""",
    # cisaillement sous 60 MPa de compression (2 pj dn, dn = -2e-8 m) : plage
    # coulomb ~21 um (solidity) / ~18,5 um (plastic) ; decharge a 12 um
    "sh": """jbMode = shear
jbNormal = -2.0e-8
jbAmp = 3.0e-5
jbUnloadAt = 0.4
""",
    # basculement : trois bras distincts, pas de decharge
    "tilt": """jbMode = tension
jbAmp = 2.0e-5
jbTilt = 0.5
""",
}

DECKS = {
    "ten_solidity": ("ten", "solidity", ""),
    "ten_plastic":  ("ten", "plastic", ""),
    "ten_origin":   ("ten", "origin", ""),
    "sh_solidity":  ("sh", "solidity", ""),
    "sh_plastic":   ("sh", "plastic", ""),
    "sh_origin":    ("sh", "origin", ""),
    "tilt_sol":     ("tilt", "solidity", ""),
    "tilt_pl_maj":  ("tilt", "plastic", ""),
    "tilt_pl_any":  ("tilt", "plastic", "jointFailRule = any\n"),
}


def main():
    for name, (path, law, extra) in DECKS.items():
        txt = BASE + "# --- loi ---\n" + LAW[law] + "# --- trajet ---\n" + PATH[path]
        if extra:
            # la derniere occurrence d une cle gagne (Config) : on retire la
            # ligne de la loi pour ne pas laisser deux valeurs dans le deck
            txt = txt.replace("jointFailRule = majority\n", "")
            txt += "# --- variante qui DOIT differer ---\n" + extra
        with io.open(os.path.join(HERE, name + ".cfg"), "w", encoding="utf-8", newline="\n") as f:
            f.write(txt)
        print("ecrit", name + ".cfg")


if __name__ == "__main__":
    main()
