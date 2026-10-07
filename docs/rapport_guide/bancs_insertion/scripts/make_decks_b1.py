# -*- coding: utf-8 -*-
"""make_decks_b1.py — banc 1 : un joint isole dans une bande de 8 CST (mesh_b1.py).

Trois trajets (traction pure theta = 0 ; cisaillement pur theta = 45 deg sous compression
uniaxiale, phi = 0 donc f_s = c et aucun frottement ; mode mixte theta = 45 deg en traction,
sigma_n = tau) x deux schemas (adaptive, intrinsic) x deux penalites (10 E/h, 100 E/h)
x deux adoucissements (linear, yan = z-curve de Munjiza / eq. 11 de Yan 2023).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "b1")

BASE = """# banc 1 (bancs_insertion) : joint isole, genere par make_decks_b1.py
mode = fdem
scenario = tension
mesh = file
meshFile = {mesh}
thickness = 1.0
T = {T}
frames = 120
rho = 2650
E = 50e9
nu = 0.0
ft = 10e6
cohesion = 25e6
frictionDeg = 0
Gf = 70
gfShearFactor = 10
pullV = {v}
pullRamp = 2e-5
gripLateralFree = true
dampingLocal = 0.7
jointXi = 0.0
contactMu = 0.0
jointResidualMu = 0
jointSoftening = {soft}
verifyFt = false
insertion = {ins}
{pen}
"""
PATHS = {  # nom : (maillage, pullV, T)
    "I":   ("b1_t0.msh", 0.02, 2.0e-3),
    "II":  ("b1_t45.msh", -0.05, 3.0e-3),
    "mix": ("b1_t45.msh", 0.02, 3.0e-3),
}

def main():
    os.makedirs(OUT, exist_ok=True)
    names = []
    for p, (mesh, v, T) in PATHS.items():
        for ins in ("adaptive", "intrinsic"):
            for pf in (10, 100):
                for soft in ("linear", "yan"):
                    key = "insertionPenaltyFactor" if ins == "adaptive" else "jointPenaltyFactor"
                    name = "b1_%s_%s_p%d_%s" % (p, ins[:5], pf, soft[:3])
                    txt = BASE.format(mesh=os.path.join(ROOT, "meshes", mesh), T=T, v=v,
                                      soft=soft, ins=ins, pen="%s = %d" % (key, pf))
                    fn = os.path.join(OUT, name + ".cfg")
                    open(fn, "w").write(txt)
                    names.append(fn)
    open(os.path.join(OUT, "liste.txt"), "w").write("\n".join(names) + "\n")
    print(len(names), "decks")

if __name__ == "__main__":
    main()
