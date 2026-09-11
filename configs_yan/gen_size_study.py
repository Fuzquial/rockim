#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# gen_size_study.py — PEUT-ON CALIBRER SUR UNE PETITE EPROUVETTE ?
#
#   python configs_yan/gen_size_study.py        # maillages + decks
#
# QUESTION. En FEM on calibre sur UN element (banc `matpoint`). En FEMDEM c'est
# impossible : la reponse vient des JOINTS entre elements, il en faut plusieurs.
# La question pratique est donc : a partir de combien d'elements la reponse
# globale cesse-t-elle de dependre de la taille ? Si une eprouvette d'une
# dizaine d'elements donne deja le bon pic et la bonne pente, la calibration
# coute 100 a 1000 fois moins cher.
#
# POURQUOI CE N'EST PAS ACQUIS. Mesure du 2026-09-11 : la reponse depend de la
# LONGUEUR COHESIVE l_ch,II = E G_fII / c^2 comparee a la taille de
# l'eprouvette (effet d'echelle de Bazant). Reduire l'eprouvette a taille de
# maille FIXE augmente donc le rapport l_ch/W et doit rendre la reponse plus
# ductile et plus resistante. L'etude mesure de combien.
#
# PROTOCOLE. Meme materiau (le cas `gfs5` du balayage : gfShearFactor = 5,
# l_ch,II = 39 mm), meme TAILLE DE MAILLE h = 1,5 mm, meme elancement 1:2,
# meme maillage NON STRUCTURE importe (Gmsh Delaunay + Netgen, `mesh = file` —
# aucune microstructure, regle de la these pour tout essai lu en facies).
# Seule la TAILLE de l'eprouvette change, d'un facteur 2 a chaque cran.
#
# LE PIEGE QUE CE SCRIPT EVITE : a `pullV` fixe, une eprouvette deux fois plus
# courte subit une vitesse de deformation deux fois plus grande. Les courbes ne
# seraient plus comparables. pullV est donc mis A L'ECHELLE de H pour tenir
# eps_point = 1,389 /s constant partout.
#
# LECTURE. Si q_pic et eps_pic sont plats de 3000 elements jusqu'a ~50, la
# calibration sur petite eprouvette est legitime et on gagne deux ordres de
# grandeur. S'ils derivent des 500, il faut calibrer a la taille d'usage.
# Les petites tailles sont jouees a TROIS GRAINES de maillage : a 12 elements
# la dispersion de realisation peut depasser l'effet cherche, et une moyenne
# sur une seule realisation n'aurait aucun sens.
# ---------------------------------------------------------------------------
import io
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

H_ELEM = 0.0015                 # taille de maille, IDENTIQUE partout
EPSDOT = 0.1 / 0.072            # 1,389 /s — la vitesse du banc de Yan

# (tag, W, H, graines)   elancement 1:2 conserve
SIZES = [
    ("L0", 0.036,  0.072,  [4211]),              # reference : ~3000 elements
    ("L1", 0.018,  0.036,  [4211]),              # ~760
    ("L2", 0.009,  0.018,  [4211, 7, 99]),       # ~190
    ("L3", 0.0045, 0.009,  [4211, 7, 99]),       # ~48
    ("L4", 0.00225, 0.0045, [4211, 7, 99]),      # ~12  <- la question posee
]

# materiau : le cas `gfs5` du balayage (le reglage retenu comme raisonnable)
MAT = dict(E=70e9, nu=0.25, rho=2600, ft=10e6, cohesion=25e6,
           frictionDeg=40, Gf=70.0, gfShearFactor=5.0)

HEADER = """# ---------------------------------------------------------------------------
# ETUDE DE TAILLE — UCS 2D FEMDEM homogene, {W:.1f} x {H:.1f} mm, graine {seed}
# Genere par configs_yan/gen_size_study.py (y lire la question posee).
#
# Taille de maille h = 1,5 mm IDENTIQUE a tous les crans ; seule l'EPROUVETTE
# change. pullV mis a l'echelle de H pour tenir eps_point = {ed:.3f} /s partout.
# Materiau : cas `gfs5` du balayage — l_ch,II = E G_fII / c^2 = 39,2 mm.
#   rapport l_ch,II / W = {ratio:5.1f}   (effet d'echelle attendu au-dela de ~1)
#
#   rockim_g1fix.exe configs_yan/{name}.cfg out_{name}
# ---------------------------------------------------------------------------
"""


def mesh_for(W, H, seed):
    name = "ucs2d_%dx%d_s%d.msh" % (round(W * 1e4), round(H * 1e4), seed)
    p = os.path.join(ROOT, "meshes", name)
    if not os.path.exists(p):
        subprocess.run([sys.executable,
                        os.path.join(ROOT, "tools", "make_unstructured_mesh.py"),
                        "box2d", "%g" % W, "%g" % H, "%g" % H_ELEM, p, str(seed)],
                       check=True, cwd=ROOT,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return "meshes/" + name


def deck(tag, W, H, seed):
    name = "size_%s_s%d" % (tag, seed)
    mesh = mesh_for(W, H, seed)
    lch2 = MAT["E"] * MAT["Gf"] * MAT["gfShearFactor"] / MAT["cohesion"] ** 2
    L = [HEADER.format(W=W * 1e3, H=H * 1e3, seed=seed, ed=EPSDOT,
                       ratio=lch2 / W, name=name)]
    a = L.append
    a("mode = fdem")
    a("scenario = tension")
    a("loading = platens")
    a("T = 3.5e-3")
    a("frames = 12")
    a("historyFlush = true")
    a("")
    a("W = %g" % W)
    a("H = %g" % H)
    a("thickness = 1.0")
    a("mesh = file")
    a("meshFile = %s" % mesh)
    a("")
    for k in ("rho", "E", "nu", "ft", "cohesion", "frictionDeg", "Gf",
              "gfShearFactor"):
        a("%s = %g" % (k, MAT[k]))
    a("")
    a("insertion = adaptive")
    a("jointSoftening = yan")
    a("insertionPenaltyFactor = 4")
    a("contactMu = 0.1")
    a("dampingLocal = 0.7")
    a("jointXi = 0.0")
    a("")
    # A L'ECHELLE : eps_point = pullV / H constant. Sans cela une eprouvette
    # 16 fois plus courte serait chargee 16 fois plus vite et la comparaison
    # melangerait effet de taille et effet de vitesse.
    a("pullV = -%g                  # eps_point = %.3f /s, comme tous les crans"
      % (EPSDOT * H, EPSDOT))
    a("pullRamp = 2e-4")
    a("gaugeLoFrac = 0.25")
    a("gaugeHiFrac = 0.75")
    a("ucsStopAfterPeak = true")
    a("ucsStopDelay = 2.0e-4")
    a("stopPeakDrop = 0.35")
    return name, "\n".join(L) + "\n"


if __name__ == "__main__":
    os.makedirs(os.path.join(ROOT, "meshes"), exist_ok=True)
    out = []
    print(" %-12s %10s %10s %10s %8s" % ("deck", "W x H (mm)", "l_ch/W", "pullV", "graine"))
    lch2 = MAT["E"] * MAT["Gf"] * MAT["gfShearFactor"] / MAT["cohesion"] ** 2
    for tag, W, H, seeds in SIZES:
        for s in seeds:
            name, txt = deck(tag, W, H, s)
            io.open(os.path.join(HERE, name + ".cfg"), "w",
                    encoding="utf-8", newline="\n").write(txt)
            out.append(name)
            print(" %-12s %5.2f x %-5.2f %9.1f %10.4f %8d"
                  % (name, W * 1e3, H * 1e3, lch2 / W, EPSDOT * H, s))
    print("\n%d decks ecrits" % len(out))
