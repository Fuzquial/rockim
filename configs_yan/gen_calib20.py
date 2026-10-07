#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# gen_calib20.py — calibration du triaxial Red Bohus a sigma_3 = 20 MPa
#
#   python configs_yan/gen_calib20.py                 # lot de depistage (8)
#   python configs_yan/gen_calib20.py c=55 phi=45 ... # un point precis
#
# CIBLE (donnees nettoyees v4, CONTINUUM/calib_bohus_triax/exp_qc/
# experimental_data_red_bohus_clean.json, essais 2_2 / 2_3 / 2_4) :
#       q_pic   = 404,8 +- 2,8 MPa        (dispersion experimentale 0,7 %)
#       eps_pic = 0,661 %
#       E       = 76,3 +- 1,3 GPa
#
# EPROUVETTE : le cran L1 de l'etude de taille du 2026-09-11 — 18 x 36 mm,
# 766 triangles, maillage NON STRUCTURE importe (Gmsh Delaunay + Netgen),
# AUCUNE microstructure. Mesure qui l'autorise : contre l'eprouvette pleine
# 36 x 72 mm (3044 triangles), L1 rend le module a -0,4 %, le pic a +1,7 % et
# la montee entiere a 0,88 % RMS — pour 8 fois moins cher (54 s contre 7,3 min
# en UCS). RESERVE : le POST-PIC differe (L0 fait un snap-back de jauge que L1
# ne fait pas), donc l'objectif ne doit porter que sur MONTEE / PIC / MODULE,
# jamais sur l'adoucissement. Les points retenus se verifient sur L0.
#
# PARAMETRES, et pourquoi ceux-la. Le balayage UCS du 2026-09-11 (8 cas, meme
# maillage) a mesure les sensibilites :
#   cohesion  -> LE levier du pic          (c 25 -> 40 MPa : q_pic +46 %)
#   frictionDeg -> gain sous confinement   (non balaye ; a mesurer ici)
#   gfShearFactor (donc l_ch,II = E G_fII / c^2) -> eps_pic et raideur de
#                 l'adoucissement          (l_ch,II 173 -> 10 mm : eps_pic
#                 0,282 -> 0,156 %)
#   ft        -> SPECTATEUR en compression (ft 10 -> 20 MPa : q_pic +2 %,
#                 eps_pic +4 %, faciès inchange) — laisse fixe.
# C'est aussi le paramétrage de l'outillage du depot (calib_quick/calib/
# design.py : ft, c, phi, l_cz), a ceci pres qu'on fixe ft.
#
# LOT DE DEPISTAGE (8 runs) : 4 points en cohesion pour la pente dq/dc, 3 en
# frottement, 3 en l_ch,II. De quoi ajuster un modele local et viser.
# ---------------------------------------------------------------------------
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

BASE = dict(E=70e9, nu=0.25, rho=2600, ft=10e6,
            cohesion=25e6, frictionDeg=40, Gf=70.0, gfShearFactor=5.0)

W, H = 0.018, 0.036             # cran L1
MESH = "meshes/ucs2d_180x360_s4211.msh"
EPSDOT = 0.1 / 0.072            # 1,389 /s, la vitesse du banc de Yan
SIG3 = 20e6

# lot de depistage : (tag, surcharges)
SCREEN = [
    ("c25",   dict(cohesion=25e6)),
    ("c45",   dict(cohesion=45e6)),
    ("c65",   dict(cohesion=65e6)),
    ("c85",   dict(cohesion=85e6)),
    ("p30",   dict(cohesion=65e6, frictionDeg=30)),
    ("p50",   dict(cohesion=65e6, frictionDeg=50)),
    ("g15",   dict(cohesion=65e6, gfShearFactor=15.0)),
    ("g1p5",  dict(cohesion=65e6, gfShearFactor=1.5)),
]

HEAD = """# ---------------------------------------------------------------------------
# CALIBRATION Red Bohus, triaxial sigma_3 = 20 MPa — point « {tag} »
# Genere par configs_yan/gen_calib20.py (y lire la cible et le protocole).
#
# Cible : q_pic 404,8 +- 2,8 MPa | eps_pic 0,661 % | E 76,3 +- 1,3 GPa
# Ce point : c = {c:.0f} MPa, phi = {phi:.0f} deg, G_fII = {gf2:.0f} J/m2
#            l_ch,II = E G_fII / c^2 = {lch:5.1f} mm   (eprouvette W = 18 mm)
#
# Eprouvette L1 (18 x 36 mm, 766 triangles, maillage non structure importe) :
# rend le module a -0,4 % et le pic a +1,7 % de l'eprouvette pleine, pour 8
# fois moins cher. L'objectif ne porte PAS sur l'adoucissement (le post-pic
# de L1 differe : pas de snap-back de jauge).
#
#   rockim_g1fix.exe configs_yan/cal20_{tag}.cfg out_cal20_{tag}
# ---------------------------------------------------------------------------
"""


def deck(tag, over):
    m = dict(BASE)
    m.update(over)
    gf2 = m["Gf"] * m["gfShearFactor"]
    lch = m["E"] * gf2 / m["cohesion"] ** 2 * 1e3
    L = [HEAD.format(tag=tag, c=m["cohesion"] / 1e6, phi=m["frictionDeg"],
                     gf2=gf2, lch=lch)]
    a = L.append
    a("mode = fdem")
    a("scenario = tension")
    a("loading = platens")
    # eps_pic vise 0,66 % a 1,389 /s -> 4,75 ms, plus pullDelay 0,6 ms, plus
    # du post-pic : T = 8 ms. stopPeakDrop coupe avant si la chute est franche.
    a("T = 8.0e-3")
    a("frames = 16")
    a("historyFlush = true")
    a("")
    a("W = %g" % W)
    a("H = %g" % H)
    a("thickness = 1.0")
    a("mesh = file")
    a("meshFile = %s" % MESH)
    a("")
    for k in ("rho", "E", "nu", "ft", "cohesion", "frictionDeg", "Gf",
              "gfShearFactor"):
        a("%s = %g" % (k, m[k]))
    a("")
    a("insertion = adaptive")
    a("jointSoftening = yan")
    a("insertionPenaltyFactor = 4")
    a("contactMu = 0.1")
    a("dampingLocal = 0.7")
    a("jointXi = 0.0")
    a("")
    a("# --- chargement : confinement etabli AVANT la phase deviatoire ------")
    a("pullV = -%g" % (EPSDOT * H))
    a("pullRamp = 2e-4")
    a("pullDelay = 6.0e-4")
    a("confiningPressure = %g" % SIG3)
    a("confiningRamp = 2e-4")
    a("confineFaces = sides")
    a("confineGaugeTime = 5.0e-4")
    a("gaugeLoFrac = 0.25")
    a("gaugeHiFrac = 0.75")
    a("ucsStopAfterPeak = true")
    a("ucsStopDelay = 2.0e-4")
    a("stopPeakDrop = 0.35")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    if len(sys.argv) > 1:                      # point unique : cle=valeur
        over, tag = {}, []
        for kv in sys.argv[1:]:
            k, v = kv.split("=")
            over[k] = float(v)
            tag.append("%s%g" % (k[:3], float(v)))
        cases = [("_".join(tag), over)]
    else:
        cases = SCREEN
    print(" %-8s %8s %6s %8s %10s" % ("point", "c(MPa)", "phi", "G_fII", "l_ch,II(mm)"))
    for tag, over in cases:
        io.open(os.path.join(HERE, "cal20_%s.cfg" % tag), "w",
                encoding="utf-8", newline="\n").write(deck(tag, over))
        m = dict(BASE); m.update(over)
        gf2 = m["Gf"] * m["gfShearFactor"]
        print(" %-8s %8.0f %6.0f %8.0f %10.1f"
              % (tag, m["cohesion"] / 1e6, m["frictionDeg"], gf2,
                 m["E"] * gf2 / m["cohesion"] ** 2 * 1e3))
