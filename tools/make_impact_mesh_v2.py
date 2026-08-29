#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_impact_mesh.py — maillage de l'essai d'impact a insert unique
# (Yang et al. 2025-2026, donnees Mines Paris ; spec 005, WP4).
#
#   python tools/make_impact_mesh.py meshes/impact_s15.msh 1.5
#
# Quatre corps, volumes physiques nommes : rock, insert (carbure), bit,
# piston (acier). L'insert et le bit sont mailles CONFORMEMENT (fragment
# OCC) — leur interface recoit des joints par groupBond.bit.insert = joints.
# Le piston vole a 0,2 mm au-dessus du bit, l'insert a 0,2 mm de la roche :
# le contact general fait le reste.
#
# Geometrie (leur fig. 5, en m) : roche cylindre R 0,125 x 0,150 ; insert
# hemisphere R 8,51 mm + fut Phi 15,88 (23,2 mm au total) ; bit Phi 30,
# 265 mm insert compris ; piston Phi 26,5 x 260. SIMPLIFICATION V1 : la
# plaque de charge et le circlip sont OMIS — l'article precise que le poids
# sur l'outil n'est pas applique dans leurs simulations ; leur role de
# distribution de masse est secondaire pour le facies de fissuration.
#
# Tailles de maille (leur fig. 6) x le facteur d'echelle s :
#   roche : 1 mm dans la boule R 12,5 mm sous l'impact, 2 mm jusqu'a
#   R 25 mm, 10 mm au bord ; insert 0,7 mm ; bit 3 mm ; piston 5 mm.
#   s = 1 reproduit l'article (~230 k tets) ; s = 1,5 est la variante
#   econome pour les runs de nuit.
# ---------------------------------------------------------------------------
import sys

import gmsh

out = sys.argv[1] if len(sys.argv) > 1 else "impact.msh"
s = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
# jeu insert/roche optionnel (3e argument, m). Defaut 0,2 mm = insert en vol ;
# ~0,02 mm = insert POSE sur la roche, la configuration de l'essai reel (le
# poids sur l'outil assied le bit) — c'est elle qui donne la courbe F-p en
# rampe de leur fig. 7b, l'onde de frappe enfoncant l'insert directement.
GAPR = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0e-4
# echelle PROPRE A LA ROCHE (4e argument, defaut = s). Les fissures vivent
# dans la roche : on peut la mailler a l'echelle de l'article (1,0) en
# laissant l'acier grossier — le dt reste commande par l'insert.
SR = float(sys.argv[4]) if len(sys.argv) > 4 else s
# mode "leger" (5e argument) : rock + insert + bit SEULEMENT — pour les runs
# de pulverisation ou seul le chargement de la roche compte, le corps
# bit+insert etant lance directement a la vitesse d indentation mesuree.
LEGER = len(sys.argv) > 5 and sys.argv[5] == "leger"
# mode "roche" (5e argument, 2026-08-28) : la ROCHE SEULE — l'outil est alors
# l'hemisphere ANALYTIQUE du solveur (toolShape = sphere, toolRadius 8,51 mm,
# toolMass, impactSpeed) : aucun element de carbure, donc le dt est commande
# par la roche. GAPR et les tailles acier deviennent sans objet.
ROCHE = len(sys.argv) > 5 and sys.argv[5] == "roche"

GAP = 2.0e-4                 # jeu piston/bit [m]
R_ROCK, H_ROCK = 0.125, 0.150
R_INS, R_SHANK, H_INS = 0.00851, 0.00794, 0.0232
R_BIT = 0.015
L_BIT = 0.265 - H_INS        # le bit fait 265 mm INSERT COMPRIS
R_PIS, L_PIS = 0.01325, 0.260
# circlip (bague carbure brasee au bit) et plaque de charge (acier) : leur
# fig. 5c-d. La plaque 119,9 x 40 x 6 mm percee a Phi 31 REPOSE sur le
# circlip (jeu 0,05 mm) : son poids charge le bit par gravite, comme au banc.
R_CLIP, H_CLIP, Z_CLIP = 0.018, 0.003, 0.030
PL_X, PL_Y, PL_H, R_HOLE = 0.1199, 0.040, 0.006, 0.0155

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("impact")
occ = gmsh.model.occ

rock = occ.addCylinder(0, 0, -H_ROCK, 0, 0, H_ROCK, R_ROCK)
zc = GAPR + R_INS                      # centre de l'hemisphere (pointe a GAPR)
if not ROCHE:
    sph = occ.addSphere(0, 0, zc, R_INS)
    shank = occ.addCylinder(0, 0, zc, 0, 0, GAPR + H_INS - zc, R_SHANK)
    ins = occ.fuse([(3, sph)], [(3, shank)])[0]
zb0 = GAPR + H_INS                     # bas du bit = haut de l'insert
bit = occ.addCylinder(0, 0, zb0, 0, 0, L_BIT, R_BIT) if not ROCHE else None
zp0 = zb0 + L_BIT + GAP
pis = None
clip = None
zpl = Z_CLIP + H_CLIP + 5.0e-5
if not LEGER and not ROCHE:
    pis = occ.addCylinder(0, 0, zp0, 0, 0, L_PIS, R_PIS)
    clip_o = occ.addCylinder(0, 0, Z_CLIP, 0, 0, H_CLIP, R_CLIP)
    clip_i = occ.addCylinder(0, 0, Z_CLIP, 0, 0, H_CLIP, R_BIT)
    clip = occ.cut([(3, clip_o)], [(3, clip_i)])[0]
    box = occ.addBox(-PL_X / 2, -PL_Y / 2, zpl, PL_X, PL_Y, PL_H)
    hole = occ.addCylinder(0, 0, zpl, 0, 0, PL_H, R_HOLE)
    plate = occ.cut([(3, box)], [(3, hole)])[0]

# conformite insert/bit ET bit/circlip (faces partagees) ; le reste au contact
if not ROCHE:
    all3 = occ.fragment(ins + (clip or []), [(3, bit)])[0]
occ.synchronize()

vols = gmsh.model.getEntities(3)
names = {}
for dim, tag in vols:
    x, y, z = occ.getCenterOfMass(dim, tag)
    # le centre de masse d'un ANNEAU est sur l'axe : on classe par cote z
    if z < 0:                                   nm = "rock"
    elif z < zb0:                               nm = "insert"
    elif Z_CLIP - 1e-6 < z < zpl - 1e-9:        nm = "circlip"
    elif zpl - 1e-9 <= z < zpl + PL_H:          nm = "plate"
    elif z < zp0 - 1e-6:                        nm = "bit"
    else:                                       nm = "piston"
    names.setdefault(nm, []).append(tag)
for nm, tags in names.items():
    p = gmsh.model.addPhysicalGroup(3, tags)
    gmsh.model.setPhysicalName(3, p, nm)
want = ({"rock"} if ROCHE
        else {"rock", "insert", "bit"} if LEGER
        else {"rock", "insert", "bit", "piston", "circlip", "plate"})
assert set(names) == want, names

# ---- tailles --------------------------------------------------------------
# raffinement de la roche : boules concentriques sous le point d'impact
RESTRICT = True
f1 = gmsh.model.mesh.field.add("Ball")
gmsh.model.mesh.field.setNumber(f1, "VIn", 0.001 * SR)
gmsh.model.mesh.field.setNumber(f1, "VOut", 1.0)
gmsh.model.mesh.field.setNumber(f1, "Radius", 0.0125)
f2 = gmsh.model.mesh.field.add("Ball")
gmsh.model.mesh.field.setNumber(f2, "VIn", 0.002 * SR)
gmsh.model.mesh.field.setNumber(f2, "VOut", 1.0)
gmsh.model.mesh.field.setNumber(f2, "Radius", 0.025)
# taille de fond : graduation 2 mm -> 10 mm entre R 25 et R 100 mm
fd = gmsh.model.mesh.field.add("Distance")
pt0 = occ.addPoint(0, 0, 0)
occ.synchronize()
gmsh.model.mesh.field.setNumbers(fd, "PointsList", [pt0])
f3 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f3, "InField", fd)
import os as _os
SRFAR = float(_os.environ.get("SRFAR", SR))   # echelle du CHAMP LOINTAIN (R>25mm)
gmsh.model.mesh.field.setNumber(f3, "SizeMin", 0.002 * SRFAR)
gmsh.model.mesh.field.setNumber(f3, "SizeMax", 0.010 * SRFAR)
gmsh.model.mesh.field.setNumber(f3, "DistMin", 0.025)
gmsh.model.mesh.field.setNumber(f3, "DistMax", 0.100)
flist = [f1, f2, f3]
if RESTRICT and not ROCHE:
    rockv = names["rock"]
    rocks = [t for (d, t) in gmsh.model.getBoundary(
        [(3, t) for t in rockv], oriented=False, recursive=False)]
    rr = []
    for fld in (f1, f2):
        g = gmsh.model.mesh.field.add("Restrict")
        gmsh.model.mesh.field.setNumber(g, "InField", fld)
        gmsh.model.mesh.field.setNumbers(g, "VolumesList", [float(t) for t in rockv])
        gmsh.model.mesh.field.setNumbers(g, "SurfacesList", [float(t) for t in rocks])
        rr.append(g)
    flist = rr + [f3]
fmin = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", flist)
gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)

# tailles par corps (points OCC) : insert fin, bit et piston grossiers
def set_pts(tags, size):
    pts = gmsh.model.getBoundary([(3, t) for t in tags], recursive=True)
    gmsh.model.mesh.setSize([p for p in pts if p[0] == 0], size)

set_pts(names["rock"], 0.010 * SRFAR)
if not ROCHE:
    set_pts(names["bit"], 0.003 * s)
    set_pts(names["insert"], 0.0007 * s)
if not LEGER and not ROCHE:
    set_pts(names["piston"], 0.005 * s)
    set_pts(names["circlip"], 0.0015 * s)
    set_pts(names["plate"], 0.004 * s)

gmsh.option.setNumber("Mesh.RandomSeed", 1)
import os
if os.environ.get("OPTI", "0") == "1":
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.5)
gmsh.model.mesh.generate(3)
if os.environ.get("OPTI", "0") == "1":
    for _ in range(3):
        gmsh.model.mesh.optimize("Netgen")
        gmsh.model.mesh.optimize("")
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(out)
ntet = len(gmsh.model.mesh.getElementsByType(4)[0])
gmsh.finalize()
print("ecrit : %s  (%d tets, echelle %.2f, roche %.2f)" % (out, ntet, s, SR))
