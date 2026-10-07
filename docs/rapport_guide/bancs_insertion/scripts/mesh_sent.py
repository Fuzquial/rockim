# -*- coding: utf-8 -*-
"""mesh_sent.py — plaque a entaille laterale (SENT) pour le banc 3, gmsh, MSH 2.2.
usage : python mesh_sent.py W H a e h out.msh   (m ; entaille rectangulaire depuis x = 0,
a mi-hauteur, longueur a, ouverture e ; maillage Delaunay uniforme de taille h)"""
import sys, gmsh
W, H, a, e, h = map(float, sys.argv[1:6]); out = sys.argv[6]
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("sent")
y0, y1 = H / 2 - e / 2, H / 2 + e / 2
P = [(0, 0), (W, 0), (W, H), (0, H), (0, y1), (a, y1), (a, y0), (0, y0)]
pt = [gmsh.model.geo.addPoint(x, y, 0, h) for x, y in P]
ln = [gmsh.model.geo.addLine(pt[i], pt[(i + 1) % len(pt)]) for i in range(len(pt))]
cl = gmsh.model.geo.addCurveLoop(ln); s = gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.synchronize(); gmsh.model.addPhysicalGroup(2, [s], 1)
gmsh.option.setNumber("Mesh.Algorithm", 5); gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.mesh.generate(2); gmsh.model.mesh.optimize("Laplace2D")
gmsh.write(out); gmsh.finalize(); print("ecrit", out)
