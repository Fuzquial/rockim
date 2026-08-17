# -*- coding: utf-8 -*-
"""Maillage 2D GRADUE pour la coupe : fin sous la surface (zone de coupe),
grossier en profondeur. Impossible en mode voronoi — la tessellation y est
uniforme par construction — d'ou ce maillage libre charge par mesh = file.

  python make_cut_mesh.py [W_mm H_mm hfin_mm hgros_mm bande_mm sortie.msh]

Le champ de taille est une rampe lineaire en y : h = hfin sur la bande
superieure, puis croissance jusqu'a hgros au fond.
"""
import os, sys
import gmsh
import numpy as np

W = float(sys.argv[1]) if len(sys.argv) > 1 else 30.0
H = float(sys.argv[2]) if len(sys.argv) > 2 else 15.0
HFIN = float(sys.argv[3]) if len(sys.argv) > 3 else 0.25
HGROS = float(sys.argv[4]) if len(sys.argv) > 4 else 1.6
BANDE = float(sys.argv[5]) if len(sys.argv) > 5 else 3.0
OUT = sys.argv[6] if len(sys.argv) > 6 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "meshes",
    "cut_graded.msh")

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("cut")
s = gmsh.model.occ.addRectangle(0, 0, 0, W * 1e-3, H * 1e-3)
gmsh.model.occ.synchronize()

# --- champ de taille : rampe en y ---------------------------------------
yTop, yBand = H * 1e-3, (H - BANDE) * 1e-3
f = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(
    f, "F", "%.6g + (%.6g - %.6g) * max(0, min(1, (%.6g - y) / %.6g))"
    % (HFIN * 1e-3, HGROS * 1e-3, HFIN * 1e-3, yBand, yBand))
gmsh.model.mesh.field.setAsBackgroundMesh(f)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.Algorithm", 6)          # frontal-Delaunay
gmsh.option.setNumber("Mesh.RandomSeed", 1)
gmsh.model.mesh.generate(2)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(OUT)

# --- bilan qualite ------------------------------------------------------
_, _, conn = gmsh.model.mesh.getElements(2)
tri = np.array(conn[0], dtype=int).reshape(-1, 3)
tags, xyz, _ = gmsh.model.mesh.getNodes()
idx = {t: i for i, t in enumerate(tags)}
P = np.array(xyz).reshape(-1, 3)[[idx[t] for t in tri.flatten()]][:, :2]
P = P.reshape(-1, 3, 2)
a = np.linalg.norm(P[:, 1] - P[:, 0], axis=1)
b = np.linalg.norm(P[:, 2] - P[:, 1], axis=1)
c = np.linalg.norm(P[:, 0] - P[:, 2], axis=1)
sp = 0.5 * (a + b + c)
A = np.sqrt(np.maximum(sp * (sp - a) * (sp - b) * (sp - c), 0))
rin = A / sp                                        # rayon inscrit
cy = P[:, :, 1].mean(axis=1)
top = cy > yBand
print("[mesh] %d triangles | %.0f %% dans la bande superieure (%g mm)"
      % (len(tri), 100 * top.mean(), BANDE))
print("[mesh] diametre inscrit : bande %.3f mm, fond %.3f mm, min global %.3f mm"
      % (2e3 * rin[top].mean(), 2e3 * rin[~top].mean(), 2e3 * rin.min()))
print("[mesh] ecrit", os.path.normpath(OUT))
gmsh.finalize()
