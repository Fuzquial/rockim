#!/usr/bin/env python3
# Maillage du banc P1.4 (problème de Lamb : force ponctuelle verticale en
# échelon sur la surface d'un demi-espace). Bloc 240 x 240 x 120 mm, source au
# centre de la face supérieure, maillage tétraédrique NON STRUCTURÉ gradué :
# h = hFin dans la demi-sphère de rayon rFin autour de la source, croissance
# linéaire jusqu'à hLoin à rLoin. La source et les stations de surface sont des
# points INCLUS dans la face (sommets exacts : aucune erreur de placement).
# Groupes physiques : volume "solid", faces xmin xmax ymin ymax bottom top,
# points src, r20, r30, r40, r50 (axe x) et d40 (diagonale, contrôle d'isotropie).
#
#   python3 docs/rapport_guide/simulations_a_lancer/configs/nouveaux_bancs/lamb_pekeris/gen_mesh.py [hFin_mm] [seed]
# Sortie : docs/rapport_guide/simulations_a_lancer/meshes/nb_lamb_h<hFin>.msh
# Recette Gmsh (Delaunay 3D + optimisation Netgen, MSH 2.2) : celle de
# tools/make_unstructured_mesh.py.
import math
import os
import sys

import gmsh

W, D, H = 0.24, 0.24, 0.12
R_FIN, R_LOIN, H_LOIN = 0.058, 0.110, 0.008
STATIONS = {"r20": (0.020, 0.0), "r30": (0.030, 0.0), "r40": (0.040, 0.0),
            "r50": (0.050, 0.0), "d40": (0.040 / math.sqrt(2), 0.040 / math.sqrt(2))}


def main():
    h_fin = float(sys.argv[1]) * 1e-3 if len(sys.argv) > 1 else 1.5e-3
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 6))
    out = os.path.join(root, "docs", "rapport_guide", "simulations_a_lancer", "meshes",
                       f"nb_lamb_h{h_fin * 1e3:g}".replace(".", "p") + ".msh")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.model.add("lamb")
    occ = gmsh.model.occ
    box = occ.addBox(0, 0, 0, W, D, H)
    cx, cy = W / 2, D / 2
    pts = {"src": occ.addPoint(cx, cy, H)}
    for g, (dx, dy) in STATIONS.items():
        pts[g] = occ.addPoint(cx + dx, cy + dy, H)
    occ.synchronize()
    # face supérieure
    top = [s for (dim, s) in gmsh.model.getEntities(2)
           if abs(gmsh.model.occ.getCenterOfMass(2, s)[2] - H) < 1e-9][0]
    gmsh.model.mesh.embed(0, list(pts.values()), 2, top)
    # champ de taille : distance à la source
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "PointsList", [pts["src"]])
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", h_fin)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", H_LOIN)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", R_FIN)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", R_LOIN)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)
    for k in ("MeshSizeExtendFromBoundary", "MeshSizeFromPoints", "MeshSizeFromCurvature"):
        gmsh.option.setNumber(f"Mesh.{k}", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    # groupes physiques
    vol = gmsh.model.getEntities(3)[0][1]
    gmsh.model.addPhysicalGroup(3, [vol], name="solid")
    tol = 1e-9
    for (dim, s) in gmsh.model.getEntities(2):
        x, y, z = gmsh.model.occ.getCenterOfMass(2, s)
        name = ("xmin" if abs(x) < tol else "xmax" if abs(x - W) < tol else
                "ymin" if abs(y) < tol else "ymax" if abs(y - D) < tol else
                "bottom" if abs(z) < tol else "top")
        gmsh.model.addPhysicalGroup(2, [s], name=name)
    for g, p in pts.items():
        gmsh.model.addPhysicalGroup(0, [p], name=g)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.write(out)
    _, etags, _ = gmsh.model.mesh.getElements(3)
    print(f"{out} : {len(etags[0])} tetraedres (hFin = {h_fin * 1e3:g} mm)")
    gmsh.finalize()


if __name__ == "__main__":
    main()
