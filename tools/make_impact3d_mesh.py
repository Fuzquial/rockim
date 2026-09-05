#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_impact3d_mesh.py — eprouvette 3D d'impact a insert unique, maillage
# GRADUE RADIALEMENT sous le point d'impact. Version 3D de
# make_impact2d_mesh.py, meme champ de taille.
#
#   python tools/make_impact3d_mesh.py W D H hFine rFine hFar out.msh [seed]
#   ex :  python tools/make_impact3d_mesh.py 0.150 0.150 0.100 0.002 0.0125 \
#              0.015 meshes/impact3d_gros.msh 1
#
# Le point d'impact est le CENTRE DE LA FACE SUPERIEURE (W/2, D/2, H), la ou
# le scenario `percussion` de rockim place l'outil.
#
# Champ de taille : Threshold sur la distance au point d'impact — hFine
# jusqu'a rFine, croissance lineaire jusqu'a hFar a 8 rFine. C'est la
# demi-boule raffinee de la fig. 6 de Yang et al. 2026 (1 mm dans R 12,5 mm,
# 2 mm jusqu'a R 25, 10 mm au bord).
#
# EN 3D LE COMPTE EXPLOSE : a taille egale, un tetraedre occupe ~1/8,5 du cube
# de son arete, alors qu'un triangle occupe ~1/2,3 du carre. Passer de 2D a 3D
# a maillage identique multiplie donc les elements par ~(L/h), et le pas de
# temps baisse encore. D'ou la gradation, indispensable ici, et la version
# GROSSIERE (hFine = 2 mm) pour valider le montage avant de raffiner.
#
# Sortie : Gmsh MSH 2.2 ASCII, tetraedres type 4 — le format `mesh = file` 3D.
# ---------------------------------------------------------------------------
import sys
import gmsh


def main():
    if len(sys.argv) < 8:
        raise SystemExit("usage: make_impact3d_mesh.py W D H hFine rFine hFar "
                         "out.msh [seed]")
    W, D, H, hFine, rFine, hFar = map(float, sys.argv[1:7])
    out = sys.argv[7]
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("impact3d")
    gmsh.model.occ.addBox(0, 0, 0, W, D, H)
    gmsh.model.occ.synchronize()

    px = gmsh.model.occ.addPoint(W / 2.0, D / 2.0, H, hFine)
    gmsh.model.occ.synchronize()
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "PointsList", [px])
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", hFine)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", hFar)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", rFine)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", 8.0 * rFine)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.option.setNumber("Mesh.Algorithm3D", 1)     # Delaunay
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(out)

    _, _, te = gmsh.model.mesh.getElements(3)
    n = len(te[0]) // 4 if te else 0
    print(f"[impact3d] eprouvette {W*1e3:.0f} x {D*1e3:.0f} x {H*1e3:.0f} mm, "
          f"impact en ({W/2*1e3:.0f}, {D/2*1e3:.0f}, {H*1e3:.0f}) mm")
    print(f"[impact3d] maille {hFine*1e3:.2f} mm jusqu'a R {rFine*1e3:.1f} mm, "
          f"{hFar*1e3:.1f} mm au-dela de R {8*rFine*1e3:.0f} mm")
    print(f"[mesh] {n} tetraedres -> {out}")
    gmsh.finalize()


if __name__ == "__main__":
    main()
