#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_impact2d_mesh.py — eprouvette 2D d'impact a insert unique, avec un
# maillage GRADUE RADIALEMENT sous le point d'impact.
#
#   python tools/make_impact2d_mesh.py W H hFine rFine hFar out.msh [seed]
#   ex :  python tools/make_impact2d_mesh.py 0.250 0.150 0.001 0.0125 0.010 \
#              meshes/impact2d_grad.msh 1
#
# POURQUOI CE MAILLAGE. Les premiers impacts 2D ont ete faits sur un bloc
# UNIFORME de 120 x 60 mm. Mesure du 2026-09-03 : AUCUN joint ne rompt sous
# le cratere (0 sur 0-20 mm de profondeur) et TOUTES les ruptures se trouvent
# entre 32 et 59 mm, sur une colonne large de 2 mm centree sur l'axe — la
# signature d'un ECAILLAGE par reflexion de l'onde sur le fond, pas d'une
# fissuration d'impact. Le bloc etait trop mince : l'onde revient en ~25 us,
# pendant le chargement.
#
# LA PARADE, telle que l'article la pose (Yang et al. 2026, IJRMMS 206,
# leur fig. 6 et leur §3.1) : une eprouvette LARGE, maillee FIN seulement
# sous l'impact et de plus en plus grossiere en s'en eloignant.
#   * leur fig. 6 : 1 mm dans la boule R 12,5 mm, 2 mm jusqu'a R 25 mm,
#     10 mm au bord ;
#   * leur §3.1 : « a minimum ratio of approximately 1:6 between the damage
#     extent and the specimen size is required to ensure that boundary
#     effects do not significantly influence crack propagation behaviour ».
#     Leurs fissures vont a 21 mm pour une eprouvette de 250 mm — soit 1:12.
#
# Le champ de taille est un Threshold sur la DISTANCE au point d'impact :
# hFine jusqu'a rFine, puis croissance lineaire jusqu'a hFar a 8 rFine.
# C'est la demi-couronne raffinee demandee, sans discontinuite de taille.
#
# Le point d'impact est le MILIEU DU BORD SUPERIEUR (x = W/2, y = H), la ou
# le scenario `percussion` de rockim place l'outil par defaut.
#
# Sortie : Gmsh MSH 2.2 ASCII, triangles type 2 — le format `mesh = file` 2D.
# ---------------------------------------------------------------------------
import sys
import gmsh


def main():
    if len(sys.argv) < 7:
        raise SystemExit("usage: make_impact2d_mesh.py W H hFine rFine hFar "
                         "out.msh [seed]")
    W, H, hFine, rFine, hFar = map(float, sys.argv[1:6])
    out = sys.argv[6]
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("impact2d")

    gmsh.model.occ.addRectangle(0, 0, 0, W, H)
    gmsh.model.occ.synchronize()

    # --- champ de taille : distance au point d'impact ----------------------
    px = gmsh.model.occ.addPoint(W / 2.0, H, 0.0, hFine)
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
    # le champ SEUL commande la taille : sinon Gmsh melange courbures et bords
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    # Delaunay + optimisation Netgen : le diametre inscrit minimal pilote la
    # CFL de rockim, l'optimisation le remonte (meme choix que les autres
    # generateurs de la maison ; le frontal est BANNI, il aligne les aretes).
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(out)

    _, _, tri = gmsh.model.mesh.getElements(2)
    n = len(tri[0]) // 3 if tri else 0
    print(f"[impact2d] eprouvette {W*1e3:.0f} x {H*1e3:.0f} mm, "
          f"impact en (x = {W/2*1e3:.0f}, y = {H*1e3:.0f}) mm")
    print(f"[impact2d] maille {hFine*1e3:.2f} mm jusqu'a R {rFine*1e3:.1f} mm, "
          f"{hFar*1e3:.1f} mm au-dela de R {8*rFine*1e3:.0f} mm")
    print(f"[mesh] {n} triangles -> {out}")
    gmsh.finalize()


if __name__ == "__main__":
    main()
