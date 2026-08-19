#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_impact_mesh3d.py — bloc 3D plein gradue autour du point d'impact
# (centre de la face superieure). Transposition directe de
# make_impact_mesh.py : meme champ Distance + Threshold, meme logique.
#
#   python tunnel_edz/tools/make_impact_mesh3d.py W D H hFine rFine hFar out.msh [seed]
#   ex :  python tunnel_edz/tools/make_impact_mesh3d.py 0.04 0.04 0.025 \
#               0.0006 0.006 0.0025 meshes/indent3d_block.msh 1
#
# Mesh.Algorithm = 5 pour la peau (Delaunay isotrope, cf. le piege du
# frontal-Delaunay du 2026-08-17) et Mesh.Algorithm3D = 1 (Delaunay) pour le
# volume, qui est le defaut de gmsh.
# ---------------------------------------------------------------------------
import sys

import gmsh
import numpy as np


def main():
    if len(sys.argv) < 8:
        raise SystemExit("usage: make_impact_mesh3d.py W D H hFine rFine hFar "
                         "out.msh [seed]")
    W, D, H, hFine, rFine, hFar = map(float, sys.argv[1:7])
    out = sys.argv[7]
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("indent3d")
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", 5)       # peau : Delaunay isotrope
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)     # volume : Delaunay
    gmsh.model.occ.addBox(0, 0, 0, W, D, H)
    gmsh.model.occ.synchronize()

    px = gmsh.model.occ.addPoint(0.5 * W, 0.5 * D, H)
    gmsh.model.occ.synchronize()
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "PointsList", [px])
    fth = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(fth, "InField", fd)
    gmsh.model.mesh.field.setNumber(fth, "SizeMin", hFine)
    gmsh.model.mesh.field.setNumber(fth, "SizeMax", hFar)
    gmsh.model.mesh.field.setNumber(fth, "DistMin", rFine)
    gmsh.model.mesh.field.setNumber(fth, "DistMax", rFine + 1.2 * rFine)
    gmsh.model.mesh.field.setAsBackgroundMesh(fth)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", hFine)
    gmsh.option.setNumber("Mesh.MeshSizeMax", hFar)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(out)

    tp, tg, cn = gmsh.model.mesh.getElements(3)
    tet = np.array(cn[0], dtype=int).reshape(-1, 4)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    idx = {t: i for i, t in enumerate(tags)}
    P = np.array(xyz).reshape(-1, 3)[[idx[t] for t in tet.flatten()], :]
    P = P.reshape(-1, 4, 3)
    a, b, c = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]
    V = np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0
    # aire totale des 4 faces -> diametre inscrit 6V/A (c'est lui qui fixe dt)
    A = np.zeros(len(tet))
    for i, j, k in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
        A += 0.5 * np.linalg.norm(np.cross(P[:, j] - P[:, i], P[:, k] - P[:, i]),
                                  axis=1)
    hin = 6.0 * V / A
    print(f"[indent3d] bloc {W*1e3:g} x {D*1e3:g} x {H*1e3:g} mm, impact au "
          f"centre de la face superieure")
    print(f"[mesh] {len(tet)} tetraedres, h inscrit min/med = "
          f"{hin.min()*1e3:.4f}/{np.median(hin)*1e3:.4f} mm")
    gmsh.finalize()


if __name__ == "__main__":
    main()
