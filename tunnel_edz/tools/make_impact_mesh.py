#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_impact_mesh.py — bloc PLEIN gradue autour du point d'impact, avec
# exactement les memes reglages de maillage que l'etude tunnel : gmsh,
# Mesh.Algorithm = 5 (Delaunay, ISOTROPE), champ de distance + seuil.
#
#   python tunnel_edz/tools/make_impact_mesh.py W H hFine rFine hFar out.msh [seed]
#   ex :  python tunnel_edz/tools/make_impact_mesh.py 40 24 0.22 12 1.5 \
#               meshes/impact_block.msh 1
#
# La seule difference avec le maillage du tunnel : il n'y a pas de cavite, et
# le raffinement est centre sur le POINT D'IMPACT (milieu du bord superieur)
# au lieu de la paroi de l'excavation.
#
# ⚠️ NE PAS remettre Mesh.Algorithm = 6 (frontal-Delaunay, le defaut du
# generateur maison) : il pave en triangles quasi equilateraux des que la
# taille est localement constante — orientation des aretes periodique a 60 deg,
# rapport pic/creux 30, mediane du plus petit angle 60,0 deg. Les fissures
# n'ont alors que trois directions disponibles. Mesure du 2026-08-17, cf.
# tunnel_edz/tools/mesh_algo_sweep.py.
# ---------------------------------------------------------------------------
import sys

import gmsh
import numpy as np


def main():
    if len(sys.argv) < 7:
        raise SystemExit("usage: make_impact_mesh.py W H hFine rFine hFar "
                         "out.msh [seed]")
    W, H, hFine, rFine, hFar = map(float, sys.argv[1:6])
    out = sys.argv[6]
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("impact")
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", 5)          # Delaunay — cf. en-tete
    gmsh.model.occ.addRectangle(0, 0, 0, W, H)
    gmsh.model.occ.synchronize()

    # raffinement autour du point d'impact : milieu du bord superieur
    px = gmsh.model.occ.addPoint(0.5 * W, H, 0.0)
    gmsh.model.occ.synchronize()
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "PointsList", [px])
    fth = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(fth, "InField", fd)
    gmsh.model.mesh.field.setNumber(fth, "SizeMin", hFine)
    gmsh.model.mesh.field.setNumber(fth, "SizeMax", hFar)
    gmsh.model.mesh.field.setNumber(fth, "DistMin", rFine)
    gmsh.model.mesh.field.setNumber(fth, "DistMax", rFine + 0.8 * rFine)
    gmsh.model.mesh.field.setAsBackgroundMesh(fth)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", hFine)
    gmsh.option.setNumber("Mesh.MeshSizeMax", hFar)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(out)

    _, _, conn = gmsh.model.mesh.getElements(2)
    tri = np.array(conn[0], dtype=int).reshape(-1, 3)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    idx = {t: i for i, t in enumerate(tags)}
    P = np.array(xyz).reshape(-1, 3)[[idx[t] for t in tri.flatten()], :2]
    P = P.reshape(-1, 3, 2)
    e = np.stack([np.linalg.norm(P[:, 2] - P[:, 1], axis=1),
                  np.linalg.norm(P[:, 0] - P[:, 2], axis=1),
                  np.linalg.norm(P[:, 1] - P[:, 0], axis=1)], axis=1)
    u, v = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]
    Ar = 0.5 * np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
    hin = 4.0 * Ar / e.sum(axis=1)
    lcz = 10e9 * 20.0 / (0.6e6 ** 2)
    print(f"[impact] bloc {W:g} x {H:g} m, impact en ({0.5*W:g}, {H:g}), "
          f"Mesh.Algorithm = 5 (Delaunay, isotrope)")
    print(f"[mesh] {len(tri)} triangles, h inscrit min/med = "
          f"{hin.min():.4f}/{np.median(hin):.4f} m")
    print(f"[mesh] l_cz = {lcz:.3f} m -> dx admissible <= {0.5*lcz:.3f} m ; "
          f"hFine = {hFine:.3f} m [{'OK' if hFine <= 0.5*lcz else 'TROP GROSSIER'}]")
    gmsh.finalize()


if __name__ == "__main__":
    main()
