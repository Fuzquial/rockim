#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_circle_mesh.py — plaque W x H percee d'un trou CIRCULAIRE R au centre,
# maillage gradue (hFine jusqu'a rFine de la paroi, hFar au loin).
#
#   python make_circle_mesh.py W H R hFine rFine hFar out.msh [seed]
#   ex :  python tunnel_edz/tools/make_circle_mesh.py 100 100 5 0.25 10 3 \
#               meshes/kirsch_r5.msh 1
#
# Sert au CONTROLE DE KIRSCH (V2) : c'est la seule geometrie pour laquelle la
# contrainte en paroi a une forme fermee, donc le seul controle quantitatif de
# la chaine pre-contrainte + relachement + rouleaux.
#
# Volontairement AUTONOME (ne touche pas tools/make_unstructured_mesh.py) : si
# le controle devient permanent, ce type se replie en une clause `circle` du
# generateur principal.
# ---------------------------------------------------------------------------
import sys

import gmsh


def main():
    if len(sys.argv) < 8:
        raise SystemExit(__doc__ or "usage: make_circle_mesh.py W H R hFine "
                                    "rFine hFar out.msh [seed]")
    W, H, R, hFine, rFine, hFar = map(float, sys.argv[1:7])
    out = sys.argv[7]
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("kirsch")
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", 6)          # frontal-Delaunay
    gmsh.option.setNumber("Mesh.Optimize", 1)

    cx, cy = 0.5 * W, 0.5 * H
    plate = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
    hole = gmsh.model.occ.addDisk(cx, cy, 0, R, R)
    gmsh.model.occ.cut([(2, plate)], [(2, hole)])
    gmsh.model.occ.synchronize()

    # courbes de la paroi = celles dont la boite englobante tient dans le trou
    wall = []
    for dim, tag in gmsh.model.getEntities(1):
        x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(dim, tag)
        if (x0 > cx - R - 1e-6 and x1 < cx + R + 1e-6
                and y0 > cy - R - 1e-6 and y1 < cy + R + 1e-6):
            wall.append(tag)
    if not wall:
        raise SystemExit("contour du trou introuvable")

    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", wall)
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 400)
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

    import numpy as np
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
    nWall = int(round(2.0 * np.pi * R / hFine))
    print(f"[mesh] {len(tri)} triangles, h inscrit min/med = "
          f"{hin.min():.4f}/{np.median(hin):.4f} m, "
          f"~{nWall} elements sur le pourtour du trou")
    gmsh.finalize()


if __name__ == "__main__":
    main()
