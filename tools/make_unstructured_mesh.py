#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_unstructured_mesh.py — genere un maillage simplexe NON STRUCTURE
# uniforme (le maillage "a la Yan et al. 2023") pour rockim `mesh = file`.
#
#   python3 tools/make_unstructured_mesh.py box3d W D H h out.msh [seed]
#   python3 tools/make_unstructured_mesh.py box2d W H   h out.msh [seed]
#   python3 tools/make_unstructured_mesh.py bench1 W D H R gap h hIns out.msh [seed]
#
# bench1 (V1) : DEUX corps — bloc W x D x H (volume physique "rock") et
# INSERT spherique de rayon R centre au-dessus (volume physique "insert"),
# separes de `gap`, mailles a h (roche) et hIns (insert). C'est le premier
# jalon de la trajectoire Solidity : l'impact insert-roche en maille.
#
# Sortie : Gmsh MSH 2.2 ASCII (tets type 4 en 3D, triangles type 2 en 2D).
# Gmsh est LE mailleur de la pratique FDEM (Akantu, OpenFDEM, litterature) ;
# pip install gmsh. Algorithmes : Delaunay + optimisation Netgen (qualite —
# le min du diametre inscrit pilote la CFL de rockim). Le seed perturbe le
# point de depart du frontal 2D pour varier la realisation.
# ---------------------------------------------------------------------------
import sys
import gmsh


def main():
    kind = sys.argv[1]
    if kind == "box3d":
        W, D, H, h = map(float, sys.argv[2:6])
        out = sys.argv[6]
        seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    elif kind == "box2d":
        W, H, h = map(float, sys.argv[2:5])
        out = sys.argv[5]
        seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    elif kind == "bench1":
        W, D, H, R, gap, h, hIns = map(float, sys.argv[2:9])
        out = sys.argv[9]
        seed = int(sys.argv[10]) if len(sys.argv) > 10 else 1
    elif kind == "tunnel":
        # plaque W x H percee d'un trou circulaire R au centre (2D) —
        # cavite pressurisee par confineFaces = bore
        W, H, R, h = map(float, sys.argv[2:6])
        out = sys.argv[6]
        seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    else:
        raise SystemExit("usage: make_unstructured_mesh.py box3d W D H h out.msh [seed]\n"
                         "       make_unstructured_mesh.py box2d W H h out.msh [seed]\n"
                         "       make_unstructured_mesh.py bench1 W D H R gap h hIns out.msh [seed]\n"
                         "       make_unstructured_mesh.py tunnel W H R h out.msh [seed]")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("rockim")
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", 6)        # frontal-Delaunay 2D
    gmsh.option.setNumber("Mesh.Optimize", 1)
    if kind == "box3d":
        gmsh.model.occ.addBox(0, 0, 0, W, D, H)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay 3D
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.model.mesh.generate(3)
    elif kind == "bench1":
        rock = gmsh.model.occ.addBox(0, 0, 0, W, D, H)
        ins = gmsh.model.occ.addSphere(0.5 * W, 0.5 * D, H + gap + R, R)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [rock], name="rock")
        gmsh.model.addPhysicalGroup(3, [ins], name="insert")
        # taille locale : fine sur l'insert, h partout ailleurs
        gmsh.option.setNumber("Mesh.MeshSizeMin", min(h, hIns))
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        insPts = gmsh.model.getBoundary([(3, ins)], recursive=True)
        for dim, tag in insPts:
            if dim == 0:
                gmsh.model.mesh.setSize([(0, tag)], hIns)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.model.mesh.generate(3)
    elif kind == "tunnel":
        plate = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
        hole = gmsh.model.occ.addDisk(0.5 * W, 0.5 * H, 0, R, R)
        gmsh.model.occ.cut([(2, plate)], [(2, hole)])
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
    else:
        gmsh.model.occ.addRectangle(0, 0, 0, W, H)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(out)
    # bilan qualite : diametre inscrit min/med (ce qui pilote le dt rockim)
    import numpy as np
    if kind in ("box3d", "bench1"):
        _, _, conn = gmsh.model.mesh.getElements(3)
        tets = np.array(conn[0], dtype=int).reshape(-1, 4)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        idx = {t: i for i, t in enumerate(tags)}
        P = np.array(xyz).reshape(-1, 3)[[idx[t] for t in tets.flatten()]]
        P = P.reshape(-1, 4, 3)
        V = np.abs(np.einsum("ij,ij->i",
                             np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                             P[:, 3] - P[:, 0])) / 6.0
        F = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
        A = sum(0.5 * np.linalg.norm(
                np.cross(P[:, f[1]] - P[:, f[0]], P[:, f[2]] - P[:, f[0]]),
                axis=1) for f in F)
        hin = 6.0 * V / A
        print(f"[mesh] {len(tets)} tets, h inscrit min/med/max = "
              f"{hin.min()*1e3:.3f}/{np.median(hin)*1e3:.3f}/{hin.max()*1e3:.3f} mm")
    gmsh.finalize()


if __name__ == "__main__":
    main()
