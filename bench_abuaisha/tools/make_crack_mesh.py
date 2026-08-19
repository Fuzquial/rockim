#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_crack_mesh.py — bloc W x H avec une DISCONTINUITE PLANE de longueur L au
# centre, d'apres AbuAisha, Eaton, Priest & Wong, J. Petrol. Sci. Eng. 154
# (2017), figure A.20 : « A discontinuity subjected to uniform fluid pressure
# of 12 MPa and embedded in high strength linear elastic medium ».
#
#   python bench_abuaisha/tools/make_crack_mesh.py W H L hFine hFar out.msh [seed]
#   ex :  python bench_abuaisha/tools/make_crack_mesh.py 8.0 8.0 1.5 0.003 0.3 \
#               meshes/parker_crack.msh 1
#
# LA DIFFICULTE, ET COMMENT ELLE EST RESOLUE. Une discontinuite d'epaisseur
# NULLE ne peut pas etre obtenue en decoupant une fente : l'ouverture attendue
# vaut 0,065 mm pour une maille de 3 mm, donc toute fente maillable serait plus
# epaisse que le resultat cherche. On maille donc le bloc PLEIN avec la ligne
# de fissure imposee comme contrainte interne (`embed`), puis on DEDOUBLE les
# noeuds de cette ligne avec le greffon `Crack` de gmsh. Les deux levres
# deviennent alors des faces EXTERIEURES du maillage, geometriquement
# confondues a t = 0, et rockim n'y insere aucun joint puisqu'elles ne
# partagent plus de noeud.
#
# CONSEQUENCE COTE ROCKIM : les levres etant des faces exterieures situees a
# moins de boreSelectR du centre, `confineFaces = bore` les selectionne — et
# elles seules, le bord exterieur du bloc etant a W/2 de la. C'est ce qui
# permet d'appliquer la pression uniforme de l'article sans une ligne de code
# nouvelle.
#
# GRADATION : fine le long de la fissure, relachee au loin. L'article annonce
# « element size range goes from 0.003 to 0.3 m ».
# ---------------------------------------------------------------------------
import sys

import gmsh
import numpy as np


def main():
    if len(sys.argv) < 7:
        raise SystemExit("usage: make_crack_mesh.py W H L hFine hFar out.msh "
                         "[seed]")
    W, H, L, hFine, hFar = map(float, sys.argv[1:6])
    out = sys.argv[6]
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1

    cx, cy = 0.5 * W, 0.5 * H
    x0, x1 = cx - 0.5 * L, cx + 0.5 * L

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("parker")
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", 5)      # Delaunay isotrope

    surf = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
    pA = gmsh.model.occ.addPoint(x0, cy, 0)
    pB = gmsh.model.occ.addPoint(x1, cy, 0)
    crk = gmsh.model.occ.addLine(pA, pB)
    gmsh.model.occ.synchronize()
    # la ligne devient une contrainte INTERNE du maillage : les aretes des
    # triangles s'alignent dessus, condition necessaire au dedoublement
    gmsh.model.mesh.embed(1, [crk], 2, surf)

    pgCrack = gmsh.model.addPhysicalGroup(1, [crk])
    gmsh.model.setPhysicalName(1, pgCrack, "crack")
    pgSurf = gmsh.model.addPhysicalGroup(2, [surf])
    gmsh.model.setPhysicalName(2, pgSurf, "rock")

    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", [crk])
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 800)
    fth = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(fth, "InField", fd)
    gmsh.model.mesh.field.setNumber(fth, "SizeMin", hFine)
    gmsh.model.mesh.field.setNumber(fth, "SizeMax", hFar)
    gmsh.model.mesh.field.setNumber(fth, "DistMin", 0.10 * L)
    gmsh.model.mesh.field.setNumber(fth, "DistMax", 2.0 * L)
    gmsh.model.mesh.field.setAsBackgroundMesh(fth)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", hFine)
    gmsh.option.setNumber("Mesh.MeshSizeMax", hFar)
    gmsh.model.mesh.generate(2)

    nBefore = len(gmsh.model.mesh.getNodes()[0])
    # LE dedoublement : les noeuds du groupe physique `crack` sont dupliques,
    # les triangles d'un cote pointent vers les copies. Les pointes restent
    # soudees (pas d'OpenBoundaryPhysicalGroup) : la fissure est interne au
    # bloc, comme dans l'article.
    gmsh.plugin.setNumber("Crack", "Dimension", 1)
    gmsh.plugin.setNumber("Crack", "PhysicalGroup", pgCrack)
    gmsh.plugin.setNumber("Crack", "NormalX", 0.0)
    gmsh.plugin.setNumber("Crack", "NormalY", 0.0)
    gmsh.plugin.setNumber("Crack", "NormalZ", 1.0)
    gmsh.plugin.run("Crack")
    nAfter = len(gmsh.model.mesh.getNodes()[0])

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
    print(f"[parker] bloc {W:g} x {H:g} m, fissure L = {L:g} m centree en "
          f"({cx:g}, {cy:g}), demi-longueur l = {0.5*L:g} m")
    print(f"[crack]  noeuds {nBefore} -> {nAfter} "
          f"(+{nAfter-nBefore} dedoubles sur les levres)")
    print(f"[mesh]   {len(tri)} triangles, h inscrit min/med = "
          f"{hin.min()*1e3:.3f}/{np.median(hin)*1e3:.3f} mm")
    if nAfter == nBefore:
        print("[ATTENTION] aucun noeud dedouble : le greffon Crack n'a pas "
              "agi, les levres ne seront PAS des faces exterieures.")
    gmsh.finalize()


if __name__ == "__main__":
    main()
