#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_cut3d_mesh.py — eprouvette de COUPE 3D au cutter PDC, avec ENTAILLE de
# depart et zone fine en COULOIR le long de la course du cutter.
# Reproduction de Heilman et al., ARMA 24-0238 (Utah FORGE, 40 x 30 x 20 mm).
#
#   python tools/make_cut3d_mesh.py W D H cutDepth jeu notchLen \
#          hFine xFineMin xFineMax bandY bandH hFar out.msh [seed]
#   ex :  python tools/make_cut3d_mesh.py 0.040 0.030 0.020 0.001016 0.000284 \
#          0.003 0.0005 0.001 0.010 0.006 0.003 0.004 meshes/cut3d_h50.msh 1
#
# L ENTAILLE. L article demarre le cutter dans une entaille dont la profondeur
# egale la passe. On y ajoute un JEU (`jeu`) : la lecon 2D de
# configs/cut2d_v3ter_rake.cfg (l. 110-123) est qu une entaille EXACTEMENT a la
# profondeur de passe fait rompre 51 joints AVANT que le cutter n atteigne la
# face verticale, contre 11 avec 0,284 mm de jeu — le plancher de l entaille et
# l arete du cutter sont sinon a la meme cote, a un arrondi pres.
# Profondeur d entaille = cutDepth + jeu. Longueur = notchLen depuis x = 0.
#
# LE COULOIR. La zone fine n est PAS une boule (make_impact3d_mesh.py) mais une
# BOITE [xFineMin, xFineMax] x [D/2 - bandY, D/2 + bandY] x [H - bandH, H] :
# la course reelle du cutter plus la profondeur fissuree. C est ce qui rend le
# 3D abordable — a hFar = 15 mm le champ lointain payait les trois quarts d un
# impact 3D pour de la roche ou rien ne se passe (mesure du 2026-09-03).
# Gradation par le champ Box de gmsh (Thickness = transition lineaire).
#
# SORTIE. MSH 2.2 ASCII, tetraedres type 4 (`mesh = file`). Le script imprime
# le compte de tetraedres et hmin = min(6V/A) — la longueur qui borne le pas de
# temps du solveur — pour dimensionner AVANT de lancer.
# ---------------------------------------------------------------------------
import sys

import gmsh
import numpy as np


def hmin_of(nodes, tets):
    P = nodes[tets]
    v = P[:, 1:] - P[:, [0]]
    V = np.abs(np.einsum("ij,ij->i", v[:, 0], np.cross(v[:, 1], v[:, 2]))) / 6.0
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    A = sum(0.5 * np.linalg.norm(np.cross(P[:, b] - P[:, a], P[:, c] - P[:, a]),
                                 axis=1) for a, b, c in faces)
    h = 6.0 * V / A
    return h.min(), np.median(h)


def main():
    if len(sys.argv) < 14:
        raise SystemExit("usage: make_cut3d_mesh.py W D H cutDepth jeu notchLen "
                         "hFine xFineMin xFineMax bandY bandH hFar out.msh [seed]")
    (W, D, H, cutDepth, jeu, notchLen, hFine, xFineMin, xFineMax,
     bandY, bandH, hFar) = map(float, sys.argv[1:13])
    out = sys.argv[13]
    seed = int(sys.argv[14]) if len(sys.argv) > 14 else 1
    depth = cutDepth + jeu

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cut3d")
    gmsh.option.setNumber("Mesh.RandomSeed", seed)

    blk = gmsh.model.occ.addBox(0, 0, 0, W, D, H)
    if notchLen > 0.0:
        ntc = gmsh.model.occ.addBox(-1e-3, -1e-3, H - depth, notchLen + 1e-3,
                                    D + 2e-3, depth + 1e-3)
        gmsh.model.occ.cut([(3, blk)], [(3, ntc)])
    gmsh.model.occ.synchronize()

    yc = 0.5 * D
    fb = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(fb, "VIn", hFine)
    gmsh.model.mesh.field.setNumber(fb, "VOut", hFar)
    gmsh.model.mesh.field.setNumber(fb, "XMin", xFineMin)
    gmsh.model.mesh.field.setNumber(fb, "XMax", xFineMax)
    gmsh.model.mesh.field.setNumber(fb, "YMin", yc - bandY)
    gmsh.model.mesh.field.setNumber(fb, "YMax", yc + bandY)
    gmsh.model.mesh.field.setNumber(fb, "ZMin", H - bandH)
    gmsh.model.mesh.field.setNumber(fb, "ZMax", H + 1.0)
    gmsh.model.mesh.field.setNumber(fb, "Thickness", 3.0 * bandH)
    gmsh.model.mesh.field.setAsBackgroundMesh(fb)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)          # Delaunay
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(out)

    tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(coords).reshape(-1, 3)
    idx = {t: i for i, t in enumerate(tags)}
    _, _, te = gmsh.model.mesh.getElements(3)
    tets = np.array([idx[t] for t in te[0]]).reshape(-1, 4) if te else np.zeros((0, 4), int)
    hmin, hmed = hmin_of(nodes, tets) if len(tets) else (0.0, 0.0)
    print(f"[cut3d] eprouvette {W*1e3:.0f} x {D*1e3:.0f} x {H*1e3:.0f} mm ; "
          f"entaille {notchLen*1e3:.1f} mm de long x {depth*1e3:.3f} mm de profond "
          f"(passe {cutDepth*1e3:.3f} + jeu {jeu*1e3:.3f})")
    print(f"[cut3d] couloir fin x in [{xFineMin*1e3:.1f}, {xFineMax*1e3:.1f}] mm, "
          f"+-{bandY*1e3:.1f} mm en y, {bandH*1e3:.1f} mm sous la surface, "
          f"hFine {hFine*1e3:.2f} mm -> hFar {hFar*1e3:.2f} mm")
    print(f"[mesh] {len(tets)} tetraedres, hmin (6V/A) {hmin*1e3:.4f} mm, "
          f"h median {hmed*1e3:.4f} mm -> {out}")
    gmsh.finalize()


if __name__ == "__main__":
    main()
