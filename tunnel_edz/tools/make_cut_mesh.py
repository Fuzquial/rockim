#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_cut_mesh.py — eprouvette de coupe 2D AVEC ENTAILLE DE DEPART, d'apres
# Heilman et al., ARMA 24-0238 (LANL/HOSS) : « we include a notch in the
# specimen where the cutter begins, the notch depth is equal to the depth of
# cut ».
#
#   python tunnel_edz/tools/make_cut_mesh.py W H depth notchLen hFine bandH hFar out.msh [seed]
#   ex :  python tunnel_edz/tools/make_cut_mesh.py 0.040 0.020 0.001016 0.003 \
#               0.00025 0.004 0.001 meshes/cut_notch.msh 1
#
# POURQUOI L'ENTAILLE. Sans elle, le cutter attaque la face verticale ENTIERE
# du bloc : toute une colonne de noeuds entre en contact au meme pas, et avec
# une penalite de E x epaisseur sur des masses nodales de 4e-5 kg, un noeud
# legerement penetre gagne des centaines de m/s en un pas. Mesure du
# 2026-08-18 : divergence des le premier contact, 794 elements retournes,
# deplacement nodal de 6,5 m sur un bloc de 24 mm. Avec l'entaille, le cutter
# n'attaque qu'une marche de la hauteur de la passe.
#
# GRADATION : fine sur une BANDE le long de la surface (la ou la coupe se
# passe), grossiere en profondeur. Mesh.Algorithm = 5 (Delaunay isotrope).
# ---------------------------------------------------------------------------
import sys

import gmsh
import numpy as np


def main():
    if len(sys.argv) < 9:
        raise SystemExit("usage: make_cut_mesh.py W H depth notchLen hFine "
                         "bandH hFar out.msh [seed]")
    W, H, depth, notchLen, hFine, bandH, hFar = map(float, sys.argv[1:8])
    out = sys.argv[8]
    seed = int(sys.argv[9]) if len(sys.argv) > 9 else 1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cut")
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", 5)

    # bloc moins l'entaille : rectangle W x H auquel on retire le rectangle
    # [0, notchLen] x [H - depth, H]
    blk = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
    # notchLen <= 0 : PAS d'entaille. Mesure du 2026-08-18 - le coin PDC 2D
    # n'a pas de face de degagement sous l'arete, donc tout fond d'entaille
    # situe DERRIERE l'arete tombe dans le coin et se fait labourer. Une
    # entaille de 3 mm a ainsi consomme tout un run avant meme la coupe.
    if notchLen > 0.0:
        ntc = gmsh.model.occ.addRectangle(0, H - depth, 0, notchLen, depth)
        gmsh.model.occ.cut([(2, blk)], [(2, ntc)])
    gmsh.model.occ.synchronize()

    # champ de taille : distance a la SURFACE DE COUPE (bord superieur du bloc
    # + fond et flanc de l'entaille), fine sur bandH puis relachement
    top = []
    for dim, tag in gmsh.model.getEntities(1):
        x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(dim, tag)
        if y0 > H - depth - 1e-9:          # tout ce qui borde la zone de coupe
            top.append(tag)
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", top)
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 600)
    fth = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(fth, "InField", fd)
    gmsh.model.mesh.field.setNumber(fth, "SizeMin", hFine)
    gmsh.model.mesh.field.setNumber(fth, "SizeMax", hFar)
    gmsh.model.mesh.field.setNumber(fth, "DistMin", bandH)
    gmsh.model.mesh.field.setNumber(fth, "DistMax", bandH + 3.0 * bandH)
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
    print(f"[cut] eprouvette {W*1e3:g} x {H*1e3:g} mm, entaille {notchLen*1e3:g} "
          f"x {depth*1e3:.3f} mm (profondeur de passe), fond a y = "
          f"{(H-depth)*1e3:.3f} mm")
    print(f"[mesh] {len(tri)} triangles, h inscrit min/med = "
          f"{hin.min()*1e3:.4f}/{np.median(hin)*1e3:.4f} mm")
    gmsh.finalize()


if __name__ == "__main__":
    main()
