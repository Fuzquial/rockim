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
#
# ---------------------------------------------------------------------------
# AJOUTS OPT-IN DU 2026-09-11 (principe VIII : le defaut reste bit-identique).
# Sans aucun drapeau, ce script produit exactement le meme maillage qu'avant.
#
#   --cylinder
#       Eprouvette CYLINDRIQUE de DIAMETRE W et de hauteur H (D est alors
#       ignore), axe vertical, impact au centre de la face superieure. C'est
#       la geometrie du banc de Yang et al. 2026 (fig. 5) : 250 mm de
#       diametre, 150 mm de haut. Le pave par defaut reste disponible.
#       Pourquoi la forme compte : une eprouvette prismatique renvoie les
#       ondes par ses ARETES et ses coins, qui refocalisent, la ou le
#       cylindre les renvoie par une paroi de courbure constante. Sur un
#       essai lu en facies de fissuration, cette difference se voit.
#
#   --mid hMid rMid
#       Gradation a DEUX crans au lieu d'un : hFine jusqu'a rFine, puis
#       hMid jusqu'a rMid, puis croissance vers hFar. C'est litteralement la
#       fig. 6 de Yang et al. 2026 — 1 mm dans la demi-boule R 12,5 mm,
#       2 mm jusqu'a R 25 mm, 10 mm au bord. Un seul cran force a choisir
#       entre payer le 1 mm trop loin ou degrader trop vite la zone ou les
#       fissures radiales se propagent.
#       Realise par un champ Min de deux Threshold : Gmsh prend en chaque
#       point la plus petite des deux tailles, ce qui enchaine les paliers
#       sans discontinuite.
#
# RAPPEL DE LA REGLE DE LA THESE : maillage NON STRUCTURE obligatoire en FDEM
# (DOCUMENTATION_rockim.md, 8.4 — un maillage structure est une CONDITION
# D'INVALIDITE : trajets de fissure biaises, divergence en phase debris).
# Delaunay + optimisation Netgen, aucune direction privilegiee.
# ---------------------------------------------------------------------------
import sys

import gmsh


def main():
    argv = list(sys.argv[1:])

    # ---- drapeaux opt-in, retires avant la lecture positionnelle ---------
    cylinder = False
    if "--cylinder" in argv:
        argv.remove("--cylinder")
        cylinder = True
    hMid = rMid = None
    if "--mid" in argv:
        i = argv.index("--mid")
        hMid, rMid = float(argv[i + 1]), float(argv[i + 2])
        del argv[i:i + 3]

    if len(argv) < 7:
        raise SystemExit("usage: make_impact3d_mesh.py W D H hFine rFine hFar "
                         "out.msh [seed] [--cylinder] [--mid hMid rMid]")
    W, D, H, hFine, rFine, hFar = map(float, argv[:6])
    out = argv[6]
    seed = int(argv[7]) if len(argv) > 7 else 1

    if hMid is not None and not (hFine <= hMid <= hFar and rFine <= rMid):
        raise SystemExit("--mid : il faut hFine <= hMid <= hFar et "
                         "rFine <= rMid")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("impact3d")
    if cylinder:
        # cylindre d'AXE Z, diametre W, hauteur H, centre en (0, 0)
        gmsh.model.occ.addCylinder(0.0, 0.0, 0.0, 0.0, 0.0, H, W / 2.0)
        xImp, yImp = 0.0, 0.0
    else:
        gmsh.model.occ.addBox(0, 0, 0, W, D, H)
        xImp, yImp = W / 2.0, D / 2.0
    gmsh.model.occ.synchronize()

    px = gmsh.model.occ.addPoint(xImp, yImp, H, hFine)
    gmsh.model.occ.synchronize()
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "PointsList", [px])
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", hFine)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", hFar)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", rFine)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", 8.0 * rFine)
    bg = ft
    if hMid is not None:
        # second palier : hMid tenu jusqu'a rMid, puis montee vers hFar.
        # Le champ Min des deux enchaine les paliers : dans R rFine le
        # premier impose hFine, entre rFine et rMid le second impose hMid,
        # au-dela les deux montent vers hFar.
        f2 = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f2, "InField", fd)
        gmsh.model.mesh.field.setNumber(f2, "SizeMin", hMid)
        gmsh.model.mesh.field.setNumber(f2, "SizeMax", hFar)
        gmsh.model.mesh.field.setNumber(f2, "DistMin", rMid)
        gmsh.model.mesh.field.setNumber(f2, "DistMax", 8.0 * rMid)
        bg = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(bg, "FieldsList", [ft, f2])
    gmsh.model.mesh.field.setAsBackgroundMesh(bg)
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
    if cylinder:
        print(f"[impact3d] eprouvette CYLINDRIQUE diam {W*1e3:.0f} x "
              f"{H*1e3:.0f} mm, impact en (0, 0, {H*1e3:.0f}) mm")
    else:
        print(f"[impact3d] eprouvette {W*1e3:.0f} x {D*1e3:.0f} x "
              f"{H*1e3:.0f} mm, impact en ({W/2*1e3:.0f}, {D/2*1e3:.0f}, "
              f"{H*1e3:.0f}) mm")
    if hMid is not None:
        print(f"[impact3d] maille {hFine*1e3:.2f} mm jusqu'a R "
              f"{rFine*1e3:.1f} mm, {hMid*1e3:.2f} mm jusqu'a R "
              f"{rMid*1e3:.1f} mm, {hFar*1e3:.1f} mm au bord")
    else:
        print(f"[impact3d] maille {hFine*1e3:.2f} mm jusqu'a R "
              f"{rFine*1e3:.1f} mm, {hFar*1e3:.1f} mm au-dela de R "
              f"{8*rFine*1e3:.0f} mm")
    print(f"[mesh] {n} tetraedres -> {out}")
    gmsh.finalize()


if __name__ == "__main__":
    main()
