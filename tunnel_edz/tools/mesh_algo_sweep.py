#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# mesh_algo_sweep.py — trouver un mailleur qui ne fabrique PAS un réseau.
#
# Constat du 2026-08-17 : `Mesh.Algorithm = 6` (frontal-Delaunay) avec un champ
# de taille CONSTANT produit un pavage quasi équilatéral — médiane du plus
# petit angle a 60,0 deg, histogramme d'orientation périodique a 60 deg,
# pic/creux = 19. Les fissures n'ont alors que trois directions disponibles et
# suivent des droites sur plusieurs mètres. C'est une condition d'invalidité
# pour un calcul de faciès de fissuration.
#
# Ce script fabrique le MEME domaine avec plusieurs réglages et mesure
# l'isotropie de chacun, pour choisir sur des nombres.
#
#   python tunnel_edz/tools/mesh_algo_sweep.py [--out dossier]
#
# Variantes essayées :
#   algo6            l'existant (référence du défaut)
#   algo5            Delaunay
#   algo1            MeshAdapt
#   algo6_nosmooth   frontal-Delaunay sans lissage laplacien
#   algo5_perturb    Delaunay + champ de taille bruité (+-12 %, incommensurable)
#   algo6_perturb    frontal-Delaunay + le meme bruit
#
# Le bruit de taille est le levier le plus sûr : il casse la commensurabilité
# qui permet au pavage de se refermer sur lui-meme, quel que soit l'algorithme.
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import gmsh
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
sys.path.insert(0, HERE)
from make_unstructured_mesh import build_tunnel_hs, TUNNEL_HS  # noqa: E402
from mesh_isotropy import read_msh22, stats  # noqa: E402

W = H = 60.0                      # domaine reduit : la mesure est locale
HFINE, RFINE, HFAR = 0.22, 8.0, 2.0

VARIANTS = [
    ("algo6",          dict(algo=6, smooth=1, perturb=0.0)),
    ("algo5",          dict(algo=5, smooth=1, perturb=0.0)),
    ("algo1",          dict(algo=1, smooth=1, perturb=0.0)),
    ("algo6_nosmooth", dict(algo=6, smooth=0, perturb=0.0)),
    ("algo5_perturb",  dict(algo=5, smooth=1, perturb=0.12)),
    ("algo6_perturb",  dict(algo=6, smooth=1, perturb=0.12)),
]


def build(path, algo, smooth, perturb, seed=1):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("sweep")
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", algo)
    gmsh.option.setNumber("Mesh.Smoothing", smooth)
    cx, cy0 = 0.5 * W, 0.5 * H - 0.5 * TUNNEL_HS["height"]
    plate = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
    loop, g = build_tunnel_hs(cx, cy0)
    cav = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.cut([(2, plate)], [(2, cav)])
    gmsh.model.occ.synchronize()
    wall = []
    for dim, tag in gmsh.model.getEntities(1):
        x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(dim, tag)
        if (x0 > cx - g["halfSpan"] - 1e-6 and x1 < cx + g["halfSpan"] + 1e-6
                and y0 > cy0 - 1e-6 and y1 < cy0 + g["height"] + 1e-6):
            wall.append(tag)
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", wall)
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 400)
    fth = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(fth, "InField", fd)
    gmsh.model.mesh.field.setNumber(fth, "SizeMin", HFINE)
    gmsh.model.mesh.field.setNumber(fth, "SizeMax", HFAR)
    gmsh.model.mesh.field.setNumber(fth, "DistMin", RFINE)
    gmsh.model.mesh.field.setNumber(fth, "DistMax", 1.8 * RFINE)
    use = fth
    if perturb > 0.0:
        # bruit multiplicatif a longueurs d'onde INCOMMENSURABLES : le pavage
        # ne peut plus se refermer, quel que soit l'algorithme.
        fp = gmsh.model.mesh.field.add("MathEval")
        gmsh.model.mesh.field.setString(
            fp, "F", f"F{fth} * (1 + {perturb} * sin(1.9319*x + 0.7)"
                     f" * sin(2.7183*y + 1.3)"
                     f" + {0.6 * perturb} * sin(4.6692*x - 2.1))")
        use = fp
    gmsh.model.mesh.field.setAsBackgroundMesh(use)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", HFINE * (1.0 - perturb))
    gmsh.option.setNumber("Mesh.MeshSizeMax", HFAR)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(path)
    _, _, conn = gmsh.model.mesh.getElements(2)
    n = len(conn[0]) // 3
    gmsh.finalize()
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "..", "meshtest"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for name, kw in VARIANTS:
        path = os.path.join(a.out, f"sweep_{name}.msh")
        n = build(path, **kw)
        P, T = read_msh22(path)
        R, ratio, amed = stats(P, T, 0.5 * W, 0.5 * H, 12.0, name)
        rows.append((name, n, R, ratio, amed))
        print()
    print("=" * 68)
    print(f"{'variante':16s} {'triangles':>10s} {'R6':>7s} {'pic/creux':>10s} "
          f"{'angle min med':>14s}")
    for name, n, R, ratio, amed in rows:
        print(f"{name:16s} {n:10d} {R:7.3f} {ratio:10.2f} {amed:13.1f} deg")
    print("\nCible : R6 et pic/creux les plus bas possible, angle median "
          "nettement sous 60 deg.")


if __name__ == "__main__":
    main()
