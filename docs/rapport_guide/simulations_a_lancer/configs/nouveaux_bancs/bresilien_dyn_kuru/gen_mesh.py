#!/usr/bin/env python3
# Maillage du disque brésilien QUASI STATIQUE du banc B11 (granite de Kuru,
# Saksala et al. 2013, D = 40,8 mm). Recette de calib_quick/make_disc_mesh.py
# (Gmsh Delaunay, champ de taille bruité seedé), petits méplats de 2 x 5 degrés
# pour asseoir les platines. Lancer depuis la racine du dépôt :
#   python3 docs/rapport_guide/simulations_a_lancer/configs/nouveaux_bancs/bresilien_dyn_kuru/gen_mesh.py
# Sortie : docs/rapport_guide/simulations_a_lancer/meshes/nb_disc41_h05_f5_s1.msh
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *[".."] * 6))
sys.path.insert(0, os.path.join(ROOT, "calib_quick"))
from make_disc_mesh import build  # noqa: E402

out = os.path.join(ROOT, "docs", "rapport_guide", "simulations_a_lancer", "meshes",
                   "nb_disc41_h05_f5_s1.msh")
nodes, tris = build(0.0408, 0.5e-3, 5.0, 1, out)
print(f"{out} : {len(tris)} triangles")
