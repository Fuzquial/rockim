"""Maillage : Gmsh pour le régime homogène, run tronqué pour l'aperçu Voronoï.

Gmsh est appelé par le script existant tools/make_unstructured_mesh.py (box2d,
Delaunay 2D, MSH 2.2), dans un processus séparé : une seule implémentation, et
un plantage de Gmsh n'emporte pas l'interface.

L'aperçu d'un maillage Voronoï (GBM) ne peut venir que du solveur, qui génère
lui-même grains et triangles. En attendant l'option `--mesh-only` (spec 007,
S3), on passe par un run tronqué à une frame (T = 2 µs), comme les bancs courts
de la campagne (gen_decks.py --smoke). C'est un lancement du solveur : il passe
par la file de calculs, donc par une action explicite de l'utilisateur.
"""
import copy
import os
import subprocess
import sys

import numpy as np

from .geometrie import elements_par_grain, nombre_elements

RACINE_G1 = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCRIPT_GMSH = os.path.join(RACINE_G1, "tools", "make_unstructured_mesh.py")


def generer_gmsh(W, H, h, graine, chemin, delai=300):
    """Écrit un maillage box2d non structuré ; rend ses statistiques."""
    os.makedirs(os.path.dirname(os.path.abspath(chemin)), exist_ok=True)
    r = subprocess.run([sys.executable, SCRIPT_GMSH, "box2d", repr(W), repr(H), repr(h),
                        chemin, str(graine)], capture_output=True, text=True, timeout=delai)
    if r.returncode != 0 or not os.path.exists(chemin):
        raise RuntimeError("Gmsh a échoué (code %d) :\n%s" % (r.returncode, (r.stderr or r.stdout)[-2000:]))
    return stats_msh(lire_msh(chemin))


def lire_msh(chemin):
    """Lecteur MSH 2.2 ASCII réduit aux triangles, comme le solveur (buildMeshFile)."""
    noeuds, tris = {}, []
    with open(chemin, encoding="utf-8", errors="replace") as f:
        lignes = iter(f.read().splitlines())
    for l in lignes:
        if l.startswith("$MeshFormat"):
            v = next(lignes).split()[0]
            if not v.startswith("2"):
                raise ValueError("format MSH %s : le solveur n'accepte que 2.2" % v)
        elif l.startswith("$Nodes"):
            for _ in range(int(next(lignes))):
                i, x, y, _z = next(lignes).split()[:4]
                noeuds[int(i)] = (float(x), float(y))
        elif l.startswith("$Elements"):
            for _ in range(int(next(lignes))):
                c = next(lignes).split()
                if c[1] == "2":                      # triangle à 3 nœuds
                    nt = int(c[2])
                    tris.append([int(v) for v in c[3 + nt:6 + nt]])
    ids = np.array(sorted(noeuds))
    rang = {v: k for k, v in enumerate(ids)}
    xy = np.array([noeuds[i] for i in ids])
    t = np.array([[rang[v] for v in tri] for tri in tris], dtype=np.int32)
    return xy, t


def stats_msh(maillage):
    xy, t = maillage
    p = xy[t]
    aires = 0.5 * np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                         - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
    cotes = np.stack([np.linalg.norm(p[:, (k + 1) % 3] - p[:, k], axis=1) for k in range(3)], 1)
    # Hauteur minimale de chaque triangle : c'est elle qui pilote le pas de temps.
    hmin = (2 * aires / cotes.max(1)).min()
    return dict(noeuds=int(len(xy)), triangles=int(len(t)), aire=float(aires.sum()),
                cote_moyen=float(cotes.mean()), hauteur_min=float(hmin),
                W=float(np.ptp(xy[:, 0])), H=float(np.ptp(xy[:, 1])))


def estimation_voronoi(essai):
    """Ce qu'on sait d'un maillage Voronoï AVANT de lancer quoi que ce soit."""
    m, e = essai.maillage, essai.eprouvette
    return dict(elements=int(round(nombre_elements(e.W, e.H, m.taille_element))),
                elements_par_grain=elements_par_grain(m.taille_grain, m.taille_element),
                grains=int(round(e.W * e.H / (0.785 * m.taille_grain ** 2))))


def essai_apercu(essai):
    """Le même essai, tronqué à une frame : sert à voir le maillage réel du solveur."""
    a = copy.deepcopy(essai)
    a.nom = essai.nom + "_apercu"
    a.sorties.T = 2.0e-6
    a.sorties.frames = 1
    return a
