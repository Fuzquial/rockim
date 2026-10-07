# -*- coding: utf-8 -*-
"""make_bar2.py — barre a DEUX COUCHES nommees, maillage Gmsh MSH 2.2 ASCII.

Sert le banc de REUSS (F3) : la barre est empilee selon z, la moitie basse
porte le volume physique « dur », la moitie haute « mou ». Les deux couches
partagent leurs noeuds a l'interface (maillage CONFORME) — c'est exactement
ce que le mode fem3d exige, puisqu'il n'a aucun contact entre corps : deux
volumes non conformes s'y traverseraient sans un mot.

Decoupage de Kuhn (6 tets par hexaedre, meme diagonale principale pour toutes
les mailles) : conforme d'une maille a l'autre, donc pas un seul noeud
dedouble. Les tets ne sont pas de belle qualite, mais ce banc mesure un module
apparent en serie, pas un trajet de fissure.

usage : python bench_phases/make_bar2.py bench_phases/bar2.msh
"""
import sys

NX, NY, NZ = 2, 2, 16          # 16 couches : l'interface tombe pile a z = H/2
W, D, H = 2.0e-3, 2.0e-3, 16.0e-3   # elancement 8 : effets de mors limites

# Decoupage de Kuhn : sommet local a = ix + 2 iy + 4 iz ; un tet par
# permutation des axes, tous issus de la diagonale 0-7.
KUHN = [(0, 1, 3, 7), (0, 1, 5, 7), (0, 2, 3, 7),
        (0, 2, 6, 7), (0, 4, 5, 7), (0, 4, 6, 7)]


def main(path, nommes=True):
    """nommes = False : le MEME maillage sans $PhysicalNames, pour le cas
    fautif « phases demandees sans source de phase » (le deck declare des
    mineraux, le maillage est un seul corps : une seule phase s'appliquerait,
    en silence — c'est refuse)."""
    nid = lambda i, j, k: 1 + i + (NX + 1) * (j + (NY + 1) * k)
    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat"]
    if nommes:
        lines += ["$PhysicalNames", "2", '3 1 "dur"', '3 2 "mou"',
                  "$EndPhysicalNames"]
    lines += ["$Nodes", str((NX + 1) * (NY + 1) * (NZ + 1))]
    for k in range(NZ + 1):
        for j in range(NY + 1):
            for i in range(NX + 1):
                lines.append("%d %.17g %.17g %.17g"
                             % (nid(i, j, k), i * W / NX, j * D / NY,
                                k * H / NZ))
    lines += ["$EndNodes", "$Elements"]
    el, e = [], 0
    for k in range(NZ):
        phys = 1 if k < NZ // 2 else 2          # bas = dur, haut = mou
        for j in range(NY):
            for i in range(NX):
                c = [nid(i + (a & 1), j + ((a >> 1) & 1), k + ((a >> 2) & 1))
                     for a in range(8)]
                for t in KUHN:
                    e += 1
                    el.append("%d 4 2 %d %d %d %d %d %d"
                              % (e, phys, phys, c[t[0]], c[t[1]],
                                 c[t[2]], c[t[3]]))
    lines.append(str(len(el)))
    lines += el
    lines += ["$EndElements", ""]
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(lines))
    print("%s : %d noeuds, %d tets (dur %d / mou %d), boite %g x %g x %g m"
          % (path, (NX + 1) * (NY + 1) * (NZ + 1), len(el), len(el) // 2,
             len(el) // 2, W, D, H))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--sans-noms"]
    main(args[0] if args else "bench_phases/bar2.msh",
         "--sans-noms" not in sys.argv)
