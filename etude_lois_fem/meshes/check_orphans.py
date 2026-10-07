# -*- coding: utf-8 -*-
"""Noeuds ORPHELINS (jamais references par un tetraedre) dans un MSH 2.2 : un point geometrique ajoute pour un
champ de taille (Distance) est maille comme noeud 0-D par Gmsh ; dans fem3d il devient un noeud de masse nulle,
epingle (FIXED), mais le contact outil par penalite lui applique quand meme une force (broche fantome sous le
pole de l'insert) et la voie Signorini divise par sa masse (NaN). Decouvert le 2026-09-05 sur Q1_c05.msh.
usage : python check_orphans.py <fichiers.msh ...>     (rapport seulement, ne modifie rien)"""
import io, os, sys
import numpy as np

def orphans(p):
    L = io.open(p, encoding="utf-8", errors="ignore").read().split("\n")
    i = L.index("$Nodes"); n = int(L[i + 1])
    nodes = np.array([[float(x) for x in l.split()[1:4]] for l in L[i + 2:i + 2 + n]])
    j = L.index("$Elements"); m = int(L[j + 1]); used = set()
    for l in L[j + 2:j + 2 + m]:
        t = l.split()
        if t[1] == "4": used.update(int(x) - 1 for x in t[-4:])
    orph = [k for k in range(n) if k not in used]
    return n, orph, nodes

if __name__ == "__main__":
    for p in sys.argv[1:]:
        n, o, nodes = orphans(p)
        print("%-26s noeuds %6d  orphelins %d  %s" % (os.path.basename(p), n, len(o), [tuple(np.round(nodes[k] * 1e3, 2)) for k in o[:4]]))
