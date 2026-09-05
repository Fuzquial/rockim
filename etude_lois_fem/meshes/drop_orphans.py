# -*- coding: utf-8 -*-
"""Retire les noeuds ORPHELINS d'un MSH 2.2 (jamais references par un tetraedre) et renumerote — sans toucher
aux tetraedres ni a leur ordre. Le fichier d'origine n'est PAS modifie : sortie <nom>_clean.msh.
usage : python drop_orphans.py <fichier.msh ...>"""
import io, os, sys

def clean(p):
    L = io.open(p, encoding="utf-8", errors="ignore").read().split("\n")
    i = L.index("$Nodes"); n = int(L[i + 1]); nodes = L[i + 2:i + 2 + n]
    j = L.index("$Elements"); m = int(L[j + 1]); elems = L[j + 2:j + 2 + m]
    used = set()
    for l in elems:
        t = l.split()
        if t[1] == "4": used.update(int(x) for x in t[-4:])
    old_ids = [int(l.split()[0]) for l in nodes]
    keep = [k for k in old_ids if k in used]
    renum = {k: r + 1 for r, k in enumerate(keep)}
    new_nodes = [l for l in nodes if int(l.split()[0]) in renum]
    new_nodes = ["%d %s" % (renum[int(l.split()[0])], " ".join(l.split()[1:])) for l in new_nodes]
    new_elems = []
    for l in elems:
        t = l.split()
        if t[1] != "4": continue                       # on ne garde que les tetraedres
        new_elems.append(" ".join(t[:-4] + [str(renum[int(x)]) for x in t[-4:]]))
    new_elems = ["%d %s" % (r + 1, " ".join(l.split()[1:])) for r, l in enumerate(new_elems)]
    out = os.path.splitext(p)[0] + "_clean.msh"
    with io.open(out, "w", encoding="ascii", newline="\n") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n" % len(new_nodes))
        f.write("\n".join(new_nodes) + "\n$EndNodes\n$Elements\n%d\n" % len(new_elems))
        f.write("\n".join(new_elems) + "\n$EndElements\n")
    print("%s : %d noeuds -> %d (%d orphelins retires), %d tetraedres -> %s" % (os.path.basename(p), n, len(new_nodes), n - len(new_nodes), len(new_elems), os.path.basename(out)))

if __name__ == "__main__":
    for p in sys.argv[1:]:
        clean(p)
