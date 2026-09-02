#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_disc_mesh.py — disque brésilien Ø D, maillage Gmsh DELAUNAY (algo 5) +
# champ de taille bruité seedé (même recette que make_box_mesh.py), avec
# deux méplats optionnels (angle 2·flatDeg vu du centre) pour le contact des
# plateaux. Sortie MSH 2.2 ASCII, disque centré en (D/2, D/2) pour que la
# boîte englobante soit [0, D] x [0, D] (rockim lit W_/H_ sur la boîte).
#
#   python calib_quick/make_disc_mesh.py --D 0.04 --h 0.8e-3 --flat 20 --seeds 1
# ---------------------------------------------------------------------------
import argparse
import math
import os
import sys
import numpy as np
import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_box_mesh import isotropy  # noqa: E402


def build(D, h, flat_deg, seed, out):
    rng = np.random.default_rng(seed)
    k = rng.uniform(3000, 12000, size=4); ph = rng.uniform(0, 2 * math.pi, size=3)
    expr = (f"{h}*(1+0.10*sin(x*{k[0]:.1f}+{ph[0]:.3f})*cos(y*{k[1]:.1f}+{ph[1]:.3f})"
            f"+0.06*sin(x*{k[2]:.1f}+y*{k[3]:.1f}+{ph[2]:.3f}))")
    R = 0.5 * D
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("disc")
    disk = gmsh.model.occ.addDisk(R, R, 0, R, R)
    if flat_deg > 0:
        # méplats : on retire deux calottes au-dessus et au-dessous
        c = R * math.cos(math.radians(flat_deg))      # distance centre-méplat
        top = gmsh.model.occ.addRectangle(0, R + c, 0, D, R)      # y > R + c
        bot = gmsh.model.occ.addRectangle(0, -c, 0, D, R)         # y < R - c
        gmsh.model.occ.cut([(2, disk)], [(2, top), (2, bot)])
    gmsh.model.occ.synchronize()
    f = gmsh.model.mesh.field.add("MathEval"); gmsh.model.mesh.field.setString(f, "F", expr)
    gmsh.model.mesh.field.setAsBackgroundMesh(f)
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    surfs = [s[1] for s in gmsh.model.getEntities(2)]
    gmsh.model.addPhysicalGroup(2, surfs, 1)
    gmsh.model.mesh.generate(2)
    ntags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(coords).reshape(-1, 3)[:, :2]
    idx = {t: i for i, t in enumerate(ntags)}
    _, etags, enodes = gmsh.model.mesh.getElements(2)
    tris = np.array([[idx[n] for n in enodes[0][3 * i:3 * i + 3]] for i in range(len(etags[0]))])
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(out); gmsh.finalize()
    return nodes, tris


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--D", type=float, default=0.04)
    ap.add_argument("--h", type=float, default=0.8e-3)
    ap.add_argument("--flat", type=float, default=20.0, help="demi-angle des meplats [deg], 0 = disque complet")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1])
    a = ap.parse_args()
    for s in a.seeds:
        name = f"disc{int(round(a.D*1e3))}_h{int(round(a.h*1e4)):02d}_f{int(a.flat)}_s{s}.msh"
        out = os.path.join(HERE, "meshes", name)
        nodes, tris = build(a.D, a.h, a.flat, s, out)
        pc, r6, am = isotropy(nodes, tris)
        print(f"{name:30s} {len(tris):6d} tri  boite x [{nodes[:,0].min()*1e3:.2f};{nodes[:,0].max()*1e3:.2f}] "
              f"y [{nodes[:,1].min()*1e3:.2f};{nodes[:,1].max()*1e3:.2f}] mm  pic/creux {pc:.2f}  R6 {r6:.3f}  angle min {am:.1f}")


if __name__ == "__main__":
    main()
