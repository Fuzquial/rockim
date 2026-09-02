#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_box_mesh.py — eprouvette rectangulaire W x H, maillage Gmsh DELAUNAY
# (Mesh.Algorithm = 5) + champ de taille bruite (+-10 %) dont les phases
# dependent d une graine : plusieurs REALISATIONS independantes du meme
# maillage isotrope, pour mesurer le bruit de realisation de la calibration.
# Le frontal (algo 6) est BANNI (pavage quasi equilateral, R6 = 0,34).
#
#   python calib_quick/make_box_mesh.py --W 0.02 --H 0.04 --h 0.8e-3 --seeds 1 2 3
#
# Sortie : calib_quick/meshes/box{W}x{H}_h{h}_s{seed}.msh (format 2.2 ASCII,
# comme box20x40_h08_algo5.msh) + metriques d isotropie sur stdout :
#   pic/creux  = max/min de l histogramme des orientations d aretes (36 classes)
#   R6         = ordre orientationnel 6-fold |<exp(6 i theta)>| (0 isotrope, 1 nid d abeille)
#   angle min  = mediane de l angle minimal des triangles
# ---------------------------------------------------------------------------
import argparse
import math
import os
import sys
import numpy as np
import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))


def isotropy(nodes, tris):
    p = nodes[tris]                      # (n, 3, 2)
    e = np.concatenate([p[:, 1] - p[:, 0], p[:, 2] - p[:, 1], p[:, 0] - p[:, 2]])
    th = np.arctan2(e[:, 1], e[:, 0]) % math.pi
    hist, _ = np.histogram(th, bins=36, range=(0, math.pi))
    r6 = abs(np.exp(6j * th).mean())
    a, b, c = (np.linalg.norm(p[:, 1] - p[:, 2], axis=1), np.linalg.norm(p[:, 2] - p[:, 0], axis=1),
               np.linalg.norm(p[:, 0] - p[:, 1], axis=1))
    def ang(x, y, z):
        return np.degrees(np.arccos(np.clip((y * y + z * z - x * x) / (2 * y * z), -1, 1)))
    amin = np.minimum(np.minimum(ang(a, b, c), ang(b, c, a)), ang(c, a, b))
    return hist.max() / max(1, hist.min()), r6, float(np.median(amin))


def build(W, H, h, seed, out):
    rng = np.random.default_rng(seed)
    k = rng.uniform(3000, 12000, size=4)          # nombres d onde des trois modes
    ph = rng.uniform(0, 2 * math.pi, size=3)      # phases
    expr = (f"{h}*(1+0.10*sin(x*{k[0]:.1f}+{ph[0]:.3f})*cos(y*{k[1]:.1f}+{ph[1]:.3f})"
            f"+0.06*sin(x*{k[2]:.1f}+y*{k[3]:.1f}+{ph[2]:.3f}))")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("box")
    gmsh.model.occ.addRectangle(0, 0, 0, W, H)
    gmsh.model.occ.synchronize()
    f = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(f, "F", expr)
    gmsh.model.mesh.field.setAsBackgroundMesh(f)
    gmsh.option.setNumber("Mesh.Algorithm", 5)              # Delaunay — JAMAIS 6 (frontal)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    try:
        gmsh.option.setNumber("Mesh.RandomSeed", float(seed))
    except Exception:
        pass
    gmsh.option.setNumber("Mesh.RandomFactor", 1e-9 * (1 + seed))
    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.mesh.generate(2)
    ntags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(coords).reshape(-1, 3)[:, :2]
    idx = {t: i for i, t in enumerate(ntags)}
    etypes, etags, enodes = gmsh.model.mesh.getElements(2)
    tris = np.array([[idx[n] for n in enodes[0][3 * i:3 * i + 3]] for i in range(len(etags[0]))])
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(out)
    gmsh.finalize()
    return nodes, tris


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--W", type=float, default=0.02)
    ap.add_argument("--H", type=float, default=0.04)
    ap.add_argument("--h", type=float, default=0.8e-3)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    a = ap.parse_args()
    os.makedirs(os.path.join(HERE, "meshes"), exist_ok=True)
    print(f"{'maillage':34s} {'tri':>6s} {'pic/creux':>9s} {'R6':>6s} {'angle min':>9s}")
    for s in a.seeds:
        name = f"box{int(round(a.W*1e3))}x{int(round(a.H*1e3))}_h{a.h*1e3:g}_s{s}.msh".replace(".", "p", 1) if False else \
            f"box{int(round(a.W*1e3))}x{int(round(a.H*1e3))}_h{int(round(a.h*1e4)):02d}_s{s}.msh"
        out = os.path.join(HERE, "meshes", name)
        nodes, tris = build(a.W, a.H, a.h, s, out)
        pc, r6, am = isotropy(nodes, tris)
        print(f"{name:34s} {len(tris):6d} {pc:9.2f} {r6:6.3f} {am:8.1f}°")


if __name__ == "__main__":
    main()
