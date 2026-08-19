#!/usr/bin/env python3
# make_impact_uni.py — bloc de roche + insert spherique, maillage UNIFORME.
#   python tunnel_edz/tools/make_impact_uni.py W HR RI h hIns out.msh [seed]
import sys, gmsh, numpy as np

W, HR, RI, h, hIns = map(float, sys.argv[1:6])
out = sys.argv[6]
seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("Mesh.RandomSeed", seed / 100.0)
rock = gmsh.model.occ.addBox(-0.5*W, -0.5*W, -HR, W, W, HR)
ins = gmsh.model.occ.addSphere(0.0, 0.0, RI + 1e-4, RI)
gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(3, [rock], name="rock")
gmsh.model.addPhysicalGroup(3, [ins], name="insert")

fm = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(fm, "F", "%g" % h)
fr = gmsh.model.mesh.field.add("Restrict")
gmsh.model.mesh.field.setNumber(fr, "InField", fm)
gmsh.model.mesh.field.setNumbers(fr, "VolumesList", [rock])
fi = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(fi, "F", "%g" % hIns)
fri = gmsh.model.mesh.field.add("Restrict")
gmsh.model.mesh.field.setNumber(fri, "InField", fi)
gmsh.model.mesh.field.setNumbers(fri, "VolumesList", [ins])
fmin = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", [fr, fri])
gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
for k in ("Mesh.MeshSizeFromPoints","Mesh.MeshSizeFromCurvature","Mesh.MeshSizeExtendFromBoundary"):
    gmsh.option.setNumber(k, 0)
gmsh.option.setNumber("Mesh.Algorithm", 5)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.mesh.generate(3)
gmsh.write(out)

_, _, npar = gmsh.model.mesh.getElements(3)
tags = np.concatenate([np.array(x).reshape(-1,4) for x in npar])
nid, nc, _ = gmsh.model.mesh.getNodes()
nc = nc.reshape(-1,3); idx = {int(t): k for k,t in enumerate(nid)}
P = np.array([nc[idx[int(t)]] for t in tags.ravel()]).reshape(-1,4,3)
cen = P.mean(axis=1); e = np.zeros(len(P))
for i in range(4):
    for j in range(i+1,4): e += np.linalg.norm(P[:,i]-P[:,j],axis=1)
e /= 6.0
haut = cen[:,2] > 0
print("[maillage] %d tetraedres : roche %d (arete med %.2f mm), insert %d (%.2f mm)"
      % (len(P), (~haut).sum(), np.median(e[~haut])*1e3, haut.sum(), np.median(e[haut])*1e3))
gmsh.finalize()
