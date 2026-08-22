"""worker.py — le worker gmsh de studio (spec 006, WP2.1), HORS processus.

Lancé par gmsh_service via QProcess :
    python -m rockim_studio.geometry.worker recette.json sortie.msh
Lit la recette JSON, construit la géométrie OCC, maille, écrit le MSH 2.2
ASCII (le format que lit rockim) plus un aperçu .vtu (pour la scène), et
imprime UNE ligne JSON de résultat sur stdout :
    {"ok": true, "triangles": N, "hMin": ..., "hMed": ..., "vtu": "..."}

Motifs du processus séparé : gmsh est mono-session, ses crashs OCC ne
doivent pas emporter la GUI, et un maillage long ne doit pas geler l'UI.

Le profil fer à cheval est IMPORTÉ de tools/make_unstructured_mesh.py
(build_tunnel_hs) — une seule implémentation, celle validée par l'étude
tunnel (constitution : chercher l'existant avant de construire).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

_STUDIO = Path(__file__).resolve().parents[2]
_ROOT = _STUDIO.parents[1]
sys.path.insert(0, str(_STUDIO))
sys.path.insert(0, str(_ROOT / "tools"))

from rockim_studio.geometry.recipes import MeshRecipe   # noqa: E402


def _fail(msg: str) -> int:
    print(json.dumps({"ok": False, "error": msg}))
    return 1


def main() -> int:
    if len(sys.argv) != 3:
        return _fail("usage: worker.py recette.json sortie.msh")
    recipe_path, out_msh = Path(sys.argv[1]), Path(sys.argv[2])
    try:
        recipe = MeshRecipe.from_dict(
            json.loads(recipe_path.read_text(encoding="utf-8")))
    except Exception as e:
        return _fail(f"recette illisible : {e}")
    errs = recipe.validate()
    if errs:
        return _fail(" ; ".join(errs))

    try:
        import gmsh
    except ImportError as e:
        return _fail(f"gmsh indisponible : {e}")

    try:
        import numpy as np
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("studio")
        gmsh.option.setNumber("Mesh.RandomSeed", recipe.seed)

        W, H, c = recipe.W, recipe.H, recipe.cavity
        plate = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
        if c.kind == "circle":
            hole = gmsh.model.occ.addDisk(c.cx, c.cy, 0, c.r, c.r)
        elif c.kind == "horseshoe":
            from make_unstructured_mesh import TUNNEL_HS, build_tunnel_hs
            prof = {k: v * c.scale if isinstance(v, (int, float)) else v
                    for k, v in TUNNEL_HS.items()}
            cy0 = c.cy - 0.5 * prof["height"]
            loop, _g = build_tunnel_hs(c.cx, cy0, prof)
            hole = gmsh.model.occ.addPlaneSurface([loop])
        else:                                   # polygon
            pts = [gmsh.model.occ.addPoint(x, y, 0) for x, y in c.points]
            lines = [gmsh.model.occ.addLine(pts[i], pts[(i + 1) % len(pts)])
                     for i in range(len(pts))]
            loop = gmsh.model.occ.addCurveLoop(lines)
            hole = gmsh.model.occ.addPlaneSurface([loop])
        gmsh.model.occ.cut([(2, plate)], [(2, hole)])
        gmsh.model.occ.synchronize()

        # Courbes de paroi de la cavité = celles qui ne bordent pas la plaque.
        wall = []
        eps = 1e-6
        for dim, tag in gmsh.model.getEntities(1):
            x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(dim, tag)
            on_border = (x0 < eps or y0 < eps
                         or x1 > W - eps or y1 > H - eps)
            if not on_border:
                wall.append(tag)
        if not wall:
            return _fail("contour de cavité introuvable après la coupe")

        fd = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fd, "CurvesList", wall)
        gmsh.model.mesh.field.setNumber(fd, "Sampling", 400)
        fth = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(fth, "InField", fd)
        gmsh.model.mesh.field.setNumber(fth, "SizeMin", recipe.hFine)
        gmsh.model.mesh.field.setNumber(fth, "SizeMax", recipe.hFar)
        gmsh.model.mesh.field.setNumber(fth, "DistMin", recipe.rFine)
        # transition étalée sur 0,8·rFine : une marche de taille dégrade la
        # qualité donc le dt (leçon bench1g)
        gmsh.model.mesh.field.setNumber(fth, "DistMax",
                                        1.8 * recipe.rFine)
        gmsh.model.mesh.field.setAsBackgroundMesh(fth)
        for opt in ("MeshSizeExtendFromBoundary", "MeshSizeFromPoints",
                    "MeshSizeFromCurvature"):
            gmsh.option.setNumber(f"Mesh.{opt}", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", recipe.hFine)
        gmsh.option.setNumber("Mesh.MeshSizeMax", recipe.hFar)
        gmsh.option.setNumber("Mesh.Algorithm", recipe.algo2d)
        gmsh.model.mesh.generate(2)

        out_msh.parent.mkdir(parents=True, exist_ok=True)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(str(out_msh))
        vtu = out_msh.with_suffix(".preview.vtk")
        gmsh.write(str(vtu))

        # qualité : diamètre inscrit (ce qui pilote le dt rockim)
        _, _, conn = gmsh.model.mesh.getElements(2)
        tris = np.array(conn[0], dtype=int).reshape(-1, 3)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        idx = {t: i for i, t in enumerate(tags)}
        P = np.array(xyz).reshape(-1, 3)[
            [idx[t] for t in tris.flatten()]].reshape(-1, 3, 3)[:, :, :2]
        a = np.linalg.norm(P[:, 1] - P[:, 0], axis=1)
        b = np.linalg.norm(P[:, 2] - P[:, 1], axis=1)
        d = np.linalg.norm(P[:, 0] - P[:, 2], axis=1)
        s = 0.5 * (a + b + d)
        area = np.sqrt(np.maximum(s * (s - a) * (s - b) * (s - d), 0.0))
        h_ins = 4.0 * area / (a + b + d)        # 2 × rayon inscrit
        print(json.dumps({
            "ok": True,
            "triangles": int(tris.shape[0]),
            "nodes": int(len(tags)),
            "hMin": float(h_ins.min()),
            "hMed": float(np.median(h_ins)),
            "hMax": float(h_ins.max()),
            "msh": str(out_msh),
            "vtu": str(vtu),
        }))
        return 0
    except Exception as e:
        return _fail(f"gmsh : {e}")
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


if __name__ == "__main__":
    sys.exit(main())
