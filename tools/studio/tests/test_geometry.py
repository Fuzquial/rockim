"""Tests M2 : recettes, worker gmsh, et l'or — comparaison au mailleur CLI."""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/studio"))

from rockim_studio.geometry.recipes import Cavity, MeshRecipe  # noqa: E402

WORKER = ROOT / "tools/studio/rockim_studio/geometry/worker.py"


def _run_worker(recipe, out_msh):
    rp = out_msh.with_suffix(".recipe.json")
    rp.write_text(json.dumps(recipe.to_dict()), encoding="utf-8")
    r = subprocess.run([sys.executable, str(WORKER), str(rp), str(out_msh)],
                       capture_output=True, text=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def test_recipe_validation():
    r = MeshRecipe()
    assert r.validate() == []
    r.cavity = Cavity(kind="circle", cx=1.0, cy=50.0, r=5.0)
    assert any("déborde" in e for e in r.validate())
    r.cavity = Cavity(kind="polygon", points=[[10, 10], [20, 10]])
    assert any("3 sommets" in e for e in r.validate())
    d = MeshRecipe(cavity=Cavity(kind="horseshoe")).to_dict()
    assert MeshRecipe.from_dict(d).cavity.kind == "horseshoe"


def test_worker_circle_and_polygon():
    try:
        import gmsh  # noqa: F401
    except ImportError:
        print("  (gmsh absent : worker sauté)")
        return
    with tempfile.TemporaryDirectory() as td:
        res = _run_worker(MeshRecipe(hFine=1.0, rFine=6.0, hFar=5.0),
                          Path(td) / "c.msh")
        assert res["ok"], res
        assert res["triangles"] > 500 and Path(res["msh"]).exists()
        assert Path(res["vtu"]).exists()
        poly = MeshRecipe(hFine=1.0, rFine=6.0, hFar=5.0, cavity=Cavity(
            kind="polygon",
            points=[[45, 45], [55, 45], [58, 52], [50, 57], [42, 52]]))
        res2 = _run_worker(poly, Path(td) / "p.msh")
        assert res2["ok"], res2
        assert res2["triangles"] > 500


def test_worker_horseshoe_matches_cli_mesher():
    """L'OR de M2 : même profil, même graine, même algo -> le worker studio
    doit produire le MÊME maillage que make_unstructured_mesh.py tunnelhs
    (même nombre de triangles, mêmes h)."""
    try:
        import gmsh  # noqa: F401
    except ImportError:
        print("  (gmsh absent : test sauté)")
        return
    with tempfile.TemporaryDirectory() as td:
        recipe = MeshRecipe(W=100, H=100, hFine=0.25, rFine=6.0, hFar=3.0,
                            seed=1, cavity=Cavity(kind="horseshoe",
                                                  cx=50.0, cy=50.0))
        res = _run_worker(recipe, Path(td) / "hs.msh")
        assert res["ok"], res
        cli = subprocess.run(
            [sys.executable, str(ROOT / "tools/make_unstructured_mesh.py"),
             "tunnelhs", "100", "100", "0.25", "6", "3.0",
             str(Path(td) / "cli.msh"), "1"],
            capture_output=True, text=True)
        import re
        m = re.search(r"(\d+) triangles", cli.stdout)
        assert m, cli.stdout + cli.stderr
        n_cli = int(m.group(1))
        assert res["triangles"] == n_cli, \
            f"studio {res['triangles']} != CLI {n_cli}"


def _main():
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                print(f"FAIL {name}: {e}")
                fails += 1
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    _main()
