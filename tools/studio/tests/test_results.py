"""Tests M1 : FrameSeries sur un dossier out_* synthétique (sans solveur).

Fabrique avec pyvista un mini-dossier de résultats au format exact du
solveur (bulk fdem_XXXX.vtu + joints fdem_joints_XXXX.vtu + frames.csv) et
vérifie découverte, temps, tableaux et cache. Sauté si pyvista absent.
"""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/studio"))


def _make_out(dirpath: Path, n_frames: int = 3):
    import numpy as np
    import pyvista as pv
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float)
    tris = np.hstack([[3, 0, 1, 2], [3, 0, 2, 3]])
    lines = np.hstack([[2, 0, 2]])
    rows = ["frame,t,toolX,toolY"]
    for i in range(n_frames):
        bulk = pv.PolyData(pts, tris).cast_to_unstructured_grid()
        bulk.cell_data["vonMises"] = np.array([1.0 + i, 2.0 + i])
        bulk.save(dirpath / f"fdem_{i:04d}.vtu")
        joints = pv.PolyData(pts, lines=lines).cast_to_unstructured_grid()
        joints.cell_data["tBreak"] = np.array(
            [0.5 if i == n_frames - 1 else -1.0])
        joints.save(dirpath / f"fdem_joints_{i:04d}.vtu")
        rows.append(f"{i},{i * 1e-5},0,0")
    (dirpath / "frames.csv").write_text("\n".join(rows) + "\n")


def test_frame_series_discovery_and_read():
    try:
        import pyvista                              # noqa: F401
    except ImportError:
        print("  (pyvista absent : test sauté)")
        return
    from rockim_studio.results.vtu_series import FrameSeries
    with tempfile.TemporaryDirectory() as td:
        out = Path(td)
        _make_out(out)
        s = FrameSeries(out)
        assert len(s) == 3
        assert s.times == [0.0, 1e-5, 2e-5]
        assert "vonMises" in s.bulk_arrays()
        assert s.bulk(2).n_cells == 2
        j = s.joints(2)
        assert j is not None and float(j.cell_data["tBreak"][0]) == 0.5
        # cache : une relecture rend le même objet
        assert s.bulk(2) is s.bulk(2)


def test_scene_broken_mask():
    try:
        import numpy as np
        import pyvista as pv
    except ImportError:
        print("  (pyvista absent : test sauté)")
        return
    from rockim_studio.views.scene import SceneView
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    lines = np.hstack([[2, 0, 1], [2, 1, 2]])
    j = pv.PolyData(pts, lines=lines).cast_to_unstructured_grid()
    j.cell_data["tBreak"] = np.array([-1.0, 3e-5])
    broken = SceneView._broken(j)
    assert broken is not None and broken.n_cells == 1


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
