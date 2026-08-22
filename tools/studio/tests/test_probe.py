"""Tests de la sonde nodale (« XY from ODB ») sur un out_* synthétique."""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/studio"))


def test_probe_extract():
    try:
        import numpy as np
        import pyvista as pv
    except ImportError:
        print("  (pyvista absent : test sauté)")
        return
    from rockim_studio.results.probe import NodeProbe
    from rockim_studio.results.vtu_series import FrameSeries
    pts0 = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float)
    tris = np.hstack([[3, 0, 1, 2], [3, 0, 2, 3]])
    with tempfile.TemporaryDirectory() as td:
        rows = ["frame,t,toolX,toolY"]
        for i in range(4):
            pts = pts0.copy()
            pts[2] += [0, 0.01 * i, 0]          # le nœud 2 monte
            m = pv.PolyData(pts, tris).cast_to_unstructured_grid()
            m.cell_data["damage"] = np.array([0.1 * i, 0.0])
            m.point_data["velocity"] = np.tile([0.0, 1.0 + i, 0.0], (4, 1))
            m.save(Path(td) / f"fdem_{i:04d}.vtu")
            rows.append(f"{i},{i * 1e-3},0,0")
        (Path(td) / "frames.csv").write_text("\n".join(rows) + "\n")
        probe = NodeProbe(FrameSeries(td))
        idx, coords = probe.nearest_node(1.02, 0.98)
        assert idx == 2 and coords[0] == 1.0
        names = probe.variables()
        assert {"u_y", "u_mag", "velocity_y", "damage"} <= set(names)
        data = probe.extract(idx, ["u_y", "velocity_y", "damage"])
        assert np.allclose(data["u_y"], [0, 0.01, 0.02, 0.03])
        assert np.allclose(data["velocity_y"], [1, 2, 3, 4])
        assert np.allclose(data["damage"], [0, 0.1, 0.2, 0.3])
        assert np.allclose(data["t"], [0, 1e-3, 2e-3, 3e-3])


def _main():
    try:
        test_probe_extract()
        print("PASS test_probe_extract")
    except AssertionError as e:
        print(f"FAIL test_probe_extract: {e}")
        sys.exit(1)


if __name__ == "__main__":
    _main()
