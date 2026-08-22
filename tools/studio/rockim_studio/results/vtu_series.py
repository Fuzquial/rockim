"""vtu_series.py — découverte et lecture des frames VTU d'un dossier out_*.

Chaque solveur écrit ses familles de fichiers ; on les regroupe par frame :
  fem/fem3d/fdem/fdem3d : <mode>_XXXX.vtu (bulk) [+ <mode>_joints_XXXX.vtu]
  dem                   : dem_XXXX.vtu (points + liaisons en cellules ligne)
  dem3d                 : dem3d_particles_XXXX.vtu + dem3d_bonds_XXXX.vtu
Les temps viennent de frames.csv (frame,t,toolX,toolY…). Lecture par les
vrais lecteurs VTK (pyvista.read) avec cache borné — jamais de parsing
maison (la leçon du rockim_gui historique).
"""
from __future__ import annotations

import csv
import re
from collections import OrderedDict
from pathlib import Path

import pyvista as pv

_FRAME_RE = re.compile(r"^(?P<stem>.+?)_(?P<num>\d{3,5})\.vtu$")
# rôles par suffixe de famille : le « secondaire » est superposé au bulk
_SECONDARY = ("joints", "bonds")


class FrameSeries:
    """Une famille de frames (bulk + éventuel secondaire) d'un dossier."""

    def __init__(self, out_dir: str | Path, cache_frames: int = 8):
        self.dir = Path(out_dir)
        self.times: list[float] = []
        self.tool_xy: list[tuple[float, float]] = []
        self.primary: list[Path] = []
        self.secondary: list[Path] = []
        self._cache: OrderedDict[Path, pv.DataSet] = OrderedDict()
        self._cache_max = cache_frames
        self._discover()

    # --- découverte -------------------------------------------------------
    def _discover(self):
        families: dict[str, dict[int, Path]] = {}
        for f in self.dir.glob("*.vtu"):
            m = _FRAME_RE.match(f.name)
            if m:
                families.setdefault(m.group("stem"), {})[
                    int(m.group("num"))] = f
        if not families:
            return
        secondary_stems = [s for s in families
                           if s.rsplit("_", 1)[-1] in _SECONDARY]
        primary_stems = [s for s in families if s not in secondary_stems]
        # bulk = la famille primaire la plus peuplée
        stem = max(primary_stems, key=lambda s: len(families[s]),
                   default=None)
        if stem is None:
            return
        nums = sorted(families[stem])
        self.primary = [families[stem][n] for n in nums]
        for sec in secondary_stems:
            if sec.startswith(stem):
                self.secondary = [families[sec].get(n) for n in nums]
                break
        self._read_frames_csv(nums)

    def _read_frames_csv(self, nums: list[int]):
        path = self.dir / "frames.csv"
        times: dict[int, float] = {}
        tool: dict[int, tuple[float, float]] = {}
        if path.exists():
            with open(path, newline="", encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    try:
                        i = int(row["frame"])
                        times[i] = float(row["t"])
                        tool[i] = (float(row.get("toolX", 0) or 0),
                                   float(row.get("toolY", 0) or 0))
                    except (KeyError, ValueError):
                        continue
        self.times = [times.get(n, float(n)) for n in nums]
        self.tool_xy = [tool.get(n, (0.0, 0.0)) for n in nums]

    # --- accès ------------------------------------------------------------
    def __len__(self):
        return len(self.primary)

    def _read(self, path: Path) -> pv.DataSet:
        if path in self._cache:
            self._cache.move_to_end(path)
            return self._cache[path]
        mesh = pv.read(path)
        self._cache[path] = mesh
        while len(self._cache) > self._cache_max:
            self._cache.popitem(last=False)
        return mesh

    def bulk(self, i: int) -> pv.DataSet:
        return self._read(self.primary[i])

    def joints(self, i: int) -> pv.DataSet | None:
        if not self.secondary or self.secondary[i] is None:
            return None
        return self._read(self.secondary[i])

    def bulk_arrays(self) -> list[str]:
        if not self.primary:
            return []
        m = self.bulk(0)
        return list(m.cell_data.keys()) + list(m.point_data.keys())
