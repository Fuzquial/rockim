"""probe.py — la sonde nodale « à la ODB » (spec 006, M1+).

Sur un run TERMINÉ (FrameSeries) : choisir un nœud (au clic ou par
coordonnées), extraire l'évolution de plusieurs variables au fil des
frames, prêt à tracer. Variables servies :

  * u_x, u_y, u_z, u_mag — déplacement du nœud, reconstruit de x(t) − x(0)
    (les VTU stockent les positions COURANTES : l'identité du nœud est son
    index, stable d'une frame à l'autre — même maillage, mêmes points) ;
  * les tableaux aux points du VTU (velocity → v_x, v_y, v_z, v_mag) ;
  * les tableaux aux CELLULES de la maille adjacente au nœud (damage,
    sigmaXX, vonMises… — convention ODB : valeur de l'élément attaché).

Résolution temporelle = les frames écrites par le run (clé `frames`) —
pour un historique dense il faut la sonde côté solveur (spec §7 S4).
"""
from __future__ import annotations

import numpy as np

from .vtu_series import FrameSeries


class NodeProbe:
    def __init__(self, series: FrameSeries):
        if not len(series):
            raise ValueError("aucune frame dans ce dossier")
        self.series = series
        self._p0 = np.asarray(series.bulk(0).points)

    # --- sélection du nœud -------------------------------------------------
    def nearest_node(self, x: float, y: float, z: float = 0.0):
        """Index et coordonnées INITIALES du nœud le plus proche."""
        d2 = ((self._p0 - np.array([x, y, z])) ** 2).sum(axis=1)
        idx = int(np.argmin(d2))
        return idx, tuple(self._p0[idx])

    def _adjacent_cell(self, node: int) -> int | None:
        """Une cellule portant ce nœud (frame 0) — pour les champs cellule."""
        mesh = self.series.bulk(0)
        try:
            cells = mesh.point_cell_ids(node)
            return int(cells[0]) if len(cells) else None
        except Exception:
            return None

    # --- extraction --------------------------------------------------------
    def variables(self) -> list[str]:
        mesh = self.series.bulk(0)
        out = ["u_x", "u_y", "u_mag"]
        if self._p0.shape[1] > 2 and np.ptp(self._p0[:, 2]) > 1e-12:
            out.insert(2, "u_z")
        for name in mesh.point_data:
            arr = mesh.point_data[name]
            if arr.ndim == 2 and arr.shape[1] >= 2:
                out += [f"{name}_x", f"{name}_y", f"{name}_mag"]
                if arr.shape[1] > 2:
                    out.insert(len(out) - 1, f"{name}_z")
            else:
                out.append(name)
        out += list(mesh.cell_data.keys())
        return out

    def extract(self, node: int, names: list[str]) -> dict[str, np.ndarray]:
        """{'t': …, nom: série} pour chaque variable demandée, au nœud."""
        n = len(self.series)
        cell = None
        first = self.series.bulk(0)
        cell_names = set(first.cell_data.keys())
        if any(nm in cell_names for nm in names):
            cell = self._adjacent_cell(node)
        out: dict[str, list] = {nm: [] for nm in names}
        for i in range(n):
            mesh = self.series.bulk(i)
            p = np.asarray(mesh.points)[node]
            u = p - self._p0[node]
            for nm in names:
                if nm == "u_x":
                    out[nm].append(u[0])
                elif nm == "u_y":
                    out[nm].append(u[1])
                elif nm == "u_z":
                    out[nm].append(u[2])
                elif nm == "u_mag":
                    out[nm].append(float(np.linalg.norm(u)))
                elif nm in cell_names:
                    out[nm].append(float(mesh.cell_data[nm][cell])
                                   if cell is not None else np.nan)
                else:
                    base, _, comp = nm.rpartition("_")
                    if base and base in mesh.point_data \
                            and comp in ("x", "y", "z", "mag"):
                        vec = np.asarray(mesh.point_data[base][node])
                        out[nm].append(
                            float(np.linalg.norm(vec)) if comp == "mag"
                            else float(vec["xyz".index(comp)]))
                    elif nm in mesh.point_data:
                        out[nm].append(float(mesh.point_data[nm][node]))
                    else:
                        out[nm].append(np.nan)
        result = {"t": np.asarray(self.series.times)}
        result.update({nm: np.asarray(v) for nm, v in out.items()})
        return result
