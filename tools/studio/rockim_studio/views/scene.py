"""scene.py — la vue 3D des résultats (spec 006, M1 / WP1.1-1.2).

QtInteractor de pyvistaqt : rendu VTK natif (GPU) dans un widget Qt. Le
widget se dégrade proprement si pyvista/pyvistaqt manquent (message +
studio utilisable sans 3D). Presets d'affichage :
  * bulk coloré par n'importe quel tableau cellule/point du VTU ;
  * « fissures » : joints ROMPUS en rouge par-dessus le bulk translucide
    (rupture = tBreak >= 0 si présent, sinon damage >= 1) ;
  * éléments érodés masqués quand le tableau `eroded` existe (fem/fem3d).
"""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLabel,
                               QSlider, QVBoxLayout, QWidget)

try:
    import numpy as np
    from pyvistaqt import QtInteractor

    from ..results.vtu_series import FrameSeries
    _HAVE_3D = True
    _IMPORT_ERROR = ""
except Exception as e:          # pyvista/pyvistaqt/OpenGL absents
    _HAVE_3D = False
    _IMPORT_ERROR = str(e)


class SceneView(QWidget):
    frame_changed = Signal(int, float)      # index, temps
    point_picked = Signal(float, float, float)   # sonde : nœud piqué

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)

        bar = QHBoxLayout()
        bar.addWidget(QLabel("Champ :"))
        self.array_box = QComboBox()
        self.array_box.currentTextChanged.connect(lambda _t: self.refresh())
        bar.addWidget(self.array_box, 1)
        self.cracks_box = QCheckBox("fissures")
        self.cracks_box.setChecked(True)
        self.cracks_box.toggled.connect(lambda _c: self.refresh())
        bar.addWidget(self.cracks_box)
        self.pick_box = QCheckBox("sonde au clic")
        self.pick_box.setToolTip("cliquer un nœud du maillage : ses "
                                 "coordonnées partent dans la Sonde nodale "
                                 "de l'onglet Courbes")
        self.pick_box.toggled.connect(self._toggle_pick)
        bar.addWidget(self.pick_box)
        lay.addLayout(bar)

        if _HAVE_3D:
            self.plotter = QtInteractor(self)
            lay.addWidget(self.plotter, 1)
        else:
            self.plotter = None
            lay.addWidget(QLabel(
                "Vue 3D indisponible — installer pyvista + pyvistaqt\n"
                f"({_IMPORT_ERROR})"), 1)

        row = QHBoxLayout()
        self.slider = QSlider(Qt.Horizontal)
        self.slider.valueChanged.connect(self._slide)
        row.addWidget(self.slider, 1)
        self.time_label = QLabel("—")
        row.addWidget(self.time_label)
        lay.addLayout(row)

        self.series = None
        self.index = 0
        self._camera_set = False

    # --- API --------------------------------------------------------------
    def load(self, out_dir: str):
        if not _HAVE_3D:
            return
        self.series = FrameSeries(out_dir)
        self.index = min(self.index, max(len(self.series) - 1, 0))
        if not len(self.series):
            self.time_label.setText("aucune frame VTU")
            return
        self.slider.blockSignals(True)
        self.slider.setRange(0, len(self.series) - 1)
        self.slider.setValue(len(self.series) - 1)   # dernière frame d'abord
        self.slider.blockSignals(False)
        self.index = len(self.series) - 1

        arrays = self.series.bulk_arrays()
        self.array_box.blockSignals(True)
        self.array_box.clear()
        self.array_box.addItems(arrays)
        for prefer in ("bulkD", "damage", "vonMises", "sigmaYY", "velocity"):
            if prefer in arrays:
                self.array_box.setCurrentText(prefer)
                break
        self.array_box.blockSignals(False)
        self._camera_set = False
        self.refresh()

    def _slide(self, value: int):
        self.index = value
        self.refresh()

    def _toggle_pick(self, on: bool):
        if self.plotter is None:
            return
        if on:
            self.plotter.enable_point_picking(
                callback=lambda p: self.point_picked.emit(
                    float(p[0]), float(p[1]), float(p[2])),
                show_message="clic gauche : piquer un nœud",
                left_clicking=True, show_point=True)
        else:
            self.plotter.disable_picking()

    def show_mesh(self, path: str):
        """Aperçu d'un maillage seul (sortie du mailleur M2) : filaire."""
        if not _HAVE_3D or self.plotter is None:
            return
        import pyvista as pv
        mesh = pv.read(path)
        self.series = None
        self.time_label.setText(f"maillage : {mesh.n_cells} cellules")
        self.plotter.clear()
        self.plotter.add_mesh(mesh, style="wireframe", color="#5a7d9a",
                              line_width=1)
        zmin, zmax = mesh.bounds[4], mesh.bounds[5]
        if zmax - zmin < 1e-12:
            self.plotter.view_xy()
        self.plotter.reset_camera()
        self._camera_set = False
        self.plotter.render()

    # --- rendu ------------------------------------------------------------
    def refresh(self):
        if self.series is None or not len(self.series):
            return
        t = self.series.times[self.index]
        self.time_label.setText(f"t = {t:.4e} s  "
                                f"[{self.index + 1}/{len(self.series)}]")
        self.frame_changed.emit(self.index, t)
        if self.plotter is None:
            return

        bulk = self.series.bulk(self.index)
        if "eroded" in bulk.cell_data:
            bulk = bulk.threshold(0.5, scalars="eroded", invert=True)

        cam = self.plotter.camera_position if self._camera_set else None
        self.plotter.clear()
        array = self.array_box.currentText() or None
        show_cracks = self.cracks_box.isChecked()
        self.plotter.add_mesh(
            bulk, scalars=array if array in
            (list(bulk.cell_data.keys()) + list(bulk.point_data.keys()))
            else None,
            cmap="viridis", opacity=0.35 if show_cracks else 1.0,
            show_edges=bulk.n_cells < 60000,
            scalar_bar_args={"title": array or ""})

        joints = self.series.joints(self.index) if show_cracks else None
        if joints is not None and joints.n_cells:
            broken = self._broken(joints)
            if broken is not None and broken.n_cells:
                self.plotter.add_mesh(broken, color="red", line_width=3,
                                      render_lines_as_tubes=joints.volume > 0)
        if cam is not None:
            self.plotter.camera_position = cam
        else:
            zmin, zmax = bulk.bounds[4], bulk.bounds[5]
            if zmax - zmin < 1e-12:      # cas 2D : vue à plat
                self.plotter.view_xy()
            self.plotter.reset_camera()
            self._camera_set = True
        self.plotter.render()

    @staticmethod
    def _broken(joints):
        if "tBreak" in joints.cell_data:
            mask = np.asarray(joints.cell_data["tBreak"]) >= 0.0
        elif "damage" in joints.cell_data:
            mask = np.asarray(joints.cell_data["damage"]) >= 1.0
        else:
            return None
        return joints.extract_cells(np.flatnonzero(mask)) \
            if mask.any() else None
