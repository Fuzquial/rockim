"""geometry.py — la vue Géométrie de M2 : dessiner la cavité, mailler.

Gauche : formulaire de la recette (plaque, type de cavité, champ de
taille). Centre : le CANEVAS — la plaque à l'échelle, la cavité dessinée ;
en mode polygone, CHAQUE CLIC AJOUTE UN SOMMET (clic droit : retirer le
dernier). Bouton « Mailler » → worker gmsh hors processus → stats +
meshFile posé dans le modèle + aperçu chargé dans la scène 3D.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import (QComboBox, QDoubleSpinBox, QFormLayout,
                               QHBoxLayout, QLabel, QPushButton, QSpinBox,
                               QVBoxLayout, QWidget)

from ..controller import Controller
from ..geometry.gmsh_service import GmshService
from ..geometry.recipes import Cavity, MeshRecipe


class _Canvas(QWidget):
    """La plaque vue de face, à l'échelle. En mode polygone les clics
    dessinent les sommets (coordonnées PHYSIQUES, converties)."""

    points_changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.recipe = MeshRecipe()
        self.setMinimumSize(320, 320)

    # physique -> pixels (marge 12 px, y vers le haut)
    def _tf(self):
        m = 12
        w, h = self.width() - 2 * m, self.height() - 2 * m
        s = min(w / self.recipe.W, h / self.recipe.H)
        ox = m + (w - s * self.recipe.W) / 2
        oy = m + (h - s * self.recipe.H) / 2
        return s, ox, oy

    def _to_px(self, x, y):
        s, ox, oy = self._tf()
        return QPointF(ox + s * x, oy + s * (self.recipe.H - y))

    def _to_phys(self, pos):
        s, ox, oy = self._tf()
        return ((pos.x() - ox) / s,
                self.recipe.H - (pos.y() - oy) / s)

    def mousePressEvent(self, ev):
        if self.recipe.cavity.kind != "polygon":
            return
        if ev.button() == Qt.LeftButton:
            x, y = self._to_phys(ev.position())
            if 0 < x < self.recipe.W and 0 < y < self.recipe.H:
                self.recipe.cavity.points.append([round(x, 3), round(y, 3)])
        elif ev.button() == Qt.RightButton and self.recipe.cavity.points:
            self.recipe.cavity.points.pop()
        self.points_changed.emit()
        self.update()

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.fillRect(self.rect(), self.palette().base())
        s, ox, oy = self._tf()
        plate = QPolygonF([self._to_px(0, 0), self._to_px(self.recipe.W, 0),
                           self._to_px(self.recipe.W, self.recipe.H),
                           self._to_px(0, self.recipe.H)])
        p.setPen(QPen(QColor("#666"), 1.5))
        p.setBrush(QBrush(QColor(160, 160, 170, 60)))
        p.drawPolygon(plate)

        c = self.recipe.cavity
        p.setBrush(self.palette().base())
        p.setPen(QPen(QColor("#c62828"), 2))
        if c.kind == "circle":
            center = self._to_px(c.cx, c.cy)
            p.drawEllipse(center, s * c.r, s * c.r)
        elif c.kind == "horseshoe":
            # silhouette approchée pour l'aperçu (voûte + piédroits) ; la
            # vraie géométrie OCC est celle du worker (build_tunnel_hs)
            hs, ht = 5.55 * c.scale, 8.85 * c.scale
            rc = 5.55 * c.scale
            y0 = c.cy - ht / 2
            path = QPolygonF()
            import math
            for i in range(33):
                a = math.pi * i / 32
                path.append(self._to_px(c.cx + hs * math.cos(a),
                                        y0 + (ht - rc) + rc * math.sin(a)))
            path.append(self._to_px(c.cx - hs, y0))
            path.append(self._to_px(c.cx + hs, y0))
            p.drawPolygon(path)
        else:
            pts = [self._to_px(x, y) for x, y in c.points]
            if len(pts) >= 2:
                p.drawPolyline(QPolygonF(pts + pts[:1] if len(pts) > 2
                                         else pts))
            p.setBrush(QBrush(QColor("#c62828")))
            for q in pts:
                p.drawEllipse(q, 3, 3)
            if not pts:
                p.setPen(QPen(QColor("#888")))
                p.drawText(self.rect(), Qt.AlignCenter,
                           "cliquer pour dessiner la cavité\n"
                           "(clic droit : retirer le dernier sommet)")
        # anneau de la zone fine
        if c.kind != "polygon" or len(c.points) >= 3:
            p.setPen(QPen(QColor("#2a7d46"), 1, Qt.DashLine))
            p.setBrush(Qt.NoBrush)
            if c.kind == "circle":
                center = self._to_px(c.cx, c.cy)
                p.drawEllipse(center, s * (c.r + self.recipe.rFine),
                              s * (c.r + self.recipe.rFine))


class GeometryPanel(QWidget):
    mesh_ready = Signal(dict)       # résultat du worker, après meshFile posé

    def __init__(self, ctrl: Controller, parent=None):
        super().__init__(parent)
        self.ctrl = ctrl
        self.service = GmshService(self)
        self.service.finished.connect(self._meshed)
        self.service.log.connect(lambda t: ctrl.log.emit(t))

        lay = QHBoxLayout(self)
        form_box = QVBoxLayout()
        form = QFormLayout()

        def dspin(val, lo, hi, step=1.0, dec=3):
            w = QDoubleSpinBox()
            w.setRange(lo, hi)
            w.setDecimals(dec)
            w.setSingleStep(step)
            w.setValue(val)
            w.valueChanged.connect(self._sync)
            return w

        self.W = dspin(100.0, 0.01, 10000)
        self.H = dspin(100.0, 0.01, 10000)
        self.kind = QComboBox()
        self.kind.addItems(["circle", "horseshoe", "polygon"])
        self.kind.currentTextChanged.connect(self._sync)
        self.cx = dspin(50.0, 0, 10000)
        self.cy = dspin(50.0, 0, 10000)
        self.r = dspin(5.0, 0.001, 5000)
        self.scale = dspin(1.0, 0.01, 100, 0.1)
        self.hFine = dspin(0.25, 1e-6, 1000, 0.05)
        self.rFine = dspin(6.0, 1e-6, 5000)
        self.hFar = dspin(3.0, 1e-6, 1000, 0.5)
        self.seed = QSpinBox()
        self.seed.setRange(1, 10 ** 9)
        self.seed.setValue(1)
        self.seed.valueChanged.connect(self._sync)

        form.addRow("W [m]", self.W)
        form.addRow("H [m]", self.H)
        form.addRow("cavité", self.kind)
        form.addRow("centre x [m]", self.cx)
        form.addRow("centre y [m]", self.cy)
        form.addRow("rayon [m]", self.r)
        form.addRow("échelle profil", self.scale)
        form.addRow("hFine [m]", self.hFine)
        form.addRow("rFine [m]", self.rFine)
        form.addRow("hFar [m]", self.hFar)
        form.addRow("graine", self.seed)
        form_box.addLayout(form)

        self.mesh_btn = QPushButton("Mailler")
        self.mesh_btn.clicked.connect(self.do_mesh)
        form_box.addWidget(self.mesh_btn)
        self.status = QLabel("—")
        self.status.setWordWrap(True)
        form_box.addWidget(self.status)
        form_box.addStretch(1)
        lay.addLayout(form_box)

        self.canvas = _Canvas()
        self.canvas.points_changed.connect(self._sync)
        lay.addWidget(self.canvas, 1)
        self._sync()

    # --- recette <- widgets ------------------------------------------------
    def recipe(self) -> MeshRecipe:
        c = self.canvas.recipe.cavity
        return MeshRecipe(
            W=self.W.value(), H=self.H.value(),
            cavity=Cavity(kind=self.kind.currentText(),
                          cx=self.cx.value(), cy=self.cy.value(),
                          r=self.r.value(), scale=self.scale.value(),
                          points=list(c.points)),
            hFine=self.hFine.value(), rFine=self.rFine.value(),
            hFar=self.hFar.value(), seed=self.seed.value())

    def _sync(self, *_a):
        kind = self.kind.currentText()
        self.cx.setEnabled(kind != "polygon")
        self.cy.setEnabled(kind != "polygon")
        self.r.setEnabled(kind == "circle")
        self.scale.setEnabled(kind == "horseshoe")
        old_pts = self.canvas.recipe.cavity.points
        self.canvas.recipe = self.recipe()
        self.canvas.recipe.cavity.points = old_pts
        errs = self.canvas.recipe.validate()
        self.status.setText("prêt" if not errs else " ; ".join(errs))
        self.mesh_btn.setEnabled(not errs and not self.service.busy)
        self.canvas.update()

    # --- maillage ----------------------------------------------------------
    def do_mesh(self, out_msh: str | None = None):
        recipe = self.canvas.recipe
        if recipe.validate():
            return
        if out_msh is None:
            base = self.ctrl.model.source_dir or Path.cwd()
            out_msh = Path(base) / "meshes/studio_mesh.msh"
        self._out_msh = Path(out_msh)
        self.status.setText("maillage en cours…")
        self.mesh_btn.setEnabled(False)
        self.service.mesh(recipe, self._out_msh)

    def _meshed(self, result: dict):
        self.mesh_btn.setEnabled(True)
        if not result.get("ok"):
            self.status.setText(f"ÉCHEC : {result.get('error')}")
            self.ctrl.log.emit(f"maillage : {result.get('error')}")
            return
        self.status.setText(
            f"{result['triangles']} triangles — h inscrit "
            f"min {result['hMin']:.4g} / méd {result['hMed']:.4g} m\n"
            f"-> {result['msh']}")
        self.ctrl.set_key("mesh", "file")
        self.ctrl.set_key("meshFile", result["msh"])
        self.ctrl.log.emit(
            f"maillé : {result['triangles']} triangles, meshFile posé")
        self.mesh_ready.emit(result)
