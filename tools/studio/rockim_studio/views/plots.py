"""plots.py — courbes : history.csv (live + post-run) et SONDE NODALE.

La sonde, c'est le « XY data from ODB » d'Abaqus : sur un run terminé,
choisir un nœud (au clic dans la scène ou par coordonnées), cocher des
variables, tracer leur évolution au fil des frames. matplotlib est
acceptable ici : des COURBES, jamais des champs maillés (la leçon de
l'ancien rockim_gui).
"""
from __future__ import annotations

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtWidgets import (QComboBox, QDoubleSpinBox, QGroupBox,
                               QHBoxLayout, QLabel, QListWidget,
                               QListWidgetItem, QPushButton, QVBoxLayout,
                               QWidget)
from PySide6.QtCore import Qt


class LivePlot(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Colonne :"))
        self.column = QComboBox()
        self.column.currentTextChanged.connect(lambda _t: self._redraw())
        bar.addWidget(self.column, 1)
        self.ref_label = QLabel("")
        bar.addWidget(self.ref_label)
        lay.addLayout(bar)

        self.canvas = FigureCanvasQTAgg(Figure(tight_layout=True))
        self.ax = self.canvas.figure.add_subplot(111)

        body = QHBoxLayout()
        body.addWidget(self.canvas, 1)
        body.addWidget(self._build_probe_box())
        lay.addLayout(body, 1)

        self.header: list[str] = []
        self.data: list[list] = []
        # run de référence superposé (tirets) : (nom, header, data)
        self.ref: tuple[str, list, list] | None = None

    def reset(self):
        self.header = []
        self.data = []
        self.column.clear()
        self.ax.clear()
        self.canvas.draw_idle()

    def set_reference(self, name: str | None, header=None, data=None):
        """Fixe (ou efface, name=None) le run de référence superposé."""
        self.ref = (name, header, data) if name else None
        self.ref_label.setText(f"réf : {name}" if name else "")
        self._redraw()

    def set_header(self, header: list):
        self.header = header
        self.column.blockSignals(True)
        self.column.clear()
        self.column.addItems(header[1:])   # colonne 0 = temps
        # présélection utile : force outil si présente
        for i, name in enumerate(header[1:]):
            if "force" in name.lower() or name in ("Fy", "Fz", "toolF"):
                self.column.setCurrentIndex(i)
                break
        self.column.blockSignals(False)
        self._redraw()

    def add_rows(self, rows: list):
        self.data.extend(rows)
        self._redraw()

    def load_csv(self, path):
        """Charge un history.csv complet (run terminé)."""
        header, data = self.read_csv(path)
        if header is None:
            return
        self.reset()
        self.data = data
        self.set_header(header)

    def _redraw(self):
        if not self.header or not self.data:
            return
        name = self.column.currentText()
        if name not in self.header:
            return
        j = self.header.index(name)
        t = [r[0] for r in self.data if len(r) > j]
        y = [r[j] for r in self.data if len(r) > j]
        # la sonde a pu passer la figure en sous-graphes : repartir d'un axe
        self.canvas.figure.clear()
        self.ax = self.canvas.figure.add_subplot(111)
        self.ax.plot(t, y, lw=1.0, label="courant")
        if self.ref is not None:
            ref_name, ref_header, ref_data = self.ref
            if name in ref_header:
                k = ref_header.index(name)
                rt = [r[0] for r in ref_data if len(r) > k]
                ry = [r[k] for r in ref_data if len(r) > k]
                self.ax.plot(rt, ry, lw=1.0, ls="--", label=ref_name)
                self.ax.legend(loc="best", fontsize=8)
        self.ax.set_xlabel(self.header[0])
        self.ax.set_ylabel(name)
        self.ax.grid(True, alpha=0.3)
        self.canvas.draw_idle()

    # --- sonde nodale (« XY data from ODB ») -------------------------------
    def _build_probe_box(self) -> QGroupBox:
        box = QGroupBox("Sonde nodale")
        v = QVBoxLayout(box)
        row = QHBoxLayout()
        self.px = QDoubleSpinBox()
        self.py = QDoubleSpinBox()
        self.pz = QDoubleSpinBox()
        for w, lbl in ((self.px, "x"), (self.py, "y"), (self.pz, "z")):
            w.setRange(-1e9, 1e9)
            w.setDecimals(4)
            row.addWidget(QLabel(lbl))
            row.addWidget(w)
        v.addLayout(row)
        self.node_label = QLabel("aucun nœud")
        self.node_label.setWordWrap(True)
        v.addWidget(self.node_label)
        self.var_list = QListWidget()
        self.var_list.setSelectionMode(QListWidget.NoSelection)
        v.addWidget(self.var_list, 1)
        self.trace_btn = QPushButton("Tracer au nœud")
        self.trace_btn.clicked.connect(self.trace_probe)
        v.addWidget(self.trace_btn)
        box.setMaximumWidth(240)
        self._probe = None
        self._probe_node = None
        return box

    def attach_series(self, series):
        """Branche la sonde sur un run chargé (FrameSeries)."""
        from ..results.probe import NodeProbe
        try:
            self._probe = NodeProbe(series)
        except Exception:
            self._probe = None
            return
        self.var_list.clear()
        for name in self._probe.variables():
            item = QListWidgetItem(name)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(
                Qt.Checked if name in ("u_y", "u_mag") else Qt.Unchecked)
            self.var_list.addItem(item)
        self.node_label.setText("entrer x/y (ou piquer dans la scène) "
                                "puis Tracer")

    def set_probe_point(self, x: float, y: float, z: float = 0.0):
        """Reçoit un point piqué dans la scène 3D."""
        self.px.setValue(x)
        self.py.setValue(y)
        self.pz.setValue(z)
        self._locate()

    def _locate(self):
        if self._probe is None:
            return None
        idx, (x, y, z) = self._probe.nearest_node(
            self.px.value(), self.py.value(), self.pz.value())
        self._probe_node = idx
        self.node_label.setText(f"nœud {idx} à ({x:.4g}, {y:.4g}"
                                + (f", {z:.4g})" if abs(z) > 1e-12 else ")"))
        return idx

    def trace_probe(self):
        if self._probe is None:
            self.node_label.setText("charger d'abord un run (Ctrl+R)")
            return
        idx = self._locate()
        names = [self.var_list.item(i).text()
                 for i in range(self.var_list.count())
                 if self.var_list.item(i).checkState() == Qt.Checked]
        if not names:
            self.node_label.setText("cocher au moins une variable")
            return
        data = self._probe.extract(idx, names)
        t = data["t"]
        fig = self.canvas.figure
        fig.clear()
        # une variable par sous-graphe, axe temps partagé (les unités
        # different : un axe commun écraserait les petites grandeurs)
        axes = fig.subplots(len(names), 1, sharex=True)
        if len(names) == 1:
            axes = [axes]
        for ax, nm in zip(axes, names):
            ax.plot(t, data[nm], marker="o", ms=3, lw=1.0)
            ax.set_ylabel(nm, fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=8)
        axes[0].set_title(f"nœud {idx}", fontsize=10)
        axes[-1].set_xlabel("t [s]")
        self.ax = axes[0]
        self.canvas.draw_idle()

    @staticmethod
    def read_csv(path):
        """Lit un history.csv -> (header, data) ou (None, None)."""
        import csv
        from pathlib import Path
        path = Path(path)
        if not path.exists():
            return None, None
        with open(path, newline="", encoding="utf-8",
                  errors="replace") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if not header:
                return None, None
            data = []
            for row in reader:
                try:
                    data.append([float(v) for v in row])
                except ValueError:
                    continue
        return header, data
