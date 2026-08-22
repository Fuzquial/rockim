"""plots.py — courbe live de history.csv (M0 ; les vraies planches sont M1).

matplotlib est acceptable ici : des COURBES, jamais des champs maillés (la
leçon de l'ancien rockim_gui). Redessin throttlé par les arrivées du
monitor (déjà cadencées à 500 ms).
"""
from __future__ import annotations

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QVBoxLayout, \
    QWidget


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
        lay.addWidget(self.canvas, 1)

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
        self.ax.clear()
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
