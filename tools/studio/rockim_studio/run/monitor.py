"""monitor.py — suivi live de history.csv pendant un run (spec 006, WP0.5).

S'appuie sur historyFlush (défaut true côté solveur) : le fichier grossit
ligne à ligne. Lecture incrémentale par offset (pas de relecture complète),
cadencée par QTimer — un watcher fichier est inutilement nerveux pour un
CSV qui s'étoffe toutes les ~n ms.
"""
from __future__ import annotations

import csv
import io
from pathlib import Path

from PySide6.QtCore import QObject, QTimer, Signal


class HistoryMonitor(QObject):
    header_ready = Signal(list)     # noms de colonnes
    rows_added = Signal(list)       # liste de lignes (listes de float/str)

    def __init__(self, parent=None, interval_ms: int = 500):
        super().__init__(parent)
        self.timer = QTimer(self)
        self.timer.setInterval(interval_ms)
        self.timer.timeout.connect(self.poll)
        self.path: Path | None = None
        self._offset = 0
        self._carry = ""
        self.header: list[str] = []

    def watch(self, out_dir: str | Path):
        self.path = Path(out_dir) / "history.csv"
        self._offset = 0
        self._carry = ""
        self.header = []
        self.timer.start()

    def stop(self):
        self.poll()                 # dernière moisson
        self.timer.stop()

    def poll(self):
        if self.path is None or not self.path.exists():
            return
        size = self.path.stat().st_size
        if size <= self._offset:
            return
        with open(self.path, "r", encoding="utf-8", errors="replace") as f:
            f.seek(self._offset)
            chunk = f.read()
            self._offset = f.tell()
        # une ligne encore en cours d'écriture reste en réserve
        text = self._carry + chunk
        complete, sep, self._carry = text.rpartition("\n")
        if not sep:
            self._carry = text
            return
        rows = []
        for row in csv.reader(io.StringIO(complete)):
            if not row:
                continue
            if not self.header:
                self.header = row
                self.header_ready.emit(row)
                continue
            rows.append([_maybe_float(v) for v in row])
        if rows:
            self.rows_added.emit(rows)


def _maybe_float(v: str):
    try:
        return float(v)
    except ValueError:
        return v
