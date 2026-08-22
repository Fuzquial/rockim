"""gmsh_service.py — pilote Qt du worker gmsh (spec 006, WP2.1).

Écrit la recette en JSON, lance le worker en QProcess (jamais de gmsh dans
le processus GUI), publie le résultat par signaux. Un seul maillage à la
fois — un maillage relancé annule le précédent.
"""
from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QObject, QProcess, Signal

from .recipes import MeshRecipe

_WORKER = Path(__file__).resolve().parent / "worker.py"


class GmshService(QObject):
    started = Signal()
    finished = Signal(dict)         # le JSON du worker, ok True/False
    log = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.proc: QProcess | None = None

    @property
    def busy(self) -> bool:
        return self.proc is not None \
            and self.proc.state() != QProcess.NotRunning

    def mesh(self, recipe: MeshRecipe, out_msh: str | Path,
             python: str = "python3"):
        if self.busy:
            self.proc.kill()
        out_msh = Path(out_msh)
        out_msh.parent.mkdir(parents=True, exist_ok=True)
        recipe_path = out_msh.with_suffix(".recipe.json")
        recipe_path.write_text(json.dumps(recipe.to_dict(), indent=1),
                               encoding="utf-8")
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.SeparateChannels)
        self.proc.finished.connect(
            lambda _c, _s: self._done())
        self.proc.errorOccurred.connect(
            lambda e: self.finished.emit(
                {"ok": False, "error": f"worker : {e}"}))
        self.proc.start(python, [str(_WORKER), str(recipe_path),
                                 str(out_msh)])
        self.started.emit()

    def _done(self):
        out = bytes(self.proc.readAllStandardOutput()).decode(
            "utf-8", errors="replace").strip()
        err = bytes(self.proc.readAllStandardError()).decode(
            "utf-8", errors="replace").strip()
        if err:
            self.log.emit(err)
        try:
            result = json.loads(out.splitlines()[-1]) if out else \
                {"ok": False, "error": "worker muet"}
        except json.JSONDecodeError:
            result = {"ok": False, "error": f"sortie illisible : {out[:200]}"}
        self.proc = None
        self.finished.emit(result)
