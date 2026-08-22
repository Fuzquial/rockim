"""runner.py — lancement du solveur via QProcess (spec 006, WP0.5).

Règles maison respectées : UN run à la fois (file d'attente séquentielle),
arrêt par terminate() (jamais kill() en première intention), env
OMP_NUM_THREADS posé explicitement. Le cfg est écrit dans le dossier de
sortie avant lancement (traçabilité : chaque out_* garde sa config exacte).
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, Signal


class Runner(QObject):
    started = Signal(str)       # dossier de sortie
    output = Signal(str)        # stdout/stderr du solveur, ligne(s)
    finished = Signal(int, str)  # code retour, dossier de sortie

    def __init__(self, parent=None):
        super().__init__(parent)
        self.proc: QProcess | None = None
        self.exe = ""
        self.threads = 0        # 0 = laisser l'environnement décider
        self._queue: list[tuple[Path, Path]] = []
        self._current_out: Path | None = None

    @property
    def busy(self) -> bool:
        return self.proc is not None \
            and self.proc.state() != QProcess.NotRunning

    def launch(self, cfg_path: str | Path, out_dir: str | Path):
        self._queue.append((Path(cfg_path), Path(out_dir)))
        if not self.busy:
            self._next()

    def stop(self):
        self._queue.clear()
        if self.busy:
            self.output.emit("--- arrêt demandé (terminate) ---")
            self.proc.terminate()
            # kill de secours si le solveur ignore terminate sous 10 s
            if not self.proc.waitForFinished(10000):
                self.output.emit("--- pas de réponse, kill ---")
                self.proc.kill()

    def _next(self):
        if not self._queue:
            return
        cfg_path, out_dir = self._queue.pop(0)
        out_dir.mkdir(parents=True, exist_ok=True)
        self._current_out = out_dir

        self.proc = QProcess(self)
        env = QProcessEnvironment.systemEnvironment()
        if self.threads > 0:
            env.insert("OMP_NUM_THREADS", str(self.threads))
        self.proc.setProcessEnvironment(env)
        self.proc.setProcessChannelMode(QProcess.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self._read)
        self.proc.finished.connect(self._finished)
        self.proc.errorOccurred.connect(
            lambda e: self.output.emit(f"--- erreur QProcess : {e} ---"))
        self.proc.start(self.exe, [str(cfg_path), str(out_dir)])
        self.started.emit(str(out_dir))

    def _read(self):
        data = bytes(self.proc.readAllStandardOutput()).decode(
            "utf-8", errors="replace")
        if data:
            self.output.emit(data.rstrip("\n"))

    def _finished(self, code: int, _status):
        out = str(self._current_out)
        self.finished.emit(code, out)
        self.proc = None
        self._next()
