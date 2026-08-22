"""controller.py — LE médiateur (spec 006, §2).

Unique point de mutation du RockimModel. Les vues émettent des intentions
(appels de méthodes), le contrôleur mute le modèle via des commandes
(pile undo/redo dès M0 — API en place, spec §2) et notifie par signaux Qt.
Aucune vue ne touche le modèle directement.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QObject, Signal

from .model.rockim_model import RockimModel


class _SetKey:
    def __init__(self, model: RockimModel, key: str, value: str):
        self.model, self.key, self.value = model, key, value
        self.had = model.is_explicit(key)
        self.old = model.cfg.pairs.get(key)

    def do(self):
        self.model.set_value(self.key, self.value)

    def undo(self):
        if self.had:
            self.model.set_value(self.key, self.old)
        else:
            self.model.unset(self.key)


class _UnsetKey:
    def __init__(self, model: RockimModel, key: str):
        self.model, self.key = model, key
        self.had = model.is_explicit(key)
        self.old = model.cfg.pairs.get(key)

    def do(self):
        self.model.unset(self.key)

    def undo(self):
        if self.had:
            self.model.set_value(self.key, self.old)


class Controller(QObject):
    model_reset = Signal()          # nouveau fichier / ouverture
    key_changed = Signal(str)       # une clé a changé (nom)
    validation_changed = Signal(list)   # [(niveau, clé, message)]
    log = Signal(str)

    def __init__(self, model: RockimModel | None = None, parent=None):
        super().__init__(parent)
        self.model = model or RockimModel()
        self._undo: list = []
        self._redo: list = []

    # --- fichier ----------------------------------------------------------
    def open(self, path: str | Path):
        self.model.open(path)
        self._undo.clear()
        self._redo.clear()
        self.log.emit(f"ouvert : {path} "
                      f"({len(self.model.cfg.pairs)} clés, mode "
                      f"{self.model.mode})")
        self.model_reset.emit()
        self._revalidate()

    def new(self):
        self.model.new()
        self._undo.clear()
        self._redo.clear()
        self.log.emit("nouveau cas")
        self.model_reset.emit()
        self._revalidate()

    def save(self, path: str | Path | None = None) -> Path:
        target = self.model.save(path)
        self.log.emit(f"enregistré : {target}")
        return target

    # --- édition ----------------------------------------------------------
    def set_key(self, key: str, value: str):
        if self.model.cfg.pairs.get(key) == value.strip():
            return
        self._apply(_SetKey(self.model, key, value))
        self._after_edit(key)

    def unset_key(self, key: str):
        if not self.model.is_explicit(key):
            return
        self._apply(_UnsetKey(self.model, key))
        self._after_edit(key)

    def undo(self):
        if self._undo:
            cmd = self._undo.pop()
            cmd.undo()
            self._redo.append(cmd)
            self._after_edit(cmd.key, structure=True)

    def redo(self):
        if self._redo:
            cmd = self._redo.pop()
            cmd.do()
            self._undo.append(cmd)
            self._after_edit(cmd.key, structure=True)

    # --- interne ----------------------------------------------------------
    def _apply(self, cmd):
        cmd.do()
        self._undo.append(cmd)
        self._redo.clear()

    def _after_edit(self, key: str, structure: bool = False):
        # 'mode' change la portée de TOUTES les clés : reconstruction.
        if key == "mode" or structure:
            self.model_reset.emit()
        else:
            self.key_changed.emit(key)
        self._revalidate()

    def _revalidate(self):
        self.validation_changed.emit(self.model.validate())
