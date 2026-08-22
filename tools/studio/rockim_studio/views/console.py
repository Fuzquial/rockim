"""console.py — journal (log solveur/studio) + panneau de validation."""
from __future__ import annotations

from PySide6.QtWidgets import (QListWidget, QListWidgetItem, QPlainTextEdit,
                               QTabWidget)
from PySide6.QtGui import QColor

from ..controller import Controller


class Console(QTabWidget):
    def __init__(self, ctrl: Controller, parent=None):
        super().__init__(parent)
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(20000)
        self.issues = QListWidget()
        self.addTab(self.log_view, "Journal")
        self.addTab(self.issues, "Validation")
        ctrl.log.connect(self.append_log)
        ctrl.validation_changed.connect(self.show_issues)

    def append_log(self, text: str):
        self.log_view.appendPlainText(text.rstrip())

    def show_issues(self, issues: list):
        self.issues.clear()
        if not issues:
            self.issues.addItem("aucun problème détecté")
            self.setTabText(1, "Validation")
            return
        n_err = 0
        for level, key, msg in issues:
            item = QListWidgetItem(f"[{level}] {key} : {msg}")
            if level == "erreur":
                item.setForeground(QColor("#c62828"))
                n_err += 1
            else:
                item.setForeground(QColor("#b26a00"))
            self.issues.addItem(item)
        self.setTabText(1, f"Validation ({n_err} err. / "
                           f"{len(issues) - n_err} al.)")
