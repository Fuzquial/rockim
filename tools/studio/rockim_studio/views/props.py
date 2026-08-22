"""props.py — formulaire de propriétés auto-généré depuis le registre.

Un widget par clé selon son type : énumération -> QComboBox, bool ->
QCheckBox, numérique -> QLineEdit à validation différée (pas de
QDoubleSpinBox : il refuse la notation 50e9 et impose la locale). Les
valeurs par DÉFAUT s'affichent en placeholder grisé ; une valeur explicite
du cfg s'affiche en gras avec un bouton « ↺ » pour revenir au défaut.
"""
from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QHBoxLayout,
                               QLabel, QLineEdit, QScrollArea, QToolButton,
                               QWidget)

from ..controller import Controller
from ..model.registry import Key


class _KeyRow(QWidget):
    def __init__(self, ctrl: Controller, key: Key, parent=None):
        super().__init__(parent)
        self.ctrl, self.key = ctrl, key
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)

        model = ctrl.model
        explicit = model.is_explicit(key.name)
        value = model.cfg.pairs.get(key.name, "")
        default = key.default_for(model.mode)

        if key.choices:
            w = QComboBox()
            w.addItem("")           # vide = défaut solveur
            w.addItems(list(key.choices))
            if explicit:
                w.setCurrentText(value)
            w.currentTextChanged.connect(self._combo_changed)
        elif key.type == "bool":
            w = QCheckBox()
            w.setTristate(True)     # PartiallyChecked = défaut solveur
            if explicit:
                w.setCheckState(Qt.Checked if value.lower() in
                                ("1", "true", "yes", "on") else Qt.Unchecked)
            else:
                w.setCheckState(Qt.PartiallyChecked)
            w.checkStateChanged.connect(self._check_changed)
        else:
            w = QLineEdit()
            if explicit:
                w.setText(value)
            if default is not None:
                w.setPlaceholderText(f"défaut : {default}")
            w.editingFinished.connect(self._line_changed)
        self.widget = w
        self._set_bold(explicit)
        lay.addWidget(w, 1)

        if key.unit:
            lay.addWidget(QLabel(key.unit))
        self.reset_btn = QToolButton()
        self.reset_btn.setText("↺")
        self.reset_btn.setToolTip("revenir au défaut du solveur")
        self.reset_btn.setVisible(explicit)
        self.reset_btn.clicked.connect(self._reset)
        lay.addWidget(self.reset_btn)

        tip = key.doc or key.name
        if key.bounds:
            tip += f"  [{key.bounds[0]} … {key.bounds[1]}]"
        if key.required:
            tip += "  (REQUISE)"
        self.setToolTip(tip)

    def _set_bold(self, on: bool):
        f = self.widget.font()
        f.setBold(on)
        self.widget.setFont(f)
        if hasattr(self, "reset_btn"):
            self.reset_btn.setVisible(on)

    def _line_changed(self):
        text = self.widget.text().strip()
        if text:
            self.ctrl.set_key(self.key.name, text)
        else:
            self.ctrl.unset_key(self.key.name)
        self._set_bold(bool(text))

    def _combo_changed(self, text: str):
        if text:
            self.ctrl.set_key(self.key.name, text)
        else:
            self.ctrl.unset_key(self.key.name)
        self._set_bold(bool(text))

    def _check_changed(self, state):
        if state == Qt.PartiallyChecked:
            self.ctrl.unset_key(self.key.name)
            self._set_bold(False)
        else:
            self.ctrl.set_key(self.key.name,
                              "true" if state == Qt.Checked else "false")
            self._set_bold(True)

    def _reset(self):
        self.ctrl.unset_key(self.key.name)
        if isinstance(self.widget, QComboBox):
            self.widget.setCurrentIndex(0)
        elif isinstance(self.widget, QCheckBox):
            self.widget.setCheckState(Qt.PartiallyChecked)
        else:
            self.widget.clear()
        self._set_bold(False)


class PropertyPanel(QScrollArea):
    """Formulaire des clés d'UN groupe UI, régénéré à la sélection."""

    def __init__(self, ctrl: Controller, parent=None):
        super().__init__(parent)
        self.ctrl = ctrl
        self.setWidgetResizable(True)
        self._show([])

    def show_group(self, group: str):
        keys = self.ctrl.model.groups().get(group, [])
        self._show(keys)

    def _show(self, keys: list[Key]):
        body = QWidget()
        form = QFormLayout(body)
        form.setLabelAlignment(Qt.AlignRight)
        for key in keys:
            label = QLabel(key.name)
            if key.required:
                label.setText(key.name + " *")
            form.addRow(label, _KeyRow(self.ctrl, key))
        self.setWidget(body)
