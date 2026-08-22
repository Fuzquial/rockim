"""tree.py — l'arbre du modèle : groupes UI -> clés explicites du cas.

Sélectionner un groupe affiche son formulaire complet dans le panneau de
propriétés ; les clés explicitement posées dans le cfg apparaissent comme
enfants (valeur en colonne 2), pour lire un cas d'un coup d'œil.
"""
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

from ..controller import Controller

# Ordre d'affichage des groupes (le reste vient après, trié).
_ORDER = ["Général", "Maillage", "Matériau", "Joints", "Contact", "Outil",
          "Conditions aux limites", "Corps et groupes", "Sorties",
          "Autres", "Inconnues"]


class ModelTree(QTreeWidget):
    group_selected = Signal(str)

    def __init__(self, ctrl: Controller, parent=None):
        super().__init__(parent)
        self.ctrl = ctrl
        self.setColumnCount(2)
        self.setHeaderLabels(["Modèle", "Valeur"])
        self.itemClicked.connect(self._clicked)
        ctrl.model_reset.connect(self.rebuild)
        ctrl.key_changed.connect(lambda _k: self.rebuild())
        self.rebuild()

    def rebuild(self):
        selected = self.currentItem().text(0) if self.currentItem() else None
        self.clear()
        model = self.ctrl.model
        groups = model.groups()
        names = [g for g in _ORDER if g in groups] + sorted(
            g for g in groups if g not in _ORDER)
        for name in names:
            explicit = [k for k in groups[name]
                        if model.is_explicit(k.name)]
            top = QTreeWidgetItem([name, f"{len(explicit)} posée(s)"
                                   if explicit else ""])
            self.addTopLevelItem(top)
            for k in explicit:
                QTreeWidgetItem(top, [k.name, model.cfg.pairs[k.name]])
            top.setExpanded(bool(explicit))
            if name == selected:
                self.setCurrentItem(top)
        self.resizeColumnToContents(0)

    def _clicked(self, item: QTreeWidgetItem, _col: int):
        top = item if item.parent() is None else item.parent()
        self.group_selected.emit(top.text(0))
