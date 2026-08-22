"""rockim_model.py — l'état complet d'un cas (spec 006, §3.2), sans Qt.

Enveloppe CfgFile + Registry : valeurs explicites du .cfg, défauts du mode
courant, groupes UI, validation. Le contrôleur (Qt) mute ce modèle et
notifie les vues ; ce module doit rester importable sans PySide6.
"""
from __future__ import annotations

from pathlib import Path

from .cfg_io import CfgFile
from .registry import DynamicFamily, Key, Registry


class RockimModel:
    def __init__(self, registry: Registry | None = None):
        self.registry = registry or Registry.load()
        self.cfg = CfgFile()
        self.path: Path | None = None
        self.dirty = False

    # --- cycle de vie -----------------------------------------------------
    def open(self, path: str | Path) -> None:
        self.cfg = CfgFile.parse(path)
        self.path = Path(path)
        self.dirty = False

    def new(self) -> None:
        self.cfg = CfgFile()
        # le studio est l'interface du FDEM : un cas neuf démarre en fdem
        self.cfg.pairs["mode"] = "fdem"
        self.path = None
        self.dirty = False

    def open_template(self, path: str | Path) -> None:
        """Ouvre une config de référence comme GABARIT : contenu chargé,
        mais aucun chemin — l'enregistrement passera par « sous… » pour ne
        jamais écraser la référence."""
        self.cfg = CfgFile.parse(path)
        self.path = None
        self.dirty = True

    def save(self, path: str | Path | None = None) -> Path:
        target = Path(path) if path else self.path
        if target is None:
            raise ValueError("aucun chemin de sauvegarde")
        self.cfg.write(target, header="écrit par rockim-studio")
        self.path = target
        self.dirty = False
        return target

    # --- accès aux clés ---------------------------------------------------
    @property
    def mode(self) -> str:
        return self.cfg.pairs.get("mode", "fem")

    def value(self, key: str) -> str | None:
        """Valeur explicite du cfg, sinon défaut du mode courant, sinon None."""
        if key in self.cfg.pairs:
            return self.cfg.pairs[key]
        k = self.registry.keys.get(key)
        if k is not None:
            return k.default_for(self.mode)
        return None

    def is_explicit(self, key: str) -> bool:
        return key in self.cfg.pairs

    def set_value(self, key: str, value: str) -> None:
        self.cfg.pairs[key] = value.strip()
        self.dirty = True

    def unset(self, key: str) -> None:
        if key in self.cfg.pairs:
            del self.cfg.pairs[key]
            self.dirty = True

    # --- structure pour l'UI ---------------------------------------------
    def groups(self) -> dict[str, list[Key]]:
        """Groupe UI -> clés visibles pour le mode courant (portée), les
        clés explicites du cfg TOUJOURS incluses même hors portée (elles
        doivent rester éditables/supprimables)."""
        mode = self.mode
        out: dict[str, list[Key]] = {}
        for k in self.registry.keys.values():
            visible = (not k.scope) or (mode in k.scope) \
                or (k.name in self.cfg.pairs)
            if visible:
                out.setdefault(k.group, []).append(k)
        for grp in out.values():
            grp.sort(key=lambda k: k.name.lower())
        # clés du cfg inconnues du registre statique : familles dynamiques
        # ou vraies inconnues — groupe dédié pour rester visibles.
        for name in self.cfg.pairs:
            hit = self.registry.lookup(name)
            if isinstance(hit, DynamicFamily):
                out.setdefault(hit.group, []).append(
                    Key(name=name, type="str", group=hit.group, doc=hit.doc))
            elif hit is None:
                out.setdefault("Inconnues", []).append(
                    Key(name=name, type="str", group="Inconnues",
                        doc="clé absente du registre — vérifier l'orthographe"))
        return out

    # --- validation (noyau M0 : types, bornes, énumérations, locale) -----
    def validate(self) -> list[tuple[str, str, str]]:
        """Retourne [(niveau, clé, message)] ; niveau = 'erreur'|'alerte'."""
        issues: list[tuple[str, str, str]] = []
        for name, raw in self.cfg.pairs.items():
            hit = self.registry.lookup(name)
            if hit is None:
                issues.append(("alerte", name, "clé inconnue du registre"))
                continue
            if isinstance(hit, DynamicFamily):
                continue
            if hit.type in ("float", "int"):
                if "," in raw:
                    issues.append(("erreur", name,
                                   f"virgule décimale dans '{raw}' — le "
                                   "solveur exige le point (locale FR)"))
                    continue
                try:
                    num = float(raw)
                except ValueError:
                    issues.append(("erreur", name,
                                   f"valeur non numérique '{raw}'"))
                    continue
                if hit.type == "int" and num != int(num):
                    issues.append(("erreur", name,
                                   f"entier attendu, reçu '{raw}'"))
                if hit.bounds:
                    lo, hi = hit.bounds
                    if lo is not None and num < lo:
                        issues.append(("erreur", name,
                                       f"{raw} < borne min {lo}"))
                    if hi is not None and num > hi:
                        issues.append(("erreur", name,
                                       f"{raw} > borne max {hi}"))
            elif hit.choices and raw not in hit.choices:
                issues.append(("erreur", name,
                               f"'{raw}' hors énumération "
                               f"{'/'.join(hit.choices)}"))
        return issues
