"""cfg_io.py — lecture/écriture des fichiers de configuration rockim.

Reproduit EXACTEMENT la sémantique de Config::load (src/Config.cpp) :
  * '#' ouvre un commentaire jusqu'à la fin de ligne ;
  * lignes vides ou sans '=' ignorées ;
  * clé et valeur trimées ; une clé répétée écrase la précédente ;
  * aucune interprétation des valeurs à la lecture (le typage vient du
    registre, à la validation — comme dans le solveur où gets/getd décident).

La classe CfgFile préserve ce qu'un dict perdrait : l'ordre des clés et les
commentaires du fichier source (archivés, réémis en tête à l'écriture —
risque R7 de la spec 006). Round-trip garanti : parse(write(parse(f))) est
sémantiquement identique à parse(f) (mêmes paires clé-valeur).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class CfgFile:
    pairs: dict[str, str] = field(default_factory=dict)  # ordre d'insertion
    comments: list[str] = field(default_factory=list)    # texte sans '#'
    source: Path | None = None

    @classmethod
    def parse(cls, path: str | Path) -> "CfgFile":
        path = Path(path)
        cfg = cls(source=path)
        cfg.parse_text(path.read_text(encoding="utf-8", errors="replace"))
        return cfg

    def parse_text(self, text: str) -> "CfgFile":
        for raw in text.splitlines():
            line, hash_, comment = raw.partition("#")
            if hash_ and comment.strip():
                self.comments.append(comment.strip())
            line = line.strip()
            if not line:
                continue
            key, eq, value = line.partition("=")
            if not eq:
                continue
            key = key.strip()
            # dict Python : l'affectation répétée écrase et CONSERVE la
            # position d'origine ; std::map écrase aussi — même sémantique.
            self.pairs[key.strip()] = value.strip()
        return self

    def write(self, path: str | Path, header: str = "") -> None:
        Path(path).write_text(self.dumps(header), encoding="utf-8")

    def dumps(self, header: str = "") -> str:
        out = []
        if header:
            out += [f"# {line}" for line in header.splitlines()]
        if self.comments:
            out.append("# --- commentaires du fichier source ---")
            out += [f"# {c}" for c in self.comments]
        if out:
            out.append("")
        width = max((len(k) for k in self.pairs), default=0)
        out += [f"{k.ljust(width)} = {v}" for k, v in self.pairs.items()]
        return "\n".join(out) + "\n"
