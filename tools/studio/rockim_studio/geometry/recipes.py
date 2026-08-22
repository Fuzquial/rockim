"""recipes.py — les recettes géométriques de M2 (spec 006, WP2.2).

Une recette = la description PARAMÉTRIQUE d'un domaine 2D à mailler :
plaque W×H percée d'une cavité (cercle, fer à cheval de l'étude tunnel,
ou polygone dessiné par l'utilisateur), avec champ de taille raffiné à la
paroi (hFine → hFar sur rFine, la transition douce des mailleurs maison).

Sérialisable en JSON : c'est le contrat entre l'UI et le worker gmsh
(processus séparé), et ce qui est archivé dans le projet .rsg pour
régénérer le maillage à l'identique (graine comprise).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field


@dataclass
class Cavity:
    kind: str = "circle"            # circle | horseshoe | polygon
    # circle : centre + rayon [m]
    cx: float = 50.0
    cy: float = 50.0
    r: float = 5.0
    # horseshoe : profil TUNNEL_HS de l'étude Wang et al., mis à l'échelle
    scale: float = 1.0
    # polygon : sommets dessinés, sens trigonométrique
    points: list = field(default_factory=list)   # [[x, y], ...]


@dataclass
class MeshRecipe:
    W: float = 100.0
    H: float = 100.0
    cavity: Cavity = field(default_factory=Cavity)
    hFine: float = 0.25             # taille à la paroi [m]
    rFine: float = 6.0              # rayon de la zone fine [m]
    hFar: float = 3.0               # taille au loin [m]
    algo2d: int = 5                 # 5 = Delaunay ISOTROPE (leçon 2026-08-17 :
    #                                 l'algo 6 fabrique un pavage quasi
    #                                 structuré qui canalise les fissures)
    seed: int = 1

    def validate(self) -> list[str]:
        errs = []
        if self.W <= 0 or self.H <= 0:
            errs.append("dimensions de plaque non positives")
        if not (0 < self.hFine <= self.hFar):
            errs.append("il faut 0 < hFine <= hFar")
        if self.rFine <= 0:
            errs.append("rFine doit être positif")
        c = self.cavity
        if c.kind == "circle":
            if c.r <= 0:
                errs.append("rayon de cavité non positif")
            if not (c.r < c.cx < self.W - c.r
                    and c.r < c.cy < self.H - c.r):
                errs.append("le cercle déborde de la plaque")
        elif c.kind == "polygon":
            if len(c.points) < 3:
                errs.append("polygone : au moins 3 sommets")
            for x, y in c.points:
                if not (0 < x < self.W and 0 < y < self.H):
                    errs.append("polygone : un sommet sort de la plaque")
                    break
        elif c.kind != "horseshoe":
            errs.append(f"cavité inconnue : {c.kind}")
        return errs

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "MeshRecipe":
        cav = Cavity(**d.get("cavity", {}))
        rest = {k: v for k, v in d.items() if k != "cavity"}
        return cls(cavity=cav, **rest)
