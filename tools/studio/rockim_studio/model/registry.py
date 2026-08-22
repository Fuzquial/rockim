"""registry.py — le registre des clés de config (spec 006, §3.1).

Fusionne deux sources :
  * keys_extracted.json — la VÉRITÉ des sites d'appel de src/*.cpp, produite
    par tools/studio/dev/extract_keys.py (types, défauts par solveur, requis) ;
  * CURATED — métadonnées d'ergonomie tenues à la main (groupe UI, doc
    courte, bornes, énumérations), enrichies progressivement depuis
    DOCUMENTATION_rockim.md §5. Une clé extraite sans entrée CURATED reste
    parfaitement utilisable (formulaire générique) ; une entrée CURATED sans
    clé extraite est une ERREUR (clé fantôme) détectée par les tests.

Aucune dépendance hors stdlib : ce module doit s'importer sans Qt ni VTK.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_EXTRACTED = _HERE / "keys_extracted.json"

# Fichier source C++ -> mode(s) solveur (portée d'une clé).
_FILE_MODES = {
    "FemSolver.cpp": ("fem",),
    "Fem3dSolver.cpp": ("fem3d",),
    "DemSolver.cpp": ("dem",),
    "Dem3dSolver.cpp": ("dem3d",),
    "FdemSolver.cpp": ("fdem",),
    "Fdem3dSolver.cpp": ("fdem3d",),
    "MatLaw.cpp": ("fem3d", "fdem", "fdem3d"),
    "Tessellation.cpp": ("fdem",),
    "Tessellation3.cpp": ("fdem3d",),
    "main.cpp": ("fem", "fem3d", "dem", "dem3d", "fdem", "fdem3d"),
    "Config.cpp": (),
    "VtkWriter.cpp": (),
}

_ALL_MODES = ("fem", "fem3d", "dem", "dem3d", "fdem", "fdem3d")


def _modes_for(filename: str) -> tuple:
    # Les .hpp (Material, Tool, RandomField…) sont des blocs partagés lus par
    # plusieurs solveurs : portée « tous » faute de plus fin — un raffinement
    # par include réel est un enrichissement CURATED, pas un blocage.
    if filename.endswith(".hpp"):
        return _ALL_MODES
    return _FILE_MODES.get(filename, ())


@dataclass
class Key:
    name: str
    type: str                       # float | int | bool | str
    required: bool = False
    # défauts par mode : {"fdem": "0.05", ...} ; None = requis/sans défaut
    defaults: dict = field(default_factory=dict)
    scope: tuple = ()               # modes qui lisent la clé
    group: str = "Autres"           # groupe de l'arbre/formulaire UI
    doc: str = ""                   # colonne « rôle » de DOCUMENTATION §5
    bounds: tuple | None = None     # (min, max) inclusifs, numériques
    choices: tuple | None = None    # énumération pour les clés str
    unit: str = ""                  # unité SI affichée

    def default_for(self, mode: str):
        return self.defaults.get(mode)


@dataclass
class DynamicFamily:
    prefix: str
    doc: str = ""
    group: str = "Corps et groupes"

    def matches(self, name: str) -> bool:
        return name.startswith(self.prefix)


# --- Métadonnées d'ergonomie (extrait de démarrage ; complété par WP0.4). ---
# Format : nom -> dict de champs Key à surcharger.
CURATED = {
    "mode": dict(group="Général", choices=("fdem", "fdem3d", "fem", "fem3d",
                 "dem", "dem3d"),
                 doc="solveur — fdem/fdem3d = le cœur du studio ; "
                     "fem/dem = modes gelés (vérification seulement)"),
    "scenario": dict(group="Général", choices=("percussion", "shear",
                 "tension", "bar_wave", "brazilian", "shpb"),
                 doc="cas de charge"),
    "geometry": dict(group="Général", choices=("box", "disc", "cylinder",
                 "shpb"), doc="forme du domaine"),
    "mesh": dict(group="Maillage", choices=("grid", "voronoi", "file"),
                 doc="générateur de maillage"),
    "meshFile": dict(group="Maillage", doc="Gmsh MSH 2.2 ASCII (mesh=file)"),
    "T": dict(group="Général", unit="s", bounds=(0.0, None),
              doc="durée physique simulée"),
    "W": dict(group="Général", unit="m", bounds=(0.0, None), doc="largeur"),
    "H": dict(group="Général", unit="m", bounds=(0.0, None), doc="hauteur"),
    "D": dict(group="Général", unit="m", bounds=(0.0, None),
              doc="profondeur (3D)"),
    "rho": dict(group="Matériau", unit="kg/m³", bounds=(0.0, None),
                doc="masse volumique"),
    "E": dict(group="Matériau", unit="Pa", bounds=(0.0, None),
              doc="module d'Young"),
    "nu": dict(group="Matériau", bounds=(0.0, 0.4999),
               doc="coefficient de Poisson"),
    "ft": dict(group="Matériau", unit="Pa", bounds=(0.0, None),
               doc="résistance en traction"),
    "cohesion": dict(group="Matériau", unit="Pa", bounds=(0.0, None),
                     doc="cohésion"),
    "frictionDeg": dict(group="Matériau", unit="°", bounds=(0.0, 89.0),
                        doc="angle de frottement"),
    "Gf": dict(group="Matériau", unit="J/m²", bounds=(0.0, None),
               doc="énergie de fissuration mode I"),
    "jointXi": dict(group="Joints", bounds=(0.0, 1.0),
                    doc="amortissement du dashpot de joint — règle maison : "
                        "0 vérification, 0.01 quasi-statique, 0.05 impact"),
    "jointSoftening": dict(group="Joints",
                           choices=("linear", "yan", "munjiza"),
                           doc="adoucissement des joints cohésifs"),
    "insertion": dict(group="Joints", choices=("intrinsic", "adaptive"),
                      doc="insertion des joints (adaptive = Yan 2023)"),
    "contact": dict(group="Contact", choices=("penalty", "potential"),
                    doc="loi de contact général"),
    "gcActivation": dict(group="Contact", choices=("full", "adaptive"),
                         doc="activation des faces de contact (Fukuda)"),
    "absorbing": dict(group="Conditions aux limites",
                      choices=("none", "sides", "all"),
                      doc="frontières absorbantes Lysmer-Kuhlemeyer"),
    "law": dict(group="Matériau", choices=("elastic", "dpr", "mc", "saksala",
                "saksala2011", "dpdfh"), doc="loi de comportement du bulk"),
    "toolShape": dict(group="Outil", choices=("disc", "flat", "pdc", "sphere",
                      "none"), doc="forme de l'outil rigide"),
    "impactSpeed": dict(group="Outil", unit="m/s", doc="vitesse d'impact"),
    "toolMass": dict(group="Outil", unit="kg", bounds=(0.0, None),
                     doc="masse de l'outil (percussion = libre)"),
    "toolRadius": dict(group="Outil", unit="m", bounds=(0.0, None),
                       doc="rayon de l'outil"),
    "outputDir": dict(group="Sorties", doc="dossier de sortie"),
    "frames": dict(group="Sorties", bounds=(0, None),
                   doc="nombre de frames VTU"),
    "historyFlush": dict(group="Sorties",
                         doc="flush de history.csv à chaque ligne"),
    "seed": dict(group="Maillage", doc="graine (jitter, Voronoï, phases)"),
    "grainSize": dict(group="Maillage", unit="m", bounds=(0.0, None),
                      doc="diamètre moyen de grain (mesh=voronoi)"),
    # --- couplage hydro-mécanique (spec 004, AbuAisha) ---------------------
    "hydro": dict(group="Hydro", doc="active le couplage hydro-mécanique "
                  "(2D fdem, scénarios à cavité)"),
    "hydroSource": dict(group="Hydro", choices=("bore", "all"),
                        doc="faces sources : paroi du forage ou toute la "
                            "frontière"),
    "hydroInjection": dict(group="Hydro", choices=("rate", "pressure"),
                           doc="pompe à débit (pression = sortie) ou "
                               "pression imposée (contrôles)"),
    "hydroRate": dict(group="Hydro", unit="m³/s/m",
                      doc="débit de pompe — 20 l/s de l'article = 0.02"),
    "hydroPressure": dict(group="Hydro", unit="Pa",
                          doc="pression imposée (mode pressure)"),
    "hydroP0": dict(group="Hydro", unit="Pa",
                    doc="pression de référence (0 = pressions effectives)"),
    "hydroRamp": dict(group="Hydro", unit="s",
                      doc="rampe cosinus de la pompe"),
    "fluidBulk": dict(group="Hydro", unit="Pa", bounds=(0.0, None),
                      doc="module du fluide K_f — 2.2 GPa = eau ; fixe "
                          "toute la chronologie (hypothèse, pas une donnée "
                          "de l'article)"),
    "fluidDensity": dict(group="Hydro", unit="kg/m³", bounds=(0.0, None),
                         doc="masse volumique du fluide"),
    "boreCX": dict(group="Hydro", unit="m", doc="centre X du forage"),
    "boreCY": dict(group="Hydro", unit="m", doc="centre Y du forage"),
    "boreSelectR": dict(group="Hydro", unit="m",
                        doc="rayon de sélection des faces de forage "
                            "(R forage + une maille de garde)"),
    # --- état in-situ et excavation (tunnel EDZ) ---------------------------
    "insituSh": dict(group="Conditions aux limites", unit="Pa",
                     doc="contrainte in-situ HORIZONTALE (x) — attention : "
                         "le sigma_H d'AbuAisha (piège de nomenclature)"),
    "insituSv": dict(group="Conditions aux limites", unit="Pa",
                     doc="contrainte in-situ VERTICALE (y)"),
    "insituSxy": dict(group="Conditions aux limites", unit="Pa",
                      doc="cisaillement in-situ"),
    "excavRelease": dict(group="Conditions aux limites",
                         doc="relâchement d'excavation (tunnel)"),
    "excavStart": dict(group="Conditions aux limites", unit="s",
                       doc="début de l'excavation"),
    "excavRamp": dict(group="Conditions aux limites", unit="s",
                      doc="durée de la rampe d'excavation"),
    "confiningPressure": dict(group="Conditions aux limites", unit="Pa",
                              doc="pression de confinement / de paroi"),
    "confineFaces": dict(group="Conditions aux limites",
                         doc="faces recevant le confinement (bore, …)"),
    "confiningRamp": dict(group="Conditions aux limites", unit="s",
                          doc="rampe du confinement"),
}


class Registry:
    def __init__(self, keys: dict[str, Key],
                 dynamic: list[DynamicFamily]):
        self.keys = keys
        self.dynamic = dynamic

    @classmethod
    def load(cls, path: Path = _EXTRACTED) -> "Registry":
        data = json.loads(path.read_text(encoding="utf-8"))
        keys: dict[str, Key] = {}
        for name, info in data["keys"].items():
            defaults: dict[str, str] = {}
            scope: set[str] = set()
            for site in info["sites"]:
                for mode in _modes_for(site["file"]):
                    scope.add(mode)
                    if "default" in site and mode not in defaults:
                        defaults[mode] = site["default"]
            k = Key(name=name, type=info["type"],
                    required=info.get("required", False),
                    defaults=defaults, scope=tuple(sorted(scope)))
            for attr, val in CURATED.get(name, {}).items():
                setattr(k, attr, val)
            keys[name] = k
        dynamic = [DynamicFamily(prefix=p, doc=d)
                   for p, d in data["dynamic_prefixes"].items()]
        return cls(keys, dynamic)

    def lookup(self, name: str) -> Key | DynamicFamily | None:
        """Clé exacte, ou famille dynamique dont name porte le préfixe."""
        if name in self.keys:
            return self.keys[name]
        for fam in self.dynamic:
            if fam.matches(name):
                return fam
        return None

    def unknown_keys(self, names) -> list[str]:
        return [n for n in names if self.lookup(n) is None]

    def phantom_curated(self) -> list[str]:
        """Entrées CURATED sans clé extraite — erreurs à corriger."""
        return [n for n in CURATED if n not in self.keys]
