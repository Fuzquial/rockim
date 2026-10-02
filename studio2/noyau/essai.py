"""Modèle d'un essai triaxial 2D (ou de traction directe) et sa traduction en deck g1.

L'utilisateur pense en essai, pas en clés : ce module porte les choix physiques
(éprouvette, chargement, maillage, matériau, phases, discontinuités) et cache la
recette du solveur. La recette triaxiale de g1 n'a pas de scénario dédié :
`scenario = tension` + `loading = platens` + `pullV < 0` + `confiningPressure`
(DOCUMENTATION_rockim.md §5.7). Elle est figée ici, une fois.

Tout objet se sérialise en dictionnaire JSON (vers_dict / depuis_dict) pour la
file de calculs et l'interface.
"""
import dataclasses as dc
import math
from typing import List, Optional, Tuple

from . import cfg
from .geometrie import segments_plans

Segment = Tuple[float, float, float, float]


@dc.dataclass
class Eprouvette:
    W: float = 0.036               # m
    H: float = 0.072               # m
    epaisseur: float = 1.0         # m (2D : contrainte plane par unité d'épaisseur)


@dc.dataclass
class Chargement:
    """Triaxial : confinement latéral établi AVANT la charge axiale par plateaux."""
    type_essai: str = "triaxial"   # triaxial | traction
    sigma3_MPa: float = 20.0
    vitesse: float = 0.1           # m/s, vitesse de fermeture des plateaux (ou de traction)
    rampe_axiale: float = 2e-4     # s
    delai_axial: float = 6e-4      # s, fin de consolidation (pullDelay)
    rampe_confinement: float = 2e-4
    temps_jauge_confinement: float = 5e-4
    arret_apres_pic: bool = True
    chute_arret: float = 0.5       # arrêt quand sigma < (1 - x) pic
    delai_arret: float = 2e-4      # s
    jauge_bas: float = 0.25
    jauge_haut: float = 0.75
    delai_arret_traction: float = 5e-5


@dc.dataclass
class Maillage:
    """`voronoi` : le solveur génère grains et maillage (GBM possible).
    `gmsh` : maillage non structuré importé, une seule phase (g1 refuse phases + file)."""
    type: str = "voronoi"
    taille_grain: float = 0.003            # m
    taille_element: float = 0.0007         # m (Voronoï : grainElemSize ; Gmsh : h)
    dispersion_tailles: Optional[float] = None   # grainSizeSpread (Laguerre)
    graine: int = 4211
    semis: str = "random"
    jitter: float = 0.5
    lloyd: int = 2
    fichier_msh: Optional[str] = None      # gmsh : chemin du .msh produit par maillage.py


@dc.dataclass
class Materiau:
    rho: float = 2624.0
    E: float = 52e9
    nu: float = 0.21584
    ft: float = 34e6
    cohesion: float = 13.6e6
    frictionDeg: float = 13.4
    Gf: float = 70.0
    gfShearFactor: float = 10.0
    crushCap: Optional[float] = 400e6
    nom: str = "bohus"


@dc.dataclass
class Phase:
    nom: str
    fraction: float                        # SURFACIQUE
    rho: float
    E: float
    nu: float
    ft: Optional[float] = None             # None : héritée du matériau de volume
    cohesion: Optional[float] = None
    frictionDeg: Optional[float] = None
    Gf: Optional[float] = None
    taille_grain: Optional[float] = None   # affinité de taille (sur-spécifiée), m


@dc.dataclass
class Paire:
    """Surcharge des propriétés d'interface entre deux phases (gb.<a>.<b>.*)."""
    a: str
    b: str
    ft: float
    cohesion: float
    Gf: float
    frictionDeg: float


@dc.dataclass
class JointsGrain:
    alpha_ten: float = 0.5
    alpha_coh: float = 0.5
    alpha_gf: float = 0.5
    alpha_e: float = 1.0
    alpha_fric: float = 1.0
    paires: List[Paire] = dc.field(default_factory=list)


@dc.dataclass
class FamillePlans:
    """Famille de plans affaiblis, ancrée sur l'origine du maillage comme dans le solveur."""
    pendage: float                         # deg
    espacement: float                      # m
    facteur: float = 0.3                   # ft et cohésion des plans x facteur
    gf: str = "follow"
    fraction_rompue: float = 0.0           # part des plans déjà rompus au départ


@dc.dataclass
class Discontinuites:
    fraction_diffuse: float = 0.0          # jointPrebrokenFrac
    graine_diffuse: Optional[int] = None   # défaut : graine du maillage + 4242
    plans: Optional[FamillePlans] = None
    segments: List[Segment] = dc.field(default_factory=list)   # dessinés à la main, m
    mu_residuel: Optional[float] = None    # défaut : tan(frictionDeg), le frottement de PIC

    def actives(self):
        return bool(self.fraction_diffuse or self.plans or self.segments)


@dc.dataclass
class Schema:
    """Schéma numérique de la campagne : insertion adaptative + loi de Yan (Yan 2023)."""
    insertion: str = "adaptive"
    adoucissement: str = "yan"
    penalite_insertion: float = 4.0
    mu_contact: float = 0.1
    amortissement_local: float = 0.7
    xi_joint: float = 0.0


@dc.dataclass
class Sorties:
    T: float = 1e-2                        # s, borne haute (l'arrêt après pic coupe avant)
    frames: int = 24
    deformations_historique: bool = True   # historyStrains : epsAx, epsLat, epsVol
    champs_deformation: bool = True        # writeStrainFields : strainXX/YY/XY + déplacement (g1 >= J3)


@dc.dataclass
class Essai:
    nom: str = "essai"
    description: str = ""
    eprouvette: Eprouvette = dc.field(default_factory=Eprouvette)
    chargement: Chargement = dc.field(default_factory=Chargement)
    maillage: Maillage = dc.field(default_factory=Maillage)
    materiau: Materiau = dc.field(default_factory=Materiau)
    phases: List[Phase] = dc.field(default_factory=list)
    joints_grain: JointsGrain = dc.field(default_factory=JointsGrain)
    discontinuites: Discontinuites = dc.field(default_factory=Discontinuites)
    schema: Schema = dc.field(default_factory=Schema)
    sorties: Sorties = dc.field(default_factory=Sorties)

    # ------------------------------------------------------------ sérialisation
    def vers_dict(self):
        return dc.asdict(self)

    @staticmethod
    def depuis_dict(d):
        d = dict(d)
        disc = dict(d.get("discontinuites", {}))
        if disc.get("plans"):
            disc["plans"] = FamillePlans(**disc["plans"])
        disc["segments"] = [tuple(s) for s in disc.get("segments", [])]
        jg = dict(d.get("joints_grain", {}))
        jg["paires"] = [Paire(**p) for p in jg.get("paires", [])]
        return Essai(
            nom=d.get("nom", "essai"), description=d.get("description", ""),
            eprouvette=Eprouvette(**d.get("eprouvette", {})),
            chargement=Chargement(**d.get("chargement", {})),
            maillage=Maillage(**d.get("maillage", {})),
            materiau=Materiau(**d.get("materiau", {})),
            phases=[Phase(**p) for p in d.get("phases", [])],
            joints_grain=JointsGrain(**jg),
            discontinuites=Discontinuites(**disc),
            schema=Schema(**d.get("schema", {})),
            sorties=Sorties(**d.get("sorties", {})),
        )

    # ------------------------------------------------------------ dérivés
    def mu_residuel(self):
        m = self.discontinuites.mu_residuel
        if m is not None:
            return m
        # Le frottement résiduel d'un joint pré-rompu vaut le frottement de PIC :
        # une pré-fissure perd sa cohésion et rien d'autre. Au-delà, le frottement
        # AUGMENTERAIT en s'endommageant (défaut du run F6, 2026-09-15).
        return round(math.tan(math.radians(self.materiau.frictionDeg)), 4)

    def segments_prerompus(self):
        segs = []
        p = self.discontinuites.plans
        if p and p.fraction_rompue > 0:
            segs += segments_plans(p.pendage, p.espacement, p.fraction_rompue,
                                   self.eprouvette.W, self.eprouvette.H)
        return segs + list(self.discontinuites.segments)

    # ------------------------------------------------------------ deck
    def vers_cfg(self):
        e, c, m, mat = self.eprouvette, self.chargement, self.maillage, self.materiau
        traction = c.type_essai == "traction"
        titre = ["%s : %s" % (self.nom, "TRACTION DIRECTE" if traction else
                              "TRIAXIAL 2D, sigma3 = %g MPa" % c.sigma3_MPa)]
        titre += [l for l in self.description.strip().splitlines() if l.strip()]
        w = cfg.Ecrivain(titre)

        w.cle("mode", "fdem")
        w.cle("scenario", "tension")
        w.cle("loading", "grips" if traction else "platens")
        w.cle("T", self.sorties.T)
        w.cle("frames", self.sorties.frames)
        if self.sorties.deformations_historique:
            w.cle("historyStrains", True)
        if self.sorties.champs_deformation:
            w.cle("writeStrainFields", True, "exige un g1 compile apres le 2026-10-02 (spec 007 J3)")

        w.section("eprouvette")
        if m.type != "gmsh":           # en Gmsh, W et H sont recalculés depuis le maillage
            w.cle("W", e.W)
            w.cle("H", e.H)
        w.cle("thickness", e.epaisseur)

        w.section("maillage")
        if m.type == "gmsh":
            w.cle("mesh", "file")
            w.cle("meshFile", m.fichier_msh or "")
        else:
            w.cle("mesh", "voronoi")
            w.cle("grainSize", m.taille_grain)
            w.cle("grainSeeding", m.semis)
            w.cle("grainJitter", m.jitter)
            w.cle("lloydIters", m.lloyd)
            w.cle("grainMesh", "delaunay")
            w.cle("grainMeshRandom", True, "OBLIGATOIRE en GBM (DOC regle 8.4)")
            w.cle("grainElemSize", m.taille_element)
            if m.dispersion_tailles:
                w.cle("grainSizeSpread", m.dispersion_tailles)
        w.cle("seed", m.graine)

        w.section("materiau de volume (%s)" % mat.nom)
        for k in ("rho", "E", "nu", "ft", "cohesion", "frictionDeg", "Gf", "gfShearFactor"):
            w.cle(k, getattr(mat, k))
        if mat.crushCap:
            w.cle("crushCap", mat.crushCap)

        if self.phases:
            w.section("phases")
            w.cle("phases", " ".join(p.nom for p in self.phases))
            for p in self.phases:
                w.cle("phase.%s.fraction" % p.nom, p.fraction)
                if p.taille_grain:
                    w.cle("phase.%s.grainSize" % p.nom, p.taille_grain,
                          "CIBLE sur-specifiee, pas la taille realisee")
                for k in ("rho", "E", "nu", "ft", "cohesion", "frictionDeg", "Gf"):
                    v = getattr(p, k)
                    if v is not None:
                        w.cle("phase.%s.%s" % (p.nom, k), v)

        jg = self.joints_grain
        w.section("joints de grain")
        w.cle("gbAlphaTen", jg.alpha_ten)
        w.cle("gbAlphaCoh", jg.alpha_coh)
        w.cle("gbAlphaGf", jg.alpha_gf)
        w.cle("gbAlphaE", jg.alpha_e)
        w.cle("gbAlphaFric", jg.alpha_fric)
        for p in jg.paires:
            for k in ("ft", "cohesion", "Gf", "frictionDeg"):
                w.cle("gb.%s.%s.%s" % (p.a, p.b, k), getattr(p, k))

        d = self.discontinuites
        if d.actives():
            w.section("discontinuites preexistantes")
            w.cle("jointResidualMu", self.mu_residuel(), "frottement de PIC : la pre-fissure perd sa cohesion seule")
            w.cle("gcSurfaceRefresh", "eager", "exige des qu'il existe des joints pre-rompus")
            if d.fraction_diffuse:
                w.cle("jointPrebrokenFrac", d.fraction_diffuse)
                w.cle("jointPrebrokenSeed", d.graine_diffuse if d.graine_diffuse is not None
                      else m.graine + 4242)
            if d.plans:
                w.cle("weakPlanes", "%s %s" % (cfg.nombre(d.plans.pendage), cfg.nombre(d.plans.espacement)))
                w.cle("weakPlaneFactor", d.plans.facteur)
                w.cle("weakPlaneGf", d.plans.gf)
            segs = self.segments_prerompus()
            if segs:
                # Une SEULE clé : répétée, la seconde écraserait la première (DOC §5.11).
                # 6 chiffres significatifs, soit 0,1 µm sur une éprouvette de 72 mm.
                w.cle("preBrokenJoints", "; ".join("%.6g %.6g %.6g %.6g" % s for s in segs))

        s = self.schema
        w.section("schema numerique")
        w.cle("insertion", s.insertion)
        w.cle("jointSoftening", s.adoucissement)
        w.cle("insertionPenaltyFactor", s.penalite_insertion)
        w.cle("contactMu", s.mu_contact)
        w.cle("dampingLocal", s.amortissement_local)
        w.cle("jointXi", s.xi_joint)

        if traction:
            w.section("traction directe")
            w.cle("pullV", abs(c.vitesse), "positif = on tire")
            w.cle("pullRamp", c.rampe_axiale)
            w.cle("gripLateralFree", True)
            w.cle("verifyFt", False)
            w.cle("gripsStopAfterPeak", c.arret_apres_pic)
            w.cle("gripsStopDelay", c.delai_arret_traction)
        else:
            w.section("chargement et confinement")
            w.cle("pullV", -abs(c.vitesse), "negatif = plateaux en compression")
            w.cle("pullRamp", c.rampe_axiale)
            w.cle("pullDelay", c.delai_axial, "le confinement s'etablit AVANT l'axial")
            w.cle("confiningPressure", c.sigma3_MPa * 1e6)
            w.cle("confiningRamp", c.rampe_confinement)
            w.cle("confineFaces", "sides")
            w.cle("confineGaugeTime", c.temps_jauge_confinement)
            w.cle("ucsStopAfterPeak", c.arret_apres_pic)
            w.cle("stopPeakDrop", c.chute_arret)
            w.cle("ucsStopDelay", c.delai_arret)
            w.cle("gaugeLoFrac", c.jauge_bas)
            w.cle("gaugeHiFrac", c.jauge_haut)
        return w.texte()
