"""Vérifications avant lancement : ce que le solveur refuserait, et les règles maison.

Chaque avis a un niveau :
  erreur  le solveur refuserait le deck, ou l'essai n'a pas de sens ; on ne lance pas ;
  alerte  le deck tourne mais une règle apprise à nos dépens est enfreinte ;
  info    une grandeur utile à lire avant de lancer.
Chaque règle cite son origine. Les tests (tests/test_validation.py) font échouer
chacune d'elles sur un cas construit pour.
"""
import dataclasses as dc
import math

from .geometrie import aire, elements_par_grain, nombre_elements

ELEMENTS_PAR_GRAIN = (35.0, 90.0)       # règle du 2026-09-07 (60 runs perdus avec `fan`)


@dc.dataclass
class Avis:
    niveau: str
    code: str
    message: str


def verifier(essai):
    A = []

    def erreur(code, msg):
        A.append(Avis("erreur", code, msg))

    def alerte(code, msg):
        A.append(Avis("alerte", code, msg))

    def info(code, msg):
        A.append(Avis("info", code, msg))

    e, c, m, mat = essai.eprouvette, essai.chargement, essai.maillage, essai.materiau
    ph, d, jg = essai.phases, essai.discontinuites, essai.joints_grain
    noms = [p.nom for p in ph]

    # ---- géométrie et maillage
    if min(e.W, e.H, e.epaisseur) <= 0:
        erreur("dimensions", "W, H et l'épaisseur doivent être positives.")
    if m.type not in ("voronoi", "gmsh"):
        erreur("maillage", "type de maillage inconnu : %s" % m.type)
    if m.taille_element <= 0:
        erreur("taille_element", "la taille d'élément doit être positive.")
    if m.type == "gmsh":
        if ph:
            erreur("gmsh_phases", "g1 refuse les phases sur un maillage importé "
                   "(FdemSolver.cpp:1618) : le GBM passe par le maillage Voronoï.")
        if not m.fichier_msh:
            erreur("gmsh_fichier", "aucun fichier .msh : générer le maillage d'abord.")
        if (jg.alpha_ten, jg.alpha_coh, jg.alpha_gf, jg.alpha_e, jg.alpha_fric) != (1, 1, 1, 1, 1) or jg.paires:
            alerte("gmsh_gb_inerte", "sur un maillage importé tous les éléments sont dans un seul grain : "
                   "l'atténuation gbAlpha et les paires gb.* sont sans effet, sans avertissement du solveur.")
        if m.taille_element > 0:
            info("elements", "environ %d triangles (Gmsh, h = %g mm)."
                 % (nombre_elements(aire(essai), 1.0, m.taille_element), m.taille_element * 1e3))
    elif m.taille_element > 0 and m.taille_grain > 0:
        if m.taille_element >= m.taille_grain:
            erreur("taille_element", "l'élément doit être plus petit que le grain.")
        n = elements_par_grain(m.taille_grain, m.taille_element)
        if not ELEMENTS_PAR_GRAIN[0] <= n <= ELEMENTS_PAR_GRAIN[1]:
            alerte("elements_par_grain", "environ %.0f éléments par grain, hors de la plage 35-90 "
                   "(règle du 2026-09-07)." % n)
        info("elements", "environ %d triangles, %.0f par grain."
             % (nombre_elements(aire(essai), 1.0, m.taille_element), n))
    if m.type == "voronoi" and not m.aleatoire:
        alerte("maillage_regle", "maillage Delaunay sans points intérieurs aléatoires (grainMeshRandom absent) : "
               "la règle du 2026-09-07 l'impose en GBM. Admis seulement pour rejouer un deck antérieur.")
    if essai.schema.penalite_joint is not None and essai.schema.insertion == "adaptive":
        alerte("penalite_inerte", "jointPenaltyFactor est inerte en insertion adaptative (le solveur l'annonce) : "
               "la pénalité effective est insertionPenaltyFactor.")
    if m.dispersion_tailles is not None:
        if not 0.0 <= m.dispersion_tailles <= 1.5:
            erreur("dispersion", "grainSizeSpread doit être dans [0 ; 1,5].")
        if m.semis != "random":
            erreur("dispersion_semis", "la dispersion des tailles exige un semis aléatoire.")

    # ---- phases
    if ph:
        s = sum(p.fraction for p in ph)
        if abs(s - 1.0) > 1e-3:
            erreur("fractions", "les fractions surfaciques somment à %.4f au lieu de 1." % s)
        if any(p.fraction <= 0 for p in ph):
            erreur("fractions", "toute phase doit avoir une fraction positive.")
        if len(set(noms)) != len(noms):
            erreur("phases_doublon", "deux phases portent le même nom.")
    if any(p.taille_grain for p in ph) and len(ph) < 2:
        erreur("taille_phase", "une taille de grain par phase n'a de sens qu'avec au moins deux phases.")
    for p in jg.paires:
        if p.a not in noms or p.b not in noms:
            erreur("paire_inconnue", "la paire %s-%s cite une phase absente." % (p.a, p.b))

    # ---- discontinuités
    if d.actives():
        mu = essai.mu_residuel()
        if mu < 0:
            erreur("mu_residuel", "jointResidualMu doit être positif ou nul (FdemSolver.cpp:3356).")
        elif mu > math.tan(math.radians(mat.frictionDeg)) + 1e-4:
            alerte("mu_residuel", "frottement résiduel %.3f supérieur au frottement de pic %.3f : "
                   "le frottement augmenterait en s'endommageant (défaut du run F6, 2026-09-15)."
                   % (mu, math.tan(math.radians(mat.frictionDeg))))
        if d.fraction_diffuse and essai.schema.insertion == "adaptive":
            info("prerompus_effectifs", "en insertion adaptative une pré-fissure isolée reste inerte : "
                 "5 %% de consigne donnaient 2,2 %% de joints glissants (F7). Lire le log à l'initialisation.")
    if not 0.0 <= d.fraction_diffuse < 1.0:
        erreur("fraction_diffuse", "jointPrebrokenFrac doit être dans [0 ; 1[.")
    if d.plans:
        if d.plans.espacement <= 0 or d.plans.facteur <= 0:
            erreur("plans", "espacement et facteur des plans affaiblis doivent être positifs.")
        if not 0.0 <= d.plans.fraction_rompue <= 1.0:
            erreur("plans", "la fraction rompue d'une famille est dans [0 ; 1].")
    for s in d.segments:
        if not all(-1e-9 <= v <= L + 1e-9 for v, L in zip(s, (e.W, e.H, e.W, e.H))):
            erreur("segment_hors", "un segment pré-rompu sort de l'éprouvette : %s." % (s,))

    # ---- chargement
    if c.type_essai == "bresilien":
        R = e.W / 2
        if abs(e.W - e.H) > 1e-12:
            erreur("disque", "le brésilien se fait sur un disque : largeur et hauteur égales (le diamètre).")
        if m.type == "gmsh":
            erreur("bresilien_gmsh", "le disque du brésilien est construit par le solveur (discMesh = native) : "
                   "pas de maillage Gmsh importé.")
        if not 0.0 <= c.aplatissement_deg < 90.0:
            erreur("aplatissement", "l'angle total du méplat doit être dans [0 ; 90[ degrés.")
        if not 0.0 < c.demi_largeur_plateau <= R:
            erreur("plateau", "la demi-largeur des plateaux doit être positive et au plus égale au rayon.")
        lo, hi = c.jauge_elastique
        if not 0.0 < lo < hi < 1.0:
            erreur("jauge_elastique", "la bande de la jauge élastique doit vérifier 0 < bas < haut < 1 (× ft).")
    if c.type_essai not in ("triaxial", "traction", "bresilien"):
        erreur("type_essai", "type d'essai inconnu : %s" % c.type_essai)
    if c.vitesse <= 0:
        erreur("vitesse", "la vitesse de chargement doit être positive.")
    if c.type_essai == "triaxial":
        if c.sigma3_MPa < 0:
            erreur("sigma3", "le confinement est une compression : sigma3 >= 0.")
        if c.sigma3_MPa > 0 and c.delai_axial < c.rampe_confinement:
            alerte("consolidation", "la charge axiale démarre avant la fin de la rampe de confinement.")
        if c.sigma3_MPa > 0 and c.temps_jauge_confinement < c.rampe_confinement:
            alerte("jauge_confinement", "la jauge de confinement est figée avant la fin de la rampe : "
                   "confAchieved sous-estimerait sigma3.")
    elif c.sigma3_MPa:
        info("traction_sigma3", "en traction directe le confinement est ignoré.")

    # ---- calcul
    if essai.sorties.T <= 0 or essai.sorties.frames < 1:
        erreur("sorties", "T doit être positif et frames au moins 1.")
    return A


def erreurs(avis):
    return [a for a in avis if a.niveau == "erreur"]
