# Guo 2014, chapitre 5 — couplage fluide-solide (Fluidity-ICOM × Y3D) et fracture hydraulique

*Fiche du 2026-08-28. Source : thèse Guo (Imperial College), p. 176-212, lues
intégralement sur PDF fourni par Fernando.* **[V]** = lu dans ce PDF ;
**[D]** = analyse personnelle. À lire avec
`specs/004-couplage-hydro-mecanique/spec.md` (AbuAisha, 2D, Draft).

## 1. L'architecture **[V]**

Couplage FAIBLE à corps immergé : Fluidity-ICOM (fluide FEM implicite,
Navier-Stokes incompressible, maillage adaptatif) × Y3D (solide FDEM
explicite). Le solide est vu du fluide comme une fraction volumique sur le
maillage fluide ; les échanges passent par un **supermesh** (Farrell) qui
garantit la conservation discrète et la 3e loi de Newton. Deux pas de temps :
Δt_fluide = 2,5e-6 s, Δt_solide = 5e-9 s → **500 pas solides par pas fluide**
(time staggering assumé). Travail propre de Guo : rendre le modèle de rupture
compatible — conversion du maillage DISCONTINU (tets + joints, aucun nœud
partagé) en maillage CONTINU (volume + surface) pour le fluide, par
regroupement des nœuds superposés à l'état initial (éq. 5.1-5.5 : moyennes
des coordonnées/vitesses du groupe).

## 2. Le suivi des fissures mouillées — la partie transférable **[V]**

Hypothèse : le solide non fissuré est IMPERMÉABLE (pas de leak-off). Quand un
joint casse, il ne devient frontière fluide-solide que selon **ES_i = nombre
de ses arêtes déjà sur la frontière solide (0-3)** :
- ES=0 (joint interne) : le fluide est « aveugle » à la fissure — cohérent
  avec l'imperméabilité ;
- ES=3 : détacher le tet, 2 triangles neufs en surface ;
- ES=2 : ouvrir un coin (1 nœud scindé) ;
- ES=1 : cas piège — traiter le joint SEUL exagérerait artificiellement
  l'aire de fissure ; l'algorithme attend un second joint cassé adjacent
  (ES=1 partageant une arête), d'où un **RETARD** du front fluide sur le
  front solide (visible dans son exemple, §5.6). Correctif proposé
  (fig. 5.14) : scinder d'abord les deux tets porteurs (nœud F → F1/F2)
  pour conserver EXACTEMENT l'aire de la fissure.

**[D]** C'est précisément la réponse structurée au manque identifié par la
spec 004 (« faces born from cracking receive nothing ») : quel que soit le
modèle de fluide retenu, il faut un REGISTRE des faces mouillées mis à jour à
chaque rupture de joint, avec une règle de connexité au réservoir — et le cas
ES=1 montre que la règle naïve retarde ou exagère. En 2D (spec 004) le
problème est plus simple : une fissure = arêtes, la connexité est un simple
parcours de graphe depuis le puits.

## 3. L'exemple numérique (§5.5) — un gabarit de validation **[V]**

100×49×50 mm, pré-fissure en coin 24,5 mm (ouverture 2,45 → élancement 1/10),
contraintes in situ 10/15/20 MPa sur les trois axes, pression fluide
équilibrée à 15 MPa puis rampe 2e4 MPa/s ; conçu pour une traction initiale
de 5 MPa = ft en tête. Maillage solide non structuré 2→5 mm. Résultats : la
fissure se propage quasi horizontalement (direction correcte vis-à-vis du
déviateur), charge critique correcte, concentration de contrainte en tête
capturée, le fluide entre dans la fissure neuve. Fluide à HAUTE viscosité
requis pour la stabilité (aveu explicite).

**[D]** Gabarit réutilisable pour valider le futur module 004 : pré-fissure +
in situ anisotrope + rampe de pression → deux critères fermés (direction de
propagation ⊥ à la contrainte la moins compressive ; charge critique) sans
avoir besoin de reproduire l'hydrodynamique fine.

## 4. Limites avouées (§5.6) et leçon d'architecture **[V]/[D]**

Pas de leak-off poro-élastique (même hypothèse qu'AbuAisha) ; fluide
incompressible → instabilités (d'où la haute viscosité) ; extension
compressible visée pour le TIR (blasting). Leçon pour rockim : la famille
« CFD immergé + supermesh » coûte un code fluide complet et un couplage
conservatif — la spec 004 a raison de partir d'AbuAisha (pression uniforme
inviscide dans le réseau connecté, chantier ~2× plus court), Lisjak (canaux,
loi cubique) en extension, Guo/Fluidity en horizon lointain si un jour le
JET (water-jet, cf. IJRMMS 2022) ou le nettoyage de fond de trou deviennent
des objectifs de simulation directs.

## 5. Découverte annexe : le chapitre 6 est la validation dynamique

La ToC du ch. 6 (p. 212, en fin de ce PDF) : « Application of the 3D fracture
model to the simulation of breakages of concrete armour units » — essais de
CHUTE (§6.3.1) et de PENDULE (§6.3.2) sur blocs Dolosse, intacts ET
pré-fissurés, COMPARÉS À DES EXPÉRIENCES PHYSIQUES (§6.3.3.1), puis
multi-corps Core-Loc sous gravité (§6.4, diversité des modes de ruine).
C'est le chapitre à demander en priorité : des benchmarks d'impact avec
cibles expérimentales, directement transposables en repères rockim.
