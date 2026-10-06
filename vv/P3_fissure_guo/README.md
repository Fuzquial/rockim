# P3.1 Fissure pressurisée, sensibilité au maillage (Guo 2014, §2.4)

Banc de la campagne V&V (`docs/VV_campagne.md`, lignes P3.1 et B1). Script : `p3_fissure_guo.py`.
Valeurs de référence de Guo : `guo_reference.json`. Résultats : `resultats.json` et figures
`fig_charge_rupture.pdf`, `fig_temps_propagation.pdf` (écrits par l'action `analyse`).

## Objectif

Reproduire en fdem3d l'étude de sensibilité au maillage de la thèse de Guo (Imperial College, 2014,
code Solidity 3D, §2.4, p. 79 à 100) et la confronter à deux références physiques indépendantes du
code : la mécanique linéaire de la rupture (LEFM) et la solution d'une fissure cohésive en plaque
infinie. Les grandeurs comparées sont la pression d'amorçage (« fracture load », fig. 2.29) et le
temps de propagation jusqu'au bord (fig. 2.30), pour cinq tailles de maille et deux énergies de
fissuration.

## Le problème de Guo

| Grandeur | Valeur | Source |
|---|---|---|
| Domaine | 120 × 120 mm dans le plan xy, épaisseur 20 mm selon z | texte p. 80, fig. 2.17 |
| Fissure | centrale, traversante, 40 mm de long (2a = 40 mm, a = 20 mm), toute l'épaisseur, plan normal à y | fig. 2.17 |
| Symétrie | demi-modèle droit, rouleau (u_x = 0) sur le plan yz central | texte p. 80 |
| Chargement | pression P = 1,0·10¹⁰ t (Pa) sur les deux lèvres, sans pesanteur | éq. 2.61 |
| Matériau | ρ = 2 340 kg/m³, E = 26 GPa, ν = 0,2, ft = 3 MPa, c = 15 MPa, φ = 30° | table 2.3 |
| Énergie de fissuration | Gf = 10 et 50 N/m | p. 81 |
| Frottement des lèvres | 0,6 | p. 81 |
| Maillages | tétraèdres structurés, h = 20 ; 10 ; 5 ; 2,5 ; 1,25 mm, soit 432 ; 3 456 ; 27 648 ; 221 184 ; 1 769 472 éléments | table 2.4 |
| Loi de joint | δp = 2fh/p0, E ≤ p0 ≤ 10E (valeur de §2.4 non publiée), δc = 3Gf/f | éq. 2.25 à 2.30 |

Le nombre d'éléments vaut 24 fois le nombre de cubes de côté h du demi-domaine (3 × 6 × 1 cubes
pour h = 20 mm) : le découpage est celui à 24 tétraèdres par cube (centre du cube et centres des
faces), reproduit ici.

### Valeurs extraites

Les dix valeurs de chaque figure ont été numérisées sur les images incluses dans le PDF de la
section (7 015 × 4 960 px, centroïdes des marqueurs, axes calés sur les graduations). Les quatre
valeurs citées dans le texte (p. 97) sont retrouvées à 0,01 MPa près, ce qui fixe l'incertitude de
lecture à ±0,02 MPa sur la charge et ±0,5 µs sur le temps.

| h (mm) | h/20 mm | P, Gf = 10 (MPa) | P, Gf = 50 (MPa) | t_prop, Gf = 10 (µs) | t_prop, Gf = 50 (µs) |
|---|---|---|---|---|---|
| 20 | 1 | 2,63 (texte) | 4,90 (texte) | 29,4 | 8,6 |
| 10 | 0,5 | 2,61 | 4,61 | 34,9 | 15,0 |
| 5 | 0,25 | 2,47 | 4,24 | 41,3 | 24,5 |
| 2,5 | 0,125 | 2,29 | 3,73 | 47,9 | 39,8 |
| 1,25 | 0,0625 | 2,10 (texte) | 3,09 (texte) | 52,3 | 63,0 |

Guo note lui-même l'absence de convergence en h (p. 99) : la charge baisse encore de 8 à 17 %
entre 2,5 et 1,25 mm. Son temps CPU par pas va de 5,0·10⁻⁴ s (h = 20 mm) à 2,62 s (h = 1,25 mm),
en série (table 2.5).

## Références physiques

Notations : a = 20 mm demi-longueur de fissure, b = W = 60 mm demi-largeur du domaine,
E' = E/(1 − ν²) = 27,08 GPa (déformation plane ; l'épaisseur de 20 mm est égale à a, l'état réel
est intermédiaire et E' = E donnerait des pressions 2 % plus basses).

### Mécanique linéaire de la rupture

Une pression P sur les lèvres d'une fissure de longueur 2a donne, par superposition, le même
facteur d'intensité qu'une contrainte lointaine P :

K_I = P √(πa) F(a/b), K_Ic = √(E' Gf), donc P_c = √(E' Gf) / (F √(πa)).

F est le facteur de largeur finie de Tada pour une fissure centrale dans une bande de largeur 2b :
F = (1 − 0,025 α² + 0,06 α⁴) √(sec(πα/2)), α = a/b = 1/3, F = 1,0724. La hauteur finie de la
plaque (120 mm) n'est pas corrigée.

| Gf (N/m) | K_Ic (MPa√m) | P_c plaque infinie (MPa) | P_c largeur finie (MPa) |
|---|---|---|---|
| 10 | 0,520 | 2,08 | 1,94 |
| 50 | 1,164 | 4,64 | 4,33 |

### Longueur caractéristique et fissure cohésive

La LEFM suppose une zone d'élaboration petite devant a. La longueur de Hillerborg vaut
l_ch = E Gf / ft² = 28,9 mm pour Gf = 10 et 144 mm pour Gf = 50, soit a/l_ch = 0,69 et 0,14 : la
zone cohésive n'est pas petite, et la charge d'amorçage d'un modèle cohésif doit dépasser la
valeur LEFM (effet d'échelle de résistance).

Le script calcule la solution exacte d'une fissure cohésive pressurisée en plaque infinie et
déformation plane (fonction `cohesive_pc`). La fissure fictive a une demi-longueur c ; les lèvres
réelles (|x| < a) portent P, la zone a < |x| < c porte la traction cohésive σ(δ). L'ouverture est
celle de Sneddon par fonction de Green :

δ(x) = (4/(πE')) ∫₀ᶜ p(s) ln| (√(c² − s²) + √(c² − x²)) / (√(c² − s²) − √(c² − x²)) | ds,

et la fermeture régulière impose K = 0 en c : P arcsin(a/c) = ∫ₐᶜ σ(s) / √(c² − s²) ds. La zone
cohésive est discrétisée en 40 panneaux constants (collocation au milieu) ; l'amorçage est
l'instant où l'ouverture en x = a atteint δc. Écart entre 40 et 80 panneaux : < 0,02 %. Pour
Dugdale (σ = ft), la condition K = 0 donne c en forme fermée, a/c = sin(π ft / (2(P + ft))).

| Loi (plaque infinie) | Gf = 10 : P (MPa), c (mm) | Gf = 50 : P (MPa), c (mm) |
|---|---|---|
| Dugdale, δc = Gf/ft | 2,72 ; 27,3 | 7,55 ; 46,3 |
| linéaire, δc = 2Gf/ft (loi par défaut de rockim) | 3,33 ; 33,4 | 10,15 ; 68,4 |
| z-curve de Munjiza, δc = 3Gf/ft (convention de Guo) | 4,04 ; 39,8 | 12,83 ; 89,5 |

Pour Gf = 50 la zone cohésive de la plaque infinie (c = 68 mm) dépasse le domaine (W = 60 mm) :
le domaine fini gouverne. L'équilibre du demi-domaine supérieur borne alors la pression : la
résultante P·a des lèvres est reprise par le ligament de longueur W − a, d'où
P ≤ ft (W − a)/a = 6,0 MPa (ligament entièrement à ft).

Lecture attendue : pour Gf = 10, la pression d'amorçage d'un modèle cohésif convergé tombe entre la
LEFM en largeur finie (1,94 MPa) et la solution cohésive linéaire en plaque infinie (3,33 MPa), la
largeur finie abaissant la seconde. Pour Gf = 50, entre 4,33 et 6,0 MPa. Les valeurs les plus fines
de Guo (2,10 et 3,09 MPa) sont au niveau de la LEFM pour Gf = 10 et en dessous pour Gf = 50 ; ce
désaccord est un des enjeux du banc. La vitesse de chargement (une durée de chargement d'environ
dix allers-retours d'onde) ajoute un effet dynamique non pris en compte par ces références.

## Mise en œuvre dans rockim

### Fissure préexistante

Les clés `pressure.<g>` (§5.21 de la documentation) exigent des faces extérieures, et refusent un
triangle intérieur à un corps ou posé sur l'interface ambiguë de deux corps non liés. La fissure est
donc modélisée par des nœuds dupliqués, comme le plugin Crack de Gmsh : sur le plan y = 60 mm, les
sommets de la partie fissurée (x < a, centres de faces compris) reçoivent un double que seuls les
tétraèdres supérieurs référencent. `buildMeshFile` ne soude pas les nœuds coïncidents : les
triangles des deux lèvres ont des sommets distincts, deviennent des faces extérieures distinctes du
même corps, sans joint entre elles, et chacun appartient à un seul tétraèdre (aucune ambiguïté). Le
front x = a n'est pas dupliqué. Chaque lèvre est un groupe physique de surface (`lip_lo`,
`lip_hi`) ; la pression suiveuse s'applique selon la normale sortante du tétraèdre porteur, donc
ouvre la fissure. Le même maillage marche en fem3d (vérifié). Le contact général traite les lèvres
si elles se referment. Aucune modification du code n'a été nécessaire.

### Maillage

Générateur propre (`ensure_mesh`), MSH 2.2 ASCII, demi-domaine x ∈ [0 ; 60] mm, 24 tétraèdres
congruents par cube, nombres d'éléments identiques à la table 2.4. Groupes : `sym` (face x = 0),
`lip_lo`, `lip_hi`, `mouth_lo`, `mouth_hi` (sommets de la bouche x = 0, y = 60 mm, à mi-épaisseur
si le maillage en a un, sinon les deux faces z = 0 et z = 20 mm, symétriques), volume `roche`.

### Decks

Clés communes : `mode`, `scenario = loads`, `mesh = file`, `meshFile`, `T`, `frames = 1`,
`historyFlush = true`, `rho`, `E`, `nu`, `ft`, `cohesion`, `frictionDeg`, `Gf`,
`gfShearFactor = 1` (un seul Gf comme chez Guo), `contactMu = 0.6`, `fix.sym = x`,
`pressure.lip_lo` et `pressure.lip_hi` égales à RATE·T avec `amplitude.<g> = 0 0 T 1` (donc
P = 10¹⁰ t exactement), `force.mouth_lo = 0 0 0` et `force.mouth_hi = 0 0 0` (colonnes
d'ouverture `U_mouth_*` dans `history.csv`). Durées : T = 500 µs pour Gf = 10 (P final 5 MPa),
800 µs pour Gf = 50 (8 MPa, au-delà de la borne de ligament).

| Variante | Clés propres | Rôle |
|---|---|---|
| `intr_pf100` (principale) | `insertion = intrinsic`, `jointPenaltyFactor = 100`, `jointXi = 0.01` | loi par défaut (adoucissement linéaire). pf = 100 : biais de célérité prédit −0,6 % par la loi de P1.1, marge sur le seuil pf ≥ 62, parce que α = 1,24 y a été mesuré avec hmin commandé par des tétraèdres aplatis alors qu'ici tous les tétraèdres sont congruents (hmin = 0,26 h) |
| `guo_pf5` | intrinsèque, `jointPenaltyFactor = 5`, `jointPenaltyLength = edge`, `jointElastic = parabolic`, `jointSoftening = munjiza`, `jointDeltaC = guo`, `jointDeath = damage`, `jointXi = 0` | la loi de la thèse (éq. 2.25 à 2.31) ; p0 = 10E, borne haute de l'éq. 2.28, soit pf = p0/(2E) = 5 |
| `adapt` | `insertion = adaptive`, `jointXi = 0.01` | continuum exact avant insertion |
| `elas_fem3d` | `mode = fem3d`, `law = elastic`, T = 100 µs | référence élastique du même maillage pour la souplesse |

`jointXi = 0.01` suit la règle maison du quasi-statique. h = 1,25 mm n'est lancé qu'avec `--fine`.

## Métriques

| Métrique | Définition |
|---|---|
| P_onset | 10¹⁰ × t_onset, t_onset = plus petit `tBreak` des joints du VTU final (recoupé par le premier nBroken ≥ 1 de `history.csv`, à un intervalle d'historique près) ; coordonnées du premier joint rompu dans `resultats.json` |
| t_reach | plus petit `tBreak` des joints de la bande |y − 60 mm| < h/2 qui touchent le bord x = 60 mm |
| t_prop | t_reach − t_onset (temps de propagation de la fig. 2.30) |
| t_split, fraction de ligament | dernier `tBreak` des joints du plan y = 60 mm au-delà du front, et fraction rompue de ces joints |
| joints hors bande | joints rompus hors de la bande, pour repérer une fissure qui dévie ou de l'endommagement diffus |
| souplesse | intégrale en temps de l'ouverture de la bouche sur [5 ; 60] µs (P ≤ 0,6 MPa, avant tout endommagement de la pointe), rapportée au fem3d du même maillage ; le rapport élimine les oscillations dynamiques, communes aux deux calculs |
| bilan | résidu B4 imprimé par le solveur, en % de l'échelle |

Les positions des joints sont lues dans la trame 0 (configuration de référence) et les `tBreak` dans
la trame finale : après la séparation, les deux moitiés s'écartent de plus d'un millimètre sous la
pression qui continue de croître.

## Critères d'acceptation

Fixés dans le script le 2026-10-04, avant le premier calcul de la campagne (les essais de fumée ne
portent que sur h = 20 et 10 mm).

| Critère | Énoncé |
|---|---|
| C1 convergence | variante principale, |P(h_min)/P(h suivant) − 1| ≤ 5 % pour chaque Gf, sur les deux maillages les plus fins |
| C2 Guo | variante `guo_pf5`, |P/P_Guo − 1| ≤ 15 % à h ≤ 10 mm ; écart seulement rapporté pour les autres variantes |
| C3 physique | variante principale au maillage le plus fin : Gf = 10, P entre 1,94 (LEFM, largeur finie) et 3,33 MPa (cohésif linéaire, plaque infinie) ; Gf = 50, P entre 4,33 et 6,0 MPa (borne de ligament) |
| C4 propagation | la fissure atteint le bord x = 60 mm avant T dans tous les runs fdem3d |
| C5 souplesse | variante principale à 2 % du fem3d sur le même maillage |
| C6 bilan | résidu B4 ≤ 0,1 % de l'échelle |

## Essais de fumée (2026-10-04, `OMP_NUM_THREADS = 1`, sorties dans `out/smoke/`)

| Run | Tétraèdres | dt (s) | P_onset (MPa) | t_prop (µs) | Guo P / t_prop | Temps (s) |
|---|---|---|---|---|---|---|
| intr_pf100, h = 20, Gf = 10 | 432 | 1,83·10⁻⁸ | 1,99 | 32,5 | 2,63 / 29,4 | 6,7 |
| intr_pf100, h = 10, Gf = 10 | 3 456 | 9,13·10⁻⁹ | 1,91 | 54,1 | 2,61 / 34,9 | 106 |
| intr_pf100, h = 20, Gf = 50 | 432 | 1,83·10⁻⁸ | 4,16 | 18,5 | 4,90 / 8,6 | 10,7 |
| guo_pf5, h = 20, Gf = 10 | 432 | 8,43·10⁻⁸ | 2,22 | 30,3 | 2,63 / 29,4 | 1,6 |
| guo_pf5, h = 20, Gf = 50 | 432 | 8,43·10⁻⁸ | 4,50 | 15,3 | 4,90 / 8,6 | 2,4 |
| adapt, h = 20, Gf = 10 et 50 | 432 | 7,80·10⁻⁸ | aucune insertion jusqu'à 5 et 8 MPa | | | 1,4 et 2,2 |
| adapt, h = 10, Gf = 10 | 3 456 | 3,90·10⁻⁸ | 4,48 | 44,6 | 2,61 / 34,9 | 23 |
| elas_fem3d, h = 20 et 10 | | 4,46·10⁻⁷ | | | | < 1 |

Le premier joint rompu est sur le plan de la fissure, juste devant le front (x = 21,7 mm à
h = 10 mm), et tout le ligament rompt. Bilans B4 entre 10⁻¹⁴ et 10⁻¹³ %. Souplesse de la bouche
rapportée au fem3d : +0,8 % pour `intr_pf100` (h = 20 et 10 mm), +32 % pour `guo_pf5`, 0 pour
`adapt`. Une première mesure par pente sur P ≤ 1 MPa donnait +6 % à h = 10 mm : l'endommagement
sous-critique de la pointe commence vers 0,6 à 1 MPa en intrinsèque, d'où la fenêtre retenue.

Deux constats déjà visibles. L'intrinsèque amorce près de la LEFM (1,9 à 2,0 MPa), bien sous la
solution cohésive (3,33 MPa), parce que les joints de la pointe lisent la séparation nodale,
singulière, et non une contrainte moyenne. L'adaptatif, au contraire, n'insère rien à h = 20 mm et
amorce à 4,5 MPa à h = 10 mm : avec des mailles grossières, l'effort passe surtout par les
tétraèdres qui ne touchent le plan qu'au front (σ1 jusqu'à 6,5 MPa à P = 8 MPa), alors que les
tétraèdres qui portent les facettes du ligament restent sous ft (2,6 MPa) ; le critère par
contrainte d'élément moyennée est aveugle à la concentration de pointe tant que h n'est pas petit
devant l_ch. `facetAverage = max` et `nodal` n'insèrent rien non plus à h = 20 mm.

## Coût estimé de la campagne

Coût mesuré de la variante principale : 5,6·10⁻⁷ s par tétraèdre et par pas sur un fil, à h = 20
comme à 10 mm. Le nombre de tétraèdres croît comme h⁻³ et le nombre de pas comme h⁻¹ : facteur 16
par division de h par deux. Estimation à 2 fils (gain supposé 1,7), hors coût du contact après la
séparation :

| h (mm) | intr_pf100, Gf = 10 / 50 | guo_pf5, Gf = 10 / 50 | adapt, Gf = 10 / 50 |
|---|---|---|---|
| 10 | 1 min / 2 min | < 1 min | < 1 min |
| 5 | 17 min / 27 min | 4 min / 7 min | 4 min / 6 min |
| 2,5 | 4,4 h / 7,1 h | 1,1 h / 1,7 h | 1,0 h / 1,6 h |
| 1,25 | 71 h / 114 h | 17 h / 28 h | 16 h / 26 h |

h = 1,25 mm n'est pas abordable pour la variante principale (environ 8 jours à 2 fils pour les deux
Gf) ; il reste possible pour `guo_pf5` et `adapt` si la queue est libre. Les VTU d'éléments à
h = 2,5 mm pèsent de l'ordre de 100 Mo par trame (deux trames par run).

## Reproduire

```bash
cd ~/rockim
python3 vv/P3_fissure_guo/p3_fissure_guo.py prepare                    # maillages et decks
python3 vv/run_queue.py P3_fissure_guo --slots 4 --threads 2            # campagne sans 1,25 mm
python3 vv/P3_fissure_guo/p3_fissure_guo.py analyse                    # resultats.json, figures
python3 vv/P3_fissure_guo/p3_fissure_guo.py run --fine --h 1.25 --var guo_pf5 --threads 2
python3 vv/P3_fissure_guo/p3_fissure_guo.py smoke --h 20 --gf 10 --threads 1   # fumée, out/smoke/
```

`--var`, `--h` et `--gf` restreignent la campagne ; `--T` raccourcit la durée en fumée.

## Réserves

- La pression de Guo est nominale ; sa loi de joint (z-curve, δc = 3Gf/f) dissipe 1,159 Gf. La
  variante `guo_pf5` la reproduit, à la quadrature près (points aux milieux d'arêtes chez Guo, aux
  nœuds dans rockim, §5.4 quater) et à la valeur de p0, non publiée.
- La définition de l'amorçage chez Guo n'est pas plus précise que « la pression à laquelle la
  fissure commence à se propager » : il peut s'agir du premier joint rompu ou d'une lecture sur les
  vues. Les instants des vues des fig. 2.19 à 2.28 sont dans `guo_reference.json` pour recouper.
- Le modèle n'est libre en y et en z que par l'équilibre des charges (rouleau seul, comme chez
  Guo) ; aucune dérive n'a été vue en fumée.
- Les références LEFM et cohésive sont planes et en plaque infinie (largeur finie corrigée
  seulement pour la LEFM) ; la hauteur finie et l'effet tridimensionnel du bord libre en z ne sont
  pas corrigés.
- `adapt` n'est pas comparable à h ≥ 10 mm (voir les constats ci-dessus).
