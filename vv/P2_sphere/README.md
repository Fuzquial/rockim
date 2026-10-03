# P2.2 Sphère élastique sur un plan : contact de Hertz

Banc de la campagne V&V de rockim (validation physique du contact général en
`fdem3d`). Script unique : `p2_sphere.py`.

## Objectif

Vérifier que le contact général par potentiel de Munjiza (`contact = potential`)
reproduit l'impact élastique de Hertz d'une sphère sur un massif : force
maximale, enfoncement maximal, durée de contact, restitution et bilan
d'énergie, et mesurer à quel point il s'en approche quand on raffine le
maillage au contact. C'est l'essai de la figure 8 de Fukuda et al. (RMRE 53,
2020), complété par la solution analytique de Hertz. La loi de potentiel n'est
pas exacte au sens de Hertz sur maillage grossier (surface facettée, pénétration
pénalisée) : le banc chiffre l'écart et sa convergence.

## Problème

- Sphère de rayon R = 20 mm, bloc de 120 x 120 x 60 mm (6R x 6R x 3R) encastré
  par sa base, même matériau : ρ = 2 650 kg/m³, E = 1 GPa, ν = 0,25.
- Vitesse d'impact v0 = 1 m/s selon -z, jeu initial 20 µm, sans frottement
  (`contactMu = 0`) ni pesanteur, sans amortissement (`scenario = loads`).
- Le choix E = 1 GPa fixe a_max/R = 0,114 : seul le rapport ρ v0²/E compte, et
  un module de roche (50 GPa) à la même vitesse donnerait a_max/R = 0,03, donc un
  maillage quatre fois plus fin au contact pour le même nombre de mailles sur le
  rayon de contact.
- Modèle quart : plans de symétrie x = 0 et y = 0 tenus en normale. Forces et
  masses du quart sont multipliées par 4 avant comparaison.
- Variante quasi statique (`qs_pf1`) : la sphère est menée à vitesse imposée
  0,2 m/s selon z (montée en cosinus sur 0,2 ms), donc rigide en z ; c'est
  l'indentation d'un demi-espace élastique par une sphère rigide.

## Solution exacte

Hertz (Johnson, Contact Mechanics, 1985, §4.2 et §11.4), bloc assimilé à un
demi-espace, masse du bloc infinie :

- module réduit E* = 1 / ((1 - ν²)/E + (1 - ν²)/E) = E / (2 (1 - ν²)) = 0,5333 GPa ;
- loi de contact F = K δ^(3/2), K = (4/3) E* √R ;
- équation du mouvement m δ'' = -K δ^(3/2), m = (4/3) π ρ R³ = 88,80 g ;
- enfoncement maximal δ_max = (15 m v0² / (16 E* √R))^(2/5) = 261,3 µm ;
- force maximale F_max = (4/3) E* √R δ_max^(3/2) = 424,8 N ;
- durée de contact t_c = 2,9432 δ_max / v0 = 769,1 µs ;
- rayon de contact maximal a_max = √(R δ_max) = 2,286 mm ;
- restitution e = 1. Dans un demi-espace, la part de l'énergie rayonnée en
  ondes élastiques vaut environ λ = 1,04 (v0/c0)^(3/5) avec c0 = √(E/ρ)
  (Hunter 1957, Johnson §11.4), soit λ = 2,2 % et e ≈ 0,989 ici. Dans le bloc
  fini, les ondes se réfléchissent (c0 t_c ≈ 0,47 m contre 0,06 m de bloc) :
  e est mesuré et documenté, avec un critère lâche.
- Indentation par sphère rigide : F = (4/3) E*_qs √R d^(3/2), E*_qs = E/(1 - ν²).

La courbe F(t) exacte est obtenue par intégration RK4 de l'équation du
mouvement (fonction `hertz_curve`).

Validité : t_c est grand devant R/c0 (33 µs), v0/c0 = 1,6e-3, a_max/R = 0,11
(petites déformations, pression maximale de Hertz environ 39 MPa, soit 4 %
de E). La souplesse du bloc fini (base encastrée à 3R) décale δ d'environ
0,2 a/H ≈ 1 %, ordre de grandeur estimé et non corrigé.

## Mise en œuvre

- Maillage Gmsh (API Python) écrit par le script : quart de sphère et quart de
  bloc disjoints, volumes nommés `sphere` et `block`, surfaces nommées `symx`,
  `symy`, `bottom`. Taille h_c = a_max / na dans une boule de rayon 1,5 a_max
  autour du point de contact, puis croissance linéaire (pente 0,3) jusqu'à 4 mm.
  Niveaux na = 2, 4, 8, 16 (h_c = 1,14 ; 0,57 ; 0,29 ; 0,14 mm ; 16 900,
  19 500, 29 900 et 90 600 tétraèdres).
- Clés rockim (toutes existantes) : `mode = fdem3d`, `scenario = loads`,
  `mesh = file`, `insertion = adaptive` avec ft et cohésion de 1e12 Pa (aucune
  insertion), `groupContinuum.block = true`, `groupContinuum.sphere = true`,
  `contact = potential`, `potPenaltyFactor`, `contactMu = 0`,
  `fix.bottom = all`, `fix.symx = x`, `fix.symy = y`,
  `groupVel.sphere = 0 0 -1` (impact) ou `velocity.sphere = free free -0.2`
  avec `amplitude.sphere = ramp 2e-4` (indentation), `trackGroup = sphere`
  (colonnes `grpZ`, `grpVz`, `grpFz`), `contactForcePairs = block:sphere`
  (colonnes `Fc_block_sphere_*`), `historyFlush = true`. La variante `pen`
  utilise `contact = penalty` et `gcRestitution = 1.0`.
- Variantes : `pot_pf1` (principale, p = E), `pot_pf10` (p = 10 E, part de la
  souplesse de contact), `pen` (pénalité nœud-face, pour mémoire, elle échoue
  le bloc glissant P2.1), `qs_pf1` (indentation quasi statique).

## Métriques

- δ(t) : rapprochement des points lointains = descente du centroïde de la
  sphère (`grpZ`) moins le jeu. Le centroïde n'est mû que par la force de
  contact, il est insensible aux vibrations propres.
- F(t) : 4 `Fc_block_sphere_z`. La force brute oscille autour de la courbe
  quasi statique (modes propres non amortis de la sphère et du bloc, facettes) :
  F_max est lu sur la force filtrée par moyenne glissante de largeur 0,05 t_c,
  la courbe exacte passant par le même filtre. Le pic brut est reporté.
- t_c : instant où δ repasse par zéro, compté depuis l'instant exact du premier
  contact (jeu / v0). La durée au seuil 1 % F_max de la force filtrée est
  reportée.
- e = vz final / v0 (`grpVz`), recoupé par l'impulsion ∫F dt = m (1 + e) v0.
- Bilan : résidu B4 du journal ; partition finale de l'énergie initiale entre
  translation de la sphère (e²), vibration de la sphère, énergie cinétique du
  bloc et énergie élastique stockée ; fermeture = somme / KE0 - 1.
- Indentation : écart L2 relatif de F(d) sur d dans [0,2 ; 1] δ_max.

## Critères d'acceptation

Fixés dans le script avant le premier calcul, appliqués au maillage le plus
fin de chaque variante au potentiel :

| grandeur | critère |
|---|---|
| F_max filtrée | écart ≤ 5 % |
| t_c | écart ≤ 5 % |
| δ_max | écart ≤ 5 % |
| e | ≥ 0,90 (documenté) |
| résidu B4 | ≤ 1e-6 % de l'échelle, chaque run |
| fermeture du bilan | ≤ 1 % |
| écart F_max | décroissant avec le raffinement |
| indentation F(d) | écart L2 ≤ 5 % |

Les seuils n'ont pas changé depuis l'essai court du 2026-10-04 (na = 2,
T = 0,6 t_c). Ce qu'il a montré a conduit à préciser deux définitions : F_max
sur la force filtrée et t_c sur le retour à zéro de δ.

## Essai court (na = 2, T = 0,6 t_c, 1 fil)

16 934 tétraèdres, dt = 4,30e-8 s, 11 205 pas, 97 s. Force filtrée maximale
443,8 N (+4,6 %), pic brut 577,8 N (+36 %), δ_max 252,0 µm (-3,6 %), écart L2
de F filtrée 6,3 %, résidu B4 7,9e-13 %. Sur ce maillage grossier le contact
est un peu trop raide (δ_max plus petit, F plus grande) : la raideur des
tétraèdres linéaires l'emporte sur la souplesse de la pénalité.

## Lancer

Depuis la racine rockim :

```
python3 vv/P2_sphere/p2_sphere.py prepare                 # maillages et decks
python3 vv/run_queue.py P2_sphere --slots 4 --threads 2   # la file commune
python3 vv/P2_sphere/p2_sphere.py analyse                 # resultats.json, figures
python3 vv/P2_sphere/p2_sphere.py all --threads 2         # ou tout en série
python3 vv/P2_sphere/p2_sphere.py smoke --tfrac 0.6       # essai court, out/_smoke
```

Sorties : `out/<variante>_na<na>/` (deck, journal, history.csv),
`meshes/sphere_plan_na<na>.msh`, `resultats.json`, `fig_impact`,
`fig_convergence`, `fig_indentation` (PDF et PNG). L'essai court vit dans
`out/_smoke/` et n'est jamais compté comme run de campagne.

Coût estimé (1 fil, mesuré à na = 2 puis extrapolé en tétraèdres x pas) :
na = 2 environ 4 min, na = 4 environ 9 min, na = 8 environ 26 min, na = 16
environ 2 h 30 par run d'impact ; l'indentation coûte 1,4 fois l'impact au
même maillage. Total des 9 runs : environ 15 000 s à 1 fil, soit 9 000 à
10 000 s cumulées à 2 fils ; en 4 slots de 2 fils, la file dure environ 2 h,
dominée par `pot_pf1_na16` (le pas de `pot_pf10` n'est pas encore mesuré :
s'il dépend de la pénalité normale, ses deux runs peuvent coûter davantage).
