# Campagne de vérification et de validation de rockim

Document vivant de la campagne V&V, commencée le 2026-10-03. La partie I fixe le cadre et fait
l'inventaire de l'existant ; la partie II décrit chaque banc de vérification, sa solution de
référence, sa mise en œuvre, ses critères et ses résultats ; la partie III est réservée à la
validation expérimentale. Les scripts et les decks sont dans `vv/`, un dossier par banc.

| Banc | Objet | Statut | Section |
|---|---|---|---|
| V1 | onde de compression dans une barre | fait : fem3d passe ; fdem3d adaptatif passe sauf le bilan ; fdem3d intrinsèque biaisé de −3 % | II.1 |
| V2 | joint seul, modes I, II et mixte, énergie Gf | prévu | |
| V3 | énergie de fissuration en traction directe | prévu | |
| V4 | contact de Hertz sphère-plan | prévu | |
| V5 | choc de deux barres | prévu | |
| V6 | patch test | prévu | |
| V7 | convergence en h et en dt | couvert en partie par V1 | II.1 |
| V8 | objectivité de la rupture | prévu | |
| V9 | bloc sur plan incliné | prévu | |

Environnement des résultats de la partie II, sauf mention contraire : macOS 27 arm64, AppleClang
21, OpenMP de Homebrew (libomp 21.1.6), `-O3 -ffp-contract=off`, binaire `build_nofma/rockim` de la
branche `claude/rockim-version-check-4ugw9e`.

## Partie I : cadre et inventaire (2026-10-03)

### 1. Cadre

Le document sépare deux démonstrations, selon ASME V&V 10 (2019) et Oberkampf et Roy (2010).

La vérification établit que le code résout correctement les équations qu'il implémente. Sa
référence est une solution analytique, une propriété mathématique (conservation, invariance,
ordre de convergence) ou une implémentation indépendante de la même loi. Elle n'utilise aucune
donnée expérimentale.

La validation établit que ces équations et ces lois représentent le granite sous percussion. Sa
référence est un essai, et cet essai ne doit pas avoir servi à caler les paramètres.

Un test de non-régression, dont la référence est une sortie antérieure de rockim, ne relève
d'aucune des deux catégories. Il garantit que le code n'a pas changé, pas qu'il est juste.

### 2. Classement des 118 tests de `tools/verify_suite.py`

| Catégorie | Nombre | Nature de la référence |
|---|---|---|
| Analytique au point matériel ou sans maillage | 11 | formule exacte, VUMAT Fortran, conservation à deux corps |
| Analytique sur une structure | 8 | barre en traction, allongement σL/E et réaction exacts |
| Invariant au repos | 26 | aucune rupture et aucune injection d'énergie sous charge nulle |
| Invariant d'équivalence | 11 | deux réglages qui doivent donner le même résultat, ou un poste d'énergie qui doit fermer le bilan |
| Pic contre ft, référence figée | 24 | écart du pic à ft, mais la tolérance porte sur l'écart mesuré un jour donné (par exemple −1,39 % ± 0,01) |
| Régression pure | 38 | comptages de joints rompus, taux médians, pics, UCS en sortie de rockim |

Les 56 tests des quatre premières lignes ont une référence indépendante du code : ce sont des
tests de vérification. Les 24 tests de pic contre ft comparent bien à une
grandeur théorique, mais leur verdict porte sur une valeur figée de l'écart, si bien qu'un écart
de 1,4 % et un écart de 15 % seraient traités de la même façon s'ils avaient été enregistrés tels
quels. Les 38 tests de régression ne démontrent pas la justesse, à quelques contrôles près qu'ils
portent en plus (bilan d'énergie de `bench1_insert_impact`, travail de contact nul de
`shpb_elastic_potential`).

Les 11 échecs de la suite `fast` sur macOS arm64 (39 tests sur 50 passent) appartiennent tous
aux catégories pic contre ft et régression. Les écarts portent sur des grandeurs post-pic que
l'arithmétique flottante de la plateforme déplace ; aucun test de vérification n'échoue.

### 3. Comparaisons existantes hors de la suite

#### 3.1 Vérification analytique

- Kirsch (λ = 1 et 0,5) : écarts de 1,7 % et 2,1 % ; Lamé −2,0 % ; Sneddon −4 %
  (`tunnel_edz/`, `bench_abuaisha/VALIDATION_hydro.md`).
- Bornes de Voigt et Reuss sur le module d'un agrégat : 69,24 GPa entre 66,77 et 71,21 GPa
  (`gbm_uniaxial/UNIAXIAL.md`, `bench_phases/`).
- Critère de Camacho-Ortiz sur un cube en traction : pic de 9,72 MPa pour ft = 10 MPa
  (journal de thèse du 2026-09-11).
- Élasticité isotrope transverse en forme fermée : écarts de 0,00 et 0,13 % (`tests_f2/t13*`,
  `t14*`).
- Points d'entrée `selftest-triax`, `selftest-cdp` (79 contrôles), `selftest-fixed`,
  `selftest-dfhplus` et `thermobench`, hors de la suite.

#### 3.2 Comparaison avec un autre code

- Abaqus, loi saksala2011, même maillage, sphère à 8 m/s (31/07) : pic filtré +7,9 %, énergie
  +0,9 %, pénétration +1,1 %, restitution 0,111 contre 0,059. Les critères avaient été fixés
  avant le calcul.
- Abaqus, loi CDP, quart de bloc à 11 m/s (05-06/09) : le pic passe de −20 % à +2,8 % après
  l'ajout de quatre clés choisies pour se rapprocher d'Abaqus. La concordance finale est donc
  un ajustement, pas une vérification indépendante.
- VUMAT Fortran saksala2011 et kstdfh : écarts de 8e-14 et 4,7e-12 au point matériel.
- Code public de Solidity : transcription de la loi de joint et cycle fermé ; Solidity crée
  lui-même 1,9e-8 J par cycle.

#### 3.3 Validation expérimentale

- Saint-Anne, insert unique à 10,66 m/s, contre Yang 2025 (essais d'Aising) : contrainte dans
  le taillant +6 %, enfoncement maximal +17 % (1,29 contre 1,10 mm simulés par Yang, 1,05 mm
  mesuré), rayon de cratère −12 %, masse de fragments +87 % (4,67 g contre 1,2 g mesuré). Les
  paramètres viennent de la table 4 de Yang, qui les a calés sur ces mêmes essais : la
  comparaison vérifie la reproduction de Solidity et ne constitue pas une validation
  indépendante. Les valeurs de référence sont lues sur des figures, à ±10 %.
- Kuru Grey à 9 m/s contre Yang 2026 : comparaison partielle ; le run s'arrête à 183 µs, la
  pulvérisation reste à 0-2 éléments contre environ 360 chez Yang.
- Red Bohus (Dumoulin 2024) : UCS, brésilien et triaxial à 20-100 MPa servent tous à la
  calibration. Le jeu bayésien a été invalidé par la correction d'amortissement.
- Courbe force-pénétration de Kuru (k₁ = 162 kN/mm, pic de 113 kN à 0,83 mm) : numérisée dans
  `percussion_bohus/fig_fp_exp_simu.py`, comparée à Abaqus seulement, jamais à rockim.
- Aucun essai SHPB n'est comparé à rockim.

### 4. Lacunes

Les lacunes de la vérification sont les suivantes.

1. La propagation d'onde n'est vérifiée nulle part dans le solveur : célérité c = √(E/ρ),
   amplitude σ = ρcv, réflexion sur un bord libre et sur un bord encastré. Le contrôle T0b de
   `selftest-toolcontact` traite la barre de Saint-Venant en forme fermée, sans maillage.
2. Le contact n'a pas de référence analytique en force : Hertz sphère-plan, durée de choc de
   deux barres égale à 2L/c, glissement d'un bloc sur un plan incliné.
3. L'énergie dissipée par une fissure n'est jamais comparée à Gf multiplié par l'aire rompue,
   ni la réponse d'un joint seul à sa loi de traction-séparation en mode I, en mode II et en
   mode mixte.
4. Aucune étude ne mesure un ordre de convergence en maillage ou en pas de temps sur un cas
   élastique. `docs/MAILLAGE_serie_2026-09-13.md` documente une série de maillages en impact,
   sans référence exacte.
5. L'objectivité de la rupture vis-à-vis du maillage (pic, énergie dissipée, longueur
   fissurée en fonction de h) n'est pas quantifiée. C'est la propriété la plus discutée du
   FDEM à insertion intrinsèque.
6. Le patch test, en élasticité linéaire et en grandes transformations, est absent.
7. Les 24 tests de pic contre ft n'ont pas de critère absolu.

Les lacunes de la validation sont les suivantes.

1. Aucune comparaison n'est faite sur des données qui n'ont pas servi à caler les paramètres.
2. Aucun essai dynamique de laboratoire (SHPB, brésilien dynamique) n'est comparé.
3. Les références de Yang ne sont pas numérisées ; elles vivent dans des tableaux Markdown et
   dans des scripts.
4. Aucune incertitude expérimentale ni critère d'acceptation n'est fixé avant le calcul, sauf
   pour la comparaison Abaqus du 31/07.

### 5. Plan proposé

#### 5.1 Vérification

| Rang | Cas | Référence | Grandeur et critère | Solveurs |
|---|---|---|---|---|
| V1 | Barre 1D sous impulsion | c = √(E/ρ), σ = ρcv, réflexion | temps d'arrivée à 1 %, amplitude à 2 %, convergence en h | fem3d, fdem3d intrinsèque et adaptatif |
| V2 | Joint seul, modes I, II et mixte | loi de traction-séparation, énergie Gf | courbe à 1e-6, énergie dissipée égale à Gf·A à 0,1 % | fdem 2D et 3D |
| V3 | Énergie de fissuration en traction directe | Gf·A, A mesurée sur la fissure | écart inférieur à 5 %, pour trois tailles de maille | fdem3d |
| V4 | Hertz sphère-plan, quasi statique | F = (4/3) E* √R δ^(3/2) | force à 3 % sur la plage élastique, convergence en h | fem3d, fdem3d avec pénalité et potentiel |
| V5 | Choc de deux barres | durée 2L/c, vitesses de sortie | durée à 2 %, quantité de mouvement exacte | fdem3d |
| V6 | Patch test | champ uniforme exact | erreur au niveau de l'arrondi | fem3d, fdem3d |
| V7 | Convergence en h et en dt sur un cas élastique | V1 ou Kirsch | ordre observé contre ordre théorique | fem3d, fdem3d |
| V8 | Objectivité de la rupture | invariance attendue de Gf | pic et énergie en fonction de h et de la pénalité, tendance documentée | fdem3d |
| V9 | Bloc sur plan incliné | seuil tan φ, accélération g(sin θ − μ cos θ) | accélération à 1 % | fdem3d |

V1, V2 et V4 viennent en premier : ce sont les trois phénomènes qui gouvernent la percussion,
et chacun coûte quelques minutes de calcul. Les 24 tests de pic contre ft reçoivent en plus un
critère absolu, par exemple |écart| inférieur à 5 %, en conservant la référence figée pour la
non-régression.

#### 5.2 Validation

| Niveau | Essai | Rôle | Statut des données |
|---|---|---|---|
| 1 | UCS, brésilien, triaxial Red Bohus | calibration | disponibles (`calibration_redbohus/targets/`) |
| 2 | Triaxial à 75 et 100 MPa, non utilisé pour le calage | validation quasi statique | disponibles ; à retirer du jeu de calage |
| 3 | SHPB ou brésilien dynamique sur granite | validation de la loi de taux | à rechercher dans la littérature |
| 4 | Force-pénétration de Kuru, impact d'Aising | validation de l'impact d'un insert | numérisation à faire |
| 5 | Cratère et fragments de Saint-Anne | validation du système | lecture de figures à ±10 %, réserve à porter |

Chaque comparaison de validation fixe avant le calcul la grandeur comparée, l'incertitude de
l'essai, le critère d'acceptation et les paramètres gelés. Les paramètres calés au niveau 1 ne
sont plus modifiés aux niveaux 2 à 5.

#### 5.3 Livrables

- `tools/verify_suite.py` : une étiquette `kind` par test (`analytique`, `invariant`,
  `regression`, `validation`) et un rapport qui compte les tests par étiquette.
- `configs/vv/` : un deck par cas V1 à V9, avec un script de dépouillement qui trace la
  simulation contre la référence.
- `data/experimental/` : les références numérisées au format CSV, avec leur source et leur
  incertitude.
- Un chapitre de méthode pour la thèse, construit sur la matrice des sections 5.1 et 5.2.


## Partie II : bancs de vérification

### II.1 V1, onde de compression dans une barre

#### Objectif

Vérifier que rockim propage une onde élastique à la bonne célérité, avec la bonne impédance et
les bonnes réflexions, et mesurer l'ordre de convergence en maillage. La percussion est un
problème d'ondes : une erreur à ce niveau se retrouverait dans tous les bancs d'impact.

#### Problème

La barre mesure 8 × 8 × 400 mm, d'axe z, en matériau élastique (ρ = 2 650 kg/m³, E = 60 GPa,
ν = 0,25). Les quatre faces latérales sont sur rouleaux (`fix.xmin = x`, `fix.xmax = x`,
`fix.ymin = y`, `fix.ymax = y`). La déformation est alors uniaxiale et le problème continu est
exactement unidimensionnel, sans la dispersion géométrique de Pochhammer qu'aurait une barre à
faces libres. La célérité est celle des ondes de compression en déformation uniaxiale :

c = √(M/ρ), avec M = E(1 − ν)/((1 + ν)(1 − 2ν)) = 72 GPa, soit c = 5 212,47 m/s.

La face `bottom` (z = 0) reçoit une vitesse imposée V(t) = v₀ a(t), avec v₀ = 0,1 m/s et a(t) un
trapèze : montée linéaire de 0 à 1 en 4 µs, palier jusqu'à 16 µs, retour à 0 à 20 µs. La
contrainte vaut ρcv₀ = 1,38 MPa, dans le domaine élastique. La face `top` (z = L) est libre. La
durée simulée, T = 250 µs, couvre 1,6 aller-retour (2L/c = 153,5 µs).

#### Solution exacte

Soit D(t) le déplacement imposé, intégrale de V de 0 à t, nul pour t < 0. La solution de
d'Alembert avec réflexions s'écrit

u(z, t) = Σₙ (−1)ⁿ [ D(t − (2nL + z)/c) + D(t − (2(n+1)L − z)/c) ].

Le premier terme est l'onde descendante, le second l'onde réfléchie par le bout libre, qui revient
avec le même signe. Le bout chargé, à déplacement imposé, renvoie l'onde avec un signe opposé,
d'où le facteur (−1)ⁿ. La réaction de l'appui sur le solide vaut

RF(t) = −M A ∂u/∂z(0, t), avec ∂u/∂z(0, t) = Σₙ (−1)ⁿ [ −V(t − 2nL/c) + V(t − 2(n+1)L/c) ] / c.

Son plateau pendant l'impulsion vaut ρcv₀A = 88,40 N. Au retour de l'onde, à t = 2L/c, la
réaction double et change de signe.

#### Mise en œuvre

Le script `vv/V1_onde/v1_onde.py` enchaîne cinq étapes.

1. Il génère les maillages tétraédriques non structurés par
   `tools/make_unstructured_mesh.py box3dbc 0.008 0.008 0.4 h`, qui nomme les six faces et les
   huit coins. Les maillages vont dans `vv/V1_onde/meshes/`.
2. Il écrit un deck par variante et par taille de maille, en `scenario = loads` (§5.21 de la
   documentation). Les stations de mesure s1, s2 et s3, à z = 100, 200 et 300 mm, sont des
   groupes `point.sK` chargés par une force nulle (`force.sK = 0 0 0`). Cette charge nulle fait
   écrire le déplacement moyen du groupe, `U_sK_z`, dans `history.csv`. Le procédé marche en
   fem3d comme en fdem3d, qui n'a pas de sondes `probes`. Le bout libre est lu de la même façon
   par `force.top = 0 0 0`, et la réaction de l'appui par `RF_bottom_z`.
3. Il lance rockim dans `vv/V1_onde/out/<variante>_h<h>/`, avec le journal `rockim.log`.
4. Il relit la position réelle de chaque station. Sur un maillage non structuré, `point.` retient
   le sommet le plus proche, à 1,8 mm du point demandé au plus sur le maillage de 2 mm ; le
   journal imprime ses coordonnées (`point.s1 : sommet 4274 a 0.00163 m du point demandé
   (… 0.100618, repère solveur)`). La solution exacte est évaluée à cette position, de sorte que
   l'erreur mesurée est celle du solveur et non celle du placement.
5. Il compare, calcule les métriques, écrit `vv/V1_onde/resultats.json` et trace les figures.

La solution exacte est évaluée sur les instants de `history.csv` (environ 2 000 lignes, un pas
d'échantillonnage de 0,125 µs). D(t) est intégré par la méthode des trapèzes sur 400 001 points ;
la vitesse étant linéaire par morceaux, la quadrature n'est inexacte que dans les quatre
intervalles qui contiennent un angle du trapèze.

#### Métriques

| Métrique | Définition |
|---|---|
| e_u(sK) | erreur L² relative du déplacement axial sur [0, T] : ‖u_rockim − u_exact‖₂ / ‖u_exact‖₂ |
| e_F | erreur L² relative de la réaction de l'appui sur [0, T] |
| erreur de célérité | c_mes/c − 1, avec c_mes = (z₃ − z₁)/(t₃ − t₁) et t_k l'instant où u franchit la moitié du déplacement final de l'impulsion (0,8 µm), interpolé linéairement, au premier passage |
| erreur de force | moyenne de RF sur le palier (5 à 15 µs) divisée par ρcv₀A, moins 1 |
| bilan d'énergie | résidu B4 imprimé par le solveur, en pourcentage de l'échelle |
| ordre observé | pente de log e_u(s2) en fonction de log h, par moindres carrés |

#### Critères d'acceptation

Ils ont été fixés dans le script avant le premier calcul et s'appliquent au maillage le plus fin
de chaque variante : e_u(s2) ≤ 2 %, |erreur de célérité| ≤ 1 %, |erreur de force| ≤ 2 %, résidu
du bilan ≤ 1e-6 % sur tous les maillages, ordre observé ≥ 1 quand au moins trois maillages sont
disponibles.

#### Variantes

| Variante | Clés propres | Maillages h (mm) | Attendu |
|---|---|---|---|
| `fem3d` | `law = elastic` | 4 ; 2,83 ; 2 ; 1,41 ; 1 | convergence vers la solution exacte, ordre voisin de 2 pour des tétraèdres linéaires |
| `fdem3d_adaptive` | `insertion = adaptive`, ft = 50 MPa | 4 ; 2,83 ; 2 | aucun joint inséré ; même réponse que le continu |
| `fdem3d_intrinsic_pf20` | `insertion = intrinsic`, `jointPenaltyFactor = 20` | 4 ; 2,83 ; 2 | onde ralentie par la souplesse des joints de pénalité, erreur qui ne décroît pas avec h |

#### Reproduire

```bash
cd ~/rockim
python3 vv/V1_onde/v1_onde.py all --exe build_nofma/rockim
```

L'action `run` lance seulement les calculs, `analyse` seulement le dépouillement ; `--variants` et
`--h` restreignent la campagne.

#### Résultats (2026-10-03)

| Variante | h (mm) | Tétraèdres | e_u(s2) | e_u(top) | e_F | Erreur de célérité | Erreur de force | Résidu B4 | Temps (s) |
|---|---|---|---|---|---|---|---|---|---|
| fem3d | 4 | 3 688 | 0,46 % | 0,64 % | 16,2 % | −0,05 % | −0,04 % | 2,6e-12 % | 0,9 |
| fem3d | 2,83 | 9 824 | 0,28 % | 0,40 % | 11,3 % | −0,03 % | −0,02 % | 2,3e-12 % | 4,2 |
| fem3d | 2 | 20 435 | 0,17 % | 0,25 % | 7,3 % | +0,04 % | < 0,01 % | 7,7e-12 % | 13 |
| fem3d | 1,41 | 48 962 | 0,09 % | 0,15 % | 4,5 % | −0,01 % | < 0,01 % | 3,1e-11 % | 42 |
| fem3d | 1 | 125 486 | 0,05 % | 0,08 % | 3,5 % | < 0,01 % | < 0,01 % | 1,2e-10 % | 167 |
| fdem3d adaptatif | 4 | 3 688 | 0,46 % | 0,65 % | 16,0 % | −0,05 % | −0,04 % | 2,2e-3 % | 26 |
| fdem3d adaptatif | 2,83 | 9 824 | 0,28 % | 0,40 % | 11,1 % | −0,03 % | −0,02 % | 1,1e-3 % | 130 |
| fdem3d adaptatif | 2 | 20 435 | 0,17 % | 0,25 % | 7,2 % | +0,04 % | < 0,01 % | 5,5e-4 % | 539 |
| fdem3d intrinsèque, pf = 20 | 4 | 3 688 | 7,7 % | 13,3 % | 59,5 % | −2,9 % | −2,9 % | 1,5e-3 % | 58 |
| fdem3d intrinsèque, pf = 20 | 2,83 | 9 824 | 7,9 % | 13,6 % | 59,9 % | −3,0 % | −3,0 % | 6,6e-4 % | 368 |
| fdem3d intrinsèque, pf = 20 | 2 | 20 435 | 8,6 % | 14,8 % | 63,8 % | −3,2 % | −3,3 % | 3,5e-4 % | 1 345 |

Valeurs exactes de référence : c = 5 212,47 m/s, plateau de réaction ρcv₀A = 88,403 N. Les temps
sont ceux d'un run avec le nombre de fils OpenMP par défaut (10 cœurs sur ce Mac), avec par moments d'autres calculs en parallèle ; ils
donnent un ordre de grandeur, pas une mesure de performance.

Figures, dans `vv/V1_onde/` : `fig_maillage.pdf` (face y = 0 de la barre et vue 3D des 40
premiers millimètres), `fig_propagation.pdf` (contrainte axiale à sept instants),
`fig_profils.pdf` (déplacement des nœuds contre la solution exacte), `fig_signaux_<variante>.pdf`
(signaux aux stations et réaction de l'appui) et `fig_convergence.pdf`. Les trois premières sont
produites par `vv/V1_onde/viz_v1.py`, sur un run fem3d à h = 2 mm et 25 trames.

#### Verdict et interprétation

fem3d passe les cinq critères. L'erreur en déplacement passe de 0,46 % à 0,05 % entre h = 4 et
1 mm ; l'ordre observé vaut 1,58. Il reste sous l'ordre 2 des tétraèdres linéaires parce que la
vitesse imposée n'est que continue : ses quatre angles limitent la régularité de la solution
exacte. La célérité et le plateau de force sont justes à 0,05 % dès h = 4 mm. L'erreur L² sur la
réaction, de 16 % à 3,5 %, vient des oscillations parasites qui suivent les angles du trapèze au
retour de l'onde (dispersion numérique du schéma explicite, figure `fig_signaux_fem3d.pdf`) ; le
plateau lui-même est exact.

fdem3d adaptatif, sans aucun joint inséré, reproduit fem3d à 1e-4 près sur toutes les métriques
mécaniques, maillage par maillage. Il échoue au seul critère du bilan d'énergie : le résidu vaut
2,2e-3 %, 1,1e-3 % puis 5,5e-4 % de l'échelle, contre 1e-12 % pour fem3d. L'en-tête du bilan
dans `src/Fdem3dSolver.cpp` (l. 7644) en donne la raison : en fdem3d, les compteurs de travail
lisent la vitesse au moment où la force est appliquée, décalée d'un demi-pas dans le schéma
saute-mouton, et le théorème travail-énergie n'y est vérifié qu'à un résidu d'ordre dt près. Le
résidu décroît avec le pas de temps, et les postes se recoupent (à h = 4 mm, travail des liaisons
1,29439e-4 J contre énergie cinétique 6,21855e-5 J plus énergie élastique 6,72442e-5 J). C'est
une erreur de comptabilité et non une création d'énergie. Le critère de 1e-6 %, calibré sur la
comptabilité exacte de fem3d, ne convient pas à fdem3d ; il n'a pas été modifié après coup, et
l'échec est consigné tel quel.

fdem3d intrinsèque avec `jointPenaltyFactor = 20`, la valeur du deck v3P, propage l'onde 2,9 à
3,2 % trop lentement, et la réaction de plateau est trop faible d'autant. Aucun joint ne se rompt :
l'écart est purement élastique. Une célérité réduite de 3 % correspond à un module oedométrique
effectif réduit d'environ 6 %, la souplesse ajoutée par les joints de pénalité, du même ordre que
les 4,5 % d'allongement supplémentaire mesurés en statique sur la barre en traction
(`loads_intrinseque_3d`). L'erreur ne décroît pas avec h (ordre observé −0,16) : le nombre de
joints par unité de longueur croît comme 1/h et compense le gain de précision des éléments. En
percussion, cela signifie qu'avec ce réglage toutes les ondes du bloc arrivent environ 3 % en
retard et que l'impédance ρc, donc la force de contact initiale, est sous-estimée d'environ 3 %,
indépendamment du maillage. Seul le facteur de pénalité règle ce biais ; un balayage de pf à h =
4 mm est le complément naturel de V1.

#### Réserves

- Le champ VTU `pressure` de fem3d vaut tr(σ)/3 (`src/Fem3dSolver.cpp`, l. 2100), positif en
  traction. Il ne s'agit pas d'une pression au sens usuel ; `viz_v1.py` en tient compte.
- Les maillages sont générés par Gmsh 4.13.1 avec la graine 1 ; un autre générateur ou une autre
  graine changerait les valeurs au-delà de la quatrième décimale, pas les conclusions.
- Le pas de temps de fdem3d est commandé par les tétraèdres les plus aplatis (hauteur inscrite
  minimale 0,22 mm à h = 2 mm), d'où des temps de calcul de 40 à 100 fois ceux de fem3d.
- La barre est confinée latéralement. La propagation dans une barre à faces libres, dispersive,
  n'est pas couverte par V1.

## Partie III : validation expérimentale

Aucun banc n'est encore ouvert. Le plan est en section I.5.2.

## Journal

- 2026-10-03 : cadre, inventaire de l'existant et plan (partie I). Banc V1 écrit, lancé et dépouillé (11 runs) ; figures du maillage et de la propagation.
