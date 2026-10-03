# Campagne de vérification et de validation de rockim

Document vivant de la campagne V&V, commencée le 2026-10-03. La démarche reprend celle des codes
FDEM publiés (Guo 2014 pour Solidity, Lisjak 2013 pour Y-Geo, Fukuda et al. 2020, Lei et al. 2016
pour HOSS ; synthèse en I.5) en deux temps, chacun avec des critères d'acceptation fixés avant le
calcul, ce que ces travaux ne font pas.

La partie II, validation physique, confronte rockim à des solutions exactes et à des résultats
théoriques sur les cas canoniques : élément et ondes, contact, rupture, essais de laboratoire,
sensibilités. La partie III, validation par benchmark, reproduit les résultats chiffrés d'articles
publiés sur les mêmes cas, pour comparer rockim aux autres codes. La partie IV confronte un jeu de
paramètres calé une fois à des essais qui n'ont pas servi au calage. Les scripts et les decks sont
dans `vv/`, un dossier par banc.

| Banc | Objet | Référence | Statut | Section |
|---|---|---|---|---|
| P1.1 | onde de compression dans une barre | d'Alembert | fait : fem3d passe ; fdem3d adaptatif passe sauf le bilan ; fdem3d intrinsèque biaisé de −3 % | II.1 |
| P1.2 | plaque trouée en 3D | Kirsch | prévu (fait en 2D, tunnel EDZ) | |
| P1.3 | cube sous chargement biaxial, patch test | élasticité linéaire | prévu | |
| P2.1 | bloc glissant jusqu'à l'arrêt | L = v²/(2µg) | prévu | |
| P2.2 | sphère sur un plan | conservation de l'énergie, Hertz | prévu | |
| P3.1 | fissure pressurisée, sensibilité au maillage | Guo §2.4 ; ténacité K_Ic | prévu | |
| P3.2 | énergie de fissuration | Gf × aire rompue | prévu | |
| P3.3 | fissure en aile | mécanique linéaire de la rupture | prévu | |
| P4.1 | flexion trois points | théorie des poutres | prévu | |
| P4.2 | essai brésilien | Hondros, σt = 2F/(πDt) | prévu | |
| P4.3 | compression polyaxiale | Mohr-Coulomb | fait (`bench_polyaxial`, suite) | |
| P5 | sensibilités : maillage, vitesse de chargement, pénalité | protocole de Guo et Lisjak | pénalité faite (balayage, loi c_eff/c = (1 + 1,24/pf)^(−1/2)) | II.1 |
| B1 | Guo 2014, §2.4 et chapitre 3 | valeurs de Solidity | prévu | |
| B2 | Lisjak 2013, fissure en aile | 7,0 MPa (Y-Geo) | prévu | |
| B3 | Fukuda et al. 2020, frottement et brésilien | Y-HFDEM 3D | prévu | |
| B4 | AbuAisha et al. 2017, hydro-mécanique | Y-Geo | fait (`bench_abuaisha/`) | |
| B5 | Wang et al. 2024, tunnel | MultiFracS | fait (`tunnel_edz/`) | |
| B6 | Yan et al. 2023, UCS adaptatif | article | fait (`ucs_yan_adaptive`) | |
| B7 | Yang et al. 2025 et 2026, impact d'insert | Solidity | partiel | |

Environnement des résultats, sauf mention contraire : macOS 27 arm64, AppleClang 21, OpenMP de
Homebrew (libomp 21.1.6), `-O3 -ffp-contract=off`, binaire `build_nofma/rockim` de la branche
`claude/rockim-version-check-4ugw9e`.

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

### 5. Ce que font les autres codes, et le plan qui en découle

Cinq démarches publiées ont été relues (2026-10-03).

| Source | Élément | Contact | Rupture | Essais de laboratoire | Critère chiffré |
|---|---|---|---|---|---|
| Guo 2014, Solidity 3D | renvoi à Munjiza | renvoi à Munjiza | fissure pressurisée sur cinq maillages, sans convergence ; flexion trois points | brésilien 1,96 MPa pour ft = 3 MPa (−35 %) ; polyaxial, angles jugés sur figure | aucun |
| Lisjak 2013, Y-Geo 2D | renvoi à Mahabadi | renvoi à Mahabadi | fissure en aile 7,0 MPa, dans l'intervalle 5,1 à 14,6 MPa de quatre modèles de LEFM | calage sur UCS et brésilien, validation à paramètres gelés sur le confinement et le litage | aucun |
| Xiang et al. 2009, Y | aucun | bloc glissant, L = v²/(2µg), accord lu sur figure | aucun | aucun | aucun |
| Fukuda et al. 2020, Y-HFDEM 3D | aucun | sphère sur un plan ; bloc glissant | effet du maillage structuré sur le mode de rupture | UCS calé, brésilien contre Hondros ; SHPB calé sur l'essai même | aucun |
| Lei et al. 2016, élément de HOSS | cube biaxial et plaque de Kirsch, convergence mesurée | aucun | aucun | cylindres entaillés, écart de 4 à 5 fois | aucun |

Ces démarches partagent un socle de cas à solution exacte (bloc glissant, sphère sur un plan,
Kirsch, cube), les essais de laboratoire classiques (flexion trois points, brésilien, UCS,
polyaxial) et des études de sensibilité au maillage, à la vitesse de chargement et à la pénalité.
Aucune ne fixe de critère d'acceptation avant le calcul, la plupart des comparaisons sont visuelles
ou portent sur les données du calage, et aucune ne montre la convergence de la rupture en maillage.
Les auteurs de HOSS (Knight et al. 2020) relèvent eux-mêmes l'absence de bancs d'essai standard
pour la rupture des milieux discontinus.

La campagne de rockim suit le même socle en deux temps : la validation physique (partie II,
blocs P1 à P5 du tableau de tête) puis la reproduction chiffrée des articles (partie III, bancs
B1 à B7). La validation sur des essais indépendants du calage (partie IV) vient en dernier : calage
unique sur Red Bohus (UCS, brésilien, triaxial à 20 et 50 MPa), puis prédiction du triaxial à 75
et 100 MPa sans retouche.

## Partie II : validation physique

### II.1 P1.1, onde de compression dans une barre (V1)

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
indépendamment du maillage. Seul le facteur de pénalité règle ce biais ; le balayage ci-dessous le
quantifie.

#### Balayage du facteur de pénalité intrinsèque

La pénalité des joints intrinsèques vaut p = pf·E/h. Un modèle de joints en série avec le continu
donne M_eff = M/(1 + α/pf), soit c_eff/c = (1 + α/pf)^(−1/2), avec α indépendant de h. La
prédiction a été écrite dans `v1_onde.py` avant le balayage, avec α = 1,21 calé sur le seul point
pf = 20, h = 4 mm.

| pf | h (mm) | Erreur de célérité | Prédite | Erreur de force | e_u(s2) | dt (s) | Temps (s) |
|---|---|---|---|---|---|---|---|
| 10 | 4 | −5,54 % | −5,55 % | −5,56 % | 14,2 % | 6,72e-9 | 41 |
| 20 | 4 | −2,92 % | −2,89 % | −2,93 % | 7,7 % | 4,89e-9 | 58 |
| 20 | 2,83 | −3,01 % | −2,89 % | −2,96 % | 7,9 % | 2,59e-9 | 368 |
| 20 | 2 | −3,23 % | −2,89 % | −3,27 % | 8,6 % | 1,47e-9 | 1 345 |
| 50 | 4 | −1,23 % | −1,19 % | −1,23 % | 3,3 % | 3,14e-9 | 87 |
| 100 | 4 | −0,65 % | −0,60 % | −0,64 % | 1,8 % | 2,24e-9 | 123 |
| 200 | 4 | −0,35 % | −0,30 % | −0,34 % | 1,0 % | 1,59e-9 | 181 |
| 200 | 2,83 | −0,35 % | −0,30 % | −0,32 % | 0,95 % | 8,29e-10 | 1 008 |

La prédiction tient à 0,05 point près de pf = 10 à 200 ; l'ajustement sur les huit points donne
α = 1,239. L'erreur de célérité est indépendante de h à pf = 200 (−0,347 et −0,345 %) et croît de
0,3 point entre h = 4 et 2 mm à pf = 20. La force de plateau suit la célérité, l'impédance ρc étant
réduite d'autant. Le pas de temps décroît en 1/√pf (rapport 0,32 entre pf = 20 et 200 pour 1/√10 =
0,316).

La règle de choix qui en découle, pour un biais de célérité visé ε, est pf ≥ α/(2ε) ≈ 0,62/ε :
pf = 62 pour 1 %, pf = 124 pour 0,5 %. Le coût en pas de temps est un facteur √(pf/20) par
rapport au réglage actuel des decks d'impact, soit 1,8 pour pf = 62 et 2,5 pour pf = 124. La figure
`fig_penalite.pdf` trace le retard de célérité en fonction de pf avec la prédiction et
l'ajustement.

#### Réserves

- Le champ VTU `pressure` de fem3d vaut tr(σ)/3 (`src/Fem3dSolver.cpp`, l. 2100), positif en
  traction. Il ne s'agit pas d'une pression au sens usuel ; `viz_v1.py` en tient compte.
- Les maillages sont générés par Gmsh 4.13.1 avec la graine 1 ; un autre générateur ou une autre
  graine changerait les valeurs au-delà de la quatrième décimale, pas les conclusions.
- Le pas de temps de fdem3d est commandé par les tétraèdres les plus aplatis (hauteur inscrite
  minimale 0,22 mm à h = 2 mm), d'où des temps de calcul de 40 à 100 fois ceux de fem3d.
- La barre est confinée latéralement. La propagation dans une barre à faces libres, dispersive,
  n'est pas couverte par V1.

## Partie III : validation par benchmark (reproduction d'articles)

Aucun banc nouveau n'est encore ouvert. Les reproductions antérieures (B4 à B7) sont documentées
dans leurs dossiers et seront résumées ici.

## Partie IV : validation sur essais indépendants du calage

Aucun banc n'est encore ouvert.

## Journal

- 2026-10-03 : balayage de la pénalité intrinsèque (pf = 10 à 200) ; prédiction tenue à 0,05 point.
- 2026-10-03 : démarche en deux temps retenue (validation physique, puis benchmark d'articles), après relecture de Guo, Lisjak, Xiang, Fukuda et Lei.
- 2026-10-03 : cadre, inventaire de l'existant et plan (partie I). Banc V1 écrit, lancé et dépouillé (11 runs) ; figures du maillage et de la propagation.
