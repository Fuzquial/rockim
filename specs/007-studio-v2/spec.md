# Spec 007 : interface graphique du code maison FDEM, deuxième version

**Créée** : 2026-10-02
**Statut** : brouillon, à valider par Fernando avant tout code
**Lignée solveur** : g1 (`rockim_g1`)
**Remplace** : spec 006 (studio PySide6 du 2026-08-22, branche `claude/rockim-last-commit-date-h9nw63`),
gelée comme référence, rien n'est supprimé.

---

## 0. Pourquoi une deuxième version

Les deux tentatives précédentes (piste Qt Quick/QML d'août, studio PySide6 + PyVista du 22/08) ont
échoué sur quatre plans à la fois : lenteur, aspect, inadéquation au travail réel, code fragile.
Le diagnostic retenu tient en trois causes.

1. **Partir du solveur au lieu du travail.** La spec 006 part du registre des 256 clés et génère un
   formulaire. Le doctorant ne pense pas en clés mais en essais. Un formulaire de 256 champs
   n'aide pas à monter un triaxial.
2. **Choisir la pile avant de mesurer.** Aucun budget de fluidité n'a été mesuré avant de coder.
   Or le cas de référence (§6) est petit : 13 000 triangles et 19 000 joints par frame. La lenteur
   venait de l'implémentation (rendu piloté depuis Python, matplotlib, VTU ASCII relus), pas du
   volume.
3. **Livrer en une journée.** Tout le code a été écrit le 22/08 sans retour d'usage, puis la branche
   n'a jamais été fusionnée ni maintenue. Le solveur a ensuite évolué sans elle.

Règles qui en découlent :

- **On construit par essais, pas par clés.** Chaque écran correspond à un type d'essai. Les clés
  restent cachées derrière des choix physiques.
- **Le budget de fluidité est mesuré sur prototype avant de choisir la pile** (§7).
- **Une tranche verticale à la fois**, validée par un usage réel avant la suivante.
- **Le solveur ignore l'interface.** g1 lit un `.cfg` et un `.msh`, écrit `out_*`. Les seules
  modifications du C++ sont des extensions de sortie listées au §5, décidées séparément.

---

## 1. Tranche 1 : l'essai triaxial 2D confiné

Récit d'usage. Fernando veut lancer une série de triaxiaux 2D sur Red Bohus pour comparer l'effet
de l'hétérogénéité. Il choisit les dimensions de l'éprouvette, le confinement et la vitesse. Il
choisit un maillage non structuré homogène ou un GBM, avec des tailles de grain uniformes ou
dispersées. Il ajoute éventuellement des joints pré-rompus. Il voit le maillage avant de lancer,
met plusieurs calculs en file, suit leur avancement, puis compare les courbes et regarde les
champs de contrainte et de déformation et les fissures au cours du temps.

La campagne `etude_triax_hetero/` (86 decks, séries A, B, C, F) fait exactement ce travail en
scripts. Elle sert de référence fonctionnelle : la tranche 1 est réussie si elle remplace
`gen_decks.py` + `lancer_lot.py` + `suivre.py` + `depouille.py` pour cette campagne, à résultats
identiques.

---

## 2. Exigences fonctionnelles

### 2.1 Préparer l'essai

| # | Exigence | Traduction g1 (cachée à l'utilisateur) |
|---|---|---|
| P1 | Éprouvette rectangulaire W × H, épaisseur | `W`, `H`, `thickness` |
| P2 | Matériau par préréglage (Bohus calibré Broyden), propriétés modifiables | bloc matériau + joints |
| P3 | Confinement σ₃ et rampe ; confinement nul = UCS | `confiningPressure`, `confiningRamp`, `confineFaces = sides`, `confineGaugeTime` |
| P4 | Chargement axial par platines, vitesse, délai après confinement | `scenario = tension`, `loading = platens`, `pullV < 0`, `pullRamp`, `pullDelay` |
| P5 | Arrêt après le pic, seuil de chute | `ucsStopAfterPeak`, `stopPeakDrop` |
| P6 | Sorties de courbes en déformation | `historyStrains = true` imposé |

La recette triaxiale est figée par l'interface. Le piège de `triax_bohus_gbm.cfg` (oubli de
`loading`, essai en mors collés) devient impossible.

### 2.2 Mailler

| # | Exigence | Traduction g1 |
|---|---|---|
| M1 | **Non structuré homogène** : taille d'élément h | Gmsh `make_unstructured_mesh.py box2d W H h out.msh seed` → `mesh = file` |
| M2 | **GBM** : taille de grain, nombre de phases, fractions, propriétés par phase | `mesh = voronoi`, `grainMesh = delaunay`, `grainMeshRandom = true`, `grainSize`, `phases`, `phase.<n>.*` |
| M3 | **Tailles de grain** uniformes ou dispersées, affinité de taille par phase | `grainSeeding = random`, `grainSizeSpread`, `phase.<n>.grainSize` |
| M4 | Propriétés des joints de grain (règle alpha ou par paire de phases) | `gbAlpha*`, `gb.<a>.<b>.*`, `gbCombine` |
| M5 | **Aperçu du maillage avant lancement**, avec nombre d'éléments par grain (cible 35-90) | rejeu de la génération (solveur en mode aperçu, voir S3) |
| M6 | Graine aléatoire visible et modifiable, pour les répétitions | `seed` |

Décision par défaut D1 : en 2D, « GBM » signifie Voronoï généré par le solveur avec maillage
Delaunay aléatoire intra-grain, et « homogène » signifie maillage Gmsh. g1 ne combine pas les deux
(`mesh = file` refuse `phases`). L'interface présente un seul choix « homogène / GBM » et masque
cette asymétrie.

### 2.3 Discontinuités

| # | Exigence | Traduction g1 |
|---|---|---|
| D1 | Aucune | — |
| D2 | Joints pré-rompus diffus, fraction et graine | `jointPrebrokenFrac`, `jointPrebrokenSeed` |
| D3 | Joints pré-rompus **dessinés** à la souris sur l'aperçu | `preBrokenJoints = "x1 y1 x2 y2; …"` |
| D4 | Familles de plans affaiblis : pendage, espacement, facteur | `weakPlanes`, `weakPlaneFactor`, … |
| D5 | Garde-fous du solveur appliqués avant lancement | `jointResidualMu >= 0`, refus de `jointContactPenalty = adaptive` et de `jointDeath = damage` |

Hors tranche 1 : **les discontinuités ouvertes** (fente avec écartement initial). g1 ne sait pas
les représenter, un joint pré-rompu est fermé et ne garde que son frottement. C'est un chantier
solveur séparé (géométrie de fente dans Gmsh, ou ouverture initiale des joints), à spécifier à part.

### 2.4 Lancer et suivre

| # | Exigence |
|---|---|
| L1 | File de calculs persistante : ajouter, retirer, réordonner, dupliquer avec une variation |
| L2 | Exécution parallèle réglable (défaut 4 jobs × 4 fils, mesure de LANCER.md) ; un run déjà terminé n'est pas relancé |
| L3 | Choix de l'exécutable g1 (`rockim_g1y19.exe` par défaut) mémorisé |
| L4 | Suivi en direct : courbe q-ε, nombre de joints rompus, avancement en temps simulé, estimation de la fin |
| L5 | Arrêt propre d'un run, reprise de la file après fermeture de l'application |
| L6 | Journal stdout conservé par run (le résumé 2D sort sur stdout, pas dans `summary.txt`) |
| L7 | Variation paramétrique simple : une clé physique × liste de valeurs (σ₃ = 0/10/20/40, par exemple) |

### 2.5 Dépouiller

| # | Exigence |
|---|---|
| R1 | Courbes q-ε axial, ε latéral, ε volumique ; convention maison q = σ − σ(fin de consolidation) |
| R2 | Superposition de plusieurs runs, légende par paramètre qui les distingue |
| R3 | Champs de contrainte par élément : σxx, σyy, σxy, von Mises, avec échelle réglable |
| R4 | Champs de déformation (voir S1) |
| R5 | Fissures : joints rompus colorés par mode (traction, cisaillement, pré-rompu), fragments |
| R6 | Navigation temporelle synchronisée entre la courbe et les champs (clic sur la courbe = frame) |
| R7 | Tableau de synthèse de la série : σ₃ atteint, q pic, E sécant, ε au pic, chute post-pic, part intergranulaire (repris de `depouille.py`) |
| R8 | Export figure PDF vectoriel, Computer Modern, fissures en rendu éléments (règle 5 du CLAUDE.md) |

### 2.6 Modeler la loi des joints (ajout du 2026-10-02, demande de Fernando)

Fernando veut « moduler la loi des joints par cliquer-glisser, varier le pic ou la courbe de
montée, pour voir comment se comporte le matériau ». Tout est déjà réglable dans g1 ; il manque
de le voir et de le toucher.

| # | Exigence | Traduction g1 |
|---|---|---|
| J1 | Courbe σ(δ) du mode I, **calculée avec les formules du solveur** (`setJointLengths`, `YanSoftening.hpp`) | — |
| J2 | Poignée du pic : verticale = `ft` ; horizontale = raideur de montée | `ft` ; `insertionPenaltyFactor` (p = k E / h) |
| J3 | Poignée de fin d'adoucissement : l'aire sous la courbe = Gf | `Gf` |
| J4 | Poignée de forme sur l'adoucissement | `jointSoftening = yan`, `yanC` (a, b fixes, champs avancés pour a, b) |
| J5 | Montée linéaire ou parabolique ; adoucissement linéaire ou de Yan | `jointElastic`, `jointSoftening` |
| J6 | Courbe τ(s) du mode II pour une contrainte normale choisie : pic c + tan φ σn, aire GfII | `cohesion`, `frictionDeg`, `gfShearFactor` |
| J7 | Courbe de référence en fantôme (le préréglage), lectures chiffrées : δe, δF, I, Gf, longueur de zone cohésive l = E Gf / ft² et son rapport à la taille d'élément | — |
| J8 | Clés modifiées listées, écrites dans le deck de l'essai | — |
| J9 | **Essai éclair** : traction, compression simple ou triaxial sur le GBM de 10 grains (banc de calibration du 2026-09-07, 9 s par évaluation sous g0, 3,5 % d'écart sur σt/ft avec 255 grains), courbes successives superposées | run g1, lancé par un clic de l'utilisateur |

La montée est 50 à 500 fois plus courte que l'adoucissement (fragile : δe = 15 nm pour une plage
de 7,6 µm) : elle s'affiche dans une loupe à part. J1-J8 ne demandent aucun calcul et réagissent
instantanément. J9 demande une mesure de coût sous g1 (insertion adaptative), soumise au GO de
Fernando, avant de fixer sa taille et sa durée.

---

## 3. Exigences non fonctionnelles, mesurables

Mesurées sur le cas de référence `etude_triax_hetero/out/F7_disc_gbm_P020` : 12 995 triangles,
19 339 joints, 26 frames, 2 004 lignes d'historique, 165 Mo de VTU ASCII.

| # | Critère | Seuil |
|---|---|---|
| N1 | Démarrage de l'application | < 2 s |
| N2 | Première ouverture du run de référence (conversion unique en cache binaire) | < 15 s, en tâche de fond, interface réactive |
| N3 | Réouverture du même run (cache) | < 1 s |
| N4 | Changement de frame | < 100 ms |
| N5 | Zoom et déplacement sur les champs | ≥ 60 images/s |
| N6 | Aucune opération sur le fil de l'interface | > 50 ms |
| N7 | Suivi d'une file de 4 runs en cours | sans gel perceptible |
| N8 | Installation par un collègue sur un poste Windows neuf | une commande documentée, < 10 min |
| N9 | Thème cohérent (clair et sombre), aucune couleur codée en dur hors du fichier de thème | revue de code |
| N10 | Point décimal partout, quelle que soit la locale Windows | test automatisé |

---

## 4. Architecture

Le principe central répond à l'échec « code fragile » : **un noyau sans interface, testable seul**,
et une interface mince par-dessus.

```
noyau (Python, sans dépendance graphique)
├── essais/triaxial2d.py   modèle de l'essai → deck .cfg (+ .msh si homogène)
│                          reprend la logique de gen_decks.py, une seule implémentation
├── maillage.py            appel Gmsh (homogène) ; aperçu GBM (S3)
├── validation.py          garde-fous du solveur + règles maison, avant lancement
├── file.py                file persistante, exécution parallèle (logique de lancer_lot.py)
├── resultats.py           lecture history.csv + VTU → cache binaire par run
└── depouillement.py       q-ε, E sécant, pic, synthèse (logique de depouille.py)

interface (pile choisie au §7)
└── écrans : Essai | Maillage | File | Résultats
```

Conséquences :

- Les scripts de la campagne deviennent des appels au noyau. On ne garde pas deux implémentations.
- Le noyau est testé sans interface : deck généré identique au deck de `gen_decks.py` pour les
  86 decks de la campagne ; synthèse identique à `RESULTATS_<lot>.csv`.
- Les VTU ASCII ne sont lus qu'une fois. Le cache binaire (tableaux numpy) sert ensuite toutes les
  vues.

---

## 5. Extensions côté solveur g1 (décidées séparément, sans effet sur les résultats)

| # | Extension | Motif | Priorité |
|---|---|---|---|
| S1 | Écrire εyy, εxy par élément et le déplacement nodal dans les VTU | R4 : g1 n'écrit aujourd'hui que εxx | tranche 1 |
| S2 | VTU binaires (`format="appended"` encodé) en option | divise la taille et le temps de lecture | tranche 1 bis |
| S3 | `rockim --mesh-only cfg out` : génère le maillage (Voronoï + Delaunay) et l'écrit sans calculer | M5 : aperçu fidèle du GBM | tranche 1 |
| S4 | Résumé de fin de run écrit aussi dans `summary.txt` en 2D | L6 : lecture directe | confort |

Sans S1, R4 se limite à εxx et aux courbes globales. Sans S3, l'aperçu GBM passe par un run tronqué
à la frame 0, ce qui reste exact mais plus lent.

---

## 6. Hors périmètre de la tranche 1

Discontinuités ouvertes ; 3D ; percussion et impact ; essai brésilien ; calibration Broyden ou PSO
intégrée ; édition libre de toutes les clés (un panneau « avancé » en lecture seule montre le deck
généré) ; modes gelés (fem, dem).

Ces sujets forment les tranches suivantes, dans l'ordre décidé à la fin de la tranche 1.

---

## 7. Choix de la pile : trois candidats, tranchés par mesure

| | A. Qt 6 C++ / QML + VTK | B. PySide6 + Qt Quick (QML) | C. Application web locale |
|---|---|---|---|
| Principe | tout natif, vue VTK dans QML | QML pour l'aspect, noyau Python direct | serveur Python local, interface dans le navigateur (TypeScript, rendu WebGL, courbes uPlot) |
| Aspect | très bon, à construire | très bon, à construire | très bon, bibliothèques de composants modernes |
| Fluidité 2D | excellente | bonne si le rendu reste côté GPU | excellente (13 000 triangles = trivial en WebGL) |
| Itération ergonomique | lente (compilation Qt + VTK sous MSVC) | rapide | rapide (rechargement à chaud) |
| Installation équipe | binaire à distribuer, chaîne de build lourde | Python + PySide6 | Python + navigateur ; empaquetage bureau possible plus tard |
| Lien au noyau Python | pont à écrire | direct | API locale |
| Passage à la 3D | VTK natif | VTK via QQuickVTKItem, fragile | vtk.js ou three.js |
| Antécédent | piste d'août abandonnée | — | — |

Recommandation provisoire : **C**, puis B en repli. A est écarté : son coût d'itération est
incompatible avec une ergonomie construite par retours d'usage.

**Protocole de décision (phase 3)** : deux prototypes jetables, C et B, chacun limité à l'écran
Résultats sur le cas de référence (courbe + champ σyy + fissures + curseur temporel). On mesure N1
à N5 et Fernando juge l'aspect sur captures. La pile retenue est celle qui tient tous les seuils ;
à égalité, l'avis de Fernando sur l'aspect tranche.

---

## 8. Jalons

| Jalon | Contenu | Validation qui le ferme |
|---|---|---|
| J0 | Cette spec | Fernando la valide |
| J1 | Prototypes C et B, mesures N1-N5, maquettes des quatre écrans | choix de la pile |
| J2 | Noyau : triaxial2d + validation + file + résultats, tests d'identité sur la campagne | 86 decks identiques, synthèse identique |
| J3 | Extensions S1 et S3 dans g1 | non-régression bit à bit des runs existants |
| J4 | Interface tranche 1 complète | Fernando monte, lance et dépouille une série réelle sans éditeur texte |
| J5 | Installation sur un second poste, guide utilisateur | un collègue l'installe seul |

Chaque jalon a de la valeur seul. Le noyau (J2) sert aux scripts même si l'interface s'arrête là.

---

## 9. Journal

**2026-10-02.** Spec rédigée après inventaire de g1 (scénario triaxial, maillage, GBM, pré-rupture,
sorties, lancement, campagne `etude_triax_hetero`). Cas de référence mesuré : F7_disc_gbm_P020.

**2026-10-02 (J1).** Prototypes C (web, WebGL2) et B (PySide6 + QtQuick3D) de l'écran Résultats,
mesurés sur F7 : voir [j1/RESULTATS_J1.md](j1/RESULTATS_J1.md). C tient N1-N5 (démarrage 0,45 s,
frame 16 ms, 60 i/s) ; B rate N1 (2,7 s) et N8 (environnement de 715 Mo). Recommandation C, en
attente de la décision de Fernando.

**2026-10-02 (décision).** Fernando retient C, l'application web locale. J2 lancé.

Précision sur le critère de J2 « 86 decks identiques » : l'identité est **sémantique**, au sens
du lecteur du solveur (`Config::load` : texte après `#` ignoré, coupure au premier `=`, la
dernière occurrence d'une clé l'emporte). Mêmes clés, mêmes valeurs numériques ; les commentaires
et l'ordre peuvent différer. Une identité au caractère près n'aurait aucun sens physique et
obligerait le noyau à recopier les notes propres à la campagne.

Constat d'inventaire utile à l'interface : dans la campagne, le témoin « homogène » n'est pas
un maillage Gmsh mais une tessellation de Voronoï à une seule phase, de même graine que le GBM.
Le choix présenté à l'utilisateur est donc à deux niveaux : maillage (Voronoï ou Gmsh), puis
nombre de phases (une ou plusieurs, Voronoï seulement).

**2026-10-02 (J2 livré).** Noyau sans interface dans `studio2/noyau/` (essai, materiaux,
geometrie, validation, maillage, file, resultats, depouillement), 267 tests en 27 s
(`python -m pytest studio2/tests`). Verdicts :

- **Decks** : les 81 runs de MATRICE.csv sont réécrits à l'identique sémantique, contre les decks
  sur disque ET contre `gen_decks.deck()` en mémoire. Le dossier `decks/` en compte 86 : les 5 de
  plus (DIAG_gf7, DIAG_phi23, GF15, GF25, GF40) sont des diagnostics écrits à la main, hors
  générateur. Contrôle falsifiant : un écart de 0,01 % sur E est détecté.
- **Dépouillement** : identique à `depouille.py` (mesures, joints, diagnostics) sur les 20 runs de
  `out/` qui ont un historique. Il n'existait aucun `RESULTATS_*.csv` : l'oracle est le script
  lui-même, appelé en mémoire.
- **Validation** : 17 règles, chacune déclenchée par un cas fautif ; la campagne passe sans
  erreur ni alerte. Le nombre d'éléments par grain est estimé avec le coefficient mesuré sur F2
  (367 grains, 12 995 éléments, soit 35,4) ; la formule théorique (33,3) aurait fait échouer la
  règle 35-90 sur une campagne qui la respecte.
- **File** : testée avec un faux solveur (parallélisme borné, ordre, échec, déjà fait, relance
  sans suppression, arrêt, rattachement par PID après fermeture). **Non vérifié** : la survie
  d'un vrai calcul g1 à la fermeture de l'interface, et tout lancement réel de g1 (règle :
  chaque lancement est validé par Fernando, même un banc court).
- **Maillage** : Gmsh `box2d` par le script existant (2,1 s, reproductible octet à octet à graine
  égale). L'aperçu Voronoï passe par un run tronqué (T = 2 µs, une frame) dans la file ;
  l'option solveur S3 reste souhaitable.
- **Cache** : `resultats.convertir` produit un cache identique octet à octet à celui de J1.

**2026-10-02 (premier lancement réel, GO de Fernando).** Aperçu tronqué de F7_disc_gbm_P020
(deck F7 sauf `T = 2e-6`, `frames = 1`), `rockim_g1y19.exe`, 1 job x 4 fils, lancé par la file
du noyau : code 0 en 2,2 s. La frame 0 et ses joints sont **identiques octet à octet** à ceux
du vrai F7 (367 grains, 12 995 éléments, 1 027 pré-rompus dont 424 libres) : le run tronqué
est un aperçu exact du maillage. Capture : [j1/apercu_maillage_F7.png](j1/apercu_maillage_F7.png).
À reprendre en J4 : la biotite et les pré-rompus partagent le violet.

**2026-10-02 (J3, S1 livré).** Clé `writeStrainFields` (défaut false) dans g1 :
`strainXX/YY/XY` par élément (Biot V − I dans le repère GLOBAL, même décomposition polaire
qu'`elementForces`, cisaillement tensoriel) et `displacement` nodal, calculés dans `writeFrame`
seulement. Exécutables : `rockim_j3ref.exe` (HEAD avant J3, même taille que g1y19) et
`rockim_j3.exe`. Vérification (GO de Fernando) : F7_disc_gbm_P020 tronqué à T = 8e-4 s
(confinement 20 MPa puis début de l'axial), 4 frames, 1 fil, trois runs en parallèle de 580 s
(1,2 µs simulée par seconde), script [j3/verif_j3.py](j3/verif_j3.py) :

- A (référence) contre B (j3 sans la clé) : 15 fichiers, tous identiques octet à octet sauf
  `config_effective.cfg`, qui liste la nouvelle clé à son défaut et nomme l'exécutable. Journaux
  identiques hors temps mur. **Le défaut est intact.**
- A contre C (j3 avec la clé) : history.csv, frames.csv, csv finaux et VTU de joints identiques
  octet à octet ; VTU d'éléments identiques sur les 10 tableaux communs, 4 tableaux ajoutés.
  **La clé ne change aucune trajectoire.**
- Contrôles physiques sur la dernière frame : déplacement = x − x0 à 1e-13 m ; strainXX − epsXX
  ≤ 1,4e-6 pour des déformations de 3e-3 (rotations faibles avant rupture) ; moyenne de strainYY
  entre les jauges −1,529e-4 contre epsGauge 1,524e-4 dans history.csv (0,3 %, signe opposé par
  convention : epsGauge est positif en compression) ; strainXX moyen −1,1e-3 sous 20 MPa latéraux.

S3 (`--mesh-only`) abandonné : l'aperçu par run tronqué (J2) est exact au bit et prend 2,2 s.

**2026-10-02 (§2.6, mesure de l'essai éclair, GO de Fernando).** Trois runs construits par le
noyau, `rockim_j3.exe`, 1 fil chacun en parallèle : Voronoï 20 × 40 mm, `grainSize = 0.010`
(10 grains), `grainElemSize = 0.0023` (367 éléments, 36 par grain), jeu fragile une phase,
insertion adaptative, `seed = 12345`. Durées murales : **traction 1,1 s, compression simple
5,5 s, triaxial 20 MPa 16,9 s** (dt = 2,4e-8 s, 33 k / 124 k / 207 k pas). Figure :
[j4/eclair_courbes.png](j4/eclair_courbes.png).

- Traction : pic 1,38 MPa = 1,06 ft (joints de grain non atténués), 7 joints rompus, arrêt
  0,24 ms après le départ.
- Compression : **les maxima bruts sont des pointes parasites** (16 lignes sur 1 103 en UCS, 71
  sur 1 896 en triaxial dépassent de plus de 5 MPa la médiane glissante). Le pic de l'enveloppe
  (médiane glissante, stable de ±3 à ±15 lignes) vaut **48,3 MPa en UCS** (campagne Yan du 11/08 :
  51 MPa, −5 %) et **89,7 MPa à σ₃ = 20 MPa** (F1_homog_P020, même matériau sur 36 × 72 mm :
  93,3 MPa filtré, −4 %). Le banc éclair reproduit donc l'essai pleine taille à 5 % près, environ
  300 fois plus vite.
- Conséquences pour J4 : (1) le dépouillement doit donner un pic FILTRÉ à côté du maximum brut
  (le brut reste, pour l'identité avec depouille.py) ; (2) en traction par mors, history.csv n'a ni
  `epsGauge` ni `nBrokTen` : le dépouillement doit lire `epsAx` et les compteurs du résumé ;
  (3) l'origine des pointes (contact des plateaux sur un maillage grossier ?) reste à établir,
  hors J4.

**2026-10-02 (J4 lancé, maquettes validées par Fernando).** Découpage, chaque lot testé avant
le suivant :

| Lot | Contenu | Test qui le ferme |
|---|---|---|
| J4.0 | noyau : pic filtré, traction par mors | tests unitaires, oracle depouille.py inchangé |
| J4.1 | `studio2/app/` : serveur stdlib (API JSON sur le noyau, file pilotée en tâche de fond), coquille à onglets | tests d'API sur un espace de travail jetable |
| J4.2 | Résultats : liste des runs (espace + dossiers de campagne), vue WebGL de J1, comparaison, synthèse | navigateur, sur les runs existants |
| J4.3 | File : état en direct, arrêt, relance, ordre, courbe et journal du run choisi | faux solveur |
| J4.4 | Essai : formulaire, schéma, validation, deck, ajout, variation de σ₃ | aller-retour formulaire -> deck = noyau |
| J4.5 | Maillage : aperçu par run tronqué, dessin de fissures | GO de Fernando pour le premier aperçu réel |
| J4.6 | Loi des joints + essai éclair | mode II confronté à tools/yan_point.cpp ; GO pour le premier éclair |

Aucun lancement de g1 pendant le développement sans GO : la file est testée avec le faux solveur.

**2026-10-02 (J4 livré).** Les cinq écrans sont en service (`python studio2/app/serveur.py`,
http://localhost:8770 ; guide : [studio2/LISEZMOI.md](../../studio2/LISEZMOI.md)), 380 tests.
Verdicts par lot :

- J4.0 : pic filtré et traction par mors ; oracle depouille.py inchangé.
- J4.1 : serveur stdlib, API testée sur espace jetable avec faux solveur.
- J4.2 Résultats : comparaison et détail WebGL sur les runs de la campagne, déformations affichées
  quand le run a `writeStrainFields`.
- J4.3 File : actions vérifiées par clics réels (ordre côté serveur, suspension des lancements).
- J4.4 Essai : **le formulaire reconstruit les 81 runs de la campagne à l'identique** à partir des
  seuls choix d'écran (tests/test_formulaire.py) ; estimation de coût calée sur les mesures du jour.
- J4.5 Maillage : aperçu réel par run tronqué, reconnu par sa clé de deck sans relancer g1 (testé
  avec le vrai aperçu de F7) ; fissures dessinées à la souris jusqu'au deck. Lien de partage d'essai.
- J4.6 Loi des joints : **la maquette se trompait sur le mode II** (cisaillement ramené à zéro) ;
  le solveur n'adoucit que la cohésion et garde le frottement, τ -> tan φ σn. Loi de référence
  `noyau/loi.py` confrontée à `tools/yan_point.cpp` (copie acceptant σn,
  `tests/donnees/yan_point_sn.cpp`) à 1e-9 près en modes I et II (σn = 0 et 20 MPa), copie
  JavaScript égale à la référence à 1e-12 (Node). Essai éclair identique aux decks validés.

**2026-10-02 (brésilien et banc, GO de Fernando).** Essai brésilien ajouté (noyau, formulaire,
validation, dépouillement, écrans). Préréglage = banc de Yan 2023 calibré, réécrit à l'identique
de `configs_yan/bd_yan_calibre.cfg` (deck du 11/08, commit 04a222f) ; ses trois écarts aux règles
de septembre sont signalés (pas de grainMeshRandom, ~123 éléments par grain, jointPenaltyFactor
inerte). Garde-fou : `historyStrains` fait échouer g1 sur un disque, il n'est jamais écrit.

Premier run lancé PAR L'INTERFACE (file de studio2/espace, rockim_j3.exe, 1 job x 4 fils) :
fini, code 0, 130 s. Confronté aux mesures notées dans le deck :

| | Deck (août) | g1 (j3) |
|---|---|---|
| jauge élastique de bande | 1,001 | 0,992 PASS ; σyy/σxx = -3,10 (théorie -3) |
| arrêt après pic | pas 96 589 / 180 471 | pas 166 623 / 172 585 |
| pic | t = 0,51 ms | t = 0,98 ms |
| σt | ft + ~9 % de surcharge dynamique | 5,30 MPa = 4,1 ft (filtré 4,88) |

**Verdict : l'élastique est exact, la rupture ne reproduit pas le banc.** Au dernier pas intact le
centre porte σxx = 3,3 MPa (2,6 ft) avec 3 joints rompus : l'insertion ne se déclenche pas sous
l'état biaxial du centre (σyy ≈ -3 σxx), alors qu'en traction directe g1 rompt à 1,06 ft (essai
éclair). La courbe oscille fortement dès le début (dampingLocal = 0,1). Les nouvelles clés
d'insertion de g1 (insertionCriterion, facetAverage, insertionHoldSteps) ont des défauts
historiques ; la cause n'est pas établie. Enquête proposée (solveur, hors interface, GO requis) :
recompiler la source du 11/08 et rejouer le même deck pour séparer régression et deck.
Capture : [j4/app_resultats_bresilien.png](j4/app_resultats_bresilien.png).

**2026-10-02 (enquête, GO de Fernando).** Source du 11/08 (commit 04a222f) extraite dans
`FDEM/rockim_g1_0811` (worktree détaché), compilée (`rockim_0811.exe`), même deck, 4 fils :
**résultats identiques à g1 au dernier chiffre** (jauge 0,991762, σt = 5,30023 MPa, arrêt au pas
166 623 / 172 585, 8 194 éléments, 63 grains). Pas de régression entre août et g1. Les valeurs
notées dans le deck viennent d'un AUTRE maillage : leur plafond de 180 471 pas suppose dt =
6,10e-9 s contre 6,37e-9 ici, cohérent avec la non-portabilité connue de la tessellation (CDP-06).
Reste ouvert : sur ce maillage-ci, l'insertion ne se déclenche pas au centre du disque avant 2,6 ft
de traction. Test suivant proposé (GO) : même deck avec grainMeshRandom = true, puis une autre
graine, pour savoir si les 4 ft tiennent au maillage.

**Non vérifié** : aucun aperçu ni essai éclair n'a encore été lancé par l'interface sur g1 (les
seuls runs réels du jour l'ont été par scripts, avec GO) ; la survie d'un calcul g1 à la fermeture
du serveur ; l'aspect en thème clair. La forme du mode II est celle du point matériel, pas d'un
joint dans un maillage.

Hors J2, à décider : remplacer `gen_decks.py` par des appels au noyau (une seule implémentation),
une fois l'interface en service.
