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

Hors J2, à décider : remplacer `gen_decks.py` par des appels au noyau (une seule implémentation),
une fois l'interface en service.
