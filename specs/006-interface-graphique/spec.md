# Feature Specification : rockim-studio — interface graphique complète (voie 2)

**Feature Branch** : `006-interface-graphique`

**Created** : 2026-08-22

**Status** : Draft — plan à valider par Fernando avant tout code

**Input** : Décision de Fernando (2026-08-22) : « faire une interface qui me
permette de faire le modèle, qui maille, qui fait les contacts et les CL — un
vrai software ». Voie 2 retenue parmi les trois voies étudiées :
application de bureau PySide6 + PyVista + API Python de Gmsh, sur le modèle
architectural de PrePoMax (contrôleur central, modèle de données, trois vues,
tout le lourd — CAO, maillage, rendu — délégué à des bibliothèques établies).

---

## 0. Principes non négociables

1. **Le solveur ne sait jamais qu'une GUI existe.** rockim continue de lire un
   `.cfg` texte et un `.msh` Gmsh 2.2 ASCII, et d'écrire son dossier `out_*`.
   Aucune modification du C++ n'est requise par ce chantier (les seuls
   correctifs côté solveur, listés au §7, sont des extensions du lecteur MSH
   déjà planifiables indépendamment). Scripts batch, suite de vérification et
   runs cluster restent intacts.
2. **Un seul modèle de données, une seule sérialisation.** Tout ce que la GUI
   affiche vient d'un objet `RockimModel` ; sauver = écrire un `.cfg` (plus,
   si géométrie native, un `.rsg` projet — §3.4) ; ouvrir un `.cfg` existant
   quelconque des 104 du dépôt doit recharger un modèle complet. Le
   **round-trip cfg → modèle → cfg est un test automatisé**, pas une intention.
3. **Ne rien réécrire qui existe.** Géométrie et maillage = API Python de Gmsh
   (noyau OCC) ; rendu 3D et picking = PyVista/VTK ; courbes = matplotlib ;
   dépouillement = les scripts `tools/` existants appelés tels quels.
4. **Périmètre fermé sur les objets de rockim — et centré FDEM** (décision
   Fernando 2026-08-22, alignée sur « rockim = FDEM » du 2026-08-14) : le
   studio est l'interface des modes **`fdem` / `fdem3d`** — joints cohésifs,
   GBM/Voronoï, Weibull, contact potentiel, essais QS (UCS/BD/triaxial),
   hydro, percussion/coupe. Les modes gelés (fem/fem3d/dem/dem3d) restent
   OUVRABLES (round-trip, lancement, courbes — c'est gratuit via le
   registre) mais ne reçoivent AUCUN développement d'interface dédié : pas
   de presets dem3d (glyphes sphères), pas de formulaire spécifique aux
   lois continues — Abaqus/VUMAT est leur maison. On ne construit PAS un
   pré-processeur généraliste : les seuls objets sont ceux que le solveur
   FDEM comprend (bloc/disque/maillage importé, corps par groupes
   physiques, outil analytique ou maillé, phases GBM, joints, CL des
   scénarios existants). Toute demande hors de cette liste passe par une
   révision de la présente spec.
5. **Windows d'abord** (machine de production du doctorant), Linux en CI.
   Python ≥ 3.10, dépendances : PySide6 (LGPL), pyvista + pyvistaqt, gmsh
   (wheel officiel), numpy, matplotlib. Pas de compilation native côté GUI.

---

## 1. Vision produit

Nom de travail : **rockim-studio** (`tools/studio/`, lancé par
`python -m rockim_studio` ou `tools/studio.cmd`).

Une fenêtre principale, trois « vues » commutables (barre d'onglets haute,
comme PrePoMax) partageant la même scène 3D :

```
┌────────────────────────────────────────────────────────────────┐
│  [Géométrie]  [Modèle]  [Calcul]  [Résultats]                  │
├──────────────┬─────────────────────────────────┬───────────────┤
│ Arbre du     │                                 │ Éditeur de    │
│ modèle       │        Vue 3D (PyVista)         │ propriétés    │
│  Géométrie   │   picking faces / cellules /    │ (formulaire   │
│  Matériaux   │        groupes physiques        │  typé de      │
│  Phases GBM  │                                 │  l'objet      │
│  Joints      │                                 │  sélectionné) │
│  Outil       ├─────────────────────────────────┤               │
│  Contact     │ Console : log gmsh / rockim,    │               │
│  CL & scénario│  validation, suivi de run      │               │
│  Sorties     │                                 │               │
├──────────────┴─────────────────────────────────┴───────────────┤
│ Barre d'état : mode, nb éléments, dt estimé, verdict validation│
└────────────────────────────────────────────────────────────────┘
```

Flux utilisateur cible (cas nominal percussion 3D « à la spec 005 ») :

1. **Géométrie** : primitive `Bloc (W×D×H)` + primitive `Insert sphérique
   (R, gap)` → deux volumes nommés `rock` / `insert` (groupes physiques).
2. **Maillage** : champ de taille (h global, hIns, raffinement cylindrique
   rFin/dFin sous l'axe) → bouton *Mailler* → aperçu du maillage, stats
   (nb tets, h_min, dt estimé par la formule du solveur).
3. **Modèle** : matériau Bohus (préréglage), Weibull/effet d'échelle,
   `groupBond.rock.insert = joints` cliqué dans l'arbre Contact,
   `groupVel.insert` saisi, jauges et trackGroups cochés.
4. **CL** : clic sur la face inférieure → *Encastrement* ; clic sur les faces
   latérales → *Absorbant* (la GUI traduit en clés `absorbing` ou groupes).
5. **Calcul** : validation (checklist §5.6), écriture du `.cfg` + `.msh`,
   lancement, suivi live (`history.csv` grâce à `historyFlush`), arrêt propre.
6. **Résultats** : frames VTU dans la même vue 3D (fissures = joints seuillés
   sur `state`, éléments érodés masqués), courbes F-δ / énergies, planche
   des critères de la spec 005 en un clic (réutilise `bench_impact/tools`).

---

## 2. Architecture logicielle

Calque assumé de PrePoMax (voir DeepWiki/PrePoMax), transposé en Python :

```
tools/studio/rockim_studio/
├── __main__.py            # point d'entrée, QApplication
├── app.py                 # MainWindow : onglets, docks, actions, settings
├── controller.py          # LE médiateur : reçoit les intentions de l'UI,
│                          #   mute le modèle, notifie les vues (signaux Qt)
├── model/
│   ├── registry.py        # REGISTRE DES CLÉS (§3.1) : nom, type, défaut,
│   │                      #   portée (modes), bornes, doc courte, groupe UI
│   ├── rockim_model.py    # RockimModel : état complet du cas (§3.2)
│   ├── objects.py         # dataclasses : Material, Phase, Joint, Tool,
│   │                      #   BodyGroup, BC, Scenario, OutputRequest…
│   ├── cfg_io.py          # sérialisation .cfg <-> RockimModel (round-trip)
│   ├── project_io.py      # projet .rsg (json) : modèle + géométrie + chemins
│   └── validate.py        # règles §5.6 (feasibilité l_cz, unités, exclusions)
├── geometry/
│   ├── gmsh_service.py    # session gmsh unique (initialize/finalize),
│   │                      #   construction OCC, groupes physiques, maillage,
│   │                      #   export MSH 2.2 ASCII, exécution HORS processus
│   │                      #   UI (QProcess vers un worker — §6 R2)
│   ├── primitives.py      # bloc, disque, cylindre, bloc+insert (bench1/1g),
│   │                      #   tunnel — paramétrés, re-générables
│   └── sizefields.py      # champs de taille : global, boule/cylindre, rampe
├── views/
│   ├── tree.py            # arbre du modèle (QTreeView + modèle Qt)
│   ├── props.py           # ProperyPanel : formulaire auto-généré depuis
│   │                      #   registry.py (spinbox borné, combo enum, unité)
│   ├── scene.py           # QtInteractor PyVista : affichage géométrie /
│   │                      #   maillage / résultats, picking, surbrillance
│   ├── console.py         # console de logs + validation
│   └── plots.py           # courbes matplotlib (history, F-δ, énergies)
├── run/
│   ├── runner.py          # QProcess rockim(.exe), env OMP_NUM_THREADS,
│   │                      #   file d'attente, arrêt (terminate), code retour
│   └── monitor.py         # tail de history.csv + summary.txt -> signaux
├── results/
│   ├── vtu_series.py      # lecture des frames (pyvista.read), cache,
│   │                      #   time-slider, seuillages standard (state, bulkD)
│   └── postpro.py         # ponts vers tools/plot_*.py, crater_metrics,
│                          #   bench_impact/tools (imports, pas de subprocess)
└── tests/                 # pytest + pytest-qt (§5)
```

Décisions structurantes :

- **Contrôleur unique + signaux Qt**, pas d'accès direct vue → modèle. C'est
  ce qui a permis à PrePoMax de grossir sans se noyer ; on garde la même
  discipline dès la première ligne.
- **gmsh dans un processus séparé** (worker lancé par QProcess qui reçoit un
  script de construction JSON et renvoie le `.msh` + un `.vtp` d'aperçu).
  Motif : gmsh est mono-session par processus, ses crashs OCC ne doivent pas
  emporter la GUI, et le GIL bloquerait l'UI pendant un maillage long.
- **La scène 3D est unique** et change de contenu selon la vue active
  (géométrie B-Rep tessellée / maillage / résultats), comme PrePoMax — pas
  trois widgets VTK.
- **Undo/redo** par pile de commandes sur le contrôleur (chaque mutation du
  modèle est une commande sérialisable). Livré en M3, mais l'API commande
  est en place dès M0 (sinon c'est irrattrapable).

---

## 3. Modèle de données

### 3.1 Le registre des clés — la pièce maîtresse

Un module `registry.py` déclarant CHAQUE clé de config sous forme :

```python
Key("jointXi", float, default=0.05, scope={"fdem", "fdem3d"},
    bounds=(0.0, 1.0), group="Joints",
    doc="ratio d'amortissement du dashpot de joint",
    advice={"vérification": 0.0, "quasi-statique": 0.01, "impact": 0.05})
```

Source : DOCUMENTATION_rockim.md §5 (≈ 150 clés) croisée avec `Config.cpp`
(vérité du parseur). Le registre pilote : le formulaire de propriétés
(généré, jamais écrit à la main), la validation, la complétion, l'infobulle
(colonne « rôle » de la doc), et le diff config (« quelles clés diffèrent des
défauts »). **T1 livre un extracteur semi-automatique** (`tools/studio/dev/
extract_keys.py`) qui parse `Config.cpp` et signale toute clé du C++ absente
du registre — c'est lui qui empêche la GUI de dériver quand le solveur bouge
(brancher au tier fast de `verify_suite.py`).

### 3.2 RockimModel

État complet d'un cas : `mode`, `scenario`, bloc matériau, listes typées
(phases, corps, CL, jauges, trackGroups, groupBond), objet outil, blocs
joints/contact/hydro, sorties, et un dict `extra` pour toute clé reconnue par
le registre mais sans objet dédié (garantie de round-trip : rien n'est
perdu, même une clé exotique). Trois invariants vérifiés en continu :

- exclusions du solveur (`law` ⊻ `phases` ; `mesh=file` + `phases` ⇒ groupes
  nommés obligatoires ; `scenario=brazilian` ⇒ `geometry=disc`…) ;
- unités SI, point décimal (le piège n°1 de la maison — locale FR) ;
- cohérence corps/groupes : tout `groupBond`/`groupVel`/`gauge` référence un
  groupe physique existant du `.msh`.

### 3.3 Deux régimes de géométrie, assumés dès le départ

- **Régime « clés »** (`mesh = grid | voronoi`) : PAS de géométrie
  explicite — le solveur génère lui-même grille ou Voronoï/GBM. La vue
  Géométrie montre un aperçu paramétrique (boîte/disque + graine de Voronoï
  en 2D via `Tessellation` rejoué en Python léger, ou simple silhouette).
  La GUI n'appelle pas gmsh dans ce régime.
- **Régime « géométrie native »** (`mesh = file`) : construction OCC dans
  gmsh (primitives §2, plus tard STEP), groupes physiques nommés, champs de
  taille, export MSH 2.2 ASCII. C'est le régime du « vrai software » et de
  la spec 005.

### 3.4 Fichiers

| extension | contenu | statut |
|---|---|---|
| `.cfg` | la config solveur — SEULE vérité côté calcul | inchangé |
| `.msh` | maillage Gmsh 2.2 ASCII | inchangé |
| `.rsg` | projet studio (JSON) : recette géométrique paramétrique, chemins cfg/msh/out, état UI (vue, sélections) | nouveau, jamais lu par le solveur |

Ouvrir un `.cfg` seul = projet dégradé sans recette géométrique (le `.msh`
est affiché mais non régénérable) — utile pour reprendre les 104 configs
existantes immédiatement.

---

## 4. Découpage en work packages

### M0 — Socle piloté par le registre (le « PrePoMax minimal »)

- **WP0.1** Squelette application : MainWindow, docks, contrôleur, signaux,
  settings persistants (géométrie fenêtre, exe rockim, threads).
- **WP0.2** `registry.py` + extracteur `extract_keys.py` + revue manuelle
  contre DOCUMENTATION §5.
- **WP0.3** `RockimModel` + `cfg_io.py` + **test de round-trip sur les 104
  configs du dépôt** (lu → écrit → relu ⇒ modèles identiques ; le fichier
  réécrit peut différer en ordre/commentaires, pas en sémantique).
- **WP0.4** Arbre du modèle + formulaire de propriétés auto-généré + console
  de validation.
- **WP0.5** `runner.py` + `monitor.py` : lancer, suivre (courbe live depuis
  `history.csv`), arrêter, rejouer ; file d'attente séquentielle (règle
  maison : UN gros job à la fois) ; bouton « suite de vérification »
  (reprend `VERIF_CONFIGS` de l'actuel `rockim_gui.py`).
- **Critère de sortie M0** : ouvrir `configs/fdem3d_percussion.cfg`, tout
  éditer au formulaire, relancer, suivre la courbe — sans toucher un éditeur
  texte. L'actuel `rockim_gui.py` est alors gelé (pas supprimé — règle n°6
  du CLAUDE.md racine).

### M1 — Vue 3D et résultats (le gain visible)

> **Budget de performance (ajout 2026-08-22, suite à la revue de Fernando —
> la lenteur et la laideur du tkinter actuel ne doivent PAS se reproduire).**
> Critères mesurés à la fin de M1, bloquants pour le jalon :
> démarrage de l'application < 3 s ; chargement du run de référence spec 005
> (~150 k tets × 50 frames) < 10 s avec navigation temporelle ensuite < 0,5 s
> par frame (cache) ; rotation/zoom fluides (> 30 fps) sur la frame la plus
> lourde ; aucune opération UI > 100 ms sur le thread principal (tout le
> long passe en worker). Rappel d'architecture qui rend ces chiffres
> atteignables : le rendu est fait par VTK (C++/GPU), le maillage par Gmsh
> (C++/OCC), les données par les lecteurs VTK binaires + numpy — Python
> n'est jamais dans une boucle par élément. Si un composant précis rate son
> budget malgré profilage, la parade est un module C++ ciblé lié par
> pybind11 (même chaîne oneAPI que le solveur), PAS une réécriture.

- **WP1.1** `scene.py` : QtInteractor PyVista, thèmes, axes, échelle.
- **WP1.2** `vtu_series.py` : slider temporel, presets d'affichage **FDEM**
  (joints rompus en surbrillance sur bulk translucide, champs bulk au
  choix, `bulkD` quand la pulvérisation est armée, fragments). Les autres
  modes s'affichent avec le rendu générique, sans preset dédié (principe
  n°4).
- **WP1.3** `plots.py` : historiques multi-colonnes, F-δ (reprend
  `plot_force_penetration.py`), énergies ; superposition de deux runs.
- **WP1.4** Aperçu du maillage d'un `.msh` (régime file) et silhouette
  paramétrique (régime clés) dans la vue Modèle.
- **Critère de sortie M1** : dépouiller un run de la spec 005 entièrement
  dans studio (fissures 3D + F-δ + énergies), sans ParaView pour le
  quotidien (ParaView reste l'outil des figures finales).

### M2 — Géométrie native et maillage (le cœur « vrai software »)

- **WP2.1** `gmsh_service.py` en worker hors processus + protocole JSON.
- **WP2.2** Primitives paramétriques : bloc, bloc+insert sphérique
  (équivalents `bench1`/`bench1g`), disque, cylindre, tunnel — chaque
  générateur de `make_unstructured_mesh.py`/`make_impact_mesh.py` devient
  une recette studio, ET les scripts CLI restent (ils appellent le même
  module — une seule implémentation).
- **WP2.3** Champs de taille interactifs (global, boule, cylindre, rampe) +
  stats de maillage + **estimation de dt et verdict de faisabilité
  `2·dx < l_cz < a`** (GUIDE §4.2 bis) affichés AVANT de lancer — c'est la
  fonctionnalité qui aurait évité les « six runs de percussion perdus » du
  2026-08-17.
- **WP2.4** Groupes physiques : nommage, couleur, liste des corps dérivée du
  `.msh`, synchronisation avec `groupBond`/`groupVel`/`gauge`/`trackGroups`.
- **Critère de sortie M2** : reconstruire le banc Yang (bloc + insert,
  maillage gradué 127 k tets) entièrement dans studio et retrouver les
  chiffres du Run 3 (même seed ⇒ même maillage, même history.csv).

### M3 — Picking, CL par sélection, confort

- **WP3.1** Picking PyVista : faces/arêtes de la géométrie B-Rep et
  faces/cellules du maillage ; surbrillance, sélection par boîte.
- **WP3.2** Affectation par sélection : encastrement, absorbant
  (sides/all), mors/platines des essais QS, vitesses initiales de corps —
  la GUI traduit sélection → clés/groupes, et matérialise l'inverse
  (cliquer une CL dans l'arbre surligne ses faces).
- **WP3.3** Undo/redo actif, duplication d'objets, bibliothèque de
  matériaux (préréglages : Bohus carte DP-DFH, Kuru Grey, St Anne…),
  modèles de cas (templates = les configs de référence du dépôt).
- **WP3.4** Études paramétriques simples : varier 1-2 clés en grille,
  générer N cfg + file d'attente, tableau comparatif des summary.txt
  (pont vers `calibration_redbohus/tools` sans les remplacer).
- **Critère de sortie M3** : monter un essai UCS/BD/triaxial complet à la
  souris, CL comprises, sans connaître le nom des clés.

### M4 — Durcissement et diffusion

- **WP4.1** Packaging : environnement `pip` verrouillé
  (`tools/studio/requirements.txt`), lanceur `.cmd`, et un exe PyInstaller
  optionnel pour la machine de bureau.
- **WP4.2** CI headless Linux (pytest-qt + VTK off-screen) branchée sur le
  même workflow que la suite du solveur.
- **WP4.3** Documentation : `GUIDE_studio.md` (dans le style de
  GUIDE_rockim.md) + section dans DOCUMENTATION_rockim.md §3.4.
- **WP4.4** Reprise du backlog d'ergonomie accumulé pendant M0-M3.

Ordre de valeur volontaire : M0 et M1 rendent service dès les premières
semaines même si M2/M3 glissent ; l'inverse (géométrie d'abord) laisserait
des mois sans bénéfice quotidien.

---

## 5. Stratégie de test

| niveau | contenu | outil |
|---|---|---|
| unitaire | registry (types/bornes), cfg_io round-trip 104 configs, validate (cas des exclusions), primitives gmsh (volumes/groupes attendus) | pytest |
| intégration | worker gmsh (recette → .msh relu par meshio/pyvista : nb éléments, groupes) ; runner sur `verify_fdem_tension.cfg` (run court réel, verdict PASS lu) | pytest, exe rockim requis |
| UI | fumée headless : ouvrir/éditer/sauver un cfg, changer de vue, picking simulé | pytest-qt, VTK off-screen |
| non-régression solveur | inchangée — la GUI n'y touche pas ; ajout d'un repère « cfg écrit par studio ≡ cfg de référence » pour 3 cas types | verify_suite.py |
| or | même seed, même maillage : le `.msh` de WP2.2 est comparé octet à octet à celui de `make_unstructured_mesh.py` pour bench1/bench1g | pytest |

Règle : **aucun WP n'est « fait » sans son test** — même discipline que les
repères du solveur.

## 5.6 Validation métier embarquée (checklist pré-run)

Portée dans `validate.py`, affichée avant chaque lancement, bloquante ou
avertissante :

- faisabilité percussion : `2·dx < l_cz = E·Gf/ft² < a` (rayon de Hertz
  calculé de l'outil et du matériau) — avertissement argumenté sinon ;
- dt estimé et durée de run extrapolée (à partir des perf notées dans la
  doc) ;
- `jointXi` conforme à la règle maison selon le scénario ;
- absorbant à ≥ une longueur de fissure attendue de la zone de process
  (règle du README) ;
- Weibull : `jointSizeEffect` sans `jointZeff` physique déclaré ⇒ alerte ;
- locale/point décimal, unités, clés inconnues du registre.

---

## 6. Risques et parades

| # | risque | parade |
|---|---|---|
| R1 | Compat VTK/PySide6 sur Windows (versions qui se marchent dessus) | versions épinglées dans requirements.txt, testées en binôme pyvistaqt ; fallback PyQt6 si une release PySide6 casse |
| R2 | gmsh instable/mono-session dans le processus GUI | worker hors processus dès WP2.1 (décision d'architecture, pas un correctif) |
| R3 | dérive registre ↔ Config.cpp | extracteur automatique + repère dans le tier fast |
| R4 | scope creep (« et si on ajoutait l'import STEP, le remaillage adaptatif… ») | principe n°4 : périmètre = objets du solveur ; toute extension passe par une révision de spec |
| R5 | gros VTU 3D lents à charger (runs de 150 k tets × 50 frames) | cache décimé pour la navigation, chargement pleine résolution à la demande ; seuils appliqués côté lecture |
| R6 | le chantier GUI cannibalise le temps de thèse | jalons M0/M1 courts et utiles seuls ; M2/M3 planifiés dans les creux (runs longs qui tournent) ; l'IA fait le gros du code, Fernando valide l'ergonomie |
| R7 | round-trip impossible sur configs exotiques (commentaires, astuces) | dict `extra` conservatif + les commentaires du cfg source archivés en tête du fichier réécrit |
| R8 | la pile Python déçoit malgré tout (perf ou rendu) | **clause de sortie** : le budget de perf M1 est bloquant ; composant fautif remplacé par un module C++/pybind11 ciblé ; au pire, le principe n°1 (solveur ignorant de la GUI) et les acquis transposables (registre des clés, format projet, recettes gmsh, mapping groupes → CL) permettent une GUI Qt/C++ ultérieure sans rien perdre du travail conceptuel. Le tout-C++ d'emblée est rejeté : itération d'ergonomie 3-5× plus lente, chaîne de build Qt+VTK+OCC lourde, et aucun gain là où ça compte (rendu, maillage, données sont déjà en C++ dans la pile retenue — PrePoMax lui-même est en C#, pas en C++) |

---

## 7. Petites extensions côté solveur (optionnelles, hors périmètre GUI)

À décider séparément, aucune n'est bloquante pour M0-M1 :

- **S1** `rockim --print-keys` : dump machine-lisible du registre de clés
  depuis `Config.cpp` (remplacerait l'extracteur R3 par la vérité même) ;
- **S2** groupes physiques 2D dans le lecteur MSH (aujourd'hui documentés
  « V1, 3D ») pour amener le régime géométrie native aux cas 2D ;
- **S3** un `--dry-run` qui parse, valide, imprime dt/l_cz et sort — la
  checklist §5.6 s'appuierait sur le solveur lui-même.

---

## 8. Estimation et jalons

Hypothèse : développement porté par Claude Code en sessions dirigées,
validation ergonomique par Fernando (~2-3 h/semaine), en parallèle des runs
de thèse.

| jalon | contenu | durée estimée | cumul |
|---|---|---|---|
| M0 | socle registre + cfg round-trip + run manager | 2 semaines | 2 sem |
| M1 | vue 3D + résultats | 2 semaines | 4 sem |
| M2 | géométrie native + maillage + faisabilité | 3 semaines | 7 sem |
| M3 | picking + CL à la souris + confort | 3 semaines | 10 sem |
| M4 | durcissement, CI, doc, packaging | 1-2 semaines | ~12 sem |

Soit **environ trois mois calendaires à mi-régime**, avec un outil déjà
utile au quotidien à partir de la fin de M0 (2 semaines). Point GO/NO-GO à
la fin de chaque jalon : chaque M est conçu pour avoir de la valeur même si
le chantier s'arrête là.

## 9. Prochaines actions

1. Validation de cette spec par Fernando (périmètre, priorités M2 vs M3,
   nom `rockim-studio`).
2. WP0.1 + WP0.2 : squelette + registre + extracteur (première PR).
3. WP0.3 : round-trip sur les 104 configs — le premier verdict chiffré du
   chantier.
