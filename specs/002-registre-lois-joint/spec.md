# Feature Specification: Registre de lois de joint cohésif

**Feature Branch**: `002-registre-lois-joint`

**Created**: 2026-08-17

**Status**: Draft

**Input**: Registre `jointLaw` calqué sur le dispatcher `law`/MatLaw, pour empiler
indéfiniment des lois de joint sans jamais en modifier ni en retirer une seule ;
préréglages par article publié ; règle d'entrée par test falsifiable ; dispatch
hors boucle chaude ; trajectoires existantes bit-identiques.

## Contexte

Le solveur fait **déjà** cohabiter six lois de volume derrière une clé unique
(`law = elastic | dpr | mc | saksala | saksala2011 | dpdfh`), chacune classe
autonome, chacune avec ses clés propres et son repère de non-régression. Ce
modèle fonctionne et n'a jamais cassé une loi en ajoutant la suivante.

La loi de **joint** n'a pas cet équivalent : elle est écrite en dur dans
`jointForces()` et pilotée par trois booléens indépendants (`jointSoftening`,
`jointShearUnload`, `jointFrictionScaled`). Cela fait **2³ = 8 combinaisons**
dont deux ou trois seulement sont documentées et couvertes par la suite. Chaque
loi supplémentaire ajouterait un booléen et **doublerait** la combinatoire :
c'est la trajectoire par laquelle la fissuration — le cœur physique du code —
deviendrait intestable.

L'étude tunnel du 2026-08-17 a rendu le défaut concret : la reproduction de
Wang et al. (2024) demande la loi de Munjiza 2004 (frottement résiduel conservé
pendant l'adoucissement, décharge par retour radial, pénalité absolue), alors
que la configuration courante suit Yan 2023 — un choix qui s'est fait par
combinaison de booléens, **sans être ni nommé ni déclaré nulle part**.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Reproduire un article par son nom (Priority: P1)

Fernando veut reproduire les résultats d'un article publié. Il désigne la loi
de cet article **par son nom** dans la configuration, lance, et retrouve dans
le journal la liste exacte des choix que ce nom implique — sans avoir à
reconstituer une combinaison de booléens ni à deviner ce qu'un défaut recouvre.

**Why this priority**: c'est la raison d'être du chantier. Sans nommage, tout
écart entre deux runs se discute à l'aveugle ; avec lui, l'écart entre deux
lois devient une mesure. C'est aussi la garantie qu'un choix silencieux ne
puisse plus s'installer.

**Independent Test**: livrable utilisable seul — deux configurations ne
différant que par le nom de la loi produisent deux runs comparables et deux
journaux qui énoncent leurs hypothèses.

**Acceptance Scenarios**:

1. **Given** une configuration posant `jointLaw = munjiza2004`, **When** le run
   démarre, **Then** le journal imprime la loi retenue, sa référence
   bibliographique et la valeur effective de chacun de ses paramètres.
2. **Given** deux configurations identiques hormis `jointLaw`, **When** les
   deux runs s'achèvent, **Then** leurs observables sont comparables terme à
   terme et l'écart est attribuable à la seule loi.
3. **Given** une configuration nommant une loi inexistante, **When** le run
   démarre, **Then** il s'arrête immédiatement en énumérant les noms admis.

---

### User Story 2 - Ajouter une loi sans toucher aux autres (Priority: P1)

Un contributeur ajoute une loi de joint (dépendance à la vitesse de
déformation, joints anisotropes, variante d'un article récent). Il écrit une
unité autonome et l'enregistre. Aucune loi existante n'est modifiée, aucune
configuration antérieure ne change de résultat.

**Why this priority**: c'est l'objectif explicite — accumuler indéfiniment. La
valeur se mesure au fait que l'ajout de la n-ième loi ne coûte pas plus cher
que la deuxième et ne met en risque aucune des n−1 précédentes.

**Independent Test**: ajouter une loi factice au registre, vérifier que la
suite de non-régression complète reste au bit près et que la loi factice est
sélectionnable.

**Acceptance Scenarios**:

1. **Given** une loi nouvelle enregistrée, **When** la suite `fast` s'exécute,
   **Then** les 15 repères rendent des valeurs **identiques au bit** à celles
   d'avant l'ajout.
2. **Given** une loi existante et un besoin de la faire évoluer, **When** le
   comportement doit changer, **Then** la modification prend un **nom nouveau**
   et la loi d'origine reste disponible et inchangée.
3. **Given** une loi candidate sans test au point matériel, **When** on tente
   de l'enregistrer, **Then** l'enregistrement est refusé (règle d'entrée).

---

### User Story 3 - Retrouver un run d'il y a six mois (Priority: P2)

Fernando rejoue une configuration écrite avant le chantier. Elle tourne sans
modification et rend exactement les mêmes nombres.

**Why this priority**: la valeur d'archive de la thèse en dépend ; c'est aussi
le principe I de la constitution. Priorité inférieure aux deux premières
seulement parce qu'elle est une contrainte à respecter, non une capacité à
livrer.

**Independent Test**: rejouer un lot de configurations du dépôt avant/après et
comparer les journaux au bit.

**Acceptance Scenarios**:

1. **Given** une configuration antérieure au chantier utilisant les trois clés
   booléennes, **When** elle est rejouée, **Then** elle produit des
   trajectoires bit-identiques et le journal indique vers quelle loi nommée
   ses clés se sont traduites.
2. **Given** une configuration ne posant aucune clé de loi, **When** elle est
   rejouée, **Then** le comportement par défaut est inchangé.

---

### Edge Cases

- Une configuration pose **à la fois** `jointLaw` et l'une des trois clés
  booléennes historiques : conflit signalé et run arrêté, plutôt qu'une
  priorité silencieuse.
- Un préréglage et une clé explicite se contredisent : la clé explicite gagne,
  et le journal signale l'écrasement en le nommant.
- Une loi ne définit pas de branche de mode II : refus à l'enregistrement, pas
  d'appel dans le vide.
- Deux lois portent le même nom de paramètre avec des sens différents : chaque
  loi ne lit que ses propres clés, préfixées par son nom.
- La loi choisie est incompatible avec l'insertion adaptative (par exemple
  parce qu'elle suppose le joint présent dès l'origine) : incompatibilité
  déclarée par la loi et vérifiée au démarrage.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Le solveur DOIT offrir une sélection de loi de joint **par nom**,
  au moins `munjiza2004`, `yan2023` et `linear`, couvrant respectivement la loi
  de l'article de référence, la loi actuellement employée et la branche
  linéaire historique.
- **FR-002**: Chaque loi DOIT être une unité autonome portant son
  adoucissement, son critère de mode II, sa règle de décharge et ses propres
  paramètres ; aucune loi ne DOIT lire les paramètres d'une autre.
- **FR-003**: Les trois clés booléennes existantes DOIVENT continuer de
  fonctionner et se traduire vers une loi nommée, la traduction étant
  imprimée au démarrage.
- **FR-004**: Toute configuration antérieure au chantier DOIT produire des
  trajectoires **bit-identiques** à 1 thread.
- **FR-005**: Le solveur DOIT refuser au démarrage, avec la liste des noms
  admis, toute loi inconnue ; et refuser tout couple de réglages
  contradictoires plutôt que d'en privilégier un silencieusement.
- **FR-006**: Le journal de démarrage DOIT énoncer la loi retenue, sa référence
  bibliographique et la valeur effective de tous ses paramètres, y compris
  ceux laissés à leur défaut.
- **FR-007**: Le solveur DOIT offrir des **préréglages nommés par article**
  (au moins celui de l'étude tunnel en cours et celui de la compression
  uniaxiale de Yan 2023), chacun se dépliant en clés explicites imprimées.
- **FR-008**: Une clé posée explicitement DOIT l'emporter sur la valeur venue
  d'un préréglage, et l'écrasement DOIT être signalé.
- **FR-009**: L'admission d'une loi au registre DOIT être conditionnée à
  l'existence (a) d'un pilote au point matériel démontrant que l'aire sous ses
  branches vaut exactement les énergies de rupture déclarées, et (b) d'au moins
  un repère dans la suite de non-régression.
- **FR-010**: La sélection de loi NE DOIT PAS introduire de coût par joint et
  par pas ; elle est résolue une fois par pas au plus.
- **FR-011**: Chaque loi DOIT déclarer ses incompatibilités éventuelles avec
  les autres capacités (insertion adaptative, mode de décharge, hétérogénéité),
  vérifiées au démarrage.
- **FR-012**: La référence des clés DOIT porter une ligne par loi, avec son
  article source, ses paramètres et son test d'admission.

### Key Entities

- **Loi de joint** : unité nommée décrivant comment un joint perd sa résistance
  — courbe d'adoucissement, critère de mode II, règle de décharge, jeu de
  paramètres propres, référence bibliographique, incompatibilités déclarées,
  test d'admission associé.
- **Registre** : collection des lois disponibles, interrogeable par nom, seule
  autorité sur les noms admis.
- **Préréglage** : ensemble nommé de choix reproduisant un article publié,
  se dépliant en clés explicites et traçables.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Les 15 repères de la suite rapide rendent des valeurs identiques
  au bit avant et après le chantier, à 1 thread.
- **SC-002**: Le contrôle à charge nulle donne 0 joint cassé et un travail de
  contact exactement nul pour **chacune** des lois du registre.
- **SC-003**: Pour chaque loi, l'aire sous les branches de mode I et de mode II
  restitue les énergies de rupture déclarées à mieux que 0,1 %.
- **SC-004**: Le surcoût de temps de calcul du mécanisme de sélection est
  inférieur à 1 % sur un cas de référence d'au moins 100 000 joints, mesuré
  machine au repos à 1 thread.
- **SC-005**: Ajouter une loi supplémentaire ne demande de modifier **aucun**
  fichier portant une loi existante (vérifiable au diff).
- **SC-006**: Deux configurations ne différant que par le nom de la loi
  produisent un écart d'observables **attribuable et chiffré**, là où l'état
  actuel ne permet que de constater une différence.
- **SC-007**: Le journal d'un run permet de reconstituer intégralement la loi
  employée sans lire le code ni la configuration.

## Assumptions

- Le périmètre couvre la loi **de joint** ; les lois de volume gardent leur
  dispatcher actuel, inchangé.
- Trois lois initiales suffisent à valider le mécanisme : celle de l'article de
  référence, celle employée aujourd'hui, et la branche linéaire historique. Les
  lois à venir (dépendance à la vitesse de déformation, joints anisotropes) ne
  font pas partie de ce chantier mais doivent pouvoir s'y ajouter sans le
  rouvrir.
- Le chantier est **bit-neutre par construction** : il ne change aucun défaut
  et ne recalibre aucune référence.
- La règle d'admission s'applique aussi aux trois lois initiales, y compris la
  branche linéaire historique, qui devra donc recevoir son pilote si elle n'en
  a pas.
- Le registre existe dans les deux dimensions (2D et 3D) dans le même chantier,
  conformément au principe III de la constitution.
- Les préréglages ne sont pas des raccourcis de commodité : chacun documente un
  article reproduit et se déplie intégralement, de sorte qu'aucun réglage ne
  reste implicite.
