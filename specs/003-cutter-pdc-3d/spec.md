# Feature Specification: Cutter PDC en 3D

**Feature Branch**: `003-cutter-pdc-3d`

**Created**: 2026-08-18

**Status**: Draft

**Input**: Porter en 3D la géométrie de cutter qui n'existe qu'en 2D — arête de
coupe, angle de coupe arrière, étendue de face, empreinte finie, chanfrein —
avec le contact et ses garde-fous, et pour critère central la reproduction du
cas 2D en déformation plane par un cutter large.

## Contexte

Un essai de coupe en trois dimensions n'est pas possible aujourd'hui. Le
solveur 3D refuse la forme « cutter » et n'accepte qu'une sphère, un poinçon
plat ou l'absence d'outil ; et dans le scénario de coupe il retombe
systématiquement sur la sphère, la forme plate n'étant honorée qu'en
percussion. On ne peut donc rayer la roche qu'à la **bille**, ce qui ne
reproduit pas la dépendance à l'**angle de coupe** — laquelle gouverne la force
de coupe réelle et fait tout l'intérêt de l'essai.

La géométrie existe pourtant, complète et éprouvée, en deux dimensions : le
cutter y est le demi-espace situé **derrière la face de coupe** et **au-dessus
de l'arête**, si bien que la roche sous la profondeur de passe n'est jamais
touchée — exactement ce que laisse un vrai cutter.

Deux acquis chèrement payés en 2D doivent voyager avec la géométrie. D'abord
l'**écrêtage de pénétration sur la taille d'élément locale** : sans lui, un
nœud du copeau qui file au fond du coin produit des pics isolés de **deux
ordres de grandeur** sur une courbe de force par ailleurs régulière. Ensuite le
**frottement régularisé**, sans lequel le contact oscille au changement de sens
du glissement.

Enfin, une anomalie à solder : les réglages de **chanfrein** sont lus et
stockés, mais n'interviennent nulle part dans le calcul du contact. Un réglage
que l'on peut écrire, que le solveur accepte, et qui ne fait rien, est pire que
son absence — il fait croire à un essai qui n'a pas eu lieu.

**Motivation scientifique.** Le banc de rayure de la thèse, monté sous un code
continu, a conclu que la loi de volume employée **n'a pas de régime de coupe**.
Disposer d'un cutter en trois dimensions permettrait de savoir si l'approche
discrète en a un, et de comparer les deux écoles sur le même essai.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Rayer la roche avec un vrai cutter en 3D (Priority: P1)

Fernando monte un essai de coupe en trois dimensions, choisit l'angle de coupe
arrière et la profondeur de passe, lance, et obtient une force de coupe qui
**dépend de l'angle** — la grandeur que l'essai est censé mesurer.

**Why this priority**: c'est la raison d'être du chantier. Sans dépendance à
l'angle, l'essai ne mesure rien de ce qui distingue un cutter d'une bille.

**Independent Test**: livrable utilisable seul — deux essais ne différant que
par l'angle de coupe produisent deux forces distinctes, dans le sens attendu
(un angle plus agressif coupe avec moins d'effort normal).

**Acceptance Scenarios**:

1. **Given** un essai de coupe 3D avec une géométrie de cutter, **When** il
   démarre, **Then** le journal énonce l'angle de coupe, l'étendue de face, la
   profondeur de passe et l'empreinte retenue.
2. **Given** deux essais identiques hormis l'angle de coupe, **When** ils
   s'achèvent, **Then** les forces de coupe diffèrent de façon monotone avec
   l'angle.
3. **Given** une géométrie de cutter demandée dans une dimension ou un scénario
   qui ne la supporte pas, **When** le run démarre, **Then** il s'arrête en
   nommant ce qui est disponible, au lieu de retomber silencieusement sur une
   autre forme.

---

### User Story 2 - Retrouver le cas plan comme limite (Priority: P1)

Fernando vérifie la géométrie et le contact d'un seul coup : il lance en 3D un
cutter **large**, à symétrie latérale, et retrouve la force par unité de
largeur mesurée en déformation plane.

**Why this priority**: c'est le seul contrôle qui valide ensemble la géométrie,
le test de contact et l'écrêtage. Sans lui, un cutter 3D qui « donne des
courbes plausibles » reste invérifiable.

**Independent Test**: comparer les deux forces par unité de largeur sur des
essais dont tout le reste — matériau, loi, maillage, profondeur, vitesse — est
identique.

**Acceptance Scenarios**:

1. **Given** un cutter large et latéralement symétrique en 3D, **When** l'essai
   atteint le régime établi, **Then** la force par unité de largeur reproduit
   celle du cas plan à quelques pour cent près.
2. **Given** le même essai avec un cutter étroit, **When** on compare,
   **Then** la force par unité de largeur s'en écarte — signe que les effets
   de bord tridimensionnels sont bien présents et non gommés.

---

### User Story 3 - Ne rien casser de l'existant (Priority: P1)

Fernando rejoue les essais du dépôt après le chantier et retrouve exactement
les mêmes nombres.

**Why this priority**: principe de bit-neutralité de la constitution. Le
chantier est purement additif — aucune configuration existante ne demande cette
géométrie en 3D — donc l'exigence est atteignable sans compromis.

**Independent Test**: la suite de non-régression, avant et après, au bit près.

**Acceptance Scenarios**:

1. **Given** la suite de non-régression rapide, **When** elle tourne après le
   chantier, **Then** tous les repères rendent des valeurs identiques au bit.
2. **Given** un essai de coupe 2D existant, **When** il est rejoué, **Then** il
   produit la même trajectoire qu'avant.

---

### User Story 4 - Solder le chanfrein (Priority: P2)

Fernando pose un chanfrein et sait, sans lire le code, s'il agit.

**Why this priority**: c'est une correction de fiabilité, pas une capacité
nouvelle ; elle peut suivre les trois premières, mais elle ne doit pas être
oubliée.

**Independent Test**: deux essais ne différant que par la longueur de chanfrein
— soit ils diffèrent, soit le journal a annoncé que le réglage est sans effet.

**Acceptance Scenarios**:

1. **Given** un chanfrein non nul, **When** le run démarre, **Then** ou bien il
   agit sur le contact, ou bien le journal déclare explicitement qu'il est sans
   effet dans cette version.
2. **Given** un chanfrein nul, **When** on compare à l'état actuel, **Then** le
   comportement est inchangé.

### Edge Cases

- Le cutter est **entièrement dégagé** de la roche : la force doit être
  rigoureusement nulle, sans force fantôme d'un outil resté à sa position par
  défaut.
- La **profondeur de passe est nulle** : le cutter effleure la surface sans
  charger.
- L'**angle de coupe est négatif** (face rejetée en arrière) : la géométrie
  reste licite et le contact cohérent, jusqu'à une borne au-delà de laquelle le
  run refuse de démarrer.
- Le cutter est **plus étroit qu'un élément** : le contact ne peut plus être
  résolu ; le run doit le signaler plutôt que produire une force quelconque.
- Un nœud du copeau **pénètre profondément** dans le cutter : la force doit
  rester bornée par l'écrêtage, pas exploser.
- Le cutter **sort du domaine** en fin de course : plus de contact, aucune
  force résiduelle.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Le solveur DOIT offrir en trois dimensions une géométrie de
  cutter décrite par la position de son **arête de coupe**, un **angle de coupe
  arrière**, une **étendue de face**, une **empreinte latérale finie** et un
  **angle de coupe latéral**.
- **FR-002**: La roche située **sous la profondeur de passe** NE DOIT JAMAIS
  être chargée par le cutter, quelle que soit la pénétration.
- **FR-003**: Le contact DOIT n'engager un point de la roche que s'il satisfait
  simultanément trois conditions : être derrière la face de coupe, se trouver
  dans l'étendue de cette face, et tomber dans l'empreinte latérale du cutter.
- **FR-004**: La force de contact DOIT rester **bornée par un écrêtage assis
  sur la taille locale des éléments**, de sorte qu'un point très enfoncé ne
  puisse pas dominer la force totale.
- **FR-005**: Le frottement DOIT être **régularisé** au changement de sens du
  glissement, sans discontinuité.
- **FR-006**: Le solveur DOIT **refuser explicitement**, en nommant les formes
  disponibles, toute géométrie d'outil non supportée par la dimension ou le
  scénario demandés — et ne JAMAIS y substituer silencieusement une autre forme.
- **FR-007**: Le journal de démarrage DOIT énoncer la géométrie retenue et
  toutes ses valeurs effectives, y compris celles laissées par défaut.
- **FR-008**: Les réglages de **chanfrein** DOIVENT soit agir sur le contact
  dans les deux dimensions, soit être déclarés sans effet au démarrage ; le
  silence actuel n'est pas acceptable.
- **FR-009**: Toute configuration antérieure au chantier DOIT produire des
  trajectoires **bit-identiques**.
- **FR-010**: La force exercée sur le cutter DOIT être **déterministe** à
  nombre de fils d'exécution fixé.
- **FR-011**: Le cutter DOIT pouvoir être **piloté en déplacement** (vitesse de
  coupe et profondeur imposées), le mode qui correspond à l'essai de rayure.
- **FR-012**: La référence des clés DOIT documenter la géométrie 3D, ses
  paramètres, ses bornes admissibles et son domaine de validité.

### Key Entities

- **Cutter** : outil rigide défini par son arête de coupe, son orientation
  (coupe arrière et coupe latérale), l'étendue de sa face, son empreinte
  latérale et son chanfrein éventuel.
- **Essai de coupe** : mouvement imposé du cutter à profondeur et vitesse
  données, dont on enregistre les composantes de force et l'endommagement
  produit.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Un cutter large et latéralement symétrique en 3D reproduit la
  force par unité de largeur du cas plan **à mieux que 10 %** en régime établi.
- **SC-002**: Cutter dégagé de la roche : force **exactement nulle** et aucun
  joint rompu.
- **SC-003**: Profondeur de passe nulle : force **exactement nulle**.
- **SC-004**: Les repères de la suite de non-régression rendent des valeurs
  **identiques au bit** avant et après le chantier.
- **SC-005**: Aucun pic isolé de force ne dépasse **cinq fois** la valeur
  médiane du régime établi — l'écrêtage tient (référence : le défaut mesuré en
  2D atteignait un facteur soixante avant correction).
- **SC-006**: La force de coupe varie de façon **monotone** avec l'angle de
  coupe sur au moins trois angles balayés.
- **SC-007**: Un lecteur du journal peut reconstituer intégralement la
  géométrie employée sans lire ni le code ni la configuration.

## Assumptions

- Le périmètre couvre la **géométrie et le contact** du cutter. La loi de
  volume et la loi de joint restent celles déjà disponibles et ne sont pas
  touchées.
- Le cutter est **rigide** : sa déformation propre n'est pas modélisée, comme
  en deux dimensions.
- L'usure du cutter est **hors périmètre**.
- Un seul cutter par essai : les configurations à plusieurs cutters, et donc
  l'interaction entre saignées voisines, ne font pas partie de ce chantier mais
  ne doivent pas être rendues impossibles par sa conception.
- Le cas plan de référence, pour le contrôle central, est un essai de coupe
  **déjà disponible** en deux dimensions, rejoué avec le même matériau, la même
  loi, la même profondeur et la même vitesse.
- La géométrie existe dans les deux dimensions à l'issue du chantier
  (principe III de la constitution) : la 2D est déjà pourvue, seul le chanfrein
  y reste à trancher.

## Amendement du 2026-08-18 — la voie du cylindre court

Remarque de Fernando en cours de rédaction : un cutter PDC réel est **un
cylindre court**, une « pièce de 1 euro », et non un coin extrudé.

Vérification faite, cette géométrie **existe déjà** dans le solveur 3D sous le
nom de poinçon cylindrique, mais elle n'est **honorée qu'en percussion** : dans
le scénario de coupe, le solveur lui substitue la sphère. Le chantier s'en
trouve considérablement réduit — il ne s'agit plus de porter le coin 2D en
trois dimensions, mais :

1. d'**honorer le poinçon cylindrique dans le scénario de coupe** (ce que
   FR-006 exige déjà : ne jamais substituer silencieusement une forme) ;
2. de lui donner une **inclinaison** — l'angle de coupe arrière est
   l'inclinaison de l'axe du cylindre — et un angle de coupe latéral ;
3. de reprendre le test de contact en trois conditions et l'écrêtage.

Cette voie est **plus fidèle au cutter réel** que l'extrusion du coin 2D, et
elle réutilise une géométrie déjà écrite et éprouvée. Les exigences et les
critères de succès de la spec restent valables tels quels : ils portent sur le
comportement observable, pas sur la construction géométrique.
