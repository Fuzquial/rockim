# Feature Specification: Pilier performance — l'architecture de calcul du banc percussif

**Feature Branch**: `001-pilier-performance-architecture`

**Created**: 2026-08-14

**Status**: Draft — en attente des clarifications C1-C4 (Fernando)

**Input**: User description: "s'attaquer au pilier absent : l'architecture de
calcul (revue biblio 2026-08-13 : rockim au niveau 2023-2025 sur 3 piliers
sur 4 — lois cohésives, contact moderne, multi-corps — le pilier manquant est
le HPC : tous les codes de référence, Y-HFDEM/Fukuda, HOSS, Irazu, ont
franchi GPU ou MPI ; rockim est OpenMP mono-nœud). Mesure P1 → décision P2
(CPU / MPI / GPU) → implémentation de l'issue choisie."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Connaître le coût réel avant de choisir (P1 — la mesure qui commande tout)

Fernando dispose d'une **charge de référence** représentative du banc de
Mines Paris (géométrie, taille de maille, durée physique) et d'un **profil de
coût mesuré** : temps mur total et répartition par module (éléments, joints,
contact, E/S, reste) sur cette charge, extrapolé à la taille réelle du banc.
La décision d'architecture se prend sur ces chiffres, pas sur des intuitions.

**Why this priority**: toute la suite (P2) est conditionnée à cette mesure —
décider avant de mesurer serait de l'architecture sur plan ; c'est la doctrine
du projet (constitution, principe IV).

**Independent Test**: lancer la charge de référence sur la machine de mesure
→ le profil (ROCKIM_PROF + wall) sort, l'extrapolation à la taille du banc
réel est documentée avec ses hypothèses (scaling observé, pas supposé).

**Acceptance Scenarios**:

1. **Given** la charge de référence définie et gelée (configs versionnées),
   **When** elle tourne machine idle à 1 thread puis aux threads disponibles,
   **Then** le rapport P1 donne : temps mur par cas, profil par module,
   scaling OpenMP observé, et l'extrapolation chiffrée au banc réel avec
   incertitude.
2. **Given** le rapport P1, **When** on applique la grille de décision P2
   (voir FR-003), **Then** une et une seule issue (a/b/c) est désignée, avec
   les chiffres qui la justifient.

---

### User Story 2 - Des runs longs qui survivent (P1 bis — le prérequis de toute issue)

Un run de plusieurs heures (banc complet, calibration) peut être interrompu
(panne, fenêtre de calcul, redémarrage) et **reprendre où il en était** sans
perte : l'état complet du solveur est sauvegardé périodiquement et le run
relancé depuis le dernier point de contrôle produit la même trajectoire que
le run ininterrompu.

**Why this priority**: indispensable quelle que soit l'issue P2 (le plan v2
le classe D7 « indispensable dès les runs de dizaines d'heures ») ; sans lui,
la calibration V4 (des dizaines de runs) est un jeu de roulette.

**Independent Test**: un run de référence coupé à mi-course puis repris
termine avec les MÊMES résultats (bit-identiques à 1 thread) que le run
d'une traite — c'est le critère, vérifiable sur un cas court.

**Acceptance Scenarios**:

1. **Given** un run avec point de contrôle activé (opt-in), **When** il est
   tué à t = T/2 puis relancé depuis le checkpoint, **Then** résumé final,
   history.csv (à partir de la reprise) et compteurs d'énergie sont
   identiques au run témoin d'une traite (1 thread).
2. **Given** une config SANS la clé checkpoint, **When** elle tourne,
   **Then** trajectoires bit-identiques à l'état actuel (constitution, I).

---

### User Story 3 - Des sorties qui ne freinent pas le calcul (P2)

Sur les gros maillages, l'écriture des frames (VTU ASCII aujourd'hui) ne
domine plus le temps de run : les sorties passent en binaire (opt-in) et
leur poids dans le temps mur devient marginal, sans rien perdre de la
lisibilité ParaView.

**Why this priority**: mesuré comme « dominant en E/S sur gros maillages »
(roadmap D1) ; devient bloquant à l'échelle du banc réel, mais ne conditionne
pas la décision P2.

**Independent Test**: sur la charge de référence, part de l'E/S dans le
temps mur avant/après ; fichiers relisibles dans ParaView, valeurs identiques.

**Acceptance Scenarios**:

1. **Given** la charge de référence avec `vtkBinary = true`, **When** elle
   tourne, **Then** l'E/S ≤ 5 % du temps mur et les champs relus (ParaView /
   scripts) sont identiques aux valeurs ASCII à la précision d'écriture.

---

### User Story 4 - Le solveur exploite la machine choisie (P2 — l'issue élue de la décision)

Selon l'issue désignée par P2, le banc complet tourne dans le budget temps
de Fernando : l'accélération est mesurée sur la charge de référence,
la physique est inchangée (mêmes trajectoires aux tolérances déclarées de
l'issue), et la suite de non-régression passe sur l'architecture cible.

**Why this priority**: c'est le pilier lui-même — mais il ne peut être
spécifié finement qu'après P1/P2 ; cette user story sera raffinée en spec
dédiée (002-...) une fois l'issue choisie.

**Independent Test**: la charge de référence sous l'architecture élue tient
le critère de budget temps (voir SC-005) ; suite fast/full verte sur la
cible ; zeroload et bilan d'énergie inchangés.

**Acceptance Scenarios**:

1. **Given** l'issue P2 implémentée, **When** la charge de référence tourne,
   **Then** le temps mur satisfait le budget déclaré en C2 et le rapport
   accélération/coût est documenté.
2. **Given** le mode accéléré OFF (défaut), **When** une config existante
   tourne, **Then** bit-identique à l'état actuel (constitution, I).

---

### Edge Cases

- Interruption PENDANT l'écriture d'un checkpoint : le point de contrôle
  précédent reste valide (écriture atomique — fichier temporaire + rename).
- Checkpoint d'un run multi-corps avec contact par potentiel : l'état des
  paires (historiques de frottement, relèves de naissance, caches) fait
  partie de l'état sauvegardé — sinon la reprise diverge silencieusement.
- Charge de référence sur une machine différente de celle de la mesure P1 :
  les chiffres P1 sont PAR MACHINE ; l'extrapolation le déclare.
- Issue GPU avec double précision indisponible/lente sur la carte cible :
  la tolérance de trajectoire de l'issue doit être déclarée AVANT (une somme
  flottante réordonnée n'est pas bit-identique — c'est une décision non
  bit-neutre au sens de la constitution).
- Suite de non-régression sur architecture parallèle : les références
  actuelles sont certifiées à 1 thread — l'issue choisie devra définir son
  équivalent (référence par configuration matérielle ?).

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: le dépôt DOIT définir une **charge de référence P1** gelée et
  versionnée : ensemble de configs représentatives du banc (impact insert
  V1, percussion longue pénalité ET potentiel, un cas à maille fine), avec
  [NEEDS CLARIFICATION C1 : dimensions et taille de maille du banc RÉEL de
  Mines Paris — éprouvette, insert/taillant, durée physique à simuler —
  pour dimensionner le cas extrapolé].
- **FR-002**: le solveur DOIT produire un profil de coût par module
  (ROCKIM_PROF existant, complété si un module n'est pas couvert) et le
  rapport P1 DOIT consigner : wall par cas, profil, scaling OpenMP mesuré
  (1, 2, N threads), mémoire crête.
- **FR-003**: la grille de décision P2 DOIT être appliquée telle quelle :
  (a) banc réel ≲ quelques heures sur la machine de Fernando → issue CPU
  (confort : D1 + D7 + OpenMP fem/dem/dem3d) ;
  (b) banc en jours ET accès cluster CPU → issue MPI (décomposition de
  domaine à la HOSS) ;
  (c) banc en semaines OU besoin paramétrique/temps réel → issue GPU
  (portage CUDA des noyaux, recettes Fukuda 2019-2024, C&G 2024b) ;
  avec [NEEDS CLARIFICATION C2 : budget temps acceptable par run de banc —
  heures ? une nuit ? — et nombre de runs paramétriques visé pour V4] et
  [NEEDS CLARIFICATION C3 : matériel accessible — GPU (modèle, où : poste
  labo, cluster école, cloud ?), cluster CPU (cœurs, MPI dispo ?), et la
  machine Windows de F1 est-elle la machine de production ?].
- **FR-004**: le restart/checkpoint DOIT être opt-in (`checkpoint = ...`),
  couvrir l'état COMPLET (nœuds, éléments, joints, historiques de contact,
  compteurs d'énergie, RNG, outil/corps), avec écriture atomique et reprise
  bit-identique à 1 thread.
- **FR-005**: la sortie VTU binaire DOIT être opt-in (`vtkBinary = true`),
  appended encodé selon le standard VTK, relisible ParaView, et couverte par
  un repère de suite (mêmes champs aux tolérances d'écriture près).
- **FR-006**: chaque livrable de ce chantier DOIT respecter la constitution
  (bit-neutralité des défauts, zeroload, suite, budget d'énergie, docs) —
  les « constitution checks » du plan les listeront un à un.
- **FR-007**: la décision P2 DOIT être documentée dans un rapport court
  (chiffres P1 → issue, coût estimé, risques) et validée par Fernando AVANT
  toute implémentation de l'issue — c'est une décision d'architecture au
  sens de la constitution (« Workflow »).
- **FR-008**: l'issue élue fera l'objet d'une spec dédiée (002-…) avec ses
  propres critères — la présente spec s'arrête à : mesure faite, décision
  prise, prérequis transverses (checkpoint, E/S) livrés.

### Key Entities

- **Charge de référence P1** : configs gelées + graines + machine déclarée ;
  l'étalon de toute mesure du pilier.
- **Rapport P1** : les chiffres (wall, profil, scaling, mémoire) +
  extrapolation au banc réel + incertitudes.
- **Décision P2** : issue (a/b/c), justification chiffrée, validation
  Fernando, date.
- **Point de contrôle** : fichier(s) d'état complet du solveur, versionné
  par un numéro de format, écrit atomiquement.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: le rapport P1 existe, reproductible (configs gelées, machine
  déclarée), et l'extrapolation au banc réel porte une incertitude
  argumentée (scaling mesuré, pas supposé).
- **SC-002**: la décision P2 est prise sur la grille FR-003, validée par
  Fernando, en connaissance des chiffres — zéro implémentation d'issue
  avant validation.
- **SC-003**: un run coupé/repris est bit-identique (1 thread) au run
  d'une traite sur le cas de démonstration ; le checkpoint coûte < 2 % du
  temps mur à la cadence par défaut.
- **SC-004**: E/S ≤ 5 % du temps mur sur la charge de référence en binaire ;
  relecture ParaView vérifiée.
- **SC-005**: (portée par la spec 002 de l'issue élue) le banc complet tient
  le budget temps déclaré en C2 sur le matériel déclaré en C3, suite verte
  sur la cible.

## Assumptions

- La machine de mesure P1 de cette session (2 cœurs cloud) sert au PROFIL et
  aux rapports de coût RELATIFS ; les chiffres ABSOLUS du banc se rapportent
  à la machine de production de Fernando [NEEDS CLARIFICATION C4 : specs de
  cette machine — CPU/cœurs/RAM/GPU éventuel] — le rapport P1 séparera
  explicitement les deux.
- Le banc V3 (piston-taillant-roche) n'étant pas encore assemblé, la charge
  de référence P1 s'appuie sur les cas existants (V1, percussion longue) +
  un cas dimensionné à l'échelle du banc réel ; elle sera re-gelée quand V3
  existera (le plan v2 plaçait P1 « dans V3 » — on le décale en amont pour
  débloquer la décision, l'extrapolation en tient lieu).
- Les données expérimentales (B7) ne sont PAS nécessaires à ce chantier.
- D0 (coût du potentiel en phase débris) est un chantier SÉPARÉ : la charge
  de référence inclut le potentiel tel quel ; si D0 aboutit avant P1, la
  charge est re-mesurée (une ligne de plus, pas un blocage).
