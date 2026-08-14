# Implementation Plan: Pilier performance — architecture de calcul

**Branch**: `001-pilier-performance-architecture` | **Date**: 2026-08-14 | **Spec**: specs/001-pilier-performance-architecture/spec.md

**Input**: spec.md (Clarified 2026-08-14, C1-C4 résolues)

## Summary

Mesurer le coût réel d'un banc percussif à l'échelle de Mines Paris
(éprouvette 10-15 cm → 0,5-2 M tets) sur la seule machine disponible
(portable 4-8 cœurs), en tirer la décision d'architecture P2 par la grille
de la spec (issue (a) CPU pressentie : ~1,4 h estimées sur 8 cœurs pour
1 M tets — à CONFIRMER par la mesure), et livrer les deux prérequis
transverses valables quelle que soit l'issue : restart/checkpoint à reprise
bit-identique (D7) et sorties VTU binaires (D1). L'issue élue sera spécifiée
à part (002) après validation Fernando du rapport P1.

## Technical Context

**Language/Version**: C++17 (g++/libstdc++ Linux, MSVC Windows), Python 3.11
pour les outils.

**Primary Dependencies**: Eigen 3.4 headers-only, OpenMP optionnel — AUCUNE
dépendance nouvelle dans ce chantier (checkpoint = binaire maison versionné ;
VTU binaire = format « appended raw » du standard VTK, écrit à la main comme
l'ASCII actuel).

**Storage**: fichiers — checkpoints binaires (`<out>/checkpoint_NNN.rkc`,
écriture atomique tmp+rename), VTU appended, configs texte inchangées.

**Testing**: `tools/verify_suite.py` (repères en dur) + contrôles
constitution (zeroload, bit-identité, budget d'énergie).

**Target Platform**: Linux (mesures de cette session) + Windows/MSVC
(machine de production — re-baseline F1 à faire).

**Project Type**: exécutable unique piloté par config (inchangé).

**Performance Goals**: rapport P1 complet ; banc extrapolé ≲ heures sur
8 cœurs sinon dossier matériel ; checkpoint < 2 % du temps mur ; E/S ≤ 5 %.

**Constraints**: constitution I-VII (opt-in bit-identique, zeroload, 2D+3D,
mesures séquentielles, suite, énergie, docs) ; RAM portable (≤ ~16 Go
supposés : 2 M tets × ~1 Ko/tet état FDEM ≈ 2 Go — tenable, à mesurer).

**Scale/Scope**: charge P1 = 4 cas gelés ; checkpoint 2D+3D ; VTU binaire
2D+3D (éléments + joints) ; rapport + décision.

## Constitution Check

*GATE avant Phase 0 — re-vérifié après Phase 1.*

| principe | conformité du plan |
|---|---|
| I. Bit-neutralité (NON NÉG.) | `checkpoint =` et `vtkBinary =` opt-in ; sans clé : trajectoires et fichiers bit-identiques. La reprise checkpoint vise la bit-identité à 1 thread (SC-003) — c'est un CRITÈRE, pas un vœu |
| II. Charge nulle | repère zeroload avec checkpoint activé (écrire/reprendre au repos : 0 cassé, gcWork = 0) |
| III. 2D+3D de front | checkpoint et VTU binaire dans FdemSolver ET Fdem3dSolver, même chantier ; la charge P1 est 3D (le banc l'est) — documenté |
| IV. La mesure fait foi | P1 séquentiel machine idle, 1 puis 2 threads ici, extrapolation déclarée par machine ; critères chiffrés posés dans la spec AVANT |
| V. Suite = contrat | nouveaux repères : `checkpoint_resume_2d/3d` (bit-identité), `vtk_binary_roundtrip` (relecture) ; réfs inchangées ailleurs |
| VI. Énergie = juge | les accumulateurs B4 font partie de l'état checkpoint (un budget faux après reprise = échec du repère) |
| VII. Docs | DOCUMENTATION (clés checkpoint/vtkBinary/resume), ROADMAP (D1/D7), même commit |

**Verdict** : PASS — aucune dérogation demandée.

## Project Structure

### Documentation (this feature)

```text
specs/001-pilier-performance-architecture/
├── spec.md          # fait (clarifié)
├── plan.md          # ce fichier
├── research.md      # Phase 0 : inventaire d'état, format VTK appended
├── p1-report.md     # le RAPPORT P1 (mesures + extrapolation + décision)
└── tasks.md         # /speckit-tasks
```

### Source Code (repository root)

```text
include/rockim/Checkpoint.hpp     # NOUVEAU : sérialisation binaire versionnée
src/FdemSolver.cpp                # +writeCheckpoint/readCheckpoint, +VTU binaire
src/Fdem3dSolver.cpp              # idem 3D
src/main.cpp (ou équiv.)          # clé resume = <fichier>
tools/verify_suite.py             # +3 repères
tools/p1_report.py                # NOUVEAU : lance la charge gelée, produit p1-report.md
configs/p1_*.cfg                  # NOUVEAU : la charge de référence GELÉE (4 cas)
meshes/p1_banc_h20.msh            # NOUVEAU : maillage banc extrapolé (~0,7 M tets, h≈2 mm)
```

## Phase 0 — Research (research.md)

1. **Inventaire d'état exhaustif** des deux solveurs (le piège n° 1 du
   checkpoint est l'état oublié) : nœuds (X0_, u_, v_, m_, flag_, kAbs/cAbs,
   lastTouch_), éléments (dN, V0/A0, phase/grain, MatState st, sigG/svm),
   joints (TOUT : D, slip, smax, omax, dn0, bonded, dead, tBreak, bmode…),
   contact (pen0_, potFt_ avec Ft/step/vRef/sepAxis, jointOfPair_, pool_,
   extOn_, deadList_, actStamp_, sweep…), outil/corps (tool_, elemGroup_,
   trackGroup_), scalaires (t_, stepCount_, nBroken_, nInserted_, TOUS les
   accumulateurs B4, gcWork_, dampWork_, peakF_, work_, keInit_…), RNG
   (std::mt19937 sérialisable par operator<<), groupes adaptatifs
   (grpsOfVert_…). Méthode : lecture systématique des headers, table
   membre → {sauver | reconstructible | transient} avec justification.
2. **Format checkpoint** : magique + version + mode + tailles, little-endian
   déclaré, doubles bruts (la bit-identité interdit toute conversion texte).
   Un fichier par point de contrôle, rotation N=2 (le précédent reste valide
   pendant l'écriture du suivant — edge case de la spec).
3. **Format VTU appended** : `format="appended" offset=…` + bloc
   `<AppendedData encoding="raw">_`, en-têtes de taille UInt64,
   header_type="UInt64" — vérifier la relecture ParaView ET par nos scripts
   regex (qui devront apprendre le binaire : `tools/vtu_read.py` partagé).
4. **Dimensionnement du cas banc** : éprouvette 12×12×12 cm, insert R 11 mm,
   h ≈ 2 mm → générer et compter ; viser 0,5-1 M tets ; T mesuré court
   (2-5 e-5 s) + extrapolation en pas — les hypothèses écrites dans le
   rapport.

## Phase 1 — Design

- **Checkpoint** : `Checkpoint.hpp` fournit un writer/reader binaire à
  primitives (`put(T)/get(T)`, vecteurs, maps) ; chaque solveur implémente
  `writeCheckpoint(path)` / `readCheckpoint(path)` en LISTANT explicitement
  chaque membre (pas de sérialisation « magique » — la liste EST la
  documentation de l'état). Clés : `checkpointEvery = <n_pas|0>` (0 = off,
  défaut), `resume = <fichier>` (charge l'état APRÈS init du maillage depuis
  la config — le maillage se reconstruit, l'état écrase ; garde : hash de
  compatibilité config/maillage dans l'en-tête). Reprise : le run témoin et
  le run coupé/repris passent par les MÊMES appels dans le même ordre —
  aucune renumérotation.
- **VTU binaire** : `vtkBinary = true|false` (défaut false) ; les writers
  actuels gagnent un backend (ASCII|appended) — mêmes champs, même ordre.
  `tools/vtu_read.py` : lecteur commun ASCII+binaire, adopté par
  crater_metrics et les scripts de figures (une seule implémentation).
- **Charge P1 gelée** (`configs/p1_*.cfg`, graines fixes) :
  1. `p1_bench1.cfg` — l'impact insert V1 tel quel (17,9 k tets, T 2e-4) ;
  2. `p1_perc_pen.cfg` — percussion 3D longue pénalité+adaptatif (19,3 k) ;
  3. `p1_perc_pot.cfg` — idem potentiel (l'état D0 assumé tel quel) ;
  4. `p1_banc.cfg` — le cas banc ~0,7 M tets, T court mesuré.
  `tools/p1_report.py` : exécute (séquentiel, idle), relève wall +
  ROCKIM_PROF + RSS crête, répète à 2 threads, écrit `p1-report.md` avec
  l'extrapolation (pas × tets, scaling observé) pour 4 et 8 cœurs.
- **Décision P2** : section finale de p1-report.md pré-remplie par la grille
  FR-003 clarifiée — validation Fernando AVANT toute spec 002.

## Ordre d'exécution (proposé pour /speckit-tasks)

1. Phase 0 complète (research.md — l'inventaire d'état d'abord : il
   conditionne tout le checkpoint).
2. US3 VTU binaire (court, indépendant, débloque les gros cas de P1 dont
   l'E/S ASCII fausserait la mesure).
3. Charge P1 + maillage banc + p1_report.py → mesures → **p1-report.md**
   → STOP décision Fernando (FR-007).
4. US2 checkpoint (le gros morceau : inventaire → impl 2D+3D → repères
   bit-identité → zeroload).
5. Docs + repères + commit par étape (constitution VII).

## Complexity Tracking

Aucune dérogation constitution demandée. Risques notés : (a) l'inventaire
d'état du checkpoint est le point de fragilité — mitigé par la liste
explicite revue membre à membre contre les headers ; (b) la RAM du cas banc
sur cette machine cloud (~2 Go attendus) — mesurée avant de lancer ;
(c) MSVC non testable d'ici (F1) — le format binaire est little-endian
explicite et sans padding implicite (écriture champ par champ).
