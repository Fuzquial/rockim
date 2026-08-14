# Tasks: Pilier performance — architecture de calcul

**Input**: plan.md (2026-08-14) | **Ordre**: research → US3 → P1 → STOP décision → US2 → docs

## Phase 0 — Research

- [ ] T001 research.md : inventaire d'état EXHAUSTIF de Fdem3dSolver (table
      membre → sauver/reconstructible/transient, revue contre le header)
- [ ] T002 research.md : idem FdemSolver (2D)
- [ ] T003 research.md : format checkpoint (magique/version/endian/rotation
      N=2/hash de compatibilité) + format VTU appended (UInt64, raw)

## Phase A — US3 : VTU binaire (D1)

- [ ] T010 Backend d'écriture appended dans le writer VTU 3D (opt-in
      `vtkBinary`, mêmes champs même ordre) [Fdem3dSolver.cpp]
- [ ] T011 Idem 2D [FdemSolver.cpp]
- [ ] T012 tools/vtu_read.py : lecteur commun ASCII+binaire ; migrer
      crater_metrics.py dessus
- [ ] T013 Repère suite `vtk_binary_roundtrip` (champs identiques aux
      tolérances d'écriture) + relecture ParaView vérifiée
- [ ] T014 Docs (clé vtkBinary) + suite fast + commit

## Phase B — Charge P1 + rapport + décision

- [ ] T020 meshes/p1_banc : éprouvette 12×12×12 cm + insert R11, h≈2 mm
      (~0,5-1 M tets) via make_unstructured_mesh.py bench1
- [ ] T021 configs/p1_{bench1,perc_pen,perc_pot,banc}.cfg gelées (graines,
      T courts mesurés pour le banc)
- [ ] T022 tools/p1_report.py : exécution séquentielle idle, wall +
      ROCKIM_PROF + RSS, 1 et 2 threads, extrapolation 4/8 cœurs
- [ ] T023 Mesures + p1-report.md (chiffres, hypothèses, incertitudes,
      section décision pré-remplie par la grille FR-003)
- [ ] T024 **STOP — validation Fernando de la décision P2** (FR-007)

## Phase C — US2 : checkpoint/restart (D7)

- [ ] T030 include/rockim/Checkpoint.hpp (writer/reader binaire à
      primitives, écriture atomique tmp+rename)
- [ ] T031 writeCheckpoint/readCheckpoint 3D (liste explicite des membres
      de T001) + clés checkpointEvery/resume + hash compat
- [ ] T032 Idem 2D (liste de T002)
- [ ] T033 Repères `checkpoint_resume_2d/3d` : run coupé à T/2 + repris ==
      run d'une traite (bit-identique 1 thread, budget B4 compris)
- [ ] T034 Repère zeroload avec checkpoint actif (0 cassé, gcWork = 0)
- [ ] T035 Coût du checkpoint mesuré (< 2 % wall à cadence défaut)
- [ ] T036 Docs + suite full + commit

## Phase D — Clôture

- [ ] T040 ROADMAP (D1/D7 faits), spec 001 status → Implemented (US1-3),
      bundle ; l'issue P2 élue → /speckit-specify 002
