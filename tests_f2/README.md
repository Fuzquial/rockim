# Bancs courts et decks de test du chantier f2 et des lots suivants

Ce dossier rassemble des tests écrits au fil des chantiers de septembre et octobre 2026. Chaque test
vérifie une capacité précise avec un critère chiffré, souvent accompagné d'une variante qui doit
échouer pour prouver que le critère discrimine. Ces tests sont courts (de la seconde à quelques
minutes) mais, à l'exception des decks de `bitid/`, ils ne sont pas branchés dans
`tools/verify_suite.py` : ils se rejouent à la main et leurs verdicts sont consignés dans les fichiers
`RESULTATS_*` du dossier, dans `CHANGELOG.md` ou dans les notes de `docs/notes/`.

Les chemins absolus Windows qui apparaissent dans certains scripts et résultats sont ceux du poste du
doctorant ; ils sont conservés tels quels.

## Vue d'ensemble

| Élément | Nature | Statut |
|---|---|---|
| `bitid/` | 9 decks courts de bit-identité (fem3d, fdem3d, fdem) | en service : joués par `tools/bitid.py` |
| `t1_*.cfg` à `t18_*.cfg`, `tg_*.cfg`, `dg_*.cfg` | 47 decks de traction directe 2D du chantier f2, une capacité par deck | validés à la main en septembre 2026 ; pas de runner |
| `check_jointtsl_header.cpp`, `check_camacho_friction.cpp`, `check_camacho_friction3d.cpp`, `bench_lotb_joints.cpp` | bancs C++ des noyaux de la loi de joint `camacho` | passés le 2026-09-11 (15/15 et 13/13 selon `CHANGELOG.md`) |
| `bench_energie_matrice.cpp`, `build_bench_energie.cmd` | banc d'énergie au point matériel de la loi de volume de la note 2026 | joué en septembre 2026 ; résultat dans `docs/notes/RETOUR_v3_2026-09-11.md` §1.4 sexies |
| `insertion_pointe/` | banc de la clé `insertionTipFactor` et de `insertion = none` | tous critères passés (`docs/notes/PORT_INSERTION_POINTE.md` §4.3) |
| `psivar/` | banc de la dilatance variable ψ(p) de DP-DFH | tous critères passés (`docs/notes/PORT_INSERTION_POINTE.md` §4.2) |
| `rang3/` | micro-banc de stabilité du basculement de raideur joint et contact | fait le 2026-10-03, verdict partiel |
| `campagne13/` | mini-tests de la campagne de correction du 13/09/2026 | tous les bilans OK |
| `_eapp.py` | dépouillement des modules apparents de `t14b_ti_*` | outil du chantier f2 |

## Decks de bit-identité, `bitid/`

Huit decks forment l'ancre de `tools/bitid_refs.json` : `fem3d_cdp_PQ_court`, `fem3d_dpr_T1_court`,
`fem3d_sk2011_cyl_court`, `fdem3d_kuru9_court`, `fdem3d_cut3d_heilman_court`, `fdem3d_visc_yan_3d`,
`fdem_toolcontact_signorini`, `fdem_ucs_yan_adaptive_court`. Le neuvième, `fdem3d_yang_v2_court`, a sa
propre ancre (`tools/bitid_refs_jointlaw.json`). Mode d'emploi : `tools/BITID.md`.

```bash
python3 tools/bitid.py --exe build/rockim
```

## Decks du chantier f2 (`t*.cfg`, `tg_*.cfg`, `dg_*.cfg`)

Tous partent du deck de vérification `configs/verify_fdem_voronoi_tension.cfg` (traction directe sur
une structure de grains de Voronoï monophasée) et n'y ajoutent que quelques clés. Les deux decks
`dg_*_grid` partent de `configs/verify_fdem_tension.cfg` (bande uniforme). Le tableau donne les clés
ajoutées. Description des capacités : `DOCUMENTATION_rockim.md` §5.10 à §5.15.

| Deck | Clés ajoutées | Ce qui est vérifié |
|---|---|---|
| `t1_seismic` | `microseismic` | catalogue microsismique des ruptures |
| `t2_prebroken`, `t5_fracpath` | `preBrokenJoints` ou `jointPrebrokenFrac`, `jointResidualMu` | joints pré-rompus, frottement résiduel |
| `t3_garde_vide`, `t4_garde_mu` | segment pré-rompu hors domaine ; frottement résiduel absent | refus attendu à l'initialisation |
| `t6_cavclose_off`, `t6_cavclose_on` | `hydro`, `hydroCavityClosure` | fermeture de la cavité hydraulique |
| `dg_press_grid`, `dg_rate_grid`, `dg_rate_voronoi` | `hydro`, injection en pression ou en débit | hydro sur grille et sur Voronoï |
| `t7_damping`, `t7_damping_off`, `t8_damping_impl` | `dampingViscous`, `dampingViscousScheme` | amortissement visqueux explicite et implicite |
| `t8_garde_orpheline` | `dampingViscousScheme` sans amortissement | refus attendu |
| `t9_thermo_sigma`, `t9b_thermo_adaptive`, `t10_thermo_cracks` | `thermal` | choc thermique de paroi, avec et sans fissuration |
| `t11_schisto*`, `t12_aniso*`, `t18_wp_bimodal` | `weakPlanes`, `weakPlaneFactor` | famille de plans de faiblesse, balayage de l'angle de 0 à 90° |
| `t13*`, `t14*` | `beddingDip`, `beddingEperp`, `beddingGperp` | élasticité isotrope transverse de Lisjak contre la forme fermée |
| `t15_gamma_neutre`, `t16_gamma_0`, `t16_gamma_90` | `beddingFtRatio`, `beddingCohRatio`, `beddingGfIRatio` | loi cohésive directionnelle ; neutre si les rapports valent 1 |
| `t17_jointstate` | `writeJointState` | état de contact des joints en sortie |
| `tg_3d`, `tg_law`, `tg_orphan`, `tg_ti_thermal` | `thermal` en 3D, avec `law`, sans `thermal`, avec litage | refus attendus (combinaisons non prises en charge) |

Résultats mesurés (extrait de `docs/notes/CHANTIER_f2.md`) :

| Essai | Attendu | Mesure |
|---|---|---|
| `t15_gamma_neutre` (rapports = 1) | bit-identique au deck sans clé | 28 fichiers, 0 différent |
| `t13c_nocund_*`, isotrope transverse à la limite isotrope, sans Cundall | écart nul | écart exactement nul sur 2 080 lignes |
| `t14c_ti_*`, modules apparents | 1,5910 et 4,3307 GPa | 1,5910 (+0,00 %) et 4,3364 (+0,13 %) |
| `t16_gamma_0` et `t16_gamma_90` | anisotropie de Lisjak, rapport voisin de 1,9 | pic de 4,35 à 8,43 MPa, rapport 1,94 |
| gardes `tg_*` | refus | refus propres |

Le plan de robustesse du 2026-09-05 (`docs/notes/PLAN_ROBUSTESSE_2026-09-05.md`, mesure B10) prévoyait
un tier `f2` dans `verify_suite.py` pour ces decks. Il n'est pas réalisé.

Lancer un deck :

```bash
./build/rockim tests_f2/t16_gamma_0.cfg out_t16_gamma_0
```

## Bancs C++ de la loi de joint et de la loi de volume

Ces fichiers se compilent seuls contre les en-têtes de `include/rockim/` et n'appellent que les
fonctions pures que le solveur appelle : ils vérifient le code livré, pas une copie. La commande de
compilation MSVC est en tête de chaque fichier.

| Fichier | Vérifie |
|---|---|
| `check_jointtsl_header.cpp` | `JointTsl.hpp` avec branche ascendante : continuité de la traction à l'insertion, raideur de charge bornée, intégrale de l'adoucissement égale à G_C |
| `check_camacho_friction.cpp` | frottement en cap du joint `camacho` en 2D (`jfric::capReturn`) : nul en traction, dissipation par cycle positive, continuité en δ_s = 0, variante discontinue qui doit échouer |
| `check_camacho_friction3d.cpp` | même banc en 3D (`camachoFrictionSlider`) |
| `bench_lotb_joints.cpp` | lot B : fonctions de joint et DIF face aux exigences de la note |
| `bench_energie_matrice.cpp` | la loi de volume de la note 2026 (`law = saksala` et options) crée-t-elle de l'énergie sur un cycle fermé ? Aucune création sur 616 cellules ; création identifiée sur un trajet pré-broyé avec rotation d'axes, avec sa formule exacte |

## `insertion_pointe/` et `psivar/`

Bancs des capacités portées de la branche insertion-pointe.

```bash
python3 tests_f2/insertion_pointe/check_tip.py --exe build/rockim --threads 2
python3 tests_f2/psivar/check_psivar.py --exe build/rockim
```

`check_tip.py` lance quatre compressions uniaxiales 2D courtes. Critères : clé absente et
`insertionTipFactor = 1` bit-identiques ; facteur 1,6 qui change la trace ; pic de résistance
inchangé (51,0807 MPa dans les trois cas). `check_psivar.py` pilote deux points matériels DP-DFH :
la clé `dfhPsiVar` agit, ψ décroît avec la pression moyenne, et `--falsify` fournit la variante qui
doit échouer.

## `rang3/` : basculement de raideur joint et contact

Deux tétraèdres liés par un joint ; A encastré, B lancé. `balayage.py` croise
`jointContactPenalty` (fixed, adaptive), `dtFactor` (0,05 à 0,9) et la vitesse initiale (12 et
13 m/s). Résultats : `RESULTATS_balayage_2026-10-03.csv`. Sur environ 1 600 oscillations,
l'amplitude ne croît dans aucun cas et le bilan d'énergie ferme à 10⁻¹² %. Réserve : les
excursions en compression sont faibles, donc le basculement est peu sollicité. Ce n'est pas une
preuve que le réglage par défaut est sûr (`docs/notes/HANDOFF_2026-10-03.md` §2.3).

```bash
python3 tests_f2/rang3/balayage.py build/rockim out_rang3
```

## `campagne13/` : campagne de correction du 13/09/2026

Mini-tests des tâches de la campagne décrite dans `docs/notes/CAMPAGNE_correction_2026-09-13.md`.
Chaque sous-dossier contient ses decks, un script `check_*.py` ou `verif_*.py` et ses résultats.

| Dossier | Tâche | Bilan consigné |
|---|---|---|
| `S1/` | instrumentation de rupture : champs `dead`, `openMax`, `pMean`, clé `jointBreakModeRef` | sorties vérifiées, bit-identité IDENTIQUE |
| `S2/` | forces de contact entre corps nommés (`contactForcePairs`) | `BILAN S2 : TOUS OK`, `BILAN S2bis : TOUS OK` |
| `S3/`, `S3bis/` | banc de joint cinématique `scenario = jointbench` : décharge, énergie de rupture, règle deux points sur trois | « schéma falsifiant entièrement tenu » ; les FAIL listés sont ceux des variantes qui doivent échouer |
| `S4/` | pulvérisation : mesures de longueur et de déformation, essais élémentaires | `BILAN : 0 verdict(s) FAIL` |
| `T1/` | vérification de `tools/crack_paths.py` | tableaux `verif_*.csv` |
| `T3/` | masses par corps d'une série de maillages | `smoke_masses.md` |
| `T4/` | fumée des decks de conformité (St Anne, série Kuru) | journaux `smoke_*.md` |
| `B/` | rejeu de S1 à S4 avec le binaire final `rockim_g1y17.exe` | `RESULTATS_check_s*_g1y17.txt` |
