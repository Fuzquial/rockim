# Outils de rockim

Scripts Python (3.9 ou plus, numpy et matplotlib ; meshio ou gmsh pour certains mailleurs), scripts
shell et petits pilotes C++. Ils se lancent depuis la racine du dépôt. Aucun ne modifie `src/`.

## Commande type de la suite de non-régression

```bash
python3 tools/verify_suite.py --exe build/rockim                     # tier fast, environ 1 min
python3 tools/verify_suite.py --exe build/rockim --tier full         # + repères longs, adaptatif, 3D grille
python3 tools/verify_suite.py --exe build/rockim --tier all          # + Voronoï 3D, environ 40 min de plus
python3 tools/verify_suite.py --exe build/rockim --only fdem3d       # filtre sur le nom
python3 tools/verify_suite.py --exe build/rockim --refs refs_macos_arm.json   # références de la plateforme
```

La suite compte 135 tests : 66 au tier `fast`, 58 au tier `full`, 11 au tier `all` (chaque tier inclut
les précédents). Elle tourne à un fil par défaut. Les invariants (charge nulle sans joint rompu,
amortissement qui ne crée pas d'énergie, compteurs, selftests) sont écrits en dur dans le script ; les
références chiffrées dépendent de la plateforme et se rebasent avec `--update-refs`. Le code de retour
vaut 0 quand tout passe.

La bit-identité d'un binaire contre une ancre se vérifie avec `bitid.py` (huit decks courts de
`tests_f2/bitid/`) :

```bash
python3 tools/bitid.py --exe build/rockim --json results/bitid.json
```

Méthode et règles : `tools/BITID.md`.

## Inventaire par famille

### Vérification et non-régression

| Script | Usage |
|---|---|
| `verify_suite.py` | suite de non-régression à trois tiers (`fast`, `full`, `all`), tolérances par nature de test, rapport `--json` |
| `bitid.py` | bit-identité de `history.csv`, `frames.csv` et des VTU de la dernière trame contre une ancre (`bitid_refs*.json`) ; `--update` reprend l'ancre |
| `bitid_refs.json`, `bitid_refs_jointlaw.json`, `bitid_refs_w18_herite.json` | ancres de bit-identité (hachages SHA-256 par deck et nombre de fils) |
| `BITID.md` | méthode de la bit-identité, options et codes de sortie de `bitid.py` |
| `ab_history.py` | deux `history.csv` du même deck (deux binaires ou deux nombres de fils) : bit-identiques sur les instants communs ? |
| `deck_smoke.py` | lancement de fumée de decks (T et trames réduits) : erreurs de clé, dt, masses, coût estimé du run complet |
| `t4_check_instr.py` | vérifie l'instrumentation des decks de conformité à partir d'un dossier de fumée de `deck_smoke.py` |
| `check_provenance.py` | vérifie que chaque extrait de code cité dans un `.tex` (`\provenance{src/Fichier.cpp:ligne}`) correspond au dépôt |
| `test_yang_estimators.py` | test sur signaux synthétiques de `yang_estimators.py` |
| `facet_cycle_camacho.py` | portage Python, terme à terme, de la branche `jointTSL = camacho` au point d'intégration : travail du joint sur un cycle |
| `yan_point.cpp` | pilote au point matériel de l'adoucissement de Yan (lié à `YanSoftening.hpp`) |
| `potvolume_partiel.cpp` | force de contact `potForce = volume` contre le potentiel de Munjiza sur deux tétraèdres en recouvrement partiel |

### Clés de configuration et binaires

| Script | Usage |
|---|---|
| `gen_keys_by_mode.py` | régénère `keys_by_mode.json` et `include/rockim/KeysByMode.hpp` à partir des lectures de clés du code |
| `keys_by_mode.json`, `obsolete_keys.json` | registre des clés par mode et table des clés renommées |
| `scan_decks.py` | balayage statique de tous les decks `*.cfg` selon la règle « clé inconnue = erreur » ; `--fix-obsolete` corrige les noms anciens |
| `scan_decks_2026-09-05.md` | rapport du balayage du 2026-09-05 (881 decks) |
| `exe_manifest.py`, `exe_manifest.json` | empreinte SHA-256 des binaires historiques `rockim_f2*.exe` et leur ligne dans le journal |
| `build.ps1` | compilation Windows par CMake et Ninja (remplace `scripts/windows/build_*.cmd`) |

### Dépouillement et post-traitement

| Script | Usage |
|---|---|
| `plot_results.py` | aperçu rapide d'un dossier de sortie : carte d'endommagement ou particules, historique de force |
| `plot_force_penetration.py` | courbes force-pénétration d'un outil disque, plusieurs runs |
| `impact2d_report.py` | impact 2D : restitution, enfoncement, joints rompus, fragments, énergie |
| `yang_report.py` | impact fdem3d de type Yang 2026 : critères et cinétique lus dans `history.csv` |
| `yang_estimators.py` | estimateurs de Yang 2025 (vitesses d'indentation et de rebond, contrainte de référence) |
| `crater_metrics.py` | métriques de cratère d'un impact fdem3d |
| `crack_paths.py` | fissures connectées (composantes de facettes rompues), longueur radiale et rayon de cratère |
| `joint_state_from_vtu.py` | reconstitue ouverture et glissement par facette à partir des VTU d'un run fdem3d |
| `law_matrix.py` | matrice fracture × conservation : une ligne par run à un instant donné |
| `bench_compare.py` | plusieurs runs d'impact côte à côte : physique à un instant commun et coût par pas |
| `export_abaqus.py` | exporte le maillage d'un run (trame 0) en `.inp` Abaqus pour une comparaison sur maillage identique |
| `morning_s1.sh` | dépouillement complet d'un run d'impact en une commande (coupes, joints, critères) |

### Figures

| Script | Usage |
|---|---|
| `fig_impact3d.py` | planche de courbes et GIF de coupe d'un impact 3D, comparés à Yang 2026 |
| `fig_kinetics.py` | cinétique de l'impact (contrainte de jauge, vitesse du bit, enfoncement) contre Yang 2026 |
| `fig_fp.py`, `fig_fp_direct.py` | force-pénétration d'un impact 3D : estimateur par quantité de mouvement, puis réaction mesurée (`contactForcePairs`) |
| `fig_retournement.py` | date le point de retour d'un impact non terminé |
| `fig_evolution_cratere.py` | mesures du cratère contre le temps |
| `fig_joints_cuts.py`, `fig_joints_only.py` | coupes exactes du réseau de joints rompus, et joints rompus seuls |
| `fig_stress_section.py` | coupe de contrainte avec le réseau de joints superposé |
| `fig_surface_compare.py`, `fig_icl_stanne_review.py`, `animate_detached_rock.py` | figures de relecture de l'impact St Anne : traces de surface, surfaces de rupture, fragments détachés |
| `fig_loads.py` | figures des charges par groupes (maillage, banc élastique, rupture) |
| `fig_mesh3d.py`, `fig_montage_impact.py` | montrer un maillage tétraédrique ou un montage d'impact avant de lancer |
| `fig_cible_dumoulin.py` | vecteur de cibles de calibration Red Bohus (Dumoulin et al. 2024) |
| `make_gif.py`, `make_gif_cut.py` | GIF animé d'un run (outil compris), et d'un essai de coupe |

### Calibration

| Script | Usage |
|---|---|
| `calibrate_bohus.py` | calibration bayésienne du modèle à grains sur les cibles du granite Red Bohus |
| `bayes_bench.py` | banc synthétique de la chaîne de calibration par émulateur, à vérité connue |

### Maillage

| Script | Usage |
|---|---|
| `make_unstructured_mesh.py` | maillage simplexe non structuré uniforme (boîte 2D ou 3D) pour `mesh = file` |
| `make_impact_mesh.py` | montage d'impact à insert unique (piston, bit, insert, roche) d'après Yang 2025-2026 |
| `make_impact2d_mesh.py`, `make_impact3d_mesh.py` | éprouvette d'impact 2D ou 3D à maillage gradué sous le point d'impact |
| `make_cut_mesh.py`, `make_cut3d_mesh.py` | éprouvette de coupe 2D ou 3D avec entaille de départ (Heilman et al. 2024) |
| `make_conformity_decks.py` | génère les decks de conformité de la campagne du 13/09, sans les lancer |
| `mesh_quality.py` | qualité d'un `.msh` : diamètre inscrit, slivers, quantiles par corps |

### Performance et files de calcul

| Script | Usage |
|---|---|
| `bench_threads.sh`, `prof_threads.sh` | gain en temps selon le nombre de fils, et phase du pas qui passe mal à l'échelle (`ROCKIM_PROF=1`) |
| `queue_benches_s25.sh`, `queue_bench_D.sh`, `queue_bench_E.sh`, `queue_bench_F.sh` | files de bancs d'impact du 13/09, un gros calcul à la fois |
| `queue_nuit_1314.sh`, `queue_nuit_1314b.sh`, `queue_nuit_1314c.sh` | files de la nuit du 13 au 14/09 (St Anne puis Kuru) |

### Interface

| Script | Usage |
|---|---|
| `rockim_gui.py` | interface graphique de pilotage : choix du deck, lancement, suivi |

Les scripts de file (`queue_*.sh`) et certains scripts de dépouillement contiennent des chemins et des
noms de binaires propres au poste Windows du doctorant. Ils sont conservés pour la traçabilité des
campagnes de septembre 2026.
