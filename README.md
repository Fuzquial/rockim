# rockim

Solveur FDEM (méthode des éléments finis et discrets combinés, lignée Munjiza) écrit pour la thèse
« forage percussif dans le granite de Red Bohus » (F. Uzquiano, Mines Paris – PSL). Il simule
l'impact d'un insert de carbure sur la roche : élasticité de volume, joints cohésifs insérés entre
éléments, contact entre fragments, frottement, effets de vitesse, bilan d'énergie.

Un encadrant qui découvre le dépôt commence par **[LIRE_EN_PREMIER.md](LIRE_EN_PREMIER.md)**.

## Ce que contient le code

Un exécutable unique, `rockim`, lit un fichier de configuration texte et écrit un dossier de
résultats. Six solveurs partagent la même interface (clé `mode`) :

| mode | modèle | fissuration | usage |
|---|---|---|---|
| `fdem` | FDEM 2D (Munjiza), grains de Voronoï | joints cohésifs, insertion intrinsèque ou adaptative | essais de laboratoire (brésilien, compression, triaxial, barres de Hopkinson), tunnel, coupe, fracturation hydraulique |
| `fdem3d` | FDEM 3D, tétraèdres | joints cohésifs triangulaires | impact d'un insert, coupe 3D |
| `fem3d` | éléments finis 3D, cinq lois de volume | endommagement et érosion | percussion 3D, comparaison avec Abaqus |
| `fem` | éléments finis 2D | endommagement et érosion | essais rapides, onde de barre |
| `dem`, `dem3d` | particules liées (BPM) | rupture des liaisons | comparaison historique |

Unités SI partout ; point décimal obligatoire dans les fichiers de configuration.

## Compiler et vérifier

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
python3 tools/verify_suite.py --exe build/rockim --tier fast
```

Eigen est récupéré automatiquement s'il manque. OpenMP est optionnel. La suite rapide de
non-régression passe 65 contrôles sur 66 sous Linux ; le seul échec, `t1_toolcontact_penalty`, est
connu et documenté. Ses valeurs de référence sont enregistrées en série (`OMP_NUM_THREADS = 1`).
Sous Windows, les scripts de compilation historiques sont dans `scripts/windows/`.

## Documentation

| document | contenu |
|---|---|
| [LIRE_EN_PREMIER.md](LIRE_EN_PREMIER.md) | parcours de lecture pour un encadrant, résultats et limites clés |
| [docs/rapport_guide/rapport_guide_rockim.pdf](docs/rapport_guide/rapport_guide_rockim.pdf) | rapport-guide : formulation, vérification, comparaisons, limites (123 pages) |
| [GUIDE_rockim.md](GUIDE_rockim.md) | guide pratique : lancer une simulation pas à pas |
| [DOCUMENTATION_rockim.md](DOCUMENTATION_rockim.md) | référence complète des commandes, des clés et des sorties |
| [CHANGELOG.md](CHANGELOG.md) | historique des changements |
| [docs/VV_campagne.md](docs/VV_campagne.md) | campagne de vérification et validation |
| [LISEZ_MOI.md](LISEZ_MOI.md) | instructions du dossier de partage (installation sous macOS) |

## Carte du dépôt

| dossier | contenu |
|---|---|
| `src/`, `include/` | code source |
| `configs/`, `configs_bench/`, `configs_yan/` | fichiers de configuration |
| `meshes/` | maillages |
| `tools/` | suite de vérification, post-traitement, figures |
| `tests_f2/`, `vv/` | campagnes de tests et bancs de vérification |
| `docs/rapport_guide/` | rapport-guide, chapitres, figures, simulations préparées, notes de la correction du contact |
| `docs/notes/` | notes de travail datées (bilans, plans, transmissions), conservées pour l'historique |
| `bench_*`, `tunnel_edz/`, `tunnel_schisto/`, `calib*`, `gbm_uniaxial/`, `etude_lois_fem/`, `exo_tunnel/` | études et bancs (README dans `bench_*`, `tunnel_edz/`, `calib*` ; pas encore dans les autres) |
| `scripts/windows/`, `scripts/historique/` | scripts de compilation Windows et scripts de lancement historiques (chemins du poste du doctorant) |
| `specs/` | spécifications des développements |

## Règle de développement

Toute capacité nouvelle est une option désactivée par défaut : sans elle, les sorties restent
identiques au bit près. Les corrections du contact du 7 octobre 2026 suivent cette règle
(`contactCandidates = vertex`, `gcBirth = offset`, `potForceExact = true`) ; voir
`docs/rapport_guide/correction_contact/`.
