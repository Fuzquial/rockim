# rockim — solveur FDEM pour le forage percussif
**Dossier de partage constitué le 2026-08-30. F. Uzquiano, Mines Paris – PSL.**

---

## Ce que c'est

**rockim** est un solveur **FDEM** (*combined finite-discrete element method*,
lignée Munjiza) écrit pour la thèse « forage percussif en granite de Red Bohus ».
Il traite l'impact d'un insert sur la roche : élasticité de volume, **joints
cohésifs** insérés entre éléments, **contact général** par potentiel de Munjiza,
frottement, effets de vitesse (DIF), pulvérisation, et un **bilan d'énergie
fermé** à sept postes.

**Ce dossier est un mini-dépôt qui fonctionne tel quel** — il compile et il
vérifie. Ce n'est pas un extrait mort.

## Comment le prendre en main, en trois commandes

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
./build/rockim configs/verify_fdem3d_tension.cfg
python3 tools/verify_suite.py --exe build/rockim --tier fast
```

La première ne demande **aucune dépendance à installer** : Eigen est récupéré
automatiquement s'il n'est pas déjà là (en-têtes seuls, pas de compilation —
il faut juste un accès réseau à la première configuration, ou
`brew install eigen` / `apt install libeigen3-dev` pour l'éviter).
La troisième rejoue **44 contrôles de non-régression** et doit imprimer
`[suite] TOUT PASSE (44/44)`.

> ✅ **Vérifié avant l'envoi, le 2026-08-30** : ce dossier, extrait dans un
> répertoire isolé, a compilé sans erreur et sa suite a rendu
> `[suite] TOUT PASSE (44/44)` — Linux x86-64, g++, `OMP_NUM_THREADS = 1`.
> Ce que vous recevez est ce qui a été testé, pas une copie de ce qui l'a été.

### Sur macOS, précisément

```sh
xcode-select --install     # compilateur Apple Clang + make + python3 (une fois)
brew install cmake         # ou l'installeur de cmake.org si pas de Homebrew
tar -xzf bundle_rockim_2026-08-30.tar.gz && cd bundle_rockim_2026-08-30
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/rockim configs/verify_fdem3d_tension.cfg
python3 tools/verify_suite.py --exe build/rockim --tier fast
```

Trois choses à savoir, toutes vérifiées dans le `CMakeLists.txt` :

* **OpenMP est OPTIONNEL** (`find_package(OpenMP QUIET)`). Apple Clang ne le
  fournit pas : le binaire compile alors **en série, sans erreur** — et c'est
  très bien ainsi, car **la suite de vérification impose de toute façon
  `OMP_NUM_THREADS = 1`** : toutes ses valeurs de référence sont enregistrées en
  série. Un build série sur Mac est donc le build *de référence*, pas un build
  dégradé. (Pour des runs longs multi-fils : `brew install libomp` puis passer
  les drapeaux `OpenMP_*` à cmake, ou `brew install gcc` et
  `-DCMAKE_CXX_COMPILER=g++-14` — utile seulement au-delà de la vérification.)
* **`verify_suite.py` n'utilise que la bibliothèque standard** de Python — le
  `python3` des Command Line Tools suffit, rien à installer.
* **Apple Silicon (M1-M4)** : les références chiffrées de la suite sont la
  baseline **Linux/g++/libstdc++** — sur une autre plateforme, une partie
  d'entre elles échoue, et pas toujours *de justesse*. **Mesuré le 2026-08-30
  sur macOS ARM (Apple Clang 21/libc++) : 10 échecs sur 44 au tier `fast`**,
  bit-identiques entre Eigen 3.4.0 et 5.0.1 (ce n'est donc pas Eigen) — dont
  les cas Voronoï, où `std::shuffle` produit un AUTRE maillage à graine égale
  (7,07 % mesuré contre 11,60 ± 0,05 attendu), et des dérives amplifiées par
  les phases de rupture (wp6_pulv : +11 %). **Les invariants, eux, passent
  tous partout** : charges nulles → 0 joint cassé, `dampWork <= 0`, selftests,
  compteurs entiers. La règle : certifier une plateforme non-Linux contre SA
  baseline — `--update-refs refs_<plateforme>.json` sur un run sain, puis
  `--refs refs_<plateforme>.json` à chaque certification. Un échec de
  *valeur* hors baseline se signale (avec les valeurs imprimées), ne se
  « corrige » pas ; un échec d'*invariant* est un bug, sur toute plateforme.

## Dans quel ordre lire

| # | fichier | pour qui |
|---|---|---|
| 1 | [`etat_de_l_art/LISEZ_MOI.md`](etat_de_l_art/LISEZ_MOI.md) | **commencer ici si vous découvrez le sujet.** Ce que fait le FDEM d'Imperial College, reconstitué sur sources primaires, et où rockim en diverge |
| 2 | [`DOCUMENTATION_rockim.md`](DOCUMENTATION_rockim.md) | **le guide de référence.** Toutes les clés de configuration, leur défaut, ce qu'elles changent, et la mesure qui le prouve |
| 3 | [`SOURCES_SOLIDITY.md`](SOURCES_SOLIDITY.md) | **la provenance.** Les citations du code d'Imperial présentes dans les sources, leur statut, et ce qui reste à vérifier |
| 4 | [`README.md`](README.md) | l'historique du projet et ses partis pris |
| 5 | `include/rockim/*.hpp` | **les en-têtes portent le raisonnement**, pas seulement les déclarations. C'est là qu'est la physique |

> **Le code fait 23 000 lignes ; la documentation en fait 26 800.** Ce n'est pas
> un accident : dans ce dépôt, une capacité non documentée est considérée comme
> non livrée.

## Trois règles du dépôt, qui expliquent sa forme

1. **Toute capacité nouvelle est opt-in, et le défaut reste bit-identique.** On
   n'améliore jamais un résultat existant par surprise.
2. **Une capacité active et muette est indiscernable d'une capacité inerte.**
   Chaque clé armée s'annonce au journal ; plusieurs refusent de démarrer plutôt
   que d'être silencieusement ignorées.
3. **Rien n'entre sans un contrôle de non-régression chiffré.** La suite en porte
   **104**, dont 44 au tier rapide, pour **237 assertions**.

## ⚠️ Ce que ce dossier NE contient PAS

* **Aucune licence.** Le dépôt d'origine n'en porte pas encore. **En l'état, ce
  code est communiqué pour lecture et discussion, et rien d'autre n'est
  concédé** — c'est une question à régler avec Mines Paris-PSL avant toute
  réutilisation.
* **Les données de campagne** (`.npz` de résultats, ~24 Mo) et les figures : hors
  périmètre, ils ne servent pas à lire le code.
* **Les maillages du banc de réplique St Anne** (`impact_fidele_r10.msh`,
  `impact_fidele_s15.msh`) : ils ne sont pas dans le dépôt d'origine non plus.
  Les deux maillages fournis (`meshes/`) suffisent à faire tourner toute la
  vérification.
* **Le code de Solidity** (Imperial College) : il est public sous LGPL-3.0 mais
  n'est pas redistribué ici. Voir [`SOURCES_SOLIDITY.md`](SOURCES_SOLIDITY.md).
* **Les articles cités** : sous droits. `etat_de_l_art/SOURCES.md` en donne la
  référence complète et le DOI.

## Où en est le travail

Un **contre-audit adversarial** a été mené le 2026-08-30 sur le bilan de rockim
contre l'état de l'art. Il a noté ce bilan **52/129**, et ses corrections font
foi : [`etat_de_l_art/chantier/CONTRE_AUDIT_corrections.md`](etat_de_l_art/chantier/CONTRE_AUDIT_corrections.md).

Trois défauts de fond en sont sortis, corrigés dans le code de ce dossier :

* **le septième poste du bilan d'énergie n'existait pas** — l'énergie
  gravitaire, que la littérature pose explicitement. Mesuré : *le résidu du
  bilan ÉTAIT ce poste manquant, à 0,02 % près*, et le garde-fou d'énergie
  coupait un run parfaitement sain à cause de lui ;
* **le diagnostic de « diffuse ratcheting » n'existait qu'en 2D** et sur un essai
  de calibration — donc nulle part sur le cas comparé à Imperial. Porté sur la
  percussion 3D, il montre **× 3,24** de joints endommagés sous insertion
  intrinsèque, et un pic d'effort inférieur de 8,1 % ;
* **la raideur tangentielle du contact n'entrait pas dans le budget de pas de
  temps** du solveur 3D.

Chacun a sa fiche dans `etat_de_l_art/chantier/`, avec sa source, sa preuve de
bit-identité, sa mesure et ses réserves.

---

*Questions, réserves et contre-exemples bienvenus — en particulier sur les points
que les fiches signalent elles-mêmes comme non tranchés.*
