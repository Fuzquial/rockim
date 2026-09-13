# Installer et faire tourner rockim sur un Mac — 14/09/2026

Le dépôt est multiplateforme : `CMakeLists.txt` a une branche non-MSVC (`-O3 -Wall -Wextra`) et
OpenMP est **facultatif**. L'exécutable `rockim_g1y19.exe` versionné est un binaire **Windows** : il
ne sert à rien sur Mac, il faut recompiler. Compter dix minutes en tout.

## 1. Les outils

```bash
xcode-select --install                 # clang, make (si pas déjà fait)
brew install cmake libomp             # cmake + OpenMP pour Apple clang
brew install gmsh                     # seulement pour régénérer les maillages
python3 -m pip install numpy matplotlib scipy gmsh
```

`libomp` est le point délicat : **Apple clang ne fournit pas OpenMP**. Sans lui le code compile et
tourne quand même, en série, mais quatre à dix fois plus lentement — et les empreintes de
`tools/bitid.py` diffèrent alors de celles du poste Windows (c'est écrit dans le `CMakeLists`).

## 2. Le clone

```bash
git clone https://github.com/Fuzquial/rockim.git
cd rockim
git checkout g1
```

Le dépôt est privé : il faut être authentifié (`gh auth login`, ou une clé SSH et l'URL
`git@github.com:Fuzquial/rockim.git`).

## 3. La compilation

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_ROOT=$(brew --prefix libomp)
cmake --build build -j 8
./build/rockim --help 2>/dev/null || ls -la build/rockim
```

Si CMake ne trouve toujours pas OpenMP (message « rockim: OpenMP introuvable — build SERIE ») :

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
  -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=$(brew --prefix libomp)/lib/libomp.dylib
```

Eigen est géré tout seul : le `CMakeLists` cherche une copie locale, puis le paquet système
(`brew install eigen`), puis télécharge les en-têtes 3.4.0 si le réseau est là.

## 4. Vérifier que la construction est saine

```bash
python3 tools/bitid.py --exe build/rockim --list          # les 8 decks de l'ancre
python3 tools/bitid.py --exe build/rockim --only fdem3d_visc_yan_3d
```

⚠️ Les empreintes de l'ancre ont été prises **sous MSVC, à 4 fils**. Sur Mac elles **différeront** :
compilateur, bibliothèque mathématique et ordre de réduction ne sont pas les mêmes. Ce n'est pas une
régression. Ce qui doit valoir sur Mac, ce sont les **critères physiques**, pas les hachages :

```bash
python3 tools/bitid.py --exe build/rockim --only fdem_ucs_yan_adaptive_court   # pic ≈ 41,7 MPa
OMP_NUM_THREADS=4 ./build/rockim tests_f2/campagne13/S3bis/cyce_plastic.cfg /tmp/cyc
# doit imprimer : W net ∝ Δt (loi saine) — critère (v) du banc de joint
```

## 5. Relancer le calcul des radiales

Le maillage n'est pas dans le dépôt (5,5 Mo), il se régénère à l'identique (graine fixée) :

```bash
python3 tools/make_impact_mesh.py meshes/impact_yang_train1_rock137_hxt.msh \
        1.0 2e-5 1.37 gap=2e-5 quality=hxt train=fixed
OMP_NUM_THREADS=8 ./build/rockim configs/stanne2025_rock137_visc0.cfg out_stanne2025_rock137
```

**Coût** : une vingtaine d'heures pour 200 µs sur 14 fils Windows. Sur un Mac portable, compter deux
à trois fois plus, et surveiller la chauffe. Pour un essai rapide, ramener `T` à 2e-5 dans le deck
(20 µs, quelques minutes) — assez pour vérifier que ça tourne, pas pour voir des fissures.

## 6. Écrire le rapport sans rien relancer

C'est le cas le plus probable. Toutes les figures sont versionnées :

```bash
git clone ... && cd rockim && git checkout g1
ls results/fig/stanne192/*.pdf output/pdf/stanne_*_192us/*.pdf
```

Onze PDF vectoriels en Computer Modern, directement inclusibles en LaTeX, plus les mesures chiffrées
dans `docs/ECARTS_guo2014_rockim_2026-09-13.md` et `docs/REPRODUIRE_stanne_radiales_2026-09-14.md`.
Aucune compilation n'est nécessaire pour cela.

## 7. Ce qui ne marchera pas sur Mac

- `tools/build.ps1` (PowerShell + Visual Studio) — remplacé par les commandes CMake ci-dessus.
- Les files d'attente `tools/queue_*.sh` : elles testent la présence d'un processus par `tasklist` et
  l'arrêtent par `Stop-Process`. Sur Mac, remplacer par `kill -0 "$PID"` et `kill "$PID"`, qui y
  fonctionnent normalement (c'est l'inverse du piège Windows du 14/09).
- `rockim_g1y19.exe` : binaire Windows.
