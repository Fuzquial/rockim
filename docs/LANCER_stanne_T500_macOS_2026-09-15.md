# Lancer le St Anne complet (500 µs) sur un Mac — 15/09/2026

Le deck est prêt : `configs/stanne2025_rock137_visc0_T500.cfg`. Il ne diffère du run terminé que par
deux clés, `T = 5.0e-4` et `frames = 33`. **Rien n'est lancé**, ce document donne la marche à suivre.

⚠️ **Je n'ai pas pu tester ces commandes sur un Mac**, faute d'en avoir un. Elles sont écrites depuis
le `CMakeLists.txt` et non de mémoire, mais l'étape 5 est là pour que vous mesuriez votre machine
avant d'engager plusieurs jours de calcul, au lieu de croire mon estimation.

## 1. Les outils, une fois

```bash
xcode-select --install
brew install cmake libomp gmsh
python3 -m pip install gmsh numpy matplotlib scipy
```

`libomp` est le point délicat. **Apple clang ne fournit pas OpenMP.** Sans lui le `CMakeLists`
affiche « rockim: OpenMP introuvable — build SERIE » et le code tourne sur **un seul cœur**. Pour un
run de plusieurs jours, ce n'est pas une option : vérifiez que ce message n'apparaît pas.

## 2. Le dépôt

```bash
git clone --branch g1 https://github.com/Fuzquial/rockim.git
cd rockim
```

Le dépôt est public, aucune authentification n'est nécessaire.

## 3. Compiler

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOpenMP_ROOT=$(brew --prefix libomp)
cmake --build build -j 8
```

Si OpenMP n'est toujours pas trouvé :

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=$(brew --prefix libomp)/lib/libomp.dylib
```

Le binaire est `build/rockim`, **sans `.exe`**. Hors MSVC le `CMakeLists` compile en `-O3 -Wall
-Wextra`. Eigen se règle seul : copie locale, puis paquet système, puis téléchargement des en-têtes.

## 4. Vérifier que la construction est saine

**Les empreintes de bit-identité ont été prises sous MSVC à 4 fils. Elles différeront sur Mac**, et
ce n'est pas une régression : compilateur, bibliothèque mathématique et ordre de réduction changent.
Ce qui doit valoir, ce sont les **critères physiques** :

```bash
python3 tools/bitid.py --exe build/rockim --only fdem_ucs_yan_adaptive_court
```

Le pic doit tomber vers 41,7 MPa. Contrôle plus proche du cas visé :

```bash
python3 tools/bitid.py --exe build/rockim --only fdem3d_yang_v2_court
```

## 5. Mesurer VOTRE machine avant d'engager trois jours

C'est l'étape à ne pas sauter. Copiez le deck en réduisant la durée :

```bash
sed 's/^T = 5.0e-4/T = 2.0e-5/; s/^frames = 33/frames = 2/' configs/stanne2025_rock137_visc0_T500.cfg > /tmp/calib20.cfg
```

Régénérez d'abord le maillage, qui n'est pas dans le dépôt mais se refait à l'identique, graine
fixée :

```bash
python3 tools/make_impact_mesh.py meshes/impact_yang_train1_rock137_hxt.msh 1.0 2e-5 1.37 gap=2e-5 quality=hxt train=fixed
```

Puis le calibrage, en remplaçant 8 par votre nombre de cœurs de performance :

```bash
OMP_NUM_THREADS=8 ./build/rockim /tmp/calib20.cfg /tmp/out_calib20
```

**Comment lire le résultat.** Notez le temps mur. Sur le poste Windows à 14 fils, les 20 premières
microsecondes ont pris **17 minutes**. Le rapport entre votre temps et ces 17 minutes est votre
facteur d'échelle, à appliquer aux totaux du tableau ci-dessous.

⚠️ Ce facteur est **optimiste** : les 20 premières microsecondes sont la partie la moins chère du
calcul, avant toute fissuration. La cadence du run complet s'est effondrée d'un facteur 18 entre le
début et la fin, de 1,20 à 0,065 µs par minute, parce que le contact coûte de plus en plus cher à
mesure que les facettes rompent.

## 6. Le run complet

```bash
caffeinate -i env OMP_NUM_THREADS=8 ./build/rockim configs/stanne2025_rock137_visc0_T500.cfg out_stanne2025_rock137_T500
```

`caffeinate -i` empêche la mise en veille. **Sans lui, un Mac portable s'endort et le calcul
s'arrête**, sans reprise possible puisque rockim n'a pas de sauvegarde d'état.

**Coût, mesuré sur le poste Windows à 14 fils** :

| Cible | Supplément sur 300 µs | Total |
|---|---:|---:|
| 450 µs, la vitesse de rebond | 38,5 h | 67,8 h, 2,8 jours |
| **500 µs, plus l'arrêt des latérales** | **51,3 h** | **80,6 h, 3,4 jours** |
| 750 µs, tout leur enregistrement | 115,4 h | 6,0 jours |

Sortie : 34 trames de 94 Mo, soit **3,2 Go**.

## 7. Trois réserves propres au portable

- **La chauffe.** Trois jours à pleine charge font descendre la fréquence. Le run finira, mais plus
  lentement que le calibrage ne le laisse croire. Une surface dure et ventilée, pas un lit.
- **L'alimentation.** Sur batterie, macOS bride les cœurs de performance. Rester branché.
- **Aucune reprise.** Si la machine s'endort, redémarre ou manque de place, tout est perdu : le deck
  n'est lu qu'au lancement et les trames ne portent pas assez d'état pour repartir. C'est la raison
  du `caffeinate`, et une raison de préférer un poste fixe pour un run de trois jours.

## 8. Dépouiller

Les scripts de figures tournent sans recompiler et lisent le dossier de sortie :

```bash
python3 tools/fig_retournement.py out_stanne2025_rock137_T500 --stem results/fig/T500_retournement
python3 tools/fig_kinetics.py     out_stanne2025_rock137_T500 --stem results/fig/T500_cinetique
python3 tools/fig_evolution_cratere.py out_stanne2025_rock137_T500 --frames 0-33 --stem results/fig/T500_evolution
```

C'est `fig_kinetics.py` qui donnera enfin la **pente de rebond** comparable aux 5,6 m/s de Yang, la
seule grandeur que le run à 300 µs n'a pas pu mesurer. Les autres commandes sont dans
`docs/FIGURES_pour_le_rapport_2026-09-13.md`.

## 9. Ce qui ne marchera pas sur Mac

- `tools/build.ps1` : PowerShell et Visual Studio. Remplacé par le CMake du §3.
- `tools/queue_*.sh` : elles testent la présence d'un processus par `tasklist` et l'arrêtent par
  `Stop-Process`. Sur Mac, `kill -0 "$PID"` et `kill "$PID"` fonctionnent normalement, c'est
  l'inverse du piège Windows.
- `rockim_g1y19.exe` : binaire Windows, inutile ici.
