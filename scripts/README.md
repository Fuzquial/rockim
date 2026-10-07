# Scripts de compilation et de lancement

Ce dossier regroupe des scripts conservés pour la traçabilité des campagnes d'août et septembre 2026.
Ils ont été déplacés de la racine du dépôt le 2026-10-07 (commit `91e922b`, `git mv`, contenu
inchangé hors la ligne `cd` décrite plus bas). Ils ne sont pas nécessaires pour compiler rockim
aujourd'hui.

La compilation de référence passe par CMake (`CMakeLists.txt` à la racine), décrite dans `README.md`
et `DOCUMENTATION_rockim.md` §2. Sous Windows, `tools/build.ps1` construit par CMake et Ninja avec les
mêmes drapeaux que les scripts ci-dessous.

## `windows/` : compilation MSVC historique

Chaque `build_<tag>.cmd` charge l'environnement de Visual Studio 2022 (`vcvars64.bat`), compile
`src\*.cpp` avec `cl` (`/std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX`, Eigen pris dans
`..\rockim\eigen-3.4.0`) et produit un binaire au nom dédié, `rockim_<tag>.exe`. Le nom dédié laissait
l'ancien binaire disponible pour comparaison et évitait de verrouiller l'exécutable d'un calcul en
cours. Les objets vont dans `obj_<tag>\`.

Avant le déplacement, ils vivaient à la racine et commençaient par `cd /d %~dp0`. Depuis le
2026-10-07, ils commencent par `cd /d %~dp0..\..` : ils remontent
de `scripts\windows\` à la racine du dépôt avant de compiler, si bien que les chemins relatifs
`src\` et `include\` restent valides quel que soit le dossier d'appel. Trois scripts gardent un
chemin absolu du poste du doctorant (`cd /d C:\Users\...\rockim_g1`) : `build_lota.cmd`,
`build_lotb.cmd` et `build_lotb_bench.cmd`, qui compilent un seul fichier ou un banc et non le
solveur complet.

| Script | Binaire ou objet | Objet du lot, d'après le commentaire d'en-tête |
|---|---|---|
| `build_a1.cmd` à `build_a8.cmd`, `build_tun.cmd`, `build_yy.cmd` | `rockim_a1.exe` à `rockim_a8.exe`, `rockim_tun.exe`, `rockim_yy.exe` | étude du tunnel (zone endommagée, excavation) ; le même en-tête est recopié d'un script à l'autre |
| `build_e1.cmd`, `build_e2.cmd`, `build_e3.cmd` | `rockim_e1.exe` | essais hydro-mécaniques AbuAisha : instrumentation, mouillage des joints (`hydroWetDamage`), pompe différée (`hydroStart`) |
| `build_hy.cmd` | `rockim_hy.exe` | correctif du signe hydro (AbuAisha 2017) |
| `build_i1.cmd` à `build_i4.cmd` | `rockim_e1.exe` | impact de Yang : pulvérisation (`bulkDamage = yang`), liaison entre corps (`groupBond`), suivi des corps et jauge (`trackGroups`) |
| `build_e.cmd` | `rockim_e.exe` | lot E : gardes (`toolShape`, NaN, traces) |
| `build_pot.cmd` | `rockim_pot.exe` | étude de coupe PDC avec contact par potentiel |
| `build_f2.cmd` | `rockim_f2.exe` | chantier f2 : joints pré-rompus, microsismicité |
| `build_fix.cmd` | `rockim_fix.exe` | corrections du guide du 2026-09-06 |
| `build_g1.cmd` | `rockim_g1.exe`, ou `rockim_g1ref.exe` avec l'argument `ref` | loi de la note « FDEM hybride à insertion adaptative » |
| `build_dev.cmd`, `build_chk.cmd`, `build_flush.cmd`, `build_size.cmd`, `build_str.cmd` | `rockim_dev.exe`, `rockim_chk.exe`, `rockim_flush.exe`, `rockim_size.exe`, `rockim_str.exe` | binaires de travail, sans commentaire d'en-tête |
| `build_lota.cmd`, `build_lotb.cmd`, `build_lotb_bench.cmd` | objets dans `obj_LOTA\`, `obj_LOTB\` | compilation séparée de `MatLaw.cpp`, de `Fdem3dSolver.cpp` et du banc du lot B |
| `run_cut.cmd` | `rockim_chk.exe` | lance la démonstration de coupe de `calibration_redbohus/` à 14 fils |

Les binaires de la série `rockim_f2*` sont décrits dans `etude_lois_fem/JOURNAL.md` et leurs empreintes
dans `tools/exe_manifest.json`.

## `historique/` : scripts de lancement du poste du doctorant

Scripts bash (Git Bash sous Windows) qui ont lancé des calculs longs en août et septembre 2026. Ils
commencent tous par un `cd` vers un chemin absolu du poste (par exemple
`C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1`) et appellent un binaire daté
(`rockim_yy.exe`, `rockim_e1.exe`, `rockim_hydro.exe`…). Ils ne tournent pas tels quels sur une autre
machine. Ils sont conservés parce qu'ils disent quel deck a été lancé avec quel binaire, dans quel
dossier de sortie.

| Script | Deck et objet |
|---|---|
| `run_3dfin.sh`, `run_3dgrad.sh`, `run_chain_3d.sh` | indentation 3D (`configs/indent3d_*.cfg`) |
| `run_b1_damping.sh` | balayage de l'amortissement de Cundall sur `indent3d_grad` |
| `run_indent2d.sh`, `run_yan2d.sh` | indentation 2D (`configs/indent2d_*.cfg`) |
| `run_impact.sh`, `run_impact_uni.sh`, `run_xi2.sh` | impact St Anne et impact d'insert de type Yang |
| `run_yang_equiv.sh`, `run_yeq_timing.sh`, `run_smoke_yy.sh` | deck équivalent Yang, mesure de temps, fumée d'une liste de decks |
| `run_e1.sh`, `run_e2.sh`, `run_e3.sh` | essais hydro-mécaniques AbuAisha 1 à 3 (`bench_abuaisha/configs/`) |
| `run_f7.sh`, `run_f7_hydro.sh`, `run_f7_prod.sh`, `run_hfs.sh` | fracturation hydraulique AbuAisha, séries f7 et pression de rupture |
| `run_essais_coupe.sh` | essais de coupe 2D du tunnel (`tunnel_edz/configs/cut2d_*.cfg`) |
