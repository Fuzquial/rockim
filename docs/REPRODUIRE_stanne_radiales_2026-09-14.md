# Reproduire le run qui donne les fissures radiales (St Anne 2025) — 14/09/2026

Tout ce qu'il faut pour relancer, depuis un clone frais de `Fuzquial/rockim`, branche `g1`, le calcul
qui a produit les radiales, la fissure médiane et la chute de contrainte sur les lèvres.

## 1. Le binaire

| | |
|---|---|
| Exécutable | `rockim_g1y19.exe` (versionné à la racine, 2,5 Mo) |
| SHA-256 | `4e5594b5e31f1899bc2ebf22aa3e84c55a97553eddb7030ee56573fb74cbce94` |
| Sources correspondantes | **exactement** celles de la branche au commit `bfb1880` et suivants — aucun commit ne touche `src/`, `include/` ni `CMakeLists.txt` après sa construction |
| Ancre de bit-identité | jouée sur `rockim_g1y18` (8/8 + 1/1 IDENTIQUE, `results/bitid_g1y18.log`). `g1y19` = `g1y18` + le seul correctif `jbDir_` du mode `jbMode = cycle`, hors du chemin par défaut |

**Le reconstruire** (Visual Studio 2022 requis, ~1 min) :

```
powershell -ExecutionPolicy Bypass -File tools\build.ps1 -Jobs 8
copy build\rockim.exe rockim_g1y19.exe
```

L'exécutable versionné est une commodité ; **les sources font foi**. Un binaire reconstruit n'aura pas
le même SHA-256 (horodatage PE) mais le même comportement, ce que vérifie `python tools/bitid.py --exe
rockim_g1y19.exe`.

## 2. Le maillage

Les `.msh` restent hors du dépôt (5,5 Mo chacun). Ils se **régénèrent à l'identique** : le mailleur pose
`Mesh.RandomSeed = 1`, il est déterministe.

```
python tools/make_impact_mesh.py meshes/impact_yang_train1_rock137_hxt.msh 1.0 2e-5 1.37 gap=2e-5 quality=hxt train=fixed
```

| Maillage | SR | Tétraèdres | Cœur r < 12,5 mm | Arête médiane du cœur |
|---|---|---|---|---|
| `impact_yang_train1_rock25_hxt.msh` | 2,5 | 24 010 | 839 | 3,59 mm |
| **`impact_yang_train1_rock137_hxt.msh`** | **1,37** | **108 667** | **14 030** | **1,42 mm** |
| `impact_yang_train1_rock073_hxt.msh` | 0,73 | 251 460 | 37 018 | 1,03 mm |

`train=fixed` fige le train de frappe à s = 1 dans les trois : la **seule** variable est la résolution
de la roche (piston 1,058 kg, bit 1,288 kg partout). Contrôle : `python tools/mesh_quality.py <maillage>`.

## 3. Le lancement

```
OMP_NUM_THREADS=14 ./rockim_g1y19.exe configs/stanne2025_rock137_visc0.cfg out_stanne2025_rock137
```

Deck : St Anne 2025 à 10,66 m/s, matériau de leur Table 4 (ρ 2731, E 57 GPa, ν 0,31, ft 7,0 MPa,
c 18,8 MPa, tanφ 1,0, GI 12, GII 800), pénalité `jointPenaltyLength = edge` et facteur 26,316 = p0/(2E)
avec les 3 000 GPa de l'ARMA 2024, loi `plastic` + ratchet, **pulvérisation désactivée** (elle est le
modèle 2026 du granite, absent de l'article 2025), sans viscosité, cap de traction moyenne neutralisé,
`writeRuptureFields = true` et `contactForcePairs = rock:insert piston:bit plate:bit`.

**Coût mesuré** : dt = 1,67 ns ; ~0,5 µs/min au départ, ~0,09 µs/min à 190 µs (le contact grossit avec
le nombre de facettes rompues). Compter **une vingtaine d'heures pour 200 µs** sur 14 fils.

Files d'attente prêtes (un seul gros job à la fois, bornes en temps mur) : `tools/queue_nuit_1314c.sh`
enchaîne ce run puis le Kuru `configs/yang2026_kuru_train1_v5.cfg`.

## 4. Les figures

```
python tools/fig_fp_direct.py out_stanne2025_rock137:"St Anne fin" --stem results/fig/x_fp
python tools/fig_kinetics.py  out_stanne2025_rock137:"St Anne fin" --stem results/fig/x_cin
python tools/fig_joints_only.py out_stanne2025_rock137 --lim 25 --depth 25 --stem results/fig/x_joints
python tools/fig_joints_cuts.py out_stanne2025_rock137 --zs=-0.5,-3,-8 --lw 1.4 --lim 25 --stem results/fig/x_coupes
python tools/fig_stress_section.py out_stanne2025_rock137 --field sigma1 --smooth 1.2 --phase rock --stem results/fig/x_sig_coupe
python tools/fig_stress_section.py out_stanne2025_rock137 --field sigma1 --plan z --z0 -2 --phase rock --stem results/fig/x_sig_dessus
python tools/crater_metrics.py out_stanne2025_rock137 --plot results/fig/x_cratere.png --sectors 16
python tools/crack_paths.py    out_stanne2025_rock137 --lw 0.7 --stem results/fig/x_branches
python tools/fig_icl_stanne_review.py --frame <k> --out output/pdf/x_ICL
python tools/fig_surface_compare.py   --frame <k> --out output/pdf/x_surface
```

Les figures de l'état à 192 µs sont versionnées dans `results/fig/stanne192/` et `output/pdf/stanne_*_192us/`
(PDF vectoriels Computer Modern, directement inclusibles en LaTeX).

## 5. Ce qui n'est PAS dans le dépôt

Les fichiers de résultats (`out_*/`, trames VTU de plusieurs dizaines de Mo chacune), les maillages et les
journaux. Le run se rejoue à l'identique avec les trois commandes ci-dessus ; les mesures chiffrées sont
consignées dans `docs/ECARTS_guo2014_rockim_2026-09-13.md` §10 et dans ce document.
