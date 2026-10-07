# Série de maillages à train figé et masses du train — T3, campagne de correction du 13/09/2026

Motif : `DIAGNOSTIC_ICL_independant_2026-09-12.md` §4 (« le jumeau s = 2,5 a un piston 26,5 % plus léger que
le s = 1 : comparer ces deux runs comme une étude de résolution est invalide ») et §6.5 (« conserver les
corps métalliques identiques entre niveaux de maille, étudier au moins trois résolutions de roche ») ;
`ECARTS_guo2014_rockim_2026-09-13.md` ligne 23 (« arête médiane 1,37 mm dans la boule, 14 722 tétras au lieu
de ~35 000 : SR ≈ 0,73 »). Cadrage : `CAMPAGNE_correction_2026-09-13.md` §T3. **Aucun run lancé** : seuls
des lancements de lecture (2 µs et 0,5 µs, `tools/deck_smoke.py`, 4 fils) ont servi à vérifier que le
solveur lit les fichiers fusionnés et à mesurer le pas de temps.

## 1. Pourquoi le train dépendait de SR (mesure)

`tools/make_impact_mesh.py` applique un champ de fond (boules 1 mm × SR à R 12,5, 2 mm × SR à R 25, seuil
2 → 10 mm × SR de R 25 à R 100) à TOUT le modèle, et `Mesh.MeshSizeExtendFromBoundary = 0` fait que les
tailles prescrites aux points OCC des corps (`set_pts` : bit 3 mm × s, insert 0,7 × s, piston 5 × s…) ne
gouvernent que les COURBES. L'intérieur et les faces latérales de l'insert, du bit et du piston suivent donc
le champ de la roche, à l'échelle SR. Mesuré avec `tools/mesh_quality.py` sur les fichiers existants :

| fichier | s | SR | piston (tétras / kg) | bit | insert |
|---|---|---|---|---|---|
| `impact_yang_s1_pose.msh` | 1 | 1 | 1 033 / 1,057028 | 5 692 / 1,288297 | 9 949 / 0,064598 |
| `impact_yang_s2.5_pose.msh` | 2,5 | 2,5 | 209 / 0,776709 | 960 / 1,080201 | 968 / 0,063134 |
| test s = 2,5, SR = 1 (défaut) | 2,5 | 1 | 822 / 1,046117 | 4 802 / 1,297868 | 9 472 / 0,064570 |

Le piston « 26,5 % plus léger » du diagnostic (0,776709 contre 1,05703 kg) est exactement cet effet : à
s = 2,5 le champ lointain vaut 25 mm et un cylindre Φ26,5 est pavé par 5-6 cordes.

## 2. Les outils (opt-in, défauts bit-identiques)

**`tools/make_impact_mesh.py`** — arguments nommés ajoutés, tous facultatifs ; sans eux le script est
inchangé (sha256 du maillage s = 2,5 régénéré = `a0047b76…` = `meshes/impact_yang_s2.5_pose.msh`, avant et
après édition).

- **`train=fixed`** : (1) le train (insert, bit, piston, circlip, plaque) est maillé SEUL, roche cachée
  (`Mesh.MeshOnlyVisible`), sous une copie du champ de fond à l'échelle **s**, par un seul `generate(3)` —
  gmsh ré-amorce son générateur aléatoire à chaque `generate()`, donc le train ne dépend de rien d'autre ;
  (2) la roche est maillée dans un **second modèle** gmsh (cylindre + point origine + champs à l'échelle
  SR) ; (3) le script écrit la fusion au format v2.2 de gmsh (CRLF, `$PhysicalNames` et tags physiques /
  élémentaires du premier modèle, seuls les nœuds des tétras — comme gmsh avec `Mesh.SaveAll = 0`, 0 nœud
  orphelin vérifié dans les fichiers du défaut).
  Deux schémas en un seul modèle ont été essayés et **réfutés** : champs `Restrict` par volume (le train
  variait avec SR — volume du piston 9,611e-5 contre 9,894e-5 m³ — parce que la perturbation aléatoire de
  gmsh est consommée par les faces de la roche, maillées avant celles du train ; et à SR = 5, roche plus
  grossière, l'insert changeait encore, 989 contre 968 tétras, ce qui exclut une fuite par `Min`) ; passes
  par visibilité et par dimension (`generate(d)` ré-exécute les dimensions inférieures et EFFACE les faces
  de tout le modèle, groupe caché compris : journal gmsh « No tetrahedra in region 1 », 8 nœuds).
  Contrepartie assumée : à SR = s le fichier n'est pas bit-identique au défaut (autre ordre de génération :
  piston 236 tétras / 1,0136e-4 m³ contre 209 / 9,894e-5 à s = 2,5).
- **`srfar=x`** : échelle du seul champ lointain (SizeMax du seuil = 10 mm × x), défaut SR. Mesuré à
  s = 2,5, SR = 1 : roche 44 936 tétras avec `srfar=2.5` contre 98 547 sans, boule identique. Non utilisé
  pour la série (le cadrage demande SR pur).
- **`verbose=1`** : journal gmsh (General.Terminal).

**`tools/mesh_quality.py`** — sortie par défaut inchangée ; options : `--ball R` (tétras dont le centroïde
est à r < R de l'origine, par corps : N, arête moyenne médiane, q10/q90, h min, h médian), `--masses`
(volume × ρ par corps ; ρ des decks `yang2026_*` : rock 2 626, bit/piston/plate 7 850, insert/circlip
15 250 ; `--rho corps=val`), `--yang` (volumes analytiques des corps tels que dessinés par le générateur,
perte de facettisation, masses de Yang 2026 — piston 1,173 kg, bit 1,509 kg — et densité corrigée).

## 3. La série (s = 1, `quality=hxt`, `train=fixed`, jeux `gap=2e-5 gapr=2e-5` comme les decks)

Commande : `python tools/make_impact_mesh.py meshes/impact_yang_train1_<nom>_hxt.msh 1.0 2e-5 <SR> gap=2e-5 quality=hxt train=fixed`

| fichier | SR | tétras | roche | boule R 12,5 (roche) : N / arête moy. médiane / h médian | roche h min | roche < 0,3 mm | masse roche | temps |
|---|---|---|---|---|---|---|---|---|
| `impact_yang_train1_rock25_hxt.msh` | 2,5 | 24 010 | 6 221 | 839 / **3,594 mm** / 1,230 | 0,7654 mm | 0 | 19,247244 kg | 5 s |
| `impact_yang_train1_rock137_hxt.msh` | 1 | 108 667 | 90 878 | 14 030 / **1,425 mm** / 0,493 | 0,2700 mm | 5 | 19,320000 kg | 13 s |
| `impact_yang_train1_rock073_hxt.msh` | 0,73 | 251 460 | 233 671 | 37 018 / **1,033 mm** / 0,360 | 0,2000 mm | 4 125 | 19,327251 kg | 33 s |

sha256 : `349a6109…` (rock25), `ec42bb6e…` (rock137), `40897167…` (rock073).

**Train, identique sur les trois fichiers** (`cmp_train.py` : mêmes nœuds, mêmes surfaces, mêmes tétras,
écart de volume 0) : insert 9 079 tétras (h min 0,2528 mm), bit 5 090 (0,3044), piston 916 (1,693),
circlip 841 (0,561), plaque 1 863 (0,871) — 17 789 tétras, arête médiane insert 1,43 mm. Masses : insert
0,064584, bit 1,288297, piston 1,058054, circlip 0,014237, plaque 0,190551 kg.

Lectures :
- La **roche du `rock137` est identique** (surface et 17 632 nœuds communs sur 17 632, 90 878 tétras) à
  celle de `impact_yang_s1_pose_hxt05.msh`, le maillage du preset hxt du 13/09 : le membre SR = 1 de la
  série est le run s = 1 actuel côté roche. Le train diffère légèrement de celui du hxt05 (insert 9 079
  contre 9 116 tétras, piston 916 contre 923 ; masses +0,1 % et −0,02 %) : c'est un train NEUF, maillé sans
  la roche avant lui.
- `rock073` atteint l'arête médiane visée (1,033 mm pour 1,0 ; 37 018 tétras dans la boule pour les
  ~35 000 attendus d'un « vrai 1 mm ») ; son total (251 460) est celui du modèle de Yang (230 788, ligne 23
  d'ECARTS), avec la réserve du diagnostic §6.5 (leur statistique de maillage n'est pas connue).
- Le coût : SR scale TOUTE la roche (champ lointain 7,3 mm au lieu de 10) : roche × 2,57, h min 0,270 →
  0,200 mm (pas de temps × 0,74 si le sliver de la roche commande). Le pas de temps mesuré par lecture est
  au §5. `srfar=1` garderait le bord à 10 mm (mesuré à s = 2,5 : roche ÷ 2,2 pour la même boule).
- `rock25` : h min de la roche 0,765 mm, mais l'insert (0,2528 mm, carbure) et les joints de la roche
  gardent la main sur le pas de temps — un maillage à roche grossière et train fin ne coûte pas ce que sa
  roche laisse croire (§5).

## 4. Les masses : maillé, analytique, Yang

`tools/mesh_quality.py <msh> --masses --yang` sur la série (train identique, roche SR = 1) :

| corps | N | V maillé [m³] | V analytique [m³] | écart | ρ | masse maillée | Yang 2026 | écart |
|---|---|---|---|---|---|---|---|---|
| rock | 90 878 | 7,357197e-3 | 7,363108e-3 | −0,08 % | 2 626 | 19,320 kg | — | — |
| insert | 9 079 | 4,235024e-6 | 4,260361e-6 | −0,59 % | 15 250 | 0,064584 | — | — |
| bit | 5 090 | 1,641143e-4 | 1,709183e-4 | **−3,98 %** | 7 850 | 1,288297 | 1,509 | **−14,6 %** (analytique −11,1 %) |
| piston | 916 | 1,347839e-4 | 1,434019e-4 | **−6,01 %** | 7 850 | 1,058054 | 1,173 | **−9,8 %** (analytique −4,0 %) |
| circlip | 841 | 9,335946e-7 | 9,330530e-7 | +0,06 % | 15 250 | 0,014237 | — | — |
| plate | 1 863 | 2,427399e-5 | 2,424739e-5 | +0,11 % | 7 850 | 0,190551 | — | — |
| bit + insert + circlip | | | | | | 1,367119 | 1,509 | −9,4 % (analytique 1,420909 : −5,8 %) |

Contrôle de l'outil : ses masses reproduisent celles du résumé du solveur à 6 chiffres sur trois maillages
(s = 2,5 `results/yang_bench_s25_v3P.log` : piston 0,776709, bit 1,0802, insert 0,0631335, rock 19,2472,
circlip 0,0142247, plate 0,191708 ; s = 1 `results/yang_adaptive_smoke.log` : 1,05703 / 1,2883 / 0,0645975 /
19,32 / 0,0142373 / 0,190551 ; maillage fusionné `train=fixed` s = 2,5 lu par `rockim_g1y16.exe` 2 µs :
0,795709 / 1,0802 / 0,0631335 / 19,2472 / 0,0142247 / 0,191708) ; la variante fausse (insert en acier)
donne 0,0325 kg et échoue comme il se doit.

**D'où vient l'écart.** Trois causes, mesurées séparément :

1. **Facettisation** (le maillage seul, −6,0 % piston, −4,0 % bit). Les faces latérales des cylindres sont
   maillées à la taille du champ de fond et non à celle des points OCC (§1) : à s = 1 le piston, à plus de
   100 mm de l'origine, reçoit 10 mm alors que ses deux cercles d'extrémité sont à 5 mm. La perte d'aire
   d'un polygone inscrit vaut ≈ (h/R)²/6 par corde : piston R 13,25, h 10 → 9,5 % à mi-longueur, 2,4 % aux
   extrémités, −6,0 % mesuré ; bit R 15, h 2 → 10 mm le long du fût → jusqu'à 7,4 %, −4,0 % mesuré ; roche
   R 125, h 10 → 0,1 % (−0,08 % mesuré). L'anneau (circlip) et la plaque percée sont à +0,1 % : les
   **évidements** (trou Φ31 de la plaque, anneau) sont exacts et ne sont pas en cause.
2. **Longueurs et géométrie du dessin** (−4,0 % piston, −11,1 % bit, même sans maillage). Les cylindres
   EXACTS de la fig. 5 tels que codés (piston Φ26,5 × 260, bit Φ30 × 265 insert compris, insert
   hémisphère R 8,51 + fût Φ15,88) pèsent 1,1257 kg et 1,3417 kg (1,4209 avec insert et circlip) : les
   masses publiées (1,173 et 1,509) demandent 4,2 % et 12,5 % (6,2 %) de volume en plus. Pour le piston
   seul cela vaut +10,9 mm de longueur (270,9) ou Φ27,05 ; pour le bit, un épaulement ou une tête plus
   large que Φ30 (gorge du circlip, portée de la plaque) que le dessin simplifié V1 du générateur omet, ou
   un **périmètre** différent de « bit 1,509 kg » (diagnostic §4 : à vérifier avec les auteurs).
3. Par conséquent, raffiner les faces latérales (récupérer 6 % et 4 %) ne suffirait pas : il manquerait
   encore 4 % et 11 %.

**Correction proposée (deck, sans code ni maillage nouveau)** : densité par corps via deux phases acier
supplémentaires — les densités sont celles du train figé de la série (`--yang`) :

```
phases = rock steel steelPiston steelBit carbide
# piston : 7850 x 1,173 / 1,058054
phase.steelPiston.rho = 8702.8
phase.steelPiston.E   = 200e9        # option (rho, E) : 221.7e9 = 200e9 x 1.10863, voir ci-dessous
phase.steelPiston.nu  = 0.29
phase.steelPiston.ft  = 1e12
phase.steelPiston.cohesion = 1e12
# bit : 7850 x (1,509 - 0,064584 - 0,014237) / 1,288297 si 1,509 kg = bit + insert + circlip
phase.steelBit.rho = 8714.5           # 9194.8 si 1,509 kg = le seul corps « bit »
phase.steelBit.E   = 200e9           # option (rho, E) : 222.0e9
phase.steelBit.nu  = 0.29
phase.steelBit.ft  = 1e12
phase.steelBit.cohesion = 1e12
groupPhase.piston = steelPiston
groupPhase.bit    = steelBit
# groupPhase.plate reste steel (sa masse n'est pas publiée)
```

Choix entre « ρ seul » et « ρ et E ensemble » : avec ρ seul, c = √(E/ρ) baisse de 5,1 % (5 047 → 4 794 m/s)
et l'impédance ρc monte de 5,3 %, la durée 2L/c de l'onde du piston s'allonge de 5 % ; avec (ρ, E) scalés du
même facteur 1,1086, c et la durée sont inchangées et ρc monte de 10,9 % — c'est ce que ferait un piston
réel de section plus grande (la masse manquante vient d'une section ou d'une longueur non dessinées ; une
section plus grande laisse c intact). Recommandation : **(ρ, E) ensemble**, à confirmer par les auteurs.
Effets collatéraux : les corps acier sont `groupContinuum` (pas de loi de joint) ; `potStiffnessByPhase =
min` fait que la raideur de contact insert/roche ne change pas et que celle piston/bit monte de 11 % ; le pas
de temps des éléments acier (∝ h/c) est inchangé sous (ρ, E). Quantité de mouvement incidente : 1,173 × 9 =
10,56 N·s contre 1,058 × 9 = 9,52 (−9,8 %) dans tous les runs faits jusqu'ici. Les `phase.<n>.fraction` sont
inopérantes avec `mesh = file` (doc §5.2) ; la lecture d'un deck à cinq phases par le solveur reste à
vérifier par un lancement de 2 µs (`tools/deck_smoke.py`), non fait ici (decks = T4).

## 5. Pas de temps et coût par pas des trois maillages (lecture 0,5 µs, deck témoin `yang2026_bench_s25_v3P.cfg`, 4 fils)

`python tools/deck_smoke.py --exe rockim_g1y16.exe --T 5e-7 --frames 1 --threads 4 --mesh meshes/impact_yang_train1_<nom>_hxt.msh configs/yang2026_bench_s25_v3P.cfg`
(deck témoin tel quel, donc `jointPenaltyLength = min` : la pénalité de TOUS les joints suit le plus petit
h inscrit du fichier, qui est ici celui de l'insert, 0,2528 mm, sur les trois maillages).

| maillage | tétras | joints | nœuds | dt | pas pour 200 µs | mur de la lecture (4 fils, init compris) | masses piston / bit / insert lues |
|---|---|---|---|---|---|---|---|
| rock25 | 24 010 | 43 772 | 96 040 | 2,983 ns | 67 043 | 37 s pour 168 pas | 1,058 / 1,288 / 0,06458 kg |
| rock137 | 108 667 | 209 380 | 434 668 | 1,112 ns | 179 890 | 116 s pour 450 pas | 1,058 / 1,288 / 0,06458 kg |
| rock073 | 251 460 | — | — | **non mesuré** (lecture encore en cours à la remise du rapport) | — | — | — |

Lecture : de rock25 à rock137 le pas de temps est divisé par 2,7 alors que h min de la roche passe de 0,765
à 0,270 mm (÷ 2,8) — c'est bien le sliver de la ROCHE qui commande, l'insert (0,2528 mm, carbure, sans
joint) non ; le coût par pas (mur/pas, init comprise) passe de ~0,22 à ~0,26 s à 4 fils pour 4,5 × plus de
tétras, l'initialisation dominant ces lectures courtes. Pour rock073 (h min 0,200 mm), extrapolation à
confirmer : dt ≈ 0,82 ns, ~243 000 pas pour 200 µs, coût ≈ 2,3 × (tétras) × 1,35 (pas) ≈ 3,1 × celui de
rock137. Les masses lues par le solveur sont celles de `--masses` (train identique).

## 6. Tests exécutés (chiffres)

| test | critère | résultat |
|---|---|---|
| défaut du mailleur après édition | sha256 de `impact_yang_s2.5_pose.msh` régénéré = fichier existant | `a0047b76…` = `a0047b76…` (3 régénérations, avant et après chaque édition) |
| `train=fixed` à SR = 2,5 / 1 / 5, s = 2,5, Delaunay | train (5 corps) : tétras, nœuds, surface, volume identiques | identiques (insert 968, bit 927, piston 236, circlip 352, plaque 562 ; volumes à 0) |
| idem sous `quality=hxt`, SR = 2,5 / 1 | identiques | identiques (808 / 658 / 119 / 171 / 288) ; HXT déterministe (deux générations, sha `0463cb81…`) |
| défaut à SR = 1, s = 2,5 (DOIT échouer) | train différent | échoue comme prévu : piston 822 tétras au lieu de 209, insert 9 472 au lieu de 968 |
| série s = 1 : rock25 / rock137 / rock073 | train identique | identique (9 079 / 5 090 / 916 / 841 / 1 863 tétras, volumes à 0) |
| format fusionné | sections, CRLF, tétras seuls, 0 orphelin, tags 1-6 | `$MeshFormat/$PhysicalNames/$Nodes/$Elements`, CRLF vrai, type 4 seul, 2 655 nœuds tous référencés, tags 1-6 |
| lecture par le solveur (`rockim_g1y16.exe`, deck témoin, maillage fusionné s = 2,5, T = 2 µs) | démarre, code 0, masses = outil | code 0, 708 pas, dt 2,828 ns, 28 s ; 6 masses = `--masses` à 6 chiffres |
| `--masses` sur `impact_yang_s2.5_pose.msh` et `_s1_pose.msh` | = résumé du solveur | 12/12 masses reproduites (§4) ; variante insert en acier : 0,0325 ≠ 0,0631 (échoue) |
| `--ball 0.0125` sur `impact_yang_s1_pose.msh` | 14 722 tétras / 1,372 mm (diagnostic §4) | 14 722 / 1,3722 mm |
| caractères de contrôle | aucun octet < 32 hors tab/LF/CR dans les .py édités | aucun |

## 7. Ce qui n'a PAS été fait, et pourquoi

- Aucun run de la série (interdit par le cadrage) ; les chiffres de coût du §5 sont des lectures de 0,5 µs.
- Pas de maillage `srfar` dans la série : le cadrage demande SR pur ; l'option existe et est mesurée à
  s = 2,5 seulement.
- Pas de correction GÉOMÉTRIQUE du train (faces latérales plus fines, épaulement du bit) : elle changerait le
  train commun à tous les runs depuis le 11/09 et ne suffirait pas (§4, cause 2) ; c'est un choix pour
  Fernando. Le périmètre de « bit 1,509 kg » n'est pas connu : deux densités candidates sont données.
- Aucun deck écrit (T4) ; le bloc de phases du §4 est une proposition, sa lecture par le solveur n'est pas
  vérifiée.
- Le train de la série n'est pas bit-à-bit celui de `impact_yang_s1_pose_hxt05.msh` (piston 916 contre 923
  tétras) : la série est cohérente en elle-même (train identique sur ses trois membres), pas avec les runs
  déjà faits.
