# Vérification indépendante de la correction du contact par potentiel (rockim)

Date : 2026-10-07. Dépôt `/home/user/rockim`, HEAD 5fae24b + modifications non commitées (src/, include/,
tools/). Aucun fichier du dépôt n'a été modifié. Tout est dans `test_contact/` (scratchpad).

Binaires recompilés ici, Release, g++ Linux, OpenMP, mêmes en-têtes Eigen 3.4.0 :
- référence : `ref/` = `git archive HEAD`, binaire `ref/bld/rockim` ;
- nouvelle version : `new/` = copie des fichiers suivis de l'arbre de travail (vérifié : seuls les 7 fichiers
  src/include du `git diff` diffèrent de `ref/`), binaire `new/bld/rockim`.

Les binaires de l'implémenteur (`build_contact/`, `correction_contact/ref_bld/`) n'ont pas été utilisés.

Conditions machine : 4 cœurs, partagés pendant toute la campagne. Une autre suite de non-régression, lancée
depuis `build_contact/rockim` par un autre processus, tournait en même temps (charge 2,5 à 7,8). Les durées
sont donc indicatives. Les deux calculs à dt/2 ont perdu leur trame VTU n° 6 (fichiers vides écrits à la
même seconde, 07:32, par les deux runs). C'est un incident d'écriture lié à l'environnement : les journaux
et les autres trames sont complets.

## Verdicts

| point | objet | verdict |
|---|---|---|
| 1 | non-régression (suite rapide, bit-identité au défaut) | confirmé |
| 2 | `selftest-potcontact2d` (gradient, boucle fermée, dérive à la naissance) | confirmé |
| 3 | effet physique sur le tunnel (contact non injectant, résidu B4 faible) | confirmé, avec une nuance sur l'interprétation de gcWork < 0 |
| 4 | sensibilité résiduelle à p | partiel : l'injection ne dépend plus de p, le faciès en dépend encore |

## 1. Non-régression

### 1.1 Suite rapide (`tools/verify_suite.py --tier fast`, OMP = 1, lancée depuis chaque copie)

| binaire | tests | verdict | échec |
|---|---|---|---|
| référence | 51 | 50/51 | `t1_toolcontact_penalty` (toolinj 2,98519 contre 4,42633 ; toolvb 10,3848 contre 8,2484 ; broken 3 contre 2) |
| nouvelle | 52 | 51/52 | le même, avec des valeurs identiques |

- Les 51 repères communs ont des champs `ok`, `detail` et `meas` strictement identiques dans les deux JSON
  (`suite_ref.json`, `suite_new.json`). Aucun écart.
- Le nouveau repère `selftest_potcontact2d` passe : potx_ramp = 1,00013, potx_off = 2,24555e-6.
- Sorties des runs de suite : 1 098 fichiers côté référence, 1 099 côté nouvelle (le fichier en plus est le
  CSV du nouveau selftest). `diff -rq` ne trouve d'écart que dans les 43 `config_effective.cfg`.

### 1.2 Calculs courts comparés en entier (OMP = 1, toutes les nouvelles clés au défaut)

Decks dans `bit/`, sorties `bit/out_<cas>_{ref,new}`, journaux `bit/<cas>_{ref,new}.log`.

| cas | deck | fichiers de sortie | fichiers différents | lignes stdout différentes |
|---|---|---|---|---|
| tunnel FDEM 2D court | `red_tip16` (frottement par défaut), T = 0,03 s | 33 | `config_effective.cfg` seulement | avertissement `gcBirth = ramp` (ajouté) ; message de signe du biais (correctif 5) ; `integration` +136 005 → +0,0683 J/m ; `residu` −136 005 J/m [CHECK] → −2,0e-11 J/m [OK] |
| essai 2D Yan | `configs_yan/bd_adap.cfg` (brésilien, adaptatif, contact pénalité), T = 6e-4 s | 89 | `config_effective.cfg` seulement | aucune hors comptage des clés et temps mural |
| essai 2D Yan sous potentiel | même deck + `contact = potential`, `gcActivation = adaptive` | 89 | `config_effective.cfg` seulement | message de signe du biais seulement |
| fdem3d | `configs/fdem3d_percussion.cfg` + `contact = potential`, T = 2e-5 s | 68 | `config_effective.cfg` seulement | aucune hors chronos `tGrid`/`tLoop`. Les compteurs de paires du potentiel 3D (85 851 043 paires, 34 799 clips avec force) sont identiques |
| tunnel `tip16_mu0`, T = 0,12 s, 2 fils | `tun/mu0_a.cfg` | 33 | `config_effective.cfg` seulement | avertissement ; message de signe ; ligne « naissances » ajoutée ; `integration` +1,03948e6 → +511 820 J/m ; `residu` −527 664 J/m [CHECK] → −1,38e-7 J/m [OK] |

Dans tous les `config_effective.cfg`, l'écart se limite au chemin de l'exécutable, au nombre de clés
consommées et aux lignes `# contactCandidates = active (defaut)` et `# potForceExact = false (defaut)`
(une seule en 3D). VTU, CSV de fin, `history.csv` et `frames.csv` sont identiques au bit près. Les seules
lignes qui changent de valeur sont celles du bilan des rouleaux (correctif 3). Les autres différences sont
du texte.

Le compteur de naissances du défaut (3 689 paires, 1,45186e6 J/m) reproduit le ΣdE_sc de l'enquête
(1,452e6 J/m), mesuré par une instrumentation indépendante.

## 2. Selftest `selftest-potcontact2d` (binaire `new`, 1,1 s, code retour 0)

| critère | annoncé | mesuré | seuil | |
|---|---|---|---|---|
| P0 sommet commun sans recouvrement | 0 | 0 (colinéaire, lu : 2,82e-7) | ≤ 1e-12 | ok |
| P0b éléments fins à \|X\| ~ 100 | −4e-21 | −4,11e-21 | > −1e-9 | ok |
| P1 gradient, Munjiza / exact | 0,157 / 2,2e-10 | 0,156586 / 2,24472e-10 | > 1e-2 / < 1e-6 | ok |
| P2 boucle fermée, Munjiza / exact | 2,85e-4 / 5,1e-17 | 2,85482e-4 / 5,14454e-17 | > 1e-5 / < 1e-12 | ok |
| P3 rampe : énergie créée / E(S0) | 1,00013 | 1,00013 (E(S0) = 0,544013) | > 0,5 | ok |
| P3 offset + exact : max \|H − H0\|/H0 | 2,2e-6 | 2,24555e-6 (H_fin − H0 = 7,9e-9) | < 1e-3 | ok |
| P3 offset + Munjiza (lu) | 0,12 | 0,121009 | — | — |

Remarque : H0 = KE0 = 0,0038, alors que E(S0) = 0,544. Rapportée à E(S0), la dérive vaut 1,6e-8, ce qui
rend le critère d'autant plus strict.

Ce selftest a été écrit par l'implémenteur. Je l'ai donc complété par un banc indépendant
(`indep/chk.cpp`), qui inclut `PotentialContact.hpp` tel quel :
- énergie `pairEnergy` contre une quadrature brute de p∫(φA+φB) sur une grille, avec φ = `Bary::phi` du
  code historique : écart maximal 1,46e-3 avec N = 3 000, 7,3e-4 avec N = 6 000, sur 22 paires. L'écart
  décroît en O(1/N), ce qui est l'erreur de la grille : l'énergie exacte est confirmée ;
- forces `pairForce` + `exactGradientForces` contre −dE/dx par différences finies centrées de `pairEnergy` :
  1 101 paires à |X| ~ 100, échelles de 1e-3 à 1. Écart maximal 4,9e-6 de |f|max, au niveau de l'arrondi
  des différences finies à |X| = 100 (pas de 1e-9 sur des coordonnées de 1e2). Pas de défaut.

## 3. Effet physique : tunnel `tunnel_hs_red.msh`, T = 0,12 s, 2 fils, nice -n 10

Decks `tun/*.cfg` : deck d'origine, T = 0,12 et meshFile absolu ; les clés sont ajoutées en fin de fichier.
(a) = défaut ; (b) = `contactCandidates = vertex` ; (c) = (b) + `gcBirth = offset` + `potForceExact = true`.
p = 1e10 N/m, dt = 5,83e-6 s.

Notations. gcWork = travail net imprimé par le contact (le signe + correspond à une injection). Travail
normal = gcWork + frottement dissipé. La ligne « contact » du B4 vaut −gcWork. E(S0) = énergie de
recouvrement à la naissance (compteur imprimé). v max = max |v| nodal, lu dans la dernière trame VTU et sur
l'ensemble des trames.

| run | gcWork (J/m) | frottement dissipé (J/m) | travail normal (J/m) | naissances ; E(S0) (J/m) | cohésif (J/m) | Cundall (J/m) | résidu B4 (J/m) | joints insérés / rompus | KE fin (J/m) | v max fin / max (m/s) | durée |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tip16_mu0 (a) | **+1,309e6** | 0 | +1,309e6 | 3 689 ; 1,452e6 (matérialisés) | 392 838 | 3,334e6 | −1,4e-7 | 6 167 / 4 889 | 139 034 | 11,6 / 24,7 | 92 s |
| tip16_mu0 (b) | −8,81e4 | 0 | −8,81e4 | 15 364 ; 16,0 | 317 944 | 1,722e6 | 3,6e-8 | 5 556 / 4 350 | 137 600 | 1,78 / 2,32 | 217 s |
| tip16_mu0 (c) | −8,27e4 | 0 | −8,27e4 | 15 560 ; 15,3 (neutralisés) | 316 551 | 1,718e6 | −1,6e-8 | 5 650 / 4 317 | 130 400 | 2,03 / 2,46 | 231 s |
| red_tip16 (a) | **+4,36e5** | 398 905 | **+8,35e5** | 3 150 ; 9,05e5 | 358 100 | 2,353e6 | 1,1e-7 | 5 965 / 4 686 | 103 600 | 2,09 / 4,35 | 90 s |
| red_tip16 (b) | −3,67e5 | 271 357 | −9,6e4 | 15 003 ; 12,5 | 301 000 | 1,424e6 | 9,8e-9 | 5 635 / 4 362 | 109 300 | 1,78 / 3,18 | 145 s |
| red_tip16 (c) | −3,60e5 | 264 551 | −9,6e4 | 15 007 ; 12,8 (neutralisés) | 304 700 | 1,409e6 | −1,4e-7 | 5 622 / 4 396 | 104 900 | 2,04 / 2,73 | 182 s |

Signes : la ligne « contact » du bilan B4 vaut −gcWork. Exemple : `red_tip16` (a) imprime contact = −436 272 J/m,
c'est-à-dire gcWork = +4,36e5 J/m, dont 398 905 J/m dissipés par frottement ; la part normale vaut donc
+8,35e5 J/m injectés.

Le travail injecté par le contact sans frottement, ou la part normale avec frottement, passe de +1,31e6
(sans frottement) et +8,35e5 J/m (avec) à −8e4 / −1e5 J/m dès `vertex`. L'énergie de naissance tombe de
1,45e6 à environ 15 J/m, soit 0,005 % du cohésif. Le pic de vitesse tombe de 24,7 à 2,5 m/s sans frottement
et de 4,35 à 2,7 m/s avec. Le résidu B4 vaut au plus 1,4e-7 J/m dans tous les cas, contre −527 664 J/m
[CHECK] avec le binaire de référence. Les chiffres de l'implémenteur sur `tip16_mu0` (défaut, vertex, les
trois clés) sont retrouvés à l'identique.

Coût (machine chargée) : sans frottement ×2,4 (b) et ×2,5 (c) ; avec frottement ×1,6 et ×2,0. Tous les runs
durent moins de 25 min.

### Contrôle supplémentaire : gcWork < 0 n'est-il qu'un biais du compteur ?

Le biais O(dt) du compteur est négatif (enquête §3.1). Un gcWork négatif pourrait donc masquer une injection
résiduelle. J'ai relancé avec dt/2 (`dtFactor = 0,1`, 41 165 pas) :

| run | gcWork (J/m) | naissances ; E(S0) (J/m) | cohésif (J/m) | joints rompus | résidu B4 (J/m) | v max (m/s) | durée |
|---|---|---|---|---|---|---|---|
| tip16_mu0 (a), dt/2 | +1,133e6 | 3 535 ; 1,200e6 | 295 300 | 4 610 | 9,2e-9 | 8,05 | 855 s* |
| tip16_mu0 (c), dt/2 | −7,68e4 | 15 243 ; 8,9 | 259 900 | 4 233 | −2,0e-10 | 2,48 | 978 s* |

\* deux runs en parallèle, charge 6 à 8.

Le cas (a) à dt/2 reproduit exactement l'enquête (1,13e6 et 1,200e6). Pour (c), gcWork ne varie que de 7 %
quand dt est divisé par 2 (−8,27e4 → −7,68e4). Ce n'est donc pas principalement le biais O(dt) : c'est un
travail réellement absorbé, très probablement l'énergie stockée dans les recouvrements en fin de run
(U_fin, de l'ordre de 3e4 J/m dans l'enquête) et la décharge Cundall. Je ne l'ai pas décomposé, faute
d'instrumentation. Conclusion : le contact n'est pas injectant, et le signe négatif ne dépend pas de dt.

## 4. Sensibilité résiduelle à p (tip16_mu0, T = 0,12 s)

p = potPenaltyFactor × E. Pour p = 1e11, dt est réduit à 3,65e-6 s (×0,625), si bien qu'une partie de
l'effet de p sur le faciès passe par dt.

| p (N/m) | clés | gcWork (J/m) | E(S0) (J/m) | cohésif (J/m) | joints insérés / rompus | Cundall (J/m) | résidu B4 (J/m) | v max (m/s) | durée |
|---|---|---|---|---|---|---|---|---|---|
| 1e9 | défaut | +1,322e5 | 1,789e5 | 353 200 | 6 551 / 4 947 | 2,252e6 | −7,2e-9 | 10,7 | 92 s |
| 1e10 | défaut | +1,309e6 | 1,452e6 | 392 838 | 6 167 / 4 889 | 3,334e6 | −1,4e-7 | 24,7 | 92 s |
| 1e11 | défaut | +7,131e6 | 8,17e6 | 471 700 | 5 551 / 4 569 | 9,476e6 | 4,1e-7 | 27,8 | 140 s |
| 1e9 | (c) | −6,44e4 | 1,28 | 345 300 | 6 412 / 4 854 | 1,994e6 | 2,9e-11 | 3,33 | 192 s |
| 1e10 | (c) | −8,27e4 | 15,3 | 316 551 | 5 650 / 4 317 | 1,718e6 | −1,6e-8 | 2,46 | 231 s |
| 1e11 | (c) | −8,03e4 | 116,5 | 240 200 | 4 974 / 3 699 | 1,564e6 | −3,0e-8 | 2,83 | 262 s |

- **Injection : ne dépend plus de p.** Avec les correctifs, gcWork vaut −6,4e4, −8,3e4 et −8,0e4 J/m sur deux
  décades de p. Au défaut, il passait de +1,3e5 à +7,1e6 J/m, soit ×54. L'énergie de naissance reste
  proportionnelle à p (1,3 → 15 → 117 J/m), mais elle est négligeable (au plus 0,05 % du cohésif) et
  neutralisée par `offset`. Cundall et v max ne suivent plus p non plus. Les chiffres de l'implémenteur
  (−6,44e4 / −8,03e4 ; 1,3 / 117 J/m) sont retrouvés.
- **Faciès : dépend encore de p, et nettement.** Entre p = 1e9 et 1e11, les joints rompus passent de 4 854
  à 3 699 (−24 %), les joints insérés de 6 412 à 4 974 (−22 %) et le cohésif de 345 à 240 kJ/m (−30 %).
  Au défaut, l'écart sur les joints rompus était plus faible (4 947 → 4 569, −8 %), parce que l'injection
  croissante avec p compensait en partie. Une part de l'écart à 1e11 vient du dt réduit. À p fixé, dt/2
  seul fait bouger (c) de −2 % sur les joints rompus et de −18 % sur le cohésif : le faciès à T = 0,12 s
  n'est donc pas convergé, ni en dt ni en p. La correction supprime l'artefact énergétique, mais elle ne
  rend pas le résultat indépendant de la pénalité : un contact plus raide limite davantage les
  interpénétrations et retient la fissuration.

## 5. Observations de lecture du code (sans incidence sur les verdicts)

- `potentialContact()` sous `vertex` : `static std::vector<long> vmark` est partagé par toutes les
  instances du solveur et n'est pas protégé entre threads. C'est sans effet avec un seul solveur par
  processus, mais fragile.
- `gcBirth = offset` laisse persister les recouvrements nés profonds, qui ne portent aucune force tant que
  E ≤ e0. C'est voulu (aucune énergie créée), et sous `vertex` ces naissances deviennent rares (E(S0) de
  15 J/m).
- `offset` et `potForceExact` ne sont pas portés en 3D, comme le documente l'implémenteur.

## Fichiers

`bit/` (decks, journaux et sorties de bit-identité) ; `suite_{ref,new}.{log,json}` (sorties de suite dans
`{ref,new}/rockim_suite_*`) ; `selftest/st.log` ; `indep/chk.cpp`, `indep/chk.log` ; `tun/*.cfg`,
`tun/*_new.log`, `tun/out_*`, `tun/times.txt` ; extracteur `scripts/extract.py`.
