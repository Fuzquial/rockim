# bitid — bit-identité générique de rockim (mesure B5/B6 du plan de robustesse, 2026-09-05)

## Pourquoi

Jusqu'au 05/09 la preuve « défaut bit-identique » de chaque build était un script `bitid_wNN.sh`
copié-collé par build (exe, decks et nombre de fils codés en dur), l'ancre a changé deux fois
(10 configs à OMP 14 pour w0…w12, un seul deck cdp à OMP 2 pour w14…w18), les hachages ne sont pas
comparables d'un nombre de fils à l'autre, et aucune liste de hachages de référence n'était au dépôt
(constats T03/T08/P17/D08 de `ETAT_DES_LIEUX_2026-09-05.md`). `tools/bitid.py` remplace ces scripts
par **un** outil, **une** liste fixe de decks, **une** ancre versionnée.

## Usage

```
python tools/bitid.py --exe rockim_f2w18.exe            # compare a l'ancre tools/bitid_refs.json
python tools/bitid.py --exe build/rockim.exe --only fem3d
python tools/bitid.py --exe rockim_f2w19.exe --json results/bitid_w19.json
python tools/bitid.py --exe rockim_f2w18.exe --update   # (re)prend l'ancre  -> CHANGELOG + commit
python tools/bitid.py --list                             # les 8 decks, mode, loi
```

Options : `--threads N` (défaut = les fils de l'ancre, deck par deck ; un N différent rend la
comparaison **NON COMPARABLE**, le script le dit), `--only <sous-chaîne>`, `--refs <json>`,
`--timeout <s>` (900), `--keep` (garde les `.vtu`), `--outroot <dossier>`, `--json <rapport>`,
`--apex-retries` (10). L'exe est cherché tel quel puis dans la racine `rockim_f2`. L'exe est lancé
**depuis la racine `rockim_f2`** : les `meshFile` des decks y sont relatifs.

Sorties : par défaut un **dossier temporaire** `rockim_bitid_<exe>_omp<N>_*` créé par
`tempfile.mkdtemp` (même convention que `verify_suite.py` ; le chemin est imprimé en tête), un
sous-dossier par deck (`history.csv`, `frames.csv`, journal `<deck>.log` à côté ; `.vtu` supprimés
après hachage — aussi sur échec ou timeout — sauf `--keep`) ; `--outroot` fixe le dossier. Rien n'est
écrit dans le dépôt hors `--update` : le rapport durable est `--json`. (Les dossiers
`tests_f2/bitid/out_*` des trois passes du 05/09 datent d'avant ce défaut ; `tests_f2/bitid/.gitignore`
les exclut du futur dépôt.) Codes de sortie : 0 = `N/N IDENTIQUE`, 1 = différence ou deck sans
référence, 2 = échec de run / exe introuvable / `--only` sans deck.

Ce qui est haché (SHA-256 du contenu, 64 hex dans le json, 16 affichés) : `history.csv`,
`frames.csv`, le ou les `.vtu` de la **dernière trame** (en fdem/fdem3d : `fdem3d_NNNN.vtu` **et**
`fdem3d_joints_NNNN.vtu`), et `probes.csv` s'il existe. Le pic de force (`peak tool force`), le
temps de calcul, `dt` et le nombre de pas sont relevés depuis la console et écrits dans le json —
une différence est annoncée avec **les deux pics de force**.

## Les 8 decks (`tests_f2/bitid/`, chacun < 5 min à OMP 4, en-tête = mode, loi, source, durée)

| # | deck | mode | loi / chemin exercé | source (raccourci) |
|---|---|---|---|---|
| 1 | `fem3d_cdp_PQ_court` | fem3d | cdp inverse, pénalité, quarterModel, confinement 20 MPa, soupapes + spall | `etude_lois_fem/bitid_w12/PQ_cdpI_P020_court.cfg` inchangé (T 130 µs) ; maillage `CONTINUUM/…/perc3d/Q1_c05.msh` (13 687 tets) |
| 2 | `fem3d_dpr_T1_court` | fem3d | dpr + dpApex + dpTension off + rankineDrive stress + power + crackband + erodeWfrac, **signorini**, **bulkViscosity**, **hencky** | `etude_lois_fem/matrice_v2/C_T1_R_P000.cfg`, T 3e-4 → 2,5e-5 (~12 µs de contact, phase de charge), frames 2 ; `T1_c05_clean.msh` (98 342 tets) |
| 3 | `fem3d_sk2011_cyl_court` | fem3d | saksala2011, cylindre interne Kuhn + jitter, pénalité | `configs/fem3d_percussion_cyl_saksala2011.cfg`, T 1e-3 → 2e-4, frames 2 |
| 4 | `fdem3d_kuru9_court` | fdem3d | 6 corps, 3 phases, groupBond, adaptive + potential + DIF yang-fig2 + bulkDamage yang + viscosité Yan graduée | `bench_impact/configs/impact_kuru9.cfg`, T 1e-3 → 1e-5, piston −9 → −30 m/s (contact piston/bit à 6,7 µs ; roche non touchée : import 6 corps + assemblage + contact entre corps), frames 2, `fragBrushV` → `fragBrushV0` (C0) ; **quatre clés du source absentes** (`energyBodyForces = on`, `budgetAbortPct = 2`, `trackGroups = bit piston insert`, `gauge.bit = 0.28 0.30` : métrologie et garde-fou d'énergie, pas la physique — constat du vérificateur, consigné dans l'en-tête ; les remettre = colonnes de `history.csv` en plus = ré-ancrage) |
| 5 | `fdem3d_cut3d_heilman_court` | fdem3d | pdc 3D, **signorini**, adaptive + potential, volume élastique | `configs/cut3d_heilman.cfg`, T 5e-4 → 6e-6, toolX 2,5 → 2,99 mm (4 014 joints insérés, 2 rompus dans la fenêtre), frames 2 |
| 6 | `fdem3d_visc_yan_3d` | fdem3d | joints linéaires + `bulkViscosityXi = 2`, intrinsèque, pénalité | test `visc_yan_3d` de `verify_suite.py`, T 4e-4 → 1,4e-4 et pullV 0,05 → 0,15 m/s (×3 : les 200 joints du plan de Kuhn cassent dans la fenêtre), frames 2 |
| 7 | `fdem_toolcontact_signorini` | fdem | pdc 2D, **signorini**, adaptive + potential (banc T1) | = test `t1_toolcontact_signorini` de `verify_suite.py` (inchangé) |
| 8 | `fdem_ucs_yan_adaptive_court` | fdem | Voronoï/Delaunay (seed), platens, yan, adaptive (fenêtre pré-pic : `ucsStopAfterPeak` présent mais inactif) | `configs_yan/ucs_adap.cfg`, T 1,5e-2 → 2e-3 (pré-pic, 41,7 MPa ; le run source s'arrête seul à 2,87 ms), frames 2 |

Le deck 7 est **le même** que le test `t1_toolcontact_signorini` de la suite (un hachage identique ici
et un PASS là-bas se corroborent ; `toolvb` 1,07878 retrouvé à 4 fils) ; le deck 6 est le test
`visc_yan_3d` de la suite, ×3 en vitesse et raccourci. Les decks 1 et 2 sont les chemins de la calibration cdp et de la matrice v2 ; 3, 4, 5
couvrent saksala2011, l'assemblage multi-corps et la coupe 3D. Les trois decks fdem3d ont été
dimensionnés par la mesure : ce mode coûte 10 à 20 fois le fem3d par tet-pas (kuru9 à T = 5e-5 et
cut3d à T = 6e-5 dépassaient 15 et 30 min), d'où des fenêtres de 6 à 140 µs et non de 50 à 500 µs.

## Règles

1. **Même nombre de fils, ou rien.** Les réductions OpenMP changent l'ordre des sommes flottantes :
   deux hachages ne se comparent qu'au même `OMP_NUM_THREADS`. L'ancre est prise à **4 fils** (la
   machine a 18 cœurs logiques : 14 pour la matrice, 4 pour la porte) ; `bitid.py` rejoue chaque
   deck aux fils de l'ancre. Un deck non déterministe à 4 fils serait listé avec `threads = 1` dans
   `DECKS` et la raison dans `note` — **aucun des 8 n'a eu besoin de cette exception** (rejeu w18 =
   8/8 identiques, voir « Résultats »).
2. **L'ancre est au dépôt** (`tools/bitid_refs.json`, `_meta.exe_sha256` = SHA-256 de l'exe qui l'a
   produite). **Toute ancre changée = `--update` + une ligne dans `CHANGELOG.md`** (« ancre bitid
   changée : <raison> ») **+ commit du json dans le même commit**. Un `--update --only` est une mise
   à jour partielle : `_meta` est conservé, la retouche est journalisée dans `_meta.partial_updates`
   (date, decks, exe), chaque deck porte son propre `exe_sha256` et sa date, et un exe différent de
   `_meta.exe_sha256` marque l'ancre `MIXTE` : même règle CHANGELOG + commit.
3. Une différence n'est **jamais** ré-ancrée en silence : soit c'est un bug (à consigner), soit
   c'est un changement d'algèbre voulu (promotion d'un défaut, réordonnancement flottant), et la
   ligne CHANGELOG le dit.
4. Un exe frais peut être refusé par Apex One (WinError 5 / rc 126) : le script réessaie toutes les
   20 s (10 essais) et écrit `apex_retries` dans le json. La détection par le texte est un motif
   **étroit** (`Permission denied`, `Access is denied`, « Accès (est) refusé » avec accents
   éventuellement remplacés) : un vrai échec de run n'est pas rejoué dix fois pour un « Acc… » fortuit.
5. Les anciens scripts (`etude_lois_fem/bitid_wNN.sh`, `baseline_hash.sh`) et leurs dossiers
   `bitid_wNN/` sont **gelés** : ils documentent la chaîne w0…w19 mais ne servent plus à prouver un
   build (le plan B6 prévoit leur déplacement dans `etude_lois_fem/bitid_archive/`, non fait ici :
   rien n'est déplacé sans accord).

## Empreinte des binaires : `tools/exe_manifest.json`

`python tools/exe_manifest.py` lit chaque `rockim_f2*.exe` de la racine (jamais modifié) et écrit
SHA-256, taille, date, la ligne du tableau « Binaires » de `etude_lois_fem/JOURNAL.md` qui le décrit
(motif `rockim_f2w1*.exe` = lettres seulement, donc w1, w1a…w1d, pas w10…w19) ou « non documente »,
plus les autres lignes du JOURNAL citant le nom complet. État du 05/09 19:12 : 47 exe, 17 décrits
par le tableau (w0, w1*, w2…w8, w11, w15, w16, w19), 30 « non documente » dont **w17, w18, w20**
(cités seulement dans la chronologie : lignes 273, 322, 440) et toute la série f2…f2u, f2p, w9,
w10, w12…w14. Le manifeste se régénère en une commande : à refaire à chaque nouveau binaire. Contrôle croisé : le sha `569d6b9a122207dc` de w16 imprimé dans le JOURNAL est retrouvé ;
`rockim_f2w11_triax.exe` et `rockim_f2w11_gif.exe` sont le **même** binaire (sha identique).

## Résultats du 2026-09-05 (ancre w18, OMP 4, machine chargée : matrice v2 à 14 fils + autres agents)

Ancre : `python tools/bitid.py --exe rockim_f2w18.exe --update` à 18:17-18:34 (passe complète), puis
`--update --only fdem3d_visc` à 18:47 (deck 6 remis à T = 1,4e-4 après un essai à 1,0e-4 sans
rupture ; le json le journalise dans `_meta.partial_updates`). `_meta.exe_sha256 =
655a87329c68b0d88e499a6eebfa73543d5869d9373685ea6430d7c5f1e8c008`, `threads = 4`, 0 reprise Apex One.

| deck | pas (dt) | ancre : durée | hachages (16 hex) history / frames / dernier vtu | pic | rejeu w18 |
|---|---|---|---|---|---|
| fem3d_cdp_PQ_court | 10 501 (1,238e-8) | 129,5 s | e43ca95062fbc58b / 6670a60ed53648e0 / bf3085a23adfad0f | 7 916,75 N | IDENTIQUE (97,7 s) |
| fem3d_dpr_T1_court | 1 968 (1,270e-8) | 221,6 s | 756bc37d5e33f286 / 14d5f7f70b97905d / a091dbf168fb6a84 | 6 775,14 N | IDENTIQUE (184,9 s) |
| fem3d_sk2011_cyl_court | 1 127 (1,775e-7) | 102,3 s | 1ea51ebd4a9da031 / 6f4b7a8006182f30 / da88a5ba839b1d62 | 30 390,2 N | IDENTIQUE (75,2 s) |
| fdem3d_kuru9_court | 2 768 (3,614e-9) | 233,8 s | 1c0d7cba74fe1c86 / b6e78ab13daadbc2 / e44b5aa3d1fd8f85 + joints 56d504dfee1e8a72 | 0 (toolShape none ; contact piston/bit −0,31 J) | IDENTIQUE (157,3 s) |
| fdem3d_cut3d_heilman_court | 3 889 (1,543e-9) | 93,6 s | abcbbdad4ef39a74 / f5604b32880fc9ff / a8922a4e3d7c02e9 + joints 1c8a348e8113179d | 19 541,3 N (4 014 joints insérés, 2 rompus) | IDENTIQUE (69,0 s) |
| fdem3d_visc_yan_3d | 19 576 (7,152e-9) | 207,5 s (269 s sous charge forte) | d5e75bf8a679913f / 37943002f9c24c38 / 3feb409709a875b5 + joints af143536820122cd | — (200 joints rompus, 11,0 mJ visqueux) | IDENTIQUE (226,0 s) |
| fdem_toolcontact_signorini | 205 405 (1,461e-9) | 99,5 s | cca67df9d8aef5fd / 13f14bd526e2f920 / e3de044d79346728 + joints 86d5b2b945f1c460 | 205 788 N/m (toolvb 1,07878 = réf. suite à 1 fil) | IDENTIQUE (77,7 s) |
| fdem_ucs_yan_adaptive_court | 136 574 (1,464e-8) | 180,0 s | abeaff1358e7a7a8 / 8d487de361fd2140 / c9ea1942ca9591f1 + joints 7af7122e40eb65ab | — (41,7 MPa pré-pic) | IDENTIQUE (164,0 s) |

- **Passe complète** : 1 268 s (21 min) pour l'ancre sous charge forte, 1 052 s (17,5 min) pour le
  rejeu ; chaque deck < 5 min à OMP 4 (le plus long : kuru9 234 s, visc 269 s sous la charge la plus
  forte de la soirée).
- **Rejeu w18 = 8/8 IDENTIQUE** (`run_w18_rejeu.log`, sorties `out_rockim_f2w18_omp4_rejeu/`) :
  les 8 decks sont **déterministes à 4 fils** — 24 hachages history/frames/vtu (+ 5 vtu de joints)
  retrouvés à l'octet. Le deck 7 l'a été trois fois (sonde de 17:31, ancre, rejeu), le deck 6 aussi
  (passe complète, ré-ancrage, rejeu). Aucun deck n'a dû être passé à 1 fil.
- **Effet du nombre de fils, mesuré** : `fem3d_cdp_PQ_court` donne un pic de 7 916,75 N à OMP 4
  contre 7 791,04 N à OMP 2 (`bitid_w18/comparaison_w18.txt`, même deck, même exe) et des hachages
  différents (history e43ca950… contre e45a7cbc…). Comparer deux hachages pris à des nombres de fils
  différents ne prouve **rien** — d'où la règle 1.
- **w17 contre l'ancre w18 = 8/8 IDENTIQUE** (`run_w17.log`, 19:00-19:18, 1 097 s, sorties
  `out_rockim_f2w17_omp4/`) : `rockim_f2w17.exe` (11:03, sha `81334bf48cc80989…`) = w18 sans la clé
  `probes` ; les 8 decks, clé absente, rendent les mêmes 24 + 5 hachages et les mêmes pics
  (7 916,75 / 6 775,14 / 30 390,2 / 19 541,3 N, 205 788 N/m) — la preuve « défaut bit-identique » de
  w18 est désormais portée par 8 decks et 3 modes au lieu d'un seul deck cdp à OMP 2.
- Fichiers de ces trois passes : `tests_f2/bitid/run_w18_update.{log,json}`,
  `run_w18_visc_update.log`, `run_w18_rejeu.{log,json}`, `run_w17.{log,json}`.
