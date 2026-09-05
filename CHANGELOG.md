# Changelog de rockim (arbre g0)

Format : [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/). Une section par tag (mesure A3/A5 du
plan de robustesse du 2026-09-05) ; en attendant la mise sous git (A1), la section `[Non publié]`
reçoit les lignes exigées par les règles déjà en vigueur — dont **toute ancre de bit-identité changée**
(`tools/bitid_refs.json`, règle de `tools/BITID.md`).

## [Non publié]

### Ajouté — port de la branche orpheline `insertion-pointe` (décision 3)

Quatre commits de `insertion-pointe` (2ead636), jamais fusionnés, dont **aucune des six clés
n'existait dans `f2`** : contenu réellement orphelin. Deux decks de `f2`
(`bench_impact/configs/impact3d_dpdfh.cfg` et `..._gros.cfg`) les posaient déjà et étaient donc
**refusés** par `rockim_f2w21` (garde C1). Tout arrive en **clé opt-in à défaut bit-identique**.
Détail complet, mesures et bancs : [`docs/PORT_INSERTION_POINTE.md`](docs/PORT_INSERTION_POINTE.md).

- **`insertionTipFactor` (défaut 1) et `insertionTipDamage` (défaut 0,5)** — insertion préférentielle
  en POINTE de fissure, 2D et 3D. Une facette dont un sommet porte déjà un joint inséré et endommagé
  (`D >= insertionTipDamage`) voit son enveloppe DIVISÉE par `insertionTipFactor` ; une facette en
  terrain vierge garde l'enveloppe nominale, donc **l'amorçage — et la résistance macroscopique
  mesurée — restent intacts** (mesuré : pic 51,0807 → 51,0807 MPa). Motivation : l'adaptatif ne
  propage que 43,7 % de ses ruptures contre 56,8 % pour l'intrinsèque, la moyenne sur deux CST
  écrasant la singularité de pointe. À `insertionTipFactor = 1` aucun test supplémentaire n'est
  évalué : **chemin d'origine au bit près**. Gardes : facteur < 1 refusé, `insertionTipDamage` hors
  [0,1] refusé, facteur > 1 sans `insertion = adaptive` refusé. Compte
  `propagations / nucléations` imprimé en fin de run.
- **`insertion = none`** — le CONTINUUM PUR (2D et 3D) : aucun joint n'existe ni ne peut naître, les
  copies de nœuds restent liées pour toujours (= éléments finis à nœuds partagés). Remplace le
  bricolage `insertion = adaptive` + `ft = 1e12`, qui a fait passer un impact 3D DP-DFH de 53 J à
  **−89 GJ en 10 µs** le 2026-08-25 (89 424 joints « inatteignables » activés par un élément
  distordu, portant `dnE = ft/pj = 8 cm`). Trois conséquences câblées : joints liés comme en
  adaptatif, groupes liés intégrant comme UN nœud, et **ressort de pénalité des joints sorti du
  budget de pas de temps** (mesuré : dt × 1,715 sur le banc 2D).
- **`dfhPsiVar` (défaut 0), `dfhPsi0`, `dfhKPsi`, `dfhPsiMax`** — dilatance VARIABLE ψ(p̄) du DP-DFH,
  la forme de `vumat_hole.f` l. 336-339 : `ψ = clamp(dfhPsi0 − dfhKPsi·p̄[MPa], 0, dfhPsiMax)`. La
  constante #5 de la carte (`dfhPsiDeg`, 15 deg) est alors morte. Clé absente ⇒ `psiDeg` fixe,
  trajectoires inchangées au bit près.
- **`dfhD` et `dfhTini` dans les `.vtu`** sous `law = dpdfh` (2D et 3D) : DMAX = max(D1,D2,D3) des
  trois endommagements directionnels du repère figé (SDV 4-6 de la VUMAT, ce que lisent les
  extracteurs du banc 6) et l'instant du premier amorçage. **Ajout de sortie pur.**
- **Decks** `tunnel_edz/configs/tunnel_tip13.cfg`, `tunnel_tip16.cfg`, `tunnel_tip20.cfg` (calibration
  du facteur 1,3 / 1,6 / 2,0 sur le tunnel EDZ), repris tels quels.
- **Bancs courts** `tests_f2/insertion_pointe/` (4 decks + `check_tip.py`, ~10 min à OMP 2) et
  `tests_f2/psivar/` (2 decks point-matériel + `check_psivar.py`, ~2 s). Tous deux impriment un
  verdict OK/ÉCHEC par critère.
- **`etude_lois_fem/meshes/drop_orphans.py` et `check_orphans.py`** : ces outils sont NOMMÉS par le
  message d'erreur du lecteur de maillage mais n'avaient pas été repris à la naissance de `g0`
  (seuls les `*.py` de la racine de `etude_lois_fem/` l'avaient été). Repris de `rockim_f2`.

### Corrigé
- `bench_impact/configs/impact3d_dpdfh.cfg` et `impact3d_dpdfh_gros.cfg` **chargent** (ils étaient
  refusés clé par clé). Trois réparations : les six clés existent maintenant ; la **commande exacte
  de régénération** du maillage est dans l'en-tête (les `.msh` n'existaient nulle part, ni dans `g0`
  ni dans `rockim_f2`, et les comptes annoncés — 98 858 et 143 451 tétraèdres — venaient de
  paramètres jamais consignés : les comptes réellement obtenus, 105 498 et 150 535, sont écrits) ;
  les decks pointent sur `*_clean.msh` car gmsh écrit le point du champ de taille (= le point
  d'impact) en nœud 0-D orphelin que le lecteur refuse.
- Message d'annonce du bloc de liaison des joints : il disait `adaptive insertion: N bonded edges`
  même sous `insertion = none`. Il dit maintenant `insertion = none: N bonded edges` dans ce cas.
  **Sortie seule**, aucun flottant touché.

### Ajouté — port de la branche orpheline `dif-intrinseque` (décision 3)
- `bench_impact/tools/fig_bilan.py` (159 lignes) : la planche du run de référence — partition de
  l'énergie contre les chiffres publiés (ARMA 24-0952 : fissuration 2,6 %, frottement 64,9 % de
  49,3 J) et les sept critères de la Table 3 de Yang et al. rapportés à leur fourchette. C'est,
  avec la coquille de `fig_fp.py`, **tout** ce que `dif-intrinseque` apportait encore : son C++ est
  un sous-ensemble strict de `joint-handoff` (`git diff --stat joint-handoff...dif-intrinseque --
  src include` est vide).
- **Corrigé** `bench_impact/tools/fig_fp.py` : coquille de chaîne non brute `label="moyenne (30 $\mu$s)"` →
  `label=r"..."` (portée de `dif-intrinseque`).

### Bit-identité
- `python tools/bitid.py --exe build/rockim.exe --threads 4` contre l'ancre de naissance
  `tools/bitid_refs.json` : **8/8 IDENTIQUE** après le port (`results/bitid_apres_port_insertion_pointe.json`).
  L'ancre est **inchangée** (aucun `--update`).
- **Trou de l'ancre comblé à la main** : aucun de ses 8 decks ne pose `law = dpdfh`, elle ne pouvait
  donc rien dire du port de ψ(p). Preuve dédiée : sur le point matériel DP-DFH
  `tests_f2/psivar/mp_dpdfh_psivar_off.cfg`, `build/rockim.exe` et le binaire d'avant le portage
  `rockim_f2w21.exe` rendent la **même trace au SHA-256 près**
  (`341fac848845623c2db4a702897f4efab6a797842d843964b6cf83274c997580`).

## [g0-0.1.0] — 2026-09-05 — naissance de g0

### Naissance
- **`g0` naît de `rockim_f2`.** Branche `g0` créée dans le dépôt
  `FDEM/rockim` à partir de `f2-2026-09-02` (commit `4d65395`, la branche la plus proche de l'arbre
  vivant `rockim_f2`), worktree `FDEM/rockim_g0`. L'**état complet** de `rockim_f2` au 2026-09-05 22 h
  (sources `src/` + `include/`, `tools/`, `tests_f2/`, `configs/`, `configs_yan/`, les documents de
  travail dont `PLAN_ROBUSTESSE_2026-09-05.md` et `ETAT_DES_LIEUX_2026-09-05.md`, les `*.py`/`*.md` de
  `etude_lois_fem/`) y est recopié. `diff -r` de `src/` et de `include/` entre `rockim_f2` et
  `rockim_g0` : **0 différence**.
- **Rien n'est supprimé** (décision 2 du plan) : `rockim_f2` n'est ni modifié ni recompilé, les
  worktrees `rockim_p1`…`rockim_p4` restent en place. Ne sont pas repris dans g0 les produits
  (`*.exe`, `*.obj`, `out_*/`, `*.vtu`, gros CSV, figures) ni les dossiers de campagne
  (`etude_lois_fem/matrice*`, `heterogeneite`, `bitid_w*`, `phaseA/B`, `bench_impact/out*`,
  `calib_quick`, `panel_contact*`, `tunnel_schisto`, `g2d`/`g3b`/`gif*`) : ils restent dans
  `rockim_f2`, gelé mais conservé. `.gitignore` réécrit en conséquence.
- **Maillages de la bit-identité versionnés** : les quatre `.msh` dont dépendent les decks
  `tests_f2/bitid/*.cfg` (`meshes/cut3d_h50.msh`, `meshes/impact_kuru_s15.msh`,
  `meshes/t1_toolcontact.msh`, `etude_lois_fem/meshes/T1_c05_clean.msh` — 7,5 Mo au total) sont
  exemptés de la règle `*.msh` du `.gitignore`, sans quoi `tools/bitid.py` ne serait pas rejouable
  depuis un clone. Le cinquième (`Q1_c05_clean.msh`, deck `fem3d_cdp_PQ_court`) reste **hors dépôt**,
  atteint par le chemin relatif inchangé `../../CONTINUUM/calib_bohus_triax/cdp_rockim/perc3d/` :
  `rockim_g0` étant, comme `rockim_f2`, un dossier de `FDEM/`, le chemin des decks n'a pas eu à être
  touché (aucun deck de bit-identité n'a été modifié).

### Ajouté
- **CMake + Ninja** (décision 9, mesure D1a) : `CMakeLists.txt` réécrit pour émettre **exactement** les
  drapeaux de `build_f2.cmd` — `FLAGS = /nologo /EHsc /O2 -std:c++17 -MT -openmp` plus
  `-D_USE_MATH_DEFINES -DNOMINMAX -I include -I ../rockim/eigen-3.4.0` — c'est-à-dire en neutralisant
  les ajouts par défaut de CMake sous MSVC : `/DWIN32 /D_WINDOWS /W3 /GR`, `/Ob2 /DNDEBUG` (qui
  désarmerait les `assert()`) et surtout `/MD` (CRT dynamique) là où `cl` sans drapeau prend `/MT`
  (vérifié : `rockim_f2w21.exe` n'importe que `KERNEL32.dll` et `VCOMP140.DLL`, l'exe CMake aussi).
  `/W4 /WX` viendront au lot D1b, séparément. `tools/build.ps1` enchaîne `vcvars64` → `cmake -G Ninja`
  → `cmake --build`, avec les `cmake.exe`/`ninja.exe` livrés par Visual Studio 2022.
- `build_f2.cmd` est **conservé** un temps comme secours (décision F du plan), il n'est plus la chaîne
  de référence.

### Bit-identité
- Ancre **inchangée** (`tools/bitid_refs.json`, héritée de `rockim_f2`, prise sur `rockim_f2w18.exe`,
  `OMP_NUM_THREADS = 4`, 8 decks). `python tools/bitid.py --exe build/rockim.exe --threads 4` sur
  l'exe produit par CMake + Ninja : **8/8 IDENTIQUE** (2026-09-05 22 h 43 ; exe CMake sha256 `92ef42b6e8168322...` ; rapport
  `results/bitid_g0-0.1.0.json`, journal `results/bitid_g0-0.1.0.log.txt`). La chaîne CMake reproduit donc l'algèbre de la
  chaîne `build_f2.cmd` : g0 naît sur l'ancre de f2, et cette ancre est **la nouvelle ancre de
  naissance** de l'arbre g0 (tag `g0-0.1.0`).


### Ajouté
- 2026-09-05 — `tools/bitid.py` (bit-identité générique, mesure B6), `tests_f2/bitid/` (8 decks
  courts : 3 fem3d, 3 fdem3d, 2 fdem), `tools/exe_manifest.py` + `tools/exe_manifest.json`
  (empreinte SHA-256 des 46 `rockim_f2*.exe` existants), `tools/BITID.md`.
- **ancre bitid changée : première ancre** — prise sur `rockim_f2w18.exe` existant (05/09 13:45,
  `build_f2.cmd`, SHA-256 `655a87329c68b0d88e499a6eebfa73543d5869d9373685ea6430d7c5f1e8c008`),
  `OMP_NUM_THREADS = 4`, 8 decks ; remplace les scripts `etude_lois_fem/bitid_wNN.sh` (ancre changée
  deux fois, hachages non comparables entre nombres de fils). Détail dans `tools/BITID.md`.
  Le dépôt git de `rockim_f2` n'existe pas encore (A1) : le json sera committé avec l'état du 05/09.
- 2026-09-05 (w21, C1/C7) — `include/rockim/KeyGuard.hpp` (audit des clés consommées après `init()` : obsolète,
  autre mode, famille dynamique, suggestion de Levenshtein ≤ 2, inconnue ; politique `unknownKeys = error | warn` ;
  écriture de `config_effective.cfg`) ; `Config::unusedKeys/lineOf/consumedKeys/effective/path/size` ;
  `tools/obsolete_keys.json` (table des clés obsolètes → nom courant, une entrée par renommage, jamais retirée) ;
  registre `tools/gen_keys_by_mode.py` étendu (`kKnown`, `kDynamicPrefix`, `kCommonPrefix`, `kObsolete` dans
  `KeysByMode.hpp` ; `prefixes_dynamic`, `obsolete` dans le json) ; `tools/scan_decks.py` (balayage statique des
  decks selon la règle, `--fix-obsolete`) + rapport `tools/scan_decks_2026-09-05.md` ; banc falsifiant
  `etude_lois_fem/bitid_w21/selftest_cles/` (17 decks, `run_cles.py`, `SELFTEST_cles.md`) ; DOC §8.11.
- **C7** : chaque run écrit `<outputDir>/config_effective.cfg` (clés consommées, valeur deck ou `(defaut)`, tri
  stable) et imprime `[rockim] cles : N consommees (n du deck, m au defaut), k du deck non lues`.
- 2026-09-05 (w22) — `KeysByMode.hpp` : table `kReaders` (lecteurs de chaque clé : solveurs / `shared`), champ
  `readers` du json enrichi (`shared` pour `ALWAYS_COMMON`) ; `Config::deck()`, `Config::seal()` ; banc
  `selftest_cles` : decks `c12`–`c15` (clés à plusieurs lecteurs hors mode, `warn`) et partie C (rejeu de
  `config_effective.cfg`, `history.csv` bit-identique sur les six modes) ; `tools/scan_decks.py` modélise les gardes
  pré-init de `main.cpp` et met à part les fichiers sans clé de solveur (matpoint, cartes matériau).

### Corrigé
- 2026-09-05 (C0) — `fragBrushV` → `fragBrushV0` dans 7 decks de `bench_impact/configs/` (`impact_kuru9*.cfg` ×4,
  `impact_pulv_a/b.cfg`, `impact_stanne_fidele_s15.cfg`) et les 7 mêmes de `FDEM/rockim/rockim_f2_wt/` (originaux
  dans `etude_lois_fem/bitid_w21/orig/`). **Les runs Kuru / pulvérisation / St Anne antérieurs tournaient SANS
  brossage des fragments** (la clé était inerte depuis l'origine) : leurs résultats changent si on les rejoue.
  Restent à corriger (worktrees d'archive, décision Fernando) : `rockim_p1` et `rockim_p2` `bench_impact/configs/impact_pulv_a/b.cfg`, `impact_stanne_fidele_s15.cfg`.

### Modifié
- 2026-09-05 (w21) — **clé inconnue = erreur (décision Fernando 2026-09-05, 20:00)** : toute clé du deck non
  consommée pendant l'initialisation et hors registre (ni commune, ni du mode courant) rend le code 1 avec un message
  qui explique (obsolète → nouveau nom ; autre mode → le mode ; faute de frappe → suggestion ; sinon « inconnue de
  rockim »), toutes les clés fautives listées d'un coup ; `unknownKeys = warn` : mêmes messages en WARNING, le run
  continue. Seul changement de défaut du chantier C. Le contrôle par mode de w20 (`keysbymode::check` avant init) est
  **fondu** dans cet audit (même sous-chaîne de message ; `warn` s'y applique) ; `keysWithPrefix("hydro")` n'est plus
  appelé en mode fdem (il marquerait consommées les clés `hydro*`). Bit-identité : ancre **inchangée** (w18, 4 fils),
  `tools/bitid.py --exe rockim_f2w21.exe --threads 4` : voir `etude_lois_fem/bitid_w21/bitid_w21.log`.
- 2026-09-05 (w22, relecture adverse de w21 — D1, D2) — **règle des lecteurs** : une clé du deck non consommée est
  légitime ssi le mode courant ou le code partagé la lit (`kReaders`) ; sinon `cle 'X' sans effet en mode Y : cle du
  mode Z seulement` / `cles des modes Z1, Z2`. w21 tenait pour « commune = légitime » toute clé lue par ≥ 2 solveurs :
  `jointSoftening`, `insertion`, `bulkDamage`, `fragBrushV0`, `gravity`… passaient en fem3d sans message (D1). 0 deck
  existant touché (957 balayés). **`config_effective.cfg` rejouable** : défauts du code en lignes commentées
  `# cle = valeur (defaut)`, clés du deck non consommées mais légitimes gardées actives ; w21 activait les défauts et le
  rejeu échouait en fem3d (`tensionShearRetention`) et fdem (`dampingLocalAfter`) sur les gardes « satellite orpheline »
  (D2). **Coût** : `Config::seal()` après l'audit — plus de verrou ni de chaîne de défaut dans les getters appelés à
  chaque pas (`confineGaugeTime`). `src/main.cpp` repassé en CRLF (original). Bit-identité : ancre **inchangée** (w18,
  4 fils), `tools/bitid.py --exe <abs>/rockim_f2w22.exe --threads 4` : **8/8 IDENTIQUE** (21:15,
  `etude_lois_fem/bitid_w22/bitid_w22.log`) ; banc `selftest_cles` 47/47 (A 21, C 6, B 20).
- 2026-09-05 (relecture de B5/B6, ancre **inchangée** : mêmes 8 × hachages) — `tools/bitid.py` :
  sorties par défaut dans un dossier temporaire (`tempfile.mkdtemp`, convention de `verify_suite.py`)
  au lieu de `tests_f2/bitid/out_*` ; `.vtu` supprimés aussi sur échec/timeout ; détection Apex One
  par motif étroit (l'ancien « Acc » aurait pris `fragBrushAccel` pour un refus) ; chemins de decks en
  séparateur posix. `tools/bitid_refs.json` : métadonnées seulement (clé orpheline `partial` retirée,
  `cfg` en posix). `tests_f2/bitid/fdem3d_kuru9_court.cfg` : en-tête complété (quatre clés du deck
  source absentes, non documentées avant). `tests_f2/bitid/.gitignore` (`out_*/`).
- 2026-09-05 19:55 — ancre bitid partielle : `fem3d_cdp_PQ_court` rebasé sur `Q1_c05_clean.msh` (w20 refuse le maillage à nœud orphelin) ; history/frames inchangés (e43ca95062fbc58b / 6670a60ed53648e0), vtu e2e42e15c8fbdf2d ; ancre toujours rockim_f2w18.exe, threads 4. Passe complète `tools/bitid.py --exe rockim_f2w20.exe` : 8/8 IDENTIQUE (w20 = gardes C3/C4/C6, sans effet sur les decks valides).
