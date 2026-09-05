# Plan de robustesse de rockim_f2 — 2026-09-05

Destinataire : Fernando, qui décide et arbitre. Rédigé à partir de `ETAT_DES_LIEUX_2026-09-05.md`
(89 constats, trois relectures en lecture seule) et d'une vérification du dépôt le 05/09 au soir
(un seul tag `avant-correctifs-2026-08-19`, cinq branches non fusionnées dans `main`, constitution
présente dans `rockim_f2_wt/.specify/memory/` seulement, `cmake.exe` et `ninja.exe` présents dans
`C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\`
mais hors PATH).

Le plan ne réinvente aucune règle : il rend **mécaniques** les règles déjà écrites dans la constitution
(bit-neutralité, charge nulle, suite, énergie, doc dans le même commit, croissance par addition). Le mode
de panne dominant relevé est le défaut **silencieux** (réglage inerte, garde aveugle, entrée acceptée sans
diagnostic, état corrompu sans erreur) ; chaque chantier vise donc d'abord à convertir un silence en
**erreur nommée** ou en **chiffre verrouillé par une commande**.

Conventions : `R` = `C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_f2` ; `G` =
`C:\Users\fuzquianoalricabi\simulations\FDEM\rockim` (dépôt git et worktrees) ; `CDP` =
`CONTINUUM\calib_bohus_triax\cdp_rockim`. Les identifiants de constats (T-xx tests/build/entrées,
P-xx physique, D-xx documentation/processus) sont attribués dans l'annexe A, dans l'ordre du tableau
de l'état des lieux. Coûts : **h-agent** (session d'agent, code + tests), **h-F** (relecture, décision
ou lancement par Fernando), **h-machine** (calcul, 14 cœurs, hors licence Abaqus). Les décisions qui
changent un défaut de comportement sont marquées **[DÉCISION F]** : elles ne sont pas prises ici.

---

## 0. Diagnostic en une page : quatre racines, un mode de panne

| racine | ce qu'elle produit aujourd'hui | chantier |
|---|---|---|
| Le code vivant est hors de tout dépôt (T01/D01, T12, T13, D11) | 5 400 lignes de diff sans commit, 43 exe à suffixes, 7 arbres, un tag, constitution absente du lieu de travail | **A** |
| La suite ne juge pas le chantier en cours (T03…T09, P17, D08, D09) | aucun test fem3d, aucun run de suite sur w0…w18, références Linux périmées, pas de refs MSVC, ancre de bit-identité changée deux fois, pas de CI | **B** |
| Les entrées n'ont pas de gardes (T02/P20/D02, T14…T18, P11, P13, P14) | clé inconnue ignorée, nœud orphelin devenu broche fantôme (28 runs à refaire), NaN détecté à l'aveugle, unités non vérifiées, `hydro` ignoré en 3D | **C** |
| Le build est manuel et muet (T10, T11, T12) | `/W1` sans `/WX`, lien partiel COMDAT (w1a-w1d invalides), aucune empreinte de version dans le binaire | **D** |

Les chantiers **E** (contact et physique), **F** (combinaisons), **G** (V&V rejouable), **H**
(documentation) et **I** (dépendance) traitent les constats qui restent une fois ces quatre racines
coupées. L'ordre A → B → C → D n'est pas négociable : sans dépôt, aucun critère d'acceptation n'est
rejouable ; sans suite qui juge, aucune garde n'est prouvée ; sans gardes, la matrice de fumée (F)
et le dossier V&V (G) mesureraient des runs faux-mais-plausibles. **Une seule exception, qui passe
avant tout** : la prise d'empreinte de départ (B5 suite `fast` + B6 ancre de bit-identité) sur le binaire
`rockim_f2w18.exe` **existant** (compilé le 05/09 13:45 par `build_f2.cmd`, vérifié présent). Elle ne
demande ni dépôt ni build, et c'est elle qui permet ensuite de prouver que la mise sous git (A1), le
passage à CMake (D1a) et la correction des avertissements (D1b) n'ont rien changé à l'algèbre.

---

## A. Le code vivant retourne sous git et prend un numéro de version

### A.1 Mesures

Le principe est de ne pas créer un nouveau dépôt mais de faire de `R` **le** worktree de `G`, en
conservant le chemin `R` que 92 scripts et les journaux référencent.

| id | constat | mesure (quoi, où) | critère d'acceptation | h-agent | h-F | h-mach. | dépend de |
|---|---|---|---|---|---|---|---|
| A1 | T01/D01 | Archiver l'état hors git : `Rename-Item R R_horsgit_2026-09-05` (lecture seule, conservé jusqu'à A9). Dans `G\rockim_f2_wt` (branche `f2-2026-09-02`) : robocopy de `src include tools configs configs_yan tests_f2 meshes calib_quick etude_lois_fem bench_impact panel_contact_* *.md *.cmd CMakeLists.txt` depuis l'archive, **sans** `*.exe *.obj *.w0 out_* build_*.log g* gif*` ; commit « w18 : état du 05/09 » ; `git worktree move G\rockim_f2_wt R`. | `git -C R rev-parse --is-inside-work-tree` rend `true` ; `diff -rq R\src R_horsgit_2026-09-05\src` et idem `include tools configs tests_f2` ne rendent que les `*.w0` et `*.obj` ; `git -C R status --porcelain` vide. | 2 | 0,5 | 0 | — |
| A2 | T12, D-m6 | `.gitignore` à la racine de `R` : `*.exe *.obj *.pdb *.ilk build/ bin/ out_*/ *.vtu build_*.log *.png g*/ gif*/ __pycache__/ *.w0`. Supprimer de `src/` et `include/` les cinq sauvegardes `.w0` (leçon COMDAT : une copie `.cpp` égarée serait compilée) **[DÉCISION F]** ; déplacer les 47 exe et 126 `out_*` dans `R_horsgit_2026-09-05` (rien n'est supprimé). | `git -C R status --porcelain --ignored` ne liste aucun `*.exe`, `*.obj`, `out_*` comme suivi ; `Get-ChildItem R\src,R\include\rockim -Filter *.w0` rend 0 fichier ; `Get-ChildItem R -Filter *.exe` rend 0 fichier hors `bin/`. | 1 | 0,5 | 0 | A1 |
| A3 | T12 | Fin des suffixes : un seul binaire `R\build\rockim.exe` (D1) ; les binaires « de résultat » sont copiés par le script de porte (B7) sous `R\bin\rockim_<git describe>.exe`. Un tag annoté par binaire ayant produit un résultat cité dans un journal : `v0.18.0` pour w18 (les w0…w17 sont reconstitués comme tags `w00`…`w17` sur les commits que le tableau de JOURNAL.md permet d'identifier, sinon un seul tag `w18-horsgit` sur l'archive). | `git -C R tag` contient `v0.18.0` ; `git -C R describe --tags` rend `v0.18.0` ou `v0.18.0-N-g<sha>[-dirty]` ; `Get-ChildItem R\bin` ne contient que des noms `rockim_v*.exe`. | 3 | 1 | 0 | A1, D1 |
| A4 | D07 | Constitution dans le dépôt de travail : elle y est dès A1 (`.specify/memory/constitution.md` est suivie sur la branche). Amendement v1.2.0 ajoutant le principe IX « toute entrée non reconnue est une erreur » (C1) et X « aucun binaire sans porte » (B7). Suppression des cinq réécritures parallèles (DOCUMENTATION §8, LISEZ_MOI, HANDOFF §6, CHANTIER_f2 §2, guide §2) au profit d'un lien. | `Test-Path R\.specify\memory\constitution.md` rend `True` ; `grep -c "IX\." constitution.md` ≥ 1 ; `grep -l "trois règles\|9 pièges" R\*.md` rend 0 fichier (les sections renvoient à la constitution). | 2 | 1 | 0 | A1 |
| A5 | T12 | `CHANGELOG.md` (format Keep a Changelog, section `[Non publié]` + une section par tag). Hook `pre-commit` (`.githooks/pre-commit`, activé par `git config core.hooksPath .githooks`) : tout commit touchant `src/` ou `include/` sans ligne ajoutée dans `CHANGELOG.md` est refusé. Rétro-remplissage des w0…w18 depuis le tableau de JOURNAL.md (une ligne par build, clé opt-in ajoutée, ancre bit-identité changée ou non). | `git -C R config core.hooksPath` rend `.githooks` ; un commit de test modifiant `src/Config.cpp` sans toucher `CHANGELOG.md` rend un code de sortie ≠ 0 avec le message « CHANGELOG.md non modifié » ; `grep -c "^## \[" CHANGELOG.md` ≥ 19. | 3 | 0,5 | 0 | A1 |
| A6 | T13 | Fusion des branches : (a) trancher les 87 lignes de `joint-handoff` absentes de `f2-2026-09-02` (`git diff joint-handoff f2-2026-09-02 -- src include` relu ligne à ligne, verdict écrit : refactor ou perte) ; (b) `git merge --no-ff f2-2026-09-02` dans `main` après passage de la porte B7 ; (c) `insertion-pointe`, `dif-intrinseque`, `claude/*` (3 distantes + 1 locale, vérifié le 05/09) : soit fusionnées, soit taguées `archive/<branche>` et supprimées **[DÉCISION F]** ; (d) la branche du worktree `studio_wt` (`claude/rockim-last-commit-date-h9nw63`, Rockim Studio — il n'existe **pas** de branche `studio`) reste, renommée `studio` **[DÉCISION F]**. | `git -C R branch --no-merged main` rend au plus la branche de `studio_wt` ; `git -C R log --oneline -1 main` = commit de fusion daté après le passage de la porte ; fichier `docs/FUSION_2026-09.md` listant les 87 lignes et leur sort. | 6 | 2 | 0,5 | B7 |
| A7 | D11 | Un seul arbre de travail : après A6, `git worktree remove` de `rockim_p2`, `rockim_p3`, `rockim_p4` (les dossiers de résultats qu'ils contiennent — `bench_abuaisha`, `tunnel_edz`, `calibration_redbohus` — sont déplacés dans la base `phd_geothermie` au même chemin relatif) ; `rockim_p1` garde le `.git` mais est ramené sur `main` et n'est plus un lieu de travail ; **[DÉCISION F]** pour chaque suppression. | `git -C R worktree list` rend au plus trois lignes (`rockim_p1` sur `main`, `R` sur la branche de travail, `studio_wt` tant que Rockim Studio est un produit séparé) ; `find G -name "*.py" \| xargs -n1 basename \| sort \| uniq -c \| awk '$1>2'` rend 0 ligne (plus de triples copies). | 3 | 2 | 0 | A6 |
| A8 | T13 | Fin des cherry-pick à l'aveugle : règle CONTRIBUTING (H5) « tout cherry-pick ou merge est suivi de `build.ps1` avant push » ; le hook `pre-push` lance `cmake --build build` et refuse le push si la compilation échoue. | Un push après une modification volontairement non compilable rend « pre-push : compilation échouée » et n'atteint pas `origin`. | 1 | 0 | 0,1/push | D1 |
| A9 | T01 | Clôture : la copie hors git est archivée dans la base (`phd_geothermie\FDEM\rockim_f2_horsgit_2026-09-05`, sans exe ni out_*) puis supprimée localement **[DÉCISION F]**, une fois A1…A5 et B7 vérifiés. | `Test-Path R_horsgit_2026-09-05` rend `False` ; la base contient le dossier archivé (commit dans `phd_geothermie`). | 0,5 | 0,5 | 0 | A1-A5, B7 |

### A.2 Ce que A change dans les habitudes

Un binaire n'a plus de nom : il a un `git describe`. Un journal ne cite plus « rockim_f2w16.exe »
mais `v0.16.0` (ou `v0.16.0-3-gabc1234-dirty`, et « dirty » est alors une anomalie à expliquer). Toute
mesure citée dans un rapport doit pouvoir être rejouée par `git checkout <tag>` + `build.ps1` + la
commande du run : c'est le critère G.

Coût du chantier A : **21,5 h-agent, 8 h-F, ~1 h-machine**.

---

## B. La suite de tests devient le juge

### B.1 Mesures

| id | constat | mesure (quoi, où) | critère d'acceptation | h-agent | h-F | h-mach. | dépend de |
|---|---|---|---|---|---|---|---|
| B1 | T07, D-m3 | Câbler dans `tools/verify_suite.py` (liste `TESTS`, tier `fast`) les bancs fem3d existants : `selftest-triax`, `selftest-cdp` (79 contrôles), `selftest-dpdfh` **avec** diff (G7), `selftest-saksala2011` **avec** diff, le selftest des sondes (`probes`, w18), `selftest-hencky` et `selftest-bv` (aujourd'hui `test_kinematics.cpp` hors build → compilé dans `rockim.exe` comme `selftest-kinematics`), Hertz court (E11), T0 fem3d ; plus trois decks fem3d courts avec runs : `configs/verify_fem3d_dp.cfg`, `verify_fem3d_rate1.cfg`, `PQ_cdpI_P020_court` (13 687 tets, 130 µs). Ajout d'une option `--list` qui imprime nom, tier, mode. | `python tools/verify_suite.py --list \| grep -c "mode=fem3d"` ≥ 12 ; `python tools/verify_suite.py --only fem3d --refs refs_msvc.json` rend `12/12`. | 8 | 1 | 1 | B5 (état de départ pris) ; C1 dépend de B1, pas l'inverse |
| B2 | T06, D09 | Références séparées par plateforme : `refs_msvc.json`, `refs_linux.json`, `refs_macos_arm.json` (existant), chacune avec `_meta = {git_describe, exe_sha256, compiler, omp_threads, date}`. Le runner **refuse** de tourner sans `--refs` ou choisit automatiquement le fichier de sa plateforme (`platform.system()` + chaîne compilateur imprimée par `rockim --version`, D3). Les références en dur dans le Python deviennent la référence Linux et ne sont plus modifiées à la main. | `Test-Path R\refs_msvc.json` et `refs_linux.json` rendent `True` ; `python tools/verify_suite.py --tier fast` sans `--refs` sur MSVC choisit `refs_msvc.json` et l'imprime en première ligne ; `git log --oneline -- tools/verify_suite.py` ne montre plus de commit « recale référence » après la date de B2 (les recalages passent par `--update-refs` + commit du json). | 6 | 1 | 2 | B1, D3 |
| B3 | T05, P08 | Certification à un fil. Fait vérifié le 05/09 : le runner écrase **déjà** l'environnement avec `OMP_NUM_THREADS = args.threads` (défaut `"1"`, `verify_suite.py:1162,1182`) — le « 35/40 » du 26/08 vient donc d'un `--threads 4` **explicite**, pas de l'environnement. Mesure : un appel avec `--threads N ≠ 1` est refusé **sauf** `--threads N --uncertified` explicite ; le résultat porte alors `"threads": N, "certified": false` et le runner n'imprime pas de compteur PASS/FAIL contre les refs (comparaison interdite) ; le nombre de fils est écrit dans chaque `results/*.json` et en tête du log. Rejeu à 1 fil des quatre tests du « 35/40 » (`fdem_voronoi_tension`, `jointdeath`, `neohooke`, `jointquad`) et correction de la note mémoire `reference-rockim-worktrees-suite`. | `python tools/verify_suite.py --tier fast --threads 4` rend le code 2 et « la suite certifie à 1 fil (ajouter --uncertified) » ; avec `--uncertified`, le json porte `"certified": false` et aucune ligne `N/N` ; `results/suite_<date>_<describe>.json` à 1 fil contient `"threads": 1` et les 4 tests en PASS ou leur écart nommé. | 2 | 0,5 | 1 | — |
| B4 | T04 | Tier `full` sous MSVC sur le binaire courant, une fois, à 1 fil ; les 8 tests à comptages entiers chaotiques (`ucs_yan_*`, `percussion_2d*`, `shpb_mini*`) reçoivent leur valeur MSVC dans `refs_msvc.json` (mécanisme `--update-refs`), et un champ `chaotic: true` qui fait comparer au médian plateforme plutôt qu'à l'entier Linux. | `python tools/verify_suite.py --tier full --refs refs_msvc.json` rend `N/N` (N ≥ 97) ; fichier `results/suite_full_<date>_<describe>.json` au dépôt. | 4 | 1 | 4 | B2, B3 |
| B5 | T03, P17 | **Première action du plan, avant A1.** Suite `fast` rejouée à 1 fil sur le binaire **existant** `R\rockim_f2w18.exe` (05/09 13:45, `build_f2.cmd`), sorties redirigées dans `results/suite_fast_2026-09-05_w18.txt` (+ `--json`), **avant** toute modification de source et avant tout build : c'est l'état de départ. Rejouée ensuite sur l'exe CMake de D1a (doit être identique) puis à chaque tag. | `results/suite_fast_2026-09-05_w18.json` existe et rend `47/47` (ou liste nommée des écarts, chacun traité dans B4 ou CHANGELOG) ; le même fichier pour l'exe D1a a les mêmes verdicts et les mêmes valeurs extraites. | 1 | 0,5 | 0,5 | — (aucune) |
| B6 | T08, D08 | Bit-identité générique : `tools/bitid.py --exe <exe> [--update]` lit `tools/bitid_refs.json` = 3 decks fem3d (`PQ_cdpI_P020_court` cdp+signorini+quarter ; `C_T1_R_P100` dpr+briques+hencky+bv ; `xval_fem3d_saksala2011` court) + 3 decks fdem3d (`impact_kuru9` court adaptive+potential+DIF+groupBond ; `cut3d_heilman` court signorini ; `visc_yan_3d`) + 2 decks fdem 2D (`verify_fdem_toolcontact`, `ucs_yan_adaptive` court), chacun ≤ 5 min à OMP 4 (la passe complète ≤ 30 min ; `PQ_cdpI_P020_court` tourne déjà en minutes, pas en secondes : les decks « court » sont raccourcis en T si besoin), `OMP_NUM_THREADS=4` fixé dans le json, hachage SHA-256 de `history.csv`, du dernier `.vtu`, de `probes.csv` ; boucle de reprise Apex One intégrée ; suppression des `.vtu` après hachage. **L'ancre initiale est prise sur `R\rockim_f2w18.exe` existant, avant A1 et avant D1** (deuxième action du plan, avec B5) ; D1a doit la reproduire à l'identique. L'ancre est **au dépôt** ; la changer = `--update` + commit + ligne CHANGELOG « ancre changée : raison ». Les scripts `bitid_wNN.sh` et `baseline_hash.sh` sont gelés (dossier `etude_lois_fem/bitid_archive/`). | `python tools/bitid.py --exe rockim_f2w18.exe --update` écrit `tools/bitid_refs.json` (8 entrées, `_meta.exe_sha256` = SHA-256 de `rockim_f2w18.exe`) ; `python tools/bitid.py --exe build\rockim.exe` rend `8/8 IDENTIQUE` sur le commit de l'ancre ; après une modification volontaire d'une constante physique, rend `7/8` et nomme le deck ; `git log -p -- tools/bitid_refs.json` montre pour chaque changement une ligne CHANGELOG correspondante (vérifié par le hook A5 étendu). | 6 | 1 | 0,5/run | — (script Python + exe existant) |
| B7 | T09 | Porte locale obligatoire : `tools/gate.ps1` = `build.ps1` (D1) → `verify_suite.py --tier fast` → `bitid.py` → écrit `results/gate_<describe>.json {ok, suite, bitid, smoke, exe_sha256, apex_retries}` → copie l'exe sous `bin/rockim_<describe>.exe` **seulement si ok**. Porte **v1** (vague 1) = build + suite + bitid ; l'étape `smoke_matrix.py --quick` (F1, sous-ensemble 5 min) s'y ajoute en vague 2 dès que F1 existe (`smoke: null` avant). Les lanceurs de campagne (`run_*.ps1/.sh`, `queue_nuit.sh`) refusent un exe sans `gate_<describe>.json` ok. Durée d'une porte ≈ 1 h (build 4 min + suite `fast` 22 min mesurées + bitid ≤ 30 min + fumée 5 min). | `tools/gate.ps1` rend le code 0 et `bin/rockim_v0.18.0.exe` ; après une modification cassant un test, rend ≠ 0 et **aucun** nouvel exe dans `bin/` ; `bash queue_nuit.sh` avec un exe non porté rend « exe sans porte » et ne lance rien. | 4 | 0,5 | 1/porte | B1-B6, D1 (F1 pour l'étape fumée) |
| B8 | T09, T-m1 | CI GitHub Actions `.github/workflows/ci.yml` sur `ubuntu-latest` : cmake `-Wall -Wextra -Werror`, `verify_suite.py --tier fast --refs refs_linux.json`, `bitid.py --refs tools/bitid_refs_linux.json` ; les étapes `gen_keys.py --check` (H1) et `check_stepmap.py` (H4) sont ajoutées au workflow **quand H1/H4 sont livrés** (vague 3), pas avant ; artefact `results/*.json` téléchargé et committé par Fernando dans `results/ci/`. Nightly hebdomadaire : tier `full` + `vv/run_all.py --tier short` (G). | `gh run list --limit 1 --json conclusion` rend `success` sur `main` ; badge dans README ; `results/ci/` contient un json par semaine ; après H1/H4, `.github/workflows/ci.yml` contient les deux étapes (grep). | 5 | 1 | 0 (GitHub) | A6, D1 (H1, H4 pour les étapes doc) |
| B9 | T09 | Le runner scanne stdout/stderr : `nan`, `inf`, `[rockim] error`, `warning`, `WARN` → FAIL avec la ligne fautive, sauf motif attendu déclaré dans le test (`expect_warn=…`). Contrôle négatif : `tests_f2/tg_nan.cfg` (dtFactor 5) attendu en FAIL. | `python tools/verify_suite.py --only tg_nan` rend `FAIL (attendu) : NaN détecté` ; un test qui imprime `[WARN]` non déclaré rend FAIL. | 2 | 0 | 0,1 | C4 |
| B10 | D-m3, P-info thermique | Les 47 decks `tests_f2/t1…t18, tg_*` (comptés le 05/09 ; l'état des lieux dit 44, 47 ou 48 selon la section) reçoivent un runner : entrées `tier="f2"` dans `TESTS` avec au minimum code 0 + pas de NaN + résidu d'énergie < 5 % + un repère chiffré tiré des `.md` qui les ont validés (t9b : 0,011 % ; t14c : 1,5910/4,3307 ; t3/t4/t8/tg_* : contrôles négatifs attendus en erreur). `_eapp.py` corrigé (lit `out_t14c_*`, seuil 1 %). | `python tools/verify_suite.py --tier f2 --refs refs_msvc.json` rend `47/47` (= `ls tests_f2/*.cfg | wc -l`) ; `results/suite_f2_*.json` au dépôt. | 8 | 1 | 3 | B2, C1 |
| B11 | P13 | Hydro entre dans la suite : deck court AbuAisha (`e1` réduit à T/10) avec `hydro_sign_check.py` intégré comme extracteur (signe de l'ouverture, pression cavité) et résidu d'énergie. | `python tools/verify_suite.py --only hydro_sign` rend PASS avec `ouverture > 0` et `residu < 10 %` (valeur réelle 2,5-10,5 % consignée, pas 0,0125 %). | 3 | 0,5 | 0,5 | B10 |

### B.2 Règle de non-régression sous Windows, écrite une fois

Sur MSVC, un test à comptage entier chaotique est jugé contre `refs_msvc.json` (`chaotic: true`, médian
± tolérance plateforme) ; un test de mécanisme (charge nulle, dampWork ≤ 0, T0/T1, selftests) est jugé
à l'identique sur toutes les plateformes. Le critère « comparaison avant/après sur la même machine »
devient donc **le** fichier de références de la machine, versionné, et non plus une lecture humaine
de deux logs.

Coût du chantier B : **49 h-agent, 8 h-F, ~13 h-machine + ≈ 1 h par porte**.

---

## C. Les entrées ont des gardes : tout silence devient une erreur nommée

### C.1 Mesures

| id | constat | mesure (quoi, où) | critère d'acceptation | h-agent | h-F | h-mach. | dépend de |
|---|---|---|---|---|---|---|---|
| C0 | T14 | Corriger `fragBrushV` → `fragBrushV0` dans les 10 decks (`bench_impact/configs/impact_kuru9*.cfg`, `impact_pulv_a/b.cfg`, et ceux de `rockim_p1`), **avant** C1 qui les refuserait. Relire chaque valeur : le brossage était inerte depuis l'origine, le résultat des runs Kuru change. | `grep -rE "^\s*fragBrushV\s*=" R G --include=*.cfg` rend 0 ligne ; CHANGELOG note que les runs Kuru antérieurs tournaient sans brossage. | 0,5 | 0,5 | 0 | — |
| C1 | T02/P20/D02 | `Config` suit les clés consommées : `kv_` devient `map<string, Entry{value, line, consumed}>`, chaque `get*/req*/has` marque `consumed` (mutable). `Config::unusedKeys()` rend la liste ; `main.cpp` l'appelle après `solver->init()` et lève `runtime_error` « clé inconnue 'bulkViscosty' (ligne 12) — vouliez-vous dire 'bulkViscosity' ? » (distance de Levenshtein ≤ 2 contre **l'ensemble des clés consommées pendant `init()`**, que `Config` connaît désormais ; **aucun registre externe n'est requis** — le registre H1, quand il existera, n'améliore que le message : « clé de mode fdem seulement » plutôt que « clé lue par aucun getter ») ; une clé lue **conditionnellement** (ex. `hydroStart` sans `hydro = on`, `jointResidualMu` sous `jointSoftening = linear`) n'est pas consommée et est donc refusée aussi : c'est voulu, c'est précisément le « réglage inerte » du lot E ; appel répété à `finalize()` pour signaler une clé lue **après** init (cas `gravity` lu dans `placeTool`, E16) comme avertissement. Échappatoire opt-in `unknownKeys = warn` **[DÉCISION F : le défaut devient l'erreur, c'est le seul changement de défaut de ce chantier]**. `getb` refuse toute chaîne hors `1/0/true/false/yes/no/on/off` ; une ligne sans `=` est une erreur avec numéro de ligne ; clé répétée = dernière gagne **avec** message « clé X redéfinie ligne a → b ». | `rockim.exe tests_f2/tg_unknownkey.cfg` (clé `kinematic = hencky`) rend le code 1 et une ligne contenant `vouliez-vous dire 'kinematics'` ; `tests_f2/tg_getb.cfg` (`cutterFloor = flat`) rend le code 1 ; la suite `fast` (B1) et la bit-identité (B6) passent sans qu'aucun deck ne déclenche la garde, et `rockim.exe --dump-keys` (C2) rend 0 clé refusée sur les decks des études actives (`matrice_v2/`, `CDP/perc3d/`, `configs/`) — chaque deck corrigé est listé dans CHANGELOG ; les decks de B4/B10 sont corrigés quand ces tiers sont câblés. | 6 | 1 | 1 | C0, B1, B6 (preuve de non-casse) ; H1 facultatif (message) |
| C2 | T02 | Traçabilité des clés : chaque run écrit `out/config_effective.cfg` (toutes les clés consommées, valeur effective, `# défaut` ou `# deck ligne n`) et `out/run_info.json` (D3) ; option `--dump-keys` sans lancer le calcul. | Tout `out_*/` produit après C2 contient `config_effective.cfg` dont le nombre de lignes = nombre de clés consommées imprimé en fin d'init ; `rockim.exe --dump-keys deck.cfg` rend la même liste et le code 0. | 2 | 0 | 0 | C1 |
| C3 | T16, T17, P11 | Import de maillage : `Fem3dSolver::buildMeshFile` et `Fdem3dSolver::buildFromTets` (et le lecteur 2D) — nœud jamais référencé par un élément → erreur « nœud orphelin id=1234 (24, 24, 32) : point de champ de taille Gmsh ? utiliser meshOrphans = drop » ; opt-in `meshOrphans = drop` supprime et compte ; « degenerate tet » nomme l'élément, ses 4 nœuds, leurs coordonnées et V0 ; qualité : dièdre min et rapport d'aspect imprimés, seuil `meshMinDihedralDeg` (défaut 0 = avertissement seul) ; `finishMesh` de **tous** les solveurs : tout nœud à `m ≤ 0` après assemblage → erreur nommée (plus jamais épinglé FIXED en silence). `drop_orphans.py`/`check_orphans.py` deviennent inutiles (conservés, marqués obsolètes). | `rockim.exe tests_f2/tg_orphan.cfg` en `mode = fem3d` **et** `fdem3d` rend le code 1 avec `nœud orphelin id=` ; le même deck avec `meshOrphans = drop` rend le code 0 et « 1 nœud orphelin supprimé » ; un maillage à tet dégénéré rend `élément id=… V0=…` ; bitid B6 inchangé (les maillages `*_clean` n'ont pas d'orphelin). | 5 | 0,5 | 0,2 | — |
| C4 | T18 | Détecteur de NaN réel et commun : `include/rockim/Guards.hpp` — `checkFinite(step, v, ke, work)` scanne **toutes** les vitesses nodales (réduction OpenMP, coût < 1 % à `nanCheckEvery = 1024`) et les énergies, nomme le premier nœud/élément non fini, écrit une dernière trame + `out/ERROR.txt`, code de sortie 2 ; appelé dans `FemSolver`, `Fem3dSolver` (remplace `u_[0]`), `DemSolver`, `Dem3dSolver`, `FdemSolver`, `Fdem3dSolver`. Piège FP opt-in `fpTrap = on` (`_controlfp_s` MSVC / `feenableexcept` gcc) pour les bancs. | `tests_f2/tg_nan_<mode>.cfg` pour les 6 modes rendent le code 2 et « NaN détecté au pas N, nœud i » ; avec `fpTrap = on` le pas d'arrêt N' vérifie N − 1024 < N' ≤ N ; `grep -c "u_\[0\]" src/Fem3dSolver.cpp` rend 0 dans le contexte garde ; `grep -l checkFinite src/*Solver.cpp` rend 6 fichiers. | 5 | 0,5 | 0,5 | — |
| C5 | T15 | Plausibilité dimensionnelle (pas un système d'unités : une carte de bornes SI) : après init, chaque solveur passe ses constantes à `Guards::plausible(name, value, lo, hi, unit)` — E ∈ [1e9, 5e11] Pa, ρ ∈ [500, 2e4] kg/m³, ft ∈ [1e5, 1e8] Pa, Gf ∈ [0,5, 1e4] J/m², c_d calculé ∈ [1e3, 1e4] m/s, vitesses imposées < 1e3 m/s, pressions < 2e9 Pa, dt ≤ dt_CFL, masse d'outil cohérente avec `quarterModel` (m/4 attendu si `toolMassQuarter` absent → avertissement), `dampingLocal` par défaut en scénario dynamique → avertissement (T21) ; dépend de C2 (`--dump-keys`). Hors bornes → erreur ; `unitsCheck = warn` opt-in. Le tableau des bornes est un fichier `tools/plausibility.json` versionné, qui déclare aussi les **valeurs sentinelles légitimes** (`ft = 0` = joint incassable en GBM, `Gf = 0`, `contactMu = 0` des bancs…) pour ne pas les refuser. **Avant activation**, les bornes sont passées sur les 954 decks recensés via `rockim.exe --dump-keys --plausibility` (C2) et chaque deck refusé est soit corrigé, soit la borne est élargie et justifiée dans le json. | `tests_f2/tg_units.cfg` (`E = 77.66`) rend le code 1 et « E = 77.66 Pa hors [1e9, 5e11] : GPa saisis en Pa ? » ; la suite entière passe sans déclencher C5 ; `results/plausibility_scan_<date>.txt` liste les 954 decks avec 0 refus non justifié ; `grep -L dampingLocal calib_quick/*.cfg CDP/perc3d/PQ*.cfg` rend 0 fichier. | 4 | 1 | 0,5 | — |
| C6 | P13 | `hydro` refusé hors `fdem` : garde **nommée** dans `main.cpp` à côté de `thermal`/`beddingDip` (mêmes lignes 147-160, vérifié : aucune occurrence de `hydro` dans main.cpp le 05/09), posée **dès la vague 1** indépendamment de C1 ; avec C1, toute clé qu'aucun getter du mode ne lit (`hydro`, `weakPlanes`, `preBrokenJoints`, `microseismic`, `grainSizeSpread`, `grainMeshRandom`, `historyStrains` en 3D ; `quarterModel`, `kinematics`, `bulkViscosity b1 b2` en fdem3d) est **déjà** refusée comme « clé lue par aucun getter » ; le registre H1 (vague 3) ne fait que remplacer ce message par « clé de mode fdem seulement ». | `rockim.exe tests_f2/tg_hydro3d.cfg` rend le code 1 et « 'hydro' n'est implémenté que pour mode = fdem » ; `rockim.exe` avec `quarterModel = true` en `fdem3d` rend le code 1 (message C1, puis message de mode après H1). | 2 | 0 | 0 | — pour la garde `hydro` ; C1 pour la généralisation ; H1 pour le message |
| C7 | P14 | Fin de la collision `bulkViscosity` : en `fdem`/`fdem3d` la clé devient `viscosityMunjiza` (μ, Pa·s ; `viscosityMunjizaXi`, `…Graded`, `viscousInInsertion` renommés de même) ; l'ancien nom en mode fdem est une **erreur** « clé renommée viscosityMunjiza ; bulkViscosity = b1 b2 désigne la viscosité de volume Abaqus (fem3d) » ; `bulkViscosity = b1 b2` reste fem3d. Decks migrés par script (`tools/migrate_keys.py`, table de renommage versionnée). | `grep -rlE "^\s*bulkViscosity" R --include=*.cfg` ne rend que des decks `mode = fem3d` ; tests `visc_yan_3d` passent avec la nouvelle clé et bitid B6 identique (le renommage ne change pas l'algèbre) ; doc H1 montre deux entrées distinctes. | 3 | 0,5 | 0,3 | C1 |
| C8 | T-m4, T18 | Arrêt propre : codes de sortie distincts (1 configuration/maillage, 2 NaN, 3 budget d'énergie, 4 licence/IO), `out/ERROR.txt` avec le message, la dernière trame écrite et `history.csv` fermé proprement ; le runner B9 lit `ERROR.txt`. | Chaque `tg_*.cfg` de contrôle négatif rend le code attendu (table dans `tests_f2/README.md`) et un `ERROR.txt`. | 2 | 0 | 0 | C4 |

### C.2 Ordre impératif à l'intérieur de C

C0 avant C1 (sinon les decks Kuru sont refusés à cause d'une coquille qu'on aurait pu corriger en
silence) ; C1 avant C6 et C7 (le registre porte les modes et les renommages) ; C3 et C4 indépendants
et lançables en parallèle de A.

Coût du chantier C : **29,5 h-agent, 4 h-F, ~2,5 h-machine**.

---

## D. Le build est reproductible, bruyant et daté

### D.1 Mesures

| id | constat | mesure (quoi, où) | critère d'acceptation | h-agent | h-F | h-mach. | dépend de |
|---|---|---|---|---|---|---|---|
| D1a | T10 | CMake + Ninja sur cette machine avec les exécutables livrés par Visual Studio (`cmake.exe` et `ninja.exe` vérifiés présents sous `…\Common7\IDE\CommonExtensions\Microsoft\CMake\`) : `tools/build.ps1` = `vcvars64` puis `cmake -S R -B R\build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=<ninja VS>` puis `cmake --build build --parallel` ; **flags strictement identiques à `build_f2.cmd`** (`/std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX`, vérifié l. 8) pour que l'exe CMake reproduise l'ancre B6 prise sur `rockim_f2w18.exe` ; si les hachages diffèrent malgré des flags identiques (ordre des .obj, `/Zi`), la cause est cherchée jusqu'à identité **ou** l'ancre est redéfinie sur l'exe CMake avec la différence consignée. `build_f2.cmd` conservé un mois comme secours puis supprimé **[DÉCISION F]**. | `tools/build.ps1` rend le code 0 et `build\rockim.exe` ; `python tools/bitid.py --exe build\rockim.exe` rend `8/8 IDENTIQUE` contre l'ancre w18 ; `docs/BUILD.md` consigne le verdict (identique, ou ancre redéfinie et pourquoi). | 4 | 0,5 | 0,5 | B6 (ancre), A1 |
| D1b | T10, T-m2 | Build bruyant : `CMakeLists.txt` passe à MSVC `/W4 /permissive- /Zi /fp:precise /openmp` avec Eigen en en-tête **externe** (`/external:I ..\rockim\eigen-3.4.0 /external:W0`, sinon `/WX` est inatteignable à cause des avertissements d'Eigen), gcc `-O3 -Wall -Wextra` ; correction de tous les avertissements `/W4` du code rockim (estimation 200-600 sur 31 k lignes ; les conversions sont corrigées, pas masquées, sauf liste explicite `#pragma warning` justifiée dans CONTRIBUTING) **par lots, chaque lot vérifié bit-identique par B6** ; `/WX` et `-Werror` activés seulement quand le compte est à 0. | `Select-String -Path build\build.log -Pattern "warning C" \| Measure-Object` rend 0 ; `grep -c "/WX" CMakeLists.txt` ≥ 1 ; `python tools/bitid.py --exe build\rockim.exe` rend `8/8 IDENTIQUE` après le dernier lot (une correction qui change l'algèbre — ex. conversion `int`→`double` mal placée — est un bug trouvé, à consigner, pas à ré-ancrer en silence). | 10 | 0,5 | 0,5 | D1a |
| D2 | T11 | Suivi de dépendances garanti par Ninja (fichiers `.d` de `cl /showIncludes`) : plus aucun `.obj` dans `src/` ; règle « rebuild complet après `.hpp` » supprimée de la mémoire et du journal, remplacée par « `build.ps1` suffit ». Contrôle : toucher `include/rockim/MatLaw.hpp` recompile exactement les `.cpp` qui l'incluent. | `ninja -C build -t deps` liste chaque `.obj` avec ses en-têtes ; après `touch include/rockim/MatLaw.hpp`, `ninja -C build -n \| grep -c "\.obj"` = `grep -l "MatLaw.hpp" src/*.cpp \| wc -l` (transitivement) ; `Get-ChildItem R -Filter *.obj -Recurse` hors `build/` rend 0. | 1 | 0 | 0,2 | D1 |
| D3 | T12 | Empreinte de version : CMake génère `include/rockim/Version.hpp` à la configuration (`git describe --tags --dirty --always`, date ISO, compilateur et version, `ROCKIM_OMP`) ; `rockim --version` l'imprime ; chaque run imprime en tête et écrit `out/run_info.json` : version, SHA-256 de l'exe (calculé sur `argv[0]` au démarrage), SHA-256 du deck, `OMP_NUM_THREADS` effectif, hôte, date, durée. | `build\rockim.exe --version` rend une ligne conforme à `^rockim v\d+\.\d+\.\d+(-\d+-g[0-9a-f]+)?(-dirty)? \d{4}-\d{2}-\d{2}T\S+ MSVC 19\.\d+ omp=\d+$` ; tout `out_*/run_info.json` a les 8 clés ; un binaire `-dirty` est refusé par la porte B7. | 3 | 0 | 0 | D1 |
| D4 | T-m3 | Apex One : Fernando demande l'exclusion de `R\build` et `R\bin` (action DSI, hors agent) ; en attendant, `gate.ps1` et `bitid.py` intègrent la boucle de reprise (20 s, 10 essais) et journalisent le nombre d'essais. | `results/gate_*.json` contient `"apex_retries": n` ; après exclusion, n = 0 sur cinq portes consécutives. | 0,5 | 0,5 | 0 | B7 |
| D5 | T10 | Parité gcc : la CI B8 compile avec `-Werror` ; localement, WSL ou le conteneur ne sont pas requis (la CI suffit). | `gh run list` vert sur le commit qui active `-Werror`. | — | — | — | B8 |

Dans le reste du plan, une dépendance notée « D1 » désigne **D1a** (le build CMake `tools/build.ps1` existe et
reproduit l'ancre) ; D1b (`/W4 /WX`) n'est prérequis de rien d'autre que du critère final n° 2.

Coût du chantier D : **18,5 h-agent, 1,5 h-F, ~1,5 h-machine** (D1a 4 + D1b 10 + D2 1 + D3 3 + D4 0,5).

---

## E. Contact et physique : chaque mécanisme a un banc qui peut échouer

### E.1 Mesures

| id | constat | mesure (quoi, où) | critère d'acceptation | h-agent | h-F | h-mach. | dépend de |
|---|---|---|---|---|---|---|---|
| E1 | P04 | Banc de frottement à deux corps, 2D et 3D, pénalité **et** potentiel : bloc élastique lancé à v0 sur un socle fixe, μ = 0,3, `dampingLocal = 0` ; solution : décélération μg, distance d = v0²/(2μg), énergie dissipée μ m g d. Entrées `friction_{penalty,potential}_{2d,3d}` tier `fast` (≤ 30 s chacune). | Les 4 entrées rendent PASS avec `mu_mesure = 0,30 ± 0,02` (pente de v(t)) et `eFric = μ m g d ± 2 %` ; le même banc avec `contactMu = 0` rend `d` ≥ 10× (contrôle qui doit « échouer » à dissiper). | 8 | 1 | 0,5 | C4, E3 |
| E2 | P04 | Banc d'empilement : N = 5 blocs sous gravité, 2D et 3D, pénalité et potentiel ; équilibre : force au socle = N m g, énergie cinétique → 0, pénétration bornée et **stationnaire** sur 1e5 pas (pas de dérive). | Entrées `stack_{penalty,potential}_{2d,3d}` : `F_socle = N m g ± 1 %`, `KE_fin < 1e-6 KE_max`, `pen(1e5) − pen(5e4) < 1e-3 pen(5e4)`. | 6 | 0,5 | 0,5 | E1 |
| E3 | P05, P07 | Bilan d'énergie fermable sur pièces : export dans `history.csv` (colonnes ajoutées **en fin**) de `toolWork`, `biasW`, `bcWork`, `KE_roche`, `W_gravité` ; `energyBodyForces` compte la pesanteur dès qu'un `gravity` est posé (le défaut `off` reste, mais un deck avec gravité et `off` reçoit un avertissement) ; script `tools/energy_check.py out_dir --tol` qui ferme le bilan à partir des colonnes et rend PASS/FAIL. Comparaison Munjiza (restitution potentiel) et Solidity (frottement E1) documentée dans `vv/contact/README.md`. | `python tools/energy_check.py out_friction_penalty_2d --tol 1` rend PASS ; sur T1b (`verify_fdem_toolcut_break`) rend un résidu < 1 % de l'apport ; `head -1 out_*/history.csv` contient `toolWork,biasW,bcWork,keRock`. | 6 | 1 | 1 | — |
| E4 | P02, P07, T25 | Indicateurs de pompe génériques dans fem3d et fdem3d : pour tout corps à vitesse imposée (outil, lame, mors) — `toolinj` = travail injecté / travail de corps rigide de référence, `toolvb` = vitesse nodale max / vitesse d'outil — calculés comme en 2D (`FdemSolver.cpp`, extracteurs `verify_suite.py:196-219`), imprimés et exportés ; entrées T1 fem3d et T1 fdem3d dans la suite, chacune avec sa variante pénalité **qui doit échouer**. | `grep -c toolinj src/Fem3dSolver.cpp src/Fdem3dSolver.cpp` ≥ 1 chacun ; `verify_suite.py --only toolcontact` rend 6 entrées PASS dont `t1_fem3d_penalty` et `t1_fdem3d_penalty` en « FAIL attendu (toolinj > 1,5) ». | 6 | 0,5 | 0,5 | E3 |
| E5 | P01, T-m5 | Rejeu du cas explosif 3D : percussion 3D homogène (grille de Kuhn), T = 4e-4 s, binaire courant (borne MOOSE, CFL réel, B4 par sous-système, potentiel) ; critères écrits **avant** le run : résidu < 2 % de KE0, vitesse nodale max < 2× vitesse d'outil, aucun NaN, cinétique des débris < 5 % de KE0 à T. Le run devient `percussion_3d_debris` tier `full` (T réduit à 3e-4 si > 1 h) ; s'il échoue, le chantier D0 du ROADMAP redevient bloquant et **aucune** étude fdem3d en phase débris n'est lancée. Validation de rebond 3D (étape E3 du plan de fiabilité) par deux blocs. | `results/suite_full_*.json` contient `percussion_3d_debris: PASS` avec les 4 chiffres ; `vv/explosif3d/VERDICT.md` daté ; `rebound_3d` tier `full` : restitution `e = 1,00 ± 0,01` (élastique, potentiel). | 6 | 2 | 6 | C4, E3 |
| E6 | P11, T17 | Garde broche fantôme dans `Fdem3dSolver::toolContact` (miroir de `if (m_[i] <= 0.0) continue;`) **et** en amont : C3 refuse tout nœud à masse nulle à l'init, donc la garde de contact ne peut plus s'exercer que sur un nœud rendu massless en cours de run (érosion) — cas alors journalisé. | `rockim.exe tests_f2/tg_orphan.cfg` en `fdem3d` rend le code 1 à l'init ; `grep -c "m_\[i\] <= 0" src/Fdem3dSolver.cpp` ≥ 1. | 1 | 0 | 0 | C3 |
| E7 | P09 | T0 sur la copie fdem3d : `Fdem3dSolver::nodeSig` (lignes 4064-4109) appelle `toolsig::impulse` (étape 9 du plan panel) ; `selftest-toolcontact3d` couvre le chemin fdem3d avec les 23 contrôles de T0 en 3D ; bitid sur les 7 decks signorini fdem3d avant/après. | `rockim.exe selftest-toolcontact3d` rend `23/23` avec écart max ≤ 1e-14 ; `grep -c "toolsig::" src/Fdem3dSolver.cpp` ≥ 1 ; bitid des 3 decks fdem3d de B6 identique ou différence expliquée (ordre des opérations flottantes) et ancre changée avec CHANGELOG. | 4 | 0,5 | 0,3 | B6 |
| E8 | P02 | Lame fem3d en Signorini **ou** verdict écrit : (a) implémenter la lame comme demi-espace de normale `n` dans le noyau `toolsig` (gap = distance signée au plan de coupe, mêmes équations, `toolShape = blade` + `toolContact = signorini` en fem3d), avec E4 ; ou (b) si (a) dépasse 12 h-agent : `main.cpp` refuse `toolShape = blade` en fem3d hors `toolContact = signorini` et `etude_lois_fem/VERDICT_lame_fem3d.md` établit que la coupe T2 (20 decks, 30-38 h) n'est pas lancée. **[DÉCISION F]** entre (a) et (b) après un devis de 2 h. | (a) : `verify_suite.py --only t1_fem3d_blade` rend PASS (`toolinj ≤ 1,1`) et la variante pénalité « FAIL attendu » ; (b) : `rockim.exe matrice_v2/C_T2_R_P000.cfg` rend le code 1 « lame fem3d : Signorini requis (VERDICT_lame_fem3d.md) ». | 12 (a) / 2 (b) | 1 | 1 | E4, E7 |
| E9 | P03, T26 | Noyau `dpr` par défaut : sans changer le défaut (bit-identité), `law = dpr` sans `dpApex` **ou** sans `rankineDrive = stress` imprime `[WARN] dpr : noyau d'origine, 4 défauts connus (JOURNAL 04/09) — utiliser law = dpr2` ; nouvelle loi-alias `dpr2` = dpr + dpApex + dpTension off + rankineDrive stress + erodeWfrac + meridian power (croissance par addition) ; les 14 decks `dpr` sans `dpApex` des études actives migrent vers `dpr2` ; les cfg `verify_fem3d_dp/rate1/rate2/tension` (4 fichiers présents dans `configs/`) entrent dans la suite (B1). Le script de balayage `keymatrix.py`, qui n'existe aujourd'hui **que dans le scratchpad d'une session d'agent**, est copié et versionné dans `tools/keymatrix.py` (prérequis aussi de F1). **[DÉCISION F]** : promotion de `dpr2` en défaut lors d'un tag majeur ultérieur, ou jamais. | `python tools/keymatrix.py --mode fem3d --law dpr --without dpApex --active` rend 0 deck dans `matrice_v2/` et `perc3d/` ; `rockim.exe <deck dpr sans apex>` imprime la ligne `[WARN] dpr` (déclarée `expect_warn` dans les tests qui l'utilisent) ; bitid B6 identique. | 4 | 1 | 0,5 | B1 |
| E10 | P10 | Hertz avec verdict : `CDP/perc3d/fig_hertz.py` → `vv/hertz/check.py` avec tolérances écrites (force −5 %/+5 %, enfoncement ±5 % vs Hertz analytique, pour Signorini et pénalité ×10) ; avertissement à l'init si `toolContact = penalty` avec `contactPenaltyFactor = 1` sur `mesh = file` (h_min = sliver) ; Hertz court (quart de bloc grossier, ≤ 2 min) dans la suite tier `full`. Les runs fem3d pénalité ×1 antérieurs (matrice v1 A-D) sont marqués « biais de contact −18 % + broche » dans le journal. | `python vv/hertz/check.py out_HZ_sig` rend `PASS force −x % enfoncement +y %` ; `verify_suite.py --only hertz_fem3d` PASS ; `grep -c "contactPenaltyFactor = 1" ` dans les decks actifs = 0 ou avec `expect_warn`. | 4 | 0,5 | 2 | B1 |
| E11 | P12 | Outil libre vs bloqué : rejeu de `PQ2_cdpH_pulseZ2_abqmesh_bvh_nrb` (`toolLockXY = true`, configuration comparable au RP bloqué d'Abaqus) et de la variante libre, sur le binaire porté ; `bilan_nuit.py --tol 5` rend PASS/FAIL par variante ; le journal et §5.20 consignent laquelle est « = Abaqus » et de combien. | `CDP/perc3d/bilan_nuit.py --tol 5` imprime une ligne `PASS/FAIL` pour `bvh` et `bvh_nrb` avec Δpic filtré ; `JOURNAL.md` de CDP a une entrée datée « verdict toolLockXY ». | 2 | 1 | 4 | B7 |
| E12 | P06 | Banc de relais joint → contact : barre pré-entaillée, joint mené à la mort en traction puis inversion en compression ; continuité de force au relais mesurée (saut ≤ 5 % de ft·A) pour `jointDeath = separation` et `damage`, `gcBirth = ramp/penalty` ; le cas `adaptive + damage + penalty` (facteur 0,01 à vie) est **refusé** par le code (F2). | Entrée `relay_joint_contact_2d` tier `fast` : `saut_relais ≤ 5 %` PASS ; `rockim.exe` avec la triplette incohérente rend le code 1 « combinaison interdite (HANDOFF §8.3) ». | 5 | 1 | 0,3 | E3 |
| E13 | P08 | `dtBudgetTangential` : reste opt-in ; entre dans la matrice F1 et dans le banc E1 3D (pénalité) où il est **attendu** nécessaire (contrôle : sans lui, E1-3D pénalité doit diverger ou dériver). | `friction_penalty_3d` PASS avec `dtBudgetTangential = on` ; variante `off` en « FAIL attendu » ou PASS documenté (alors la clé est inutile et le journal le dit). | 1 | 0,5 | 0,3 | E1 |
| E14 | T21 | Amortissement : la recalibration Red Bohus (ROADMAP E2) est une **étude**, hors plan ; le plan impose seulement que tout deck de calibration ou dynamique pose `dampingLocal` explicitement (C5 avertit) et que `docs/REFERENCE.md` dise en une phrase quand l'amortissement local est légitime (quasi-statique) et quand il fausse (dynamique). | `grep -L "^dampingLocal" calib_quick/*.cfg CDP/*.cfg bench_impact/configs/*.cfg` rend 0 ; C5 en place. | 1 | 0,5 | 0 | C5 |
| E15 | T19, T20, T22, T23, T24, T27 | Les 25 défauts datés d'août-septembre sont **déjà corrigés** ; ce que le plan ajoute est le test qui les aurait vus : lot E → E1 `ftScale` inerte (test « ftScale change le pic » en dpr/saksala), E2 `toolShape` (test forme), E3 `chamferLen` (test 3D), E4 `pairKey` dem3d (test tangentiel), E6-E8 (colonnes exportées) ; signe hydro → B11 ; `jointResidualMu` inerte → test « μ résiduel change eFric » ; `nVert_` → F1 (hydro + intrinsic) ; pompe → E4 ; COMDAT → D2 ; broche → C3. Chacun devient une entrée de suite ou de matrice F ; une dizaine d'entrées sont **nouvelles** (ftScale, toolShape, chamferLen 3D, pairKey dem3d, colonnes E6-E8, jointResidualMu, crushCap, critère tBreak, signe force-pénétration, μ par phase, dt parabolique), à ~1 h chacune. | Table `docs/DEFAUTS_2026-08.md` : 25 lignes, colonne « test qui le verrait » remplie d'un nom d'entrée `verify_suite`/`smoke_matrix` existant (`verify_suite.py --list` les contient tous) et, pour chaque entrée nouvelle, le contrôle négatif qui **échoue** quand le défaut est réintroduit (revert local du correctif → FAIL). | 12 | 1 | 1 | B1, F1 |
| E16 | P07 | `gravity` lu à l'init dans les 3D pour tous les scénarios (et non dans `placeTool` seulement) ; parité 2D/3D n° 6 fermée. | Deck de traction 3D avec `gravity = -9.81` : colonne `W_gravité` non nulle ; C1 ne signale plus `gravity` comme lue tardivement. | 1 | 0 | 0,1 | C1 |

### E.2 Ce que E ne décide pas

E ne promeut aucun défaut (dpr, pénalité, `energyBodyForces`, `dtBudgetTangential`) : il rend chaque
défaut **bruyant** quand il est connu pour être faux et fournit à Fernando le banc qui permettra de
décider la promotion à un tag majeur, avec ré-ancrage B6 explicite.

Coût du chantier E : **79 h-agent (69 si E8-b), 12 h-F, ~18 h-machine** (somme des lignes recomptée ; la
version initiale annonçait 76/66 pour une somme de 75/65).

---

## F. Combinaisons : chaque paire de clés principales a tourné une fois

| id | constat | mesure (quoi, où) | critère d'acceptation | h-agent | h-F | h-mach. | dépend de |
|---|---|---|---|---|---|---|---|
| F1 | P18, P19, T24 | `tools/smoke_matrix.py` : pour chaque mode, une liste de clés principales versionnée (`tools/smoke_keys.json`) — fem3d : `law` × {`elastic, dpr, dpr2, mc, saksala, saksala2011, dpdfh, cdp`} et {`kinematics=hencky`, `bulkViscosity`, `toolContact=signorini/penalty`, `quarterModel`, `confiningPressure=50e6`, `cdpCompLength`, `matWeibullM`, `probes`, `toolPulseForce`, `toolLockXY`, `symmetryY+blade`, `absorbing`} ; fdem3d : `law` × {`contact=potential`, `gcActivation`, `jointDeath`, `bulkDamage`, `strainRateDIF`, `groupBond`, `toolContact=signorini`, `confiningPressure`, `absorbing`, `viscosityMunjiza`, `dtBudgetTangential`, `energyBodyForces`} ; fdem 2D : idem + `hydro`, `thermal`, `bedding*`, `insertion=intrinsic/adaptive`, `phases`. Chaque **paire** sur un mini-deck (≤ 2 000 éléments, ≤ 20 µs, ≤ 30 s, OMP 4) ; par run : code 0, pas de NaN, résidu d'énergie < 5 %, hachage stable sur deux exécutions ; résultat `results/smoke_<describe>.json` (paire → ok/fail/interdit). La liste des paires « jamais exercées » est produite par `tools/keymatrix.py` (versionné en E9, aujourd'hui dans un scratchpad de session) et non recopiée à la main. Mode `--quick` (30 paires, 5 min) pour la porte B7 ; complet (≈ 300 paires, ≈ 2,5 h) à chaque tag. | `python tools/smoke_matrix.py --exe build\rockim.exe` rend `N paires : N−K ok, 0 échec, K interdites refusées` ; `results/smoke_v0.19.0.json` au dépôt ; les paires « jamais exercées » de l'état des lieux (hencky+saksala/dpdfh/mc, cdp+hencky+confinement, quarter+dpr, fdem3d+confinement, hydro+intrinsic, dpdfh fdem3d + potentiel/adaptatif/DIF) apparaissent toutes avec `ok`. | 12 | 1 | 2,5/passe | C1, C4, E3 |
| F2 | P06, P18 | Combinaisons interdites refusées par le code : `tools/forbidden_combinations.json` (paire ou triplet, mode, raison, source) chargé par `main.cpp` après init ; liste initiale : `adaptive + jointDeath=damage + gcBirth=penalty` ; `toolSignoriniRelax > 0` ; `thermal + TI` ; `law + phases` ; `blade fem3d + penalty` (E8-b) ; `quarterModel + symmetryY` ; `hydro` hors fdem (C6). La matrice F1 vérifie que chaque interdite rend le code 1 avec le message. | `python tools/smoke_matrix.py --forbidden-only` rend `K/K refusées` ; ajouter une interdite au json sans test rend la CI rouge (`smoke_matrix --check-json`). | 3 | 1 | 0,2 | F1 |
| F3 | P18 | Trois combinaisons « longues » (au-delà de la fumée) tournent une fois avec bilan complet : `cdp + hencky + bv + confiningPressure 100 MPa` (quart de bloc) ; `saksala2011 + hencky + bv` (xval court) ; `dpdfh + signorini + bv` ; verdict dans `docs/COMBINAISONS_2026-09.md`. | Trois `out_*` avec `run_info.json`, `energy_check.py` PASS (< 2 %), et une ligne de verdict chacune. | 3 | 1 | 6 | F1 |

Coût du chantier F : **18 h-agent, 3 h-F, ~9 h-machine + 2,5 h par passe complète**.

---

## G. Un dossier V&V rejouable en une commande

Chaque validation citée dans un rapport a un dossier `R\vv\<nom>\` avec `README.md` (quoi, pourquoi,
source de la cible, date), `expected.json` (cibles, tolérances, provenance), `run.py` (lance le ou les
runs, ou lit les sorties existantes) et `check.py` (une ligne `PASS`/`FAIL` + chiffres). `vv/run_all.py
--tier short` (< 1 h) et `--tier long` (nuit). Résultats dans `vv/results/<date>_<describe>.json`.

| id | constat | validation | commande | tolérance écrite | h-agent | h-F | h-mach. |
|---|---|---|---|---|---|---|---|
| G1 | P-m2 | Abaqus GBM/continuum saksala2011 (xval, maillage identique) : déplacer `G\xval_abaqus\compare_xval.py` et les CSV Abaqus (petits) dans `vv/abaqus_saksala2011/`. | `python vv/run.py abaqus_saksala2011` | pic filtré 10 µs ±15 %, perte KE outil ±10 %, restitution ±0,05 (critères existants) | 2 | 0 | 0,5 |
| G2 | P16 | Abaqus CDP impact (P_cdpQ_v11, 42,5 kN) : `bilan_nuit.py` → `check.py` ; CSV Abaqus extraits stockés dans `vv/abaqus_cdp_impact/abq/` ; variante `toolLockXY` (E11) comme cas de référence. **Plus** le point matériel cdp vs EF Abaqus (`CDP/compare_abaqus.py`, aujourd'hui sans seuil : RMS pré-pic et rapport pic EF/point 0,855-0,896 imprimés) → `vv/abaqus_cdp_matpoint/check.py`, tier `short`. | `python vv/run.py abaqus_cdp_impact --tier long` ; `python vv/run.py abaqus_cdp_matpoint` | impact : pic filtré ±5 %, travail ±5 %, pénétration ±10 % ; point matériel : RMS pré-pic ≤ 2 % de σ_pic, rapport pic EF/point ∈ [0,85, 0,90] (valeurs consignées) **[DÉCISION F sur les seuils]** | 5 | 1 | 4 |
| G3 | P10 | Hertz (E10). | `python vv/run.py hertz` | ±5 % force et enfoncement | 1 | 0 | 0,5 |
| G4 | P21 | Yan 2023 UCS (adaptatif) contre **l'article** (≈ 51 MPa), non contre le compteur Linux. | `python vv/run.py yan2023_ucs` | UCS 51 ± 2 MPa, E 99 ± 2 %, une bande à 65 ± 10° | 2 | 0,5 | 0,5 |
| G5 | P13 | AbuAisha : un cas court (`e1`, T/10) en tier `short` (signe + ouverture Parker ±10 %) ; les 11 decks longs en tier `long` avec résidu réel 2,5-10,5 % consigné comme **cible** (non 0,0125 %). | `python vv/run.py abuaisha [--tier long]` | signe > 0 ; ouverture ±10 % Parker ; résidu ≤ 11 % | 4 | 1 | 9 (long) |
| G6 | P21 | Lisjak TI (t14c, deux angles). | `python vv/run.py lisjak_ti` | 1,5910 et 4,3307 ± 1 % | 2 | 0 | 0,5 |
| G7 | P15 | DP-DFH point matériel contre le pilote Fortran : `selftest-dpdfh` écrit le CSV, `check.py` le diffuse contre `VUMATS\dfh\dfh_check.csv` (copié dans `vv/dpdfh_fortran/ref/`) ; idem `saksala2011` contre `test_saksala_2011.f90` (référence régénérée une fois par `ifort`, CSV stocké). Les deux entrent dans la suite `fast` (B1). | `python vv/run.py dpdfh_fortran` ; `saksala2011_fortran` | écart max ≤ 1e-10 (4,7e-12 et 8e-14 mesurés) | 3 | 0,5 | 0 |
| G8 | P-info thermique | Thermique t9b : forme fermée Eα∆T/(1−ν). | `python vv/run.py thermique_t9b` | ±0,1 % (0,011 % mesuré) | 1 | 0 | 0,2 |
| G9 | P21 | Signorini T0/T1 : déjà dans la suite ; `vv/signorini/` ne fait que pointer sur `verify_suite --only toolcontact`. | `python vv/run.py signorini` | tolérances existantes | 0,5 | 0 | 0,1 |
| G10 | P21 | `vv/run_all.py` + rapport `vv/RAPPORT_<date>.md` généré (table nom / cible / mesuré / écart / verdict / commit). | `python vv/run_all.py --tier short` | tous PASS | 3 | 0,5 | 1 |

Critère global : `python vv/run_all.py --tier short` rend `11/11 PASS` (10 validations + le point matériel
cdp de G2) ; `--tier long` a un résultat daté de moins de 30 jours dans `vv/results/`. Coût du chantier G :
**23,5 h-agent, 3,5 h-F, ~17 h-machine (long compris)**.

---

## H. La documentation est générée là où elle peut l'être, et vérifiée partout

| id | constat | mesure (quoi, où) | critère d'acceptation | h-agent | h-F | dépend de |
|---|---|---|---|---|---|---|
| H1 | D05, D04 | Registre des clés **dans le code** : chaque appel `cfg_.getX("clé", défaut)` est précédé d'un commentaire structuré `// @key clé [unité] modes=fem3d,fdem3d : une phrase` ; `tools/gen_keys.py` scanne `src/` et `include/`, produit `docs/CLES.md` (clé, type, défaut, modes, unité, phrase, fichier:ligne) et `tools/keys_registry.json` (lu par C1 pour la suggestion et par C6 pour les modes) ; `--check` échoue sur toute clé lue sans `@key` ou tout `@key` sans lecture. Les 31 clés absentes et les ~40 tokens orphelins de la doc sont résolus à cette occasion. | `python tools/gen_keys.py --check` rend `391 clés lues, 391 documentées, 0 orpheline` (nombre exact au jour du run) et le code 0 ; en CI (B8) ; ajouter un `getd("nouvelleCle")` sans `@key` rend le code 1. | 14 | 2 | — |
| H2 | D04 | Séparation journal / référence : `DOCUMENTATION_rockim.md` (2 351 l.) scindé en `docs/REFERENCE.md` (sections 1-4, 6, 8 stables + `CLES.md` inclus) et `docs/journal/2026-08.md`, `2026-09.md` (les §5.x datés, §5.20, les 73 horodatages) ; `DOCUMENTATION_rockim.md` devient un index de 30 lignes. Règle CONTRIBUTING : la référence ne contient ni date ni nom de binaire. | `grep -cE "2026-0[89]-[0-9]{2}" docs/REFERENCE.md` rend 0 ; `grep -c "rockim_f2w" docs/REFERENCE.md` rend 0 ; `wc -l DOCUMENTATION_rockim.md` ≤ 40. | 8 | 2 | H1 |
| H3 | D06 | Guide Overleaf `guide_rockim.tex` : table des fichiers générée par `tools/gen_arch.py --tex` (wc -l, rôle, en-tête) et incluse (`\input{arch_table.tex}`) ; nouvelle section « Un pas de temps de Fem3dSolver » et « … de Fdem3dSolver » (H4) ; suppression des références à `rockim_gui.py`/`export_abaqus.py` ; `\date` = tag. Boucle Overleaf (pull → MIR → push) documentée dans la mémoire, appliquée. | `python tools/gen_arch.py --check guide_rockim.tex` rend « table conforme au commit v0.19.0 » ; `grep -c "Fem3dSolver" guide_rockim.tex` ≥ 10 ; PDF compilé sans erreur (MiKTeX) et poussé sur Overleaf. | 8 | 2 | H4 |
| H4 | D10 | Carte d'un pas de temps pour les 4 solveurs principaux (`FdemSolver::step`, `Fdem3dSolver::step`, `Fem3dSolver::step`, `FemSolver::step`) : `docs/PAS_DE_TEMPS.md` (ordre des appels, ce que chaque étape lit et écrit, où les gardes C4 s'exercent) ; `tools/check_stepmap.py` vérifie que les fonctions citées existent et apparaissent dans cet ordre dans le corps de `step()`. | `python tools/check_stepmap.py` rend `4 solveurs : ordre conforme` ; en CI. | 6 | 1 | — |
| H5 | D03 | `CONTRIBUTING.md` : comment ajouter une clé (H1 + F1 + test), une loi (selftest + registre + `smoke_keys.json`), un test (référence par plateforme), une validation (G) ; checklist de PR (porte B7 verte, CHANGELOG, bitid attendu ou ré-ancré, doc dans le même commit) ; règle des cherry-pick (A8) ; convention de langue (français, accents, messages d'erreur en français). | `Test-Path R\CONTRIBUTING.md` ; exercice I2 réussi en suivant le fichier seul. | 4 | 1 | H1 |
| H6 | D-m4, D-m5 | Schéma des sorties : `out/history_schema.json` (colonnes, unités, version) écrit à côté de `history.csv` ; `docs/REFERENCE.md` §6 liste les 10 en-têtes par mode × scénario (générés par `gen_keys.py --history`). Les 4 scripts lisant par index sont convertis en lecture par nom. | Tout `out_*/history_schema.json` a autant d'entrées que de colonnes ; `grep -l "\[:, *[0-9]" G\**\*.py` rend 0 fichier actif. | 3 | 0 | — |
| H7 | D-m1, D-m2 | Hiérarchie des documents d'entrée : `README.md` (anglais, théorie) renvoie à `LISEZ_MOI.md` (ordre de lecture en 10 lignes : constitution → CONTRIBUTING → REFERENCE → PAS_DE_TEMPS → journal) ; `GUIDE_rockim.md` de `G` marqué obsolète ; `FICHE_rockim.md` de `G` mis à jour et index des 6 journaux. | `LISEZ_MOI.md` ≤ 60 lignes ; `FICHE_rockim.md` a une entrée datée ≥ 2026-09 ; `grep -c "rockim_gui\|export_abaqus" docs/*.md README.md` rend 0. | 3 | 1 | H2 |

Coût du chantier H : **46 h-agent, 9 h-F**.

---

## I. Réduire la dépendance à une seule personne et à un seul arbre

| id | constat | mesure (quoi, où) | critère d'acceptation | h-agent | h-F | dépend de |
|---|---|---|---|---|---|---|
| I1 | D03 | Deux relecteurs par changement : tout changement de `src/`, `include/`, `tools/verify_suite.py`, `tools/bitid_refs.json` passe par une branche + PR GitHub (même seul) ; relecture 1 = agent (`/code-review` ou session dédiée, commentaire dans la PR : « défaut silencieux possible ? test qui le verrait ? ») ; relecture 2 = Fernando (approbation) ; fusion seulement CI verte + porte locale jointe (`results/gate_*.json` dans la PR). | `gh pr list --state merged --limit 20 --json reviews,statusCheckRollup` : chaque PR a ≥ 1 revue et `SUCCESS` ; `git log --merges main` : chaque fusion cite un numéro de PR. | 2 | 1/PR | A6, B8 |
| I2 | D03 | Walkthroughs enregistrés : `docs/walkthroughs/{contact_outil, contact_fragments, lois_MatLaw, suite_et_refs, build_et_porte, entrees_Config}.md`, chacun = parcours de lecture de 30 min (fichiers, lignes, invariants, banc qui le prouve) + un exercice final ; enregistrement vidéo facultatif (OBS, 20 min) déposé dans la base. Test d'acceptation : une session d'agent **neuve** (sans mémoire) ajoute une clé factice `demoKey` de bout en bout (registre, garde, test, doc, PR) en suivant CONTRIBUTING + walkthrough, en < 1 h. | 6 fichiers présents ; la PR de l'exercice est fusionnée puis révertie (`git log --grep demoKey` montre les deux commits) ; durée consignée. | 12 | 3 | H5 |
| I3 | D11 | Un seul arbre (A7) et un seul lieu de résultats validés (la base `phd_geothermie`). | Critère A7. | — | — | A7 |
| I4 | D03 | Un second humain (encadrant, co-doctorant ou collègue du laboratoire) fait le walkthrough `suite_et_refs` et lance la porte une fois **[DÉCISION F : qui]**. Hors du pouvoir d'un agent ; le plan le nomme parce que la robustesse d'un code à un seul mainteneur a un plafond. | `results/gate_*.json` avec un `hostname` différent de la machine de Fernando. | 0 | 2 | I2 |
| I5 | D14, D13 | Duplication 2D/3D : pas de refonte (voir §L) ; seulement (a) le contact outil fem3d/fdem3d partagé via `ToolSignorini.hpp` (E7, E8-a) et une brique `ToolContact3d.hpp` pour `placeTool`/`setupConfinement`/`setupBoundaries` **si** bitid reste identique ; (b) gel de la croissance : `tools/funclen.py --max 150 --baseline tools/funclen_baseline.json` en CI, échoue seulement si une fonction **dépasse sa taille de référence** ou si une nouvelle fonction > 150 l. apparaît. | `python tools/funclen.py --check` rend 0 en CI ; lignes identiques `comm` entre `Fem3dSolver.cpp` et `Fdem3dSolver.cpp` ≤ 40 (126 aujourd'hui) après (a), bitid B6 identique. | 10 | 1 | E7 |
| I6 | D12 | VTU : opt-in `vtuFormat = binary` (XML appended raw, Float32 pour les champs, Float64 pour les points) ; hachage B6 calculé sur `history.csv` et sur les champs relus par `pyvista` (pas sur les octets) pour rester comparable entre formats. | Deck avec `vtuFormat = binary` : taille < 30 % de l'ASCII, `pyvista.read` rend des champs égaux à 1e-6 relatif à l'ASCII ; défaut ASCII inchangé (bitid). | 4 | 0,5 | B6 |

Coût du chantier I : **28 h-agent, 8 h-F + 1 h-F par PR**.

---

## J. Constats mineurs et informatifs : traitement groupé

| constat | traitement | où |
|---|---|---|
| T-m1 résultats Linux non persistés | B8 (CI produit `results/ci/*.json`) | B8 |
| T-m2 logs de build sans avertissement | D1b (`/W4 /WX`, `build.log` conservé par la porte) | D1b |
| T-m3 antivirus | D4 | D4 |
| T-m4 parsing (`getb`, ligne sans `=`, clé répétée) | C1 | C1 |
| T-m5 instabilité 3D gelée | E5 | E5 |
| T-info inventaire des tests, commentaire « ≈ 1 min » périmé | B1 met à jour l'en-tête (`--list` imprime la durée mesurée de la dernière passe) | B1 |
| T-info exceptions sous OpenMP (points faibles : « degenerate tet » sans id) | C3 | C3 |
| P-m1 limites 3D documentées (confinement scalaire, `cutterFloor` lu FALSE, chamferLen inerte en 2D) | C1 (`getb` strict), H1 (`@key` porte « inerte en 2D » et C1 en fait un avertissement « clé lue mais sans effet dans ce mode ») | C1, H1 |
| P-m2 saksala2011 xval hors suite | G1 | G1 |
| P-m3 mc/elastic/briques sans suite | B1 (`selftest-mc`, `selftest-kinematics`, `selftest-bv`), E9 (briques dans `verify_fem3d_*`) | B1, E9 |
| P-info algorithmes de contact, bancs existants, thermique 2D, ce qui existe en 3D, inventaire des lois, piles exercées, mécanique de la suite | Reversés dans `docs/REFERENCE.md` (H2) et `docs/PAS_DE_TEMPS.md` (H4) ; la matrice F1 remplace la liste manuelle des piles exercées | H2, H4, F1 |
| D-m1 outils cités absents, D-m2 cinq documents d'entrée | H7 | H7 |
| D-m3 bancs non branchés | B1, B10 | B1, B10 |
| D-m4 en-têtes history.csv, D-m5 lecteurs par index | H6 | H6 |
| D-m6 hygiène (47 exe, 13 obj, 126 out_*, 5 `.w0`) | A2 | A2 |
| D-m7 journal éclaté sans index | H7 (index dans `FICHE_rockim.md`) ; les journaux de chantier restent datés et ne sont pas réécrits | H7 |
| D-info post-traitement (13 outils, 190 uniques) | A7 réduit les copies ; pas d'autre action | A7 |

---

## K. Calendrier en trois vagues et coût total

### K.1 Vagues

| vague | contenu | durée indicative | ce qui **doit attendre** la fin de la vague | ce qui peut tourner **en parallèle** |
|---|---|---|---|---|
| **1 — couper les racines** | **Jour 0 (avant tout) : B5 + B6 sur `rockim_f2w18.exe` existant** (empreinte de départ, aucun build) ; puis A1-A5, D1a (bit-identique à l'ancre), D1b, D2-D3, C0, C1-C4, C6 (garde `hydro`), C7-C8, B1-B4, B7 (porte v1), E6, E16 | 5 jours ouvrés d'agent (≈ 95 h-agent), 3-4 h-F/jour de relecture et décisions | **Tout nouveau run de campagne** : la matrice v2 (phases C1 + D déjà à refaire pour la broche), l'étude hétérogénéité fem3d, le bayésien. Raison : tant que C1/C3/C5 ne sont pas en place, un run peut être faux-mais-plausible et devra être refait. Les 42 decks lancés le 04/09 07:00 sur w10 sont concernés par la broche fantôme (maillages T1/Q1) : ils sont à rejouer de toute façon. | Dépouillement et rédaction des résultats déjà acquis (Abaqus CDP 42,5 kN, phase B des briques) ; préparation des decks bayésiens et hétérogènes (ils seront re-validés par C1 avant lancement) ; lecture des walkthroughs à mesure qu'ils sortent. |
| **2 — la physique a ses bancs** | A6-A9, B8-B11, C5, E1-E5, E7-E13, E15, F1-F2, G1, G3, G4, G6-G10, D4-D5 | 10 jours ouvrés d'agent (≈ 130 h-agent), 2 h-F/jour | Les études **fdem3d** en phase débris (attendent E5) ; la coupe fem3d T2 (attend E8) ; toute comparaison Abaqus/FDEM 3D (attend C7 puis une capacité bv/hencky fdem3d, hors plan) ; les études utilisant `dpr` sans `dpr2` (E9). | Les études **fem3d** cdp/dpr2 + signorini + hencky + bv (hétérogénéité, bayésien, matrice v2 C1/D rejouée) sur le binaire porté `bin/rockim_v0.19.x.exe` : elles sont couvertes par C1-C5, B1, B6 dès la vague 1 ; chaque run porte son `run_info.json`. |
| **3 — documenter et déléguer** | H1-H7, I1-I6, F3, G2, G5 (long), E14 | 10 jours ouvrés d'agent (≈ 100 h-agent), 1-2 h-F/jour | Rien : les études continuent ; seule la rédaction du manuscrit gagne à attendre H3 (guide Overleaf à jour). | Toutes les études ; la nightly (B8) et `vv --tier long` tournent seules. |

### K.2 Coût total

| chantier | h-agent | h-F | h-machine |
|---|---|---|---|
| A git et versions | 21,5 | 8 | 1 |
| B suite juge | 49 | 8 | 13 (+1/porte) |
| C entrées et gardes | 29,5 | 4 | 2,5 |
| D build | 18,5 | 1,5 | 1,5 |
| E contact et physique | 79 (69) | 12 | 18 |
| F combinaisons | 18 | 3 | 9 (+2,5/passe) |
| G V&V | 23,5 | 3,5 | 17 |
| H documentation | 46 | 9 | 0 |
| I dépendance | 28 | 8 (+1/PR) | 0 |
| **total** | **≈ 315 h-agent** (305 si E8-b) | **≈ 57 h-F + relectures** | **≈ 62 h + porte/passes** |

Ordre de grandeur cohérent avec les chantiers de la semaine (une clé opt-in avec bit-identité et banc =
1-2 h-agent ; une matrice de fumée = 1 journée ; CMake + CI = 1 journée) : les deux postes qui dépassent
franchement ces étalons sont D1b (10 h : 200-600 avertissements à corriger sans changer l'algèbre) et
E8-a (12 h : lame Signorini fem3d), tous deux marqués incertains ; H1 (14 h) est linéaire dans le nombre
de clés (391 × ~1 min + générateur).

Cinq semaines calendaires si l'agent travaille cinq jours par semaine et Fernando relit chaque soir ; le
temps machine est absorbé la nuit (aucune licence Abaqus n'est consommée sauf G2 et E11, qui réutilisent
des extractions existantes). Les 62 h-machine s'entendent à 14 cœurs ; les portes (≈ 1 h) et passes de
matrice (2,5 h) s'ajoutent à chaque tag.

### K.3 Ordre de dépendance (résumé)

```
B5, B6 sur rockim_f2w18.exe (jour 0, aucune dépendance) → A1 → A2..A5 → D1a (= ancre) → D1b → D2, D3
C0 → C1 → C2 → C5 ; C1 → C7 ; C6 (garde hydro) et C3, C4 indépendants ; H1 facultatif pour C1/C6 (message)
B1..B4 (+ B5, B6) → B7 (porte v1) → A6 (fusion) → A7 (un arbre) → B8 (CI) → I1 ; F1 → B7-quick (v2)
E3 → E1 → E2, E13 ; E4 → E8 ; E7 → E8-a ; C3 → E6 ; C4 → E5 ; B1 → E9, E10 ; E9 (keymatrix.py) → F1
C1, C4, E3 → F1 → F2 ; G7 → B1 ; H1 → H2 → H3, H7 ; H4 → H3 ; H5 → I2 → I4
```

Règle d'ordre appliquée : ce qui **protège** des régressions (empreinte B5/B6, dépôt A1, porte B7) précède
ce qui **ajoute** (gardes C, bancs E, matrice F) ; ce qui **bloque les études en cours** est nommé dans la
colonne « doit attendre » de K.1 (matrice v2, hétérogénéité fem3d, bayésien en vague 1 ; fdem3d débris et
coupe T2 en vague 2).

---

## L. Définition finale : « la fragilité a disparu quand… »

Dix critères binaires, tous exécutables par une commande, tous rejouables par la CI ou la porte locale.
Le plan est **terminé** le jour où les dix rendent vrai sur le même commit tagué.

| n° | critère | commande qui le vérifie | rend |
|---|---|---|---|
| 1 | Le code vivant est sous git, propre, tagué, dans un seul arbre de travail | `git -C R status --porcelain ; git -C R describe --tags --exact-match ; git -C R worktree list \| Measure-Object -Line` | vide ; un tag ; ≤ 3 lignes (`rockim_p1`, `R`, `studio_wt`) |
| 2 | Le build est reproductible et muet à `/W4 /WX` avec suivi de dépendances | `tools\build.ps1 ; Select-String build\build.log -Pattern "warning C" \| Measure-Object` | code 0 ; 0 |
| 3 | Chaque binaire et chaque run portent leur empreinte | `build\rockim.exe --version ; Get-ChildItem out_* -Filter run_info.json -Recurse \| Measure-Object` = nombre de `out_*` | ligne de version conforme ; égalité |
| 4 | La suite juge sur chaque plateforme avec ses propres références, à 1 fil, résultats persistés | `python tools\verify_suite.py --tier full` (refs auto) ; `gh run list --limit 1 --json conclusion` | `N/N` ; `success` |
| 5 | La bit-identité est générique et ancrée au dépôt | `python tools\bitid.py --exe build\rockim.exe` | `8/8 IDENTIQUE` (ou ancre changée avec ligne CHANGELOG dans le même commit) |
| 6 | Aucune entrée n'est acceptée en silence | `python tools\verify_suite.py --only "tg_"` (clé inconnue, getb, orphelin ×2 modes, tet dégénéré, NaN ×6 modes, unités, hydro3d, clé hors mode, combinaison interdite) | tous en « FAIL attendu » = PASS du contrôle négatif |
| 7 | Chaque défaut silencieux d'août-septembre a le test qui l'aurait vu | `python tools\verify_suite.py --list` contient les 25 noms de `docs/DEFAUTS_2026-08.md` | 25/25 |
| 8 | Toute paire de clés principales a tourné sans NaN, avec bilan fermé et hachage stable ; les interdites sont refusées | `python tools\smoke_matrix.py --exe build\rockim.exe` | `0 échec, K/K interdites refusées` |
| 9 | Chaque validation citée se rejoue en une commande avec tolérance écrite | `python vv\run_all.py --tier short` ; `Get-ChildItem vv\results` (long < 30 j) | `11/11 PASS` ; fichier daté |
| 10 | La documentation ne peut plus dériver du code | `python tools\gen_keys.py --check ; python tools\check_stepmap.py ; python tools\gen_arch.py --check guide_rockim.tex ; python tools\funclen.py --check` | quatre codes 0 |

Ces dix critères sont exactement le contenu de `tools/gate.ps1 --full` : la porte complète **est** la
définition de la robustesse, et elle tourne à chaque tag.

---

## M. Ce que le plan ne fait pas, et pourquoi

| hors plan | raison |
|---|---|
| **MPI, GPU, refonte de performance** | La fragilité relevée est une fragilité de **correction**, pas de débit : 25 défauts silencieux, zéro plainte de temps de calcul non expliquée (le facteur 14 du maillage était une erreur d'attribution). Une machine à 14 cœurs suffit aux études prévues ; MPI casserait le déterminisme à fil fixe sur lequel reposent B6 et la bit-identité, et doublerait la surface de code (principe III déjà coûteux). À reconsidérer seulement si une étude chiffre un besoin > 10× le débit actuel. |
| **Hydro et thermique en 3D** | Ce sont des **capacités**, pas de la robustesse. L'hydro 2D est une pression uniforme sans loi cubique ni leak-off, validée sur un seul benchmark 2D avec un résidu de 2,5-10,5 % ; la porter en 3D (mouillage de faces de tets, `updateWetBoundary` 3D, cavité 3D) coûte 2-3 semaines et n'aurait aucune cible de validation 3D. Le plan se limite à **refuser** `hydro` et `thermal` hors 2D (C6), ce qui suffit à supprimer le mode de panne. |
| **bulkViscosity Abaqus et hencky en fdem3d** | Nécessaires pour une comparaison Abaqus/FDEM 3D à l'identique (P14), mais ce sont des capacités nouvelles ; le plan lève seulement l'ambiguïté de nom (C7). À planifier comme spec après la vague 2 si une étude FDEM 3D vs Abaqus est décidée. |
| **Refonte de la duplication 2D/3D (5 800 lignes homonymes)** | Le principe III (2D et 3D de front) est un choix de conception de Fernando ; le refondre en briques partagées est un chantier de 100+ h-agent à risque élevé de changer l'algèbre flottante (ancre B6). Le plan gèle la croissance (I5-b) et ne partage que le contact outil (E7, E8-a), où la copie manuelle a déjà produit un défaut. |
| **Recalibration Red Bohus (ROADMAP E2), promotion de `dpr2`/Signorini/`energyBodyForces` en défaut** | Ce sont des **décisions scientifiques** de Fernando, avec ré-ancrage de bit-identité et impact sur les résultats publiés. Le plan fournit les bancs (E9, E10, E3) et les avertissements ; il ne tranche pas. |
| **Réécriture de `verify_suite.py` (77 Ko, 60 regex sur stdout)** | Les regex sont fragiles (tout `print` modifié casse un repère) mais fonctionnent ; les remplacer par une lecture de `history.csv`/`run_info.json` est souhaitable et se fera **test par test** au fil de B1/B10, sans réécriture globale qui ferait perdre la provenance des 111 références. |
| **Traduction intégrale des commentaires en une langue** | 9 367 lignes de commentaires mixtes ; le plan fixe la convention pour le code nouveau (H5) et n'entreprend pas de réécrire l'existant. |
| **Achèvement de la réplique Solidity (Yang/Guo, P05)** | Le déficit de frottement (32 J publiés contre 0,58-0,69 J) est mesuré sur un run arrêté à 19 % de la passe, avec la réserve écrite que le régime de broyage n'est pas atteint. Mener ce run à 100 % (plusieurs heures) et rendre un verdict est une **étude**, pas de la robustesse ; le plan fournit ce qui lui manquait pour être concluante : le bilan fermable sur pièces (E3), le banc de frottement à deux corps (E1) et le relais joint → contact (E12). |
| **Suppression de quoi que ce soit sans accord** | Toutes les suppressions (worktrees, `.w0`, exe, archive hors git, branches) sont marquées **[DÉCISION F]** ; le plan déplace et archive, il ne supprime pas. |

---

## Annexe A. Table de correspondance constat → mesure

Chaque constat **bloquant** ou **majeur** de l'état des lieux (ordre du tableau source) et la ligne du plan
qui le traite. Les constats mineurs/info sont en §J.

### Section « Tests, build, robustesse des entrées »

| id | gravité | thème (état des lieux) | mesure(s) |
|---|---|---|---|
| T01 | bloquant | dépôt git de rockim_f2 | A1, A2, A9 |
| T02 | bloquant | clés inconnues | C1, C2, H1 |
| T03 | majeur | dernier run MSVC (aucune suite sur w0…w18) | B5, B7 |
| T04 | majeur | tier full MSVC 42/50 | B4 |
| T05 | majeur | « 35/40 » à OMP 4 | B3 |
| T06 | majeur | références de plateforme (pas de refs MSVC/Linux) | B2 |
| T07 | majeur | couverture fem3d nulle | B1, B10 |
| T08 | majeur | bit-identité fem3d, ancre changée | B6 |
| T09 | majeur | CI absente | B7, B8, B9 |
| T10 | majeur | deux chaînes de build divergentes | D1a, D1b, D5 |
| T11 | majeur | lien partiel COMDAT | D2 |
| T12 | majeur | versionnage par suffixe d'exe | A3, A5, D3 |
| T13 | majeur | fragmentation des branches, cherry-pick | A6, A8 |
| T14 | majeur | fragBrushV inerte dans 10 decks | C0 |
| T15 | majeur | unités non vérifiées | C5 |
| T16 | majeur | import de maillage fem3d (orphelin, tet dégénéré, slivers) | C3 |
| T17 | majeur | masse nulle = broche fantôme | C3, E6 |
| T18 | majeur | détection des NaN aveugle/absente | C4, C8 |
| T19 | majeur | lot E dette silencieuse | E15 |
| T20 | majeur | signe hydro, contrôle en norme | B11, E15 |
| T21 | majeur | amortissement et calibration | E14, C5 |
| T22 | majeur | défauts 22-29/08 | E15 |
| T23 | majeur | suite et outillage 28-30/08 | B2, E15, G5 |
| T24 | majeur | crash nVert_ (hydro + intrinsic) | F1, E15 |
| T25 | majeur | contact outil, pompe ×408 | E4 |
| T26 | majeur | noyau d'origine fem3d, 4 défauts | E9 |
| T27 | majeur | lien partiel + broche (03-05/09) | D2, C3 |

### Section « Physique et couverture »

| id | gravité | thème | mesure(s) |
|---|---|---|---|
| P01 | bloquant | instabilité latente 3D jamais rejouée | E5 |
| P02 | bloquant | coupe fem3d aveugle (lame en pénalité) | E4, E8 |
| P03 | bloquant | dpr par défaut = 4 défauts | E9 |
| P04 | majeur | bancs de contact fragments absents | E1, E2 |
| P05 | majeur | fiabilité Yang/Guo, bilan non fermable | E3, E1, E12 (bancs) ; achèvement du run = étude, §M |
| P06 | majeur | relais joint → contact | E12, F2 |
| P07 | majeur | bilan d'énergie aveugle aux pompes, gravité 3D | E3, E4, E16 |
| P08 | majeur | chaos et plateforme, dtBudgetTangential | B3, B4, E13 |
| P09 | majeur | trois implémentations du contact outil | E7 |
| P10 | majeur | Hertz sans verdict, pénalité ×1 par défaut | E10, G3 |
| P11 | majeur | broche fantôme, garde partielle | C3, E6 |
| P12 | majeur | outil libre vs bloqué (42,5 kN) | E11, G2 |
| P13 | majeur | hydro 2D seulement, non gardé | C6, B11, G5 |
| P14 | majeur | collision de nom bulkViscosity | C7 |
| P15 | majeur | dpdfh/saksala2011 sans script de diff | G7, B1 |
| P16 | majeur | cdp vs Abaqus sans seuil (impact **et** point matériel) | G2 |
| P17 | majeur | suite jamais rejouée sur w0-w18 | B5 (= T03) |
| P18 | majeur | matrice fem3d jamais exercée | F1, F3 |
| P19 | majeur | matrice fdem3d | F1 |
| P20 | majeur | clés inconnues silencieuses | C1 (= T02) |
| P21 | majeur | validations non rejouables | G1-G10 |

### Section « Documentation, processus, dépendance »

| id | gravité | thème | mesure(s) |
|---|---|---|---|
| D01 | bloquant | pas un dépôt git | A1 (= T01) |
| D02 | bloquant | clés inconnues | C1 (= T02) |
| D03 | bloquant | une seule personne connaît le code | H5, I1, I2, I4 |
| D04 | majeur | documentation dérivée en journal | H2 |
| D05 | majeur | dérive code → doc (391/342 clés) | H1 |
| D06 | majeur | guide Overleaf périmé sur Fem3dSolver | H3 |
| D07 | majeur | constitution absente du dossier de travail | A4 |
| D08 | majeur | bit-identité non automatisée | B6 (= T08) |
| D09 | majeur | suite : refs en dur, regex, pas de refs MSVC | B2, B1 (= T06) |
| D10 | majeur | pas de carte d'architecture ni de pas de temps | H4, H3 |
| D11 | majeur | sept copies de l'arbre | A7, I3 |
| D12 | majeur | VTU ASCII 5-10 Mo/image | I6 |
| D13 | majeur | 27 fonctions > 150 lignes | I5-b |
| D14 | majeur | duplication 2D/3D institutionnalisée | I5-a, E7 |

---

## Annexe B. Décisions attendues de Fernando avant la vague 1

1. **C1** : la clé inconnue devient une erreur par défaut (`unknownKeys = warn` pour s'en sortir). C'est
   le seul changement de défaut du plan ; sans lui, C6, C7, F2 et le critère 6 n'existent pas.
2. **A2/A7/A9** : suppression des `.w0`, des worktrees p2-p4 après archivage, de la copie hors git après
   vérification.
3. **A6** : sort des 87 lignes de `joint-handoff` et des branches `insertion-pointe`, `dif-intrinseque`.
4. **E8** : lame fem3d en Signorini (12 h) ou interdite avec verdict (2 h), après devis.
5. **E9** : `dpr2` comme alias bruyant maintenant ; promotion en défaut, ou jamais, à un tag majeur.
6. **G2** : seuils de la comparaison Abaqus CDP (proposés ±5 % pic, ±5 % travail, ±10 % pénétration).
7. **I4** : le second humain.
8. **K.1** : arrêt des campagnes pendant la vague 1 (5 jours), ou poursuite en acceptant de rejouer.
9. **D1a** : conserver `build_f2.cmd` un mois comme secours, puis le supprimer.

Ne demande **aucune** décision et peut commencer immédiatement : B5 + B6 sur `rockim_f2w18.exe` (empreinte de
départ), C6 (garde `hydro` dans `main.cpp`), C3, C4 (gardes de maillage et de NaN, opt-out impossible car
il n'y a pas de comportement légitime à préserver).

*Rédigé le 2026-09-05 pour arbitrage ; ne modifie rien d'autre que ce fichier.*

*Relecture adverse du 2026-09-05 (soir) : 62/62 constats bloquants/majeurs retrouvés dans l'annexe A avec
mesure et critère ; 14 références fichier:ligne de l'état des lieux revérifiées dans la source (toutes
exactes) ; corrections apportées : ordre B5/B6 avant A1/D1 sur l'exe w18 existant, scission D1a/D1b,
C1/C6 sans dépendance au registre H1, B3 reformulé (le runner force déjà `--threads`, défaut 1), B7 sans
dépendance à F1 en vague 1, B8 sans dépendance à H1/H4, A6/A7/L1 corrigés (`studio_wt` sur une branche
`claude/*`, pas de branche `studio`), `keymatrix.py` à versionner (E9, F1), tests_f2 = 47, C5 avec valeurs
sentinelles et balayage préalable des 954 decks, G2 étendu au point matériel, P05 explicité en §M, E15
recoûté 8 → 12 h, sommes E/G/total recomptées. Copie antérieure : `PLAN_ROBUSTESSE_2026-09-05.md.bak`.*

## Annexe C. Décisions prises par Fernando le 2026-09-05 (20:00-21:10)

1. **C1** — la clé inconnue devient une **erreur** par défaut, avec un message qui **explique la raison** (faute de frappe → suggestion ;
   obsolète → nouveau nom ; autre mode → le mode où elle agit ; inconnue de rockim) ; `unknownKeys = warn` pour la transition ;
   configuration effective écrite par run. En cours : w21.
2. **A2/A7/A9** — **rien n'est supprimé**. Nouvel arbre `FDEM/rockim_g0/` né sous git (branche `g0`, état de `rockim_f2` à la validation
   de w21) ; `rockim_f2` gelé (fichier GELE.md, plus aucun build) à la fin des campagnes en cours ; worktrees p2-p4 archivés en bundle et
   laissés en place.
3. **A6** — « fusionne tout dans g0 sans rien toucher dans les dossiers respectifs » : `insertion-pointe` porté en clés opt-in (contenu
   orphelin), les 87 lignes manquantes de `joint-handoff` triées et portées, `dif-intrinseque` = ses 2 scripts Python ; f2 et worktrees intacts.
4. **E8** — **Signorini pour la lame fem3d** (le choix physiquement et numériquement pertinent) ET protocole de coupe T2 refait (bloc
   60 × 20 mm, Lysmer avec ressorts, passe 1-1,5 mm, images toutes les 50 µs, témoins) ; développement dans g0.
5. **E9** — dans g0, `dpr` = le noyau corrigé, **sans alias ni variante caduque** ; principe étendu : croissance par addition à l'intérieur
   d'une version, **élagage permis à une version majeure** (CHANGELOG des clés retirées, erreurs nommées, arbres gelés conservés) ; nouvelle
   ancre de bit-identité à la naissance de g0.
6. **G2** — validation Abaqus CDP : trois états, **PASS** sous la tolérance proposée (5 % pic et travail, 10 % enfoncement et raideur de
   décharge), **DOUTE** entre la tolérance et 10 % (remonté à Fernando), **FAIL** au-delà de 10 % ; le volume broyé au pic reste
   **indicatif** (zone de doute jusqu'à 20 %) sauf avis contraire.
7. **I4** — second humain : **en suspens** (personne disponible pour l'instant) ; walkthroughs préparés et exercice « session neuve » joué par
   un agent sans contexte.
8. **K.1** — les campagnes **ne sont pas arrêtées** pendant la vague 1 (gardes w20/w21 et bit-identité générique en place ; rejeu si g0
   révèle un défaut sur un chemin exercé, identifié par l'empreinte).
9. **D1a** — `build_f2.cmd` : **sans objet** avec g0 (f2 gelé le conserve, g0 naît avec CMake).

Points ouverts hors plan : règle de verdict de la matrice (A 20 % vs B 10 %), critère d'érosion du modèle figé (max_i par défaut),
protocole T2 refait, marches 3-4 de l'hétérogénéité, bayésien bloc 1.
