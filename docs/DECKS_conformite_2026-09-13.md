# Decks de conformité (T4, campagne de correction du 13/09 après-midi) — à valider par Fernando

Cadrage : `docs/CAMPAGNE_correction_2026-09-13.md` §T4. Motifs : COMPLEMENT §2 (« commencer par St Anne »,
impact 3D 2025 sans pulvérisation), DIAGNOSTIC §4 (`meanTensionCapFactor = 0` dans tout essai de conformité),
ECARTS §5 (cap caché actif dans tous les runs, étiquette rn/rs à corriger). **Aucun deck n'a été lancé au-delà
de 2 µs** (règle 6) ; tout ce qui suit est soit lu dans les sources, soit mesuré par le lancement de fumée
`tools/deck_smoke.py` (T = 2 µs, 1 trame, 4 fils, `rockim_g1y16.exe`, journaux dans `tests_f2/campagne13/T4/`).

## 1. Les huit decks

| deck | base | ce qui diffère du témoin A (`yang2026_bench_s25_v3P_300.cfg`), liste exacte | maillage | T | exe requis |
|---|---|---|---|---|---|
| `configs/stanne2025_bench_s25_visc0.cfg` | témoin A, matériau St Anne | `rho 2731`, `E 57e9`, `nu 0.31`, `ft 7.0e6`, `cohesion 18.8e6`, `frictionDeg 45` (tanφ 1,0), `Gf 12`, `gfShearFactor 66.667` (800/12), `contactMu.rock 0.6`, `jointPenaltyLength edge`, `jointPenaltyFactor 26.316` (p0/(2E), p0 = 3 000 GPa ARMA 24-0952), `bulkDamage off` (+ `bulkDamagePhase/Delta0/DeltaF/Dmax` et `contactDamageCoupling` retirées : le code exige la seconde avec `bulkDamage = yang`, `Fdem3dSolver.cpp:1353`), `meanTensionCapFactor 0`, `groupVel.piston 0 0 -10.66`, `T 4.5e-4`, `frames 18`, `meshFile` T3 | `meshes/impact_yang_train1_rock25_hxt.msh` (T3, **absent au moment du test**) | 450 µs | g1y16 suffit |
| `configs/stanne2025_bench_s25_visc.cfg` | `_visc0` (généré par `tools/make_conformity_decks.py`) | idem + `bulkViscosity 2000` (**exploratoire** : ARMA « Mass Damping Coefficient 4000 » lu comme η D de Guo éq. 2.6 = 2 μ D, convention du banc C ; unité et opérateur non publiés) | idem | 450 µs | g1y16 suffit |
| `configs/yang2026_bench_s25_v4_A.cfg` | témoin A (v3P_300) | `meanTensionCapFactor 0`, `jointBreakModeRef slipRef` | `impact_yang_s2.5_pose.msh` (celui de la série v3, à dessein) | 300 µs | **g1y17** (S1) |
| `configs/yang2026_bench_s25_v4_B.cfg` | banc B (`_solidity.cfg`) | `jointPenaltyFactor 25`, `jointShearUnload solidity`, `jointSecantRatchet off`, `jointPenaltyLength edge`, `potPenaltyFactor 0.25`, `gcBirth penalty` + les deux clés v4 | idem | 300 µs | g1y17 |
| `configs/yang2026_bench_s25_v4_B1.cfg` | banc B1 (`_law.cfg`) | `jointShearUnload solidity`, `jointSecantRatchet off` + les deux clés v4 | idem | 300 µs | g1y17 |
| `configs/yang2026_bench_s25_v4_B2.cfg` | banc B2 (`_pen.cfg`) | `jointPenaltyFactor 25`, `jointPenaltyLength edge` + les deux clés v4 | idem | 300 µs | g1y17 |
| `configs/yang2026_bench_s25_v4_C.cfg` | banc C (`_solidity_visc.cfg`) | banc B + `bulkViscosity 2000` + les deux clés v4 | idem | 300 µs | g1y17 |
| `configs/yang2026_bench_s25_v4_D.cfg` | banc D (`_v3P_contact.cfg`) | `potPenaltyFactor 0.25`, `gcBirth penalty` + les deux clés v4 | idem | 300 µs | g1y17 |

Les listes « ce qui diffère » des v4 sont **calculées** par `tools/make_conformity_decks.py` (diff des paires
clé = valeur du deck source contre le témoin, commentaires ignorés) et recopiées dans l'en-tête de chaque deck ;
le corps du deck source est recopié octet pour octet, les deux clés sont ajoutées en fin de deck. Le fichier
`stanne2025_bench_s25.cfg` sans suffixe n'existe pas : l'amortissement étant inconnu, il n'y a pas de variante
« par défaut » (cadrage : deux variantes `_visc0` et `_visc`).

Choix faits sans données publiées (à contester) : la loi de joint et le contact du témoin A pour St Anne (la loi
`solidity` et le contact p0/200 sont mesurés à part par B/B1/B2/D ; si Fernando préfère le paquet Solidity pour
St Anne, ajouter `jointShearUnload = solidity`, `jointSecantRatchet = off`, `potPenaltyFactor = 0.263`
(15 GPa / 57 GPa), `gcBirth = penalty`) ; le DIF 2025 = `strainRateDIF yang`, `difExpT 0.17`, `difExpS 0.07`,
`continuous`, `none` comme le témoin (le journal confirme « ft et Gf reçoivent DIF_traction, cohésion et GfII
DIF_compression, l'angle de frottement est inchangé » — c'est la règle 2025 du COMPLEMENT §2) ; la vitesse
10,66 m/s (cas d'analyse des fig. 14-15 de Yang 2025 ; autres cas St Anne publiés : 2,475 ; 6,235 ; 9,41 ;
11,84 m/s) ; T = 450 µs (arrêt des radiales à 388 µs, fig. 14). La série v4 garde le maillage s = 2,5 historique
(piston 0,777 kg) pour que la seule différence v3 → v4 soit les deux clés ; la série à train figé est T3.

## 2. Coût estimé (mesuré au lancement de 2 µs, pas de temps réel de chaque deck)

Base : `results/yang_bench_s25_v3P.log` — 70 733 pas pour 200 µs en 5 498 s à 14 fils sur la machine partagée,
soit **77,7 ms/pas** ; « ~20 min seul » du cadrage → **17 ms/pas** (estimation, non mesurée).

| deck | dt mesuré (ns) | pas pour T | coût partagé 14 fils | coût seul (estim.) |
|---|---|---|---|---|
| `stanne2025_bench_s25_visc0` / `_visc` (sur le maillage s = 2,5 historique, voir §4) | 4,165 | 108 052 (450 µs) | 8 399 s ≈ 2,3 h | 1 837 s ≈ 31 min |
| `v4_A`, `v4_B1`, `v4_D` (pénalité 20 sur hmin_) | 2,828 (= v3P, 2,82757e-9) | 106 098 (300 µs) | 8 247 s ≈ 2,3 h | 1 804 s ≈ 30 min |
| `v4_B`, `v4_B2`, `v4_C` (pénalité `edge` 25) | 4,079 | 73 556 (300 µs) | 5 718 s ≈ 1,6 h | 1 250 s ≈ 21 min |

Ces coûts supposent le même ms/pas que le v3P (10 563 tets, 19 623 joints) ; le maillage T3 à train figé n'a pas le
même nombre de tétras (roche SR 2,5 + train s = 1) — à remesurer par `tools/deck_smoke.py` dès qu'il existe.
Total de la série v4 : ≈ 42 000 s (11,6 h) partagé, ≈ 2,5 h seul ; les deux St Anne : ≈ 4,7 h partagé, ≈ 1 h seul.

## 3. Critères de succès

- **St Anne** (Yang 2025 fig. 14 à 10,66 m/s ; Table 3 pour les ordres de grandeur des sept critères, cas non
  précisé dans le texte extrait) : zone broyée sous l'insert avant 130 µs ; radiales, latérales et **médiane**
  de 130 µs à la fin de la charge (254 µs) ; arrêt des médianes à 291 µs, des radiales à 388 µs ; longueur de
  fissure radiale 20-25 mm (`tools/crack_paths.py`, T1 : composante connectée au noyau) ; rayon de cratère
  10-12 mm ; vitesse d'indentation 9,4-9,9 m/s et rebond 6,9-7,1 m/s (**pentes** du déplacement, T2) ;
  enfoncement ≈ 1,5 mm. **Échec** = réseau diffus sans radiale connectée au-delà de ≈ 12 mm, ou bit non freiné.
  `_visc` contre `_visc0` : mesure l'effet de l'amortissement inconnu ; si les deux échouent de la même façon, le
  moteur joint/contact est en cause indépendamment de l'amortissement.
- **Série v4** : le critère du banc source (plateau de réaction ≥ 45 kN avant p = 0,5 mm, v_bit max ≤ 6,3 m/s,
  nPulv > 0, traces de joints au-delà de r = 9 mm) ; en plus (a) **compteur d'écrêtage** du journal (S1) = 0
  attendu — s'il est > 0, le cap changeait la physique de la v3 correspondante ; (b) proportions
  traction/cisaillement à la rupture à comparer à la v3 (`jointBreakModeRef = slipRef` ne change que
  l'étiquette : `history.csv` doit être identique à la v3 si le cap n'a jamais mordu — c'est le contrôle de
  la bit-identité de S1 sur un vrai deck).

## 4. Ce qui a été vérifié (lancements de 2 µs, `rockim_g1y16.exe`, 4 fils, `tests_f2/campagne13/T4/`)

| fumée | résultat |
|---|---|
| `smoke_g1y16_meshs25.md` — St Anne `_visc0` et `_visc`, **meshFile surchargé** vers `impact_yang_s2.5_pose.msh` (le maillage T3 n'existait pas) | 2/2 démarrent et finissent, code 0 ; dt 4,16468 ns, 481 pas ; 81 clés du deck consommées + 1 non lue à l'init mais légitime (`budgetAbortPct`, lue après init) ; journal : `jointPenaltyLength = edge`, « pénalité de joint 26.316 E/h », frottement de contact rock/steel/carbide = 0,6, DIF a_t 0,17 / a_s 0,07 sur ft-Gf / c-GfII, `groupVel.piston = -10.66`, `bulkDamage` absent du résumé (rien à ventiler) ; résidu B4 −8,0e-14 J (`_visc0`), +2,0e-13 J (`_visc`) ; masses **du maillage historique** : piston 0,7767 kg, bit 1,0802, insert 0,0631, roche 20,017 kg (= 19,247 × 2731/2626) ; KE piston 44,13 J (½ · 0,7767 · 10,66²). Avertissements : potKt_ hors budget de dt (le même que le témoin) ; `_visc` en a un de plus : « la viscosité est armée et entre dans sigG » (`viscousInInsertion = 1`, défaut, identique au banc C) |
| `smoke_g1y16_refus.md` — les six v4 **tels qu'écrits**, g1y16 | 6/6 **refusés** avec exactement **une** clé chacun : `cle 'jointBreakModeRef' (ligne N du deck) inconnue de rockim` — preuve que toutes les autres clés sont valides et que la garde `unknownKeys = error` refuse bien un deck v4 sur un binaire sans S1 (le refus est le comportement voulu). Le dt est imprimé avant le refus : 2,828 ns (A, B1, D) / 4,079 ns (B, B2, C) |
| `smoke_g1y16_sansS1.md` — les six v4 avec la ligne `jointBreakModeRef` **retirée** (`--drop-key`) | 6/6 démarrent et finissent, code 0 ; `meanTensionCapFactor = 0` consommée du deck (`config_effective.cfg`) ; résidus B4 entre −6,7e-14 et +5,8e-14 J ; mêmes masses ; dt du v4_A = celui du v3P (2,82757e-9 s) |

Contrôle sur les fichiers : 0 caractère de code < 32 hors tab/LF/CR dans les 8 decks et les 2 outils ;
`python tools/make_conformity_decks.py --check` → 7/7 identiques (génération idempotente).

## 5. Ce qui n'a PAS été fait, et pourquoi

- **Le maillage T3** `meshes/impact_yang_train1_rock25_hxt.msh` n'existait pas pendant T4 (agent T3 en
  parallèle) : les decks St Anne le référencent, la fumée a tourné sur le maillage s = 2,5 historique. À refaire
  sans `--mesh` dès qu'il existe (`python tools/deck_smoke.py --exe rockim_g1y17.exe --out tests_f2/campagne13/T4
  --tag g1y17 configs/stanne2025_bench_s25_visc0.cfg configs/stanne2025_bench_s25_visc.cfg`) ; les masses du
  train doivent alors valoir celles de la série s = 1 (piston 1,057 kg, bit 1,288 kg) et non 0,777 / 1,080.
- **Les v4 n'ont pas été démarrés avec `jointBreakModeRef`** : la clé est celle de S1, absente de g1y16 et du
  source au moment du test (`grep jointBreakModeRef src/` vide). L'agent Build rejoue la fumée avec g1y17
  (même commande, `--tag g1y17`, sans `--drop-key`) : attendu 6/6 démarrent.
- **Aucun run complet** (règle 6) : les critères du §3 ne sont pas mesurés.
- La conversion `jointPenaltyFactor = p0/(2E)` (26,316) reste **conditionnelle** (COMPLEMENT §3) : le banc S3 doit
  la vérifier sur la courbe σ(dn) d'un joint unique avant de lire le St Anne comme une réplique de leur pénalité.
- `viscousInInsertion` n'a pas été posé à 0 dans `_visc` (aligné sur le banc C, qui porte le même avertissement).

## 6. Fichiers

Decks : `configs/stanne2025_bench_s25_visc0.cfg` (écrit à la main, source), `configs/stanne2025_bench_s25_visc.cfg`
et `configs/yang2026_bench_s25_v4_{A,B,B1,B2,C,D}.cfg` (générés). Outils : `tools/make_conformity_decks.py`
(génération + `--check`), `tools/deck_smoke.py` (fumée avec surcharge T/frames/meshFile, retrait d'une clé,
tableau de coût). Journaux : `tests_f2/campagne13/T4/smoke_*.md` et un dossier par run (journal, deck de
travail, `config_effective.cfg`, `history.csv` ; les `.log` et `.vtu` sont ignorés par git, `*.log` et `*.vtu`
du `.gitignore` — les résumés `smoke_*.md` portent les chiffres). Documentation : §9 de `DOCUMENTATION_rockim.md`
(outils) et `CHANGELOG.md` ; la ligne `meanTensionCapFactor` (défaut réel 3, pas « 0 = off » ; compteur) a été
corrigée par S1 pendant T4, T4 n'y touche pas. Aucune clé de solveur ajoutée : registre `tools/keys_by_mode.json`
inchangé par T4.
