# Correction du contact par potentiel de rockim : implémentation

Date : 2026-10-07. Dépôt `/home/user/rockim`, branche `claude/adoring-bardeen-qc4ovq`, rien de commité.
Binaire de la nouvelle version : `/home/user/rockim/build_contact/rockim` (cmake Release). Binaire de
référence (HEAD 5fae24b, sources extraites par `git archive`) : `correction_contact/ref_bld/rockim`.
Sources : `docs/rapport_guide/ENQUETE_CONTACT.md` (numéros des correctifs ci-dessous = ceux de la
demande), `correction_contact/LITTERATURE.md`.

Fichiers modifiés : `src/FdemSolver.cpp`, `include/rockim/FdemSolver.hpp`,
`include/rockim/PotentialContact.hpp`, `src/Fdem3dSolver.cpp`, `include/rockim/Fdem3dSolver.hpp`,
`src/main.cpp` (fins de ligne CRLF conservées), `tools/verify_suite.py`, `DOCUMENTATION_rockim.md`,
`CHANGELOG.md`, registre des clés régénéré (`tools/gen_keys_by_mode.py` → `include/rockim/KeysByMode.hpp`,
`tools/keys_by_mode.json` ; séparateurs de chemin Windows remis à la main pour garder un diff minimal).

## Résumé des clés

| correctif | clé | valeurs | défaut | change une trajectoire ? |
|---|---|---|---|---|
| 1 | `contactCandidates` | `active` \| `vertex` | `active` | oui, opt-in (2D et 3D) |
| 2 | `gcBirth` | `ramp` \| `penalty` \| `offset` (nouveau) | `ramp` | oui, opt-in (2D) |
| 2 | (avertissement si `gcBirthTau < dt` sous `ramp`, ligne de bilan) | — | actif | non (stdout seulement) |
| 3 | (correction leapfrog des rouleaux) | — | actif | non (lignes de bilan seulement) |
| 4 | `potForceExact` | `false` \| `true` | `false` | oui, opt-in (2D) |
| 5 | (signe des commentaires / message) | — | — | non (texte) |
| 6 | `selftest-potcontact2d` + 3 repères de suite | — | — | — |

## Correctif 3 (bilan des rouleaux), actif par défaut

Identité du saute-mouton par ddl : ½M(v⁺)² − ½M(v⁻)² = dt·F·v⁻ + (dt·F)²/2M. Sur le ddl x d'un nœud
ROLLERX, la vitesse est remise à 0 après le kick : la part F_x²dt²/2M était inscrite dans `biasW_` sans
jamais atteindre l'énergie cinétique. La correction leapfrog se calcule désormais sur la force CONTRAINTE
(F_x = 0 sur le rouleau), et l'amortisseur Lysmer x n'est plus compté sur ce ddl. F et f_ ne sont pas
modifiés (trajectoire et réactions intactes). Exact dès que v_x = 0 sur le rouleau (vrai dès le 2e pas).

- `src/FdemSolver.cpp:9283-9296` (chemin groupé, copie `Fb` avec `Fb.x() = 0`), `:9305-9308`
  (amortisseur cX), `:9409-9416` (chemin par nœud : `fy²` seul), `:9425-9426` (cAbsX).
- 3D : déjà juste (les axes imposés `uMask_` sont exclus du biais et de l'amortisseur), pas de ROLLERX en
  3D. Rien à faire.
- Les autres vitesses imposées (FIXED, DRIVEX, PRESCRIBED) ne passent pas par la ligne de biais : rien à
  retrancher. Observation non corrigée : les ddl LIBRES de DRIVEX et de `gripFree` n'ajoutent pas leur
  terme (dt·F)²/2M (biais manquant, pas en trop) ; hors du périmètre demandé.

Preuve (stdout complet, binaire ref contre nouveau, deck smoke de l'enquête T = 0,03 s, 2 fils) : seules
changent les lignes `integration` (+136 005 → +0,068 J/m) et `residu` (−136 005 J/m [CHECK] →
−3,27e-11 J/m [OK]), plus l'avertissement et le message du correctif 5 ; `history.csv`, `frames.csv`,
VTU et CSV de fin identiques au bit près (`diff -rq`). Sur `tip16_mu0` T = 0,12 s : résidu −527 664 →
−1,38e-7 J/m, trajectoire identique (`diff -rq` des sorties : seul `config_effective.cfg` diffère, par
les deux nouvelles clés au défaut).

## Correctif 5 (signe du biais du compteur)

Σ f·v_old·dt = ΔKE − Σ|f|²dt²/2m : sur un canal conservatif le biais est NÉGATIF. Corrigés :
`src/FdemSolver.cpp:10073-10083` (commentaire + message imprimé du résumé, qui disait « petit résidu
positif = biais ») et `:10809-10819` (commentaire du selftest-potential2d, « biais POSITIF »). Même
correction dans la ligne `contact` de `DOCUMENTATION_rockim.md` (§5.4). Changement de texte seulement.

## Correctif 4 (forces nodales gradient exact), `potForceExact` (défaut false)

E = p∫_S(φ_A + φ_B)dA, S = A∩B, φ = 3 min λ. La répartition de Munjiza vaut
f_a = −∂E/∂x_a + p·I_A·∇λ_a (nœuds de A), I_A = ∫_S φ_A dA, et symétriquement pour B (LITTERATURE §2.2,
dérivation par transport de Reynolds). Avec la clé : f_a ← f_a − p·I_A·∇λ_a, f_b ← f_b − p·I_B·∇λ_b.
Terme auto-équilibré (Σ∇λ = 0, moment nul) : résultante, moment et 3e loi inchangés.

- `include/rockim/PotentialContact.hpp:210-389` : `clipHalfPlane`, `AffBary`, `overlapIntegrals`
  (`:275`, découpe de S en 9 régions argmin λ_A = i, argmin λ_B = j ; ∫3λ_i = 3λ_i(centroïde)·aire),
  `pairEnergy` (`:332`), `exactGradientForces` (`:344`), `birthOffsetScale` (`:374`).
- `src/FdemSolver.cpp:7361-7365` (appel après `pairForce`, avant le facteur de naissance `sc`), lecture
  de clé `:1514-1527`, membre `potExact_` (`FdemSolver.hpp:977-995`).
- Piège trouvé et corrigé en cours de route : en coordonnées absolues (|X| ~ 100 m sur le tunnel), λ = g·X + c
  et l'aire des sous-polygones perdent tous leurs chiffres pour un élément fin : I_A = −1,46 pour une aire
  de 2e-13 m². `overlapIntegrals` travaille désormais dans un repère local (origine au 1er sommet de A),
  ce qui ramène le minimum de (I_A + I_B)/h² à −4e-21 sur 200 000 tirages (P0b du selftest). L'énergie
  rendue est en outre bornée à ≥ 0 (arrondi de recouvrements rasants). Les premières campagnes tunnel
  faites avant ce correctif ont été rejouées.
- 3D NON porté : il faudrait I_A = ∫_V φ_A dV par découpe du polyèdre de recouvrement selon les 12 plans
  de médiane (16 régions), ∇λ_a = −𝒜_a ν_a/(3V). Exception documentée (constitution III).

Validation (`selftest-potcontact2d`, configuration de `loop2.cpp`, critères posés avant calcul) :
- P1 gradient : |f − (−dE/dq)|/|f|max = 0,157 (Munjiza, défaut attendu > 1e-2) ; 2,2e-10 (exact, < 1e-6).
- P2 boucle fermée sur ddl nodaux (66 paires de ddl, r = 0,02, 4 000 points) : max |W| = 2,85e-4
  (Munjiza, attendu > 1e-5) ; 5,1e-17 (exact, < 1e-12). C'est le test `loop.cpp/loop2.cpp` de l'enquête.

## Correctif 1 (candidats), `contactCandidates = vertex` (défaut active)

Littérature (LITTERATURE §1) : Y/Y-Geo voient toutes les paires d'éléments ; Guo 2014 (éqs. 2.39-2.46)
active à la rupture d'un joint l'anneau de TOUS les éléments contenant ses 6 nœuds ; Fukuda et al. 2021
(semi-ACAA) montre que les contacts manquants pendant l'adoucissement en cisaillement des cohésifs vivants
produisent une fragmentation parasite, et active dès le début de l'adoucissement. Retenu :

  𝓔_cand = 𝓔_actif ∪ Ring(sommets des éléments actifs) ∪ Ring(sommets de tout joint INSÉRÉ),
  Ring(v) = {E | v ∈ E} sur les sommets d'origine `vOf_` (topologie initiale), exclusion inchangée des
  paires reliées par un joint VIVANT.

Activation dès l'insertion (D* = 0, l'équivalent du « seuil ≈ 1 » de Fukuda pour un CZM extrinsèque).
Recalcul à chaque appel par balayage de `jt_` (résultat identique à l'ajout événementiel proposé par la
note, puisque l'insertion est définitive). Le contact est évalué même quand `act_` est vide. Le retrait
de la grille « % 8 » de `deadList_` n'a pas été ajouté : sous `vertex`, les éléments d'un joint mort sont
déjà candidats depuis son insertion, la latence de `deadList_` ne joue plus pour le potentiel (et
`gcSurfaceRefresh = eager` existe pour le reste).

- `src/FdemSolver.cpp:7212-7251` (bloc 1b), `:7564-7566` (`act_` vide), clé `:1528-1548`, compteur au
  bilan `:10165-10184`, membres `FdemSolver.hpp:977-995`.
- 3D, à l'identique : `src/Fdem3dSolver.cpp:1826-1844` (clé), `:5935-5967` (bloc 1b, 4 sommets par tet, 6
  nœuds par joint), `:6413` (`act_` vide), `:7690-7692` (compteur), `Fdem3dSolver.hpp:1365-1368`.
- Aucune paire reliée par un joint vivant n'est ajoutée : l'exclusion `jointOfPair_`/`!dead` est appliquée
  après la recherche, sans changement.
- P0 du selftest : 2 000 paires à sommet commun sans recouvrement (secteurs disjoints, tournées, échelles
  1e-3 à 10, |X| ~ 100) → aire de recouvrement exactement 0 (aucune force parasite). Arêtes colinéaires
  (deux triangles qui se touchent le long d'une droite, cas absent d'un éventail conforme) : 2,8e-7 de
  l'aire, comportement ordinaire de `pairForce`, lu sans critère.
- Charge nulle : 2D (`verify_fdem_tension`, 3 clés) 0 casse, travail de contact −1,5e-23 J/m ; 3D
  (`verify_fdem3d_tension`, pullV = 1e-12, T = 1e-4) 0 casse, travail de contact 0 J exactement, sorties
  au défaut identiques entre binaires.
- Coût : tunnel `tip16_mu0` (insertion adaptative), séquentiel 2 fils, machine calme : 80 s (défaut),
  137 s (`vertex`, ×1,7), 177 s (les trois clés, ×2,2). Sous insertion intrinsèque tous les joints sont
  insérés donc tous les éléments deviennent candidats : ×2 à ×4 en 2D (traction 89 s → 200-330 s sous
  charge), ×13 sur la charge nulle 3D (46 s → 607 s).

## Correctif 2 (naissance), `gcBirth = offset` (défaut ramp)

Ni la rampe ni la pénalité ne sont neutres : faire passer le facteur sc de 0 à 1 à géométrie fixée crée
∫e·dsc = E(S0) sans travail (la note confirme ; « re-scaling Y3D » et « contact birth » de Munjiza non
trouvés dans les sources accessibles). La note recommande l'offset avec cliquet de LS-DYNA (`IGNORE = 1`,
F = k(d − d0), d0 mis à jour quand la pénétration diminue) transposé au potentiel. Formulation retenue,
en énergie (fonction de E seule) :

  à la naissance : e0 = E(S0) ;
  U(E) = (E − e0)²/E si E > e0, 0 sinon ; force = U′(E)·(−∇E), U′(E) = 1 − (e0/E)² ;
  cliquet : si E ≤ e0, e0 ← max(0, E) (n'agit que sur la branche U = 0).

Propriétés : U(S0) = 0 (aucune énergie créée) ; U′(e0) = 0 (force continue, pas de saut au franchissement,
contrairement à la « version simple » max(0, e − e0) de la note) ; U′ → 1 quand E ≫ e0 ; e0 = 0 redonne
Munjiza ; le champ est un gradient exact avec `potForceExact = true`, puisque U ne dépend que de E. Écart à
la note : celle-ci propose un décalage de NIVEAU c0 (E_eff = p∫(φ_A+φ_B−c0)₊), qui a les mêmes propriétés
mais exige un clip de plus par paire et une force nodale à dériver pour le domaine {φ_A+φ_B > c0} ; la forme
en E seule réutilise le gradient exact du correctif 4. Je l'ai donc préférée ; les deux sont défendables.

- `include/rockim/PotentialContact.hpp:374` (`birthOffsetScale`) ; `src/FdemSolver.cpp:7370-7394`
  (compteur pur de naissance pour TOUS les modes, puis branche `offset`) ; clé et refus de `gcBirthTau`
  avec `offset` `:1433-1477` ; champ `PotHist::e0` (`FdemSolver.hpp:1298-1300`).
- AVERTISSEMENT (aucun changement) quand `gcBirthTau < dt` sous `ramp` (`:1460-1477`) : sur le tunnel
  « gcBirthTau = 1e-06 s < dt = 5.83025e-06 s -> relax = 0.00293734 : la rampe degenere en MARCHE ».
- Nouvelle ligne de bilan, tous modes : « contact, naissances : N paires, energie de recouvrement a la
  naissance Σ E(S0) J/m (materialisee sans travail par la rampe | ... penalty | neutralisee par offset) ».
  Sur `tip16_mu0` défaut : 3 689 paires, 1,45186e6 J/m, soit exactement le ΣdE_sc de l'enquête (1,452e6,
  instrumentation indépendante).
- 3D non porté (même raison que le correctif 4 : E exacte d'un recouvrement tet-tet à intégrer).

Validation P3 du selftest (deux blocs CST élastiques SANS cohésion, une copie du sommet commun chacun,
nés en recouvrement profond S0 = coin de 15° à 27°, E = 10, ν = 0,25, p = 10, pressés l'un vers l'autre à
±0,1, 30 000 pas, H = KE (synchronisée) + U_el + U_contact, critères fixés avant) :
- `ramp` avec τ = dt/5,8 (rapport du tunnel) et forces exactes : énergie créée H_fin − H0 = 1,00013 E(S0)
  (E(S0) = 0,544, KE0 = 0,0038) : le défaut est reproduit (critère > 0,5).
- `offset` + forces exactes : max_t |H − H0|/H0 = 2,2e-6 (critère < 1e-3) ; H_fin − H0 = 7,9e-9.
- `offset` + forces de Munjiza (lu) : 0,12 de dérive relative à KE0 — le terme non gradient pèse ici
  parce que le recouvrement est profond. Les deux correctifs vont ensemble.

## Correctif 6 (tests)

- `rockim selftest-potcontact2d` (`src/FdemSolver.cpp:11147-11460`, déclaration
  `PotentialContact.hpp:939-943`, `src/main.cpp:56` et `:192-198`) : P0, P0b, P1, P2, P3 ci-dessus,
  exit 0 = PASS, 1 s.
- `tools/verify_suite.py` : regex `potx_ramp`, `potx_off`, `birthe` ; repère `selftest_potcontact2d`
  (tier fast : pass_tag, potx_ramp = 1,00013 ± 1e-3, potx_off = 0 ± 1e-3) ; `zeroload_contactfix_2d`
  (tier full : 0 casse, dampWork ≤ 0, gcwork = 0 ± 1e-12 → −1,5e-23) ; `contactfix_tension_2d` (tier full,
  témoin `gcbirth_ramp_2d` −1,24549 % : err_pct −1,55359 % ± 0,01, 24 morts, résidu B4 9e-14). Les deux
  repères full passent avec le binaire final (la référence −1,46909 % mesurée avant la correction du repère
  local a été recalée à −1,55359 %, c'est une clé opt-in).

## Suite rapide et bit-identité

`python3 tools/verify_suite.py --tier fast`, OMP = 1, Linux g++ :
- référence (HEAD) : 50/51, échec connu `t1_toolcontact_penalty` (toolinj 2,98519 contre 4,42633) ;
- nouvelle version (binaire final) : 51/52 = les 51 mêmes repères (50 PASS, même échec connu, valeurs
  identiques) + `selftest_potcontact2d` PASS. Les valeurs mesurées des 51 repères sont identiques
  (comparaison des JSON) et les 1 098 fichiers de sortie des 51 runs sont identiques au bit près
  (`diff -rq` des dossiers de suite), hors `config_effective.cfg` qui liste les nouvelles clés au défaut
  (`contactCandidates = active`, `potForceExact = false`).
- Deux calculs courts comparés en entier, stdout et sorties : smoke tunnel T = 0,03 s et `tip16_mu0`
  T = 0,12 s (2 fils) : trajectoires identiques au bit près ; seules changent les lignes de bilan (correctif
  3), le message (correctif 5), l'avertissement et la ligne « naissances » (correctif 2). 3D charge nulle :
  `history.csv` identique entre binaires.

## Résultats sur le tunnel (`tip16_mu0`, `tunnel_hs_red.msh`, T = 0,12 s, 2 fils, binaire final)

| p | clés | gcWork (J/m) | naissances, Σ E(S0) (J/m) | cohésif (J/m) | résidu B4 (J/m) |
|---|---|---|---|---|---|
| 1e10 | défaut | +1,309e6 | 3 689 ; 1,452e6 créés | 392 838 | −1,4e-7 |
| 1e10 | potForceExact | +1,534e6 | 3 428 ; 1,668e6 créés | 405 605 | −7,1e-8 |
| 1e10 | offset | −2,79e4 | 2 720 ; 2,345e6 neutralisés | 374 680 | −1,8e-8 |
| 1e10 | offset + exact | −3,14e4 | 2 859 ; 1,396e6 neutralisés | 360 254 | −1,2e-8 |
| 1e10 | vertex | −8,81e4 | 15 364 ; 16 | 317 944 | 3,6e-8 |
| 1e10 | vertex + exact | −8,27e4 | 15 101 ; 19 | 308 294 | −4,8e-9 |
| 1e10 | vertex + offset | −7,55e4 | 15 267 ; 11 neutralisés | 312 016 | −2,0e-8 |
| 1e10 | les trois | −8,27e4 | 15 560 ; 15 neutralisés | 316 551 | −1,6e-8 |
| 1e9 | défaut | +1,322e5 | 3 106 ; 1,79e5 | 353 241 | −7,2e-9 |
| 1e9 | les trois | −6,44e4 | 17 516 ; 1,3 | 345 282 | 2,9e-11 |
| 1e11 | défaut | +7,131e6 | 4 065 ; 8,17e6 | 471 658 | 4,1e-7 |
| 1e11 | les trois | −8,03e4 | 13 973 ; 117 | 240 202 | −3,0e-8 |

Critères posés par l'enquête (§4-1) : gcWork ≤ 0 pour p = 1e9, 1e10, 1e11 avec les trois clés : TENU.
Énergie de naissance < 1 % de l'énergie cohésive : TENU (0,0004 %, 0,005 %, 0,05 %). `vertex` seul tarit
la source (16 J/m au lieu de 1,45e6) ; `offset` seul la neutralise sans la tarir (les recouvrements
profonds restent, mais ne créent rien) ; `potForceExact` seul ne change pas l'ordre de grandeur (le terme
non gradient pèse 0,3 % ici, comme mesuré par l'enquête). Le faciès change avec les clés (cohésif 317 kJ/m
contre 393 kJ/m, 4 350 joints rompus contre 4 889 à p = 1e10) : c'est attendu, une partie de la
fragmentation du défaut était entretenue par l'énergie créée. Le critère « éventail de 6 triangles avec
mort d'un joint » de l'enquête est couvert par P3 (deux blocs à sommet commun nés en recouvrement), pas
par un éventail complet avec joints.

Fichiers : decks et journaux `correction_contact/runs/` (`f_*.log` = binaire final, `smoke_*`, `t12_ref`
/ `t12_new`, `tm_*` = chronométrage), cas de suite `correction_contact/suite_new_cases/`, suites
`suite_ref*.{log,json}`, `suite_final*`, `suite_full_new*`, bancs autonomes `correction_contact/tests/`.

## Points ouverts

- 3D : `gcBirth = offset` et `potForceExact` à porter (intégrale exacte de φ sur le polyèdre).
- Comptage du travail au point milieu (v̄ = ½(v⁻ + v⁺), LITTERATURE §5), qui supprimerait `biasW_` et
  rendrait chaque poste lisible : non fait (changerait les postes affichés de tous les runs).
- Biais (dt·F)²/2M manquant sur les ddl libres de DRIVEX et `gripFree` (observation, non corrigé).
- Repère de suite pour le résidu B4 des rouleaux : aucun deck à rouleaux court dans la suite ; le smoke
  tunnel (30 s) pourrait en servir.

## Corrections après revue (REVUE.md, 2026-10-07)

Diff d'avant revue sauvegardé : `correction_contact/fix/avant_revue.diff`. Fichiers touchés en plus :
`src/Config.cpp` (M5). Rien de commité.

- **C1 (e0 périmé)** — `PotHist::lastEval` (`FdemSolver.hpp`), sous `offset` seulement : une paire non
  évaluée au pas précédent (pairForce faux, sortie des candidats ou de la boîte, élément inversé) RENAÎT
  à son retour, e0 = E courant (`FdemSolver.cpp`, branche `birthOffset_`). Compteurs `nPotRebirth_`,
  `potRebirthE_`, ligne de bilan « contact, renaissances ». Précision sur le scénario de la revue : un retour
  SOUS l'ancien e0 était déjà rattrapé par le cliquet (e0 ← E, identique à une renaissance) ; le vrai
  défaut était le retour en un pas AU-DESSUS de l'ancien e0 : force 1 − (e0/E)² d'emblée et
  U = (E − e0)²/E matérialisé sans travail. Test : P4 du selftest (ci-dessous), qui échoue sans le correctif
  (force 0,64 de la pleine au retour, 0 renaissance).
- **C2** — texte corrigé (`PotentialContact.hpp` au-dessus de `birthOffsetScale`, bannière d'init,
  DOCUMENTATION, CHANGELOG) : e0 = E(S1) > 0 pour toute paire, y compris par approche ; aucune paire n'est
  en Munjiza pur sous offset.
- **C3** — DOCUMENTATION (moniteur d'énergie + ligne `lateralRollers`) et CHANGELOG (titre de section
  « bit-identiques SAUF sous budgetAbortPct > 0 + lateralRollers »).
- **C4** — `python3 tools/gen_keys_by_mode.py` relancé ; séparateurs Windows remis dans `files_scanned`.
  `contactCandidates` : lecteurs « fdem fdem3d », retirée de `kTable` (plus propre à un mode) ;
  450 clés, 285 communes, 165 propres à un mode.
- **C5 (tests)** — sonde `PotContactProbe` (friend de `FdemSolver`, définie dans `FdemSolver.cpp`) :
  petit deck en dur (bande 8 × 8, h = 1, E = 10 donc potP_ = 10, contactMu = 0), init silencieuse, puis
  positions et vitesses IMPOSÉES sur les nœuds d'une paire et appel de `generalContact()` →
  `potentialContact()` du solveur.
  - P3 réécrit : mêmes blocs CST, forces de contact du SOLVEUR (deux éléments de bord sans joint ni sommet
    commun, déplacés hors du corps). Résultats identiques à l'ancienne réimplémentation : rampe
    1,00013 E(S0), offset + exact 2,2e-6, offset + Munjiza 0,121.
  - P4 (nouveau, C1) : cinématique imposée, force du solveur comparée à la force attendue recalculée
    indépendamment (pairForce, exactGradientForces, birthOffsetScale) : naissance (0), plus profond, cliquet
    (0), séparation en un pas, retour direct à 0,5 E1 (0 : renaissance), 0,7 E1 (> 0,05 de la pleine) ;
    écart 0, 1 naissance, 1 renaissance.
  - P5 (nouveau, complétude de `vertex` + C6) : éventail du sommet intérieur (4, 4), eA tourne autour du
    sommet jusqu'à recouvrir eB (40 pas). E(naissance)/E_plein : active 1,0 (née profonde, à la mort de
    joints de eA et eB) ; vertex + joint de l'éventail mort 0,0055 ; vertex + joint seulement endommagé
    (D = 0,5) 0,0055 ; vertex sans aucun joint endommagé : jamais candidate (déclencheur C6).
  - Mutations vérifiées (binaire reconstruit à chaque fois, puis restauré) : renaissance retirée → P4
    [ECHEC] ; `H.e0` laissé à 0 à la naissance → [FAIL] ; appel `exactGradientForces` du solveur retiré →
    P3 (0,121) et P4 (0,044) [ECHEC].
  - Bout en bout, tier fast : `contactfix_perc_{ramp,vertex,three}_2d` — `fdem_percussion.cfg` réduit
    (30 × 25, T = 60 µs, jeu nul, `absorbing = none`, `dampingLocal = 0`, `lateralRollers = true`), où 4-5
    joints meurent EN COMPRESSION. birthe : 2,91e-3 J/m (ramp, critère ≥ 2e-3), 8,88e-4 (vertex, ≤ 1,5e-3),
    8,91e-4 (les trois) ; gcwork −1,19 / −1,96 / −1,94 J/m (critère ≤ 0) ; résidu B4 à rouleaux −1,2e-13 /
    −1,3e-12 / −1,6e-12 J/m (critère 0 ± 1e-9 ; binaire HEAD sur le même deck : −1,73 J/m, trajectoire
    identique). C'est le repère des rouleaux demandé (résidu petit) et il échoue avec l'ancien bilan.
  - Refus (tier fast, init seule) : `refuse_contactcandidates_value`, `refuse_gcbirth_value`,
    `refuse_offset_tau`, `refuse_vertex_penalty`, `refuse_offset_penalty`, `refuse_potexact_penalty`,
    `refuse_contactcandidates_3d`, `refuse_potexact_3d` ; avertissement `warn_offset_sans_exact` ; M5 :
    `getb_casse_potexact`, `getb_valeur_inconnue`.
  - `verify_suite.py` : contrôles `max:`/`min:` (bornes, propriétés), `expect_error` (code ≠ 0 + message),
    `expect_out` (ligne exigée). La regex `birthe` sert désormais (repères `contactfix_perc_*`).
- **C6** — déclencheur de l'anneau (2D bloc 1b et 3D) : `!J.bonded && (adaptive_ || J.D > 0 || J.dead)`.
  En adaptatif : inchangé (joint inséré). En intrinsèque : joint endommagé ou mort. Bannières 2D/3D et
  message d'erreur mis à jour. Coût mesuré (OMP = 1, machine partagée avec le testeur, à lire comme ordre
  de grandeur) : petite percussion intrinsèque T = 30 µs (`contactfix_perc` à 30 µs) : défaut 3,4 s ;
  `vertex` 28,6 s avant (×8,5) → 7,7 s après (×2,3) ; à T = 60 µs : 5,3 s / 17,7 s / 31,8 s (défaut,
  vertex, les trois) ; charge nulle 3D (`zl3d_vtx`, T = 1e-4) 607 s → 370 s (défaut 46 s ; le reste vient
  des voisins par sommet des éléments ACTIFS — toute la peau extérieure — que C6 ne touche pas) ;
  `contactfix_tension_2d` (intrinsèque, full) 200-330 s → 58 s, `zeroload_contactfix_2d` → 67 s.
  Insertion adaptative (tunnel `tip16_mu0`) : déclencheur inchangé.
- **M1** — ligne de bilan : « materialisable sans travail par la rampe : en totalite a geometrie figee,
  quelle que soit gcBirthTau » ; commentaire de l'avertissement τ < dt (seul le cas le plus brutal est
  signalé ; au-delà la rampe étale sans annuler).
- **M2** — ligne `contact` de la DOCUMENTATION : la relève par aire ne protège pas une paire née en
  recouvrement ; renvoi vers `offset` et `vertex`.
- **M3** — avertissement à l'init si `gcBirth = offset` sans `potForceExact = true` (repère
  `warn_offset_sans_exact`).
- **M5** — `Config::getb` insensible à la casse : balayage des 3 918 `.cfg` de `/home/user` (hors
  `test_contact`) et des générateurs (`.py`, `.sh`, `.cmd`, `.ps1`) : toutes les valeurs des 48 clés lues
  par `getb` sont en minuscules (true 4 173, false 1 378, on 268, off 40, 1 69, 0 33) ; aucun deck changé.
  Valeur non reconnue : lue false comme avant, signalée une fois par clé sur stderr.
- **M6** — `AffBary::set` rend false (sans diviser) pour den ≤ 1e-300 ; `overlapIntegrals` et
  `exactGradientForces` sortent à vide.
- **M7** — commentaire de P3 et texte de la DOCUMENTATION : le recouvrement initial se détend d'abord
  (E 0,544 → 0, cliquet), puis le contact se réengage (Uc max ≈ 0,5 H0).
- Non traités (hors demande) : M4 (messages 3D), M8 (coût du compteur `pairEnergy` au défaut).

### Vérifications
- Build `build_contact` (Release), `selftest-potcontact2d` : [PASS] (P0-P5), 1,7 s.
- `python3 tools/verify_suite.py --exe build_contact/rockim --tier fast` (OMP = 1) : **65/66**, seul
  échec connu `t1_toolcontact_penalty` (2,98519, inchangé). Les 52 repères d'avant : valeurs mesurées
  identiques à `suite_final.json` (seul `selftest_potcontact2d` gagne la mesure `potx_reb`), et
  `diff -rq` des dossiers de sortie contre `suite_final/` : aucun fichier différent. 14 repères nouveaux,
  tous PASS. Journal : `fix/suite_fast.{log,json}`.
- Tier full : `zeroload_contactfix_2d` PASS (0 casse, gcwork +8,4e-25 J/m) ; `contactfix_tension_2d`
  recalé −1,55359 → **−2,07991 %** (C6 + C1, clés opt-in ; 24 morts, résidu −1,3e-13 J/m).
- Tunnel `tip16_mu0` (p = 1e10, T = 0,12 s, 2 fils, binaire corrigé ; `fix/f_*.log`) — critère gcWork ≤ 0
  toujours tenu avec la renaissance C1 : les trois clés −8,27e4 → −7,89e4 J/m (15 539 naissances, 13 J/m
  neutralisés ; 224 112 renaissances, 15 J/m neutralisés), résidu −1,7e-8 J/m, 4 350 casses ;
  offset + exact −3,14e4 → −3,04e4 J/m (21 057 renaissances, 24 J/m). Les chiffres du tableau ci-dessus
  sont donc ceux d'avant C1 pour les lignes `offset`.
