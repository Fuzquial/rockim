# SPEC — portage de la loi « FDEM hybride à insertion adaptative » dans rockim_g1

Contrat d'implémentation. Version de départ : copie conforme de `rockim_g0`
(`rockim_g1ref.exe`, 2 263 040 octets, identique à `rockim_fix.exe`).
Source : note de travail de septembre 2026, sections 1 à 3.
Audit de départ : `..\rockim_g0\AUDIT_loi_adaptative_2026-09-11.md`.

---

## 0. Règle absolue — croissance par addition (principe VIII)

**Aucun comportement existant ne change.** Toute capacité nouvelle arrive par une clé
opt-in dont le défaut reproduit *exactement* rockim_g0. Le critère de recette est la
**bit-identité** : `tools/bitid.py` doit rendre le même résultat qu'avec `rockim_g1ref.exe`
sur tous les decks de référence, clé absente.

Corollaires :
- jamais de `throw` nouveau sur un deck qui passait avant, **sauf** les trois gardes
  explicitement listées en §5 (elles refusent des clés aujourd'hui **inertes**, donc
  aucun deck correct ne les porte) ;
- jamais de changement de défaut, même « évidemment meilleur » (`capP0`, `jointXi`,
  `dampingLocal`, `erodeD` restent où ils sont — on les éteint **dans le deck**) ;
- toute clé nouvelle doit être ajoutée à `tools/keys_by_mode.json` **et**
  `include/rockim/KeysByMode.hpp` régénéré, sinon `KeyGuard` la refuse à l'exécution.
  **Ce fichier est géré centralement — ne le modifiez pas, déclarez vos clés dans votre
  rapport de fin.**

---

## 1. Fichiers partagés — déjà écrits, NE PAS MODIFIER

### `include/rockim/JointTsl.hpp` (nouveau)
Noyau physique des joints extrinsèques, fonctions **pures**, identiques 2D/3D :

| fonction | § de la note | rôle |
|---|---|---|
| `jtsl::phiInsert(tn, ts, tn0, ts0)` | 2.2 eq. 12 | critère elliptique, rend Φ_F |
| `jtsl::dif(rate, rate0, expo)` | 2.2 eq. 13 | DIF_k(x) = max[1, (x/ε̇₀)^a_k] |
| `jtsl::effSep(dn, ds, beta)` | 2.4 eq. 16 | δ_m |
| `jtsl::effTrac(tn, ts, beta)` | 2.4 eq. 16 | t_m |
| `jtsl::stampInsertion(tn, ts, tn0, ts0, GIc, GIIc, eta)` | 2.4 | rend `Stamp{tmIns, beta, Gc, dmF}` gelé |
| `jtsl::traction(S, dm, dmMax)` | 2.4 eq. 17 | charge + décharge sécante |
| `jtsl::damage(S, dmMax)` | 2.4 eq. 17 | D = δ_m^max/δ_m^f |
| `jtsl::split(tm, dm, dn, beta, tn, tsScale)` | 2.4 eq. 19 | partition t_n / t_s |
| `jtsl::shearCap(coh, mu, tn, D, mobilised)` | 2.4 | (1−D)·cohésion + **D**·μ⟨−t_n⟩ |
| `jtsl::viscous(eta, rate)` | 2.5 | option B |
| `jtsl::tetEdgeLength(V)` | 1.3 eq. 6 | h_e = (12V/√2)^(1/3) |
| `jtsl::compBandLimit(E, Gc, sigC0)` | 1.3 eq. 7 | borne de snap-back |

**Convention de signe : traction normale POSITIVE à l'ouverture** (section 2 de la note).

### `MatState` (déjà existant, `include/rockim/MatLaw.hpp:344-361`)
Porte déjà `wPlas`, `wDamT`, `wDamC`, `epvEq`, `D`, `Dc`, `kappa`, `ftScale`, `eroded`.
**Rien à ajouter pour l'instrumentation énergétique** : il suffit de les récolter.

### Contrat MatLaw ↔ solveurs — UN SEUL point de contact
La signature `stress(eps, s, dt, lc)` **ne change pas** (elle est virtuelle pure et
réimplémentée par une dizaine de classes dérivées). La longueur de bande de
**compression** passe par un champ nouveau de `MatState` :

```cpp
double lcComp = -1.0;   // §1.3 : longueur de bande de la COMPRESSION.
                        // <= 0  -> utiliser lc (comportement de rockim_g0).
                        // >  0  -> b_c = fc0 * lcComp / compGIIc.
```

C'est le **seul** champ à ajouter à `MatState`. Le solveur le renseigne une fois
à l'initialisation ; `MatLaw` s'en sert pour β_c uniquement, jamais pour la traction.

---

## 2. Clés — liste exhaustive et figée

### Matrice (lues dans `MatLaw::make`, `src/MatLaw.cpp`)

| clé | valeurs | défaut | § | effet |
|---|---|---|---|---|
| `bulkTensionDamage` | `on` \| `off` | `on` | 1.4 | `off` : le bloc de Rankine (`MatLaw.cpp:483-507`) ne s'exécute pas, `s.D` et `s.wDamT` restent 0, et les deux gardes de bande **en traction** (`kf <= 0.05*k0` et la garde de construction) sont levées. La contrainte nominale devient `σ = σ̄⁺ + (1−ω_c)σ̄⁻`, exactement l'eq. 8. |
| `mhForm` | `dp` \| `principal` | `dp` | 1.2 | `principal` : critère `f = σ̄₁ − σ̄₃ − B⟨σ̄₃⟩ⁿ − σ_c` en contraintes **principales**, retour plane/edge/apex. Exige `meridian = power`. |
| `dpFlowForm` | `dp` \| `mc` | `dp` | 1.2 | `mc` : potentiel `g = σ̄₁ − m_ψ σ̄₃`, `m_ψ = (1+sinψ)/(1−sinψ)`. |

### Solveurs FDEM — **noms et sémantique strictement identiques en 2D et en 3D**

| clé | valeurs | défaut | § | effet |
|---|---|---|---|---|
| `facetAverage` | `arith` \| `volume` | `arith` | 2.1 eq. 10 | `volume` : σ_F = (V⁺σ⁺ + V⁻σ⁻)/(V⁺+V⁻) au lieu de 0,5/0,5. S'applique à la contrainte **et** au taux. |
| `facetRate` | `scalar` \| `tensor` | `scalar` | 2.1 eq. 11 | `tensor` : stocke le tenseur taux par élément et forme `ε̇_eq = n·ε̇_F·n` et `γ̇ = 2‖ε̇_F·n − (n·ε̇_F·n)n‖`. Sans lui, DIF_t et DIF_s partagent le même scalaire. |
| `insertionCriterion` | `or` \| `elliptic` | `or` | 2.2 eq. 12 | `elliptic` : `jtsl::phiInsert(...) >= 1` remplace le OU logique. |
| `insertionHoldSteps` | entier ≥ 1 | `1` | 2.2 | n_h : exige Φ_F ≥ 1 sur n pas **consécutifs**. Compteur remis à 0 dès que Φ < 1. |
| `difExpT` | réel ≥ 0 | = exposant actuel | 2.2 eq. 13 | a_t, exposant du DIF de **traction**. |
| `difExpS` | réel ≥ 0 | **absente = chemin actuel (0,07)** | 2.2 eq. 13 | a_s, exposant du DIF de **cisaillement**. La note impose a_s < a_t ; **avertir** si a_s ≥ a_t, ne pas refuser. **Amendé après les lots** : le contrat écrivait « = `difExpT` », ce qui aurait changé le DIF de cisaillement d'un deck d'ancre (`fdem3d_kuru9_court`, `strainRateDIF = yang-fig2` où `difExpT` = 0,1707 alors que le cisaillement passe par le 0,07 de l'eq. 2 de Yang). Les trois lots ont, indépendamment, préféré le principe VIII — à raison. Le défaut qui en résulte (a_t = 0,1707, a_s = 0,07) satisfait en prime a_s < a_t. |
| `jointTSL` | `penalty` \| `camacho` | `penalty` | 2.4 | `camacho` : loi **initialement rigide**, t_m^ins = traction transmise, δ_f = 2G_C/t_m^ins, décharge sécante. Exige `insertion = adaptive`. |
| `jointTSLRise` | réel ≥ 0 | **`1e-3`** sous `camacho` | 2.4 | **Ajoutée en seconde passe (2026-09-11).** Branche ascendante courte de Papoulia-Sam-Vavasis, en fraction de δ_m^f : le joint naît au sommet d'une branche élastique fictive de longueur δ_0 = rise·δ_f, posée dans la direction de la traction tamponnée. Trois effets, tous nécessaires : raideur de charge/décharge **bornée** à k₀ = t_ins²/(2·rise·G_C) ; **continuité de traction** à l'insertion (à rise = 0, `split()` rendait 0 là où la facette liée transmettait t_ins — un Dirac par insertion) ; ∫t dδ sur la branche adoucissante **inchangé** = G_C (banc `tests_f2/check_jointtsl_header.cpp`). `rise = 0` explicite est accepté avec avertissement : **prouvé inconditionnellement instable** le 2026-09-11 (trois runs, voir CHANGELOG). Refusée hors `camacho`. |
| `jointMixLaw` | `none` \| `bk` | `none` | 2.4 eq. 18 | `bk` : G_C de Benzeggagh-Kenane, **gelé au ratio de mode à l'insertion**. |
| `jointBKEta` | réel | `2.0` | 2.4 eq. 18 | η, plage recommandée [1,5 ; 2,5]. |
| `jointFrictionMobilised` | `off` \| `damage` | `off` | 2.4 | `damage` : le terme frottant est multiplié par **D** (`jtsl::shearCap(..., true)`). |
| `jointEtaN` | Pa·s/m | `0` | 2.5 | option B, terme visqueux **normal** après insertion. |
| `jointEtaS` | Pa·s/m | `0` | 2.5 | option B, terme visqueux **tangentiel** (n'existe pas dans g0). |
| `jointViscousInCriterion` | `on` \| `off` | `on` | 2.5 | `off` : `tauLim` évalué sur `sigEl` et non sur `sig`, pour que l'amortisseur ne crée pas d'effet de vitesse sur le seuil. |
| `gbCombine` | `mean` \| `min` | `mean` | 2.6 eq. 21 | `min` : t_n0^F = χ_t·min(t_n0^p, t_n0^q) sur les facettes inter-granulaires. |
| `jointWeibullXu` | réel ≥ 0 | `0` | 2.6 eq. 22 | seuil x_u du Weibull à 3 paramètres, en **fraction** du seuil de facette. |
| `jointWeibullScale` | `volume` \| `area` | `volume` | 2.6 | `area` : facteur d'échelle x₀(A_F/A_ref)^(−1/m_w) sur l'**aire** de la facette. |
| `compBandLength` | `solver` \| `tetEdge` | `solver` | 1.3 eq. 6 | `tetEdge` : renseigne `s.lcComp = jtsl::tetEdgeLength(V_e)`. |
| `energyBreakdown` | `off` \| `on` | `off` | 3.2 eq. 26 | `on` : colonnes `eVp,eDamT,eDamC` ajoutées **en fin** de l'en-tête de `history.csv`. |

### Exclusions mutuelles à vérifier à la lecture (throw, message explicite)
1. `jointTSL = camacho` **exige** `insertion = adaptive` (une loi sans K_n n'a aucun
   sens sur un joint intrinsèque qui doit coller le continuum).
2. `strainRateDIF` armé **et** (`jointEtaN > 0` ou `jointEtaS > 0`) → refus.
   La note : « **une seule des deux options, jamais les deux** » (§2.5).
3. `mhForm = principal` exige `meridian = power`.
4. `jointMixLaw = bk` exige `jointTSL = camacho` (BK n'a de sens que sur une énergie
   de rupture mixte unique, pas sur deux longueurs critiques séparées).

---

## 3. Travaux par fichier

### LOT A — `src/MatLaw.cpp` + `include/rockim/MatLaw.hpp`
1. `MatState::lcComp` (+ doc). `b_c = fc0 * (s.lcComp > 0 ? s.lcComp : lc) / compGIIc`.
2. `bulkTensionDamage = off` : encadrer `MatLaw.cpp:483-507` par `if (!br_.tensionNone)`.
   **Attention** : `GfLoc` (l. 484) est réutilisé l. ~610 par `erodeWfrac` — le sortir du
   bloc ou désarmer aussi `erodeWfrac`. Lever les deux gardes de bande en traction.
   Refuser la clé sur les lois qui ignorent `BrickOpts` (`mc`, `saksala2011`, `dpdfh`,
   `dfhplus`, `cdp`) sur le modèle exact de la garde de `tensionDamage = fixed`.
3. `mhForm = principal` : décomposition spectrale du prédicteur, retour en contraintes
   principales sur `f = σ̄₁ − σ̄₃ − B⟨σ̄₃⟩ⁿ − σ_c`, traitement des arêtes (σ₁=σ₂, σ₂=σ₃) et
   de l'apex, recomposition. Le squelette existe : `MohrCoulombLaw`
   (`MatLaw.cpp:83-198`, `planeReturn`/`edgeReturn`/`apexClamp`) — l'envelope y est
   linéaire, il faut la remplacer par la loi puissance (Newton borné local, `yieldQPower`
   donne le patron). Conserver la viscosité : `dlam = F/(H + η/dt)`.
   **Le crochet de Macaulay ⟨σ̄₃⟩ doit être respecté** : pour σ̄₃ ≤ 0 la note donne
   `f = σ̄₁ − σ̄₃ − σ_c`, pas un repli sur le cône linéaire.
4. `dpFlowForm = mc` : direction d'écoulement `(1, 0, −m_ψ)` en base principale.
5. Avertissement (pas un refus) si `meridian = power` et `capP0` non posé explicitement :
   `law = saksala` allume un cap **écrouissant** par défaut, contraire à « viscoplasticité
   parfaite » (§1.2), et invisible en triaxial.

### LOT B — `src/Fdem3dSolver.cpp` + `include/rockim/Fdem3dSolver.hpp`  *(solveur cible)*
1. **Instrumentation (étape bloquante, à faire en premier).**
   Récolter `e.st.wPlas / wDamT / wDamC` dans `Fdem3dSolver::finalize` sur le modèle de
   `src/Fem3dSolver.cpp:2597`. Colonnes `eVp,eDamT,eDamC` dans `history.csv`
   (`Fdem3dSolver.cpp:4685`) sous `energyBreakdown = on`. Champs VTU `damage`, `omegaC`,
   `epvEq` **quand `law_` existe** (pas de clé — ajouter des tableaux nommés au VTU ne
   peut rien casser). `Fdem3dSolver.cpp:4618-4629`.
2. `compBandLength = tetEdge` → `e.st.lcComp = jtsl::tetEdgeLength(V0)`.
3. `facetAverage = volume` — 1 site pour la contrainte (`~2003`), 5 pour le taux
   (`~2015`, `~2066`). Les volumes sont déjà là (`Elem::V0`).
4. `facetRate = tensor` — stocker le tenseur taux par élément dans le repère **global**.
5. `insertionCriterion = elliptic` + `insertionHoldSteps` — lambda `testJoint` (`~2014-2024`).
   Compteur `nPhi` dans `struct Joint`, remis à 0 dès que Φ < 1. Sortir le compteur à côté
   de `nInserted_`.
6. `difExpT` / `difExpS` séparés, alimentés par `ε̇_eq` et `γ̇` quand `facetRate = tensor`.
7. `jointTSL = camacho` — champs `Stamp` dans `struct Joint` (`tmIns`, `beta`, `Gc`, `dmF`,
   `dmMax`), tamponnés dans `activateJoint()` (`~2051-2078`) à partir des `sig`/`tau`
   **déjà passés en argument**, branche nouvelle dans `jointForces()` (`~3000-3150`).
   La compression reste reprise par la pénalité de contact existante.
8. `jointMixLaw = bk` + `jointBKEta` — via `jtsl::stampInsertion`.
9. `jointFrictionMobilised = damage` — `jtsl::shearCap(..., true)` (`~3079`).
   **Mettre en cohérence le clamp d'activation** `fsNow` (`~2062-2066`) avec le même cap,
   sinon un joint inséré en compression subit une chute instantanée de μ|σ_n| au pas suivant.
10. `jointEtaN` / `jointEtaS` / `jointViscousInCriterion` (`tauLim` sur `sigEl`, `~3095`).
11. `gbCombine = min` (`~1641-1647`), `jointWeibullXu`, `jointWeibullScale = area`
    (`applyJointStatistics`, `~1735-1775`).
12. **Garde** : `erodeD` / `erodeEpv` / `erodeDc` posés en fdem3d → `throw`, avec le
    message « l'érosion n'est pas implémentée en FDEM (`grep -c eroded` = 0) : l'élément
    garderait masse, joints et nœuds, et continuerait d'entrer avec un poids 0,5 dans la
    moyenne de facette ».

### LOT C — `src/FdemSolver.cpp` + `include/rockim/FdemSolver.hpp`
**Miroir strict du LOT B**, mêmes noms de clés, même sémantique, mêmes messages.
Attention : le solveur 2D **duplique** la boucle de facette en une branche OpenMP et une
branche série (`~3446-3472` et `~3490-3510`) — toute modification doit être faite **deux
fois**, et les deux branches doivent rester bit-identiques entre elles.
`loading = grips | platens` reste 2D : c'est ici que les joints se calibrent.

---

## 4. Fichiers interdits aux lots

`include/rockim/KeysByMode.hpp`, `tools/keys_by_mode.json`, `include/rockim/JointTsl.hpp`,
`include/rockim/KeyGuard.hpp`, `include/rockim/Guards.hpp`, `include/rockim/Material.hpp`,
`include/rockim/YanSoftening.hpp`, `include/rockim/YangDif.hpp`, `src/main.cpp`,
`CMakeLists.txt`, `build_g1.cmd`, et **tout fichier d'un autre lot**.

Un besoin sur l'un de ces fichiers se **déclare dans le rapport de fin**, il ne s'applique pas.

---

## 5. Les trois gardes autorisées à refuser

Elles ne refusent que des clés **aujourd'hui inertes**, donc aucun deck correct ne les porte :
1. `erodeD` / `erodeEpv` / `erodeDc` en `fdem` et `fdem3d` (§ LOT B.12 / LOT C).
2. `jointTSL = camacho` sans `insertion = adaptive`.
3. `strainRateDIF` **et** `jointEta*` simultanés.

---

## 6. Recette

1. **Bit-identité** : `tools/bitid.py` inchangé sur les decks de référence.
2. **Falsifiabilité** : chaque clé nouvelle doit avoir un banc court où elle **change**
   le résultat dans le sens prévu, et une variante qui **doit échouer**.
3. **Bilan d'énergie** : sur le deck de la note, `W_bit = E_el + E_kin + D_vp + D_ωc +
   D_coh + D_fric + E_abs + E_art` doit fermer, et `D_coh / Σ A_e ≈ G_Ic` en mode I dominant.
4. **Objectivité** : `∫ t_m dδ_m = G_C` indépendamment du maillage et du dépassement
   à l'insertion — c'est la propriété que la loi initialement rigide apporte.
5. **Pas de temps** : mesurer le gain attendu (×1,96 selon l'audit) sous `jointTSL = camacho`.
   **Résultat (2026-09-11) : ce gain n'existe pas, et l'audit s'était trompé de comparaison.** Le
   ×1,96 comparait `insertion = adaptive` à `insertion = none` — un schéma **sans aucun joint**. Une
   loi extrinsèque stable a besoin (a) d'une raideur de charge bornée k₀ = t_ins²/(2·rise·G_C) qui,
   à rise = 1e-3, est du même ordre que la pénalité `10 E/h` de la note, et (b) de toute façon de la
   pénalité de contact en compression. Le budget CFL doit porter **toutes** les facettes, liées
   comprises (elles peuvent s'insérer à tout pas) — le lot B les en avait exclues, d'où un dt ×3,49
   « saturé »… et un run qui explose à 91 µs quel que soit dt. Règle : sous `camacho`, chaque
   facette pèse `max(k₀, kPara·pj)·A0/3` au budget.
