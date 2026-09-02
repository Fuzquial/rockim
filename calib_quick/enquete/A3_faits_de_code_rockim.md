# Tâche A3 — Faits de code `rockim_f2` pour les trois calibrations Red Bohus

Sources lues intégralement ou par sections : `src/FdemSolver.cpp` (8548 l.), `include/rockim/FdemSolver.hpp`, `include/rockim/Material.hpp`, `src/MatLaw.cpp` + `MatLaw.hpp`, `src/Config.cpp`, `src/main.cpp`, `include/rockim/YanSoftening.hpp`, `include/rockim/YangDif.hpp`, `src/Fdem3dSolver.cpp` (l. 600-630), `DOCUMENTATION_rockim.md` (§3-§6, §8, §5.16), `calib_quick/README.md`, `CHANTIER_f2.md` §7, `HANDOFF_2026-09-02.md`, `README.md` (racine, EN, l. 478-495), et le journal `out_q1_homog_P050/run.log` + `history.csv`. Chemins absolus sous `C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_f2/`. Aucune simulation n'a été lancée ; les seuls chiffres « mesurés » viennent des journaux existants.

Rappel structurel qui conditionne tout : le lecteur de configuration **ignore en silence toute clé inconnue** (`src/Config.cpp:15-30`, doc §4 l. 99) ; les nombres sont parsés strictement (`Config.cpp:38-56`).

---

## (i) Clés de joint : sémantique, défaut, formule

Tous les défauts ci-dessous sont lus dans `Material::from` (`include/rockim/Material.hpp:77-88`) ou dans `FdemSolver::init` ; les décks `calib_quick/*.cfg` tournent en **`insertion = adaptive` + `jointSoftening = yan`**, ce qui change la clé de pénalité effective (voir `jointPenaltyFactor`).

| clé | défaut | sémantique / formule | fichier:ligne |
|---|---|---|---|
| `ft` | 10e6 Pa | résistance en traction du **joint** (et seuil d'insertion adaptative σ_n ≥ ft·DIF sur la contrainte moyenne des deux éléments adjacents projetée sur la normale de l'arête). Élastique : dnE = ft/pj. Validation ft > 0 | `Material.hpp:28,80` ; insertion `FdemSolver.cpp:3312-3339` ; `setJointLengths` 3227-3228 |
| `cohesion` | 25e6 Pa | cohésion du joint : cap de cisaillement τ_lim = c_eff + μ_eff·T(σ_n) avec T = max(0, −σ_n) (`jointShearEnvelope = yan`, défaut) ou −min(σ_n, ft) (`yang`) ; c_eff = f(D)·c (yan) ou (1−D)·c (linear). Seuil d'insertion \|τ\| ≥ c·DIF + tanφ·T(σ_n) | `Material.hpp:29,81` ; `FdemSolver.cpp:4997-5020` ; `YangDif.hpp:134-136` ; insertion 3335-3339 |
| `frictionDeg` | 40° | angle de frottement de **pic** du joint : J.tanPhi = tan φ ; validé dans [0, 89) | `Material.hpp:30,82,175` ; `FdemSolver.cpp:2366` |
| `Gf` | 70 J/m² | énergie de fissuration mode I ; ouverture finale dnF = dnE + kI·Gf/ft avec **kI = 1/∫f dD = 2,5886** (yan), 2 (linear), 3 (`jointDeltaC = guo`) ; sous Weibull, Gf n'est PAS tiré par défaut | `Material.hpp:33,85` ; `FdemSolver.cpp:3237-3253` ; `YanSoftening.hpp:73-81` |
| `gfShearFactor` | 10 | GfII = facteur × Gf ; glissement final slipF = kI·GfII/c (même kI). Validé > 0 | `Material.hpp:34,86` ; `FdemSolver.cpp:2313,3252` |
| `jointPenaltyFactor` | 20 | pénalité **intrinsèque seulement** : pj = facteur·E_joint/h, h = ½(h_A + h_B) (taille inscrite 4A/P des deux éléments) dès que `voronoi_` est vrai — ce qui est le cas en `mesh = file` (l. 2025) et Voronoï ; hmin global seulement en `mesh = grid`. **En `insertion = adaptive` cette clé est INERTE** : pj = `insertionPenaltyFactor` (défaut 4)·E/h. Le journal du cas 1 le confirme : « activation penalty 4 E/h » (`out_q1_homog_P050/run.log:4`). Le solveur 3D avertit (`Fdem3dSolver.cpp:623-626`), le 2D **non** | `FdemSolver.cpp:2306-2307, 2359-2360, 2025` |
| `insertionPenaltyFactor` | 4 | pénalité des joints ACTIVÉS en adaptatif ; à l'activation, dn0 = min(σ_n, ft)/pj et glissement initial −τ0/pj assurent la continuité de contrainte | `FdemSolver.cpp:2306, 3408-3413` |
| `jointXi` | 0.05 | amortisseur de joint : c_d = 2ξ√(pj·L_trib·m_eff), m_eff = ½·min(m_a, m_b), borné à m_eff/dt ; bilatéral sur joint intact (D < 1), somme écrêtée ≤ 0 en contact (dn < 0). Règle maison : 0,01 en quasi-statique (doc l. 172) | `FdemSolver.cpp:903, 4967-4987` |
| `jointSoftening` | `linear` | `linear` : c_eff = (1−D)·c, traction linéaire de ft à 0 sur dnF − dnE = 2Gf/ft ; `yan` ≡ `munjiza` : facteur f(D) = [1 − (a+b−1)/(a+b)·exp(D(a+cb)/((a+b)(1−a−b)))]·[a(1−D) + b(1−D)^c], a = 0,63, b = 1,8, c = 6 (`yanA/B/C`), ∫₀¹f dD = 0,386307 ; D mixte = √(rn² + rs²) avec rn = (dn−dnE)/(dnF−dnE), rs = \|s\|/slipF | `FdemSolver.cpp:770-791` ; `YanSoftening.hpp:55-69` ; `FdemSolver.hpp:214-216` |
| `jointFrictionScaled` | 0 | 1 = le terme de Coulomb est aussi multiplié par f(D) : μ_eff = f(D)·tanφ (éq. 10 littérale de Yan) ; **lu seulement sous `jointSoftening = yan`** (inerte sous `linear`, piège noté dans le code l. 815-820) | `FdemSolver.cpp:784, 4999-5001` |
| `jointResidualMu` | −1 (= non posé) | frottement RÉSIDUEL : μ_eff = μ_res + (tanφ − μ_res)·g, g = f(D) (yan) ou max(0, 1−D) (linear). Généralise la précédente : μ_res = tanφ ⇒ défaut historique ; μ_res = 0 ⇒ `jointFrictionScaled = 1`. **Exclusives** : le run s'arrête si les deux sont posées (avec μ_res ≥ 0). Fonctionne aussi sous `linear`, contrairement à `jointFrictionScaled` | `FdemSolver.cpp:835-853, 5013-5016` ; `FdemSolver.hpp:437-450` |
| `jointWeibullM` | 0 (off) | m > 1 : chaque joint reçoit J.stat = (−ln(1−u))^(1/m)/Γ(1+1/m) (moyenne 1) ; **ft et cohésion × stat** ; Gf/GfII × stat seulement si `weibullScope = strengthGf` ; dnF/slipF recalculés. Champ VTU `ftScale` = stat. Se compose avec `jointSizeEffect` | `FdemSolver.cpp:2406-2473` (facteur 2424-2427, application 2450-2457) |
| `strengthCorrLength` | 0 | 0 = tirages i.i.d. par joint dans l'ordre des joints, `mt19937(fieldSeed)` ; > 0 = un champ gaussien `RandomField(W, H, ℓ, ℓ_B, angle, fieldSeed)` évalué au milieu du joint, copule u = Φ(g) (+ `strengthCorrLengthB`, `strengthCorrAngleDeg`) | `FdemSolver.cpp:2420-2446` ; `RandomField.hpp` |
| `fieldSeed` | seed + 777 | graine du champ/tirage, indépendante du maillage | `FdemSolver.cpp:2421-2422` |
| `seed` | 12345 | graine du **maillage** (Voronoï, jitter, phases) via `mt19937(seed)` ; en `mesh = file` le maillage ne l'utilise pas — seules les graines dérivées comptent : `fieldSeed` = seed+777, `weakPlaneSeed` = seed+5151, `jointPrebrokenSeed` = seed+4242. Tirages dépendants de la plateforme (doc §8.9) | `FdemSolver.cpp:1375, 1445, 1676, 2044, 2422, 2724, 3011` |
| `weibullScope` | `strength` | `strengthGf` = Gf et GfII suivent le facteur (Gf ∝ stat ⇒ ℓ_cz ∝ 1/stat) ; **aucune option ne garde ℓ_cz constant** (il faudrait Gf ∝ stat²) | `FdemSolver.cpp:2407-2408, 2453-2456` ; `MatLaw.cpp:1211-1214` |

Compléments utiles à la calibration : `jointDeath` (`separation` par défaut : le joint ne meurt qu'à dnMax > 3·dnF, sinon il reste le contact frottant de ses lèvres — `FdemSolver.cpp:371-380`) ; `jointShearUnload` (`plastic` par défaut, l. 859-865) ; `jointSizeEffect`/`jointZeff`/`jointSizeEffectM` (effet d'échelle (Zeff/V_J)^(1/m), doc l. 191-193, code 2852+).

---

## (ii) Phases : clés acceptées et joints de grain

### `phase.<nom>.*` — la liste exhaustive (`Material.hpp:206-224`)

Chaque phase part du bloc matériau global puis surcharge : **`rho`, `E`, `nu`, `ft`, `cohesion`, `frictionDeg`, `Gf`, `gfShearFactor`** (l. 214-221), **`fraction`** (obligatoire, `reqd`, l. 223 ; normalisée à somme 1, l. 234) et, côté tessellation, **`grainSize`** (`FdemSolver.cpp:2067-2074`, affinité de taille, exige ≥ 2 phases ; doc §5.16). Rien d'autre n'est lu sous ce préfixe (pas de `phase.<nom>.crushCap`, ni `phase.<nom>.mc*`, ni `phase.<nom>.jointWeibullM`). Validation stricte par phase (`Material.hpp:169-177`). `contactMu.<phase>` existe pour le contact général (`FdemSolver.cpp:1038-1055`).

Avec un bulk élastique (pas de `law`), les clés `phase.<nom>.ft/cohesion/frictionDeg/Gf/gfShearFactor` ne servent qu'aux **joints INTRA-granulaires** de la phase (type 0 : « intra-grain: bulk », `FdemSolver.cpp:2311-2315`) ; `E`, `nu`, `rho` servent au bulk (tables par phase `DmP_`, `nuP_`, `rhoP_`, l. 679-688, lues dans `elementForces` l. 4286) et à la pénalité des joints.

### Joints INTER-granulaires : `gbAlpha*` et `gbHeteroFactor` (`FdemSolver.cpp:2316-2329`)

Condition : `voronoi_` et `A.grain ≠ B.grain`. Avec s = `gbHeteroFactor` si `A.phase ≠ B.phase` (type 2), s = 1 sinon (type 1) :

```
E      = gbAlphaE   · ½(E_A + E_B)
ft     = s · gbAlphaTen · ½(ft_A + ft_B)
coh    = s · gbAlphaCoh · ½(c_A + c_B)
Gf     = s · gbAlphaGf  · ½(Gf_A + Gf_B)
GfII   = s · gbAlphaGf  · ½(gfs_A·Gf_A + gfs_B·Gf_B)      # même alpha que Gf
phiDeg = gbAlphaFric · ½(φ_A + φ_B)                        # pas de s
```

Noms exacts et défauts (1,0, tous > 0 exigés) : `gbAlphaTen`, `gbAlphaCoh`, `gbAlphaGf`, `gbAlphaE`, `gbAlphaFric`, `gbHeteroFactor` (`Material.hpp:186-198`). **Écart doc/code** : la doc §5.2 (l. 147) dit que `gbHeteroFactor` s'applique « sur les résistances » ; le code l'applique aussi à **Gf et GfII** (l. 2326-2328), pas à E ni φ.

**Capacité déjà présente et décisive pour la calibration 3** : surcharges **par paire de phases** `gb.<a>.<b>.{ft, cohesion, Gf, gfShearFactor, frictionDeg, E}` (`Material.hpp:104-135, 236-262` ; application `FdemSolver.cpp:2333-2348`). Une valeur posée **remplace** le résultat alpha/hetero pour cette propriété ; l'ordre des deux noms est indifférent ; `gb.<a>.<a>.*` fonctionne pour les frontières homophases (boucle j ≥ i, l. 244) ; `gfShearFactor` exige `Gf` posé (l. 254-258). C'est exactement la table à six paires d'Aboayanah et al. 2024 citée dans le commentaire. Pénalité de ces joints : pj = pf·E_joint/h avec E_joint la moyenne (ou la surcharge). Résumé imprimé au journal (l. 2371-2375).

Sortie : `type` (0/1/2) dans `fdem_joints_*.vtu` et `fdem_final_joints.csv` ; le résumé donne la fraction intergranulaire de la casse (`run.log:33`).

---

## (iii) `law = mc` et phases ; crushCap ; ce qu'exigerait une plasticité par phase

**Refus** : `src/FdemSolver.cpp:693-697` — `if (cfg_.has("law")) { if (phases_.n() > 1) throw "'law' (bulk constitutive law) is a SINGLE material model: it cannot be combined with mineral 'phases'..." }`. Doc §5.5 l. 638-639. Autre garde : `law` + litage TI refusés (l. 2511-2514). `main.cpp:91-100` n'autorise `law` qu'en fem3d/fdem/fdem3d.

**Pourquoi** : une seule instance `std::unique_ptr<MatLaw> law_` (`FdemSolver.hpp:400`), construite depuis le **bloc matériau global** `mat_` (surchargeable par `bulkFt/bulkCohesion/bulkFrictionDeg/bulkGf`, l. 711-724) via `MatLaw::make(kind, mBulk, cfg_, lcMax)` (l. 725) ; la loi `mc` lit ses clés globalement : `mcCohesion` (déf. = cohesion), `mcFrictionDeg` (déf. = frictionDeg), `mcDilationDeg` (déf. 0, ψ ≤ φ) dans `MatLaw.cpp:1220-1233`. `MohrCoulombLaw` (`MatLaw.cpp:54-173`) : critère à arêtes en contraintes principales, retour de Clausen (plan / arête / apex), écoulement non associé par ψ, **parfaitement plastique** (c, φ constants — ni adoucissement ni durcissement). Appel en déformation plane : `E3` avec ε_zz = 0, `law_->stress(E3, e.st, dt_, hEl_[eI])` (l. 4302-4312). L'état plastique est déjà **par élément** (`Elem::st`, `FdemSolver.hpp:143` ; `MatState::epsP`, `MatLaw.hpp:95`).

**crushCap** : `crushCapP_[phase] = cfg crushCap` sinon **8 × cohésion de la phase** (l. 679-688) → la clé est **globale** (une seule valeur pour toutes les phases) mais son **défaut est par phase**. Appliqué comme plafond radial du déviateur de von Mises (l. 4341-4361) ; **désactivé sous `law`** (cap = 1e300, l. 4354). Note in situ : le cap ne voit pas σ0 (l. 535-541). `meanTensionCapFactor` (déf. 3·ft de la phase, l. 419, 4363-4366) est un second garde-fou actif hors `law`.

**Ampleur d'une plasticité par phase (évaluation, non implémentée)** :

| pièce | ce qu'il faut | lignes touchées |
|---|---|---|
| conteneur | `std::vector<std::unique_ptr<MatLaw>> lawP_` indexé par `e.phase` + booléen `lawOn_` | `FdemSolver.hpp:400` |
| construction | boucle sur `phases_.mat[p]` avec `MatLaw::make(kind, phases_.mat[p], cfg, lcMax)` ; lever le garde 693-697 ; surcharges `bulk*` par phase | `FdemSolver.cpp:693-745` |
| clés | `phase.<nom>.mcCohesion/mcFrictionDeg/mcDilationDeg` : `MatLaw::make` lit `c.getd("mcCohesion")` global → passer un préfixe de clé ou un `Config` filtré | `MatLaw.cpp:1220-1233` |
| point chaud | `lawP_[e.phase]->stress(...)` ; 9 occurrences de `law_` en 2D (gardes l. 4302, 4337, 4354, 4363, 4374 + init) | `FdemSolver.cpp:4302-4312` |
| parité 3D | 8 occurrences de `law_` dans `Fdem3dSolver.cpp` (règle de parité du dépôt) | `Fdem3dSolver.cpp` |
| preuve | bit-identité clés absentes (deck GBM), `selftest-mc` intact, 1 repère de suite | `tools/verify_suite.py` |

Estimation : **~60-100 lignes en 2D** (+ autant pour la parité 3D), demi-journée avec les preuves — **ampleur moyenne**. Limites physiques à garder en tête : `mc` n'a ni adoucissement ni ψ(p) ; le bilan d'énergie « stocké élastique » est approché sous `law` (doc l. 1359-1361) ; `bulkDamage` et `bulkModel = neohookean` sont exclusifs de `law` (l. 353-358, 413-418).

---

## (iv) `history.csv` en scénario `tension`

Cadence : une ligne tous les ⌈nSteps/2000⌉ pas (`main.cpp:122`), vidage après chaque ligne (`historyFlush`, l. 140).

**`loading = grips` (défaut)** — en-tête `t,gripFy,sigma,sigmaPeak,nBroken` + `,nInserted,nDamaging` si `insertion = adaptive` (+ colonnes hydro/thermique/bulkDamage si armées) — `FdemSolver.cpp:7626-7641`, confirmé par `out_q1_homog_P050/history.csv:1`. Valeurs (l. 7705-7708, 4148-4169) :
- `gripFy` = Σ f_y sur les nœuds PRESCRIBED (rangée y = H, l. 2267-2273), en N par mètre d'épaisseur ; le bas est FIXED ;
- `sigma` = \|gripFy\|/(W·thk), **contrainte axiale TOTALE** (Pa) ;
- `sigmaPeak` = max glissant, verrouillé à la première chute sous 30 % du pic après nBroken > 0 (l. 4163-4169) ;
- `nBroken` = joints ayant atteint D ≥ 1 (compté une fois, `tBreak`, l. 5121-5122 ; les pré-cassés exclus) ;
- `nInserted` = joints activés (bonded = false, hors pré-cassés), `nDamaging` = joints avec 0 < D < 1 (l. 7687-7694) — **c'est le proxy disponible de σ_ci** (premier `nInserted`/`nDamaging` > 0 ; la doc §6.1 l. 1345-1347 dit à tort qu'aucun compte d'insertion n'est écrit en tension — le code l'écrit en adaptatif).
- **Aucune colonne de déformation** en mode grips : `plot_quick.py:38-60` reconstruit ε par l'intégrale analytique de la rampe cosinus des mors / H.

**`loading = platens`** — colonnes ajoutées : `epsPlaten,epsSpec,epsGauge,nBrokTen,nBrokShear,nFrag,confAchieved,peakLocked` (l. 7634-7636, 7709-7717). Définitions (`gaugeStrain`, l. 3871-3888 ; `setupStrainGauge` l. 3841-3869) :
- `epsPlaten` = (gap0 − gap)/gap0 (déformation MACHINE, contient la complaisance de la pénalité de platine) ;
- `epsSpec` = (ū_y bas − ū_y haut)/H (faces de l'éprouvette, positif en compression) ;
- `epsGauge` = (ū_y bande basse − ū_y bande haute)/L0, bandes de demi-largeur hmin/2 à `gaugeLoFrac`·H et `gaugeHiFrac`·H (0,25/0,75) — l'extensomètre intérieur ;
- `nBrokTen`/`nBrokShear` = nombre de joints avec `bmode` = 1 / 2, attribué **une fois** à l'instant D = 1 : rn = (dnMax − dnE)/(dnF − dnE), rs = max\|slip\|/slipF (ou le moteur origine), bmode = 1 si rn ≥ rs sinon 2 (l. 5126-5143) ; les pré-cassés portent bmode = 4 et ne comptent dans aucun des deux ;
- `confAchieved` = jauge de confinement figée (voir vi) ; `peakLocked` = marqueur 0/1.
`setupStrainGauge` sort immédiatement hors platens (l. 3842), donc `epsSpec`/`epsGauge` **n'existent pas** en grips.

**Déformation LATÉRALE ou VOLUMIQUE : NON, dans aucun mode.** Où la lire sans code :
1. `fdem_XXXX.vtu` : points = X0 + u (`writeFrame`, l. 7512), champ **`epsXX` = ε_xx co-rotée par élément** (`e.exx = eps(0)`, l. 4484 ; VTU l. 7548) — pour l'éprouvette chargée selon y, c'est directement la déformation latérale par élément ; moyenne pondérée par l'aire = ε_lat de l'éprouvette, à la cadence des `frames` (22 dans les decks) ; les faces x = 0 / x = W s'identifient par les coordonnées de la frame 0 (indices de points stables ; nœuds dupliqués par élément, 3 par triangle) ;
2. `fdem_final_elements.csv` : centroïdes DÉFORMÉS (X0 + u), fin de run seulement (l. 7772-7780) ;
3. `fdem_final_joints.csv` : extrémités déformées d'une lèvre (l. 7799-7807) ;
4. `fdem_nodal_displacement.csv` (x0, y0, ux, uy) **seulement si `gravity > 0`** (l. 7783-7790).
En déformation plane ε_v = ε_ax + ε_lat (ε_zz = 0), donc σ_cd = inversion de ε_v demande ε_lat à la cadence de l'historique → ajout de code (voir liste finale, item 1).

---

## (v) Coût

**D'où vient dt** (`computeStableDt`, l. 3991-4072) :

```
K_i = 2·E_phase·thk                                   (élément, l. 3993-3995)
    + Σ_joints  kPara·pj·(L0/2)·thk  sur les 4 nœuds  (l. 4000-4003 ; kPara = 2 si jointElastic = parabolic)
dt_crit,i = 2·sqrt(m_i / (K_i + extraContacts·kContact))     (l. 4028-4031 ; extraContacts = 2, kContact = max(kp_, kpPlaten_) = E_max·thk en tension, l. 917, 4014)
dt = dtFactor · min( min_i dt_crit,i , hCfl/c_P , borne visqueuse )   (l. 4036-4038, 4071-4072 ; dtFactor = 0,2)
```

avec m_i = ρ·A0·thk/3 par copie de nœud (l. 2264), pj = pf·E/h, h = ½(h_A + h_B) (4A/P), hmin = min des h (l. 2030-2031). **En adaptatif la boucle de joints est inconditionnelle : les joints encore liés comptent** (à pf = 4 au lieu de 20, d'où le « dt ×2 » de l'adaptatif, CHANTIER §7.1 l. 348-350 ; doc l. 179). Mesuré : cas 1 (3388 triangles, hmin 0,2795 mm) dt = 4,06 ns, 541 571 pas ; cas 3 (5814, hmin 0,199 mm) dt = 3,11 ns, 706 964 pas (`run.log:8-9` des deux sorties). Estimation d'ordre de grandeur pour un triangle de 0,8 mm (h ≈ 0,4 mm, E = 80 GPa) : élément 1,6e11, deux joints à 4E/h ≈ 6,4e11, contacts 1,6e11 N/m → les ressorts de joint font ~2/3 de K_i, et la borne CFL (≈ 0,28 mm/6050 m/s = 4,6e-8 s → 9e-9 après dtFactor) ne mord pas [estimation, NON MESURÉE].

**Passer `jointPenaltyFactor` de 20 à 10** : en `insertion = adaptive` (tous les decks `calib_quick`) **gain nul, la clé est inerte** (l. 2306-2307). En intrinsèque : dt ∝ 1/√K ≈ ×1,3-1,4 (les termes élément + contact ne bougent pas) [estimation] et la complaisance en série doublerait : doc l. 171 « ~4-5 % sur E au facteur 20 », mesurée −4,8 % (`t9`, doc l. 1041) et −2,0 %/−4,8 % (`t14b`, l. 1189) → ~8-10 % à 10 [extrapolation, NON VERIFIE]. En adaptatif, le levier équivalent est `insertionPenaltyFactor` (4 → 2 : dt ≈ ×1,2 d'après la décomposition ci-dessus, complaisance nulle avant insertion par construction, doc l. 179 ; effet sur la décharge/contact des joints activés non documenté [NON VERIFIE]).

**OpenMP** : 13 régions `#pragma omp` en 2D — `elementForces` (l. 4283), `bodyForces` (si gravité), `insertionSweep` (3295), `jointForces` (5186, par fils avec réduction en ordre de fil), recherche de contact **seulement si ≥ 4096 nœuds actifs** (5931, sinon sérielle), `toolContact`, `integrate` (7229, 7360). Aucun `omp_set_num_threads` dans le code (grep) : le nombre de fils est celui du runtime. Doc §3.2 l. 68 ne donne que la sémantique (1 fil = bit-identique). **Seule mesure d'efficacité existante** : `README.md:492-495` (racine), démo percussion Voronoï, 18 cœurs logiques hybrides : 267 s → 180 s, ×1,5. **Aucune mesure à ~3000 éléments** ; le cas 1 fait 146 s/541 571 pas = 0,27 ms/pas (`run.log:36`, nombre de fils non journalisé). La phrase du README de calib « sans la variable rockim se bride à ~3 cœurs » (`calib_quick/README.md:5`) n'a aucune contrepartie dans le code ni mesure jointe [NON VERIFIE]. Le tableau de coût du README (2930 éléments, 5,13 ns, 429 000 pas, 77 s, l. 63-64) **ne correspond pas** au journal de la sortie citée (3388 éléments, 4,06 ns, 541 571 pas, 146 s) [incohérence à trancher]. Protocole minimal pour mesurer : deck cas 1 avec `T = 2e-4`, `OMP_NUM_THREADS` ∈ {1, 2, 4, 7, 14}, ~10 min au total.

**Processus indépendants en parallèle : sûrs.** Toutes les sorties sont écrites sous le dossier de sortie passé en argument (`history.csv` `main.cpp:124`, VTU/`frames.csv` l. 7605, CSV finaux 7772-7826) ; aucun fichier temporaire, aucun état partagé, lectures seules (`cfg`, `meshFile`) ; états statiques limités à des lectures d'environnement (l. 57, 4697, 5727). Conditions : dossiers de sortie distincts, `OMP_NUM_THREADS` fixé par processus pour que la somme ≤ cœurs (deux runs à 14 fils se sur-souscrivent), et la règle de la mémoire de session (un job 14 cpus à la fois, OneDrive en pause) à faire valider par Fernando.

---

## (vi) Confinement et colonne `sigma`

| clé | défaut | sémantique | fichier:ligne |
|---|---|---|---|
| `confiningPressure` | 0 (off) | p > 0 [Pa] ; pression SUIVEUSE sur les faces extérieures **d'origine** sélectionnées à l'init ; force par nœud −½·p·L·thk·n avec n la normale sortante **courante** (Q − P) ; travail compté (`confWork_`) | `setupConfinement` 3907-3985 ; `confiningForces` 6411-6429 |
| `confiningRamp` | 0 | p(t) = p·½(1 − cos(πt/rampe)) ; 0 = échelon (WARNING) | 6413-6415, 3974-3979 |
| `confineFaces` | `sides` | `sides` : arêtes dont les deux nœuds sont à x < 1e-9 ou x > W − 1e-9 ; `all` : + toutes les faces sauf haut/bas en TENSION (mors), sauf la face outil hors tension, sauf le fond FIXED ; `bore` : faces à moins de `boreSelectR` de (`boreCX`, `boreCY`) | 3950-3961 |
| `confineGaugeTime` | 3·rampe (plancher 20 dt) | instant unique où `confAchieved_` = moyenne pondérée par l'aire de `e.sxx` sur les éléments dont le centroïde est dans le **cœur central** (\|cx − W/2\| ≤ W/4, \|cy − H/2\| ≤ H/4) ; figée ensuite ; colonne `confAchieved` en platens seulement, résumé dans tous les cas (−50 MPa à 5e-10 % dans le cas 1, `run.log:32`) | 4140-4145, 7181-7192, 8147-8158 |

**`thk`** = clé `thickness`, défaut **1,0 m** (l. 112 ; `FdemSolver.hpp:404`) : `sigma` = \|gripFy\|/(W·thk) est en Pa, `gripFy` en N par mètre d'épaisseur ; la masse porte aussi thk (l. 2264).

**Point à trancher pour le dépouillement (0 ligne de C++)** : `sigma` est la contrainte axiale **totale**. Pendant `pullDelay`, les mors sont bloqués (vg = 0, l. 7256-7259) et seuls les flancs sont pressés ; en déformation plane avec ε_y = 0 la réaction axiale vaut ν/(1−ν)·σ₃ = 16,67 MPa — c'est exactement l'offset mesuré dans `out_q1_homog_P050/history.csv` (moyenne 16,6667 MPa sur 2,5-3,0e-4 s). Le README définit q = sigma − cet offset (l. 150-151), soit 726,6 − 17,2 = 709,4 MPa ; le q expérimental des cibles est σ₁ − σ₃ = sigma − 50 = **676,6 MPa (+12,9 % vs 599,2, et non +18 %)** ; cas 3 : 514,7 + 17 − 50 ≈ 482 MPa (−20 %, et non −14 %). L'état de consolidation simulé n'est pas isotrope (σ_y = 16,7, σ_x = 50, σ_z = 16,7 MPa au départ de la rampe axiale) : une consolidation isotrope demanderait un ajout (liste ci-dessous, item 6), les mors excluant toute pression axiale (`confineFaces = all` saute haut/bas en TENSION, l. 3957).

---

## Écarts doc / README / code relevés

1. `gbHeteroFactor` multiplie Gf et GfII (code 2326-2328) — la doc §5.2 ne parle que des résistances.
2. Doc §6.1 l. 1345-1347 : « ni le nombre de joints insérés … n'est écrit » — faux en adaptatif (`nInserted,nDamaging`, l. 7632-7634).
3. CHANTIER §7.1 l. 345-346 : « en `mesh = file`, pj = 4E/hmin avec hmin GLOBAL » — le code met `voronoi_ = true` en `mesh = file` (l. 2025) et prend donc h **local** (l. 2359) ; ce qui reste vrai est que le plus petit triangle fixe dt par sa masse nodale. À vérifier sur le binaire de l'époque [NON VERIFIE].
4. `calib_quick/README.md` l. 63-64 vs `out_q1_homog_P050/run.log:8` : 2930/5,13 ns/77 s contre 3388/4,06 ns/146 s.
5. Le 2D n'avertit pas qu'un `jointPenaltyFactor` posé est inerte en adaptatif (le 3D le fait, `Fdem3dSolver.cpp:623-626`) — le commentaire « pénalité 20 » des decks est trompeur.
6. Doc l. 1077 (`weakPlaneGf = follow` : « l_cz inchangée ») contredit la règle ℓ_cz = E·Gf/ft² du README (Gf ∝ f donne ℓ_cz ∝ 1/f).

---

## À ajouter au code pour la calibration — classé par ampleur

**Petit (≤ 30 lignes, sortie seule ou garde, bit-neutre)**
1. Colonnes `epsLat` (ū_x face x = W − ū_x face x = 0)/W et `epsVol = epsAx + epsLat` dans `history.csv`, dans les deux modes (`setupStrainGauge` 3841-3869 à armer aussi en grips ; `gaugeStrain` 3871-3888 ; en-tête/ligne 7626-7641, 7705-7717) → σ_cd par inversion de ε_v. Variante sans code : `frames` ↑ et moyenne d'`epsXX` des VTU.
2. `epsSpec`/`epsGauge` en mode grips (même chantier que 1) : ε mesurée au lieu de l'intégrale analytique de `plot_quick.py`.
3. Colonnes `Dmax`/`Dmean` (endommagement des joints) — lacune notée par la doc l. 1346 ; ~10 lignes dans `historyRow`.
4. `weibullScope = lcz` : `J.Gf *= stat²`, `J.GfII *= stat²` (l. 2453-2456 + validation `MatLaw.cpp:1211-1214`) pour disperser (ft, c) à ℓ_cz constant — la règle du README l. 101-104.
5. Avertissement 2D « `jointPenaltyFactor` posé mais schéma adaptatif » (miroir de `Fdem3dSolver.cpp:604-626`, à placer vers l. 898-901).
6. `plot_quick.py` : q = sigma − `confiningPressure` (et documenter l'écart de 33 MPa) — 0 ligne C++.

**Moyen (50-150 lignes, touche la physique, preuve de bit-identité requise)**
7. **Consolidation isotrope** : pression axiale égale à p sur la rangée des mors (ou platines) pendant `pullDelay`, puis bascule en vitesse imposée — `integrate` 7256-7275 / 7385-7397 et sélection des faces 3955-3959 ; clé opt-in (`confineAxial = true`). Fait partir la courbe à q = 0 comme l'essai.
8. **Plasticité par phase** (`law` par phase, voir iii) : ~60-100 lignes 2D + parité 3D + `phase.<nom>.mc*`.
9. Clé `phase.<nom>.crushCap` si l'on garde le bulk élastique plafonné : 5 lignes (l. 686) — petit en code, moyen en interprétation (le cap est un garde-fou, pas une loi).

**Grand (> 150 lignes ou changement de structure)**
10. Weibull de VOLUME en 2D (`matWeibullM` n'existe qu'en fem3d et n'est consommé que par saksala2011/dpdfh via `MatState::ftScale`, doc l. 649-652 ; `MatLaw.hpp:106`) : exigerait des tables élastiques par élément (`DmP_` est par phase, l. 679-688, 4286) — pas nécessaire tant que l'hétérogénéité est portée par les joints.
11. Adoucissement/dilatance variable dans `mc` (loi parfaitement plastique aujourd'hui, `MatLaw.cpp:54-173`).

**Sans code, mais à budgéter** : réserve 1 du README (ε̇ ≈ 5 s⁻¹, pullV ÷ 10 → coût ×10) ; polydispersité : hmin suit le plus petit grain (dt ÷ 3,6 à σ = 0,8, doc §5.16 l. 1311-1315) — poser `phase.biotite.grainSize` pour éviter le grain unique.

---

## [NON VERIFIE] — récapitulatif
- Nombre de fils réellement utilisés par les runs `calib_quick` (non journalisé) et la bride « ~3 cœurs » sans `OMP_NUM_THREADS`.
- Efficacité OpenMP à 3000 éléments (aucune mesure ; seule mesure : ×1,5 sur une percussion Voronoï 18 cœurs).
- Les gains chiffrés ×1,3-1,4 (pénalité 20 → 10, intrinsèque) et ×1,2 (4 → 2, adaptatif) sont des estimations à partir de la formule de dt, non mesurées ; la complaisance ~8-10 % au facteur 10 est une extrapolation linéaire de la doc.
- L'affirmation de CHANTIER §7.1 sur hmin global en `mesh = file` (contredite par le code actuel).
- L'origine du tableau de coût du README (2930 éléments) — aucune sortie ne lui correspond.