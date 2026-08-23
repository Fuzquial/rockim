# Volet 1 — L'insertion adaptative de rockim face à Yan, Zheng & Wang (IJRMMS 169, 2023, 105439)

Budget web épuisé (8/8). ScienceDirect (403) et ADS (405) ont refusé la récupération : le contenu de l'article est donc reconstitué à partir des **résumés de moteur de recherche portant sur la page éditeur** (marqués `[VERIFIE-résumé]`) et de la mémoire projet (`[MEMOIRE]`). Tout ce qui concerne rockim est marqué `[CODE]` et vérifié ligne à ligne.

---

## A. La règle réellement implémentée dans rockim `[CODE]`

Fichier : `C:\Users\fuzquianoalricabi\simulations\FDEM\rockim\rockim_p1\src\FdemSolver.cpp`, lignes 1855-2067 (`buildBindingTables`, `rebindVertex`, `insertionSweep`, `activateJoint`), plus `include/rockim/YangDif.hpp:119-121`.

### A.1 Quand — cadence et périmètre

`insertionSweep()` est appelée **une fois par pas de temps**, exactement une fois : les deux occurrences (l. 2692 et 2703) sont les deux branches du `if (fProf.on)` de `step()`. Position dans le pas : **après `elementForces()`/`bodyForces()`, avant `jointForces()`** — donc un joint né au pas *n* transmet déjà sa traction au pas *n* (l. 2703, commentaire explicite).

Le balayage est **exhaustif et non filtré** : boucle sur `jt_` entier, seul filtre `if (!J.bonded) continue`. Pas d'intervalle de vérification (contrairement à Pandolfi & Ortiz qui vérifient tous les *N* pas), pas de sous-liste de candidats, pas de restriction au voisinage des fissures existantes. Parallélisation OpenMP `schedule(static)` avec vecteurs `Hit` par thread fusionnés en `critical`.

### A.2 Où — évaluation de la contrainte

```cpp
double sxx = 0.5 * (A.sxx + B.sxx);   // A = el_[J.eA], B = el_[J.eB]
double sig = n·σ·n;   double tau = e·σ·n;
```

- **Moyenne arithmétique non pondérée des deux tenseurs élémentaires voisins** (CST → un tenseur constant par élément). Ni partition nodale des forces (Camacho & Ortiz), ni limite stockée par facette (Pandolfi & Ortiz), ni pondération par aire.
- Le repère de l'arête est **actualisé en configuration déformée** : `P` et `Q` sont les milieux des copies co-localisées (`X0_+u_`) aux deux extrémités ; `e = (Q−P)/L`, `n = (e_y, −e_x)`.
- Garde numérique unique : `if (L < 1e-14) continue` (arête dégénérée).
- Convention : `sig > 0` en traction.

### A.3 Quel seuil — l'enveloppe

```cpp
fs = dC*J.coh + J.tanPhi * mcFrictionTerm(sig, J.ft, yangEnv_);   if (fs<0) fs=0;
if (sig >= dT*J.ft || std::abs(tau) >= fs) → hit
```

avec `mcFrictionTerm(sig, ft, yang) = yang ? -min(sig, ft) : max(0, -sig)`.

Deux enveloppes distinctes selon la clé `shearEnvelope` :

| `shearEnvelope` | terme de frottement | en compression | en traction (0 < σn < ft) |
|---|---|---|---|
| **défaut** (`yangEnv_ = false`) | `max(0, −σn)` | `fs = c + tanφ·|σn|` | `fs = c` (aucune réduction) |
| `= yang` | `−min(σn, ft)` | idem | `fs = c − tanφ·σn`, saturé à `c − tanφ·ft` |

**Point de conformité à trancher** : la formule de l'article telle que tu la cites (`|τ| ≥ c − σn tanφ`, sans réserve) correspond à `shearEnvelope = yang`. **Le défaut de rockim tronque la branche en traction.** Vérifie la clé dans le .cfg du tunnel : c'est un écart de fidélité silencieux (avec c = 25 MPa, φ = 40°, ft = 10 MPa, l'écart sur `fs` atteint 8,4 MPa près de la coupure en traction).

Le critère est un **OU** de deux branches indépendantes (traction OU cisaillement), sans mode mixte, sans norme effective. À comparer : la littérature CZM extrinsèque récente utilise une contrainte effective `σ_eff = √(σn² + σt²/β²)` avec compression annulée `[VERIFIE]` (arXiv 2511.14323) — un critère **plus sélectif** que le OU, qui lui déclenche sur la plus permissive des deux branches.

**Facteurs dynamiques (DIF de Yang)** : si `difOn_`, le critère voit la résistance **dynamique** (`dT·ft`, `dC·c`, `ė` = moyenne des deux éléments) ; le terme de frottement n'est **pas** amplifié (choix documenté l. 1958-1962).

**Aucune hétérogénéité de seuil par défaut** : `J.stat` (facteur de Weibull replié dans `ft` et `coh`) vaut **1 pour tous les joints** tant que `jointWeibullM ≤ 0`, valeur par défaut (`FdemSolver.cpp:1745`). Sur un maillage régulier et un champ lisse, **toutes les arêtes de la zone plastique ont rigoureusement le même seuil**.

### A.4 Ordonnancement et débit

```cpp
std::sort(hits.begin(), hits.end(), [](x,y){ return x.jI < y.jI; });
for (const Hit& h : hits) activateJoint(h.jI, h.sig, h.tau);
```

- Tri **par indice de joint**, motivé explicitement par le déterminisme inter-threads — **pas par sévérité**.
- **Toutes** les arêtes candidates du pas sont activées, **sans plafond, sans quota, sans relaxation intermédiaire** : la contrainte est gelée à l'état évalué en début de balayage et n'est pas recalculée entre deux activations. Aucune clé de bridage n'existe (`grep` : pas de `maxInsert*`, `insertEvery`, `insertionInterval`).
- Conséquence directe : si 500 arêtes franchissent l'enveloppe au même pas, les 500 sont coupées — alors que physiquement la première à s'ouvrir déchargerait ses voisines. **C'est le contraire de Camacho & Ortiz, qui avancent arête par arête** `[MEMOIRE, acquis]`.

### A.5 Gardes à l'insertion — `activateJoint`

1. **Idempotence** : `if (!J.bonded) return;`.
2. **DIF figé** : `ft *= dT`, `Gf *= dT`, `coh *= dC`, `GfII *= dC`, puis recalcul de `dnE = ft/pj`, `dnF = dnE + kI·Gf/ft`, `slipF = kI·GfII/coh` (`kI = 2` en linéaire, `1/yanI_` en loi Yan/Munjiza). Multiplicatif, se compose avec Weibull.
3. **Continuité de contrainte normale** : `J.dn0 = min(sig, ft)/pj`. Le joint transmet à ouverture géométrique nulle **exactement** la traction que portait le continuum. Arête déclenchée en traction (`sig ≥ ft`) → `dn0 = ft/pj = dnE` → le joint démarre **pile au pic**, D = 0 mais adoucissement immédiat. Arête déclenchée en cisaillement → `dn0 = sig/pj` sous-critique, **négatif si σn < 0** (joint pré-comprimé).
4. **Continuité de contrainte tangentielle** : `tau0 = clamp(tau, ±fsNow)` avec `fsNow` calculé sur les résistances *post-DIF*, puis `slip[0] = slip[1] = −tau0/pj` — décalage de glissement plastique tel que la traction d'essai `pj·(dtg − slip)` vaille `tau0` à `dtg = 0`.
5. **Topologie** : `rebindVertex()` sur les deux sommets. Union-find sur l'éventail d'éléments, arêtes encore liées = arêtes connectantes ; les groupes ne peuvent que **se scinder**. Reproduit la fig. 7 de l'article (pointe de fissure = sommet entier, sommet traversé = éclatement en composantes).
6. `++nInserted_` — simple compteur, **pas un plafond**.

### A.6 Pénalité

`pj = insertionPenaltyFactor · E/h`, **défaut 4,0** (`FdemSolver.cpp:1668`, `:505`) contre 20 E/h pour la base intrinsèque de rockim `[MEMOIRE]`. Ordre de grandeur : la littérature CZM extrinsèque place `k⁻ ≈ α E/h` avec `α ∼ 10¹–10⁴`, ligne de base 10 E/h `[VERIFIE]` (arXiv 2511.14323) ; Yan compare à 100 E `[MEMOIRE]`. **rockim est donc en bas de la fourchette d'un facteur 2,5 à 25.**

Avec les chiffres du .cfg tunnel (E = 50 GPa, h = 2 mm, ft = 10 MPa, Gf = 70 J/m²) : `pj = 1,0e14 Pa/m`, soit **4× la raideur d'un élément** (`E/h = 2,5e13`) ; `dnE = 0,1 µm` contre `dnF − dnE = 2Gf/ft = 14 µm`. La pénalité ne pilote donc **pas** l'énergie de mode I (la branche adoucissante est 140× plus longue) : elle pilote la **raideur de contact/compression** et le glissement élastique. Chaque arête coupée ajoute ~25 % de la complaisance d'un élément dans la direction normale.

### A.7 Récapitulatif : ce qui n'existe pas

| Dispositif | rockim | Référence |
|---|---|---|
| Limite de traction tirée **par facette** | non (sauf `jointWeibullM` opt-in) | Pandolfi & Ortiz 2002 `[MEMOIRE]` |
| Weibull par facette comme frein physique | opt-in, **désactivé par défaut** | Zhou & Molinari 2004 `[MEMOIRE]` |
| Avancement arête par arête, relaxation entre insertions | non | Camacho & Ortiz 1996 `[MEMOIRE]` |
| Vérification tous les N pas (throttling temporel) | non (tous les pas) | Pandolfi & Ortiz `[MEMOIRE]` |
| Ordonnancement par sévérité / dépassement | non (par indice) | — |
| Plafond d'insertions par pas ou par sommet | non | — |
| Distinction nucléation / propagation | non | — |
| Contrainte effective mixte plutôt qu'un OU | non | arXiv 2511.14323 `[VERIFIE]` |

---

## B. Ce que dit l'article, et verdict de conformité

### B.1 Identification `[VERIFIE]`

C. Yan, Y. Zheng, G. Wang, *A 2D adaptive finite-discrete element method for simulating fracture and fragmentation in geomaterials*, **Int. J. Rock Mech. Min. Sci. 169 (2023) 105439**, DOI `10.1016/j.ijrmms.2023.105439`. Confirmé par ADS (`2023IJRMM.16905439Y`) et par la page ScienceDirect `S1365160923001132`.

### B.2 Contenu annoncé `[VERIFIE-résumé]`

Résumé de la page éditeur : FDEM 2D à **éléments cohésifs insérés dynamiquement**, insertion adaptative « quand la contrainte dépasse une valeur critique » ; la méthode est **identique au FDEM conventionnel à deux différences près** : (i) un **modèle constitutif purement post-pic** pour les cohésifs, (ii) un **nouveau schéma de mise à jour nodale**. Deux points de conformité importants pour rockim :

- « post-pic » = pas de branche élastique dans le cohésif inséré. C'est exactement ce que produit `dn0 = min(sig, ft)/pj` : le joint né en traction démarre au pic, l'adoucissement commence immédiatement. **Conforme.**
- « nouveau schéma de mise à jour nodale » = le pendant du dédoublement de nœuds de leur fig. 7. rockim l'implémente **en dual** (nœuds déjà dupliqués, liaison par union-find, groupes qui ne font que se scinder). Ce mécanisme de *node binding* fait par ailleurs l'objet d'un article compagnon du même groupe : *Implementation of extrinsic cohesive zone model (ECZM) in 2D FDEM using node binding scheme*, **Computers and Geotechnics** (2023, `S0266352X23002276`) `[VERIFIE — titre/revue via résultats de recherche]`. **À récupérer : c'est probablement là que se trouvent les détails d'implémentation manquants dans l'IJRMMS**, y compris d'éventuelles gardes de débit.

### B.3 Cas de validation `[VERIFIE-résumé]`

- un **exemple de mécanique des milieux continus à solution analytique** (précision + efficacité) ;
- **disque brésilien** ;
- **compression triaxiale**.

À noter : ce sont trois essais **de laboratoire, chargement quasi-statique monotone, champ à forte localisation imposée par la géométrie**. Aucun cas de type tunnel profond, aucun cas d'impact/percussion, aucun cas à **zone plastique étendue au seuil**. `[MEMOIRE]` complète : UCS ~51 MPa reproduit à 47,8-51,1, φ 22,8° vs 22,87, gain dt ×8-18 vs 100 E.

**Le domaine de validation publié ne couvre pas ton cas d'usage.** C'est le point central : rien dans l'article ne démontre le comportement du critère quand la zone critique est *aréale* et non *linéaire*.

### B.4 Limites — ce que dit la littérature `[VERIFIE-résumé]`

Aucune limite explicite n'est extractible de l'abstract. Les limites génériques de la famille sont documentées :

- les cohésifs extrinsèques restent **dépendants du maillage** en l'absence de stratégie d'adaptation, la fissure ne pouvant se propager que le long des arêtes ;
- pathologies connues : **propagation parasite des ondes élastiques**, **effets de vitesse de pointe de fissure**, effets de maillage ;
- l'avantage revendiqué est la suppression de la **complaisance artificielle** pré-rupture — ce que rockim confirme quantitativement (E apparent 99,1 % vs 95,6 %) `[MEMOIRE]` ;
- pour la stabilité : la pénalité `k⁻` grande rétrécit directement Δt ; trois mécanismes d'instabilité identifiés — raideur cohésive initiale élevée, **saut de raideur à la transition cohésif↔contact**, adoucissement (numériquement bénin) ; pathologie observée : **le nombre de fragments croît sans borne après l'événement de rupture principal au lieu de plafonner**, avec dérive exponentielle d'énergie `[VERIFIE]` (arXiv 2511.14323). **C'est très exactement ta signature.** Remèdes proposés là-bas : surveillance locale des « points chauds » où `k⁺` dépasse un seuil admissible ; pénalité adaptative `k⁻ = k⁺(D)` pour supprimer la discontinuité de raideur.

### B.5 Verdict

**rockim est fidèle au schéma publié sur tout ce qui est publié, avec deux réserves et un angle mort.**

Fidèle : critère en moyenne des deux tenseurs voisins projetée sur l'arête ; OU traction/Coulomb ; balayage à chaque pas sans limite d'insertions ; loi purement post-pic ; continuité de contrainte à l'insertion (aire restante = Gf exactement) ; scission progressive des sommets par union-find ; aire de joint et pénalité `α E/h`.

Réserve 1 — **enveloppe** : le défaut `shearEnvelope ≠ yang` supprime la réduction de `fs` en traction, ce qui n'est pas la formule que tu attribues à l'article. Une clé à poser explicitement dans tous les .cfg concernés.

Réserve 2 — **pénalité** : `insertionPenaltyFactor = 4` est un choix propre à rockim, hors de la fourchette usuelle. Non fautif, mais non justifié par l'article.

Angle mort — **l'article ne dit rien du débit d'insertion** parce que ses trois cas de validation ne le sollicitent jamais. rockim hérite donc d'un « pas de règle » qui n'a jamais été un « pas de règle validé ».

---

## C. Diagnostic du tapis d'insertion, et remèdes

### C.0 Contrôle préalable obligatoire — 28 000 insérés ≠ 28 000 cassés

`nInserted_` compte les **activations**, pas les ruptures. Dans le schéma intrinsèque, « joint cassé » signifie **D = 1** ; dans l'adaptatif, un joint inséré à `dn0 = dnE` a **D = 0**. Comparer 28 000 insertions à ~10 000 ruptures MultiFracS est un faux rapprochement. **Refais la comparaison sur un observable commun** : nombre de joints à `D ≥ 0,99`, ou mieux **énergie de rupture dissipée cumulée** (∑ Gf·A des joints, seul invariant physique). Il est parfaitement possible que les 18 000 « en trop » soient des arêtes **déverrouillées mais non ouvertes**, cosmétiquement visibles en nuage dans les VTU. Ce test coûte un dépouillement, pas un run.

### C.1 Six mécanismes, du plus au moins probable

**M1 — Seuils rigoureusement uniformes (`jointWeibullM = 0` par défaut).** Champ lisse + seuil unique = franchissement synchrone de toute la couronne plastique. C'est le « frein physique » de Zhou & Molinari qui manque. **Probabilité : très élevée.** Signature attendue : l'instant de première insertion est net et le compteur explose en quelques pas.

**M2 — La garde de continuité de contrainte supprime le blindage.** C'est le mécanisme le plus contre-intuitif et il découle d'une implémentation *correcte*. Dans un CZM extrinsèque classique (Camacho-Ortiz), l'insertion **relâche** la traction d'un coup : une onde de décharge part et **protège les voisines**. Ici `dn0 = min(σ,ft)/pj` fait transmettre **exactement** la traction précédente : l'insertion est **invisible** pour le voisinage à l'instant où elle a lieu. La relaxation ne peut venir que de l'adoucissement, qui exige une ouverture `dnF − dnE = 2Gf/ft ≈ 14 µm` — soit **des dizaines à des centaines de pas**. Pendant tout ce temps, la couronne continue de franchir le seuil. **Le tapis est la conséquence logique de la garde de continuité combinée à M1.** Corollaire méthodologique : ne casse **pas** la garde (elle protège Gf) — mets le frein ailleurs.

**M3 — Complaisance normale ajoutée par une pénalité basse en compression.** Avec `pj = 4 E/h`, chaque arête coupée insère un ressort de 4× la raideur d'un élément dans la direction normale. Sous 40 MPa de pression de cavité et sur une couronne entièrement coupée, cela **assouplit l'anneau**, augmente la convergence, donc la déformation, donc les insertions : boucle de rétroaction positive. En impact, même effet dans la zone broyée. **Probabilité : moyenne-élevée, effet cumulatif.**

**M4 — Le lissage par moyenne des deux CST détruit le pouvoir de sélection.** Un CST n'a aucun gradient interne ; la moyenne de deux CST est un opérateur de moyenne sur un patch de 2 éléments. Toutes les arêtes de l'éventail autour d'un sommet voient une valeur quasi identique → **elles s'insèrent ensemble** → le sommet éclate en 6 copies libres d'un coup. Motif géométrique attendu : des **étoiles** de fissures courtes autour des sommets, pas des lignes. Compare ce motif à tes VTU : s'il est là, M4 est confirmé.

**M5 — La branche de cisaillement en compression.** Dans une couronne plastique de Mohr-Coulomb, `|τ| = c + tanφ|σn|` est atteint **partout simultanément** : c'est la définition même de la zone plastique. Le critère ne peut alors **rien sélectionner** — il « peint » la zone plastique. Le nuage diffus serait donc majoritairement composé d'insertions en **cisaillement**, pas en traction. **Test décisif, à faire absolument.**

**M6 — Aucune distinction nucléation / propagation.** Rien ne privilégie l'arête en pointe de fissure existante. Sur les cas de laboratoire de l'article la singularité de pointe suffit à sélectionner ; sur un champ déjà saturé, elle ne suffit plus.

### C.2 Instrumentation à ajouter (pure addition, principe VIII)

Avant tout remède, **mesure**. Trois champs à ajouter dans `Joint` et à écrire dans les VTU de joints à côté de `bonded`, remplis dans `activateJoint` (qui reçoit déjà `sig` et `tau`) :

- `J.modeIns` ∈ {traction, cisaillement, les deux} — quelle branche a déclenché ;
- `J.sigIns`, `J.tauIns` — l'état au déclenchement ;
- `J.tIns` — l'instant, pour tracer l'histogramme temporel des insertions.

Plus, dans le résumé de fin de run : **nombre d'insertions par pas** (max et médiane), et **répartition traction/cisaillement**. Un seul run instrumenté tranche entre M1, M4 et M5.

### C.3 Remèdes, par rapport efficacité/coût

**R1 — Hétérogénéité de seuil corrélée. Coût nul, déjà implémenté.**
`jointWeibullM = 6` (voire 4-8 pour un granite) **plus** `strengthCorrLength` = une longueur physique de grain/microfissuration (`FdemSolver.cpp:1756`, avec `strengthCorrLengthB` et `fieldSeed`). Le champ **corrélé** est très supérieur au Weibull indépendant ici : il crée des **zones faibles connexes** de la taille de la longueur de corrélation, qui sont les germes de macro-fissures, là où un Weibull i.i.d. ne fait que bruiter l'ordre de franchissement sans créer de chemin préférentiel. Le fichier `configs/tunnel_bore_corr.cfg` fait déjà exactement ça (`m = 6`, `ℓ = 6 mm`, `fieldSeed = 555`) — **mais ne pose pas `insertion = adaptive`** : la comparaison corrélé/adaptatif n'a apparemment jamais été faite. C'est le premier run à lancer. Pose aussi `weibullScope = strengthGf` pour que Gf suive.

**R2 — Ordonner par sévérité et brider le débit. Coût : ~15 lignes, opt-in.**
Remplacer le tri par `jI` par un tri sur le **dépassement relatif** `max(σn/ft, |τ|/fs)` décroissant, puis n'activer qu'un sous-ensemble : soit `maxInsertionsPerStep` (plafond absolu), soit `insertionQuantile` (n'activer que les q % les plus critiques). Le déterminisme inter-threads est préservé en départageant les ex æquo par `jI`. Effet physique : la première arête ouverte a le temps de commencer à adoucir avant que ses voisines ne soient coupées — c'est **exactement** l'avancement arête par arête de Camacho & Ortiz, transposé à un balayage global. **Attention** : un plafond trop bas introduit une dépendance au pas de temps ; préférer le quantile ou un seuil de dépassement (voir R3).

**R3 — Hystérésis nucléation / propagation. Coût : ~10 lignes, opt-in.**
Deux seuils : une arête **adjacente à un joint déjà inséré** s'insère à l'enveloppe nominale ; une arête isolée exige un **surdépassement** `insertionNucleationFactor` ∈ [1,05 ; 1,3]. L'information est déjà disponible gratuitement (`jointsOfVert_` donne le voisinage de chaque sommet, et `J.bonded` l'état). C'est le remède le plus ciblé contre M6 et il transforme mécaniquement un nuage en lignes.

**R4 — Plafond topologique par sommet. Coût : ~10 lignes, opt-in.**
Interdire plus de *k* insertions par sommet et par pas (k = 1 ou 2). Empêche directement l'éclatement en étoile de M4, et force les fissures à progresser comme des chemins. Utilise `jointsOfVert_`, déjà construit.

**R5 — Relever la pénalité. Coût : une clé, mais Δt en baisse.**
Balayage `insertionPenaltyFactor` ∈ {4, 20, 100} sur le cas tunnel. Si le compteur d'insertions chute nettement, M3 est confirmé. Δt varie en 1/√pj : ×5 sur pj coûte ×2,2 sur le mur. Variante fine, dans l'esprit du remède « pénalité adaptative » de la littérature `[VERIFIE]` : **pénalité asymétrique**, élevée en compression (là où elle remplace le continuum) et modérée en traction (là où l'adoucissement domine de toute façon).

**R6 — Fermer l'écart d'enveloppe.** Poser `shearEnvelope = yang` explicitement pour être sur la formule de l'article, et le documenter dans le .cfg. À faire dans le même run que R1 pour ne pas mélanger les effets.

### C.4 Plan de runs proposé (4 runs, à valider avant lancement)

| # | Config | Question tranchée |
|---|---|---|
| 0 | rerun instrumenté du cas actuel (C.2) | Traction ou cisaillement ? Combien d'insertions par pas ? Histogramme de D ? |
| 1 | + `jointWeibullM = 6`, `strengthCorrLength`, `weibullScope = strengthGf` | R1 suffit-il à localiser ? |
| 2 | + `insertionPenaltyFactor = 20` puis `100` | M3 pèse-t-il ? |
| 3 | + R3 (hystérésis) et/ou R4 (plafond par sommet) | Le frein topologique fait-il apparaître les blocs ? |

Le run 0 est indispensable : il peut à lui seul annuler le diagnostic (cf. C.0) et éviter les trois autres.

---

## Sources

- [A 2D adaptive finite-discrete element method for simulating fracture and fragmentation in geomaterials — ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1365160923001132) `[VERIFIE-résumé — page 403 en accès direct]`
- [A 2D adaptive finite-discrete element method… — NASA ADS](https://ui.adsabs.harvard.edu/abs/2023IJRMM.16905439Y/abstract) `[identification VERIFIE ; abstract inaccessible, 405]`
- [Implementation of extrinsic cohesive zone model (ECZM) in 2D FDEM using node binding scheme — Computers and Geotechnics](https://www.sciencedirect.com/science/article/abs/pii/S0266352X23002276) `[VERIFIE — titre/revue ; contenu non lu — À RÉCUPÉRER EN PRIORITÉ]`
- [Stability of Extrinsic Cohesive-Zone Model with Penalty-Based Contact in Explicit Dynamic Fragmentation Simulations — arXiv 2511.14323](https://arxiv.org/html/2511.14323) `[VERIFIE — page HTML lue intégralement]`
- [Numerical modeling of crack propagation with dynamic insertion of cohesive elements — Engineering Fracture Mechanics](https://www.sciencedirect.com/science/article/abs/pii/S0013794419309324) `[non lu]`
- [Chengzeng Yan — profil ScienceDirect](https://www.sciencedirect.com/author/56237479900/chengzeng-yan) `[non lu]`
- Code : `C:\Users\fuzquianoalricabi\simulations\FDEM\rockim\rockim_p1\src\FdemSolver.cpp` (l. 1855-2067, 2682-2712, 1668, 1745-1770), `include\rockim\YangDif.hpp` (l. 119-121), `configs\tunnel_bore*.cfg` `[CODE — lu]`


## Lacunes

LACUNES ASSUMÉES (budget 8 appels épuisé, ScienceDirect 403 / ADS 405 / colab.ws 403 / ouci sans résultat) :

1. Le TEXTE de Yan 2023 n'a pas pu être lu. Les équations 7-8 citées dans le commentaire du code n'ont PAS été revérifiées : je m'appuie sur ce que le commentaire de FdemSolver.cpp en dit et sur les résumés de moteur de recherche de la page éditeur. En particulier, je n'ai PAS pu vérifier : (a) si l'article réduit vraiment fs en traction (question shearEnvelope = yang vs défaut) ; (b) si l'article prescrit une valeur de pénalité pour le cohésif inséré ; (c) si l'article mentionne un ordonnancement, un plafond ou un quota d'insertions par pas ; (d) si l'article utilise une contrainte effective mixte plutôt qu'un OU de deux branches ; (e) si l'article discute une hétérogénéité de seuil.

2. L'article compagnon Computers and Geotechnics S0266352X23002276 (node binding scheme, même groupe, 2023) n'a été identifié que par son titre. C'est très probablement là que se trouvent les détails d'implémentation absents de l'IJRMMS — À RÉCUPÉRER EN PRIORITÉ au prochain budget.

3. Aucun papier citant Yan 2023 n'a pu être ouvert : les critiques éventuelles de la méthode n'ont pas été recensées.

4. Wang et al. 2024 / MultiFracS : aucune vérification faite (hors périmètre du volet 1 tel que formulé, mais la comparaison 28 000 vs 10 000 en dépend).

5. Aucune config de tunnel profond type Wang 2024 n'existe dans rockim_p1/configs — les tunnel_bore*.cfg sont l'analogie du banc 6 Abaqus et ne posent PAS insertion = adaptive. Le cas dont parle la question a donc été lancé ailleurs : je n'ai pas pu lire ses clés réelles (shearEnvelope, jointWeibullM, insertionPenaltyFactor, jointSoftening), et les réserves de conformité de la partie B.5 sont donc formulées comme des CONTRÔLES À FAIRE, pas comme des constats.

6. Les chiffres de compliance de M3 (25 % de la complaisance d'un élément par arête coupée) sont un ordre de grandeur calculé sur les paramètres de tunnel_bore.cfg, pas une mesure.