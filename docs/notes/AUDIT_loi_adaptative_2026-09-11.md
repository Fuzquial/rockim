# rockim_g0 face à la note « Lois constitutives pour un FDEM hybride à insertion adaptative »

Audit du 2026-09-11, sur `FDEM\rockim_g0` (exe `rockim_fix.exe`, 2026-09-07, 30 463 lignes de C++).
Document source : *Lois constitutives proposées pour un FDEM hybride à insertion adaptative —
matrice viscoplastique-endommageable et joints cohésifs extrinsèques dépendants de la vitesse*,
note de travail, septembre 2026, application à l'indentation dynamique du granite de Kuru.

Méthode : 14 briques de la note = 14 audits du source, chacun suivi d'un contre-audit adverse
chargé de le réfuter, puis une critique de complétude. 30 agents, 1 040 lectures/greps/runs.
Les mesures citées ont été refaites dans le dépôt, pas reprises de la documentation.
La brique B14 (bilan d'énergie) a échoué sur sa sortie structurée ; son contenu est couvert
par la critique §1.8 et par les vérifications directes de la section 5 ci-dessous.

---

## 0. Réponse en un paragraphe

**Oui, l'architecture de la note s'assemble et tourne aujourd'hui**, en `mode = fdem3d` avec
`insertion = adaptive`, maillage Voronoï et microstructure granulaire : run de fumée 1 854 tétraèdres
/ 3 528 facettes liées, **zéro clé de deck non lue**, résidu d'énergie 1,3e-18 J, +13 % de CPU
seulement par rapport à l'élastique. Sur 13 briques auditées : **2 conformes, 11 partielles,
0 absente**. Mais trois de ces « partielles » sont des absences sur l'équation centrale, et le
vrai problème n'est pas la physique manquante — c'est **six mécanismes dissipatifs parasites armés
par défaut** (le double comptage que la note interdit) et **l'absence totale d'instrumentation
énergétique en fdem3d**, qui rend l'affirmation « 70-90 % de l'énergie dans la zone broyée »
invérifiable dans le solveur cible.

---

## 1. Brique par brique

| § de la note | Brique | rockim_g0 | Clés |
|---|---|---|---|
| 1.2 | Enveloppe MH en loi puissance, viscoplasticité **parfaite** | **partiel** — enveloppe méridienne `q = fc0 + B·σ3^n` exacte (restituée à 1e-7 %), aucun écrouissage, return mapping sur la contrainte **effective** (Grassl-Jirásek). Mais cône **circulaire** dans le plan déviatorique, pas l'hexagone σ1−σ3 : les deux ne coïncident qu'en compression triaxiale | `law=saksala`, `meridian=power`, `merB`, `merN`, `merFc0`, **`capP0=0`** |
| 1.2 | Écoulement non associé, angle ψ | **partiel** — ψ réglable depuis le 2026-09-07, les deux formes portées (β linéaire et `sin ψ·B·n·p^(n−1)`). Mais potentiel **lisse en (p,√J2)**, pas `g = σ̄1 − m_ψ σ̄3` ; κ n'intègre que la part déviatorique | `dpDilationDeg` (défaut 0) |
| 1.3 | ω_c = Ac[1−exp(−βc ε^vp)], crack-band | **partiel** — forme exacte, régularisée par la longueur d'élément, défauts Ac = 0,98 et G_IIc = 1e4 J/m². Mais `h_e` n'est pas `(12V/√2)^(1/3)` : `V^(1/3)` en fem3d, diamètre inscrit `6V/A` en fdem3d → βc 2,0× et 5,0× trop petit. Aucune garde de snap-back | `compDamage=crackband`, `compAc`, `compGIIc` |
| 1.4 | Split spectral, **compression seule** dégradée | **partiel, et c'est le point noir** — `nominal = (1−Dc)·σ⁻ + (1−D)·σ⁺` est écrit mot pour mot (MatLaw.cpp:581), sur les directions principales effectives courantes. **Mais aucune clé ne met D à zéro** : le bloc d'endommagement de Rankine s'exécute inconditionnellement (MatLaw.cpp:483-507). La matrice dissipe **toujours** en traction | contournement `bulkFt=1e12` (qui éteint au passage `erodeD`) |
| 1.5 | Effet de vitesse `s_MH` | **conforme** — sur-contrainte linéaire `F_fin = η·dλ/dt` = la consistance de la note ; ×10 de vitesse → ×10 de sur-contrainte, mesuré. Bancs `matpoint` et `verify_fem3d_rate1/2.cfg` | `saksalaEta` [Pa·s], **pas** en unités de s_MH |
| 2.1 | Traction de facette pondérée **en volume** + taux à l'interface | **partiel/absent** — `t_n = n·σ·n`, `t_s` corrects, sur les contraintes **nominales** des deux éléments. Mais moyenne **arithmétique 0,5/0,5** (les volumes sont pourtant déjà là). Et **aucun tenseur taux** : rockim ne garde qu'un scalaire équivalent filtré par élément, donc ni `n·ε̇·n` ni `γ̇` ne sont formables | — |
| 2.2 | Critère elliptique Φ_F, DIF **gelé**, hystérésis n_h | **DIF gelé : conforme** (`strainRateDIFArm = insertion \| envelope \| continuous`, tamponné à l'activation). **Ellipse : absente** — le critère est un **OU logique** (`σn ≥ DIF·ft` **OU** `\|τ\| ≥ fs`), aux trois mêmes sites 2D-OpenMP / 2D-série / 3D. **Hystérésis : absente** (zéro occurrence de `holdSteps`). **σ*_n : jamais calculé** | `strainRateDIF`, `strainRateDIFArm` |
| 2.3 | Dédoublement des nœuds ssi l'étoile se scinde | **conforme, la seule brique pleinement conforme** — implémenté en dual : nœuds dupliqués au maillage, « scinder » = relâcher une liaison, groupes = composantes connexes par union-find sur les joints encore `bonded`. Masse au prorata des volumes (ρV/4 par copie), u et v hérités sans saut. 2D **et** 3D | automatique sous `insertion=adaptive` |
| 2.4 | TSL **initialement rigide**, `t_ins` = traction transmise | **absent** — la loi reste **intrinsèque** (raideur `pj = pf·E/h`, branche élastique `dnE = ft/pj`). Seule l'*existence* du joint est extrinsèque. À l'activation la traction est **écrêtée au seuil nominal** (`J.dn0 = min(sig, J.ft)/J.pj`), donc le joint naît au pic et le dépassement reste un saut. `δ_f` jamais recalculé depuis `t_ins` | `insertionPenaltyFactor` (défaut 4 E/h) |
| 2.5 | Option B : η_n δ̇_n après insertion | **partiel** — la forme additive normale existe et le garde `bonded` la limite bien aux joints insérés. Mais **pas de terme tangentiel**, η_n n'est pas un paramètre matériau (reconstruit depuis un taux d'amortissement adimensionnel, la pénalité, l'aire et la masse nodale), **actif par défaut** (ξ = 0,05), et rien n'interdit de l'armer en même temps que le DIF | `jointXi` (défaut 0,05) |
| 2.6 | χ_t, χ_G inter/intra + Weibull **par facette** | **partiel** — seuil par facette, avant insertion, distinguant inter/intra, Weibull par facette, **actif en fdem3d avec insertion adaptative** (vérifié : 4 786 facettes, 3 744/502/540 intra/homophase/hétérophase, 4 786 valeurs distinctes). Mais la règle inter est une **moyenne** et non un **minimum**, le Weibull est à 2 paramètres normalisé à moyenne 1 (ni x_u ni x_0 libre), et l'échelle porte sur le **volume** des deux voisins, pas sur l'**aire** de la facette | `gbAlphaTen/Coh/Gf`, `gb.<a>.<b>`, `jointWeibullM` |
| — | Mode mixte **Benzeggagh-Kenane** | **absent** — zéro occurrence de « benzeggagh » / « kenane » dans src/, include/, docs et decks. Ce qui existe est une ellipse sur les rapports de déplacement normalisés `D = √(rn²+rs²)` (Munjiza/Yan), le ratio de mode est **réévalué à chaque pas** et non gelé, et G_Ic/G_IIc pilotent **deux longueurs critiques séparées** au lieu d'une énergie mixte | `jointSoftening=yan` |
| — | Contact + frottement sur joint inséré | **partiel** — pénalité présente et réglable (`insertionPenaltyFactor`, 10 E/h atteignable par une clé), frottement de Coulomb présent. **Mais le mélange est faux** : la cohésion porte bien (1−D), le terme frottant **n'est pas multiplié par D** — il vaut μ⟨−t_n⟩ à pleine valeur **dès D = 0** | `insertionPenaltyFactor`, `frictionDeg` |

---

## 2. Ce qu'il faut écrire de zéro

1. **Le critère elliptique Φ_F et l'hystérésis n_h** (§2.2). ~250-320 lignes, 6 fichiers, à écrire
   trois fois (la branche 2D est dupliquée OpenMP/série). Aucune refonte.
2. **La TSL initialement rigide** (§2.4). ~450-550 lignes. **Mais voir §4 : c'est ×1,96 sur dt.**
3. **Benzeggagh-Kenane** et le gel du ratio de mode à l'insertion.
4. **Le tenseur taux de déformation par élément** dans le repère global, sans lequel `DIF_t(ε̇_eq)`
   et `DIF_s(γ̇)` séparés de §2.2 sont impossibles.

Tout le reste est de l'addition locale (< 100 lignes par item) ou une simple carte de deck.

---

## 3. Le double comptage : six mécanismes parasites armés par défaut

La note pose : *« un mécanisme par physique, et jamais deux fois la même dissipation »*.
Si l'on assemblait la loi aujourd'hui avec les défauts, seraient simultanément armés :

| # | Mécanisme | Doublon avec |
|---|---|---|
| 1 | **Endommagement de Rankine de la matrice**, inextinguible (MatLaw.cpp:483-507) | les joints cohésifs — **le grief central de la note** |
| 2 | **Cap volumique écrouissant** `capP0 = 8c` armé par `law = saksala` (MatLaw.cpp:4816, `s.pc += capH·dev`) | contredit « viscoplasticité parfaite », invisible en triaxial |
| 3 | **Frottement de joint pleine valeur dès D = 0** (μS = 1,0, FdemSolver.cpp:5243-5244) | le cisaillement viscoplastique de la matrice |
| 4 | **Amortisseur de joint `jointXi = 0,05`**, dont la traction visqueuse entre dans le cap de Coulomb **et dans `fnSum`** | pollue la courbe force-pénétration elle-même |
| 5 | **Amortissement local de Cundall `dampingLocal = 0,05`** en percussion 3D | non physique ; le code exige lui-même `dampingLocal = 0` sur tout banc de force |
| 6 | **Érosion fantôme** : `erodeD = 0,98` armé sans clé dans `MatLaw::make`, mais `grep -c eroded` = **0** dans les deux solveurs FDEM (29 en fem3d) | ni dissipation ni suppression : l'élément garde masse, joints et nœuds, et continue d'entrer avec un poids 0,5 dans la moyenne de facette — un voisin survivant voit sa traction **divisée par deux** au moment où il porte seul la charge |

Point rassurant vérifié : `crushCap`, `meanTensionCapFactor` et `bulkDamage` sont neutralisés dès
que `law_` existe (Fdem3dSolver.cpp:2721, 2727, 2736).

---

## 4. Trois faits que la note ne pouvait pas anticiper

**(a) Une TSL vraiment rigide vaut ×1,96 sur le pas de temps.** `computeStableDt` charge la pénalité
de **tous** les joints, y compris les 3 528 encore liés et jamais évalués :

| schéma | dt |
|---|---|
| `insertion = intrinsic` (pf = 20) | 1,82e-9 s |
| `insertion = adaptive` (pf = 4) | 4,02e-9 s |
| `insertion = none` (≈ TSL rigide) | **7,86e-9 s** |

C'est la moitié du temps CPU de chaque percussion 3D. Le plus gros levier économique du dossier.

**(b) La viscosité de la matrice s'éteint exactement là où la note en a besoin.** Sur-contrainte
linéaire **sans aucune saturation**. Mesures (Red Bohus, `law = saksala`, `meridian = power`,
`capP0 = 0`, Rankine neutralisé, σ3 = 20 MPa, η = 0,05e6 Pa·s, à ε_ax = 1 %) :

| ε̇ | q | ε_vp équivalente | travail plastique |
|---|---|---|---|
| 1 /s | 406,6 MPa | 4,76e-3 | 1,94 MJ/m³ |
| 10³ /s | 663,4 MPa | 1,46e-3 | 0,86 MJ/m³ |
| 10⁴ /s | 762,2 MPa | 1,85e-4 | 0,12 MJ/m³ |
| 10⁵ /s | 775,1 MPa | 1,9e-5 | **0,0125 MJ/m³** |

À 10⁵ /s le travail plastique est divisé par **155**. Comme ω_c est piloté par `epvEq`,
**l'endommagement de compression s'éteint avec lui**. Sous un insert (ε̇ local 10⁴-10⁶ /s), le
mécanisme censé porter 70-90 % de l'énergie devient quasi élastique-rigide. C'est une objection de
physique, pas de code : il faut soit recalibrer `s_MH` sur 10³-10⁶ /s, soit ajouter une saturation.

**(c) La séquence de calibration §3.3 est impossible dans le solveur cible.** `fdem3d scenario must
be percussion | shear | tension` (Fdem3dSolver.cpp:386) — ni brésilien, ni SHPB, ni UCS ; et
`loading = grips | platens` n'est enregistré que `{"loading", "fdem"}`, donc **2D seulement**.
Les paramètres de **joint** (ft, Gf, c, GfII, gbAlpha*) ne peuvent être calibrés qu'en 2D ou en
traction directe 3D. À trancher explicitement, pas par omission.

---

## 5. L'instrumentation : le verrou structurant

Vérifié directement :

- **fdem3d, champs d'élément VTU** : `vonMises, sigma1, tauMax, fragment, phase, grain` (+ `bulkD`).
  **Ni D, ni ω_c, ni ε^vp, ni travail plastique.**
- **fdem3d, `history.csv`** : `…, eEl, eJnt, eGc, eFric, eCund, eLys` — **un seul poste `eEl`
  pour tout le volume**, et le code avoue lui-même que l'énergie élastique stockée est « exacte en
  élastique, approx sous law/caps » (Fdem3dSolver.cpp:4820-4821).
- **fem3d écrit tout** (`damage, eroded, erodedBy, epvEq, omegaC`) — mais n'a **aucun joint**.

> Le solveur bien instrumenté n'a pas la physique de la note ; le solveur qui a la physique n'a pas
> l'instrumentation. On ne peut aujourd'hui **ni démontrer les 70-90 %, ni prouver l'absence de
> double comptage**.

Le bilan `W_bit = E_el + E_kin + D_vp + D_ωc + D_coh + D_fric + E_abs + E_art` de §3.2 est donc
**hors de portée en fdem3d** tant que `wPlas / wDamT / wDamC` ne sont pas récoltés dans
`Fdem3dSolver::finalize` (le patron existe déjà en fem3d, Fem3dSolver.cpp:2597).

Ce qui existe déjà côté scénario : `toolShape = sphere | flat | pdc | none`,
`toolContact = penalty | signorini`, outil libre avec masse et `impactSpeed`, `scenario = percussion`
— l'indentation par bouton sphérique visée par la note est disponible. Mais `toolSignoriniGroup`,
le correctif qui fait dialoguer l'impulsion de Signorini avec la liaison de nœuds par groupes,
est enregistré `{"toolSignoriniGroup", "fdem"}` : **2D seulement**. En fdem3d l'impulsion est
dimensionnée par copie de nœud alors que l'intégrateur travaille par groupe.

---

## 6. Ordre de travail proposé

0. **Instrumentation d'abord** (~60-80 lignes) : récolter `wPlas / wDamT / wDamC` dans
   `Fdem3dSolver::finalize`, colonnes dans `history.csv`, champs `damage / omegaC / epvEq` au VTU.
   Sans cela **aucune** étape suivante n'est falsifiable.
1. **Hygiène des défauts** (0 ligne de C++) : `capP0 = 0`, `dampingLocal = 0`, `jointXi = 0`,
   `erodeD = 2`, `bulkFt = bulkGf = 1e12`, `absorbing = all`, `toolContact = signorini`.
   Plus trois refus explicites (~5 lignes chacun), dont `erodeD/erodeEpv` en FDEM.
2. **`bulkTensionDamage = off`** (~60-100 lignes, un seul fichier) : tuer le double comptage n°1.
3. **Facteur D sur le frottement de joint** (~15 lignes/solveur + le clamp d'activation).
4. **Sortir le terme visqueux du critère** (`tauLim` sur `sigEl`, ~2 lignes opt-in).
5. **Trancher la viscosité** (§4b) : recalibrer `s_MH` ou ajouter une saturation.
6. **TSL initialement rigide** — reclassée : pas du formalisme, **×1,96 sur dt**, avant toute
   campagne longue.
7. **Calibration** : porter `loading = platens` en fdem3d, ou assumer la calibration 2D des joints.

**Vraiment plus tard** (aucun ne bloque un run ni ne fausse un bilan d'énergie) : forme principale
hexagonale de MH, potentiel `g = σ̄1 − m_ψ σ̄3`, convention `h_e`, moyenne pondérée en volume à la
facette, critère elliptique et hystérésis, Benzeggagh-Kenane, Weibull à 3 paramètres.

---

## 7. Réserve de méthode

Cinq erreurs ont été corrigées par les contre-audits, dont une affirmation répétée par **quatre**
briques sur treize : « `law` et `phases` sont exclusifs, donc la matrice MH ne peut pas cohabiter
avec une microstructure GBM ». **C'est faux** : la garde est `if (phases_.n() > 1) throw`
(Fdem3dSolver.cpp:423), c'est-à-dire le **multiphase** qui est refusé, pas le Voronoï monophase.
Le run de fumée le prouve : `bulk law = saksala, 1854 tets` + `adaptive insertion: 3528 bonded faces`
+ `voronoi: 53 grains, 1 phase(s), 747 grain-boundary joints`, zéro clé non lue.

Une hétérogénéité **entre phases minérales** (quartz / feldspath / mica) reste donc hors d'atteinte
avec `law = saksala` ; l'hétérogénéité **inter/intra-granulaire** de §2.6, elle, fonctionne.

Enfin : la note vise le granite de **Kuru** (E = 60 GPa, σc0 = 235 MPa, Hokka 2016) ; la campagne
du projet porte sur **Red Bohus** (E ≈ 77,7 GPa). Les valeurs de départ du tableau §3.4 sont à
retransposer.
