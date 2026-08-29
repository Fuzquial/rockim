# LOT 5 — Le FDEM d'Imperial et rockim : synthèse et plan de travail
*Livrable final de [MISSION_etat_de_l_art_2026-08-29.md](MISSION_etat_de_l_art_2026-08-29.md).
Rédigé le 2026-08-29. Commanditaire : F. Uzquiano.*

---

## 0. CE QUE CE DOCUMENT VAUT

### 0.1 Ce sur quoi il repose

**Neuf sources primaires**, lues de première main sur PDF fournis par
F. Uzquiano — la totalité de la littérature d'impact du groupe d'Imperial, plus
sa source de formulation et deux livrables ORCHYD :

| code | référence |
|---|---|
| **UCL** | Guo L., Xiang J., Latham J.-P., Izzuddin B., « A generic computational model for three-dimensional fracture and fragmentation problems of quasi-brittle materials », manuscrit déposé (UCL Discovery) — **la formulation, éq. 1-20** |
| **XLF** | Xiang J., Latham J.-P., Farsi A., « Algorithms and Capabilities of Solidity… », DEM7 (Dalian, 2016), *Springer Proc. Phys.* **188** ch. 16 — **le frottement tangentiel, éq. 3-5** |
| **IJ191** | Yang X., Xiang J., Naderi S., Wang Y., Aising J., Ugarte I., Latham J.-P., *IJRMMS* **191** (2025) 106125 — **St Anne + Rhune, 7 critères** |
| **IJ206** | Yang X., Xiang J., Naderi S., Wang Y., Aising J., Ugarte I., Latham J.-P., « High-fidelity modelling of fragmentation and pulverisation in hard granite under percussion loading: a FDEM-based approach », *IJRMMS* **206** (2026) 106660 — **la pulvérisation** |
| **JR** | Yang X., Xiang J., Latham J.-P., Naderi S., Wang Y., *JRMGE* (2025) — **Kuru Grey** |
| **ANN** | Naderi S. *et al.*, *JRMGE* (2025) — **réseau de neurones sur données FDEM, modèle 2D** |
| **A952** | Yang X. *et al.*, **ARMA 24-0952**, « Where does the energy go in percussion drilling? FDEM's answer » — **le bilan d'énergie sur St Anne** |
| **A788** | Gerbaud L. *et al.*, **ARMA 24-0788**, « Can DTH Hammer Drilling Deliver? » — **confinement jusqu'à 130 MPa ; nomme « SOLIDITY »** |
| **ORCHYD** | livrables **D6.1** et **D6.4**, plus Dumoulin *et al.*, *Rock Mechanics Bulletin*, DOI 10.1016/j.rockmb.2024.100169 — **Red Bohus, modèle continu CDP** |

### 0.2 Discipline de preuve, et son contrôle

Trois étiquettes, jamais confondues : **[ÉTABLI]** lu dans la source, citation et
page à l'appui · **[INFÉRÉ]** déduit, raisonnement donné · **[OUVERT]** non
tranché par les sources.

Le dépouillement a été **contre-vérifié de façon adversariale** : un second
lecteur, ignorant du premier, retournait à la page et cherchait la
surinterprétation. **673 verdicts rendus : 614 confirmés (91,2 %), 54 surinterprétations,
3 faux, 2 introuvables.** Les 59 défauts ont été corrigés avant d'entrer ici.
Les taux ne sont pas décoratifs : ils disent qu'environ **une affirmation sur
onze** d'une première lecture allait plus loin que sa citation. C'est la mesure
du piège que la CORRECTION 2 avait identifié. *(Chiffres définitifs : le
dépouillement a terminé à 28 agents sur 28, sans erreur.)*

### 0.3 Ce que ce document ne couvre pas

* **L'audit du code de rockim est partiel** : 4 éléments sur 11 audités et
  vérifiés ligne à ligne, 7 en première lecture, **aucun contre-audit** — une
  limite de session les a coupés. Le [lot 4](biblio_insertion/2026-08-29_lot4_bilan_rockim.md) §0 et §7 le disent.
* **Le premier lot a été bâti sans accès web** au-delà d'un moteur de recherche :
  toute la bibliographie du [lot 1](biblio_insertion/2026-08-29_lot1_bibliographie_imperial.md) est en **[MÉTA]**, jamais en [LU].
  Les neuf sources ci-dessus l'ont depuis largement dépassé.

---

## 1. LE PAYSAGE

**[ÉTABLI]** Le FDEM n'est pas un code, c'est une famille à **trois lignées** :
**Imperial + QMUL** (Y2D/Y3D → VGW → VGeST → **Solidity**, relancé sous ce nom en
2016), **Los Alamos** (HOSS : Rougier, Knight), **Toronto** (Y-Geo, Irazu). Munjiza
co-signe les trois à des époques différentes : **sa signature n'établit pas qu'un
article décrit Solidity.**

**[ÉTABLI]** Le code utilisé pour l'impact est nommé une seule fois dans tout le
corpus — A788 p. 1 : « our in-house rock fracture software, **SOLDITY** [sic]
that leverages a hybrid FDEM method ».

---

## 2. LA FORMULATION D'IMPERIAL — ce qui est ÉTABLI

| élément | contenu | source |
|---|---|---|
| **loi de joint** | σ à trois branches : compression linéaire à pente **2·pj**, durcissement **parabolique** à tangente nulle au pic, adoucissement **z·ft** | UCL éq. 11 p. 13 |
| **fonction z** | heuristique à **a = 0,63, b = 1,8, c = 6,0**, « derived originally for concrete », reprise faute de données spécifiques | UCL éq. 12 p. 13 |
| **D mixte** | déplacement **normalisé**, combinaison quadratique des modes I et II ; D = 1 dès qu'un mode atteint son δ critique | UCL éq. 13 p. 13 |
| **résistance au cisaillement** | f_s = c − σ_n tan φ, tension cut-off automatique | UCL éq. 5 p. 11 |
| **ouverture au pic** | **δ_np = 2 f_t h / p₀**, h = **longueur moyenne des arêtes du joint** | UCL éq. 6-7 p. 12 |
| **énergie de rupture** | G_f ≈ ⅓ f δ_c (aire de la parabole approchant l'exponentielle) | UCL éq. 10 p. 12 |
| **pénalité** | **p₀ est une CONTRAINTE**. Recommandation **E ≤ p₀ ≤ 10E**, choisie pour le pas de temps ; **souplesse artificielle assumée, module d'Young NON corrigé** | UCL éq. 8-9 p. 12 |
| **pénalité pratiquée** | **3000 GPa pour E = 57 GPa (St Anne) = 52,6 E** ; Rhune 1800/37 = 48,6 E — **5× au-dessus de leur propre borne** | A952 Table 1 p. 4 |
| **rupture du joint** | « labelled as failed when **at least two integration points** have zero stress components » (2 sur 3) | UCL p. 14 |
| **contact** | détection **NBS** ; interaction par **potentiel de Munjiza**, force distribuée | UCL éq. 17 p. 16 |
| **détection paresseuse** | aucun couple dans le continu ; ils naissent **à la rupture d'un joint**, par six groupes nodaux ; les couples à joint **actif** sont exclus | UCL éq. 14-16 pp. 14-15 |
| **rampe de naissance** | après une fissure de **cisaillement sous compression**, la force monte linéairement sur **n_total ≈ 10 pas** | UCL éq. 18 p. 17 |
| **frottement** | **f_t = −k_t δ_t − η v_t**, plafonné à −µ f_n dès que f_t ≥ µ f_n | XLF éq. 4-5 p. 4 |
| **DIF** | deux lois par morceaux ; **compression → cohésion et G_II**, **traction → f_t et G_I** ; frottement interne intact | IJ206 éq. 1-2 p. 3 |
| **pulvérisation** | σ = (1−D)σ̄ au-delà de ε_d, D piecewise en δ_m (une **longueur**), plafonné à D_max ; **assignée aux ÉLÉMENTS, « not intended for joint elements »** | IJ206 éq. 3-4 p. 4 |
| **retrait des fragments** | identification par **connectivité**, puis **post-traitement anti-gravité** (v₀ = 2,5 mm/s) ; **aucune érosion en cours de calcul** | IJ191 §2.3 ; IJ206 §2.1 |
| **bilan d'énergie** | six postes ; **amortissement et erreur obtenus PAR SOUSTRACTION** | A952 éq. 3-7 p. 3 |
| **partage mesuré** | St Anne, piston 9,41 m/s, 49,3 J : **2,6 % à la fissuration, 64,9 % au frottement** ; Rhune 2,46 % et 62,2 % | A952 §4 pp. 5-7 |
| **insertion** | **intrinsèque à 100 %**, nœuds dédoublés **avant chargement** | UCL p. 15 |
| **périmètre** | **la roche seule** ; « except in the insert, which is treated as a **fracture-free component** » | ANN p. 6870 |
| **maillage** | **fixe** ; fissures uniquement le long des faces de tétraèdres ; raffinement décidé **a priori** | UCL p. 14 ; IJ191 p. 4 |

### 2.1 Les trois aveux qui valent plus que les équations

1. **Leur 65 % de frottement est gonflé par leur propre maillage.** A952 p. 5 :
   « for the model used in this study, **rock fragments smaller than 1 mm cannot
   be further fractured** […] the **friction energy output by FDEM simulation
   should include some of the fracture energy** ». Plus l'angularité : « Because
   of the relatively **sharp tetrahedral elements** […] a **large friction
   force** will be generated ». **Le 2,6 % de fissuration est un plancher, pas
   une mesure.**
2. **Le modèle de pulvérisation est superflu sur le calcaire.** IJ206 §2.2 p. 4 :
   pour St Anne et Rhune, « previous validated FDEM studies were able to
   reproduce the main fragmentation characteristics **without introducing the
   additional damage model or modified sliding friction treatment** ».
3. **Leur coefficient de frottement glissant joue DEUX rôles, dont un numérique.**
   JR p. 5 : « Since FDEM simulations rely on mesh boundaries for cracking, the
   generated rock fragments **always have sharp edges** […] **it is necessary to
   reduce the sliding friction coefficient** » — c'est la compensation de
   l'angularité des tétraèdres. Mais JR p. 10 ajoute un rôle **physique** :
   « The **abrupt change in the friction coefficient from intact internal to crack
   wall sliding (i.e. from 1.96 to 0.39)** leads to a **fast energy accumulation**
   as would be required to **initiate radial cracks**. » Le **contraste** entre
   frottement interne et glissement est un moteur d'amorçage. Dans les deux cas,
   la valeur reste non transférable — un contraste dépend des deux termes.
4. **Dans JR, l'effet de vitesse est absorbé dans G_I et G_II**, pas appliqué par
   le DIF : « The energy release rates, G_I and G_II, were **artificially
   increased** to consider the loading rate effect through the validation
   process » (p. 4). **Leurs G ne sont donc pas des propriétés quasi-statiques**,
   et ne se comparent pas à ceux d'IJ206, qui accompagnent un DIF explicite. Même
   symbole, deux sens.

### 2.2 Et la preuve de la troisième : le même granite, calibré deux fois

**[ÉTABLI]** Kuru Grey, même groupe, un an d'écart :

| | JR (2025) | IJ206 (2026) |
|---|---|---|
| E (GPa) | 67 | **60** |
| G_I (J/m²) | 20 | **50** |
| G_II (J/m²) | 1500 | **1000** |
| cohésion (MPa) | 46,49 | **29,84** |
| **µ glissant** | **0,39** | **0,18** |

**Aucun paramètre d'Imperial n'est transférable** — ni entre roches, ni entre
codes, ni entre deux de leurs propres articles. Leur propre avertissement le dit :
ils calibrent « as an **integrated parameter set** » (IJ206 p. 6).

---

## 3. CE QUI EST INFÉRÉ

* **[INFÉRÉ]** L'insertion est intrinsèque : le mot n'est jamais écrit, mais les
  nœuds sont dédoublés « before loading starts » (UCL p. 15) — un schéma
  extrinsèque ne pourrait pas définir la connectivité ainsi.
* **[INFÉRÉ]** Les joints sont absents de l'acier et du carbure : leurs tables
  omettent **exactement les cinq paramètres** qui définissent un joint (G_I,
  G_II, f_t, c, tan φ).
* **[INFÉRÉ]** Il n'y a **aucun maillage adaptatif** : UCL p. 41 range le
  raffinement localisé en **perspective de recherche**. On ne met pas en
  perspective ce qu'on a implémenté.
* **[INFÉRÉ]** L'écart entre la recommandation (E-10E) et la pratique (≈50 E)
  s'explique par le régime : la borne basse vise le quasi-statique, l'impact
  exige davantage. **Non écrit.**
* **[INFÉRÉ]** L'algorithme tangentiel remonte à **Xiang, Munjiza & Latham
  (2009), *IJNME* 79(8), 946-978** : XLF écrit « Xiang et al (2009) developed
  further the FEMDEM method by taking account of the sliding friction force », et
  la thèse de Guo attribuait l'implémentation à « Dr Jiansheng Xiang ».

---

## 4. CE QUI RESTE OUVERT

| question | pourquoi elle reste ouverte |
|---|---|
| **le couplage endommagement → frottement** | IJ206 décrit l'effet (« reducing **their** sliding friction coefficient », pp. 4 et 11) mais **ne publie aucune équation**. Sa Table 1 ne porte qu'une valeur par matériau. Impossible de trancher entre µ(D) par élément et µ calibré bas globalement. |
| **la valeur de k_t et de η** | non publiées, nulle part |
| **la règle pour une paire de matériaux différents** | non publiée sur **six** sources. JR Table 2 donne pourtant granite 0,39 et acier 0,1 : une règle est **nécessairement** exercée par leur code, ils ne l'écrivent pas. |
| **le référent du « Penalty Number »** | joint, contact, ou les deux ? A952 n'en donne qu'un |
| **l'unité du Mass Damping Coefficient** | non imprimée |
| **C_d et ε_d** du modèle d'endommagement | « a material-dependent constant », sans valeur |
| **le pas de temps d'IJ191** | jamais publié — l'article de validation le plus détaillé du corpus n'est pas reproductible sur ce point |
| **l'objectivité au maillage** | leur Table 2 de balayage ne publie **que** le nombre d'éléments et le temps CPU. Aucune courbe de convergence. Les mots *objectivity*, *regularisation*, *characteristic length* : zéro occurrence. |

---

## 5. LE BILAN DE ROCKIM

**[ÉTABLI par lecture du code]** — détail au [lot 4](biblio_insertion/2026-08-29_lot4_bilan_rockim.md).

### 5.1 Ce qui est conforme

La loi de joint est **complète et fidèle**, opt-in par clés : z-curve aux mêmes
constantes (`YanSoftening.hpp:53-68`, `munjiza` étant un **alias** de `yan`),
trois branches (`FdemSolver.cpp:3987-4005`), D mixte (`:3958, 3917, 3968`),
irréversibilité (`:3970`), **règle des 2 points sur 3** (`Fdem3dSolver.cpp:2777-2783`).
Le potentiel de Munjiza est là, conservation vérifiée à 3,7e-12 (2D) et 2,0e-8
(3D). Le DIF est implémenté **au-delà** de l'article. L'exclusion des couples à
joint vivant est là.

### 5.2 Les deux seuls blocages

1. **La longueur de référence h.** rockim mesure le **diamètre inscrit** 6V/A
   (`Fdem3dSolver.cpp:1392`) ; Imperial la **longueur moyenne des arêtes**
   (UCL p. 12). Rapport **2,4495** pour un tétraèdre régulier. **Le deck de
   réplique, à 26,32, est ≈ 2,45 fois trop raide** ; l'équivalence correcte est
   26,32 × 0,4082 = **10,74**.
2. **L'injection d'énergie par le contact.** Le dépôt a mesuré **+3,66 J sur
   l'impact 3D P1 (6,9 % de KE₀) et 11,1 J (20 %) en insertion intrinsèque**.
   Tant qu'un canal injecte un cinquième de l'énergie d'entrée, **le partage
   2,6 / 64,9 % est hors d'atteinte**.

### 5.3 Le problème d'intégrité

> **⚠️ PRÉMISSE FAUSSE — voir `chantier_imperial_2026-08-29/A03_resourcer_attributions.md`.**
> `solidity-solver-open` est le dépôt public d'Imperial College London (LGPL-3.0),
> cloné et lu le 2026-08-26, provenance documentée en quatre endroits du dépôt.
> Ce n'est pas « un code qui n'est pas le leur ». Ce qui reste vrai : ce n'est pas
> la version qui a produit l'article de 2026. Le paragraphe ci-dessous est
> conservé tel quel et doit être réécrit sur la prémisse corrigée.

**117 attributions à `Solidity` / `Y3D*.c` subsistent dans le code, les en-têtes,
la suite de vérification et les decks.** La CORRECTION 2 avait purgé les bilans,
pas le code. Deux clés promettent une réplication qu'elles ne font pas :
`contactMu.<phase>` (règle de paire = minimum) et `contactDamageCoupling = solidity`
(raideur de contact × (1−D)) — **or ni la règle de paire ni le couplage sur la
pénalité ne sont publiés nulle part.** Table de rachat au lot 4 §1.1.

### 5.4 Ce que rockim a et qu'Imperial n'a pas

Treize entrées au lot 4 §3. **Quatre sont des résultats scientifiques**, pas de
l'ingénierie :

* **le garde-fou crack-band** (`MatLaw.cpp:1304-1314`) — Imperial n'a **aucun
  critère d'objectivité publié** ;
* **la séparation mesurée** du levier pénalité et du levier schéma d'insertion :
  à pénalité égale, l'écart de schéma ne vaut que **+1,5 point** ;
* **la provenance du 3000 GPa** retrouvée chez Turon, Dávila, Camanho & Costa
  (2007), règle K = α·E/t avec α ≈ 50 — **et la contradiction interne de Guo mise
  au jour** : il recommande E ≤ p₀ ≤ 10E deux phrases après avoir cité Turon, et
  les auteurs de l'article ont suivi Turon ;
* **une prédiction confirmée** : l'exposant **0,1707** dérivé de la figure 2(b) le
  **2026-08-18**, imprimé **0,17** par IJ206 un an plus tard.

Plus une **suite de non-régression de 98 contrôles** (et non 42), un **bilan
d'énergie fermé** là où le leur est un résidu, et le refus plutôt que le silence
quand une clé serait inerte.

---

## 6. LE PLAN DE TRAVAIL

*Ordonné par ce qui débloque. Chaque étape porte son critère de réussite
mesurable — si le critère n'est pas atteint, l'étape a échoué et il faut le dire.*

### Étape 1 — Fermer le canal d'injection du contact · **1 j** · **BLOQUANT**

Tester `gcBirth = penalty` sur l'impact 3D (jamais fait). Exige
`contact = potential`, que les decks d'impact posent déjà.

> **Critère** : l'injection de la branche normale tombe **sous 1 % de KE₀**
> (contre 6,9 % en adaptatif et 20 % en intrinsèque). Si elle reste au-dessus de
> 3 %, le correctif ne suffit pas et il faut instruire la naissance conditionnelle
> (étape 7) avant d'aller plus loin.

### Étape 2 — Aligner la longueur de référence · **1 j** · **BLOQUANT**

Calculer h comme la **longueur moyenne des arêtes de la facette de joint**, et
non le diamètre inscrit. Opt-in et bannière — la clé change les résultats.

> **Critère** : à `jointPenaltyFactor = 26,32` sous la nouvelle mesure de h,
> l'**ouverture au pic** δ_nE d'un joint doit valoir **2 f_t h_arête / p₀** à
> mieux que 1 %, contrôlé au point matériel par `tools/yan_point.cpp`. Repli
> immédiat si l'étape est reportée : poser **10,74** au deck, en documentant que
> c'est 26,32 corrigé du rapport de longueurs.

### Étape 3 — Re-sourcer les 117 attributions · **1 j** · intégrité

Remplacer les citations rachetables par leur vraie source ; renommer
`contactDamageCoupling = solidity` et `jointDeltaC = solidity` ; porter
l'avertissement de la CORRECTION 2 en tête des deux solveurs.

> **Critère** : `grep -rc "Y3D[a-z]*\.c" src include tools` rend **zéro**, et
> chaque clé survivante cite un article, une page et une équation.

### Étape 4 — Rejouer la réplication St Anne · **2 j** · validation

Avec les étapes 1 et 2, **et sans le modèle de pulvérisation** — Imperial dit
qu'il est superflu sur le calcaire (§2.1). Frottement à **0,6**, pas 0,18.

> **Critères, les leurs** (A952, piston 9,41 m/s, 49,3 J) : part de fissuration
> **2,0-3,5 %**, part de frottement **55-75 %**. **En sachant que leur 65 % est
> un plancher gonflé par leur maillage** : un frottement plus bas chez nous n'est
> pas nécessairement un échec, c'est peut-être un maillage moins anguleux. À
> discuter, pas à corriger d'office.

### Étape 5 — Restreindre les joints à la roche · **2 j** · raffinement

**Aucun filtre par phase n'existe** (vérifié : `jointPhase`, `jointsIn`,
`noJoint`, `jointBodies`, `skipJoint` → zéro occurrence). Il faut l'écrire, dans
les deux solveurs, avec sa clé, sa bannière et son contrôle de non-régression.
*(Estimation révisée le 2026-08-29 : elle était de 0,5 à 1 j, sous l'hypothèse
qu'un filtre existait.)*

> **Critère** : nombre de joints divisé par ≈ 2 ; **vitesse de rebond du taillant
> rapprochée de la mesure** ; champ de fissuration dans la roche **inchangé** à
> la tolérance de la suite de non-régression.

### Étape 6 — Rendre visible la pénalité · **0,5 j** · hygiène

Imprimer la pénalité sous `intrinsic` en **2D** (fait en 3D seulement, commit
`d4be57b`), imprimer la raideur `pj` réelle (min/moy/max), valider `pf > 0`.

> **Critère** : un run 2D intrinsèque imprime sa raideur, sa provenance
> (deck/défaut) et avertit si le deck pose la clé de l'autre schéma ; `pf = 0`
> lève une exception au lieu de produire des NaN.

### Étape 7 — Armer la naissance de contact sur le mode de rupture · **0,5 j**

Conditionner le relevé de naissance à `J.bmode == 2` (cisaillement) **et**
`J.fDeath < 0` (compression), comme l'éq. 18. Les deux existent déjà, mais
`bmode` est déclaré sortie seule.

> **Critère** : plus aucun relevé sur une naissance en traction franche ; le
> résidu du bilan d'énergie ne se dégrade pas.

### Étape 8 — Balayer `gcBirthTau` · **0,5 j**

1e-6 s couvre ≈ 518 pas au dt de St Anne, ~50× le n_total ≈ 10 d'Imperial. Aucune
ligne du dépôt ne justifie cette valeur.

> **Critère** : une courbe résidu d'énergie contre τ, et une valeur choisie sur
> elle plutôt que par défaut.

### Étape 9 — Porter k_t dans le budget de pas de temps du 3D · **0,5 j**

Le 2D compte `potKt_` (`FdemSolver.cpp:3134`) ; **le 3D ne le compte pas**
(`Fdem3dSolver.cpp:2115-2122`). Or Xiang, Munjiza, Latham & Guises (2009) p. 677
avertissent que le calcul des forces tangentielles exige un pas plus petit que le
cas sans frottement — « somewhat alarming », écrivent-ils.

> **Critère** : à jeu de paramètres constant, le `dt` annoncé par le solveur 3D
> baisse dès que `potTangentFactor` monte. S'il ne bouge pas, le correctif n'est
> pas branché.

### Étape 10 — Ajouter le banc analytique du frottement · ~~0,5 j~~ **1-2 j**

> **⚠️ ESTIMATION CORRIGÉE LE 2026-08-29.** rockim n'a **aucun scénario** capable
> d'exprimer ce banc (`percussion | shear | tension | brazilian | shpb`), ni de
> clé de vitesse initiale de corps. Un scénario neuf touche **43 points de
> branchement dans 17 fonctions**. Trois voies chiffrées, et une recommandation,
> dans `chantier_imperial_2026-08-29/A12_banc_frottement.md`.

Rectangle lancé sur un plan, `L = v_i²/(2µg)`, configuration publiée complète
(lot 2c §3ter). `verify_suite.py` n'a **aucun** contrôle du chemin tangentiel de
contact — ses contrôles de frottement portent tous sur le joint.

> **Critère** : distance d'arrêt simulée à moins de 1 % de `v²/(2µg)` au pas fin,
> et l'écart doit **croître** au pas grossier — c'est le comportement que les
> auteurs rapportent, et le reproduire valide le chemin en plus de la valeur.

### Ce qu'il ne faut PAS faire

* **Ne pas implémenter de maillage adaptatif.** Imperial n'en a pas et le range
  en perspective. Ce serait **dépasser** l'état de l'art, pas le rejoindre.
* **Ne pas monter la pénalité à 50 E.** C'était ma recommandation du lot 3, elle
  était fausse d'un facteur ~5 (lot 4 §4).
* **Ne pas passer `contactMu` à 0,18 sur St Anne.** Imperial y met **0,6**.
* **Ne pas attendre grand-chose de la détection par événement** : le dépôt a
  mesuré que le poste dominant est ailleurs (10-15 % du mur, pas un facteur 7).
* **Ne pas chercher à aligner `k_t` sur Imperial** : aucun nombre n'est publié
  (huit sources). Le rapport `k_t/k_n = 2/7` des decks vient du code
  non-Imperial et doit être requalifié, pas corrigé.
* **Ne pas activer la pulvérisation sur un cas St Anne.** Si elle change les
  résultats, c'est un artefact.

---

## 7. CE QUE CETTE SESSION A CORRIGÉ

Du dépôt :

* le frottement **0,18 ne vaut que pour le granite** ; St Anne est à **0,6**, ce
  que le deck faisait déjà — la fiche 2026 §3 concluait à tort que le deck était
  3,3 fois trop élevé ;
* la pénalité de **3000 GPa = 52,6 E est réelle** : la CORRECTION 1 du 29/08 qui
  la réfutait était elle-même fausse, ayant lu un `mat.txt` non-Imperial ;
* la « moitié manquante » du couplage sur la pénalité de contact **n'existe pas
  chez Imperial** : le §5 de la fiche 2026 est sans objet.

De moi-même :

* la coquille du DIF **avait été trouvée par le dépôt le 2026-08-18**, avant moi
  et mieux — crédit porté ;
* le **3000 GPa avait déjà été analysé** par le dépôt, jusqu'à sa provenance ;
* ma recommandation **R2 était fausse** (facteur ~5), **R4 surévaluée** (10-15 %
  et non ×7), et mon **facteur 14** de coût de maillage mal attribué ;
* le papier ARMA d'orchyd.eu **n'était pas celui que j'annonçais** au lot 1.

---

## 8. CE QUI RESTE À FAIRE DE LA MISSION

1. **Les sept audits de code non faits et les onze contre-audits** (lot 4 §0).
   En particulier : le couplage D → frottement de rockim est-il continu ou en
   tout-ou-rien ? existe-t-il déjà un filtre de joints par phase ?
2. **Xiang, Munjiza & Latham (2009), *IJNME* 79(8), 946-978** — dernière pièce
   susceptible de donner **k_t et η**. Si elle ne les donne pas, le `2/7` de
   rockim devient un choix documenté en propre, ce qui est tenable.
3. **Guo, Xiang, Latham & Izzuddin (2016), *Eng. Fract. Mech.* 151, 70-91** — la
   seule étude de sensibilité au maillage du modèle, hors dossier.
4. **Les livrables WP6 d'ORCHYD** : vous êtes chez le coordinateur, et un
   livrable contient les decks qu'un article comprime.
