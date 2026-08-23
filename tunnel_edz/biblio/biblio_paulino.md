« VOLET 3 — École Paulino / Celes : le maillage comme règle de déclenchement cachée »

Convention de marquage : **[V]** = page consultée dans ce run ; **[M]** = connaissance propre non revérifiée ici ; **[D]** = dérivation géométrique élémentaire faite ici (vérifiable au crayon).

---

## 1. Ce que cette école dit vraiment : le critère n'est qu'un tiers du problème

Chez Paulino/Celes, l'insertion extrinsèque repose sur un **triptyque indissociable**, et la littérature du groupe montre que corriger le seul critère de contrainte ne suffit jamais :

| Pilier | Question | Référence pivot |
|---|---|---|
| **Critère** | *quand* ? seuil, contrainte évaluée où, garde temporelle | Camacho–Ortiz 1996 ; Papoulia–Sam–Vavasis 2003/2005 |
| **Topologie** | *où* exactement ? quels nœuds dupliquer, comment garder le modèle cohérent | TopS : Celes–Paulino–Espinha 2005 ; Paulino–Celes–Espinha–Zhang 2008 **[V]** |
| **Maillage** | *quelles directions sont offertes* ? quel surcoût de longueur/ténacité | Papoulia–Vavasis–Ganguly 2006 ; Rimoli–Rojas 2015 **[V]** ; Spring/Leon–Paulino 2014 **[V]** |

Le diagnostic « nuage diffus » de rockim est très exactement le symptôme d'un pilier 1 + pilier 2 non traités, **amplifié** par le pilier 3.

---

## 2. Règles d'activation exactes

### 2.1 Le seuil (forme canonique, reprise telle quelle par toute l'école)

Critère de traction effective de Camacho–Ortiz, confirmé mot pour mot sur une source récente **[V]** :

```
σ_eff = sqrt( σ_n²  +  σ_t² / β² )        avec σ_n < 0  →  σ_n := 0
insertion sur la facette  ⟺  σ_eff ≥ σ_c
```

où β = rapport résistance cisaillement/traction. **[V]** (arXiv 2511.14323, éq. 9, attribution explicite à Camacho & Ortiz).

Deux écarts importants avec l'enveloppe de Yan 2023 utilisée dans rockim :
- **pas de branche frottante en compression** : σ_n comprimé est *écrêté à zéro*, il ne « paye » pas le cisaillement. Chez Yan, l'enveloppe `|τ| ≥ c − σ_n tan φ` **rend le seuil plus facile à atteindre quand la compression augmente** dès que tan φ dépasse la pente effective du champ. Dans une zone plastique confinée de tunnel profond ou sous un insert de percussion, σ_n est fortement compressif partout : l'enveloppe frottante y **abaisse** le seuil sur un domaine étendu, alors que le critère Camacho–Ortiz écrêté n'y déclenche presque rien. C'est un candidat direct au « tapis ». **[V]** pour la forme Camacho–Ortiz, **[M]** pour la lecture comparative.
- **un seul scalaire**, pas de test séparé mode I / mode II : l'ellipse σ_eff mélange les modes et ne peut pas déclencher deux populations de fissures indépendantes.

### 2.2 Où la contrainte est évaluée — le point le plus sous-estimé

Deux familles, et elles ne donnent **pas** le même nombre d'insertions :

- **(A) Contraintes élémentaires moyennées.** σ pris aux points de Gauss des deux éléments voisins, moyenné, projeté sur la normale de la facette. C'est ce que fait la lignée Zhang–Paulino–Celes **[M]**, ce que confirme la source récente (« computed from **B** and **D** of adjacent continuum elements ; no nodal-force-based computation ») **[V]**, et c'est aussi la formulation de Yan 2023 dans rockim.
- **(B) Partition nodale des forces internes (Camacho–Ortiz 1996).** La traction de la facette est reconstruite à partir des forces internes nodales effectivement transmises à travers elle. **[M]**

La différence est fondamentale et n'est pas cosmétique : avec (B), le cohésif inséré **reprend exactement la force qu'il vient de retirer** au continuum → aucun saut de force à l'instant d'insertion. Avec (A), la contrainte moyennée n'est pas la force réellement transmise : chaque insertion émet une **petite onde parasite**. À 28 000 insertions, ces ondes forment un bruit de fond qui re-déclenche le critère ailleurs — mécanisme d'avalanche auto-entretenue produisant précisément un nuage de fissures courtes. **[M] — hypothèse mécaniste, à tester (voir §6, test T2).**

### 2.3 La garde temporelle : continuité en temps (Papoulia–Vavasis)

C'est *la* règle de déclenchement propre à cette école, et elle est presque toujours omise dans les implémentations maison.

- Papoulia, Sam & Vavasis, *Time continuity in cohesive finite element modeling*, IJNME 58 (2003) : introduction de la notion de **continuité en temps** pour les modèles « initially rigid » (= extrinsèques) : l'interface est inactive jusqu'à ce que la traction atteigne un niveau critique. **[V]** (page Wiley + résumé)
- Sam, Papoulia & Vavasis, *Obtaining initially rigid cohesive finite element models that are temporally convergent*, Eng. Fract. Mech. 72 (2005) : cadre fondé sur une **condition de continuité en temps** produisant des modèles **temporellement convergents** en dynamique explicite. **[V]** (fiche UWaterloo)

Conséquence opérationnelle : un balayage naïf « à chaque pas, si σ_eff ≥ σ_c alors insère avec T(0) = σ_c » **n'est pas temporellement convergent**. Entre deux pas la contrainte dépasse le seuil d'un Δσ = σ_eff − σ_c > 0 ; on force ensuite T = σ_c : le saut −Δσ est un dépôt d'énergie parasite, et surtout **le nombre de facettes insérées dépend de Δt**. Sur un champ lisse au seuil (le cas du tunnel), tout le plateau franchit le seuil dans la même fenêtre Δt → insertion « en tapis » dont l'étendue est un artefact de pas de temps. La correction est simple : **initialiser la loi cohésive à la traction courante réellement portée** (donc σ_max local = σ_eff au moment de l'insertion, pas σ_c), ce qui rend T(t) continue.

### 2.4 La garde de stabilité à l'insertion (résultat récent, très pertinent)

Pour une loi cohésive extrinsèque linéaire, la raideur tangente d'un élément fraîchement inséré s'écrit **[V]** :

```
k⁺(d) = (σ_c / δ_c) · (1 − d) / d       →  diverge quand d → 0
```

Un cohésif neuf (d ≈ 0) est donc **quasi rigide**, et il fait s'effondrer le pas de temps critique. L'analyse de rayon spectral donne : à Δt = 0,1·Δt_c,b, **l'instabilité apparaît pour d < d̃ ≈ 0,02** **[V]**. La borne de Gershgorin doit être calculée **en incluant toutes les facettes cohésives potentielles**, pas seulement celles déjà insérées **[V]**. Le papier observe une **croissance non physique du nombre de fragments** qu'il attribue à cette instabilité numérique et non au modèle **[V]**.

> À retenir pour rockim : **un excès de fragments courts peut être un symptôme d'instabilité de pas de temps à l'insertion, pas un symptôme de critère.** C'est le premier test à faire, il est gratuit (§6, T1).

---

## 3. TopS : *où* et *comment* insérer — la duplication de nœud conditionnelle

Réf. **[V]** : Paulino, Celes (Filho), Espinha & Zhang, *A general topology-based framework for adaptive insertion of cohesive elements in finite element meshes*, **Engineering with Computers** (2008), Springer — « algorithme fondé sur la topologie qui **classifie les facettes fracturées** et effectue les **changements topologiques** nécessaires pour mettre à jour le modèle », 2D et 3D, éléments d'ordre quelconque. Base : la structure **TopS**, représentation réduite (éléments + nœuds seulement, arêtes/facettes comme **entités implicites** reconstruites à la volée par adjacence, accès O(1)) **[M]** ; extension parallèle **ParTopS** (Espinha, Celes, Rodriguez, Paulino, *Engineering with Computers* 2009) **[M]**.

La règle d'insertion, en trois temps **[M]** :

1. **Marquer** les facettes fracturées par le critère (§2).
2. **Classifier chaque nœud touché** : on retire les facettes fracturées de l'étoile du nœud et on compte les **composantes connexes** d'éléments restantes.
3. **Dupliquer le nœud en autant de copies qu'il y a de composantes connexes**, puis créer les éléments cohésifs.

**C'est la garde anti-sur-insertion la plus efficace de cette école, et elle est purement topologique.** Une facette insérée **isolée** ne déconnecte l'étoile d'aucun de ses deux nœuds (on peut toujours en faire le tour par l'autre côté) : **aucun nœud n'est dupliqué**, le cohésif partage ses deux nœuds avec le continuum et, en 2D linéaire, **son ouverture est identiquement nulle**. Il est **inerte** tant qu'il ne se raccorde pas à une surface libre ou à une autre fissure. La topologie fait donc gratuitement le tri entre « tapis d'insertions parasites » (inerte, sans dissipation, sans ouverture) et « chaîne connexe » (vraie macro-fissure qui découpe un bloc).

> **Diagnostic n°1 pour rockim** : si rockim **duplique les nœuds inconditionnellement** à chaque insertion (ce que fait le schéma FDEM natif, où les nœuds sont pré-dupliqués), alors **chaque** facette insérée devient immédiatement une micro-fissure ouvrante et dissipante → nuage diffus, énergie pompée partout, plus de localisation possible. Les 28 000 joints de rockim contre ~10 000 de MultiFracS s'expliquent en partie par là : ce ne sont pas 28 000 fissures, ce sont peut-être 8 000 fissures + 20 000 facettes qui, dans un cadre TopS, seraient restées inertes. **[M] — hypothèse forte, test T3 en §6.**

Un second usage de la topologie : la classification par composantes connexes gère **automatiquement** les jonctions, bifurcations et fragments détachés (un nœud entouré de k clusters devient k nœuds), ce qui est la condition pour que des **blocs** émergent au lieu d'un champ de micro-coupures.

---

## 4. Zhang, Paulino & Celes 2007 — microbranchement et biais de maillage

Réf. **[V]** : *Extrinsic cohesive modelling of dynamic fracture and microbranching instability in brittle materials*, **IJNME 72(8) 893–923 (2007)**. Contenu confirmé **[V]** : évaluation critique **intrinsèque vs extrinsèque** concluant à la **nécessité de l'extrinsèque** ; structure de données topologique nouvelle pour manipuler rapidement et robustement le maillage évolutif lors de l'insertion adaptative ; vitesse limite de fissure et **résistance croissante avec la vitesse**.

Éléments complémentaires **[M]** (non revérifiés ici, à confirmer sur le PDF si vous y tenez) :
- Motif de l'extrinsèque : la loi intrinsèque introduit une **compliance artificielle** qui abaisse la vitesse d'onde effective et fausse l'angle et l'espacement des microbranches → l'extrinsèque est imposé par la physique du branchement, pas par le coût.
- Résolution requise : plusieurs éléments (ordre 3–5) dans la longueur de zone cohésive `ℓ_cz ≈ E·G_c/σ_c²` ; sinon le critère est atteint **simultanément** sur toute une bande.
- L'angle de branchement et l'espacement des branches sont **sensibles à l'orientation du maillage** — le motif est co-déterminé par le maillage.

---

## 5. Le maillage : preuves chiffrées du biais

### 5.1 Le fait géométrique de base **[D]**

Une fissure contrainte de suivre les arêtes, entre deux directions disponibles α et β, pour une direction physique θ :

```
L_discret / L_vrai = [ sin(β−θ) + sin(θ−α) ] / sin(β−α)
```

- Maillage **quad pur** (directions {0°, 90°}) et θ = 45° : ratio = **√2 = 1,414 → +41,4 %**, et ce ratio **ne converge pas** quand h → 0 (paradoxe de l'escalier).
- Maillage **4k** (quad + les deux diagonales : directions {0°, 45°, 90°, 135°}, **4 directions**) et θ = 22,5° : ratio = 2·sin22,5°/sin45° = **1,0824 → +8,2 %**.
- Maillage **triangles rectangles structurés** (une seule diagonale : **3 directions** {0°, 45°, 90°}) : même borne de +8,2 % *d'un côté*, mais l'écart est **asymétrique** → biais directionnel net vers la diagonale unique.

### 5.2 Rimoli & Rojas 2015 — la mesure **[V]**

*Meshing strategies for the alleviation of mesh-induced effects in cohesive element models*, **Int. J. Fracture** (2015) / arXiv:1302.1161. Ils formalisent l'**anisotropie induite par le maillage** et la **ténacité induite par le maillage**, quantifiées par des **diagrammes polaires du path deviation ratio** **[V]** :

```
η = L_g / L_e   (chemin le plus court le long des arêtes / distance euclidienne)
ε = η − 1
```

Chiffres relevés sur le texte **[V]** :

| Maillage | ε mesuré | commentaire |
|---|---|---|
| **4k** (λ ≈ 1/200) | **ε_max ≈ 0,08** dans les directions intermédiaires | colle exactement à la dérivation +8,2 % **[D]** |
| **4k + NP + ES** (perturbation nodale + edge-swap, λ ≈ 1/250) | **moyenne(ε) ≈ 0,037** | |
| **Maillage aléatoire** (Delaunay, λ ≈ 1/250) | **moyenne(ε) ≈ 0,043** | |
| **K-means** (λ ≈ 1/250) | saturation ≈ **0,04** | comme 4k+NP |
| **Conjugate-directions** (λ ≈ 1/250) | **moyenne(ε) ≈ 0,018**, η ≈ **1,015** à spn = 512 | meilleur du lot |

Point capital **[V]** : les maillages 4k présentent une **« rugosité intrinsèque » qui empêche de descendre sous 3,6–4 % d'erreur, quel que soit le raffinement**. Autrement dit : **raffiner ne guérit pas le biais de maillage**, seule la *nature* du maillage le fait. Idem pour le Delaunay aléatoire, qui **sature** à ~4,3 %.

**Construction du maillage à directions conjuguées** **[V]** : (1) mailler par **K-means** — clusteriser spn·n points aléatoires en n nœuds, puis trianguler par Delaunay ; (2) **subdivision barycentrique** de chaque simplexe (barycentres + connexions selon les relations de faces) ; (3) l'enrichissement vient de ce que, **pour un triangle équilatéral, la subdivision barycentrique fournit une direction orthogonale pour chaque direction portée par les arêtes du triangle initial** — d'où le nom « directions conjuguées ».

Interprétation mécanique de η, et c'est là qu'est le lien avec le **motif** : la fissure paye G_c par unité de **longueur d'arête** parcourue, donc la ténacité apparente vaut `G_eff(θ) ≈ η(θ)·G_c`. Un diagramme polaire non circulaire est une **ténacité anisotrope numérique** : la fissure est **attirée vers les minima de η**, c'est-à-dire vers les directions du maillage. Corollaire pour votre cas : sur un champ isotrope et lisse (tunnel profond), si η(θ) est presque plat mais bruité — cas du Delaunay aléatoire, ε ≈ 4,3 % sans direction privilégiée — **il n'existe aucune direction « moins chère » pour organiser une macro-fissure** : le milieu casse partout au même prix. Le nuage diffus est la réponse rationnelle d'un solveur à un maillage sans direction préférentielle *et* sans autre brisure de symétrie (ni bruit matériau type Weibull, ni contrainte topologique).

### 5.3 Papoulia–Vavasis–Ganguly — les maillages pinwheel

- Motivation, telle que formulée par le groupe **[V]** : dans les approches classiques, « les fissures sont forcées de se propager le long des frontières d'éléments, suivant des chemins qui **requièrent plus d'énergie par unité d'extension de fissure** que dans le continuum d'origine, ce qui conduit à des solutions erronées ».
- Ganguly, Vavasis & Papoulia, *An algorithm for two-dimensional mesh generation based on the pinwheel tiling*, **SIAM J. Sci. Comput.** (algorithme **PINW**) **[V]** : extension des pavages pinwheel de **Radin et Conway**, produisant des **triangles de rapport d'aspect borné**, « utile en modélisation par éléments d'interface cohésifs **quand le chemin de fissure est un résultat de la simulation** ».
- Papoulia, Vavasis & Ganguly (IJNME 2006), *Spatial convergence of crack **nucleation** using a cohesive finite-element model on a pinwheel-based mesh* **[V pour l'existence/titre]**. Le mot **nucleation** est décisif : ce qu'ils démontrent converge, ce n'est pas le trajet d'une fissure préexistante, **c'est l'endroit où la fissure naît** — exactement votre problème (28 000 nucléations dispersées vs des macro-fissures).
- Mécanisme **[M]** : le triangle pinwheel (1, 2, √5) se subdivise en 5 copies semblables tournées de ±arctan(1/2) ≈ 26,57°, **multiple irrationnel de π** → sous raffinement, les orientations d'arêtes deviennent **denses et équidistribuées** sur [0, π) (équidistribution prouvée par Radin & Sadun). Le pavage possède alors la propriété quasi-isopérimétrique : **η → 1 quand h → 0**. C'est le seul des maillages cités dont le biais **s'annule** par raffinement, au lieu de saturer.

### 5.4 Spring, Leon & Paulino — maillages polygonaux + découpe adaptative

- Spring, Leon & Paulino, *Unstructured polygonal meshes with adaptive refinement for the numerical simulation of dynamic cohesive fracture*, **Int. J. Fract.** (2014) **[V]** : les **maillages polygonaux à germes aléatoires fournissent une discrétisation isotrope qui ne biaise pas les motifs de fissuration** **[V]** ; maillage par **tessellations de Voronoï centroïdales** (CVT) pour la qualité d'élément **[V]** ; raffinement adaptatif pour mieux capter les motifs **[V]**.
- Leon, Spring & Paulino, *Reduction in mesh bias for dynamic fracture using **adaptive splitting** of polygonal finite elements*, **IJNME** (2014) **[V]** : constat clé **[V]** — « la discrétisation polygonale **limite les directions de fissure possibles à chaque nœud** », problème traité par **découpe adaptative d'élément, qui augmente le nombre de directions potentielles à chaque pointe de fissure** ; des **études géométriques** démontrent le bénéfice du polygonal + découpe sur les discrétisations structurées et non structurées.
- Complément **[M]** : un sommet de Voronoï générique est **trivalent** → seulement **3 directions disponibles au nœud** (contre 6–8 en Delaunay), donc une erreur de direction locale pouvant atteindre ~60°, mais **sans anisotropie globale** puisque ces 3 directions sont aléatoires. La découpe (split) crée à la volée une arête **dans la direction voulue** : c'est un **enrichissement directionnel à la demande**, autrement dit **la seule stratégie de l'école qui rende la direction de fissure indépendante du maillage sans avoir à raffiner**.
- Adaptation de maillage (Park, Paulino, Celes, Espinha, IJNME 2012, raffinement/déraffinement pour CZM dynamique) **[M]** : effet secondaire précieux, c'est un **throttle spatial de fait** — loin de la zone de process le maillage reste grossier, donc peu de facettes candidates, donc peu d'insertions possibles, sans qu'aucun compteur ne soit nécessaire.

---

## 6. Synthèse opérationnelle pour rockim : gardes et remèdes, par ordre de coût croissant

**Ce que cette école ne fait PAS** : elle n'ordonne pas les insertions, ne limite pas leur nombre par pas, n'introduit pas de compteur ni de délai. **Aucun throttling algorithmique.** Ses quatre gardes sont ailleurs — temporelle, topologique, géométrique, résolutive — et elles sont toutes absentes par défaut d'une implémentation « Yan 2023 » naïve.

**T1 — Garde de stabilité (gratuit, à faire en premier).** Relancer le tunnel à Δt/4. Si le nombre de joints cassés chute nettement, une partie du nuage est de l'**instabilité d'insertion** (k⁺ → ∞ pour d → 0, instable pour d < 0,02 à Δt = 0,1 Δt_c) **[V]**, pas de la physique. Remède : calculer Δt_c par borne de Gershgorin **en incluant toutes les facettes cohésives potentielles** **[V]**, ou plafonner k⁺ (d_min plancher ≈ 0,02).

**T2 — Continuité en temps (peu coûteux, fort impact attendu).** Initialiser le cohésif inséré avec **la traction courante** (σ_max local := σ_eff à l'insertion) plutôt qu'avec σ_c, conformément à Papoulia–Sam–Vavasis **[V]**. Test de convergence temporelle : halver Δt ; le nombre d'insertions doit être **stable**. S'il ne l'est pas, le « tapis » est un artefact de Δt.

**T3 — Duplication nodale conditionnelle (remède structurel).** Adopter la règle TopS : ne dupliquer un nœud **que si** son étoile se déconnecte, en comptant les composantes connexes après retrait des facettes fracturées **[V]** pour le cadre, **[M]** pour l'algorithme détaillé. Effet attendu : les facettes isolées deviennent **inertes** (ouverture nulle, dissipation nulle) et seules les **chaînes connexes** ouvrent → passage mécanique du nuage aux blocs. Métrique de contrôle à ajouter dès maintenant, même sans changer le code : **fraction des joints insérés dont l'ouverture dépasse 10 % de δ_c**. Si elle est très faible, vos 28 000 joints sont majoritairement du bruit comptable ; si elle est élevée, le problème est bien énergétique.

**T4 — Audit géométrique du maillage (gratuit, sans simulation).** Reproduire le protocole Rimoli–Rojas sur votre maillage de tunnel : plus court chemin le long des arêtes entre paires de nœuds, **diagramme polaire de η(θ)**, moyenne et max de ε **[V]**. Cibles : moyenne(ε) < 2 % (niveau conjugate-directions **[V]**) ; si vous mesurez ~4 % avec un polaire quasi circulaire, vous avez la **signature exacte du milieu sans direction privilégiée** décrite en §5.2, et le nuage est prévisible avant même de lancer le calcul.

**T5 — Changer la nature du maillage, pas sa finesse.** Par ordre décroissant de gain mesuré : **conjugate-directions** (K-means + subdivision barycentrique, ε ≈ 1,8 %) **[V]** > **NP+ES** ou aléatoire (≈ 3,7–4,3 %, **saturants**) **[V]** ; **pinwheel** si vous voulez un biais **qui s'annule** par raffinement plutôt que qui sature **[M]** ; **polygonal CVT + splitting adaptatif** si vous acceptez de toucher au mailleur en cours de calcul, seule option qui donne la direction physique exacte à la pointe **[V]**. Retenir de **[V]** : sur 4k, raffiner ne fera **jamais** descendre sous ~3,6–4 %.

**T6 — Résolution de la zone cohésive.** Vérifier `ℓ_cz ≈ E·G_c/σ_c²` et exiger 3–5 éléments dedans **[M]**. Si h > ℓ_cz, toute la bande plastique franchit le seuil dans le même pas : le tapis est alors une **conséquence arithmétique** du sous-maillage, indépendamment de tout le reste.

**T7 — Enveloppe de rupture.** Comparer, sur votre champ de tunnel, la carte des facettes candidates sous l'enveloppe frottante de Yan (`|τ| ≥ c − σ_n tan φ`) et sous le critère écrêté de Camacho–Ortiz (`σ_n < 0 → 0`) **[V]**. Si la première illumine une aire massivement plus grande en compression, la branche frottante est un contributeur de premier ordre au tapis — ce qui rejoint le point de votre volet acquis sur la branche compression-frottement que les codes modernes suppriment.

---

## 7. Ce que ce volet démontre en une phrase

Sur un champ lisse au seuil, l'école Paulino/Celes ne compte pas sur le critère pour localiser : elle laisse le critère insérer largement et **fait porter la sélection du motif par la topologie (duplication nodale conditionnelle) et par la géométrie du maillage (η(θ))** ; Papoulia–Vavasis ajoutent la seule garde qui rende le résultat indépendant de la discrétisation temporelle. Un solveur qui duplique les nœuds inconditionnellement, initialise à σ_c et maille en Delaunay aléatoire cumule les trois défauts — et produit exactement un nuage de fissures courtes.

---

**Sources consultées dans ce run [V]**
- [Paulino, Celes, Espinha, Zhang — A general topology-based framework for adaptive insertion of cohesive elements (Engineering with Computers, 2008)](https://link.springer.com/article/10.1007/s00366-007-0069-7)
- [Rimoli & Rojas — Meshing strategies for the alleviation of mesh-induced effects in cohesive element models (Int. J. Fract., 2015)](https://link.springer.com/article/10.1007/s10704-015-0013-6) · [texte intégral arXiv:1302.1161](https://arxiv.org/abs/1302.1161)
- [Spring, Leon & Paulino — Unstructured polygonal meshes with adaptive refinement for dynamic cohesive fracture (Int. J. Fract., 2014)](https://link.springer.com/article/10.1007/s10704-014-9961-5)
- [Leon, Spring & Paulino — Reduction in mesh bias for dynamic fracture using adaptive splitting of polygonal finite elements (IJNME, 2014)](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.4744)
- [Ganguly & Vavasis — An algorithm for two-dimensional mesh generation based on the pinwheel tiling (Semantic Scholar)](https://www.semanticscholar.org/paper/An-Algorithm-for-Two-Dimensional-Mesh-Generation-on-Ganguly-Vavasis/76a1e6273482de9222c9943bed9ed1c9c898dd13)
- [Papoulia, Sam & Vavasis — Time continuity in cohesive finite element modeling (IJNME, 2003)](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.778)
- [Sam, Papoulia & Vavasis — Obtaining initially rigid cohesive finite element models that are temporally convergent](https://uwaterloo.ca/scholar/vavasis/publications/obtaining-initially-rigid-cohesive-finite-element-models-are-temporally)
- [Zhang, Paulino & Celes — Extrinsic cohesive modelling of dynamic fracture and microbranching instability (IJNME 72, 2007)](https://www.academia.edu/736760/Extrinsic_cohesive_modelling_of_dynamic_fracture_and_microbranching_instability_in_brittle_materials)
- [Stability of Extrinsic Cohesive-Zone Model with Penalty-Based Contact in Explicit Dynamic Fragmentation Simulations (arXiv:2511.14323)](https://arxiv.org/html/2511.14323v1)
- [Advances in Cohesive Zone Modeling of Dynamic Fracture (chapitre Springer)](https://link.springer.com/chapter/10.1007/978-1-4419-0446-1_12)

## Lacunes

Budget de 8 appels épuisé (4 WebSearch groupés + 1 WebFetch ar5iv Rimoli-Rojas + 2 WebSearch + 1 WebFetch arXiv 2511.14323). Lacunes assumées, à combler dans un volet ultérieur si nécessaire :

1. Zhang-Paulino-Celes 2007 : seuls le titre, la revue, la pagination et les conclusions générales sont [VERIFIE] (résumé de recherche). La formule exacte d'évaluation de la traction de facette (moyenne des contraintes élémentaires vs partition nodale), le nombre d'éléments exigé dans la zone cohésive, et la fraction de Δt_c utilisée restent [MEMOIRE] — non revérifiés. Le PDF n'a pas été ouvert (interdit par le budget).

2. TopS (Celes-Paulino-Espinha IJNME 2005) : la référence bibliographique et le principe « entités implicites / représentation réduite » sont [MEMOIRE]. Seul le papier-cadre de 2008 (Engineering with Computers) a été confirmé, et uniquement au niveau du résumé (« classifie les facettes fracturées et effectue les changements topologiques »). L'algorithme précis de duplication nodale par composantes connexes — pilier de mon diagnostic T3 — est reconstruit de mémoire et devrait être vérifié sur le PDF de 2008 avant d'être cité dans le manuscrit.

3. Le point le plus important à vérifier : ma conclusion « une facette insérée isolée reste inerte car aucun nœud n'est dupliqué » est une déduction logique [MEMOIRE], pas une citation. Si elle est fausse, tout le remède T3 tombe.

4. Papoulia-Vavasis-Ganguly IJNME 2006 : titre et objet ([VERIFIE] via résumé indirect), mais AUCUN chiffre de convergence de nucléation obtenu. L'équidistribution des orientations pinwheel (Radin & Sadun) et la géométrie du triangle (1,2,√5 ; rotation arctan(1/2)) sont [MEMOIRE].

5. Rimoli-Rojas : les chiffres sont solides [VERIFIE via ar5iv] mais l'extraction indique explicitement qu'aucun tableau comparatif exhaustif n'existe et qu'aucun comptage de directions par nœud n'est donné. Le paramètre λ (rapport taille d'élément / domaine) et la définition de « spn » n'ont pas été confirmés. Le maillage pinwheel ne figure pas dans les chiffres extraits — je ne peux donc pas comparer pinwheel et conjugate-directions sur la même métrique.

6. Spring/Leon-Paulino : aucun chiffre. Les « geometric studies » sont mentionnées comme existantes [VERIFIE] mais leurs valeurs (déviation angulaire moyenne par type de maillage, gain du splitting) n'ont pas été récupérées. La trivalence des sommets de Voronoï (3 directions par nœud) est [MEMOIRE].

7. Paulino-Park-Celes-Espinha 2010 (nodal perturbation + edge-swap) : n'a pas été recherché directement ; son existence est attestée indirectement par les labels « 4k avec NP+ES » du tableau de Rimoli-Rojas [VERIFIE]. L'amplitude de perturbation recommandée n'est pas connue.

8. Park-Paulino-Celes-Espinha 2012 (raffinement/déraffinement adaptatif) : cité de mémoire, non vérifié.

9. Non traité faute de budget : comment MultiFracS (schéma intrinsèque, Wang et al. 2024) produit ses ~10 000 joints et ses blocs — la comparaison quantitative rockim/MultiFracS reste asymétrique.