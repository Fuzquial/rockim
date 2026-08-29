# LOT 3 — Schémas d'insertion et maillage chez Imperial : mémo comparatif
# et recommandation pour rockim

*Fiche du 2026-08-29. Livrable 3 de [MISSION_etat_de_l_art_2026-08-29.md](../MISSION_etat_de_l_art_2026-08-29.md) §4.
Sources lues de première main dans cette session. Chaque affirmation porte sa
page. Dépouillement thématique conduit par agents, puis **vérifié
adversarialement** ; les points capitaux ont en outre été recontrôlés à la main
sur le texte source par le rédacteur.*

**Codes sources** (mêmes qu'aux lots 2a-2c) :
**D1** = Guo, Xiang, Latham & Izzuddin, manuscrit UCL (formulation) ·
**D2** = Yang et al., *IJRMMS* **191** (2025) 106125 (St Anne, Rhune) ·
**D3** = Yang et al., *JRMGE* (2025) (Kuru Grey) ·
**D4** = Naderi et al., *JRMGE* (2025) (réseau de neurones) ·
**D6** = Yang et al., *IJRMMS* **206** (2026) 106660 (pulvérisation).

---

## 1. RÉPONSE COURTE AUX QUATRE QUESTIONS DU BRIEF

| question du brief §4 | réponse |
|---|---|
| que fait l'insertion **intrinsèque** ? | joints partout dès le maillage, nœuds dédoublés avant chargement, pénalité p₀ en contrainte, **souplesse artificielle assumée et non corrigée** |
| que fait l'insertion **adaptative / extrinsèque** ? | **rien : elle n'existe pas chez eux.** Aucune des cinq sources n'en décrit, n'en cite ni n'en utilise |
| le **maillage adaptatif** de la thèse de Guo | **il n'existe pas non plus.** Le maillage est fixe, et D1 range le raffinement local en **perspective de recherche** |
| le **périmètre des joints** | **la roche seule.** L'outil est non fracturable — dit explicitement une seule fois, dans D4 ; déduit ailleurs de l'absence des paramètres de joint dans les tables |

**La prémisse du brief sur le maillage adaptatif était fausse.** Je l'avais
signalé comme suspicion au lot 1 ; c'est maintenant établi.

---

## 2. L'INSERTION EST INTRINSÈQUE — et personne ne le dit avec ce mot

**[LU]** D1 p. 8, la phrase de discrétisation :

> « The computational model works in a domain discretised by **4-node tetrahedral
> elements and special 6-node joint elements**. »

**[LU]** D1 p. 9, l'acte d'insertion :

> « To achieve the objective of separating tetrahedral elements, **6-node joint
> elements are inserted between 4-node tetrahedral elements** (Figure 2) so the
> failure criterion can be applied to the joint elements. »

**[LU]** D1 p. 14, la portée :

> « **Because the entire domain is discretised by tetrahedral elements and joint
> elements**, the contact means contact between tetrahedral elements. »

**[LU]** D1 p. 16, l'état avant fracture :

> « **before fracture formation the tetrahedral elements are purely connected by
> joint elements**, which contribute to the nodal forces f_joint »

**[INFÉRÉ — et c'est la preuve la plus solide]** D1 p. 15 :

> « the tetrahedral elements in group G1 **do not really share node N1** because
> the space discretisation scheme inserts joint elements between tetrahedral
> elements […] the definition of 'sharing node' […] is that the nodes have the
> **same initial coordinates before loading starts**. »

*Raisonnement* : les nœuds sont **déjà dédoublés**, et ce dédoublement est
**antérieur au chargement**. Un schéma extrinsèque, qui crée les joints en cours
de calcul, ne pourrait pas définir la connectivité par des coordonnées initiales
identiques : il n'y aurait rien à dédoubler avant que le critère ne se déclenche.
**L'insertion est donc intrinsèque, à 100 % du domaine fracturable, dès le
maillage.**

**[ABSENT — recherche exhaustive sur cinq sources]** Les mots `extrinsic`,
`adaptive insertion`, `on-the-fly`, `activated joint`, `insertion criterion` :
**zéro occurrence**. Imperial ne discute même pas l'alternative. Le débat
intrinsèque/extrinsèque qui occupe la littérature Camacho-Ortiz-Pandolfi
(cf. [`biblio_volet2_fondations.md`](biblio_volet2_fondations.md)) **n'a pas
lieu chez eux** : la question est tranchée par défaut depuis Munjiza.

---

## 3. LE PÉRIMÈTRE DES JOINTS : la roche seule

### 3.1 La seule phrase explicite de tout le corpus

**[LU]** **D4, p. 6870** — et c'est un modèle **2D**, à éléments triangulaires :

> « **Joint elements are embedded among the triangular elements, except in the
> insert, which is treated as a fracture-free component.** »

Une seule phrase, dans un article sur un réseau de neurones, en 2D. C'est tout ce
que la littérature dit franchement sur la question.

### 3.2 Ailleurs, il faut déduire — et la déduction est solide

**[INFÉRÉ]** D2 p. 6 et D3 p. 5 : les tables de propriétés des métaux ne portent
**que** densité, module d'Young, coefficient de Poisson et frottement glissant.
Il leur manque **exactement les cinq paramètres qui définissent un joint** :
G_I, G_II, résistance à la traction, cohésion, coefficient de frottement interne.

| | D2 Table 4 (roches) | D2 Table 5 (acier, carbure) | D3 Table 2 (acier) |
|---|---|---|---|
| G_I, G_II | oui | **absents** | **absents** |
| f_t, c, tan φ | oui | **absents** | **absents** |

*Raisonnement* : sans résistance ni énergie de rupture, un joint n'a pas de loi.
**Aucun joint fracturable n'est paramétré dans l'acier ni dans le carbure.**

**[LU]** D3 p. 3 le dit d'ailleurs par la bande, en restreignant son propos :

> « In the FDEM model, **the intact rock** is discretised into tetrahedral
> elements in three dimensions. The elements next to each other are connected
> with the so-called joint elements. »

### 3.3 Un troisième traitement, encore différent

**[LU]** D1 p. 19, essai brésilien :

> « the fracture model is **only applied to the disc specimen**, and the **steel
> platens are assumed to be rigid**, which means there is no deformation in the
> platens, so **material properties are not needed for them**. The platens are
> only meshed for convenience reasons in the mesh creation; **the actual mesh is
> not involved in the stress calculations** »

Donc Imperial dispose d'au moins **trois régimes** pour un corps non fracturable :
**rigide** (D1, plateaux), **élastique sans joints** (D2/D3/D4, outil), et
**déformable avec joints** (la roche). Le lot 2c a montré que Solidity possède
formellement les deux formulations, déformable et corps rigide.

### 3.4 La conséquence pour rockim, chiffrée

Le dépôt pose **34 507 joints**, sur toute face interne des **trois** corps.
Imperial n'en pose que dans la roche. Deux effets :

1. **On amollit l'outil.** Un insert en carbure à 600 GPa tapissé de joints à
   pénalité finie n'a plus 600 GPa de raideur effective. Or c'est l'outil qui
   transmet l'impulsion : sa raideur conditionne la force de contact et le rebond,
   deux des sept critères de validation.
2. **On paie le pas de temps du carbure pour rien.** La condition de stabilité
   est gouvernée par min(h·√(ρ/E)) ; le carbure a le plus grand E et la plus
   petite maille (0,7 mm de surface d'insert). Les joints qu'on y met n'apportent
   aucune physique — le carbure ne casse pas dans ces essais — et coûtent en
   raideur locale.

**C'est le levier le moins cher du dépôt : supprimer des joints ne demande aucun
algorithme nouveau.**

---

## 4. LA PÉNALITÉ ET LA SOUPLESSE ARTIFICIELLE : assumée, jamais corrigée

**[LU]** D1 p. 12, éq. 6-7 — la pénalité entre dans le **durcissement pré-pic** :

    delta_np = 2 f_t h / p0        delta_sp = 2 f_s h / p0

où « **h is the mean length of the edges of the joint element** ; **p0 is the
penalty term characterising the stiffness of the joint element** », avec
(éq. 8) `lim_{p0→∞} delta_p = 0`.

**[LU]** D1 p. 12, le raisonnement complet — c'est la réponse frontale à la
question ouverte au dépôt depuis le 28/08 :

> « ideally the value of the penalty term p₀ should be **large enough so that the
> extra elasticity introduced into the domain by the joint elements can be
> negligible** (Klein et al., 2001; Turon et al., 2007). However, a **larger
> penalty term may cause numerical stability problems** (Schellekens and de
> Borst, 1993), which usually **requires smaller time-steps** in the explicit
> time integration scheme. For the current serial numerical code, the
> computational time would be unbearably long if the time-step is too small.
> Therefore, to maintain a balance between accuracy and computational efficiency,
> the value of the penalty term p₀ is usually chosen as » **E ≤ p₀ ≤ 10E** (éq. 9)

**Verdict sans ambiguïté** : ils **connaissent** la souplesse parasite, ils ne la
**corrigent pas** en relevant le module d'Young, ils la **bornent** en montant la
pénalité aussi haut que le pas de temps le permet. La compensation est
numérique, pas constitutive.

**Et rappel du [lot 2a](2026-08-29_lot2a_parametres_stanne.md) §3** : leur
pratique sur l'impact est **≈ 50 E** (3000 GPa pour E = 57 GPa sur St Anne),
soit **cinq fois au-dessus de leur propre borne haute**. La recommandation
E-10E vaut pour leurs cas quasi-statiques ; l'impact exige davantage.

**[LU]** D4 reprend le même arbitrage en 2D avec `k = 10E`, « avoiding excessive
computational costs ».

---

## 5. LE MAILLAGE ADAPTATIF N'EXISTE PAS — et voici la preuve

**[LU]** D1 p. 14, la phrase qui ferme le sujet :

> « the computational model proposed in this paper is **based on a fixed mesh**,
> so fractures can **only initiate and propagate along faces of tetrahedral
> elements**. »

**[LU]** D1 p. 39, en discussion :

> « The computational model described in this paper **works on fixed meshes**,
> which means it only allows fractures to propagate along tetrahedral element
> boundaries. This mesh dependency of fracture patterns is not significant if the
> element size is small enough compared with the scale of the domain (Guo et al.,
> 2016). »

**[LU] LA PREUVE DÉCISIVE** — D1 p. 41, section perspectives. Restreindre le
modèle de fracture au sous-domaine autour des fissures est présenté comme un
**travail futur** :

> « Similar approaches can be developed for the FEMDEM method, in which the
> fracture model is only applied [au sous-domaine autour des fissures] »

**On ne range pas en perspective ce que l'on a déjà implémenté.** Le raffinement
localisé n'existe pas dans ce code.

**[LU]** D2 p. 4 confirme que le raffinement est **a priori**, décidé au
maillage :

> « The mesh near the centre of the rock impact point was refined. This is
> because, in the FDEM simulation, **cracks are only allowed to propagate along
> the element boundaries**. »

### Le faux ami à ne pas confondre

**[LU]** D1 p. 15 emploie le mot « adaptively » **une seule fois**, et il porte
sur la **détection de contact**, jamais sur le maillage :

> « The **contact detection algorithm** is capable of **adaptively searching**
> tetrahedral elements to form and update contact couples »

---

## 6. CE QUI EST VRAIMENT « À LA DEMANDE » CHEZ EUX : la détection de contact

**[LU]** D1 p. 14, éq. 14-16 — l'optimisation réelle, et elle est élégante :

> « To improve the computational efficiency, **instead of finding contact couples
> everywhere in the domain, the contact detection in the continuum region (no
> fractures) is only activated after new fractures are formed**. »

et p. 15, le critère d'exclusion :

> « contact couples that **still have active joint elements between them are
> excluded** from this detection process. »

Après rupture d'un joint, le couple `C_ff = {tet+, tet−}` est ajouté (éq. 14),
puis six groupes G1…G6 sont formés autour des six nœuds du joint rompu (éq. 15)
et la détection ne tourne qu'à l'intérieur et entre ces groupes (éq. 16).

**[INFÉRÉ]** L'architecture d'Imperial déplace donc le coût de l'adaptativité de
la **géométrie** (remaillage, coûteux et fragile) vers la **topologie de
contact** (liste de couples, incrémentale). C'est un choix d'ingénierie, pas un
renoncement — et c'est une piste directement transposable à rockim, à examiner
au lot 4.

---

## 7. CE QUE COÛTE LE MAILLAGE FIXE — les chiffres publiés

**[LU]** D2 p. 5, Table 2 — le prix de la finesse, sur leur banc d'impact :

| maille fine (mm) | maille grossière (mm) | éléments | temps CPU (h) |
|---|---|---|---|
| 0,5 | 0,5 | 3 146 234 | **70,78** |
| 0,5 | 1 | 830 235 | 32,28 |
| **1** | **2** | **230 788** | **5,00** |
| 2 | 2 | 179 513 | 3,12 |
| 3 | 3 | 147 896 | 2,24 |
| 4 | 4 | 134 193 | 2,06 |

**Passer du maillage de production (1/2 mm, 5,00 h) au plus fin testé (0,5/0,5 mm,
70,78 h) coûte un facteur 14,2** — mais ce passage divise par deux **les deux**
tailles, fine et grossière, pas seulement la fine. *(Formulation corrigée le
2026-08-29 : la rédaction initiale — « diviser par deux la maille la plus fine
multiplie le coût par ≈ 14 » — attribuait à la seule maille fine un facteur qui
vient des deux. Le tableau ne contient aucun couple ne différant que par la
maille fine, donc l'effet propre de celle-ci n'est pas mesurable sur ces
données.)* Le
maillage de production (1/2 mm) est le troisième : **230 788 éléments, 5 h**.

**[LU]** D1 p. 40 impute la lourdeur au couple maillage fixe + code série :

> « the computational time can become **unaffordable** if the number of elements
> is significantly large, which is a **common issue for fixed-mesh-based
> numerical models in serial codes**. »

**[LU]** D2 p. 16 en tire la limite d'ambition :

> « this study simulates the impact behaviour of a **single insert**. Due to the
> high computational cost, **full-scale bit drilling tests would be extremely
> expensive**. »

**Pas de temps publiés** : D3, 10⁻⁹ s pour 373 066 tétraèdres ; D1, 2×10⁻⁹ s
(brésilien) et 5×10⁻⁸ s (polyaxial). **D2 n'en publie aucun** — le mot
« time step » n'y apparaît pas une fois.

---

## 8. L'OBJECTIVITÉ AU MAILLAGE : ce que le corpus démontre, et ce qu'il ne
## démontre pas

C'est le point où il faut être le plus honnête, dans les deux sens.

### 8.1 Ce que la formulation régularise, et ce qu'elle ne régularise pas

**[INFÉRÉ, dimensionnel, à partir de D1]** La loi cohésive est régularisée **en
énergie** — δ_c découle de G_f par l'éq. 10 (`G_f ≈ ⅓ f δ_c`) et **ne contient
pas h**. Mais la branche **pré-pic ne l'est pas** : `δ_np = 2 f_t h / p₀`
(éq. 6) fait entrer la longueur d'arête moyenne dans l'écrouissage.

### 8.2 Ce qu'ils assument, en le disant

* **D1 p. 39** : dépendance des faciès « **not significant** if the element size
  is small enough compared with the scale of the domain ».
* **D2 p. 4-5** : « a **small inherent dependence** of cracking resistance on the
  **chance orientations of the joint element directions** associated with
  unstructured mesh generation, **which is appropriate** » — coefficient de
  variation **≈ 2 %** sur les grandeurs de transfert d'énergie et **≈ 10 %** sur
  celles de fissuration.
* **D2 p. 13** : les surfaces de fissure « en zigzag » **raccourcissent** la
  fissure latérale simulée.
* **D3 p. 5** : « FDEM simulations rely on mesh boundaries for cracking, the
  generated rock fragments **always have sharp edges** ».
* **D6 (2026) §2.2 p. 4** : l'élément minimal de 1 mm est bien plus gros que la
  poudre mesurée à 0,035 mm ; descendre à cette taille multiplierait le nombre
  d'éléments par **plus de 23 000**.

### 8.3 Ce que le corpus ne démontre PAS — et c'est une vraie faiblesse

1. **D2 ne publie aucun résultat en fonction de la taille de maille.** Sa Table 2
   ne donne **que** le nombre d'éléments et le temps CPU. **Aucune grandeur
   physique, donc aucune courbe de convergence.** L'affirmation que le maillage
   retenu ne compromet pas la précision **n'est appuyée par aucune donnée
   publiée**.
2. **D3 conclut plus large que sa preuve.** Il écrit que le modèle n'est pas
   dépendant du maillage, sur la foi de **trois réalisations non structurées à
   taille de maille constante** — il n'a jamais fait varier la taille.
3. **D2 reconnaît lui-même la faiblesse** de son chiffre d'objectivité :
   « based on **just three different meshes** ».
4. **Aucune des trois sources** n'emploie les mots `objectivity`,
   `regularisation`, `characteristic length` ou `process zone length`. **Il
   n'existe dans ce corpus aucun critère chiffré comparable au `dx < ℓc/2` que la
   branche CONTINUUM de la thèse s'impose** (cf. [CONTINUUM.md](../../phd_geothermie/phd/CONTINUUM.md) §4).
5. **Aucune quantification** de la dépendance au maillage du **nombre de
   fragments** ni de l'**énergie dissipée**, alors que D2 admet qu'un fragment
   « may be a single tetrahedral element » — grandeur structurellement liée à la
   discrétisation.
6. **Aucun mass scaling** mentionné, **aucune frontière absorbante** : les effets
   de bord sont traités par la **taille de l'échantillon** seule.
7. **La seule vraie étude de sensibilité du modèle est hors du dossier** : D1 la
   renvoie à **Guo, Xiang, Latham & Izzuddin (2016)**, *Eng. Fract. Mech.*
   **151**, 70-91 — l'article dont le §2.4 de la thèse est déjà en fiche au dépôt
   ([`guo2014_s24_maillage.md`](guo2014_s24_maillage.md)).

> **À porter au crédit de rockim (lot 4).** Le dépôt a un garde-fou crack-band
> (`MatLaw.cpp:1304-1314`, `throw` si le plus gros élément dépasse E·G_f/f_t²) et
> une checklist d'objectivité écrite. **Imperial n'a rien de tel de publié.**
> Sur ce point le dépôt n'est pas en retard — il est plus discipliné.

---

## 9. LA DÉCOUVERTE DU LOT 3 : le granite Kuru est calibré DEUX FOIS,
## DIFFÉREMMENT, par la même équipe

**[LU]** Comparaison de D3 Table 2 (*JRMGE* 2025) et de D6 Table 1
(*IJRMMS* 206, 2026) — **même roche, même groupe, un an d'écart** :

| paramètre | D3 (2025) | D6 (2026) | écart |
|---|---|---|---|
| masse volumique (kg/m³) | 2630 | 2626 | ≈ 0 |
| **module d'Young (GPa)** | **67** | **60** | **−10 %** |
| Poisson | 0,26 | 0,24 | |
| **G_I (J/m²)** | **20** | **50** | **× 2,5** |
| **G_II (J/m²)** | **1500** | **1000** | **−33 %** |
| résistance en traction (MPa) | 11,4 | 10,98 | −4 % |
| **cohésion (MPa)** | **46,49** | **29,84** | **−36 %** |
| frottement interne | 1,96 | 1,85 | −6 % |
| **frottement glissant** | **0,39** | **0,18** | **−54 %** |

D3 crédite sa table à « (Saksala et al., 2014) » ; D6 recalibre tout.

### Et D3 dit pourquoi le frottement est bas — c'est un correctif numérique

**[LU]** D3 p. 5, verbatim :

> « the **sliding friction coefficient between rock fragments was also adjusted**
> to control the mutual friction between fragments. **Since FDEM simulations rely
> on mesh boundaries for cracking, the generated rock fragments always have sharp
> edges. When simulating the mutual friction between rock powders, it is
> necessary to reduce the sliding friction coefficient** to better match the
> ejection and splashing of rock powders. »

**Le coefficient de frottement glissant n'est pas une propriété physique de la
roche : c'est un bouton de compensation de l'angularité des tétraèdres.** Ils
l'écrivent. Et la valeur de ce bouton a été divisée par deux entre deux articles
sur la même roche, quand le modèle d'endommagement est arrivé.

### DEUX NUANCES apportées par la contre-vérification (2026-08-29, soir)

**1. Le µ bas n'est pas SEULEMENT un correctif de maillage — c'est aussi un
contraste voulu.** **[LU]** D3 p. 10 (imprimée 6104), verbatim :

> « The **abrupt change in the friction coefficient from intact internal to crack
> wall sliding (i.e. from 1.96 to 0.39)** leads to a **fast energy accumulation**
> as would be required to **initiate radial cracks**. The lower sliding friction
> coefficient also **facilitates the movement of fragments**, allowing rapid
> energy release. »

Le paragraphe ci-dessus était donc trop catégorique. Le coefficient joue **deux
rôles à la fois** : compenser l'angularité des tétraèdres (§9, cité plus haut)
**et** créer un saut 1,96 → 0,39 entre frottement interne de la roche intacte et
glissement sur lèvre de fissure, saut qui **concentre l'énergie nécessaire à
l'amorçage des radiales**. Le second rôle est physique, pas numérique. Cela ne
change rien à la non-transférabilité — au contraire, un contraste dépend des
**deux** valeurs, donc de la roche.

**2. Dans D3, l'effet de vitesse est absorbé dans G_I et G_II, PAS appliqué par
le DIF.** **[LU]** D3 p. 4 (imprimée 6098), verbatim :

> « the loading rate effect is also taken into consideration. **The energy release
> rates, G_I and G_II, were artificially increased to consider the loading rate
> effect through the validation process**, ensuring good agreement with the
> experimental results. In previous work, this approach has also been
> successfully applied to simulate dynamic Brazilian splitting tests and dynamic
> three-point bending experiments (Farsi, 2017), and rock blasting (Yang et al.,
> 2017). »

**Conséquence directe sur le tableau ci-dessus.** Les G_I = 20 et G_II = 1500 de
D3 ne sont **pas** des propriétés quasi-statiques : ce sont des valeurs
**gonflées par l'effet de vitesse**, calibrées. Celles de l'article de 2026
(G_I = 50, G_II = 1000) accompagnent un DIF appliqué **explicitement** à G_I et
G_II. **Les deux jeux ne désignent donc pas la même grandeur** — même symbole,
deux sens.

**[INFÉRÉ]** L'écart entre les deux calibrations du granite Kuru n'est donc pas
seulement une recalibration : c'est en partie un **changement de convention**,
l'effet de vitesse passant d'un pré-gonflement des énergies à une loi DIF
explicite. Je ne prétends pas expliquer le sens de l'écart (G_I monte de 20 à 50,
ce qui n'est pas la direction attendue d'un dégonflement) — seulement établir que
**comparer ces deux tableaux ligne à ligne n'a pas de sens.**

> **Conséquence, et elle est lourde.** Reprendre un coefficient de frottement
> d'Imperial, c'est reprendre une compensation calibrée pour **leur** maillage,
> **leur** géométrie d'élément et **leur** version du modèle. Ce n'est
> transférable ni entre roches (fiche 2026 §4), ni entre codes, ni même **entre
> deux articles du même groupe**. Cela renforce, s'il en était besoin, leur
> propre avertissement : les paramètres sont calibrés « as an **integrated
> parameter set** ».

---

## 10. RECOMMANDATION MOTIVÉE POUR ROCKIM

Le brief demande une recommandation, pas un inventaire. La voici, par ordre de
rapport gain/effort.

### R1 — Restreindre les joints à la roche. **À faire.**

*Ce que fait Imperial* : joints dans la roche seule ; outil élastique sans joints
(D4 p. 6870), voire rigide (D1 p. 19).
*Ce que fait rockim* : 34 507 joints sur les trois corps.
*Gain attendu* : raideur d'outil correcte, donc force de contact et rebond
corrects — deux des sept critères de validation ; et un pas de temps qui cesse
d'être payé pour du carbure qui ne casse pas.
*Effort* : faible. C'est un filtre au moment de la pose, pas un algorithme.
*Critère de réussite mesurable* : à deck constant, la vitesse de rebond du
taillant doit se rapprocher de la mesure, et le nombre de joints tomber d'un
facteur ≈ 2 sans changer le champ de fissuration dans la roche.

### R2 — Porter la pénalité de joint à ≈ 50 E et la mesurer. **À faire, en premier.**

> **⚠️ CETTE RECOMMANDATION EST FAUSSE. Voir [lot 4](2026-08-29_lot4_bilan_rockim.md) §4.**
> Deux erreurs. (1) Le dépôt avait déjà établi que sous la convention parabolique
> de Guo l'équivalent de leur p₀ = 3000 GPa vaut **p₀/(2E) = 26,32**, et le deck
> de réplique y est déjà — le « 20 » comparé ici est le DÉFAUT, pas la valeur de
> réplique. Monter à 50 aurait été deux fois trop raide. (2) L'équivalence 26,32
> annule h des deux côtés, ce qu'elle n'a pas le droit de faire : rockim mesure h
> comme le **diamètre inscrit** (6V/A) et Imperial comme la **longueur moyenne
> des arêtes**, soit un rapport 2,4495. L'équivalence correcte est
> **26,32 × 0,4082 = 10,74**. Il faut donc **descendre** la pénalité, ou changer
> la mesure de h. Le texte ci-dessous est conservé tel quel — on ne réécrit pas
> l'historique.

*Ce que fait Imperial* : **3000 GPa pour E = 57 GPa** sur St Anne (lot 2a §3),
soit 52,6 E — bien au-dessus de leur propre recommandation E-10E, qui vaut pour
le quasi-statique.
*Ce que fait rockim* : `jointPenaltyFactor = 20`, soit 2,6 fois moins.
*Pourquoi en premier* : la mesure du dépôt du 29/08 donne le levier pénalité
**4,4 fois** plus fort que le levier schéma d'insertion. C'est le paramètre le
plus sensible du problème et nous sommes du mauvais côté.
*Critère de réussite* : rejouer le banc St Anne à 20, 50 et 100 ; la courbe
force-pénétration et le bilan d'énergie doivent converger en montant, et
l'écart 20→50 doit être plus grand que l'écart 50→100.
*Risque à surveiller* : le pas de temps. C'est exactement l'arbitrage que D1
p. 12 décrit — mesurer le coût, pas le supposer.

### R3 — Abandonner la piste du maillage adaptatif. **Ne pas faire.**

Elle n'existe pas chez Imperial, et D1 la range en perspective de recherche
(p. 41). L'implémenter serait **dépasser** l'état de l'art, pas le rejoindre.
Ce n'est pas interdit — mais ce n'est plus une réplication, et le coût est sans
commune mesure avec R1 et R2. **À écarter du chemin critique.**

### R4 — Reprendre leur détection de contact paresseuse. **À évaluer.**

> **⚠️ À RELIRE AVEC UN CHIFFRE QUE LE DÉPÔT AVAIT DÉJÀ.** `DOCUMENTATION_rockim.md:183`
> mesure que sur la percussion longue le poste dominant est l'intégration exacte
> des 33 M de clips AVEC force (≈ 2/3 du run) ; les clips vides ne pèsent que
> ~12 %. R4 attaquerait donc **10-15 % du mur, pas un facteur 7**. C'est une
> optimisation légitime, pas la clé du CPU. Voir [lot 4](2026-08-29_lot4_bilan_rockim.md) §2.3.

*Ce que fait Imperial* : aucun couple de contact dans la région continue ; ils
naissent à la rupture d'un joint, par groupes autour de ses six nœuds
(D1 éq. 14-16), et les couples ayant encore un joint actif sont exclus.
*Gain attendu* : CPU, potentiellement important sur un domaine majoritairement
intact — ce qui est le cas d'un impact unique.
*Effort* : moyen. C'est une restructuration de la liste de couples, pas de la
physique.
*Critère de réussite* : temps par pas réduit **à résultats bit-identiques** —
c'est une optimisation, elle ne doit rien changer d'autre.

### R5 — Ne pas importer leurs coefficients de frottement. **Règle, pas tâche.**

§9 ci-dessus. Le frottement glissant est chez eux un correctif d'angularité de
maillage, recalibré à chaque changement de modèle. **Sur St Anne, garder 0,6**
(lot 2a §2). Ne jamais transférer un µ d'une roche, d'un maillage ou d'un article
à l'autre.

### R6 — Publier notre objectivité, puisqu'ils ne publient pas la leur.

Le corpus Imperial ne contient aucune courbe de convergence en taille de maille,
aucun critère de longueur caractéristique, et une étude d'objectivité que ses
propres auteurs qualifient de « based on just three different meshes ». Le dépôt
a un garde-fou crack-band et une checklist. **C'est un différenciateur du
manuscrit, pas une dette.** À écrire comme tel.

---

## 11. CE QUI RESTE OUVERT APRÈS LE LOT 3

| question | statut |
|---|---|
| l'étude de sensibilité au maillage du modèle de fracture 3D | **hors dossier** : Guo et al. (2016), *EFM* **151**, 70-91. À demander si l'objectivité devient un enjeu du manuscrit |
| la valeur du paramètre visqueux η de l'éq. 1 de D1 | **non publiée**, nulle part, et jamais mentionnée dans les cas dynamiques |
| le pas de temps de D2 | **non publié** — le modèle de validation le plus détaillé du corpus n'est pas reproductible sur ce point |
| le nombre de cœurs des temps CPU de la Table 2 de D2 | **non publié** : les 70,78 h et 2,06 h ne sont pas comparables entre eux sans lui |
| la taille de grain du calcaire St Anne | **non publiée** (celle du grès Rhune l'est : 0,2 mm ; celle du granite Kuru : 0,27-1,5 mm) |
| existence d'un mass scaling dans Solidity | **inconnue.** Zéro occurrence dans la littérature — ce qui ne prouve pas son absence du code |
