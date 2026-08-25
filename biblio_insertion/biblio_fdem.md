# VOLET 4 — L'insertion cohésive dans les codes FDEM : règles de déclenchement, manière d'insérer, gardes et remèdes

*Convention de marquage : **[V]** = page/notice consultée dans ce run ; **[M]** = connaissance propre non revérifiée ici ; **[D]** = dérivation ou raisonnement personnel explicite (à vérifier avant citation).*
*Budget de 8 appels web épuisé ; rédaction avec l'acquis. Aucun PDF ouvert.*

---

## 1. Cadre : trois familles de règles, pas deux

La littérature FDEM oppose rituellement « intrinsèque » et « extrinsèque », mais du point de vue des **règles de déclenchement** il y a en réalité trois régimes distincts, et c'est le troisième qui explique le problème observé sur rockim.

| Régime | Quand la facette est créée | Qui décide | Effet sur le champ |
|---|---|---|---|
| **I — Intrinsèque (Munjiza 2004)** | à t=0, sur **toutes** les arêtes | personne (topologie) | réseau de ressorts souples partout : complaisance + bruit spatial permanents |
| **II — Extrinsèque ordonné (Camacho–Ortiz 1996)** | à la volée, **arête par arête**, la plus chargée d'abord, champ réactualisé | l'ordonnancement | front de fissure qui avance ; décharge propagée avant le test suivant |
| **III — Extrinsèque « en balayage »** (Yan 2023, rockim) | à la volée, **toutes les arêtes au seuil au même pas**, testées contre le **même** champ non relaxé | le pas de temps | tapis d'insertions simultanées |

Le régime III n'est pas « l'extrinsèque » : c'est une **variante Jacobi** de l'extrinsèque, alors que Camacho–Ortiz et Pandolfi–Ortiz sont des variantes **Gauss–Seidel** (acquis du volet précédent). **[D]** C'est, à mon sens, le diagnostic central de ce volet : le nuage diffus n'est pas une conséquence de « l'extrinsèque » en soi mais de l'absence de réactualisation du champ à l'intérieur d'une passe de balayage. Toutes les facettes d'une bande plastique testent un tenseur de contrainte qui n'a encore été relaxé par aucune des insertions de la même passe ; elles franchissent donc le seuil **ensemble**, et la décharge élastique qui aurait dû protéger les voisines arrive un pas trop tard.

---

## 2. Munjiza 2004 — l'intrinsèque historique et sa pénalité

**[M]** Munjiza, *The Combined Finite-Discrete Element Method*, Wiley 2004 ; loi de fissure issue de Munjiza, Andrews & White 1999, IJNME 44:41–57 (« combined single and smeared crack model »).

Règles :
- **Où** : élément cohésif d'épaisseur nulle à 4 nœuds sur **chaque** arête interne du maillage triangulaire, inséré au préprocesseur **[V]** (confirmé par la lignée Y-Geo/Irazu : *« zero-thickness, four-noded cohesive crack elements are inserted a priori along all interior edges »*).
- **Quand** : jamais — il n'y a pas de critère de déclenchement, seulement une variable d'endommagement D qui démarre.
- **Seuil** : la loi a une **branche montante élastique** σ = (2δ/δ_p − (δ/δ_p)²)·f_t jusqu'à δ_p = 2·h·f_t/p, puis une branche adoucissante z(D) calée sur béton (paramètres a≈0,63 ; b≈1,8 ; c≈6,0). **[M]**
- **Garde** : la pénalité p (distincte de la pénalité de contact p_c) est le seul garde-fou ; recommandée « grande devant E ».

Le point qui compte pour vous : **la branche montante n'est pas physique, c'est un artefact de régularisation**, et sa raideur p/h est un paramètre libre. Deux recommandations contradictoires coexistent dans la littérature courante, ce qui prouve à soi seul qu'il s'agit d'un bouton numérique et non d'une propriété matériau :
- *« penalty values two orders of magnitude larger than the input elastic moduli »* (p ≈ 100 E) — lignée Lisjak/Grasselli **[V]** ;
- *« penalty values are set equal to 10 times the Young's modulus value »* — autre calibration FDEM **[V]**.

---

## 3. Lignée Yan Chengzeng / MultiFracS (2016 → 2026)

**[V]** Yan, Zheng et coll. ; logiciel **MultiFracS**, FDEM multiphysique parallélisé GPU, appliqué à la dessiccation des sols, la rupture des roches tendres, l'injection, le tir de mine, la fracturation hydraulique, le transfert thermique de contact, la fissuration thermo-induite, la géothermie et le couplage THM. L'adaptatif de 2023 y est implémenté.

Trois jalons à distinguer soigneusement — ils ne portent **pas** sur le même objet :

1. **Yan, Zheng & Wang 2023, IJRMMS 169:105439** — *A 2D adaptive finite-discrete element method for simulating fracture and fragmentation in geomaterials* **[V]**. Formulation retenue par rockim : *« the new FDEM adaptively inserts the cohesive elements when the stress exceeds a critical value »*, par opposition explicite au FDEM conventionnel qui pose les cohésifs partout d'avance. **[M]** Le nouveau modèle constitutif est introduit conjointement à l'insertion.

2. **Yan et coll. 2023, Computers and Geotechnics** — *Implementation of extrinsic cohesive zone model (ECZM) in 2D FDEM using **node binding scheme*** (S0266352X23002276) **[V]**. C'est le papier de **mécanique d'insertion** : les nœuds sont dupliqués d'emblée mais **liés** (contrainte de liaison nodale) ; « insérer » = relâcher la liaison. Formulation clé relevée dans ce run **[V]** : *« extrinsic cohesive elements directly enter the strain-softening stage after insertion and do not bear material elastic deformation, which guarantees the continuum behavior of numerical models prior to fracture onset and thus can readily avoid the artificial compliance problem »*. **[D]** Lue à l'envers, cette phrase est exactement votre problème : le continuum est **parfait** avant l'insertion, donc parfaitement lisse, donc sans germe de localisation.

3. **Yan et coll. 2024, Computers and Geotechnics** — *Modelling dynamic fracture and fragmentation of rocks under multiaxial coupled static and dynamic loads with a parallelised 3D FDEM* (S0266352X24004191) **[V]**, et **Yan et coll. 2023, Eng. Fract. Mech.** sur les géomatériaux à forte teneur en inclusions avec stratégie de reconstruction **[V]**. Ce dernier est significatif pour votre diagnostic : la lignée Yan **réintroduit de l'hétérogénéité explicite** (inclusions, GBM) dès qu'elle veut des motifs de rupture réalistes.

**[M]** Rappel de vos propres notes : la viscosité μ de Yan est le critique de Munjiza 2h√(Eρ) — c'est un amortissement de l'interface, pas un throttle d'insertion. Aucun mécanisme d'ordonnancement ou de limitation du nombre d'insertions par pas n'est décrit dans la lignée Yan à ma connaissance.

---

## 4. Lisjak & Grasselli — Y-Geo / Irazu

**[V]** Mahabadi, Lisjak, Munjiza & Grasselli, *Y-Geo: New Combined Finite-Discrete Element Numerical Code for Geomechanical Applications* (Int. J. Geomech., 2012) ; **Irazu** (Geomechanica Inc.), présenté à l'ARMA 2016 comme *« a new fully-parallel finite-discrete element code »*, avec couplage thermique et hydro-mécanique.

**Position sur l'insertion : strictement intrinsèque.** **[V]** Aucun des deux codes ne fait d'insertion adaptative ; les CCE sont a priori sur toutes les arêtes, et la complaisance artificielle est traitée par le seul réglage de p (≈100 E, cf. §2). Leur argument de mésindépendance porte sur la **trajectoire** (*« no prior assumptions for crack trajectory are necessary ; any arbitrary fracture trajectory can be captured within the limitations of the mesh topology »* **[V]**), pas sur l'énergie.

**[M]** Compléments de la même lignée, à vérifier : Tatone & Grasselli (2015, IJRMMS) proposent une procédure de calibration 2D qui montre que E macroscopique, ν et l'UCS **émergents** dépendent conjointement de p et de la taille de maille ; Lisjak & Grasselli (2014, JRMGE) en font la revue. Y-Geo apporte en outre un critère de Mohr–Coulomb avec coupure en traction sur les CCE et l'initialisation des contraintes in situ — pertinent pour votre tunnel profond.

**[V] À lire en priorité** : *Investigation of the influence of model control parameters on fracture characteristics of GPU parallel FDEM*, Geomech. Geophys. Geo-energ. Geo-resour. 2023 (Springer, 10.1007/s40948-023-00651-y). C'est, dans ce que j'ai vu, le titre qui vise le plus directement le lien **paramètres numériques (dont la pénalité) → caractéristiques de fissuration**.

---

## 5. Fukuda et coll. — FDEM GPGPU : la seule lignée qui a documenté un « mode de fissuration parasite »

Chronologie **[V]** (titres et revues vérifiés, contenus détaillés **[M]**) :

- **2019, IJNAMG** — *Development of a GPGPU-parallelized hybrid FDEM for modeling rock fracture* ; et **2019** — *GPGPU-parallelized 3D combined FDEM modelling of rock fracture with **adaptive contact activation approach*** ; **2019, RMRE** — simulateur 3D quasi-statique/dynamique.
- **2021, IJRMMS** — *Modelling of dynamic rock fracture process using FDEM with a novel and efficient **contact activation scheme*** (S1365160921000332). **Règle exacte relevée [V]** : le *semi-ACAA* *« adaptively activates contact calculations for continuum solid elements around the cohesive element **which has just been subjected to shear softening while its softening function has just satisfied a prescribed threshold** »*. Appliqué au marbre en BTS et UCS dynamiques.
- **2024, Computers and Geotechnics** — *Development of a GPGPU-parallelized FDEM based on **extrinsic** cohesive zone model with **master–slave algorithm*** (S0266352X23006997) **[V]**.
- **2024/2025, Computers and Geotechnics** — *Development of a GPGPU-parallelized **3-D** FDEM with a novel and simple implementation of extrinsic cohesive zone model* (S0266352X24005822) **[V]** — fetch refusé (403), abstract non lu.

**Deux enseignements pour vous.**

(i) **Distinguer activation de contact et insertion de cohésif.** Le semi-ACAA n'insère pas de cohésif : il décide **quand allumer le calcul de contact** autour d'une facette déjà en adoucissement. La règle de déclenchement est **relative à un événement voisin**, pas à un seuil de contrainte absolu évalué partout. **[D]** C'est structurellement un throttle « de proche en proche » : rien ne s'active loin d'une facette déjà endommagée. Transposée à l'insertion, cette règle interdirait mécaniquement le tapis.

(ii) **La motivation explicite est un mode de fissuration parasite.** La synthèse du moteur de recherche attribue à cette lignée l'idée que l'approche adaptative *« overcomes spurious fracturing mode associated with FDEM simulations »* **[V, mais synthèse composite — à revérifier sur l'article]**. C'est la formulation la plus proche de votre observation que j'ai trouvée, mais elle va dans le sens **inverse** : chez Fukuda l'adaptatif *corrige* le mode parasite, chez vous il le crée.

---

## 6. HOSS / LANL

**[V]** *HOSS: an implementation of the combined finite-discrete element method*, Comput. Particle Mech. 2020 ; Rougier, Munjiza et coll., *The combined plastic and discrete fracture deformation framework for finite-discrete element methods*, IJNME 2020 ; rapports LANL Rougier–Knight–Munjiza 2013 (HOSS–MUNROU) et Knight–Rougier–Lei 2015 (version éducative). Applications : coalescence de fissures dans le granite (RMRE 2019), impact de plaque volante dans le granite, gouge de faille.

**Règle d'insertion** : **[V]** la description publique reste intrinsèque — *« specimens are modeled as an assembly of elastic bulk elements interconnected by cohesive elements, simulating **inherent material flaws** that evolve into potential cracks »*, l'élément étant supprimé quand D=1. Noter le vocabulaire : les cohésifs intrinsèques sont explicitement présentés comme **des défauts matériels préexistants**, c'est-à-dire comme de l'hétérogénéité — voir §8.

**[V] Piste distincte et récente** : brevet US 12 135 925, *Libraries-based explicit fracture and fragmentation framework*, qui formalise un **critère d'insertion de fracture discrète** sous la forme d'un champ scalaire φ dépendant d'un **paramètre de transition α**, du tenseur de Cauchy C et du tenseur de déformation ε. **[D]** Le « paramètre de transition » est précisément le degré de liberté qui manque à rockim : un critère qui n'est pas un simple seuil binaire sur σ_n, mais une fonction pondérée qui peut porter une mémoire ou une hystérésis.

**[V]** Rougier et coll. 2020 (plasticité + fracture discrète) est le remède structurel côté LANL : si le continuum **plastifie** avant d'atteindre le seuil de traction, la zone plastique se décharge et cesse d'être au seuil partout simultanément. C'est exactement ce que fait votre VUMAT DP-DFH côté Abaqus — et ce que rockim, avec un continuum élastique, ne fait pas.

---

## 7. (b) Comparaisons publiées des motifs de fissuration intrinsèque vs extrinsèque

Ce que j'ai pu vérifier :

- **[V]** *Comparative analysis of extrinsic and intrinsic cohesive models of dynamic fracture*, Int. J. Solids Struct. 2003 (S0020768303001719, indexé côté Illinois/groupe Paulino — **auteurs non revérifiés**). Conclusions relevées : *« intrinsic models are less numerically stable than the extrinsic ones »* ; *« under steady-state propagation conditions, some intrinsic cohesive models lead to unrealistic results as the crack opening velocity becomes negative at the cohesive zone tip »* ; *« the intrinsic model introduces artificial compliance depending on the area of cohesive element surfaces introduced and the cohesive element property »*.
- **[V]** Yan et coll. 2023 (node binding) : l'extrinsèque *« guarantees the continuum behavior prior to fracture onset »* et *« readily avoids the artificial compliance problem »*.
- **[V]** Un travail 2026 de FDEM **grain-based** pour la fracturation hydraulique en roche cristalline rapporte que l'ECZM *« avoided the artificial stress oscillations characteristic of the ICZM »* (référence exacte non capturée).
- **[V]** *Unified cohesive zone model (UCZM) for fracturing and fragmenting solids*, Eng. Fract. Mech. 2024 (S0013794424007616) — non lu, mais le titre promet précisément l'unification intrinsèque/extrinsèque. **À lire.**

**Constat honnête et important : toutes les comparaisons publiées que j'ai vues sont à charge contre l'intrinsèque, et sur des critères de complaisance, d'oscillations et de stabilité — jamais sur la topologie du réseau de fissures (blocs vs nuage).** Votre observation (extrinsèque → nuage diffus, intrinsèque → blocs) n'a pas d'équivalent publié dans la lignée FDEM d'après ce run. Les seuls travaux de tonalité voisine sont côté fragmentation dynamique : **[M]** Molinari et coll. 2007 (IJNME, *« the question of energy convergence »*) sur la non-convergence de la taille de fragment en cohésif ; et **[M, confiance moyenne]** Zhang, Paulino & Celes 2007 (IJNME) sur le micro-branchement excessif en extrinsèque. **C'est donc un espace de publication ouvert pour vous.**

---

## 8. (c) La souplesse de pénalité intrinsèque comme hétérogénéité / complaisance numérique implicite

Personne, dans ce que j'ai vu, ne l'écrit ainsi. Ce qui est écrit :
- **[V]** la complaisance artificielle est reconnue, quantifiée par « l'aire des surfaces cohésives introduites et les propriétés de l'élément », et **combattue** (p = 10 E à 100 E selon les auteurs) ;
- **[V]** HOSS présente ses cohésifs intrinsèques comme *« simulating inherent material flaws »* — c'est-à-dire, en toutes lettres, comme un **modèle de défauts**, donc d'hétérogénéité ;
- **[V]** en FDEM intrinsèque, *« macroscopic deformation prior to fracture is governed by the **coupled** stiffness of the solid bulk and the cohesive interfaces »* (FDEM CAES 2026) : le module macroscopique est un objet composite, pas le module d'entrée.

**[D] Estimation d'ordre de grandeur (ma dérivation, à vérifier).** Avec la loi de Munjiza, la raideur initiale d'un joint vaut k = p/h. Le supplément de déformation apporté par un joint traversé est σ/p ; sur une longueur L on traverse N = L/h joints, d'où un supplément total Lσ/p : la **complaisance moyenne est indépendante du maillage** et vaut 1/E_eff = 1/E + 1/p, soit −9 % de module pour p = 10 E et −1 % pour p = 100 E. En revanche, la **fluctuation spatiale** ne l'est pas : dans un maillage non structuré, la densité locale d'arêtes normales à une direction donnée fluctue, et la fluctuation relative de raideur dans une fenêtre contenant n arêtes est de l'ordre de (E/p)·n^(−1/2). Pour p = 10 E et n ≈ 4–10 arêtes dans une fenêtre de taille ℓc, cela donne **3 à 5 % de dispersion du module effectif, spatialement corrélée à la topologie du maillage** — soit l'équivalent d'un Weibull de module m ≈ 15–30 sur la contrainte locale, exactement l'ordre de grandeur du « throttle physique » de Zhou & Molinari (acquis).

Autrement dit : **l'intrinsèque contient gratuitement un bruit de champ de l'ordre du pourcent, et c'est probablement lui qui localise chez MultiFracS-intrinsèque ; l'adaptatif de Yan, en supprimant proprement la complaisance, supprime aussi ce bruit et livre un continuum trop parfait.** Le prix de l'exactitude est la perte du germe. **Cette formulation ne se trouve dans aucune référence consultée — c'est une hypothèse à tester, et c'est le cœur potentiel de votre contribution.**

---

## 9. Remèdes, classés par rapport coût/effet pour rockim

**Rang 1 — Passer le balayage de Jacobi à Gauss–Seidel (coût nul, effet attendu maximal).** Dans une passe, trier les facettes candidates par excès décroissant (σ_n − f_t, ou |τ| − (c − σ_n tanφ)), insérer la première, **réactualiser les forces internes des deux éléments voisins**, retester. Filiation directe Camacho–Ortiz. Variante bon marché si la réactualisation complète est trop coûteuse : **non-maximum suppression** — n'insérer, dans un voisinage de rayon ℓc, que la facette d'excès maximal, reporter les autres au pas suivant. **[D]** Non publié en FDEM à ma connaissance.

**Rang 2 — Verrou de voisinage temporisé.** Après insertion d'une facette, geler ses arêtes adjacentes pendant Δt ≈ h/c_d (temps d'arrivée de l'onde de décharge). Justification physique explicite ; analogue de l'activation de proche en proche de Fukuda (semi-ACAA) **[V]**.

**Rang 3 — Diagnostic avant remède : l'histogramme d'insertions par pas.** Compter n_ins(t). Une signature en rafales (10³ insertions sur un pas, puis rien) prouve le tapis ; une signature lissée l'infirme et renvoie le problème vers le critère lui-même. Corollaire : mesurer l'**overshoot** σ_n/f_t au moment de l'insertion — si la médiane dépasse nettement 1, le dépassement dissipé instantanément injecte une énergie parasite qui alimente la cascade.

**Rang 4 — Insertion temporellement continue.** **[M]** Papoulia, Sam & Vavasis (IJNME 2003, *Time continuity in cohesive finite element modeling*) et Sam, Papoulia & Vavasis (EFM 2005) montrent qu'insérer sur un test au pas de temps crée une **discontinuité temporelle de la traction**, source d'oscillations et de non-convergence quand Δt → 0 ; le remède est d'insérer à l'instant exact du franchissement (interpolation sous-pas). Référence de premier plan pour justifier proprement un throttle.

**Rang 5 — Réinjecter le bruit supprimé.** Weibull par facette (m ≈ 15–30, cf. §8) sur f_t et c. Coût nul, opt-in, conforme à votre règle de croissance par addition. **C'est le test falsifiant de l'hypothèse du §8** : si le motif redevient blocky à m ≈ 20, la complaisance intrinsèque est bien le germe.

**Rang 6 — Le test croisé décisif.** rockim sait faire les deux schémas. Faire tourner le **même** tunnel en intrinsèque avec p = 10 E, 30 E, 100 E, 300 E. Si le motif se re-diffuse quand p croît (complaisance → 0), l'hypothèse est démontrée expérimentalement dans votre propre code, sans dépendre d'aucune référence. **C'est la figure qui vend l'article.**

**Rang 7 — Rendre le continuum non lisse.** Loi bulk adoucissante ou plastique (DP-DFH côté rockim), à la manière du cadre plasticité + fracture discrète de Rougier et coll. 2020 **[V]**. Coûteux mais c'est le remède physique, pas numérique.

**Rang 8 — Non-localité du critère.** Vous moyennez la traction sur les **2** éléments voisins : c'est une longueur non locale minimale. Élargir la moyenne à un disque de rayon ℓc/2 favorise la sélection d'un maximum unique par bande.

**Rang 9 — Restaurer la branche compression–cisaillement–friction** de Camacho–Ortiz (acquis), essentielle en tunnel profond et en percussion où σ_n < 0 domine.

---

## 10. Références (état de vérification)

**[V] — notices consultées dans ce run**
- Yan, Zheng & Wang, IJRMMS 169 (2023) 105439 — adaptive 2D FDEM. `S1365160923001132`
- Yan et coll., Comput. Geotech. (2023) — ECZM in 2D FDEM, **node binding scheme**. `S0266352X23002276`
- Yan et coll., Comput. Geotech. (2024) — parallelised 3D FDEM, multiaxial static+dynamic. `S0266352X24004191`
- Yan et coll., Eng. Fract. Mech. (2023) — géomatériaux à fortes inclusions. `S0013794423001297`
- Fukuda et coll., IJRMMS (2021) — contact activation scheme (semi-ACAA). `S1365160921000332`
- Fukuda et coll., Comput. Geotech. (2024) — ECZM GPGPU, master–slave. `S0266352X23006997`
- Fukuda et coll., Comput. Geotech. (2024/25) — 3-D FDEM, simple ECZM. `S0266352X24005822` *(403, non lu)*
- Fukuda et coll., IJNAMG 2019 (10.1002/nag.2934) ; RMRE 2019 (10.1007/s00603-019-01960-z)
- Mahabadi, Lisjak, Munjiza & Grasselli — **Y-Geo**, Int. J. Geomech. 2012 ; **Irazu**, ARMA-2016-516 ; geomechanica.com/software/irazu
- *Influence of model control parameters on fracture characteristics of GPU parallel FDEM*, Geomech. Geophys. Geo-energ. 2023, 10.1007/s40948-023-00651-y — **prioritaire**
- **HOSS**, Comput. Particle Mech. 2020, 10.1007/s40571-020-00349-y ; Rougier et coll., IJNME 2020, 10.1002/nme.6255 ; coalescence granite RMRE 2019 10.1007/s00603-019-01773-0
- Brevet US 12 135 925 — critère d'insertion φ(α, C, ε)
- *Comparative analysis of extrinsic and intrinsic cohesive models of dynamic fracture*, IJSS 2003, `S0020768303001719` *(auteurs non revérifiés)*
- *Unified cohesive zone model (UCZM)*, Eng. Fract. Mech. 2024, `S0013794424007616` — **à lire**
- Paulino, Celes, Espinha & Zhang — *A general topology-based framework for adaptive insertion of cohesive elements in FE meshes*, Eng. with Computers 2008, 10.1007/s00366-007-0069-7
- *Adaptive insertion of cohesive elements for delamination in laminated composites*, Eng. Fract. Mech. 2015, `S0013794415004208`

**[M] — à revérifier avant citation**
- Munjiza, *The Combined Finite-Discrete Element Method*, Wiley 2004 ; Munjiza, Andrews & White, IJNME 44 (1999) 41–57
- Tatone & Grasselli, IJRMMS 2015 — procédure de calibration 2D FDEM (pénalité, taille de maille)
- Lisjak & Grasselli, J. Rock Mech. Geotech. Eng. 2014 — revue
- Papoulia, Sam & Vavasis, IJNME 2003 ; Sam, Papoulia & Vavasis, EFM 2005 — continuité temporelle de l'insertion
- Molinari et coll., IJNME 2007 — convergence énergétique en fragmentation cohésive
- Zhang, Paulino & Celes, IJNME 2007 — micro-branchement en extrinsèque *(confiance moyenne)*
- Radovitzky, Seagraves, Tupek & Noels, CMAME 2011 — Galerkine discontinu + cohésif extrinsèque **sans changement topologique** : élimine la complaisance intrinsèque *sans* insertion à la volée. Alternative architecturale sérieuse pour rockim.


## Lacunes

LACUNES ET POINTS NON VERIFIES (budget 8 appels epuise, aucun PDF ouvert)

1. Aucun abstract d'article lu integralement. Le WebFetch sur Fukuda 2024/25 (S0266352X24005822) a renvoye HTTP 403 ; ScienceDirect bloque. Toutes les citations de contenu proviennent de snippets de moteur de recherche, potentiellement composites (le moteur agrege plusieurs sources dans une meme phrase). A revalider avant toute citation dans le manuscrit, en particulier la phrase attribuee a la lignee Fukuda sur le "spurious fracturing mode".

2. Auteurs de "Comparative analysis of extrinsic and intrinsic cohesive models of dynamic fracture" (Int. J. Solids Struct. 2003, S0020768303001719) NON identifies. Indexe cote Illinois Experts, donc probablement groupe Paulino, mais non confirme. Ne pas citer avec des noms.

3. Regle d'insertion exacte de Yan 2023 IJRMMS non lue : je n'ai pas verifie si le balayage est Jacobi (tout au meme pas contre le champ non relaxe) ou s'il comporte un ordonnancement. Le diagnostic central du rapport (regime III = variante Jacobi) repose sur la description que VOUS donnez de rockim, pas sur une lecture de l'article. A verifier dans l'article original avant de le presenter comme une critique de Yan.

4. Aucun papier POSTERIEUR a 2023 critiquant explicitement l'adaptatif de Yan n'a ete trouve. Les papiers 2024-2026 de la lignee (3D multiaxial, GBM hydraulique) l'etendent sans le remettre en cause. Il est possible qu'une telle critique n'existe pas encore — ce qui est en soi un resultat exploitable.

5. Le point (b) est la plus grosse lacune : AUCUNE comparaison publiee des MOTIFS de fissuration (blocs vs nuage diffus) entre intrinseque et extrinseque n'a ete trouvee en FDEM. Toutes les comparaisons vues portent sur la complaisance, les oscillations de contrainte et la stabilite numerique. Les rapprochements avec Molinari 2007 (convergence energetique) et Zhang-Paulino-Celes 2007 (micro-branchement) sont [MEMOIRE] et le second est de confiance moyenne.

6. Le point (c) n'est etaye par AUCUNE reference explicite. Personne, dans les sources consultees, n'ecrit que la penalite intrinseque agit comme une heterogeneite implicite qui localise. Les elements verifies (complaisance artificielle reconnue ; HOSS decrivant ses cohesifs comme "inherent material flaws" ; raideur macro "couplee" bulk+interface) sont des indices convergents, pas une affirmation publiee. L'estimation quantitative du bruit (E/p)·n^(-1/2) ~ 3-5 % equivalent Weibull m~15-30 est MA derivation a partir de la loi de Munjiza (delta_p = 2 h ft / p) : la formule de delta_p elle-meme est de memoire et doit etre reverifiee dans Munjiza 2004 avant usage.

7. Papoulia/Sam/Vavasis (continuite temporelle de l'insertion) : references de memoire, non verifiees dans ce run, mais ce sont les references les plus importantes du volet pour justifier un throttle. A verifier en priorite lors du prochain run web.

8. Non couvert faute de budget : le detail des regles d'insertion dans HOSS (le brevet US 12135925 mentionne un critere phi(alpha, C, epsilon) mais je n'ai pas verifie qui en est le titulaire ni s'il correspond a HOSS) ; la position exacte d'Irazu commercial (une capacite d'insertion adaptative a-t-elle ete ajoutee depuis 2020 ?) ; le contenu du UCZM 2024 (Eng. Fract. Mech.) qui promet l'unification intrinseque/extrinseque et est probablement LA reference manquante du volet.

9. Prochain run web recommande (6-8 appels) : (a) UCZM Eng. Fract. Mech. 2024 ; (b) Papoulia/Sam/Vavasis ; (c) Geomech. Geophys. Geo-energ. 2023 sur les parametres de controle FDEM ; (d) le papier GBM-FDEM hydraulique 2026 citant ICZM vs ECZM ; (e) verification des auteurs de l'IJSS 2003.