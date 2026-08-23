« VOLET 2 — Pathologies documentées du schéma cohésif extrinsèque et remèdes publiés »

Note de méthode : budget épuisé (8 appels : 6 recherches, 2 fetch dont 1 bloqué par le portail Springer). Chaque énoncé est marqué **[VERIFIE]** (page/abstract consulté dans ce run), **[TITRE VERIFIE]** (référence vue dans une liste de résultats, contenu non lu) ou **[MEMOIRE]** (connaissance propre, non revérifiée). Les paragraphes « lecture pour rockim » sont mon analyse, pas une citation.

---

## 1. La question posée, reformulée en termes de littérature

Le nuage diffus de rockim n'est pas un bug isolé : c'est la conjonction de **trois pathologies indépendantes** que la littérature cohésive a identifiées séparément entre 1996 et 2025, et qui s'additionnent exactement dans votre configuration (extrinsèque + explicite + contact par pénalité + champ lisse + seuils homogènes) :

1. **discontinuité temporelle à l'activation** (Papoulia–Sam–Vavasis) → oscillations et non-convergence en temps ;
2. **absence de « throttle » physique** : un champ lisse + un seuil unique = un ensemble d'activation de mesure non nulle, donc insertion « en tapis » (Zhou–Molinari, Molinari et al., Levy–Molinari) ;
3. **dérive énergétique du couple cohésif↔contact pénalisé** en explicite, qui produit de la **fragmentation artificielle** (arXiv 2511.14323, 2025).

Le point (2) est votre hypothèse, et elle est **étayée quantitativement** par la littérature. Le point (3) est une piste que vous n'aviez pas listée et qui est, en FDEM, au moins aussi probable que (2).

---

## 2. Pathologie A — Discontinuité temporelle à l'activation (le péché originel de l'extrinsèque)

**Papoulia, Sam & Vavasis (2003), *Time continuity in cohesive finite element modeling*, IJNME (doi 10.1002/nme.778)** et **Sam, Papoulia & Vavasis (2005), *Obtaining initially rigid cohesive finite element models that are temporally convergent*, Engineering Fracture Mechanics 72:2247–2267** [VERIFIE — noter la revue : *Eng. Fract. Mech.*, pas CMAME].

Contenu établi [VERIFIE, d'après l'abstract/résumé consultés] :
- ils introduisent la notion de **continuité temporelle** pour les modèles dits *initially rigid* — « une interface est inactive jusqu'à ce que la traction qui la traverse atteigne un niveau critique », c'est-à-dire **exactement le schéma de Yan et al. 2023 et de rockim** ;
- **thèse centrale : les méthodes de cette classe sont temporellement discontinues sauf disposition particulière** ;
- **conséquences nommées** : comportement **oscillatoire**, **non-convergence en temps** (raffiner Δt ne stabilise pas la solution), et **dépendance à des paramètres de régularisation non physiques** ;
- le papier de 2005 construit un cadre général, **exempt de régularisation**, produisant des modèles initialement rigides **temporellement convergents en dynamique explicite**.

Mécanisme (reconstitution, [MEMOIRE]) : au pas où l'on teste, la traction a **dépassé** ft (overshoot O(Δt·σ̇)). Si l'on insère en initialisant la loi cohésive à ft (ou pire, à la traction courante avec un ressort de rigidité finie), la force nodale **saute** ; ce saut est un choc numérique qui rayonne. Le remède PSV consiste à traiter l'instant d'activation comme une inconnue : interpoler l'instant exact de franchissement dans le pas, activer la facette **avec la traction qu'elle porte à cet instant**, et rendre la traction continue en temps à travers l'activation.

**Lecture pour rockim.** Votre balayage « à chaque pas, sans limite du nombre d'insertions » est le cas le plus exposé : c'est le cas où l'overshoot est minimal (bien) mais où le nombre d'événements de saut par pas est maximal (mal). Chaque insertion émet une micro-onde de décharge ; ces ondes vont **abaisser localement la contrainte chez certaines voisines mais l'augmenter chez d'autres** (concentration en pointe), ce qui alimente la cascade. Deux contrôles falsifiants peu coûteux : (i) faire varier Δt d'un facteur 4 à maillage fixe — si le compte de joints cassés bouge fortement, vous êtes dans la non-convergence temporelle décrite par PSV ; (ii) tracer l'histogramme de l'overshoot σ/ft − 1 au moment de l'insertion : s'il est large, l'activation n'est pas temporellement continue.

**Prolongements de la même école** [TITRE VERIFIE, contenu non lu] : *Non-differentiable energy minimization for cohesive fracture*, Int. J. Fract. (2017) ; *Energy minimization versus criteria-based methods in discrete cohesive fracture simulations*, Comput. Mech. (2021). Message [MEMOIRE] : choisir **quelles** facettes ouvrir par **minimisation d'énergie** plutôt que par franchissement de critère supprime la non-unicité de l'ensemble d'activation — c'est le remède « propre » au problème de l'ordonnancement, et il est directement pertinent quand des dizaines de facettes franchissent le seuil au même pas.

Voir aussi [TITRE VERIFIE] : *A two-field modified Lagrangian formulation for robust simulations of extrinsic cohesive zone models*, Comput. Mech. (2013) — remplacer la pénalité de « rigidité initiale infinie » par un multiplicateur de Lagrange, ce qui élimine la raideur artificielle avant activation.

---

## 3. Pathologie B — Intrinsèque : la localisation que vous enviez est peut-être un artefact

**Falk, Needleman & Rice (2001), *A critical evaluation of cohesive zone models of dynamic fracture*, J. Phys. IV Proc. 11(Pr5):43–50** [VERIFIE].
Résultats retenus :
- les prédictions sont **sensibles à la loi cohésive employée** ;
- les lois à **réponse élastique initiale** (= intrinsèques) produisent un **branchement spontané à grande vitesse**, **mais modifient les propriétés élastiques linéaires du corps** ;
- conséquence explicite : **l'espacement des surfaces cohésives ne peut pas être raffiné arbitrairement et devient une longueur caractéristique de la simulation**.

**Kubair & Geubelle (2003), *Comparative analysis of extrinsic and intrinsic cohesive models of dynamic fracture*, Int. J. Solids Struct. 40(15):3853–3868** [VERIFIE, niveau abstract] : les modèles **extrinsèques** (raideur initiale infinie, traction = résistance) sont **plus stables** que les intrinsèques (pente initiale finie) ; en propagation spontanée les intrinsèques sont **moins stables numériquement** ; en régime stationnaire certains modèles intrinsèques donnent des résultats **non physiques** (vitesse d'ouverture négative en pointe de zone cohésive).

**Lecture pour rockim — point de méthode important.** Comparer 28 000 (rockim, extrinsèque) à ~10 000 (MultiFracS, intrinsèque) suppose que les deux comptent la même chose. Ce n'est pas le cas :
- en **intrinsèque**, la population de cohésifs est **fixée à l'avance** ; le « nombre de joints cassés » est le nombre de cohésifs ayant atteint la séparation complète parmi cette population ;
- en **extrinsèque**, votre compteur est le nombre d'**activations**, dont une large fraction peut stagner à D ≪ 1.
**Avant toute modification du schéma, refaites le comptage sur l'observable commune** : facettes avec D ≥ 0,99 (ou ouverture ≥ δc, ou énergie dissipée ≥ 0,9 Gf). Si le 3× s'effondre, le problème n'est pas « trop d'insertions » mais « insertions qui ne coalescent pas » — un diagnostic entièrement différent, qui pointe vers la pathologie D (biais de maillage) plutôt que vers le seuil.

Et le corollaire inconfortable : d'après FNR, la localisation nette de l'intrinsèque **est en partie produite par la compliance artificielle des interfaces** (bandes molles préexistantes qui concentrent la déformation) et par le fait que l'espacement cohésif y agit comme longueur physique. Le schéma de référence n'est donc pas nécessairement « plus juste » : il est **régularisé par un artefact**. Votre remède ne doit pas être « imiter l'intrinsèque » mais « introduire l'hétérogénéité que l'intrinsèque introduit par accident, de façon contrôlée et physiquement calibrée ».

---

## 4. Pathologie C — Dérive énergétique cohésif↔contact pénalisé et **fragmentation artificielle** (le suspect n° 1 en FDEM)

**arXiv:2511.14323 (2025), *Stability of Extrinsic Cohesive-Zone Model with Penalty-Based Contact in Explicit Dynamic Fragmentation Simulations*** [VERIFIE — abstract intégral lu].

Résultats :
- sur un benchmark 2D, avec des réglages de pénalité de contact et de pas de temps **standards**, on observe une **croissance exponentielle de l'énergie** et une **fragmentation artificielle** en résultat ;
- trois mécanismes isolés et quantifiés : **(i)** la **raideur cohésive initiale divergente**, qui contraint le pas de temps stable ; **(ii)** les **sauts discontinus de raideur à l'interface cohésif↔contact** ; **(iii)** la **discontinuité introduite par l'adoucissement cohésif** ;
- estimations d'erreur analytiques, diagnostics en espace des phases et métriques de croissance d'énergie montrent que le **basculement répété cohésif↔contact accumule de petites erreurs d'énergie par pas en dérive à long terme** ;
- dans l'espace de paramètres exploré, la stabilité exige **des pas de temps très en dessous de la limite usuelle** ;
- remède évalué : **pénalité adaptative** liant la raideur de contact à la raideur cohésive évolutive → supprime la discontinuité et **restaure la conservation d'énergie**, mais autorise plus d'interpénétration → **outil de diagnostic, pas remède définitif** ;
- conclusion des auteurs : le **contact par pénalité n'est pas viable** pour des simulations de fragmentation énergétiquement cohérentes sur le long terme avec des statistiques de fragments physiquement signifiantes.

**Lecture pour rockim — c'est presque votre configuration exacte.** FDEM = contact par pénalité + cohésif adoucissant + explicite + insertion extrinsèque. Un tunnel profond est précisément un cas « long terme » (beaucoup de pas, beaucoup de basculements ouverture/fermeture sous confinement). **Un nuage diffus de fissures courtes est la signature attendue d'une injection d'énergie parasite** : de l'énergie ajoutée uniformément casse partout un peu, au lieu de nourrir une macro-fissure. Contrôle immédiat, sans rien changer au code : **bilan énergétique global** (cinétique + élastique + dissipée cohésive + frottement + travail des forces extérieures) en fonction du temps, et **énergie dissipée totale / (Gf × nombre d'insertions)**. Si le bilan dérive à la hausse et si l'énergie par insertion est très inférieure à Gf, la cause est ici et non dans le critère.

---

## 5. Pathologie D — Dépendance au maillage du nombre de fragments et du trajet de fissure

**Molinari, Gazonas, Raghupathy, Rusinek & Zhou (2007), IJNME 69:484–503** — titre exact : *The cohesive element approach to dynamic fragmentation: **the question of energy convergence*** [VERIFIE ; l'intitulé « mesh sensitivity » correspond au proceeding AIP Conf. Proc. 845:654, *Numerical Convergence of the Cohesive Element Approach in Dynamic Fragmentation Simulations* [TITRE VERIFIE]].

Résultats [VERIFIE] :
- question traitée : **convergence en énergie** de la solution EF pour la fragmentation à grande vitesse de chargement, en élasticité linéaire petites déformations ;
- fournit une règle pour **l'espacement propre des zones cohésives en fonction de la vitesse de chargement**, et une « feuille de route » pour choisir **les tailles de mailles ET les distributions de tailles de mailles** en 2D/3D ;
- **résultat clé : introduire un léger degré d'aléa de maillage (mesh randomness) améliore la convergence du problème de fragmentation jusqu'à DEUX ORDRES DE GRANDEUR.**

C'est, à ma connaissance, l'énoncé publié le plus direct de votre hypothèse : **sur un maillage régulier avec des propriétés homogènes, le problème de fragmentation ne converge pas ; le désordre (même faible) est ce qui sélectionne les sites de fissuration.**

**Zhou & Molinari (2004), *Dynamic crack propagation with cohesive elements: a methodology to address mesh dependency*, IJNME 59** [TITRE VERIFIE ; doi 10.1002/nme.857 ; contenu [MEMOIRE]] : le trajet et la vitesse de fissure dépendent de l'orientation des arêtes ; le remède passe par l'aléa (maillage et/ou propriétés).

Compléments [MEMOIRE] :
- **Papoulia, Vavasis & Ganguly (2006), IJNME** — maillages **pinwheel** : ils possèdent la propriété isopérimétrique, donc toute courbe est approchable par des chemins d'arêtes → **convergence spatiale de la nucléation/du trajet de fissure**, ce que ne garantit aucun maillage triangulaire quelconque.
- **Rimoli & Rojas (2015), Int. J. Fract.**, *Meshing strategies for the alleviation of mesh-induced effects in cohesive element models* — quantifie l'**anisotropie artificielle de ténacité effective** induite par le maillage : la ténacité apparente dépend de l'angle entre la direction physique de propagation et les orientations d'arêtes disponibles.
- **Vocialta, Richart & Molinari (2017), IJNME**, *3D dynamic fragmentation with parallel dynamic insertion of cohesive elements* [TITRE VERIFIE] — l'insertion dynamique en parallèle et ses coûts/artefacts.
- **Radovitzky, Seagraves, Tupek & Noels (2011), CMAME** [MEMOIRE] — hybride **Galerkin discontinu / éléments cohésifs** : les interfaces existent dès le départ, la continuité est imposée par pénalité intérieure **consistante** (donc sans compliance artificielle) et l'activation ne change plus la topologie en cours de calcul. C'est le remède structurel à la fois à la discontinuité d'activation (§2) et au saut de raideur (§4).

**Lecture pour rockim.** Un nuage de fissures **courtes** est aussi la signature d'un maillage qui ne sait pas représenter la direction de la macro-fissure : la fissure « paie » moins cher en se fragmentant en segments alignés sur des arêtes disponibles qu'en suivant la direction physique. Cela rejoint votre règle maison (mailler **non structuré**, cf. fiche ℓc DP-DFH). Contrôle : distribution des orientations des joints cassés — si elle épouse l'histogramme des orientations d'arêtes du maillage plutôt qu'une direction mécanique, c'est du biais de maillage.

---

## 6. Le résultat central demandé (c) : homogène ⇒ diffus ; il faut de l'hétérogénéité pour localiser

**Zhou & Molinari (2004), *Stochastic fracture of ceramics under dynamic tensile loading*, Int. J. Solids Struct. 41:6573–6596** [VERIFIE] :
- modèle de micro-fissuration avec **distribution stochastique des défauts internes** : **distribution de Weibull des résistances locales**, les résistances **des facettes** sont dispersées et suivent une Weibull ;
- **les structures fragiles sont vues comme des corps contenant des défauts initiaux, modélisés comme les facettes partagées par deux éléments voisins** ; chaque facette possède **sa** résistance et **sa** énergie de rupture ; quand la contrainte sur la facette dépasse **sa** résistance, la facette est **activée comme micro-fissure traitée en élément cohésif** — c'est mot pour mot le schéma extrinsèque de rockim, mais **avec un seuil par facette** ;
- **« grâce au concept de volume effectif de la théorie de Weibull, la dépendance au maillage indésirable du calcul numérique est significativement réduite »** ;
- application : cube de **SiC sous impact**, avec un **critère d'initiation en contrainte de cisaillement effective** et un **algorithme de contact frottant intégré à la procédure d'élément cohésif** ; balayage en vitesse de chargement pour obtenir les résistances en compression.

**Levy & Molinari (2010), *Dynamic fragmentation of ceramics, signature of defects and scaling of fragment sizes*, JMPS 58(1):12–26** [VERIFIE, niveau résumé] : les **défauts jouent un rôle critique** ; les fissures s'initient **à des positions apparemment aléatoires**, se propagent et coalescent pour former les fragments, avec des ondes de décharge ; l'étude **relie la distribution de défauts au nombre de fragments obtenu**, sur une large gamme de vitesses de déformation, comparée à la théorie, au numérique et à l'expérience existants. [MEMOIRE] : la conclusion complémentaire habituellement retenue est une **transition de régime** — à basse vitesse la population de défauts fixe la taille de fragment, à très haute vitesse l'inertie domine et la sensibilité aux défauts s'estompe.

Autres jalons [MEMOIRE] : Espinosa & Zavattieri (2003) — propriétés cohésives aléatoires au niveau des grains dans un polycristal ; Tijssens, Sluys & van der Giessen (2000) — surfaces cohésives extrinsèques pour le béton avec dispersion de résistance ; Grady–Kipp / Glenn–Chudnovsky — tailles de fragments par bilan d'énergie, référence à laquelle Molinari confronte ses simulations. Vu aussi en liste [TITRE VERIFIE, auteurs non vérifiés] : *Disorder effects in dynamic fragmentation of brittle materials*, JMPS (2003) ; *Mesh Objective Stochastic Simulations of Quasibrittle Fracture* (rapport OSTI).

**Énoncé analytique correspondant** (mon analyse, à assumer comme telle) : dans un champ **lisse** avec un seuil **uniforme**, l'ensemble {x : σ_eff(x,t) = f(x)} n'est pas un ensemble de points isolés — c'est **le bord d'une région entière qui franchit le seuil dans une fenêtre de temps de largeur ~ (gradient de contrainte)⁻¹·σ̇·Δt**. Le balayage à chaque pas convertit cette fenêtre en un **tapis d'insertions simultanées**, et le tapis rend impossible la sélection d'un mode : toutes les facettes s'adoucissent un peu, aucune ne draine l'énergie élastique du voisinage, la localisation ne s'amorce jamais. Introduire une dispersion de seuil de coefficient de variation CV transforme la fenêtre temporelle en **séquence** de franchissements isolés ; chaque insertion isolée décharge son voisinage (ondes de décharge de Levy–Molinari) et **inhibe** ses voisines : c'est le mécanisme de sélection. Le désordre n'est pas un décor, c'est **l'opérateur de sélection de mode**.

---

## 7. Catalogue des remèdes publiés, traduits en règles de déclenchement pour rockim

À implanter en **capacités opt-in additives** (principe VIII : ne rien retirer, défaut inchangé).

**R1 — Seuil par facette, Weibull, normalisé en volume/aire effectif(ve)** [Zhou–Molinari 2004, VERIFIE ; Pandolfi–Ortiz 2002 stockent déjà une limite par facette, acquis].
`insertion.strength_scatter = weibull`, module m, graine. Tirer ft_i (et c_i, corrélé, pour ne pas déformer l'enveloppe) par facette **à la création du maillage**, pas à l'exécution. **Point non négociable** : la **normalisation par l'aire effective** — ft_i = ft₀·(A₀/A_i)^{1/m}·(−ln U_i)^{1/m}. Sans elle, raffiner le maillage multiplie les maillons faibles et **augmente mécaniquement** le nombre d'insertions ; c'est exactement ce que la phrase « le concept de volume effectif réduit significativement la dépendance au maillage » sanctionne. Pour Red Bohus, m est à calibrer sur la dispersion mesurée de résistance en traction (typiquement m ≈ 5–15 pour un granite ; CV ≈ 1,2/m à 0,1 près [MEMOIRE]).

**R2 — Aléa géométrique de maillage** [Molinari et al. 2007, VERIFIE] : jitter des nœuds / distribution de tailles d'éléments, **jusqu'à 2 ordres de grandeur sur la convergence**. Bon marché, sans toucher la loi. Permet aussi de **séparer expérimentalement** le rôle du bruit géométrique de celui du bruit de résistance.

**R3 — Activation temporellement continue** [PSV 2003/2005, VERIFIE] : interpoler l'instant de franchissement dans le pas et initialiser la loi cohésive **à la traction effectivement portée**, de sorte que la traction soit continue à travers l'activation. À défaut, au minimum : borner l'overshoot admissible et sous-cycler le pas où il est dépassé.

**R4 — Ordonnancement / limitation du débit d'insertion** [Sam et al. 2005 ; Papoulia 2017/2021, cadre par minimisation d'énergie] : classer les candidates par **marge relative** (σ_eff/f − 1), n'insérer que les k plus critiques par pas, recalculer, itérer. C'est la contrepartie algorithmique du « throttle physique » de Weibull ; les deux se cumulent.

**R5 — Critère non local** [MEMOIRE ; convergent avec la partition nodale de Camacho–Ortiz, acquis] : évaluer la traction motrice comme moyenne pondérée sur un voisinage de rayon ~ ℓc plutôt que par facette isolée. Supprime la réponse point-à-point au bruit de contrainte et empêche le tapis d'être déclenché par des pics d'intégration.

**R6 — Retard / hystérésis** [MEMOIRE] : exiger que le critère soit satisfait pendant une durée ≥ t_dwell (ordre de ℓc/c) ou avec une marge, plutôt qu'à un instant. Variante physique et déjà présente chez vous : le **DIF** élève ft là où le chargement est rapide — c'est un throttle dépendant du taux, à préférer à un throttle purement numérique.

**R7 — Vérification espacée** [Pandolfi–Ortiz 2002, acquis] : tester tous les N pas. Attention : **R7 sans R3 aggrave** l'overshoot et donc la discontinuité temporelle. Les deux vont ensemble.

**R8 — Gardes énergie/contact** [arXiv 2511.14323, VERIFIE] : (i) plafonner la raideur cohésive initiale ; (ii) **pénalité de contact adaptative** asservie à la raideur cohésive courante pour supprimer le saut au basculement (diagnostic, cf. interpénétration accrue) ; (iii) réduire Δt sous la limite usuelle ; (iv) **auditer le bilan d'énergie en continu** et refuser toute conclusion sur des statistiques de fragments si l'énergie dérive.

**R9 — Critère d'initiation en cisaillement effectif avec frottement, sous confinement** [Zhou–Molinari 2004, VERIFIE pour l'existence du critère et du contact frottant intégré ; forme exacte MEMOIRE, cf. Camacho–Ortiz 1996, acquis] : en compression confinée (tunnel profond), un critère σ_n ≥ ft OU |τ| ≥ c − σ_n tanφ appliqué à une contrainte moyennée sur deux éléments **surestime le cisaillement moteur** partout dans la zone plastique. Utiliser une traction effective unique σ_eff = √(τ²/β² + ⟨σ_n⟩²) en traction et σ_eff = |τ| − μ|σ_n| en compression, avec le **même** β/μ que celui de la loi cohésive post-insertion, garantit la continuité entre critère d'insertion et loi de comportement — l'incohérence entre les deux est une cause classique d'insertions qui ne s'ouvrent jamais.

**R10 — Remède structurel** [MEMOIRE] : formulation **DG/cohésif hybride** (Radovitzky et al. 2011) ou multiplicateur de Lagrange à deux champs (Comput. Mech. 2013, [TITRE VERIFIE]) — supprime la raideur artificielle avant activation et le changement de topologie en cours de calcul. Lourd, mais c'est la seule voie qui traite A et C à la racine.

---

## 8. Plan de falsification proposé (peu coûteux, avant toute réécriture)

1. **Recompter** sur observable commune (D ≥ 0,99) et non sur les activations — cf. §3. Peut à lui seul dissoudre le facteur 3.
2. **Histogramme de D** des joints insérés. Un mode piqué à D ≪ 1 = tapis d'insertions avortées → pathologie B/§6. Un mode bimodal = localisation partielle.
3. **Bilan d'énergie** et **dissipation par insertion / Gf** → teste §4 (arXiv 2511) sans rien modifier.
4. **Δt ÷ 4** à maillage fixe → teste la non-convergence temporelle de PSV.
5. **Balayage m** : m = ∞ (homogène actuel) → m = 15 → m = 8, maillage inchangé, normalisation par aire active. Si les blocs apparaissent, l'hypothèse « le seuil homogène met tout au seuil simultanément » est confirmée expérimentalement dans votre code.
6. **Contrôle croisé** : jitter de maillage seul, seuils homogènes (Molinari 2007). Sépare bruit géométrique et bruit matériau.
7. **Orientations** des joints cassés vs orientations d'arêtes → teste le biais de maillage (§5).

Les points 1 à 4 ne coûtent que du dépouillement sur les runs existants.


## Lacunes

Budget de 8 appels épuisé (6 WebSearch + 2 WebFetch, dont un bloqué par le portail d'authentification Springer). N'ont PAS pu être vérifiés dans ce run et restent marqués [MEMOIRE] ou [TITRE VERIFIE] :

1. Seagraves & Radovitzky 2010, « Advances in Cohesive Zone Modeling of Dynamic Fracture », chapitre 12 de « Dynamic Failure of Materials and Structures » (Springer, ISBN 978-1-4419-0446-1) : la page Springer redirige vers idp.springer.com (303). Le chapitre existe en accès libre sur academia.edu (URL vue dans 4 des 6 recherches) — à lire en priorité au prochain run : c'est la revue qui synthétise intrinsèque/extrinsèque + dépendance au maillage + remède DG.
2. Radovitzky, Seagraves, Tupek & Noels 2011 CMAME (hybride DG/cohésif) : non vérifié, référence donnée de mémoire (année et liste d'auteurs à confirmer).
3. Papoulia, Vavasis & Ganguly 2006 IJNME (maillages pinwheel, propriété isopérimétrique) : non vérifié.
4. Rimoli & Rojas 2015 Int. J. Fract. (anisotropie artificielle de ténacité induite par le maillage) : non vérifié.
5. Zhou & Molinari 2004 IJNME 59 « Dynamic crack propagation with cohesive elements: a methodology to address mesh dependency » : titre et doi (10.1002/nme.857) vus en liste, contenu non lu — la « méthodologie » exacte proposée reste de mémoire.
6. Levy & Molinari 2010 JMPS : seul le résumé de recherche a été obtenu ; la transition de régime « défauts dominants à basse vitesse / inertie dominante à haute vitesse » est de mémoire.
7. Kubair & Geubelle 2003 : lu au niveau d'un résumé de moteur de recherche, pas de l'abstract éditeur ; la condition quantitative sur la pente initiale intrinsèque (rapport delta_c/delta_e admissible) n'a pas été retrouvée.
8. Sam, Papoulia & Vavasis 2005 : le mécanisme algorithmique précis de l'activation temporellement continue (interpolation de l'instant de franchissement, itération sur l'ensemble d'activation) est reconstitué de mémoire, pas lu.
9. arXiv:2511.14323 : abstract intégral lu et fiable, mais aucune valeur numérique (facteur de réduction de Δt requis, seuils de raideur) n'a pu être extraite — la lecture du PDF était interdite.
10. Non couvert faute d'appels : Espinosa & Zavattieri 2003, Tijssens–Sluys–van der Giessen 2000, Klein et al. 2001 (pénalité de raideur intrinsèque sur la vitesse d'onde), Molinari–Zhou–Ramesh sur l'effet des propriétés matériau, et toute la littérature FDEM/Munjiza spécifique sur l'insertion adaptative postérieure à Yan 2023.
11. Deux références vues en liste mais dont les auteurs n'ont pas été vérifiés : « Disorder effects in dynamic fragmentation of brittle materials » (JMPS 2003) et « Mesh Objective Stochastic Simulations of Quasibrittle Fracture » (OSTI) — ne pas citer sans vérification d'auteurs.
12. Correction utile issue de ce run : le titre exact de Molinari et al. 2007 IJNME 69:484-503 est « ... the question of ENERGY CONVERGENCE » (et non « mesh sensitivity », qui correspond au proceeding AIP Conf. Proc. 845:654) ; et Sam-Papoulia-Vavasis 2005 est paru dans Engineering Fracture Mechanics 72:2247-2267 (et non CMAME).
