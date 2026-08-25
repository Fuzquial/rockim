# VOLET 5 — Côté roches : hétérogénéité et localisation

*Question posée à ce volet : que dit la littérature « mécanique des roches » sur le fait qu'un milieu homogène diffuse et qu'un milieu hétérogène localise ? Avec quelles valeurs de m de Weibull, et avec quelles longueurs de corrélation quand le champ est corrélé ?*

**Marquage** : [VÉRIFIÉ] = page/résumé consulté dans ce run ; [MÉMOIRE] = connaissance propre non revérifiée (les 8 appels autorisés ont été consommés : 5 recherches, 3 fetches dont 1 redirect et 1 refus 403 ScienceDirect).

---

## 1. Le diagnostic mécanique, avant les références

Le problème de rockim n'est pas « l'insertion adaptative est fausse ». Il est que **le critère d'insertion est appliqué à un champ lisse et, dans la zone plastique, quasiment confondu avec la surface de charge du bulk**. Dans une zone où le matériau est à la limite d'écoulement *par construction*, toutes les arêtes vérifient le critère au même pas : l'insertion devient dégénérée (« tapis »), et le nombre de joints cassés (28 000 vs 10 000) mesure la taille de la zone plastique, pas un motif de fracturation.

La littérature roches apporte deux remèdes distincts, souvent confondus :
- **(R1) briser la simultanéité** par un champ de résistance qui n'est pas constant (Weibull, grains, litage) — le « throttle physique » de Zhou & Molinari déjà acquis, mais ici justifié par la physique de la roche et non par la numérique du cohésif ;
- **(R2) créer de la traction locale sous compression globale**, que le continuum lisse ne produit pas — c'est l'apport spécifique des modèles à grains (Lan, Martin & Hu).

Un schéma intrinsèque (MultiFracS, Wang et al. 2024) obtient R1 *gratuitement* : les interfaces cohésives pré-insérées introduisent une compliance discrète, orientée par le maillage, qui perturbe le champ avant même toute rupture. L'insertion extrinsèque adaptative, elle, préserve l'exactitude du continuum — et perd le bruit qui localisait. [MÉMOIRE, argument mécanique]

---

## 2. (a) Preuves publiées : homogène → diffus, hétérogène → localisé

### 2.1 Tang 1997 / RFPA — le pionnier, et l'argument le plus direct

- **Tang, C.A. (1997), « Numerical simulation of progressive rock failure and associated seismicity », IJRMMS 34(2), 249–261** [VÉRIFIÉ pour le titre/la revue ; pagination MÉMOIRE]. Le milieu est un maillage d'éléments élastiques-endommageables ; un élément « casse » (module réduit, résistance résiduelle) quand un critère de Mohr-Coulomb à *tension cut-off* est atteint. **La seule source de localisation du modèle est la distribution de Weibull affectée élément par élément** au module d'Young et à la résistance. [VÉRIFIÉ : « rock heterogeneity is considered by assuming that certain mechanical properties, such as Young's modulus and strength of the elements within a model, conform to a Weibull distribution »]
- Le corollaire est mécaniquement inévitable et c'est exactement le cas rockim : **en compression uniaxiale d'un bloc homogène, le champ est uniforme ; tous les éléments franchissent le seuil au même incrément.** Il n'y a alors ni nucléation ni bande — la rupture est soit uniforme (nuage), soit entièrement dictée par les conditions aux limites (frettage des plateaux, singularité de coin). Le RFPA existe *parce que* le continuum homogène ne localise pas. [MÉMOIRE — c'est le raisonnement, la formulation explicite « un modèle homogène ne localise pas » n'a pas été retrouvée mot pour mot dans les pages consultées ; le papier RRFPA 2025 consulté n'affirme que « l'hypothèse de propriétés uniformes surestime substantiellement la résistance et la déformation élastique » [VÉRIFIÉ]]
- Suite directe : **Tang & Kaiser (1998), « Numerical simulation of cumulative damage and seismic energy release during brittle rock failure — Part I: Fundamentals », IJRMMS** [VÉRIFIÉ existence via ScienceDirect] — l'émission acoustique simulée est simplement le comptage des éléments cassés ; le passage « beaucoup de petits événements dispersés » → « peu de gros événements corrélés » est **piloté par m**. C'est la métrique jumelle de votre « nuage diffus vs macro-fissures ».
- **Tang et al. (2000), « Numerical studies of the influence of microstructure on rock failure in uniaxial compression — Part I: effect of heterogeneity », IJRMMS 37, 555–569** [MÉMOIRE] : l'étude paramétrique de référence sur m ; transition d'un comportement quasi-ductile/diffus (m faible) à un comportement fragile avec une macro-bande unique (m élevé).
- **Tang et al., fracturation progressive autour d'une excavation sous compression biaxiale** [VÉRIFIÉ existence] : le même moteur reproduit les *breakouts* en V et le zonage de la zone endommagée — analogue direct de votre tunnel profond.

### 2.2 FDEM argilites — Lisjak, Grasselli & Vietor 2014

- **Lisjak, A., Grasselli, G., Vietor, T. (2014), « Continuum–discontinuum analysis of failure mechanisms around unsupported circular excavations in anisotropic clay shales », IJRMMS 65, 96–115** [VÉRIFIÉ titre/volume/pages]. Le point capital pour rockim : **l'hétérogénéité y est *orientationnelle et déterministe*, pas statistique** — la résistance des éléments cohésifs (ft, c) est modulée selon l'angle entre l'arête et le litage. Cela suffit à briser la simultanéité et à produire des macro-fractures organisées (rupture par flexion/délaminage côté parallèle au litage, cisaillement ailleurs), donc une EDZ en chevrons conforme aux observations du Mont Terri.
- Suite : **Lisjak et al. (2015), « Hybrid Finite-Discrete Element Simulation of the EDZ Formation and Mechanical Sealing Process Around a Microtunnel in Opalinus Clay », RMRE** [VÉRIFIÉ existence] ; et **Lisjak et al. (2015), TUST, tunnel circulaire en Opalinus** [VÉRIFIÉ existence].
- **Leçon transposable au granite Red Bohus** : le granite est isotrope à l'échelle du bloc, donc la modulation par le litage n'existe pas — il faut la remplacer par une hétérogénéité de **grain** (§2.4), sinon rockim n'a strictement aucune source de dissymétrie.

### 2.3 BPM — Potyondy & Cundall 2004

- **Potyondy, D.O. & Cundall, P.A. (2004), « A bonded-particle model for rock », IJRMMS 41(8), 1329–1364** [MÉMOIRE]. Le BPM ne contient *aucun* tirage de résistance : l'hétérogénéité est **géométrique** (distribution de rayons, coordination variable, orientations de contact aléatoires). Elle produit un champ de forces de contact fortement fluctuant → nucléation sur les contacts les plus sollicités → coalescence en bandes. C'est la démonstration qu'une hétérogénéité *de structure* remplace une hétérogénéité *de résistance*.
- Défauts connus et documentés du BPM disque/sphère : **rapport UCS/résistance en traction trop bas (≈3–5 au lieu de 10–25 pour un granite) et enveloppe de rupture trop linéaire** [MÉMOIRE]. Ces deux défauts ont motivé les GBM.

### 2.4 Grain-based models — la source de traction locale

- **Cho, Martin & Sego (2007), « A clumped particle model for rock », IJRMMS 44** [MÉMOIRE] : agréger les particules en clumps restaure le ratio UCS/T.
- **Potyondy (2010), grain-based model (Voronoi + joints inter-grains)** [MÉMOIRE] : deux jeux de propriétés, intra-grain vs joint de grain, avec le joint plus faible.
- **Lan, Martin & Hu (2010), JGR 115, B01202, « Effect of heterogeneity of brittle rock on micromechanical extensile behavior during compression loading »** [MÉMOIRE] — **la référence la plus pertinente pour votre problème** : sous compression *purement* compressive à l'échelle macro, le contraste d'élasticité et la géométrie polygonale des grains génèrent des **contraintes de traction locales** aux joints de grains. Dans un continuum homogène, ces tractions n'existent pas ; seul le critère de cisaillement peut déclencher, et il se déclenche partout à la fois. Cela explique mécaniquement pourquoi votre insertion adaptative « en tapis » est *cisaillante et diffuse* plutôt que *tensile et localisée*.
- **Modèles Voronoi 2D continuum pour le Mine-By tunnel** [VÉRIFIÉ existence, RMRE 2023] : la tessellation seule (sans champ aléatoire) suffit à faire apparaître le notch de spalling.

### 2.5 FDEM tunnels durs — Vazaios, Diederichs, Vlachopoulos

- **Vazaios, Vlachopoulos & Diederichs (2019), « Assessing fracturing mechanisms and evolution of the excavation damaged zone of tunnels in interlocked rock masses at high stresses using a finite-discrete element approach », JRMGE** [VÉRIFIÉ titre ; contenu non lu — ScienceDirect a renvoyé 403].
- **Vlachopoulos & Vazaios (2018), Advances in Civil Engineering, « The Numerical Simulation of Hard Rocks for Tunnelling Purposes at Great Depths: A Comparison between the Hybrid FDEM Method and Continuous Techniques »** [VÉRIFIÉ existence].
- Résultat rapporté : le FDEM (Irazu) reproduit la **position et la profondeur du notch** du Mine-By Experiment de l'URL, en accord avec la microsismique enregistrée [VÉRIFIÉ via résumé de recherche]. Le mécanisme de spalling y est explicitement tensile/extensile, pas cisaillant — cohérent avec §2.4. Les travaux du groupe combinent DFN (échelle massif) **et** hétérogénéité à l'échelle du grain [MÉMOIRE].

---

## 3. (b) Valeurs de m employées et effet sur nombre / longueur des fissures

**Définition RFPA** : `f(x) = (m/x₀)·(x/x₀)^(m−1)·exp(−(x/x₀)^m)`, avec x = E ou σ_c d'un élément, x₀ le paramètre d'échelle et **m l'*indice d'homogénéité*** — **m grand = homogène** [VÉRIFIÉ : « a larger value of m (m>1) means the parameters are more concentrated around their mean, the material is relatively more homogeneous »]. ⚠️ Piège de vocabulaire à surveiller dans la littérature : certains auteurs parlent d'*heterogeneity index* pour le même symbole, ce qui inverse le sens des phrases — un des extraits recueillis dans ce run est ambigu sur ce point.

| Plage de m | Régime | Effet sur le motif |
|---|---|---|
| m ≈ 1,1–2 | très hétérogène (roche altérée, béton, charbon) | endommagement **diffus**, pic très abaissé et arrondi, quasi-ductile, très nombreux événements de faible magnitude |
| m ≈ 3–6 | granite/marbre « courant » en RFPA | **régime cible** : nucléation multiple *puis* coalescence en une ou deux macro-bandes |
| m ≈ 10–20 | roche « homogène » | pic élevé et brutal, très peu d'événements précurseurs, une macro-fissure unique |
| m ≳ 50–100 | homogène numérique | dégénérescence : rupture quasi simultanée pilotée par les conditions aux limites — **le cas rockim actuel** |

[valeurs et bornes : MÉMOIRE ; le fetch RRFPA 2025 ne donne **aucune** valeur de m et le confirme explicitement [VÉRIFIÉ]]

Éléments corroborés dans ce run [VÉRIFIÉ, extraits de résumés] : « la distribution des résistances se concentre quand m augmente… le matériau est plus homogène » ; « le motif de fissuration au pic est diagonal pour les indices faibles, en Y renversé pour les indices élevés » ; « l'angle de rupture par cisaillement augmente avec l'indice d'homogénéité ».

**Conséquences pratiques, et le contre-sens à éviter** :
1. **L'hétérogénéité n'est pas monotone.** Trop peu → tapis simultané ; **trop → re-diffusion**, car chaque élément faible casse pour son propre compte sans jamais entraîner ses voisins. Le motif « blocs découpés » vit dans une fenêtre intermédiaire de m. Ne pas répondre au problème rockim par « le plus de dispersion possible ».
2. **La longueur des fissures est gouvernée par la corrélation, pas par la variance.** Une pointe de fissure ne progresse que si l'arête suivante est faible *aussi*. Avec un tirage indépendant par arête (bruit blanc), la probabilité de trouver une arête au-dessus de la moyenne est ≈ 1/2 à chaque pas → **arrêt après quelques arêtes → exactement le nuage de fissures courtes observé**. D'où le §4.
3. **Effet d'échelle Weibull** : σ ∝ V^(−1/m). Avec m = 5, doubler le volume élémentaire abaisse la résistance apparente de ≈ 13 %. Introduire de l'hétérogénéité **abaisse la résistance macro** : il faudra **recalibrer la moyenne** pour conserver le calage GBM Red Bohus (ft = 34, c = 13,6, φ_joint = 13,4°), sinon la validation Abaqus (+7,9 % de pic filtré) est perdue.
4. **m méso ≠ m mesuré en labo.** Les dispersions expérimentales d'UCS sur éprouvettes donnent des m apparents élevés (≈ 10–30) ; le m méso du RFPA est un **paramètre de calage**, à ajuster sur la forme de la courbe post-pic et sur le taux d'AE, pas à mesurer directement [MÉMOIRE].
5. **Répartir la dispersion sur ft et c, peu sur E.** Disperser E perturbe la rigidité, la vitesse d'onde et donc le pas critique — rédhibitoire en percussion. Ordre de grandeur recommandé : CoV(ft) ≈ 0,2–0,4 (m ≈ 3–6), CoV(E) ≤ 0,10–0,15. [MÉMOIRE]
6. **Critique moderne du tirage RFPA** : tirer E et σ_c du même aléa impose une corrélation parfaite entre les deux — le papier RRFPA 2025 reproche justement à la démarche Weibull classique de « ne pas capturer la corrélation intrinsèque entre ces propriétés » [VÉRIFIÉ].

---

## 4. (c) Longueurs de corrélation spatiale — le levier décisif

**Le point clé, et probablement la correction la plus rentable pour rockim** : un tirage indépendant par élément/arête est un **bruit blanc dont la longueur de corrélation vaut la taille de maille**. Conséquences : (i) l'hétérogénéité n'est **pas objective au maillage** (raffiner change la statistique du milieu, pas seulement la discrétisation) ; (ii) elle ne crée aucun *défaut* de taille physique capable de piloter une macro-fissure. Il faut un **champ aléatoire corrélé** (gaussien ou lognormal), caractérisé par une **longueur de corrélation θ** issue de la microstructure.

Valeurs publiées relevées :

| Source | Champ | Longueur de corrélation | Statut |
|---|---|---|---|
| **RRFPA (RMRE 2025), « The Random RFPA Method for Modelling Rock Failure »** | E **et** UCS, fonction d'autocorrélation **exponentielle-cosinus**, δx et δy | **5 mm** (isotrope, δx = δy) | [VÉRIFIÉ — page lue] |
| **Modèle de champ aléatoire basé nanoindentation (RMRE 2025)** — granite | E, champ gaussien, moyenne **65,64 GPa**, **CoV 0,18** | **0,57 mm** | [VÉRIFIÉ via extrait de résumé] |
| **« Using correlated random fields for modeling the spatial heterogeneity of rock »** | E, résistance | facteur de longueur de corrélation introduit **dans la loi de Weibull** pour corréler les propriétés avec la distance entre éléments | [VÉRIFIÉ : existence + résumé ; valeurs non lues] |
| Modélisation stochastique en longue taille (IJGE 2025) | E et UCS, « scale of fluctuation » | échelle de fluctuation prédéfinie | [VÉRIFIÉ via extrait] |
| Ordre de grandeur pour un granite à grain moyen 1–3 mm (type Red Bohus) | ft, c | **θ ≈ 1–5 mm**, soit 1 à 2 tailles de grain | [MÉMOIRE, extrapolation] |
| **Blair & Cook (1998), IJRMMS, analyse statistique de la fracturation en compression** | résistance locale corrélée | montre que c'est la corrélation, pas la variance seule, qui contrôle la localisation | [MÉMOIRE, non revérifié] |

**Règles d'usage** [MÉMOIRE, pratique standard des champs aléatoires] :
- résoudre le champ : **dx ≤ θ/3 à θ/5**. Au-delà, le champ dégénère en bruit blanc et on retombe sur le cas 4.2 ;
- attention à la convention : pour une ACF exponentielle, l'**échelle de fluctuation de Vanmarcke θ_V = 2δ**, où δ est la longueur de corrélation ; les articles mélangent les deux ;
- **effet attendu et testable** : θ ↗ ⇒ moins de sites de nucléation *indépendants*, chacun pilotant une zone de taille θ ⇒ **fissures plus longues, moins nombreuses**. C'est le paramètre qui doit convertir vos 28 000 fissures courtes en un réseau de macro-fissures découpant des blocs, **à nombre total d'insertions décroissant**.

---

## 5. Traduction en règles de déclenchement / d'insertion pour rockim

| # | Règle | Justification |
|---|---|---|
| **A** | **Un seuil par ARÊTE, tiré une fois et gelé** (ft_e, c_e), jamais un seuil global. Cohérent avec la limite de traction stockée par facette de Pandolfi & Ortiz et le Weibull par facette de Zhou & Molinari (volet acquis). | §2.1, §3 |
| **B** | **Tirer un champ CORRÉLÉ**, pas un bruit blanc : bruit blanc convolué par un noyau gaussien de portée θ (ou Karhunen-Loève / moyennes mobiles), évalué aux centres d'éléments, puis **affecté à l'arête par la moyenne des 2 éléments adjacents** — la même moyenne que rockim utilise déjà pour projeter la traction. Cohérence exacte entre l'évaluateur et le seuil. | §4 |
| **C** | Alternative structurelle (isotrope, granite) : **partition Voronoi + `groupBond`** (déjà implémenté d'après la mémoire projet) → ft/c distincts pour arêtes intra-grain vs inter-grain. Reproduit R2 (traction locale) et pas seulement R1. | §2.3–2.4 |
| **D** | **Disperser ft et c, très peu E** (CoV(E) ≤ 0,15) : en percussion, disperser E déplace le pas critique et disperse la vitesse d'onde. | §3.5 |
| **E** | **Re-calibrer la moyenne** après introduction de l'hétérogénéité (effet d'échelle Weibull V^(−1/m)) pour retrouver ft = 34 / c = 13,6 macroscopiques. | §3.3 |
| **F** | **Découpler le critère d'insertion de la surface de charge du bulk.** Si l'insertion se déclenche exactement là où le matériau plastifie, elle est dégénérée par construction. Un sur-seuil aléatoire (règle A) ou un critère additionnel (déformation/ouverture critique) restaure une hiérarchie temporelle entre arêtes. | §1 |
| **G** | **θ est physique, dx est numérique.** Le seul test d'objectivité valide : raffiner dx **à θ constant** et vérifier l'invariance du motif (nombre de macro-fissures, taille des blocs). Un motif qui change avec dx à θ fixe signale que le bruit vient du maillage, pas du matériau. | §4 |
| **H** | **Diagnostic quantitatif à instrumenter** : (i) nombre d'insertions ; (ii) **longueur moyenne des chaînes connexes d'arêtes cassées** ; (iii) **fraction d'arêtes cassées isolées** (aucune voisine cassée). Signature « nuage diffus » = fraction d'isolées élevée + longueur de chaîne ≈ 1–3 arêtes. Balayer (m, θ) et tracer longueur de chaîne vs θ : c'est la courbe qui prouvera ou réfutera l'hypothèse. | §3.2, §4 |

**Contre-indications** : (1) l'hétérogénéité ne remplace pas une longueur interne — si l'énergie de rupture du cohésif est mal dimensionnée, on déplace le problème ; (2) m trop faible re-diffuse ; (3) θ plus grand que la zone de process masque tout effet ; (4) toute comparaison chiffrée à Wang et al. 2024 doit être refaite après recalibrage (règle E).

---

## 6. Références (avec statut)

- Tang (1997) IJRMMS 34(2), *Numerical simulation of progressive rock failure and associated seismicity* — [VÉRIFIÉ titre/revue]
- Tang & Kaiser (1998) IJRMMS, *Cumulative damage and seismic energy release, Part I: Fundamentals* — [VÉRIFIÉ existence]
- Tang et al. (2000) IJRMMS 37, 555–569, *Influence of microstructure, Part I: effect of heterogeneity* — [MÉMOIRE]
- Zhu & Tang (2004) RMRE 37(1), *Micromechanical model for simulating the fracture process of rock* — [MÉMOIRE]
- *The Random RFPA Method for Modelling Rock Failure*, RMRE (2025), PMC12048433 — [VÉRIFIÉ, page lue : θ = 5 mm, ACF exponentielle-cosinus]
- *Nanoindentation-Based Random Field Model for Fracture of Heterogeneous Rock*, RMRE (2025) — [VÉRIFIÉ via résumé : granite, E = 65,64 GPa, CoV 0,18, θ = 0,57 mm]
- Lisjak, Grasselli & Vietor (2014) IJRMMS 65, 96–115 — [VÉRIFIÉ titre/volume/pages]
- Lisjak et al. (2015) RMRE, EDZ microtunnel Opalinus ; Lisjak et al. (2015) TUST, tunnel circulaire — [VÉRIFIÉ existence]
- Potyondy & Cundall (2004) IJRMMS 41(8), 1329–1364, BPM — [MÉMOIRE]
- Cho, Martin & Sego (2007) IJRMMS 44, clumped particle model ; Potyondy (2010) GBM — [MÉMOIRE]
- Lan, Martin & Hu (2010) JGR 115 B01202 — [MÉMOIRE]
- Vazaios, Vlachopoulos & Diederichs (2019) JRMGE, EDZ interlocked rock masses, FDEM — [VÉRIFIÉ titre, contenu non lu — 403]
- Vlachopoulos & Vazaios (2018) Adv. Civ. Eng., FDEM vs continu à grande profondeur — [VÉRIFIÉ existence]
- Modèles Voronoi 2D continuum, Mine-By Experiment, RMRE (2023) — [VÉRIFIÉ existence]
- Blair & Cook (1998) IJRMMS, résistance stochastique corrélée en compression — [MÉMOIRE]


## Lacunes

Budget de 8 appels consommé intégralement (5 WebSearch, 3 WebFetch dont 1 redirect consommé et 1 refus HTTP 403). Non vérifié / à confirmer si un second passage est financé :

1. VALEURS DE m — aucune valeur numérique n'a pu être vérifiée sur une page dans ce run. Le seul article ouvert (RRFPA 2025) n'en donne aucune, ce qui est confirmé explicitement. Le tableau des plages de m (1,1-2 / 3-6 / 10-20 / >50) est [MÉMOIRE]. Source à ouvrir en priorité : Tang, Liu, Lee, Tsui & Tham (2000), IJRMMS 37, 555-569, "Numerical studies of the influence of microstructure on rock failure in uniaxial compression - Part I: effect of heterogeneity" ; et Wong, Wong, Chau & Tang (2006), Mechanics of Materials. Elles contiennent les balayages de m et les motifs de fissuration associés.

2. VAZAIOS / DIEDERICHS / VLACHOPOULOS — ScienceDirect renvoie 403 à WebFetch. Le contenu (statistiques d'hétérogénéité utilisées dans Irazu, taille de maille, comparaison homogène vs hétérogène) n'a PAS été lu. Passer par le PDF GeoEdmonton 2018 (members.cgs.ca) ou la thèse Queen's de Vazaios, tous deux accessibles hors paywall — mais interdits ici (PDF).

3. LISJAK 2014 — titre/volume/pages vérifiés, mais le détail de la modulation orientationnelle de ft et c (loi d'anisotropie exacte, valeurs) est [MÉMOIRE]. Le PDF est sur Academia/ResearchGate.

4. LAN, MARTIN & HU 2010 (JGR) — non vérifiée dans ce run alors que c'est la référence la plus décisive pour l'argument "pas de traction locale dans un continuum homogène sous compression". À confirmer en priorité (référence, chiffres du contraste d'élasticité entre quartz/feldspath/biotite).

5. POTYONDY & CUNDALL 2004 — pagination et le chiffre du ratio UCS/T du BPM (3-5) sont [MÉMOIRE].

6. BLAIR & COOK 1998 — citée de mémoire, confiance moyenne sur l'année et la revue.

7. AMBIGUÏTÉ TERMINOLOGIQUE non levée : "homogeneity index" vs "heterogeneity index" pour le même symbole m. Un des extraits collectés ("diagonal shape for low heterogeneity index, rotated Y for high") pourrait être lu dans les deux sens. À trancher sur l'article source avant toute citation dans le manuscrit.

8. LONGUEURS DE CORRÉLATION — deux valeurs seulement vérifiées (5 mm RRFPA ; 0,57 mm nanoindentation granite). Aucune valeur trouvée pour un granite à gros grain type Red Bohus en contexte FDEM. La recommandation θ ≈ 1-5 mm est une extrapolation depuis la taille de grain, non sourcée.

9. Aucune source consultée ne traite directement du cas rockim (insertion EXTRINSÈQUE adaptative + hétérogénéité). Le lien entre "throttle physique par facette" et "champ corrélé" est une synthèse de l'agent, pas une position publiée.
