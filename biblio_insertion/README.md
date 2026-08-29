# Revue bibliographique — insertion cohésive extrinsèque / adaptative

Six volets produits les 2026-08-24/25 pour instruire la question de Fernando :
« le critère d'insertion et la manière d'insérer sont-ils incorrects ? ».
Sources marquées **[VERIFIE]** (contenu lu ou récupéré en ligne) ou
**[MEMOIRE]** (connaissance non revérifiée) — la distinction est portée dans
chaque volet, elle n'est pas décorative.

| fichier | contenu |
|---|---|
| `biblio_volet2_fondations.md` | Camacho & Ortiz 1996, Ortiz-Pandolfi 1999, Pandolfi-Ortiz 2002, Ruiz et al. — la généalogie du schéma extrinsèque : critère effectif, évaluation nodale, absence d'ordonnancement et de plafond, résistances par facette dès 2002 |
| `biblio_impl.md` | Yan, Zheng & Wang IJRMMS 169 (2023) 105439 confronté ligne à ligne à l'implémentation rockim (`insertionSweep`, `activateJoint`, gardes de continuité, pénalité) |
| `biblio_patho.md` | pathologies reconnues : discontinuité temporelle à l'activation (Papoulia et al.), dépendance au maillage du nombre de fragments (Molinari), remèdes publiés (Weibull par facette, Zhou & Molinari 2004) |
| `biblio_paulino.md` | école Paulino/Celes : TopS, PPR, isotropie des chemins (pinwheel, Rimoli & Rojas), biais de maillage |
| `biblio_fdem.md` | les codes FDEM : Munjiza 2004 (intrinsèque), MultiFracS, Y-Geo/Irazu, Fukuda, HOSS ; comparaisons publiées intrinsèque vs extrinsèque |
| `biblio_roches.md` | applications roches : Tang/RFPA (l'hétérogénéité de Weibull nécessaire à la localisation), Lisjak EDZ, BPM, GBM, tunnels profonds |

## Volet Imperial College / Solidity (mission état de l'art, 2026-08-29)

| fichier | contenu |
|---|---|
| `2026-08-29_lot1_bibliographie_imperial.md` | **LOT 1** — bibliographie annotée du FDEM d'Imperial (Y → VGW/VGeST → Solidity) : 20 pièces, généalogie du code, statut d'accès de chacune, **liste de téléchargement ordonnée par valeur**, trous identifiés, et demandes d'extraits de la thèse de Guo. ⚠️ Aucune source n'y est lue en plein texte : le conteneur de la session n'avait pas d'accès sortant au-delà du moteur de recherche. Étiquettes [MÉTA]/[RÉSUMÉ]/[INFÉRÉ]/[SUPPOSÉ] portées partout. |
| `2026-08-29_lot4_bilan_rockim.md` | **LOT 4** — bilan de rockim contre l'état de l'art, avec plan d'action ordonné. Résultat principal : **117 attributions à un code qui n'est pas celui d'Imperial subsistent dans le code, les en-têtes, la suite de vérification et les decks** — table de rachat fournie. Deux blocages seulement : la longueur de référence h (facteur 2,45) et l'injection de 11,1 J par la branche normale du contact. Corrige les recommandations R2 et R4 du lot 3, toutes deux fausses. Treize entrées où rockim devance Imperial. |
| `2026-08-29_lot3_insertion_maillage.md` | **LOT 3** — mémo comparatif des schémas d'insertion et du maillage, avec recommandation motivée (R1 à R6). L'insertion d'Imperial est **intrinsèque**, les joints sont dans **la roche seule**, la souplesse artificielle est **assumée et non corrigée**, et **le maillage adaptatif n'existe pas** — D1 le range en perspective de recherche. Contient la découverte que le **granite Kuru est calibré deux fois différemment** par la même équipe, et que leur coefficient de frottement glissant est un **correctif d'angularité de maillage**, dit par eux. |
| `2026-08-29_lot2c_frottement_tangentiel.md` | **LOT 2c** — l'algorithme tangentiel de frottement, trouvé là où personne ne le cherchait : `f_t = −k_t δ_t − η v_t` plafonné par Coulomb (Xiang, Latham & Farsi 2017, éq. 4-5 p. 4). Il y a bien un ressort tangentiel chez Imperial. Clôt le lot 2. |
| `2026-08-29_lot2b_couplage_endommagement_contact.md` | **LOT 2b** — l'article de pulvérisation 2026 dépouillé : le seul couplage (1−D) publié porte sur la **contrainte d'élément** ; le mot « penalty » n'apparaît pas une fois. Le DIF, avec une **coquille de l'article 2025 démontrée par continuité** (exposant 0,17 et non 0,07). Et l'aveu d'Imperial : le modèle de pulvérisation **n'est pas nécessaire pour St Anne**. |
| `2026-08-29_lot2a_parametres_stanne.md` | **LOT 2a** — les paramètres de St Anne sur sources primaires (ARMA 24-0952 + IJRMMS 191), le bilan d'énergie chiffré, et **DEUX CORRECTIONS** : le frottement 0,18 ne vaut que pour le granite (St Anne est à 0,6, comme le deck) ; la pénalité de 3000 GPa = 52,6 E est réelle, la « correction » du 29/08 qui la réfutait était fausse. |
| `yang2026_pulverisation.md` | l'article IJRMMS 206 (2026) dépouillé sur PDF — source primaire. §5 porte depuis le 29/08 un avertissement : son dernier paragraphe attribuait à Solidity un couplage lu dans `/home/user/solidity`. |
| `guo2014_*.md` | quatre sections de la thèse de Guo (Imperial 2014), lues sur PDF : maillage §2.4, contact et intégration §2.3.4-5, couplage fluide ch. 5, Dolosse ch. 6 |

Conclusions opérationnelles : voir [BILAN_insertion_adaptative.md](../BILAN_insertion_adaptative.md)
§2 (ce que dit la littérature) et §7 (idées écartées avec la raison).

Réserve de méthode : la première exécution de la revue a été fauchée par une
limite de session ; seul le volet 2 en a réchappé. Les cinq autres ont été
relancés avec une discipline de lecture PDF stricte (tranches bornées) — d'où
le nombre de sources marquées [MEMOIRE] dans certains volets, à revérifier
avant toute citation dans le manuscrit.
