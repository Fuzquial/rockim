# Arrêt des fissures sur les plans de litage délaminés : théorie, mesures, observations in situ et solutions pour le code FDEM

*Revue bibliographique et propositions, état des vérifications au 2026-09-01.*

**Convention de citation.** [C] : référence dont l'existence et le contenu cité ont été confirmés par une seconde lecture indépendante du texte intégral. [L] : référence lue par la première passe (texte intégral ou résumé, DOI vérifié) mais non relue par la seconde. [NV] : référence non vérifiée, citée de seconde main, ou dont les valeurs sont de mémoire. Les chiffres marqués « déduit » sont des calculs faits ici à partir des sources, pas des valeurs publiées.

---

## 1. La question et ce que le run montre

Le run 2D (déformations planes, λ = 1, maillage à arêtes alignées sur le litage, paramètres cohésifs de la table 5.1 de Lisjak 2013 [L] : f_t 0,65/0,16 MPa, c 9/1 MPa, G_Ic 7,0/0,4 J/m², G_IIc 35/10 J/m², φ 22°) produit une zone fracturée en losange dont les bords sont des plans de litage délaminés. À l'intérieur du losange, environ 30 % des joints rompus font plus de 60° avec le litage : le code sait donc produire des fissures qui traversent les couches, l'« offre » de traversée existe. Mais aucune fissure ne franchit un plan une fois celui-ci délaminé : la frontière est un plan ouvert, pas un front de fissuration. La question est de savoir si cet arrêt est (a) la physique de l'argilite ou (b) un artefact du rapport G_Ic,min/G_Ic,max = 0,057 (b1), de la 2D (b2), du maillage aligné (b3) ou de λ = 1 (b4).

Traduits dans le vocabulaire des critères de bifurcation, les paramètres du run valent : Γ_i/Γ_b = 0,057 (mode I) et 0,286 (mode II) ; σ_b/σ_i = 4,07 (traction) et 9,0 (cohésion) ; contraste élastique de Dundurs α = 0 (même matériau de part et d'autre du plan). Gupta, Argon et Suo 1992 [L] montrent que l'orthotropie ne modifie pratiquement pas les cartes de délamination, qui ne dépendent que de α : l'isotropie transverse de l'argilite ne change pas les seuils ci-dessous.

---

## 2. La théorie prédit la déviation puis le piégeage avec ces paramètres, à toute incidence

### 2.1 Le critère en ténacité place le run 4 à 15 fois sous le seuil de traversée

He et Hutchinson 1989, IJSS [C], donnent le critère de référence : une fissure arrivant sur une interface dévie si Γ_i/Γ_b < G_d/G_p et traverse sinon (éq. 12), le rapport G_d/G_p étant indépendant de la longueur d'avancée (éq. 11). Valeurs vérifiées dans le texte (β = 0) :

| Incidence sur le plan | Seuil G_d/G_p à α = 0 | Source dans l'article |
|---|---|---|
| 90° | ≈ 1/4 (déviation simple, cas contrôlant) | fig. 3, p. 1056-1057 : « the critical ratio is approximately 1/4 » |
| 60° | ≈ 0,55-0,6 | fig. 11, a/l = 0,1 (lecture ± 0,05) |
| 45° | ≈ 0,7-0,75 | idem |
| 30° | ≈ 0,85-0,9 | idem |

La conclusion de l'article est explicite : pour α entre −0,5 et 0,25, « the toughness of the interface must be less than about one quarter of the toughness of the material across the interface if all cracks are to be deflected », et « the competition becomes more favorable to deflection the more oblique is the crack impinging the interface » (p. 1064). Martínez et Gupta 1994 [NV, cités par Chen 2019] confirment que la déviation double a un rapport plus bas que la déviation simple : c'est bien la simple qui contrôle.

Avec Γ_i/Γ_b = 0,057, le run est 4,4 fois sous le seuil à 90° et environ 15 fois sous le seuil à 30°. Leguillon, Lacroix et Martin 2000 [L] montrent en outre que leur critère asymptotique, qui inclut le décollement en avant de la pointe, est « en général plus favorable à la déviation que celui de He et Hutchinson » et tend vers lui : 1/4 est une borne basse du domaine de déviation. Kendall 1975 [L] obtient, pour une fissure courte de Griffith dans une plaque homogène, R_ad < R_co/[4π(1−ν²)], soit ≈ 0,09 pour ν = 0,3 (déduit) : plus sévère que He-Hutchinson, mais 0,057 reste en dessous. Chandler et al. 2016 [C], appliquant un critère énergétique aux mesures sur shale, situent le seuil à G_c,A/G_c,ST ≈ 3,8, au-dessus duquel « the crack should always deflect into the Short-Transverse orientation regardless of the loading conditions » ; 3,8 correspond à Γ_i/Γ_b ≈ 0,26, cohérent avec le 1/4 de He-Hutchinson (rapprochement fait ici). Le run, à 17,5, est dans ce régime « quel que soit le chargement ».

Deux réserves écrites par He et Hutchinson (p. 1058) s'appliquent au FDEM. La première, « condition (12) implicitly assumes these intrinsic flaws of comparable size », est satisfaite : les joints cohésifs des deux chemins ont la même taille. La seconde, « dynamic effects may alter the conclusions somewhat when the impinging crack is traveling at a significant fraction of the elastic wave speed », vise directement un code explicite, et va dans le sens du piégeage (voir 2.6).

### 2.2 Les critères en résistance placent le run du côté de la déviation, de peu

Cook et Gordon 1964 [C] : pour une pointe de rayon fini, la traction σ_x en avant de la pointe, parallèle au plan de fissure, atteint « approximately constant fraction (~1/5) of the peak stress concentration σ_y,max » ; un plan de faiblesse s'ouvre en avant de la pointe si σ_i < σ_b/5. Avec σ_i/σ_b = 0,246 > 0,2, la pré-délamination à la Cook-Gordon n'est marginalement pas prédite (déduit). Ce ~1/5 vaut pour la pointe elliptique d'Inglis ; les critères en champ singulier sont moins exigeants : Gupta et al. 1992 [L, via Alam 2017] placent la transition à σ_b/σ_i ≈ 3,5, et Parmigiani-Thouless (2.3) à ≈ 3,2. Avec 4,07 le run est au-dessus de ces deux bornes : la résistance seule est déjà du côté de la déviation, mais avec une marge de 15 à 25 % seulement, à comparer au facteur 4,4 en ténacité. La ténacité est le paramètre décisif.

### 2.3 Le critère couplé confirme que la ténacité décide ici

Parmigiani et Thouless 2006 [C] (preprint auteur en libre accès, relu) modélisent les deux chemins par des zones cohésives. Résultats vérifiés dans le texte : (i) « apparent absence of any lower bound for the ratio of the substrate to interface toughness to guarantee crack penetration » ; (ii) « no matter how tough an interface is, crack deflection can always be induced if the strength of the interface is low enough compared to the strength of the substrate » ; (iii) il existe « a lower bound for the ratio of the substrate strength to interfacial strength, below which penetration is guaranteed no matter how brittle the interface », asymptote σ̂_s/σ̂_i ≈ 3,2 (fig. 7). Le régime est fixé par la longueur cohésive ℓ_cz = EΓ/σ̂² rapportée à la géométrie : petite, la ténacité contrôle et l'on retrouve le 1/4 de He-Hutchinson ; grande, la résistance contrôle. Deux compléments importants pour le FDEM : le critère de He-Hutchinson n'est retrouvé exactement que si un branchement préexistant est inclus (appendice), et avec des lois auto-similaires (σ̂_s/σ̂_i = Γ_s/Γ_i) la transition tombe à ≈ 4, coïncidence jugée « probably coincidental » par les auteurs.

Foulk et al. 2008 [L] précisent : si un seul chemin est modélisé, cohésif et LEFM coïncident ; si plusieurs chemins coexistent et se développent simultanément, ce qui est le cas d'un maillage FDEM à joints partout, la zone cohésive qui s'ouvre sur un chemin écrante l'autre, et résistance et ténacité décident ensemble. Pro et al. 2018 [L] donnent ℓ_cz = M·EΓ/σ̂² avec M entre 0,21 et 1 selon la loi. Pour le run (déduit) : massif, E_∥ 3,8 GPa, Γ 7 J/m², σ 0,65 MPa, ℓ_cz ≈ 13-63 mm ; interface, E_⊥ 1,3 GPa, Γ 0,4 J/m², σ 0,16 MPa, ℓ_cz ≈ 4-20 mm. Ces longueurs sont petites devant la galerie, donc régime « ténacité », mais elles doivent être comparées à la taille des éléments : si l'arête dépasse quelques centimètres, la zone cohésive n'est pas résolue et le rapport 17,5 agit comme un rapport résistance × ouverture critique plus que comme un rapport d'énergies (voir 5.2).

L'école Leguillon-Martin ajoute que « a low toughness interface is not systematically a sufficient condition to promote the initiation of deflection » (Martin et al. 2008 [L]), le critère couplé exigeant énergie et contrainte (Leguillon 2002 [L] ; Martin et al. 2001 [L] : forte sensibilité au rapport des extensions, que He-Hutchinson supposent égales). Au run, les deux conditions sont satisfaites.

### 2.4 Les expériences montrent que la LEFM surestime la déviation près du seuil, pas à 0,057

Alam, Parmigiani et Kruzic 2017 [L] (PMMA, σ_b/σ_i = 4,1 et 3,15, Γ_b/Γ_i = 1,68 et 1,42) observent la pénétration à 90° et une transition à 85° et 80°, là où He-Hutchinson prédisaient 57° et 47° et le CZM 80° et 73°. Chen et al. 2019 [L] rapportent, pour Γ_i/Γ_b = 0,25, une transition mesurée entre 60° et 70° (LEFM : 90°), et pour 0,37 aucune déviation à 40° ; « the LEFM criterion might overestimate the tendency of the crack deflection ». Paggi et Reinoso 2017 [L] retrouvent les asymptotes LEFM avec une interface fragile et voient apparaître pénétration et déviation simultanées quand la zone de process s'allonge. Toutes ces corrections jouent près du seuil ; à 0,057 la marge est telle que tous les cadres convergent.

### 2.5 Une fissure entrée dans l'interface y reste

He et Hutchinson 1989, J. Appl. Mech. [C] : le taux de restitution d'un branchement hors de l'interface est accru vers le matériau souple, réduit vers le raide, et « the results suggest a tendency for a crack to be trapped in the interface irrespective of the loading when the compliant material is tough and the stiff material is at least as tough as the interface ». À α = 0, le rapport G_branchement/G_interface est d'ordre 1 ; avec Γ_b = 17,5 Γ_i, la fissure de litage est piégée quel que soit le chargement (déduit). La frontière-plan-délaminé que rien ne franchit est exactement ce régime.

### 2.6 Mode mixte et dynamique n'inversent pas le verdict

La fissure déviée simple est en mode mixte, ψ = tan⁻¹(K₂/K₁) ≈ −40° à α = 0 (fig. 4 de He-Hutchinson [C]) ; en interpolant entre 0,057 (mode I) et 0,286 (mode II) on reste sous 0,25 (déduit). Chalivendra et Rosakis 2008 [L] (Homalite, interface adhésive faible) : à 30° aucune tentative de pénétration même à +50 % de charge ; à 45° la fissure lente (0,35-0,4 c_R) pénètre en plus de dévier, la rapide (0,75-0,8 c_R) « only deflects » ; à 60° pénétration et déviation aux deux vitesses ; la pénétration se fait « approximately at right angle to the inclined interface ». La vitesse favorise donc le piégeage : dans un explicite, un déconfinement trop rapide renforce l'arrêt.

### 2.7 Verdict théorique

Avec Γ_i/Γ_b = 0,057 et σ_b/σ_i = 4,07, l'arrêt sur le plan de litage puis le piégeage sont ce que tous les critères prédisent, à toute incidence et quel que soit le chargement. Ce n'est pas une anomalie numérique. Mais la théorie dit aussi ce que vaut un rapport réaliste (§ 3.1) : entre 0,10 et 0,30 en G. À 0,10 le run reste dans le régime « déviation quel que soit le chargement » (1/0,10 = 10 > 3,8) ; à 0,30 il passe dans le régime dépendant du chargement (1/0,30 = 3,3 < 3,8, et 0,30 > 0,25 à 90°), où la déviation reste prédite à toute incidence oblique (seuils 0,55-0,9) mais devient sensible à K_II/K_I à 90°. En résistance, si l'on prend les rapports éprouvette mesurés (≈ 2-2,3, § 3.2) au lieu des 4,07 du joint, on tombe sous la borne 3,2 de Parmigiani-Thouless où la pénétration est garantie. Un jeu de paramètres réaliste ne supprime donc pas la déviation : il la rend fréquente au lieu de certaine, ce qui est exactement ce que les laboratoires observent (§ 3.3).

---

## 3. La roche : anisotropie de ténacité mesurée et EDZ observées

### 3.1 Le rapport de ténacité mesuré est de 2 à 4 en K, de 3 à 10 en G, jamais 17

Conventions ISRM : Arrester (A) et Divider (D) = fissure en travers du litage ; Short-Transverse (ST) = fissure le long du litage.

| Roche, état | Source | K travers (A / D), MPa·m^0,5 | K long (ST) | Rapport travers/long |
|---|---|---|---|---|
| Mancos, sec, K_Ic^c | Chandler 2016 [C] | 0,65 / 0,72 | 0,21 (lits faibles, 5/7) ; 0,52 (lits forts, 2/7) | 3,1-3,4 ; 1,25-1,4 |
| Mancos, 22 °C après 150 °C | Chandler 2017 [C] | 0,49 / 0,56 (K_Ic) | 0,22 | 2,2-2,5 |
| Nash Point, SCB | Forbes Inskip 2018 [C] | 0,74 / 0,71 | 0,24 | ≈ 3 |
| Anvil Points, 20 et 40 gal/t | Schmidt-Huddle 1977 [C] | 0,9 / 1,05 ; 0,58 / 0,6 | 0,72 ; 0,31-0,40 | ≈ 1,4 ; ≈ 1,7 (lecture d'histogramme) |
| Longmaxi | Luo 2018 [L] | 1,5 (β = 0°) | 0,6 (β = 90°) | 2,1 |
| Mudstone lité, sec / 8,1 % d'eau | Yang 2020 [L] | 0,82-0,85 / 0,45-0,54 | 0,73 / 0,21 | 1,2 / 2,5 |
| Marcellus | Lee 2015 [NV] | 0,73 max | 0,18 min | ≤ 4 |

Chandler 2016 [C] précise que la correction de ductilité m est maximale en ST (1,83), donc que l'anisotropie de K_Ic^c « should potentially be regarded as a minimum ». Yang 2020 [L] montre que l'anisotropie croît avec la teneur en eau (1,18 → 2,5), point pertinent pour une argilite saturée. Forbes Inskip 2018 [C] fournit aussi la fonction angulaire complète : K_Ic = 0,24 (0°), 0,34 (15°), 0,38 (30°), 0,45 (45°), 0,53 (60°), 0,69 (75°), 0,74 (90°), et les résistances en traction ST 3,69, A 7,63, D 8,65 MPa, soit un rapport de résistance travers/long de 2,1-2,3. Synthèse : K_travers/K_long ≈ 3 (± 1) pour les shales les plus anisotropes, 1,2-2 pour les autres ; aucune mesure ne dépasse ≈ 4.

Le passage à G, la grandeur que le code utilise, dépend du traitement de l'anisotropie élastique :

| Conversion | G travers / G long | Γ_i/Γ_b |
|---|---|---|
| K ≈ 3-3,4, E isotrope (Chandler [C] : E = 35,65 GPa) | 9,58 (lits faibles), 1,56 (lits forts) | 0,10 ; 0,64 |
| E mesuré en flexion par orientation (Chandler [C] : 21/11/8-12 GPa) : G_c 27 (D), 38 (A), 6 (ST faible), 19 (ST fort) J/m² | 4,5-6,3 ; travail de rupture R_SR : 6,5-9 | 0,16-0,22 |
| Nash Point [C], (0,24/0,74)², E isotrope | 9,5 | 0,105 |
| K ≈ 3 converti avec E_∥/E_⊥ ≈ 3 (Opalinus : Lisjak 3,8/1,3 [L] ; Liu 2024 [L] : « up to three times ») | ≈ 3 | ≈ 0,3 |
| Effet d'échelle direct, Marcellus (Li, Jin, Cusatis 2019 [L], arXiv) : G_f 29,0 (A), 37,9 (D), 44,8 (ST) N/m | 0,65-0,85, inversé | > 1 |
| Table 5.1 de Lisjak [L] : 7,0 / 0,4 J/m² | 17,5 | 0,057 |

Le 17,5 vaut donc 1,8 fois le plus grand rapport mesuré (9,6, E isotrope), 2 à 4 fois les rapports avec E anisotrope ou en travail de rupture (4,5-9), et il est à l'opposé de la seule mesure directe d'énergie (Cusatis, < 1). La fourchette réaliste pour Γ_i/Γ_b est 0,10-0,30.

### 3.2 Le 17,5 de Lisjak est une constante de calibration, pas une mesure

La thèse (p. 124 [L]) écrit : « the values of the strength parameters were obtained as final result of the calibration process », les cibles étant les essais de Bock 2009 [NV, cité par Lisjak] : UCS 11,6 (P), 14,9 (S), 4,1 MPa (45°) ; brésilien 1,30 (P) et 0,67 MPa (S). Le rapport de traction éprouvette reproduit vaut donc 0,52, comparable au 0,43-0,48 de Nash Point [C], alors que le rapport imposé au joint est 0,246, deux fois plus sévère. L'explication vraisemblable (déduit, chapitre 5 non relu) est la dilution par le chemin en zigzag : dans l'approche « smeared » (Lisjak, Grasselli, Vietor 2014 [L]), la fissure macroscopique en travers du litage emprunte des joints obliques dont les propriétés sont interpolées, et l'anisotropie au niveau du joint doit être exagérée pour que l'anisotropie de l'éprouvette sorte juste. Le 17,5 hérite de cette exagération et dépend en principe du maillage qui l'a produit (h = 0,30 mm dans la calibration, une autre taille dans la galerie).

L'interpolation elle-même mérite un contrôle. Avec une interpolation linéaire en G entre 0,057 et 1, les joints à 30°, 45° et 60° du litage reçoivent 0,37, 0,53 et 0,69 de G_max ; les mesures de Nash Point [C] converties en G (E isotrope) donnent 0,26, 0,37 et 0,51. Les chemins obliques du run sont donc relativement plus tenaces que dans la roche, mais le rapport interface/chemin oblique reste très en dessous des seuils obliques de He-Hutchinson dans les deux cas (0,083 contre 0,55-0,6 à 60° pour le run ; 0,21 pour Nash Point).

Aucune valeur de K_Ic par orientation n'a pu être extraite pour l'Opalinus (Valente 2012 [L] : SCB parallèle au litage seulement ; Nejati 2020 [L] : valeurs non accessibles) ni pour le Callovo-Oxfordien (thèse Abdulmajid 2020 [NV], PDF verrouillé). C'est un trou réel : la fourchette 0,10-0,30 vient des shales, pas de l'argilite cible.

### 3.3 Dans les galeries, les fissures d'extension ne traversent pas le litage ; la zone, elle, le traverse

**Bure, Callovo-Oxfordien** (litage quasi horizontal, σ_H N150-155°E, σ_v ≈ σ_h, σ_H/σ_h ≈ 1,3 à −490 m ; Wileveau 2007 [L], Djizanne 2019 [L], Armand 2013 [L]). Extensions de la zone fracturée mesurées, Table 2 d'Armand et al. 2013 [L] reprenant Armand et al. 2014 [L], en diamètres (min/moy/max) :

| Galerie | Contraintes dans la section | Voûte : extension / cisaillement | Paroi : extension / cisaillement | Radier : extension / cisaillement |
|---|---|---|---|---|
| // σ_H | σ_v ≈ σ_h, λ ≈ 1 | 0,1-0,15 D / aucun | 0,01-0,4 D / 0,7-1,0 D | 0,1-0,15 D / aucun |
| // σ_h | σ_H = 1,3 σ_v | 0,2-0,4 D / 0,5-0,8 D | 0,1-0,2 D / aucun | 0,2-0,5 D / 0,8-1,1 D |

Dans la galerie // σ_H, dont la section est le cas le plus proche de λ = 1 avec litage horizontal, la zone est allongée le long du litage avec un rapport paroi/voûte de 5 à 8 pour la zone cisaillée. Une anisotropie de contrainte de 30 % suffit à basculer la zone à la verticale, en travers du litage : « the in situ horizontal stress is about 1.3 times the vertical one… main reason of the anisotropic damage pattern » ; « in all the drifts, the major convergence is measured where the fracture zone is located ». Mécanisme : 75 % des fractures induites sont en mode II/III, 25 % en mode I ; « spalling is not the prevailing mechanism… shear failure seems to occur first from the excavation front face ». Les fractures de cisaillement (chevrons), inclinées sur un litage sub-horizontal, le traversent (déduit) ; elles sont profondes et peu transmissives, les fractures d'extension transmissives restant près de la paroi (Armand 2014 [L]). Souley, Vu, Armand 2022 [L] reproduisent ces formes avec des plans de faiblesse ubiquistes plus une matrice élastoplastique en 3D sous champ de contraintes anisotrope. Point capital pour un run 2D : les fractures qui traversent le litage à Bure naissent au front, en cisaillement, c'est-à-dire par un effet 3D absent d'une section en déformations planes.

**Mont Terri, Opalinus** (litage ≈ 45° SE ; σ_1 ≈ 7 MPa sub-vertical ; E ⊥ 4 GPa contre 10 GPa // ; Bossart 2002 [L], Yong 2010 [L] ; les valeurs de σ_2 et σ_3 de Martin et Lanyon 2003 [NV] n'ont pas été relues). Pour les ouvrages parallèles à la direction du litage :

- Labiouse et Vietor 2014 [L], cylindres creux déchargés : carottes parallèles au litage, « cracks sub-parallel to the bedding planes open and lead to a buckling failure in two regions that extend from the borehole in the direction normal to bedding » ; carottes perpendiculaires, aucune rupture ; « striking similarity » avec l'overcarottage in situ ; sens inverse de l'argile de Boom (cisaillement conjugué).
- Kupferschmied et al. 2015 [L] et Amann et al. 2017 [L], forage BHM-3 : à court terme, cisaillement tangentiel le long du litage, wing cracks et horsetails, réseau limité à un quart de diamètre ; à 6 jours, réseau « chimney-like » d'au moins un diamètre développé « in both lateral and radial extent perpendicular to the bedding plane orientation » ; à 30 jours, plus de deux diamètres de plaques flambées ; « buckling is associated with the formation of extensional fractures normal to bedding in the center and lateral to the buckling zone » ; la dissipation des surpressions interstitielles est le moteur.
- Nussbaum et al. 2011 [L] : « extensional EDZ fractures often do not propagate through pre-existing natural discontinuities such as well developed inclined bedding planes and/or fault planes » ; « whatever their orientation, pre-existing faults spatially limited the propagation of EDZ fractures ».
- Yong et al. 2017 [L], EZ-B : « the physical manifestation of the relatively weak bedding plane strength is dominated by bedding-perpendicular displacement as opposed to bedding-parallel shear ».
- Marschall et al. 2017 [L], HG-A : breakouts « mainly at locations where bedding is oriented tangential to the tunnel circumference » ; ouverture des plans de litage, spalling, flambage.

Le mécanisme de progression en travers du litage n'est donc pas la pénétration d'une fissure isolée : c'est une succession de délaminations plus profondes, couche après couche, dont les plaques flambent et se fissurent en extension normalement au litage à l'intérieur des plaques. Ce processus demande une face libre, des plaques élancées, et du temps (jours à semaines). À l'instant de l'excavation, l'observation (BHM-3 à quelques heures) est une zone petite et parallèle au litage.

**Tournemire** (litage sub-horizontal, σ_v 3,8, σ_h 2,1, σ_H 4,0 MPa ; Matray et al. 2007 [L], Rejeb et Cabrera 2006 [L]) : EDZ = « unloading joints, mimicking the drift shape » plus fissures de désaturation parallèles au litage ; extension 0,16-0,22 R, indépendante de l'âge ; 30 cm en forage parallèle au litage, 50 cm en forage incliné à 45°.

### 3.4 Réponse à la question posée

Une fissure d'extension isolée arrivant sur un plan de litage bien développé s'y arrête ou s'y dévie : c'est ce que disent les laboratoires (Chandler [C] : la moitié des brésiliens Arrester rejetés pour déviation dans le litage, les short-rod Arrester rompant par traction transversale et « become trapped there » ; Schmidt-Huddle [C] : « load increases occurred when the crack was seen to begin propagating in the bedding planes causing crack tip blunting and delamination » ; Forbes Inskip [C] : les fissures obliques ou normales « commonly deflect toward the weaker Short-transverse orientation » ; Nejati 2020 [L], 124 essais dont Opalinus : les fissures « kink towards the bedding planes ») et ce que dit le terrain (Nussbaum 2011). L'hypothèse (a) est fondée pour la fissure. Mais in situ la ZONE traverse le litage sous trois conditions absentes du run : une anisotropie de contrainte (Bure, 1,3 suffit), un front d'excavation 3D (chevrons en cisaillement), et une cinématique de flambage de plaques avec le temps (Mont Terri). Le losange borné est cohérent avec l'observation à t = 0 d'une galerie à λ ≈ 1 ; il n'est pas cohérent avec l'EDZ mature.

---

## 4. Ce que les codes font pour obtenir des fissures traversantes

**Lisjak et al. 2015** [L] (Y-Geo, tunnel FE de Mont Terri, élasticité transversalement isotrope, résistance des éléments cohésifs fonction de l'orientation, calibration sur les convergences) décrivent leur mécanisme : « failure initiates due to shearing of bedding planes critically oriented with respect to the compressive circumferential stress… Slippage-induced rock mass deconfinement then promotes extensional fracturing in the direction perpendicular to the bedding orientation » ; « the shape and extent of the EDZ around tunnels and boreholes are strongly dependent upon the relative orientation between excavation axis and bedding planes ». La thèse [L] conclut de même : « delamination along bedding planes and subsequent extensional fracturing as key mechanisms of the damage process potentially leading to buckling and spalling phenomena ». Les fissures « en travers » de Lisjak ne sont donc pas des fissures qui pénètrent une interface intacte : ce sont des fissures d'extension nées dans les plaques délaminées, alimentées par le déconfinement que crée le glissement du litage aux orientations critiques. Deux ingrédients les rendent possibles : un litage incliné par rapport à la contrainte tangentielle (à 45° dans le tunnel FE), et un rapport de ténacité en mode II (0,286) et de cohésion (0,111) qui rend le glissement du litage facile. Avec un litage horizontal et λ = 1, les orientations critiques de cisaillement existent aux épaules de la galerie, mais leur mobilisation dépend du chemin de déconfinement.

**Lisjak et al. 2016** [L] (HG-A) montrent que « fracture termination is simulated at the intersection with a pre-existing discontinuity » : leur code, comme le vôtre, arrête les fissures sur les discontinuités, et ils le présentent comme un résultat conforme aux observations.

**Souley, Vu, Armand 2022** [L] obtiennent les formes de Bure avec des joints ubiquistes et une matrice élastoplastique sous champ anisotrope en 3D : c'est le champ de contraintes et la 3D qui font traverser, pas un critère d'interface.

**Renshaw et Pollard 1995** [L] fournissent le critère de traversée d'une interface frictionnelle non liée, validé sur six roches : la traversée exige que « the compressive stress required to prevent slip along the interface at the moment when the stress on the opposite side of the interface is sufficient to initiate a fracture » soit disponible. La forme usuelle −σ'_xx/(T₀−σ'_yy) > (0,35+0,35/μ)/1,06 est de mémoire et non vérifiée ici [NV]. Sarmadivaleh et Rasouli 2014 [L] l'étendent à une interface cohésive non orthogonale ; Zeng et Wei 2017 [NV] ajoutent frottement et anisotropie de contrainte. Dans un FDEM à contact frottant, ce mécanisme est natif : une fissure ré-amorce de l'autre côté d'un plan délaminé si le plan est fermé et bloqué en frottement, de sorte que la traction en avant de la pointe se transmet ; il ne peut rien se passer si le plan est ouvert ou glisse. Le critère n'est donc pas à implémenter mais à vérifier dans le run (état de contact des plans-frontières à l'arrivée des fissures).

**Paggi et Reinoso 2017** [L], **Chen et al. 2019** [L], **Aranda et al. 2020** [L] montrent que champ de phase et zones cohésives retrouvent les asymptotes LEFM et produisent pénétration et déviation simultanées quand la zone de process grandit. Aucun article de cette revue ne traite du contournement 3D d'un plan délaminé de taille finie : c'est un trou de la littérature, pas seulement de cette recherche.

---

## 5. Solutions classées

Chaque solution est décrite par ce qu'elle change, pourquoi, son coût en unités de « un run 2D de la configuration actuelle », et ce qu'elle mesure. Grandeurs à extraire systématiquement : aire et élancement du losange (grand axe/petit axe, orientation), extension maximale normale et parallèle au litage en diamètres, fraction de joints rompus à plus de 60° du litage, nombre d'événements de traversée (fissure continue de part et d'autre d'un plan délaminé), état de contact des plans-frontières (ouvert, fermé-glissant, fermé-bloqué), vitesse de pointe maximale.

### 5.1 Paramétrique immédiat

**S1. Rapport G_Ic,min/G_Ic,max de 0,057 à 0,10, 0,25, 0,30 et 0,50, résistances inchangées.** Pourquoi : 0,10 est la valeur mesurée à E isotrope (Chandler [C], Nash Point [C]) ; 0,25 est le seuil de He-Hutchinson à 90° [C] ; 0,30 est la valeur K ≈ 3 convertie avec E_∥/E_⊥ = 3 ; 0,50 est au-dessus du seuil à 90° mais sous les seuils obliques. Coût : 4 runs. Mesure : si le losange persiste à 0,30 et 0,50, l'arrêt ne tient pas au 0,057 ; s'il se dissout entre 0,10 et 0,30, (b1) est confirmée ; 0,25 doit être le point de bascule à incidence normale si le code respecte la LEFM.

**S2. Découpler résistance et ténacité.** Deux variantes : f_t 0,246 conservé avec G à 0,30 (test « ténacité seule ») ; G 0,057 conservé avec f_t à 0,45 (rapport éprouvette mesuré 0,43-0,52 [C]) puis 0,31 (borne 1/3,2 de Parmigiani-Thouless [C]). Pourquoi : Parmigiani-Thouless montrent que sous σ_b/σ_i ≈ 3,2 la pénétration est garantie quelle que soit la ténacité ; Foulk 2008 [L] que dans un système multi-chemins les deux comptent. Coût : 3 runs. Mesure : nombre de traversées ; si la variante f_t = 0,31 traverse avec G = 0,057, le code est dans le régime « résistance » et le rapport d'énergie est secondaire, ce qui signale une zone cohésive non résolue (voir S5).

**S3. λ = 1,3 et λ = 1/1,3.** Pourquoi : à Bure, 1,3 renverse l'orientation de la zone ; c'est le levier in situ le plus fort, et il active le mécanisme de Lisjak (cisaillement du litage aux orientations critiques, déconfinement, extension normale au litage). Coût : 2 runs. Mesure : orientation du losange ; apparition de fissures d'extension normales au litage aux épaules ; comparaison aux extensions de la Table 2 d'Armand (0,15 D en voûte contre 0,7-1,0 D en paroi pour λ ≈ 1 ; 0,5-1,1 D en voûte/radier pour 1,3). C'est aussi la seule solution qui rapproche le run d'une donnée in situ chiffrée.

### 5.2 Numérique

**S4. Témoin isotrope, rapports = 1 sur le même maillage aligné.** Pourquoi : si le losange persiste sans anisotropie, il est un artefact du maillage ; s'il disparaît, il est porté par les paramètres. Coût : 1 run. C'est le contrôle falsifiant de (b3) et il doit passer avant tout le reste.

**S5. Maillage à deux familles.** Rangées de joints alignés sur le litage à l'espacement réel des plans (les 0,6 m du descriptif, ou l'espacement stratigraphique cible), triangles non structurés et aléatoires entre les rangées avec propriétés isotropes (pas d'interpolation), et distribution bimodale des rangées à la Chandler [C] : cinq septièmes à 1/9,6, deux septièmes à 1/1,56. Variante : maillage tourné de 10-15° par rapport au litage, anisotropie portée uniquement par la fonction d'interpolation. Pourquoi : dans le run, chaque arête alignée est « le lit le plus faible » ; la roche a des lits faibles et des lits forts, et des fissures d'extension qui traversent les lits forts dans les plaques (Kupferschmied [L]). Coût : 2 runs plus le maillage. Mesure : le losange doit devenir irrégulier, borné par les rangées faibles, avec des traversées dans les rangées fortes ; si le maillage tourné reproduit un losange en zigzag, l'alignement n'est pas la cause.

**S6. Contrôle de résolution de la zone cohésive.** Calculer ℓ_cz = EΓ/σ̂² (13-63 mm massif, 4-20 mm interface, § 2.3) et le comparer à l'arête du maillage de la galerie. Si l'arête dépasse ℓ_cz, raffiner localement ou rescaler G pour que l'énergie dissipée par élément corresponde ; sinon le rapport 17,5 n'est pas un rapport d'énergies au sens des critères (Pro 2018 [L] : la LEFM ne s'applique que si ℓ_cz est petite devant la structure ET résolue). Coût : un calcul, puis éventuellement 1 run raffiné.

**S7. Chemin de chargement et vitesse.** Relever la vitesse de pointe maximale et la comparer à c_R ; ralentir le déconfinement (relaxation progressive du noyau) si les fissures dépassent 0,3-0,4 c_R. Pourquoi : He-Hutchinson [C] réservent leur critère aux fissures lentes ; Chalivendra-Rosakis [L] montrent que la vitesse favorise le piégeage. Coût : 1 run. Mesure : nombre de traversées à vitesse réduite.

### 5.3 Modèle

**S8. Diagnostic Renshaw-Pollard sur les plans-frontières.** Extraire, aux instants d'arrivée des fissures, σ_n et τ sur le plan délaminé et l'état de contact. Pourquoi : si le plan est ouvert ou glisse (τ > σ_n tan 22°), aucune traction ne se transmet et l'arrêt est physique au sens de Renshaw-Pollard [L] ; s'il est fermé et bloqué et que la fissure s'arrête quand même, c'est le rapport d'énergies qui bloque. Coût : post-traitement seulement. Ce diagnostic tranche entre « frontière mécanique » et « frontière énergétique » sans changer le code.

**S9. Frottement et reprise de cohésion partielle sur les joints délaminés (clé opt-in, défaut intact).** Deux tests de sensibilité : φ résiduel augmenté (aucune valeur vérifiée pour le litage d'Opalinus dans cette revue, à documenter), et reprise d'une fraction de f_t et de G après fermeture. Pourquoi : la reprise de cohésion n'a de base physique qu'à long terme et par voie hydraulique (auto-colmatage, Blümling 2007 [L]) ; c'est un test de sensibilité, pas un modèle. Coût : développement plus 2 runs. Mesure : à partir de quelle fraction de cohésion reprise les traversées apparaissent.

**S10. Mécanisme de flambage.** Vérifier si des plaques délaminées existent près de la paroi là où le litage est tangentiel et si elles portent des fissures d'extension normales au litage. Pourquoi : c'est le mécanisme par lequel la zone traverse le litage à Mont Terri (Labiouse-Vietor [L], Kupferschmied [L]) et dans les runs de Lisjak [L]. Coût : post-traitement, puis éventuellement S3 et S7 qui le favorisent. Ce qu'on ne peut pas attendre d'un FDEM mécanique sec : la progression sur des jours pilotée par les pressions interstitielles.

**S11. 3D limitée.** Une tranche mince avec un plan délaminé de taille finie (patch) pour tester le contournement ; le front d'excavation, source des chevrons de Bure, est hors de portée. Coût : élevé. Mesure : une fissure contourne-t-elle un patch délaminé de largeur finie ; aucune référence de cette revue ne répond à cette question, ce qui en fait un résultat publiable en soi.

### 5.4 Protocole de validation et arbre de décision

Ordre recommandé : S4 (témoin), S8 (diagnostic), S6 (résolution), puis S1, S2, S3, puis S5, S7, enfin S9-S11.

| Observation | Conclusion |
|---|---|
| Losange persiste à rapports = 1 (S4) | Artefact de maillage (b3) ; passer à S5 avant tout |
| Losange disparaît à rapports = 1, persiste à 0,30 et 0,50 (S1) et sur maillage tourné (S5) | Physique (a) au sens des critères : la frontière est robuste au paramètre |
| Losange se dissout entre 0,10 et 0,30 (S1) | Calibration (b1) : 0,057 exagère un comportement réel |
| Traversées apparaissent avec f_t = 0,31 à G = 0,057 (S2) | Régime « résistance » ; vérifier S6, le rapport d'énergie n'est pas ce que le code voit |
| Plans-frontières ouverts ou glissants à l'arrivée des fissures (S8) | Arrêt mécanique conforme à Renshaw-Pollard ; seule la contrainte normale (λ, profondeur) peut le lever |
| λ = 1,3 réoriente la zone et crée des fissures normales au litage aux épaules (S3) | (b4) : le run à λ = 1 est correct pour une galerie // σ_H, incomplet pour les autres ; comparer aux extensions d'Armand |
| Losange persiste à λ = 1,3, maillage aléatoire, rapport 0,30 | Chercher dans la 2D et l'absence de front (b2), seule hypothèse restante |

Le résultat qui tranche entre « physique » et « artefact » est la combinaison S4 + S1 à 0,30 + S5 tourné : si le losange survit aux trois, c'est la physique des critères ; s'il tombe à l'un d'eux, c'est l'ingrédient correspondant. Le résultat qui tranche entre « run juste » et « EDZ juste » est S3 comparé à la Table 2 d'Armand.

---

## 6. Bibliographie

### Confirmées (existence et contenu relus par la seconde passe)

- He M.-Y., Hutchinson J.W. (1989). Crack deflection at an interface between dissimilar elastic materials. Int. J. Solids Structures 25(9), 1053-1067. DOI 10.1016/0020-7683(89)90021-8.
- He M.-Y., Hutchinson J.W. (1989). Kinking of a crack out of an interface. J. Appl. Mech. 56(2), 270-278. DOI 10.1115/1.3176078.
- Cook J., Gordon J.E. (1964). A mechanism for the control of crack propagation in all-brittle systems. Proc. R. Soc. Lond. A 282(1391), 508-520. DOI 10.1098/rspa.1964.0248.
- Parmigiani J.P., Thouless M.D. (2006). The roles of toughness and cohesive strength on crack deflection at interfaces. J. Mech. Phys. Solids 54(2), 266-287. DOI 10.1016/j.jmps.2005.09.002.
- Chandler M.R., Meredith P.G., Brantut N., Crawford B.R. (2016). Fracture toughness anisotropy in shale. J. Geophys. Res. Solid Earth 121, 1706-1729. DOI 10.1002/2015JB012756.
- Chandler M.R., Meredith P.G., Brantut N., Crawford B.R. (2017). Effect of temperature on the fracture toughness of anisotropic shale and other rocks. Geol. Soc. London Spec. Publ. 454, 295-303. DOI 10.1144/SP454.6.
- Schmidt R.A., Huddle C.W. (1977). Fracture mechanics of oil shale, some preliminary results. Sandia Laboratories SAND 76-0727, OSTI 7119762.
- Forbes Inskip N.D., Meredith P.G., Chandler M.R., Gudmundsson A. (2018). Fracture properties of Nash Point shale as a function of orientation to bedding. J. Geophys. Res. Solid Earth 123, 8428-8444. DOI 10.1029/2018JB015943.

### Lues par la première passe (DOI vérifié), non relues par la seconde

Théorie : Alam M., Parmigiani J.P., Kruzic J.J. (2017) Eng. Fract. Mech. 181, 116-129, DOI 10.1016/j.engfracmech.2017.05.013. Gupta V., Argon A.S., Suo Z. (1992) J. Appl. Mech. 59(2S), S79-S87, DOI 10.1115/1.2899511. Kendall K. (1975) Proc. R. Soc. Lond. A 344, 287-302, DOI 10.1098/rspa.1975.0102. Leguillon D., Lacroix C., Martin E. (2000) J. Mech. Phys. Solids 48(10), 2137-2161, DOI 10.1016/S0022-5096(99)00101-5. Martin E., Leguillon D., Lacroix C. (2001) Compos. Sci. Technol. 61(12), 1671-1679, DOI 10.1016/S0266-3538(01)00067-7. Leguillon D. (2002) Eur. J. Mech. A/Solids 21(1), 61-72, DOI 10.1016/S0997-7538(01)01184-6. Martin E., Leguillon D. (2004) Int. J. Solids Struct. 41(24-25), 6937-6948, DOI 10.1016/j.ijsolstr.2004.05.044. Martin E., Poitou B., Leguillon D., Gatt J.M. (2008) Int. J. Fract. 151(2), 247-268, DOI 10.1007/s10704-008-9228-0. Foulk J.W., Johnson G.C., Klein P.A., Ritchie R.O. (2008) J. Mech. Phys. Solids 56(6), 2381-2400, DOI 10.1016/j.jmps.2007.12.006. Pro J.W., Sehr S., Lim R.K., Petzold L.R., Begley M.R. (2018) J. Mech. Phys. Solids 121, 480-495, DOI 10.1016/j.jmps.2018.08.015. Chen H. et al. (2019) Comput. Methods Appl. Mech. Eng. 347, 1085-1104, DOI 10.1016/j.cma.2019.01.014. Chalivendra V.B., Rosakis A.J. (2008) Eng. Fract. Mech. 75(8), 2385-2397, DOI 10.1016/j.engfracmech.2007.08.005. Renshaw C.E., Pollard D.D. (1995) Int. J. Rock Mech. Min. Sci. 32(3), 237-249, DOI 10.1016/0148-9062(94)00037-4. Sarmadivaleh M., Rasouli V. (2014) Rock Mech. Rock Eng. 47, 2107-2115, DOI 10.1007/s00603-013-0509-1. Kim J.-B., Shin H., Lee W., Rhee K.Y. (2008) Arch. Appl. Mech. 78, 811-819, DOI 10.1007/s00419-007-0195-0. Paggi M., Reinoso J. (2017) Comput. Methods Appl. Mech. Eng. 321, 145-172, DOI 10.1016/j.cma.2017.04.004. Aranda M.T. et al. (2020) Theor. Appl. Fract. Mech. 105, 102389, DOI 10.1016/j.tafmec.2019.102389.

Ténacité : Forbes Inskip N.D., Meredith P.G. (2021) Rock Mech. Rock Eng., DOI 10.1007/s00603-021-02403-4. Li W., Jin Z., Cusatis G. (2019) Rock Mech. Rock Eng., arXiv 1710.07366 (DOI de la revue non vérifié). Luo Y. et al. (2018) Sci. Rep. 8, DOI 10.1038/s41598-018-26846-y. Yang J., Li L., Lian H. (2020) PLoS ONE 15(8), e0237909, DOI 10.1371/journal.pone.0237909. Zhou Q. et al. (2023) Rock Mech. Rock Eng., DOI 10.1007/s00603-023-03454-5. Suo Y. et al. (2020) Rock Mech. Rock Eng., DOI 10.1007/s00603-020-02131-1. Valente S. et al. (2012) Rock Mech. Rock Eng. 45, 767-779, DOI 10.1007/s00603-012-0225-2. Nejati M. et al. (2020) Int. J. Solids Struct. 195, 74-90, DOI 10.1016/j.ijsolstr.2020.03.004. Gehne S. et al. (2020) J. Geophys. Res. Solid Earth 125, DOI 10.1029/2019JB018971. Lisjak A. (2013) PhD thesis, University of Toronto (TSpace). Lisjak A., Grasselli G., Vietor T. (2014) Int. J. Rock Mech. Min. Sci. 65, 96-115, DOI 10.1016/j.ijrmms.2013.10.006. Liu L. et al. (2024) Rock Mech. Rock Eng., DOI 10.1007/s00603-024-04030-1.

EDZ : Armand G. et al. (2014) Rock Mech. Rock Eng. 47(1), 21-41, DOI 10.1007/s00603-012-0339-6. Armand G., Noiret A., Zghondi J., Seyedi D.M. (2013) J. Rock Mech. Geotech. Eng. 5(3), 221-230, DOI 10.1016/j.jrmge.2013.05.005. Djizanne H. et al. (2019) Geomech. Energy Environ. 17, 3-15, DOI 10.1016/j.gete.2018.11.003. Wileveau Y. et al. (2007) Phys. Chem. Earth 32, 866-878, DOI 10.1016/j.pce.2006.03.018. Souley M., Vu M.-N., Armand G. (2022) Rock Mech. Rock Eng. 55, 4183-4207, DOI 10.1007/s00603-022-02841-8. Bossart P. et al. (2002) Eng. Geol. 66, 19-38, DOI 10.1016/S0013-7952(01)00140-5. Nussbaum C. et al. (2011) Swiss J. Geosci. 104, 187-210, DOI 10.1007/s00015-011-0070-4. Amann F. et al. (2017) Swiss J. Geosci. 110, 151-171, DOI 10.1007/s00015-016-0245-0. Marschall P. et al. (2017) Swiss J. Geosci. 110, 173-194, DOI 10.1007/s00015-016-0246-z. Yong S. et al. (2017) Rock Mech. Rock Eng. 50, 1959-1985, DOI 10.1007/s00603-017-1212-4. Yong S., Kaiser P.K., Loew S. (2010) Int. J. Rock Mech. Min. Sci. 47, 894-907, DOI 10.1016/j.ijrmms.2010.05.009. Labiouse V., Vietor T. (2014) Rock Mech. Rock Eng. 47, 57-70, DOI 10.1007/s00603-013-0389-4. Kupferschmied N. et al. (2015) Int. J. Rock Mech. Min. Sci. 77, 105-114, DOI 10.1016/j.ijrmms.2015.03.027. Lisjak A. et al. (2015) Tunn. Undergr. Space Technol. 45, 227-248, DOI 10.1016/j.tust.2014.09.014. Lisjak A. et al. (2016) Rock Mech. Rock Eng. 49, 1849-1873, DOI 10.1007/s00603-015-0847-2. Matray J.-M., Savoye S., Cabrera J. (2007) Eng. Geol. 90, 1-16, DOI 10.1016/j.enggeo.2006.09.021. Rejeb A., Cabrera J. (2006) Eurosafe 2006, INIS x1r0h-ph217. Blümling P. et al. (2007) Phys. Chem. Earth 32, 588-599, DOI 10.1016/j.pce.2006.04.034. Hale S. et al. (2021) Solid Earth 12, 1581-1600, DOI 10.5194/se-12-1581-2021.

### Non vérifiées, déduites ou citées de seconde main

- Martínez D., Gupta V. (1994) J. Mech. Phys. Solids 42(8), 1247-1271, DOI 10.1016/0022-5096(94)90034-5 : non lu, contenu via Chen 2019.
- He M.-Y., Evans A.G., Hutchinson J.W. (1994) Int. J. Solids Struct. 31(24), 3443-3455, DOI 10.1016/0020-7683(94)90025-6 : métadonnées seulement.
- Hutchinson J.W., Suo Z. (1992) Adv. Appl. Mech. 29, 63-191, DOI 10.1016/S0065-2156(08)70164-9 : accès refusé.
- Xu L.R., Huang Y.Y., Rosakis A.J. (2003) J. Mech. Phys. Solids 51(3), 461-486, DOI 10.1016/S0022-5096(02)00080-7 : via Alam 2017.
- Zeng X., Wei Y. (2017) J. Mech. Phys. Solids 101, 235-249, DOI 10.1016/j.jmps.2016.12.012 : DOI vérifié, contenu via extraits.
- Lee H.P. et al. (2015) J. Geophys. Res. Solid Earth 120, 169-181, DOI 10.1002/2014JB011358 : valeurs via Chandler 2016 et Li-Jin-Cusatis.
- Abdulmajid M. (2020) Thèse Sorbonne Université, theses.fr 2020SORUS064 : résumé seulement, PDF verrouillé ; seule source identifiée de K_Ic par orientation pour le COx.
- Martin C.D., Lanyon G.W. (2003) Int. J. Rock Mech. Min. Sci. 40, 1077-1088, DOI 10.1016/S1365-1609(03)00113-8 : existence vérifiée, valeurs de contraintes de mémoire.
- Bock H. (2009) : cibles de calibration citées par Lisjak 2013, non consultées.
- Mirkhalaf et al., Kadin et al. : cités par Chen 2019, non consultés.
- Formule de Renshaw-Pollard à constantes 0,35/1,06 : de mémoire, non vérifiée.
- Extrait de moteur de recherche attribué au corps d'Armand 2014 (chevrons « modérément inclinés », rôle des « quasi-horizontal bedding planes ») : non vérifié dans le texte.