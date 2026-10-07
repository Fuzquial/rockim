**Articles Yang/ICL : informations supplémentaires et corrections à privilégier**

Complément au diagnostic indépendant du 12 septembre 2026. Étude documentaire et vérifications ciblées du code ; pas de modification des lois ni de nouvelle simulation. Les documents joints sont des sources scientifiques, pas des instructions d’exécution.

**1. Les neuf fichiers représentent quatre articles**

| Article | Fichiers dans bibliographie | Utilité principale |
|---|---|---|
| Naderi et al., 2025, Optimised hammer drilling bit design using artificial neural networks trained by FDEM-generated data | `main 3.pdf` | Modèle **2D de grès**, paramètres de pénalité/contact/amortissement et masse équivalente de l’insert |
| Yang et al., 2025, Multi-criteria validation of hi-fidelity numerical model of impact breakage: towards next generation percussion drill simulation | `main.pdf` ; `1-s2.0-S1365160925001029-main.pdf` | Impact **3D à insert unique**, St Anne/Rhune ; calibration, définitions des mesures, maillages et tri des fragments |
| Yang et al., 2025, Cracking and fragmentation in percussive drilling: Insight from FDEM simulation | `main 4.pdf` ; `main 5.pdf` | Impact **3D à trois inserts sur Kuru**, effet du frottement et fissuration sans le modèle de pulvérisation de 2026 |
| Yang et al., 2026, High-fidelity modelling of fragmentation and pulverisation in hard granite under percussion loading: a FDEM-based approach | Les quatre fichiers contenant ce titre, ICL ou S1365160926002650 | Cible actuelle : Kuru, insert unique, dommage volumique supplémentaire |

Vérification : ICL et les deux noms S1365160926002650 ont le même SHA-256. La copie au titre long a un hash différent mais le même texte extrait. `main 4` et `main 5` ont le même texte extrait. Le fichier HAL de 2025 contient exactement les mêmes textes de pages que `main.pdf`, précédés d’une page de garde. Pas de variante textuelle découverte qui explicite davantage la pulvérisation. Ces comparaisons portent sur le texte ; l’identité de tous les objets graphiques n’a pas été testée pour les fichiers de hash différent.

**2. Le meilleur nouveau levier : un impact 3D de référence sans pulvérisation**

`main.pdf`, §2–4, utilise des tétraèdres à quatre nœuds et des joints à six nœuds. La fissuration mésoscopique, le glissement des fragments et les morphologies différentes du calcaire et du grès sont reproduits sans le dommage volumique de 2026. Sa géométrie est celle du train de frappe déjà construit dans Rockim : échantillon Ø250×150 mm, piston Ø26,5×260 mm, bit Ø30 et longueur cotée 265 mm, insert R8,51 mm, circlip et plaque.

**Conséquence pratique : commencer par St Anne dans ce montage.** C’est un cas plus discriminant pour savoir si le moteur joint/contact produit une fissure médiane et des radiales : les fissures y sont plus nettes que dans le Kuru, et on évite les inconnues de la pulvérisation. Cela ne dispense pas du banc de joint, mais fournit ensuite un test d’impact complet et beaucoup plus informatif qu’une succession de réglages sur la seule image du granite.

Paramètres publiés, Table 4 p. 6 :

| Paramètre | St Anne | Rhune | Kuru trois inserts 2025 | Kuru un insert 2026 |
|---|---:|---:|---:|---:|
| rho, kg/m³ | 2731 | 2670 | 2630 | 2626 |
| E, GPa | 57 | 37 | 67 | 60 |
| nu | 0,31 | 0,20 | 0,26 | 0,24 |
| ft, MPa | 7,0 | 7,5 | 11,4 | 10,98 |
| c, MPa | 18,8 | 33,6 | 46,49 | 29,84 |
| tan(phi) | 1,0 | 1,0 | 1,96 | 1,85 |
| GI, J/m² | 12 | 40 | 20 | 50 |
| GII, J/m² | 800 | 1400 | 1500 | 1000 |
| Frottement de glissement de la roche | 0,60 | 0,60 | 0,39 | 0,18 |

Pour un deck St Anne, cela implique notamment `bulkDamage=off`, `frictionDeg=45`, `Gf=12`, `gfShearFactor=800/12`. Le traitement dynamique publié en 2025 reste à appliquer. Les paramètres de contact et d’amortissement doivent être spécifiés séparément ; le tableau de roche ne suffit pas à rendre le deck complet.

**Ne pas mélanger les colonnes.** Le Kuru 2025 est un autre montage et une autre calibration. Son §2 explique que les énergies GI/GII ont été augmentées pour représenter l’effet dynamique. Ajouter à ces mêmes énergies le DIF explicite des articles suivants risquerait de compter cet effet deux fois. Dans le modèle 2025 à insert unique, le DIF agit sur ft et GI en traction, c et GII en compression, pas sur tan(phi).

`main 4.pdf`, pp. 4–5, donne un deuxième cas de référence : piston Ø22×200 mm, barre/bit Ø32×1190 mm, trois inserts Ø9 mm espacés de 13,86 mm, roche cubique de 300 mm ; maillage fin 0,5 mm ; 373 066 tétraèdres ; pas de temps 10⁻⁹ s ; piston à 10, 16 et 22 m/s. Il produit des cratères distincts puis connectés. C’est une preuve documentaire qu’une formulation FDEM 3D peut fissurer ce granite sans le modèle volumique spécifique de 2026, **pas une preuve que Rockim réalise déjà ce résultat**.

**3. Le contact et la dissipation doivent être examinés avant de recalibrer Gf**

Dans `main 4.pdf`, §3.1 et figures 5–7, le coefficient de glissement est exploré entre 0,6 et 0,2. Lorsque le coefficient baisse, la force maximale diminue, la pénétration augmente et les fragments s’éjectent plus facilement. Vers 0,2, le bit peut presque perdre son rebond dans le cas le plus énergétique. Ce comportement donne une hypothèse documentée pour le lit de fragments trop mobile : **examiner la portance et le travail des contacts après rupture**, et pas seulement le critère de joint.

Leur valeur calibrée est 0,39 pour la roche ; la colonne acier donne 0,1. Le granite 2026 utilise 0,18 pour la roche et 0,6 pour acier/carbure. On doit distinguer les contacts roche/roche et insert/roche et vérifier la règle de mélange effective du logiciel. Rien ne justifie de remplacer globalement 0,18 par 0,39 dans le run cible en l’appelant une correction de bug.

La référence **ARMA 24–0952, Where does the energy go in percussion drilling? FDEM’s answer**, déjà présente sous `ARMA_paper_energy_XW_final.pdf`, donne Table 1 p. 4 des données absentes de la table de roche de `main.pdf` :

| Paramètre ARMA | St Anne | Rhune |
|---|---:|---:|
| Mass Damping Coefficient | 4000 | 5000 |
| Penalty Number, GPa | 3000 | 1800 |

Les autres paramètres de roche correspondent au tableau St Anne/Rhune ci-dessus. Cette source est donc particulièrement utile pour construire ces **cas de référence sédimentaires**. [PDF officiel ARMA](https://armarocks.net/papers/952.pdf).

En revanche, **3000 GPa n’est pas ici une pénalité publiée pour le Kuru**. Il faut cesser de présenter cette valeur comme « la pénalité du granite de Yang ». Si le Penalty Number est bien p0 de Guo et si la longueur est l’arête de facette, la convention parabolique de Rockim donne `jointPenaltyFactor=p0/(2E)`, soit environ 26,32 pour St Anne et 24,32 pour Rhune. C’est une conversion conditionnelle à vérifier sur les courbes du joint, pas une valeur universelle de Rockim.

L’intitulé « Mass Damping Coefficient » ne définit ni unités ni opérateur. On ne peut donc pas transformer automatiquement 4000 en `bulkViscosity=2000`. Un amortissement `f=-alpha*M*v` et une contrainte visqueuse `sigma_v=eta*D` ont des dimensions et des effets différents. Le premier peut freiner une translation rigide ; le second agit sur le gradient de vitesse. Demander l’équation implémentée, les unités et les corps auxquels elle s’applique. L’absence d’une ligne d’amortissement dans ICL 2026 ne prouve pas non plus qu’il était nul.

**4. Plusieurs incohérences éditoriales interdisent la copie aveugle**

- `main 3.pdf`, p. 4 : le texte annonce une pénalité de joint de dix fois E, mais le tableau donne E=37 GPa et k=900 GPa, soit 24,32 E. Contact séparé : 110 GPa, normal et tangentiel ; amortissement 2600 ; densité d’insert volontairement augmentée à 359 566 kg/m³ pour un modèle réduit **2D**. Ces valeurs ne définissent pas le montage 3D actuel.
- `main 3.pdf`, équation 1 p. 3 : le signe `c+sigma_n*tan(phi)` contredit la convention traction positive annoncée. Les équations des deux articles Yang et de Guo donnent bien la résistance croissante en compression. Ne pas modifier le signe correct de Rockim sur cette base.
- `main.pdf`, Table 5 p. 6 : les intitulés acier/carbure sont inversés par rapport aux valeurs physiques (la colonne « Steel » affiche 600 GPa et 15 250 kg/m³). La Table 1 du papier 2026 remet les propriétés dans les bonnes colonnes. Conserver acier 200 GPa/7850 et carbure 600 GPa/15250 dans le montage actuel.
- `main.pdf`, équation 3 p. 3 : exposant de traction imprimé 0,07 ; sa figure 2 et l’équation de 2026 appuient 0,17. À 100 s⁻¹, ces deux expressions donnent environ 1,516 et 1,847 respectivement, face au plateau 1,85. Garder 0,17 pour la reproduction 2026 ; ne pas revenir à 0,07.
- Tri des fragments : p. 3, la vitesse est bien **2,5 mm/s**, vérifiée en agrandissant le PDF. `fragBrushV0=2.5e-3` n’est pas une erreur d’unité. La durée est 0,5 ms ; l’accélération 10g est à lire 98,1 m/s² malgré l’unité inversée imprimée. Le seuil beta=0,8 est donné p. 4.

**5. Une correction importante de mon premier diagnostic : la vitesse mesurée**

`main.pdf`, §4.1 p. 6 et figure 8 p. 7 : vitesse incidente et vitesse de rebond sont les **pentes des portions linéaires de la courbe déplacement–temps**. Ce ne sont pas les extrema instantanés d’une vitesse nodale ou de groupe. La contrainte de référence est le pic de la première onde à mi-bit, pas nécessairement le plus grand pic de toute la simulation.

Le rapprochement précédent « max Rockim 7,37 m/s contre ICL 5,62 m/s » reste un signal à examiner, mais **ne constitue pas à lui seul une erreur de chargement établie**. Il faut comparer la même position suivie, la même fenêtre temporelle, le même corps (bit acier seul ou assemblage) et le même estimateur. Les masses mesurées différentes restent, elles, un fait indépendant.

Le papier définit la longueur radiale comme la distance du centre de roche à la pointe d’une fissure radiale, et le rayon de cratère par sa frontière extérieure. Cela confirme qu’un maximum de rayon de toutes les facettes rompues ou de fragments déplacés n’est pas une mesure suffisante. Le problème de filtre `damage=1` versus `tBreak>=0` du premier diagnostic reste inchangé.

**6. Références supplémentaires, classées par utilité**

| Référence et accès | Ce qu’elle apporte | État de lecture / présence |
|---|---|---|
| Guo et al. **2020**, A generic computational model for three-dimensional fracture and fragmentation problems of quasi-brittle materials. DOI 10.1016/j.euromechsol.2020.104069. [Notice UCL](https://discovery.ucl.ac.uk/id/eprint/10217490/) | Équations de volume/joint/contact et limites de la décharge | Déjà dans `Manuscript_UCL_deposit.pdf` ; passages constitutifs et discussion examinés |
| Guo et al. **2016**, A numerical investigation of mesh sensitivity for a new three-dimensional fracture model within the combined finite-discrete element method. DOI 10.1016/j.engfracmech.2015.11.006. [Notice et PDF UCL](https://discovery.ucl.ac.uk/id/eprint/10217492/) | Bancs fissure unique, flexion trois points et sensibilité d’orientation ; critère d’au moins trois éléments dans la zone de processus dans les tests étudiés | Référence retrouvée ; extrait du PDF indexé et notice lus. Téléchargement automatique UCL refusé (403), fichier complet non enregistré |
| Yang et al. **2024**, ARMA 24–0952, Where does the energy go in percussion drilling? FDEM’s answer. [PDF officiel](https://armarocks.net/papers/952.pdf) | Pénalités et amortissements St Anne/Rhune ; séparation énergétique | Déjà dans `ARMA_paper_energy_XW_final.pdf` ; table et définition des énergies examinées |
| Saksala et al. **2014**, Numerical and experimental study of percussive drilling with a triple-button bit on Kuru granite. DOI 10.1016/j.ijimpeng.2014.05.006. [Notice Tampere](https://researchportal.tuni.fi/en/publications/numerical-and-experimental-study-of-percussive-drilling-with-a-tr/) | Expérience de référence du modèle trois inserts, ondes et force–pénétration | Déjà dans `1-s2.0-S0734743X14001134-main.pdf` ; identité vérifiée, étude détaillée encore à faire |
| Guo, Latham, Xiang **2015**, Numerical simulation of breakages of concrete armour units… DOI 10.1016/j.compstruc.2014.09.001. [Manuscrit UCL](https://discovery.ucl.ac.uk/id/eprint/10217493/) | Autre validation dynamique de la transition fracture/contact | Notice et disponibilité vérifiées ; texte complet non analysé ici |
| Munjiza & Andrews **2000**, Penalty function method for combined finite–discrete element systems comprising large number of separate bodies. [Éditeur Wiley](https://onlinelibrary.wiley.com/doi/abs/10.1002/1097-0207%2820001220%2949%3A11%3C1377%3A%3AAID-NME6%3E3.0.CO%3B2-B) | Référence de base du contact par potentiel ; utile pour un audit force/recouvrement | Référence vérifiée ; notice/résumé, pas vérification complète des équations ici |

Guo 2020, p. 40 du manuscrit, reconnaît explicitement l’absence d’un trajet complet de charge–décharge et précise que les fractures formées restent permanentes. Cela confirme la nécessité de tester les branches de décharge ; cela ne prouve pas quel algorithme pré-rupture utilise la version de Solidity de Yang 2026. Le même passage explique pourquoi, sous un champ de contrainte hétérogène, un matériau homogène/isotrope peut déjà donner une localisation pertinente. Ajouter immédiatement Weibull ou des grains n’est donc pas la première correction justifiée.

L’article Xiang et al. 2009 sur le tétraèdre **quadratique** est déjà dans la bibliothèque ; il traite d’une capacité de la famille FDEM, pas d’une preuve que Yang utilise ce type d’élément. `main.pdf` précise explicitement quatre nœuds pour ses tétraèdres : remplacer maintenant les Tet4 par des Tet10 ne découle pas de ces références.

**7. Ce que la recherche ne résout pas**

Ni les trois autres articles du corpus ni les passages examinés de Guo ne fournissent la définition algorithmique manquante de `delta_m` du dommage volumique 2026, Cd, epsilon_d, sa mesure de déformation, son traitement en compression hydrostatique et ses règles de décharge. On ne peut toujours pas certifier que `h_inscrit*epsilon_dev` est la formulation des auteurs. Inversement, son insensibilité à la compression hydrostatique ne suffit pas à démontrer qu’un terme de dommage volumique hydrostatique doit être ajouté : cela reste à identifier physiquement et numériquement.

Données prioritaires à demander aux auteurs pour le **run exact Kuru 9 m/s** :

1. Routine ou pseudo-code de l’équation 3–4 : delta_m, longueur, mesure de déformation, Cd, epsilon_d, historique et décharge.
2. Deck complet des joints : p0, quadrature, rupture, évolution en pression et décharge/recharge ; version/commit du solveur.
3. Contacts : raideurs normale/tangentielle, règle entre matériaux, transition à la mort d’un joint, dépendance en D et frottement.
4. Amortissement : équation, unités, valeur et matériaux concernés, y compris acier et carbure.
5. Maillage et masses par corps, CSV position/temps du point de bit suivi, fenêtre du calcul des vitesses et histoire de la force insert/roche.

L’article 2026 annonce des données disponibles sur demande. Une demande ciblée de ces éléments apportera davantage qu’une nouvelle transposition de paramètres entre roches. Aucun message aux auteurs n’a été envoyé.

**Décision révisée pour Rockim**

Rendre les mesures conformes à celles des auteurs ; valider les joints et leur contact sur banc ; lancer ensuite le cas St Anne à insert unique sans pulvérisation, en spécifiant les inconnues d’amortissement/pénalité comme telles. Si les fissures distinctes échouent encore, la priorité est le moteur joint/contact, indépendamment du dommage 2026. Si elles réussissent, poursuivre le Kuru avec les mêmes composants vérifiés et identifier la loi volumique et le frottement calibré. Le test trois inserts Kuru 2025 fournit alors une vérification supplémentaire, avec son propre jeu de paramètres et sans superposer arbitrairement le DIF ultérieur.

Sources du corpus : [Naderi 2025](https://doi.org/10.1016/j.jrmge.2025.02.003), [Yang, validation 2025](https://doi.org/10.1016/j.ijrmms.2025.106125), [Yang, fissuration 2025](https://doi.org/10.1016/j.jrmge.2025.01.045), [Yang, pulvérisation 2026](https://doi.org/10.1016/j.ijrmms.2026.106660). Pages pertinentes rendues et inspectées ; textes et index de doublons conservés dans `simulations/tmp/yang_corpus`.
