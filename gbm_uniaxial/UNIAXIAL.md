# Essai uniaxial : microstructure GBM contre temoin homogene

*Meme maillage (1001 grains de Voronoi, 40 724 tetraedres, eprouvette 16 x 16 x 32 mm), meme chargement, loi elastique. Seul le materiau change : trois phases minerales contre un materiau unique aux moyennes de Voigt. Toute difference vient donc du CONTRASTE DE RAIDEUR et de rien d'autre. Image depouillee : fem3d_0021.vtu.*


> **Reserve emise par le solveur lui-meme.** Le coefficient de Poisson varie de 0,17 a 0,36 entre les phases. Les tetraedres lineaires a un point d'integration se raidissent d'autant plus que ce coefficient est grand : l'artefact est CORRELE A LA PHASE et SOUS-ESTIME le contraste. Les ecarts ci-dessous sont des bornes basses.


> **Seconde reserve.** L'eprouvette est un parallelepipede et non un cylindre, et le maillage intragranulaire est l'eventail par defaut de la tessellation interne. Ces deux points sont a corriger avant tout usage en rapport.


## 1. Module apparent

| grandeur | valeur |
|---|---|
| borne de Voigt, moyenne des modules | 71.21 GPa |
| borne de Reuss, moyenne des souplesses | 66.77 GPa |
| module mesure, GBM | **69.24 GPa** |
| module mesure, temoin homogene | **71.32 GPa** |
| ecart GBM / temoin | **-2.91 %** |

Le temoin doit retomber sur Voigt, puisqu'on lui a donne cette moyenne. Le GBM doit tomber EN DESSOUS : un assemblage reel est toujours plus souple que la moyenne des modules, parce que les grains souples se deforment davantage.


## 2. Concentration de contrainte

Contrainte equivalente de von Mises, ponderee par le volume des elements.

| | mediane | quantile 99 % | quantile 99,9 % | maximum | q99 / mediane |
|---|---|---|---|---|---|
| GBM | 1.3 MPa | 1.8 | 1.9 | 2.1 | **1.33** |
| temoin homogene | 1.3 MPa | 1.4 | 1.4 | 1.4 | **1.04** |
| rapport GBM / temoin | 1.00 | 1.27 | 1.39 | 1.53 | 1.27 |

La derniere colonne est la mesure directe de la concentration : de combien la queue haute depasse le coeur de la distribution. Le temoin homogene donne la dispersion due au seul maillage et aux bords ; tout ce que le GBM a en plus vient du contraste.


## 3. Repartition par phase

| phase | fraction volumique | module | sigma_eq moyen, GBM | idem temoin | rapport |
|---|---|---|---|---|---|
| feldspath | 62.0 % | 70.0 GPa | 1.3 MPa | 1.3 MPa | **0.97** |
| quartz | 31.0 % | 83.1 GPa | 1.5 MPa | 1.3 MPa | **1.16** |
| biotite | 7.0 % | 29.3 GPa | 0.7 MPa | 1.3 MPa | **0.50** |

Un contraste de raideur charge les grains RAIDES et decharge les souples. Si les rapports suivent l'ordre des modules (quartz > feldspath > biotite), l'effet est bien elastique et non un artefact.


## 4. Les maxima se logent-ils aux frontieres de grain ?

La definition naive, *partager une face avec un autre grain*, attrape 91 % du volume quand un grain ne compte qu'une quarantaine de tetraedres : elle ne discrimine rien. On mesure donc la DISTANCE du centre de chaque element au joint de grain le plus proche, rapportee au rayon de grain equivalent.

Rayon de grain equivalent : 1.25 mm. Distance mediane au joint : 0.32 mm, soit 0.26 rayon.

| population | distance mediane au joint, GBM | idem temoin | ecart |
|---|---|---|---|
| les 10 % les plus charges | 0.254 rayon | 0.253 rayon | **+1 %** |
| les 1 % les plus charges | 0.252 rayon | 0.254 rayon | **-1 %** |
| le 0,1 % le plus charge | 0.256 rayon | 0.276 rayon | **-7 %** |

Une distance PLUS COURTE que celle du temoin veut dire que les maxima se rapprochent des joints. Le temoin, materiau uniforme sur le meme maillage, donne la reference : ce qu'on lit chez lui vient du seul maillage.

| tranche de distance au joint | part du volume | sigma_eq moyen, GBM | idem temoin | rapport |
|---|---|---|---|---|
| 0.15 a 0.30 rayon | 86.3 % | 1.33 MPa | 1.33 MPa | **0.998** |
| 0.30 a 0.50 rayon | 10.8 % | 1.32 MPa | 1.33 MPa | **0.989** |
| 0.50 a 1.00 rayon | 2.7 % | 1.31 MPa | 1.34 MPa | **0.980** |
