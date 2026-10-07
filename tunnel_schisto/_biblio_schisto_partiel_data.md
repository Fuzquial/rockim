# Angle « None » — 34 references (NON VERIFIEES : la verification adverse n a pas tourne)

METHODE ET RESERVES. Le budget WebSearch de la session etait deja epuise ; j'ai donc interroge directement les API OpenAlex, Crossref, Semantic Scholar, Europe PMC et OSTI, puis lu les resumes (jamais le texte integral) via ces API ou via les pages editeur quand elles etaient accessibles. Tous les DOI ci-dessous ont ete verifies par Crossref/OpenAlex. « lu » = resume lu ; « deduit » = titre/metadonnees verifies mais contenu infere ; aucune reference n'est inventee. Les valeurs chiffrees sont uniquement celles presentes dans les resumes : pour calibrer, il faudra ouvrir les PDF (la plupart sont en acces libre).

CE QUE DISENT LES DONNEES 2018-2026 (roches metamorphiques foliees uniquement, shales exclus).
1) Courbe en U/V quasi universelle en compression (UCS, E, cohesion, resistance a long terme, resistance dynamique) : minimum a 30-60 deg selon la roche. Ardoise (Alejano 2021, Yang 2024, Li 2021, Zhu 2024/2025) : rupture par le plan de clivage sur 15-75 deg, minimum 45-60 deg. Phyllite (Xu 2021, Wen 2022, Wu 2024, Deng 2025, Liu 2025) : minimum 30-60 deg, en W a sec et en U sature (Wu 2024). Schiste (Condon 2020, Liu 2026, Peng 2024, Yan 2020) : minimum 45-60 deg ; le schiste chloriteux est peu anisotrope. Gneiss (Acosta-Violay 2020, Young 2020, Gupta 2025) : minimum 30-45 deg, parabole du second ordre. Le confinement attenue l'anisotropie non lineairement jusqu'a une pression critique (Liu 2026 ; Wu 2024 : +20 a +100 % de resistance).
2) Elasticite : Condon 2020 est LA reference pour calibrer une elasticite TI realiste : E statique de 21 a 117 GPa selon l'orientation dans le schiste de Poorman, avec un module de cisaillement G13 anormalement bas attribue au glissement sur la foliation des qu'elle est cisaillee ; l'anisotropie statique est bien plus forte que l'anisotropie dynamique (Thomsen eps=0,133, gamma=0,119). Alejano 2021 et Li 2021 donnent les 5 constantes TI de l'ardoise (seule E parallele au plan est dependante de la taille).
3) Traction et tenacite : BTS et traction directe croissent de 0 a 90 deg (en U pour la BTS « apparente », en S pour la traction directe : Liu 2022, Ou 2020, Liu 2020). K_Ic de l'ardoise en V (Li 2023, K_Ic = 0,094 sigma_t + 0,036 en geometrie divider). Schiste : K_divider > K_arrester > K_short-transverse, propagation intergranulaire en short-transverse et transgranulaire sinon (Jahnke 2022). Gneiss migmatitique : anisotropie de K_Ic par SCB (Bercakova 2017, sans chiffres dans le resume). Granodiorite foliee du Grimsel : les fissures de mode I deviennent (kink) vers la foliation, critere MTS + T-stress le meilleur ; FPZ proportionnelle a (K_Ic/sigma_t)^2 (Nejati 2020, Dutler 2018). L'anisotropie s'efface en dynamique : coefficient N = 3,25 (statique) contre 1,35 (SHPB) sur l'ardoise (Liu 2020).
4) Eau : UCS de l'ardoise saturee -24,3 %, fragilite -17,1 %, glissement precoce sur les feuillets (Zheng 2022) ; l'eau favorise le cisaillement sur les plans et reduit la fragilite (Yang 2024) ; phyllites tres sensibles (Wen 2022, Deng 2025, Zhou 2021), l'eau attaquant biotite et argiles.
5) Mecanismes observes : glissement frictionnel inter-feuillets a 45 deg et flambage des feuillets verticaux a 90 deg (Deng 2025) ; fendage axial a 0 deg passant a un kinking plastique sous fort confinement (Liu 2026) ; microfissures de traction dominantes (Wen 2022) ; fissures initiees dans les zones a micas (Young 2020) ; trajet de fissure lisse, intergranulaire et plus rapide (x2) dans le gneiss que dans le granite (Jacobsson 2024).
6) Parametres de joint pour un modele cohesif : Garcia-Fernandez 2018 deduit c et phi du plan de schistosite a partir de bresiliens a angle variable, valide par cisaillement direct ; Alejano 2021 et Deng 2025 fournissent des criteres de Jaeger ameliores (plage d'angles ou le plan gouverne).

RIEN TROUVE / LACUNES HONNETES : aucun article experimental recent dedie aux kink bands dans les roches foliees (seules des mentions dans Liu 2026 et Deng 2025 ; les references restent Donath 1961 et Anderson 1974) ; aucune courbe K_Ic(angle) chiffree pour un gneiss ou un micaschiste depuis 2018 en dehors de Jahnke 2022 (schiste) ; pas de donnees BTS/K_Ic recentes pour micaschistes stricto sensu ; les resumes de Xu 2018 (phyllite triaxiale, reference la plus citee), Zhu 2025 (ardoise carbonee) et Gonzalez-Fernandez 2025 n'ont pas pu etre lus (paywall ou champ absent des API). Pour situer, les fondateurs pre-2016 non listes ici : Cho et al. 2012 (schiste de Yeoncheon, gneiss d'Asan) et Tan et al. 2015 (bresiliens sur gneiss de Freiberg et ardoise de Moselle).

POUR LE CODE FDEM : calibrer (i) l'elasticite TI avec un G13 faible (Condon 2020), (ii) la resistance cohesive des joints de foliation via c-phi de Garcia-Fernandez 2018 ou le critere de Jaeger d'Alejano 2021, (iii) la tenacite avec le rapport divider/arrester/short-transverse de Jahnke 2022 et la relation K_Ic-sigma_t de Li 2023, puis verifier que le modele reproduit la courbe en U et la transition de modes (fendage 0 deg, glissement 30-60 deg, traversant 90 deg). Liu 2022 est le seul FDEM experimental-numerique recent sur ardoise : il n'aligne pas le maillage mais introduit un critere de traction anisotrope phenomenologique dans les elements cohesifs.

- **Alejano L.R., Gonzalez-Fernandez M.A., Estevez-Ventosa X., Song F., Delgado-Martin J., Munoz-Ibanez A., Gonzalez-Molano N., Alvarellos J. (2021) Anisotropic deformability and strength of slate from NW-Spain. Int. J. Rock Mech. Min. Sci. 148:104923** (2021) — https://doi.org/10.1016/j.ijrmms.2021.104923
  - roche : Ardoise de toiture, Galice (NW Espagne) ; methode : Compression uniaxiale et triaxiale a plusieurs pendages du clivage, jauges pour les 5 constantes elastiques TI
  - anisotropie : Elasticite transversalement isotrope (5 constantes) + critere composite : plan de faiblesse de Jaeger la ou le clivage gouverne, Hoek-Brown pour la matrice
  - resultat : Rupture par le plan de clivage observee pour les pendages de 15 a 75 deg ; criteres distincts normal/parallele a la foliation ameliorant la prediction de la resistance
  - limite : Quasi-statique seulement, pas de tenacite ; heterogeneite de l'ardoise traitee dans Gonzalez-Fernandez 2025 ; valeurs chiffrees dans le corps du texte, pas dans le resume ; acces libre : oui ; confiance : lu

- **Gonzalez-Fernandez M.A., Perez-Rey I., Song F., Muralha J., Day J.J., Giacomini A., Alejano L.R. (2025, en ligne ; vol. 2026) Understanding the anisotropic stress-strain behavior of heterogeneous slate in uniaxial compressive strength testing. J. Rock Mech. Geotech. Eng.** (2025) — https://doi.org/10.1016/j.jrmge.2025.07.032
  - roche : Ardoise heterogene (meme groupe que Alejano 2021, Galice) ; methode : Essais UCS avec analyse des courbes contrainte-deformation selon l'orientation
  - anisotropie : Interpretation des courbes contrainte-deformation en fonction de l'angle et de l'heterogeneite (deduit du titre)
  - resultat : Prolonge Alejano 2021 sur la forme des courbes contrainte-deformation ; contenu non lu
  - limite : Resume non accessible via les API ; a lire directement (acces libre) ; acces libre : oui ; confiance : deduit

- **Li K., Yin Z.-Y., Han D., Fan X., Cao R., Lin H. (2021) Size Effect and Anisotropy in a Transversely Isotropic Rock Under Compressive Conditions. Rock Mech. Rock Eng. 54:4639-4662** (2021) — https://doi.org/10.1007/s00603-021-02558-0
  - roche : Ardoise (origine non precisee dans le resume) ; methode : UCS, triaxial et resistance residuelle sur plusieurs diametres et orientations beta
  - anisotropie : 5 constantes elastiques TI mesurees (seule E parallele au plan isotrope depend de la taille) + criteres de rupture dependant taille/orientation/confinement
  - resultat : Relation resistance-taille-orientation-confinement capturee par un critere unique ; resistance residuelle anisotrope sans effet de taille
  - limite : Criteres empiriques ; origine de l'ardoise non indiquee ; pas de valeurs dans le resume ; acces libre : inconnu ; confiance : lu

- **Li K., Cheng Y., Yin Z.-Y., Han D., Meng J. (2020) Size Effects in a Transversely Isotropic Rock Under Brazilian Tests: Laboratory Testing. Rock Mech. Rock Eng. 53:2623-2642** (2020) — https://doi.org/10.1007/s00603-020-02058-7
  - roche : Meme ardoise que Li 2021 ; methode : Bresiliens a plusieurs diametres et angles de foliation
  - anisotropie : BTS en fonction de l'angle et de la taille (deduit du titre)
  - resultat : Compagnon du precedent cote traction indirecte ; resume non lu
  - limite : Contenu deduit du titre ; acces libre : inconnu ; confiance : deduit

- **Liu P., Liu Q., Huang X., Hu M., Bo Y., Yuan D., Xie X. (2022) Direct Tensile Test and FDEM Numerical Study on Anisotropic Tensile Strength of Kangding Slate. Rock Mech. Rock Eng. 55:7765-7789** (2022) — https://doi.org/10.1007/s00603-022-03036-x
  - roche : Ardoise de Kangding, tunnel du chemin de fer Sichuan-Tibet ; methode : Traction directe a 5 pendages avec emission acoustique + FDEM (elements finis-discrets combines)
  - anisotropie : Critere phenomenologique de rupture en traction anisotrope (coefficient d'anisotropie, base Nova-Zaninetti) implante dans les elements cohesifs FDEM ; valide aussi sur bresiliens de 6 autres roches litees
  - resultat : Resistance en traction directe croissant en S avec l'angle de foliation ; 3 modes : activation des feuillets, mixte, matrice ; EA en 3 phases (calme, sauts, explosion)
  - limite : Details (elasticite TI ? maillage aligne ?) non visibles dans le resume ; traction seulement ; pas d'acces libre ; acces libre : non ; confiance : lu

- **Yang X., Li J., Zhang Y., Lei J., Li X., Huang X., Xu C. (2024) Experimental Study on Mechanical Properties of Anisotropic Slate under Different Water Contents. Applied Sciences 14(4):1473** (2024) — https://doi.org/10.3390/app14041473
  - roche : Ardoise de Changsha (Chine) ; methode : XRD + UCS + triaxial a 0, 30, 45, 60, 90 deg et 3 etats hydriques (sec, naturel, sature)
  - anisotropie : Hoek-Brown modifie pour la resistance, loi de Hooke generalisee (TI) pour E
  - resultat : UCS, E et cohesion en U avec l'angle ; nu et phi peu sensibles a l'angle et a l'eau ; l'eau favorise le cisaillement sur les plans et reduit la fragilite
  - limite : Pas de valeurs chiffrees dans le resume ; revue MDPI ; acces libre : oui ; confiance : lu

- **Zheng L., Xie H., Xu Z., Deng J., Wang D., Zhang G., Li C., Zhang R., Feng G. (2022) A Comparison of Mechanical Properties and Failure Processes of Saturated and Unsaturated Slate from Sichuan-Tibet Plateau Area, China. Lithosphere 2022:4503366** (2022) — https://doi.org/10.2113/2022/4503366
  - roche : Ardoise, roche encaissante d'un tunnel du plateau Sichuan-Tibet ; methode : UCS avec emission acoustique, sature vs non sature
  - anisotropie : Pas de balayage d'angle ; effet de la foliation lu via une chute precoce de sigma1 (glissement sur les feuillets)
  - resultat : UCS saturee -24,3 % ; indice de fragilite -17,1 % ; EA diffuse et de faible energie en sature, concentree et violente en non sature
  - limite : Une seule orientation ; angle non precise ; acces libre : oui ; confiance : lu

- **Ou X., Zhang X., Feng H., Zhang C., Zhou X., Wang L. (2020) Static and Dynamic Brazilian Tests on Layered Slate considering the Bedding Directivity. Adv. Civ. Eng. 2020:8860558** (2020) — https://doi.org/10.1155/2020/8860558
  - roche : Ardoise litee a anisotropie moyenne (Chine) ; methode : Bresiliens statiques et SHPB a 5 pendages, camera rapide
  - anisotropie : Courbe BTS apparente vs angle
  - resultat : BTS statique et dynamique en U ; renforcement dynamique maximal a 45 deg ; le mode de rupture depend du degre d'anisotropie
  - limite : BTS apparente (rupture non centrale) ; origine non precisee ; acces libre : oui ; confiance : lu

- **Liu Y., He C., Wang S., Peng Y., Lei Y. (2020) Dynamic Splitting Tensile Properties and Failure Mechanism of Layered Slate. Adv. Civ. Eng. 2020:1073608** (2020) — https://doi.org/10.1155/2020/1073608
  - roche : Ardoise litee (Chine) ; methode : SHPB bresilien a 7 angles, camera rapide, comparaison statique
  - anisotropie : Coefficient d'anisotropie N et loi charge-vitesse par angle
  - resultat : Charge de rupture croissante de 0 a 90 deg ; N = 3,25 en statique contre 1,35 en dynamique ; rupture le long des feuillets pour <= 45 deg, centrale pour > 60 deg
  - limite : Charge de rupture plutot que contrainte ; origine non precisee ; acces libre : oui ; confiance : lu

- **Li E., Wei Y., Chen Z., Zhang L. (2023) Experimental and Numerical Investigations of Fracture Behavior for Transversely Isotropic Slate Using Semi-Circular Bend Method. Applied Sciences 13(4):2418** (2023) — https://doi.org/10.3390/app13042418
  - roche : Ardoise litee ; methode : Bresiliens + SCB a angles 0-90 deg (geometrie divider) + methode hybride elements finis-cohesifs (FCEM)
  - anisotropie : Elements cohesifs avec proprietes distinctes selon le litage (FCEM) ; mode II explore numeriquement
  - resultat : sigma_t et K_Ic en V de 0 a 90 deg ; relation K_Ic = 0,094 sigma_t + 0,036 ; bons accords experience-simulation
  - limite : K_IIc uniquement numerique ; origine non precisee ; acces libre : oui ; confiance : lu

- **Zhu Y., Wang X., Liu B., Liu X., Xue H. (2025) Effects of foliation angle on mechanical characteristics of carbonaceous slate using uniaxial compression tests. J. Rock Mech. Geotech. Eng. 17(4)** (2025) — https://doi.org/10.1016/j.jrmge.2024.08.012
  - roche : Ardoise carbonee, tunnels de l'ouest de la Chine ; methode : UCS a differents angles de foliation
  - anisotropie : Courbe UCS/E vs angle (deduit)
  - resultat : Resume non accessible ; l'article compagnon Zhu et al. 2024 (Applied Sciences 15(1):236, doi 10.3390/app15010236, lu) montre une resistance a long terme en U avec l'angle et un modele de fluage fractionnaire
  - limite : Contenu deduit ; a lire (acces libre) ; acces libre : oui ; confiance : deduit

- **Wang Z., Feng G., Liu X., Zhou Y. (2023) An Experimental Investigation on the Foliation Strike-Angle Effect of Layered Hard Rock under Engineering Triaxial Stress Path. Materials 16(17):5987** (2023) — https://doi.org/10.3390/ma16175987
  - roche : Ardoise mince a pendage raide ; methode : Vrai triaxial (sigma1 croissante, sigma2 constante, sigma3 decroissante) a 2 angles de direction (strike) avec EA
  - anisotropie : Effet de l'angle de direction de la foliation (et non du pendage)
  - resultat : Resistances similaires a 0 et 90 deg en vrai triaxial classique, mais capacite portante plus faible a 0 deg sous chemin de deconfinement ; fractures de surface uniquement le long des feuillets
  - limite : Petits echantillons ; 2 orientations ; acces libre : oui ; confiance : lu

- **Garcia-Fernandez C.C., Gonzalez-Nicieza C., Alvarez-Fernandez M.I., Gutierrez-Moizant R.A. (2019, en ligne 2018) New methodology for estimating the shear strength of layering in slate by using the Brazilian test. Bull. Eng. Geol. Environ. 78:2283-2297** (2018) — https://doi.org/10.1007/s10064-018-1297-3
  - roche : Ardoise du NW de l'Espagne ; methode : Bresiliens a angle variable + champ de contrainte analytique + Mohr-Coulomb sur le plan de schistosite ; validation par cisaillement direct
  - anisotropie : Cohesion et angle de frottement du plan de schistosite deduits des seuils de rupture des bresiliens
  - resultat : c et phi du plan de foliation obtenus par bresilien tres proches de ceux du cisaillement direct
  - limite : Valeurs numeriques absentes du resume ; acces libre : inconnu ; confiance : lu

- **Xu G., He C., Su A., Chen Z. (2018) Experimental investigation of the anisotropic mechanical behavior of phyllite under triaxial compression. Int. J. Rock Mech. Min. Sci. 104:100-112** (2018) — https://doi.org/10.1016/j.ijrmms.2018.02.017
  - roche : Phyllite (tunnel, Chine) ; methode : Triaxial a plusieurs angles de foliation et confinements
  - anisotropie : Courbes resistance/E vs angle et confinement, modes de rupture (deduit)
  - resultat : Reference la plus citee sur la phyllite ; contenu non verifie au-dela du titre
  - limite : Resume non accessible via les API ; paywall ; acces libre : non ; confiance : deduit

- **Xu J., Fei D., Yu Y., Cui Y., Yan C., Bao H., Lan H. (2021) Research on crack evolution law and macroscopic failure mode of joint Phyllite under uniaxial compression. Sci. Rep. 11** (2021) — https://doi.org/10.1038/s41598-021-83571-9
  - roche : Phyllite jointee (foliation traitee comme joint) ; methode : UCS a plusieurs inclinaisons + simulation numerique
  - anisotropie : Coefficient d'effet de joint ; seuils sigma_ci et sigma_cd
  - resultat : UCS en U, minimum a 60 deg ; rupture traction+cisaillement le long du joint pour 30-75 deg ; sigma_ci = 0,30-0,59 sigma_f, sigma_cd = 0,44-0,86 sigma_f ; sigma_cd a 90 deg = borne fiable de l'UCS
  - limite : Origine non precisee ; joint assimile a la foliation ; acces libre : oui ; confiance : lu

- **Wen G., Hu J., Wu Y., Zhang Z.-X., Xu X., Xiang R. (2022) Mechanical Properties and Failure Behavior of Dry and Water-Saturated Foliated Phyllite under Uniaxial Compression. Materials 15(24):8962** (2022) — https://doi.org/10.3390/ma15248962
  - roche : Phyllite foliee ; methode : UCS sec et sature a plusieurs angles, EA + DIC + MEB
  - anisotropie : Courbes pic/deformation vs angle, 4 modes de rupture par plage d'angle
  - resultat : Pic en U ou V avec l'angle ; l'eau affaiblit fortement ; microfissures de traction dominantes ; l'eau favorise la microfissuration dans biotite et argiles
  - limite : Pas de chiffres dans le resume ; acces libre : oui ; confiance : lu

- **Wu Y., Yang D. (2024) Study on mechanical characteristics of phyllite by conventional triaxial compression under multi-factors coupling action. Sci. Rep. 14** (2024) — https://doi.org/10.1038/s41598-024-81802-3
  - roche : Phyllite (roche tendre metamorphique de tunnel) ; methode : Triaxial conventionnel croisant confinement, etat hydrique et angle de litage
  - anisotropie : Courbes pic/E vs angle selon l'etat hydrique
  - resultat : Pic et E en W (sec) et en U (sature) ; confinement : +20,6 a +55,4 % (sature) et +23,3 a +102,2 % (sec) sur le pic ; modes : glissement-cisaillement, compression-traction-cisaillement, traction, composite
  - limite : Origine non precisee ; acces libre : oui ; confiance : lu

- **Deng T., Liu B., Shi X., Chu Z., Zhang X., Yu M. (2025) Anisotropic Mechanical Properties of Sericite Phyllite Under Dry and Saturated Conditions. Rock Mech. Rock Eng.** (2025) — https://doi.org/10.1007/s00603-025-04467-y
  - roche : Phyllite sericiteuse ; methode : Triaxial a 0, 45 et 90 deg, sec et sature
  - anisotropie : Modele de module lie a la foliation (FRM, fonction angle et confinement) + critere de resistance par morceaux EPSC (Jaeger ameliore)
  - resultat : Resistance et deformation en U ; a 45 deg glissement frictionnel inter-feuillets avec chute brutale du module ; a 90 deg flambage des feuillets verticaux ; deformation en A pour les echantillons satures confines
  - limite : 3 angles seulement ; paywall ; acces libre : non ; confiance : lu

- **Liu Q., Zhang C., Hou J., Li L. (2025) Experimental and numerical research on the anisotropic compressive mechanical behavior of layered rock. Sci. Rep. 15** (2025) — https://doi.org/10.1038/s41598-025-32026-6
  - roche : Phyllite naturelle ; methode : Triaxial en chargement et dechargement a 0, 30, 45, 90 deg + PFC2D
  - anisotropie : Contacts DEM distincts pour la matrice et le plan de litage ; le mode composite depend du rapport resistance matrice/plan
  - resultat : E et resistance decroissent puis croissent avec l'angle ; cohesion residuelle toujours inferieure au pic ; PFC2D coherent avec l'experience
  - limite : 4 angles ; DEM et non FDEM ; acces libre : oui ; confiance : lu

- **Zhou Y., Su S., Li P. (2021) Mechanical Behavior, Energy Release, and Crack Distribution Characteristics of Water-Saturated Phyllite under Triaxial Cyclic Loading. Adv. Civ. Eng. 2021:3681439** (2021) — https://doi.org/10.1155/2021/3681439
  - roche : Phyllite ; methode : Triaxial cyclique MTS 815 + EA, sec et sature
  - anisotropie : Pas de balayage d'angle indique ; effet eau
  - resultat : Phyllite tres sensible a l'eau ; variable d'endommagement energetique critique environ 0,80
  - limite : Angle de foliation non renseigne dans le resume ; acces libre : oui ; confiance : lu

- **Condon K.J., Sone H., Wang H.F., EGS Collab Team (2020) Low Static Shear Modulus Along Foliation and Its Influence on the Elastic and Strength Anisotropy of Poorman Schist Rocks, Homestake Mine, South Dakota. Rock Mech. Rock Eng. 53** (2020) — https://doi.org/10.1007/s00603-020-02182-4
  - roche : Schiste de Poorman (Homestake, site EGS Collab, Dakota du Sud) ; methode : UCS a plusieurs orientations de l'axe de symetrie, modules statiques et dynamiques
  - anisotropie : Elasticite TI complete ; anisotropie statique dominee par un G13 anormalement bas
  - resultat : E statique de 21 a 117 GPa selon l'orientation ; UCS de 21,9 a 194,6 MPa, minimum a 45-60 deg ; G13 bas attribue au glissement sur la foliation des qu'elle est cisaillee ; Thomsen dynamiques eps=0,133, gamma=0,119 bien plus faibles
  - limite : Roche heterogene de site, forte dispersion ; acces libre : oui ; confiance : lu

- **Jahnke B., Ruplinger C., Bate C.E., Trzeciak M., Sone H., Wang H.F. (2022) Fracture toughness of schist, amphibolite, and rhyolite from the Sanford Underground Research Facility (SURF), Lead, South Dakota. Sci. Rep. 12** (2022) — https://doi.org/10.1038/s41598-022-20031-y
  - roche : Schiste de Poorman (plus amphibolite et rhyolite) ; methode : CCNBD (disque bresilien a chevron) en 3 geometries : divider, arrester, short transverse ; effet de la vitesse de chargement ; microscopie
  - anisotropie : K_Ic selon l'orientation de la fissure par rapport a la foliation
  - resultat : K_divider > K_arrester > K_short-transverse ; propagation intergranulaire en short transverse, transgranulaire sinon ; K_Ic diminue quand la vitesse de chargement diminue ; tortuosite dependante de l'orientation
  - limite : Valeurs chiffrees absentes du resume ; acces libre : oui ; confiance : lu

- **Liu Q., Xiang J., Yin X., Song K. (2026) Effects of confining pressure and loading direction on the mechanical behavior of schist with high content and aggregation degree of mica. PLoS ONE 21** (2026) — https://doi.org/10.1371/journal.pone.0344580
  - roche : Micaschiste a forte teneur et agregation de micas (tunnel) ; methode : Micro-essais + compression a sigma3 de 0 a 20 MPa et plusieurs angles alpha
  - anisotropie : Degre d'anisotropie fonction du confinement ; confinement critique de transition anisotrope-isotrope
  - resultat : Pic en epaule et E en U avec alpha ; le confinement attenue non lineairement l'anisotropie ; a 90 deg traction+cisaillement puis cisaillement pur ; a 0 deg fendage puis kinking plastique ; a 30 deg glissement ; fissures directionnelles le long des clivages de mica
  - limite : Microstructure specifique (micas agreges), comportement differe des schistes pauvres en mica ; acces libre : oui ; confiance : lu

- **Peng Y., Du Z., Chen P., Yao Y., Liu G., Wu L. (2024) Study on Dynamic Mechanical Properties and Failure Pattern of Thin-Layered Schist. Applied Sciences 14(19):9101** (2024) — https://doi.org/10.3390/app14199101
  - roche : Schiste mince du groupe de Wudang (Chine) ; methode : SHPB en compression a plusieurs angles et espacements de schistosite + simulation
  - anisotropie : Angle et espacement des plans de schistosite
  - resultat : Resistance dynamique en U ; 0 et 90 deg : fendage axial ou selon la schistosite ; 30-60 deg : cisaillement ; espacement 22 -> 7 mm : -20,3 % de resistance
  - limite : Dynamique uniquement ; acces libre : oui ; confiance : lu

- **Yan S., Wang Q., Wang H., Qiu S., Zeng Z., Fang Y. (2020) Strength Control Factors of Chlorite Schist under Schistose Structure. Int. J. Design & Nature and Ecodynamics 15(5)** (2020) — https://doi.org/10.18280/ijdne.150503
  - roche : Schiste chloriteux (tunnel) ; methode : UCS a angles 0-90 deg, etats hydriques et d'alteration
  - anisotropie : Courbe UCS vs angle
  - resultat : UCS en V avec l'angle ; anisotropie faible ; eau et alteration degradent la structure schisteuse
  - limite : Revue peu connue ; pas de chiffres dans le resume ; acces libre : oui ; confiance : lu

- **Ozbek A., Gul M., Karacan E., Alca O. (2018) Anisotropy effect on strengths of metamorphic rocks. J. Rock Mech. Geotech. Eng. 10(1):164-175** (2018) — https://doi.org/10.1016/j.jrmge.2017.09.006
  - roche : Phyllite, schiste, gneiss et calcschiste du massif du Menderes (Turquie) ; methode : Marteau de Schmidt et essais de laboratoire perpendiculairement et parallelement a la foliation
  - anisotropie : Deux orientations seulement (perpendiculaire vs parallele)
  - resultat : Resistances perpendiculaires generalement superieures aux paralleles ; gneiss et calcschiste forts, phyllite et schiste poreux et faibles
  - limite : Pas de courbe complete en angle ; acces libre : oui ; confiance : lu

- **Acosta M., Violay M. (2020) Mechanical and hydraulic transport properties of transverse-isotropic Gneiss deformed under deep reservoir stress and pressure conditions. Int. J. Rock Mech. Min. Sci. 130:104235** (2020) — https://doi.org/10.1016/j.ijrmms.2020.104235
  - roche : Gneiss de Cresciano (Suisse) ; methode : Triaxial en conditions de reservoir profond a 0, 30, 45, 60, 90 deg, porosite, EA, permeabilite et Vp
  - anisotropie : Deux classes de proprietes : mecaniques en U, transport en decroissance monotone
  - resultat : Dilatance et pic en U (max a 0 et 90 deg, min a 30-45 deg) ; permeabilite et Vp maximales a 0 deg et minimales a 90 deg
  - limite : Conditions geothermiques specifiques ; acces libre : oui ; confiance : lu

- **Young R.P., Nasseri M.H.B., Sehizadeh M. (2020) Mechanical and seismic anisotropy of rocks from the ONKALO underground rock characterization facility. Int. J. Rock Mech. Min. Sci. 126:104190** (2020) — https://doi.org/10.1016/j.ijrmms.2019.104190
  - roche : Gneiss veine et pegmatite granitique d'Olkiluoto (Finlande) ; methode : Vrai triaxial sur chemin de contrainte hybride, EA, tomographie X, lames minces
  - anisotropie : Angle de foliation ; fissures initiees dans les zones a micas
  - resultat : Gneiss veine a 30 deg le plus faible : rupture en cisaillement a 55 MPa, coherente avec le terrain ; la pegmatite emet sans fracture macroscopique
  - limite : Peu d'orientations ; acces libre : oui ; confiance : lu

- **Gupta A., Panthee S., Bhandari A., Shrestha S. (2025) Analysis of strength anisotropy and correlation between UCS and point load test in augen gneiss at varying anisotropic angles. Discover Geoscience** (2025) — https://doi.org/10.1007/s44288-025-00315-2
  - roche : Gneiss oeille (Nepal) ; methode : UCS et essai de charge ponctuelle de 0 a 90 deg
  - anisotropie : Parabole du second ordre resistance vs beta
  - resultat : Forte correlation resistance-angle en parabole ; relation lineaire UCS-PLT
  - limite : Revue nouvelle ; chiffres non dans le resume ; acces libre : oui ; confiance : lu

- **Jacobsson L., Brander L. (2024) Tensile Fracture Initiation and Propagation of Granite and Gneiss at Wedge Splitting Tests: Part 2. Rock Mech. Rock Eng. (Part 1 : Int. J. Fract. 2025, doi 10.1007/s10704-025-00857-z)** (2024) — https://doi.org/10.1007/s00603-024-04257-y
  - roche : Gneiss et granite (Suede) ; methode : Essais de fendage au coin avec DIC et lames minces
  - anisotropie : Pas de balayage d'angle ; microstructure et joints de grains
  - resultat : Gneiss : trajet de fissure lisse et intergranulaire, vitesse plus de 2 fois celle du granite ; longueur critique et FPZ mesurees ; taux de restitution d'energie le long du trajet
  - limite : Pas d'effet d'orientation etudie ; acces libre : oui ; confiance : lu

- **Bercakova A., Melichar R., Obara Y., Ptacek J., Soucek K. (2017) Evaluation of Anisotropy of Fracture Toughness in Brittle Rock, Migmatized Gneiss. Procedia Engineering 191** (2017) — https://doi.org/10.1016/j.proeng.2017.05.260
  - roche : Gneiss migmatitique folie (Rep. tcheque) ; methode : Flexion semi-circulaire (SCB) a plusieurs orientations de la foliation
  - anisotropie : K_Ic selon l'orientation de la foliation
  - resultat : Anisotropie de K_Ic mise en evidence ; methode jugee efficace
  - limite : Resume sans chiffres ; acte de conference ; acces libre : oui ; confiance : lu

- **Perez-Rey I., Masoumi H., Mas Ivars D., Suikkanen J., Alejano L.R. (2025) Size-dependent strength behavior of heterogeneous veined gneiss under uniaxial compression. Rock Mechanics Letters** (2025) — https://doi.org/10.70425/rml.202502.12
  - roche : Gneiss veine (Finlande) ; methode : 114 UCS sur diametres de 14 a 100 mm
  - anisotropie : Aucune ; effet de taille et heterogeneite
  - resultat : Resistance croissante de 14 a 50 mm puis legere baisse jusqu'a 100 mm : loi unifiee d'effet de taille (USEL) applicable aux roches metamorphiques
  - limite : Pas d'anisotropie ; variabilite due aux veines ; acces libre : oui ; confiance : lu

- **Nejati M., Aminzadeh A., Amann F., Saar M.O., Driesner T. (2020) Mode I fracture growth in anisotropic rocks: Theory and experiment. Int. J. Solids Struct. (compagnon : Nejati et al. 2020, Theor. Appl. Fract. Mech. 107:102494, doi 10.1016/j.tafmec.2020.102494, non lu)** (2020) — https://doi.org/10.1016/j.ijsolstr.2020.03.004
  - roche : Granite (granodiorite) folie du Grimsel et argile a Opalinus ; methode : 124 essais de tenacite mode I ; comparaison aux criteres MTS, MERR, MSED
  - anisotropie : Tenacite apparente et angle de deviation selon l'orientation de la foliation
  - resultat : Les fissures de mode I deviennent (kink) vers la foliation plus faible ; tous les criteres surestiment l'angle de deviation et sous-estiment la tenacite apparente ; MTS avec T-stress est le meilleur
  - limite : Roche gneissique granitique, pas un schiste ; paywall ; acces libre : non ; confiance : lu

- **Dutler N., Nejati M., Valley B., Amann F., Molinari G. (2018) On the link between fracture toughness, tensile strength, and fracture process zone in anisotropic rocks. Eng. Fract. Mech. 201:56-79** (2018) — https://doi.org/10.1016/j.engfracmech.2018.08.017
  - roche : Granodiorite foliee du Grimsel ; methode : SCB pour K_Ic, bresiliens pour sigma_t, DIC pour la zone d'elaboration (FPZ)
  - anisotropie : Rapports d'anisotropie de K_Ic, sigma_t et longueur de FPZ
  - resultat : FPZ semi-elliptique, rapport longueur/largeur environ 2 ; longueur de FPZ proportionnelle a (K_Ic/sigma_t)^2, rapport d'anisotropie coherent avec les modeles de type Irwin/strip-yield
  - limite : Roche gneissique granitique ; valeurs dans le corps du texte ; acces libre : oui ; confiance : lu
