# Rapport A2 — Pratique de calibration triaxiale en FDEM / GBM pour les granites (revue de littérature)

Date : 2026-09-02. Contexte : trois calibrations Red Bohus (homogène, Weibull, GBM) sur `rockim_f2`, σ₃ = 20 et 50 MPa à ajuster, 75 et 100 MPa à prédire. Conventions : **[V]** = DOI et contenu vérifiés (résumé éditeur, texte intégral ou API Crossref/Semantic Scholar/OpenAlex) ; **[V-DOI]** = DOI vérifié mais contenu lu seulement en extrait ; **[NON VERIFIE]** = affirmation ou chiffre que je n'ai pas pu confirmer. Les fichiers locaux sont cités `fichier:ligne`.

---

## 0. Cibles locales (rappel chiffré, source fichiers)

| σ₃ [MPa] | pic q [MPa] | ε_pic [%] | chute post-pic (fraction moyenne) | q final [MPa] | CI/q_pic | CD/q_pic | E [GPa] |
|---|---|---|---|---|---|---|---|
| 20 | 404,0 | 0,652 | 0,300 | 294,8 | 0,566–0,569 | 0,615–0,686 | 75–78 |
| 50 | 599,2 | 0,949 | 0,259 | 344,9 | 0,547–0,558 | 0,724–0,728 | 77–78 |
| 75 | 703,0 | 1,076 | 0,247 | 576,4 | 0,519–0,543 | 0,725–0,751 | 77–81 |
| 100 | 798,2 | 1,204 | 0,074 | 774,4 | 0,568–0,601 | 0,592–0,747 | 79–90 |

Sources : `FDEM/rockim_f1/calib_triax3d/targets_triax_bohus.json` (clés `q_peak_mean_MPa`, `eps_peak_microstrain`, `chute_fraction_moyenne`) et `CONTINUUM/calib_bohus_triax/exp_qc/seuils_sbm_bohus.json` (`CI_frac`, `CD_frac`, `E_GPa`, méthode SBM). Attention : les fractions CI/CD sont rapportées au **déviateur q**, pas à σ₁ ; c'est pourquoi elles dépassent les 0,42–0,47 classiques de l'uniaxial (voir §2).

État du code au 2026-09-02 (`calib_quick/README.md:69-73`) : cas homogène 709 MPa (+18 %), Weibull m = 8 sans effet, GBM α = 1 515 MPa (−14 %) ; chute simulée 52–64 % contre 26 % ; ε_pic 0,65–0,87 % contre 0,95 %.

---

## 1. (i) Le bulk des grains : élastique ou élasto-plastique ?

### 1.1 Tableau des codes et des études

| Code / étude | Bulk | Ce qui porte la non-linéarité | Source |
|---|---|---|---|
| **Y-Geo** (Mahabadi, Lisjak, Munjiza, Grasselli 2012) | triangles **élastiques** (isotrope ou transverse-isotrope) | éléments cohésifs 4 nœuds, critère Mohr-Coulomb + coupure en traction, loi de frottement quasi-statique | Int. J. Geomech. 12(6), DOI 10.1061/(ASCE)GM.1943-5622.0000216 **[V]** ; Lisjak et al. 2013 GJI 195:423 (« the bulk material behaves elastically; failure is modeled through embedded cohesive crack elements »), DOI 10.1093/gji/ggt221 **[V]** |
| **Y-Geo GBM Stanstead** (Zhao et al. 2015) | phases quartz / biotite / feldspath **élastiques** (E 83,1 / 17,2 / 56,4 GPa) | joints à c = 24,2 MPa, φ = 51,8° identiques, ft par phase | GJI 203:1246, DOI 10.1093/gji/ggv370 **[V]** (table 1 lue) |
| **Irazu GB-FDEM** (Aboayanah, Abdelaziz, Haile, Zhao, Grasselli 2024) | phases **élastiques** (« principles of non-linear elastic fracture mechanics », aucune plasticité de grain) | joints intra / homophase / hétérophase distincts | RMRE 57:4679–4706, DOI 10.1007/s00603-024-03789-7 **[V]** (texte intégral OA) |
| **Irazu GBM** (Abdelaziz, Zhao, Grasselli 2018) | élastique (idem cadre Irazu) | « explicit modelling of grain boundaries », transition joints de grain → bande de cisaillement intra-phase | Comput. Geotech. 103:73–81, DOI 10.1016/j.compgeo.2018.07.003 **[V]** (résumé) |
| **Irazu GBM** (Zhang et al. 2021) | « elastic strain is simulated by triangular elements based on linear elastic continuum theory » | hiérarchie intra > homophase > hétérophase | Materials 14:3969, DOI 10.3390/ma14143969 **[V]** (OA) |
| **Irazu + FEM-FDEM hybride** (Lisjak, He, Ha, Tatone, Mahabadi 2024) | **élasto-plastique Drucker-Prager / Mohr-Coulomb** ajouté aux éléments finis | motivation explicite : matériaux à « irreversible damage processes under load without breaking » (pentes, roches tendres) | ARMA24, DOI 10.56952/arma-2024-0620 **[V]** |
| **Solidity** (Imperial ; Yang, Xiang, Latham) | **néo-hookéen élastique + viscosité** (`CauchyTet4`, `Y3Dfd.c` l. 716) | joints cohésifs z-curve ; frottement résiduel de joint **nul** (`dpefm = 0`, l. 1091), le résiduel vient du contact | code public LGPL (mémoire `reference-solidity-source-publique.md`) ; Lei, Latham, Xiang 2016 JCM-FEMDEM, RMRE, DOI 10.1007/s00603-016-1064-3 **[V]** |
| **HOSS** (LANL) sur granite Hwangdeung (Euser et al. 2019) | table 1 = E 55 GPa, ν 0,15, ρ 2650 : bulk **élastique** dans cette étude ; le cadre Munjiza admet « non-softening material nonlinearity (i.e., plasticity) » dans le régime durcissant | UCS exp 209 / sim 212,1 MPa ; BD 9,2 / 9,9 MPa | arXiv 1805.06032 (PDF lu) ; RMRE DOI 10.1007/s00603-019-01773-0 **[V-DOI]** |
| **HOSS – cadre plastique** (Rougier et al. 2020) | **plasticité + fracture discrète** en grandes déformations (décomposition multiplicative) | motivation : grandes déformations, matériaux ductiles | IJNME, DOI 10.1002/nme.6255 **[V]** ; HOSS général : Knight et al. 2020 Comp. Part. Mech. 7:765, DOI 10.1007/s40571-020-00349-y **[V]** |
| **Y-HFDEM** (Fukuda et al. 2019, 2020 ; Min, Fukuda et al. 2020) | « the continuous, i.e., elastic, deformations of rock … are modeled » | CE6 cohésifs, pénalités P_n = 10E, P_open = 100E, P_overlap = 1000E | IJNAMG 43, DOI 10.1002/nag.2934 **[V]** ; RMRE 53:1079, DOI 10.1007/s00603-019-01960-z **[V]** ; Min et al. 2020 Appl. Sci. (PDF lu) [DOI NON VERIFIE] |
| **MC-FDEM** (Ye, Zhang, Chen, Li 2025) | **Mohr-Coulomb élasto-plastique** dans les éléments solides + joints cohésifs | motivation : « inadequate characterization of material plastic damage » du FDEM classique ; granite altéré | IJRMMS 194:106233, DOI 10.1016/j.ijrmms.2025.106233 **[V]** ; porté dans rockim `law = mc` (`DOCUMENTATION_rockim.md:644`) |
| **FDEM-Weibull** (Deng, Liu, Lu 2023 ; Deng et al. 2024) | élastique + tirage Weibull des paramètres | hétérogénéité statistique | Comput. Geotech. 154:105138, DOI 10.1016/j.compgeo.2022.105138 **[V-DOI]** ; EABE 168:105924, DOI 10.1016/j.enganabound.2024.105924 **[V-DOI]** (résumés non accessibles) |
| **UDEC-GBM** (Lan, Martin, Hu 2010) | grains polygonaux **élastiques**, hétérogénéité d'E par phase | contacts : CI piloté par l'hétérogénéité géométrique, résistance par les contacts | JGR 115:B01202, DOI 10.1029/2009JB006496 **[V]** |
| **3DEC-GBM** (Wang & Cai 2018, 2019) | grains **rigides, élastiques ou cassables** (contacts intra-grain) comparés | verdict : grains élastiques → « very low residual friction angle of contact » nécessaire et « brittle rock deformation behavior cannot be captured » ; grains cassables : bon post-pic, dilatance correcte | Comput. Geotech. 101:224, DOI 10.1016/j.compgeo.2018.04.016 **[V]** ; IJRMMS 115:60–76, DOI 10.1016/j.ijrmms.2019.01.008 **[V]** (résumé) |
| **UDEC-BBM inélastique** (Sinha, Shirole, Walton 2020) | blocs **inélastiques** (plasticité de bloc) multi-minéraux | ciblent explicitement post-pic, résiduel, dilatance dépendant du confinement, comportement CWFS, 0–60 MPa | JGR 125, DOI 10.1029/2019JB018844 **[V]** |
| **UDEC-GBM** (Chen & Konietzky 2014) | grains élastiques, **contacts élasto-plastiques**, fluage | durée de vie, endommagement sous-critique | Tectonophysics 633:164, DOI 10.1016/j.tecto.2014.06.033 **[V]** |
| **Bohus, continuum** (Saadati et al. 2014 ; Shariati et al. 2018, 2022) | **plasticité KST/DP à dilatance variable** + endommagement anisotrope ; surface de charge jusqu'à 750 MPa de pression | conclusion des essais quasi-œdométriques + indentation + DVC : « compressible elastoplasticity should be accounted for » | IJNAMG 38, DOI 10.1002/nag.2235 **[V]** ; RMRE, DOI 10.1007/s00603-018-1646-3 **[V]** ; IJNAMG, DOI 10.1002/nag.3303 **[V]** ; RMRE 2022 DVC, DOI 10.1007/s00603-022-02991-9 **[V-DOI]** |

### 1.2 Réponse

- **Le standard FDEM et FDEM-GBM (Y-Geo, Irazu, Solidity, Y-HFDEM, HOSS-granite) est un bulk élastique** ; toute la non-linéarité est portée par les joints cohésifs (Mohr-Coulomb + coupure en traction) et, en GBM, par le contraste de raideur entre phases et l'affaiblissement des joints de grain.
- **Qui a mis de la plasticité dans le bulk, et pourquoi** : (a) LANL/HOSS pour les grandes déformations et les matériaux ductiles (Rougier 2020) ; (b) Geomechanica en 2024, hybride FEM-FDEM à DP/MC, pour les pentes et roches qui s'endommagent sans se rompre ; (c) Ye et al. 2025 pour un granite altéré (le FDEM classique décrit mal l'endommagement plastique) ; (d) en DEM-GBM, Sinha 2020 (blocs inélastiques) pour reproduire simultanément post-pic, résiduel et dilatance, et Wang & Cai 2019 qui montrent que des grains purement élastiques ne rendent pas la fragilité et exigent un frottement résiduel de contact quasi nul ; (e) pour le Bohus lui-même, l'école KTH/Grenoble (Saadati, Shariati, Forquin, Hild) n'a jamais pu se passer d'une plasticité compressible à dilatance variable, mais à des pressions (indentation, œdométrique) bien supérieures à 50 MPa.
- Aucun des travaux FDEM-GBM de granite listés (Mahabadi 2014, Abdelaziz 2018, Aboayanah 2024, Fukuda 2019/2020, Zhang 2021) n'utilise de plasticité de grain **[V pour ceux dont le texte a été lu ; Abdelaziz 2018 et Fukuda : déduit du cadre de code, NON VERIFIE ligne à ligne]**. « Liu/Wong » de la question : je n'ai pas identifié de couple Liu–Wong en FDEM ; les travaux pertinents sont Wong & Peng (PFC-GBM, §6) et Liu, Cai & Huang 2018 (UDEC-GBM, hétérogénéité géométrique, DOI 10.1016/j.compgeo.2017.11.013 **[V]**).

---

## 2. (ii) Reproduction de σ_ci et σ_cd

### 2.1 Références expérimentales de cadrage

- Nicksiar & Martin 2013 (376 essais) : CI/σ_pic = **0,42–0,47** en uniaxial « regardless of the material properties », **0,50–0,54** en confiné. Eng. Geol. 154:64–76, DOI 10.1016/j.enggeo.2012.12.007 **[V-DOI, chiffres issus d'un extrait de recherche]**.
- Bohus (fichier `seuils_sbm_bohus.json`) : CI/q_pic = 0,55–0,57 à 20–50 MPa, CD/q_pic = 0,62–0,73 — cohérent avec la montée du rapport sous confinement.
- Aboayanah et al. 2024, granite gris (échantillon S-1 expérimental) : CI = 44,5 % UCS (69 MPa), CD = 84 % UCS (131 MPa) **[V]**.

### 2.2 Ce qui pilote CI et CD dans les modèles (chiffres)

| Levier | Effet mesuré dans la littérature | Source |
|---|---|---|
| **Contraste de raideur entre phases** (biotite 17–29 GPa vs quartz 83 GPa) | Aboayanah 2024 : « macrocracks initiated and propagated along grain boundaries … driven by the high stiffness contrast between the compliant biotite and the stiffer feldspar/quartz » ; CI et CD moyens simulés = **46 % et 88 %** du pic (SBM et AEM) | DOI 10.1007/s00603-024-03789-7 **[V]** |
| **Hétérogénéité géométrique** (forme, taille, distribution de taille) | Lan 2010 : « crack-initiation stress was found to be controlled primarily by microscale geometric heterogeneity », la résistance par l'hétérogénéité des contacts ; Liu, Cai, Huang 2018 : indice d'hétérogénéité géométrique → CI, CD, pic ; Nicksiar & Martin 2014 : CI = mécanisme de traction, cisaillement de joints de grain seulement près du pic, la distribution de taille compte | DOI 10.1029/2009JB006496 **[V]** ; 10.1016/j.compgeo.2017.11.013 **[V]** ; RMRE 47:1165, DOI 10.1007/s00603-013-0451-2 **[V-DOI]** |
| **Résistance des joints hétérophases** | Aboayanah 2024 : ft(Qz-Fsp) = **2,0 MPa** contre 14 (quartz) et 8 (feldspath) « to mimic the higher density of pre-existing cracks » ; G_fI des joints Fsp-Fsp et Bt-Fsp = **1,8–2,0 J/m²** (très fragiles) ; Zhang 2021 : hétérophase ft 10 / intra 20 MPa, G_fI 50 / 900 J/m² | **[V]** |
| **Microfissures préexistantes** | Mahabadi, Tatone, Grasselli 2014 : « microscale heterogeneity and microcracks should be considered to accurately predict the tensile strength and failure behavior » (μCT + densité de microfissures en lame mince) | JGR 119, DOI 10.1002/2014JB011064 **[V]** |
| **Taille de grain** | Aboayanah 2024 : CI, CD et UCS décroissent **log-linéairement** avec d ; granite gris (d ≈ 1,9–3,1 mm) CI 71 / CD 162,7 / UCS 198,8 MPa ; pegmatite (12,5–20 mm) 56 / 128,8 / 139,2 MPa | **[V]** |
| **Biotite** | Aboayanah 2024 : la teneur en biotite abaisse CI et CD ; rapport quartz/feldspath corrélé positivement aux propriétés | **[V]** |
| **Weibull (dispersion des joints)** | Deng et al. 2023/2024 en FDEM-Weibull : existent, chiffres **[NON VERIFIE]** ; localement, m = 8 ne change ni le pic ni ε_pic (`calib_quick/README.md:72,81-83`) et le théorème du « plancher de frottement » (mémoire `project-calib-triax3d-rockim.md`) montre que disperser ft et c sans toucher tan φ ne peut pas produire l'arc σ_cc → σ_ci sous σ₁ = σ₃ tan²(45 + φ/2) | local |
| **Plasticité de bloc** | Sinha 2020 : blocs inélastiques nécessaires pour la non-linéarité axiale pré-pic (le point de non-linéarité σ–ε est retenu comme CD) | **[V]** |

### 2.3 Lecture pour rockim

L'arc σ_ci → σ_cd des granites vient dans la littérature de trois sources : (1) le contraste de raideur (GBM à modules par phase — le cas 3 en a déjà l'effet : rupture à 0,65 % ; `README.md:84-86`), (2) des joints de grain hétérophases beaucoup plus faibles que les joints intra-grain (facteur 0,14–0,5 en ft, 0,05–0,3 en G_fI), (3) des défauts préexistants (Mahabadi 2014 ; `jointPrebrokenFrac` dans rockim). La dispersion Weibull seule ne le fait pas (local, vérifié) ; la plasticité de bulk le fait (Sinha 2020) mais elle n'est pas le mécanisme physique dominant à 20–50 MPa.

---

## 3. (iii) Post-pic et résiduel en triaxial

### 3.1 Frottement résiduel des joints rompus (valeurs)

| Étude | Frottement de pic (joint) | Frottement résiduel / de fracture | Remarque |
|---|---|---|---|
| Lisjak et al. 2013 GJI, table 1 (essai à fissure préexistante, E = 3 GPa — matériau modèle, **[NON VERIFIE : nature exacte]**) | φ interne 35° | **« fracture friction angle » 35°** (paramètre distinct) ; frottement mors–éprouvette 0,1 | Y-Geo distingue φ_f de φ (cf. AbuAisha et al. 2015 éq. 7.5 citée dans `DOCUMENTATION_rockim.md:177`) **[V pour la table ; AbuAisha NON VERIFIE]** |
| Zhao et al. 2015 Stanstead | φ 51,8°, c 24,2 MPa | non extrait **[NON VERIFIE]** | 0–30 MPa de confinement |
| Aboayanah et al. 2024 | μ intra = **1,27** (φ ≈ 52°) | μ homophase **1,14**, hétérophase **0,81** ; mors 0,1 | c 27–35 MPa partout ; le résiduel est le Coulomb du joint rompu + contact |
| Zhang et al. 2021 | μ intra 1,2 | homophase 1,1 ; hétérophase 0,9 | |
| Euser et al. 2019 (HOSS) | φ 35°, c 55,4 MPa | frottement mors 0,5 (explicitement responsable de la transition traction → cisaillement) | |
| Solidity (source) | joint : c, φ | **μ_res joint = 0** (`dpefm = 0`) ; frottement de **contact** 0,6 (calcaire), **0,18** (granite de Kuru) contre 1,85 de pic | `DOCUMENTATION_rockim.md:177` ; chiffres du papier Yang **[NON VERIFIE contre l'article]** |
| Yang et al. 2022 (Abaqus 3D-FDEM, granite) | — | μ = 0,3 | Frontiers Earth Sci., DOI 10.3389/feart.2022.998521 **[V]** |
| Wang & Cai 2019 (3DEC-GBM) | — | grains rigides ou élastiques : « a near-zero residual friction angle of contact is needed to capture post-peak deformation behavior due to the grain interlocking effect » ; grains cassables : « does not require a very low contact residual friction angle » | **[V]** |
| Sci. Rep. 2024 (granite, 0–20 MPa, exp + PFC) | c_pic 21,4 MPa, φ_pic 57,7° | **c_res 4,7 MPa, φ_res 52,6°** ; PFC μ = 0,2 | DOI 10.1038/s41598-024-72834-w **[V-URL]** |

### 3.2 Rôle des blocs polygonaux

- Bahrani, Kaiser, Valley 2014 : marbre « granulé » (cohésion des joints de grain détruite) → UCS < **50 %** de l'intact mais **≈ 80 %** de l'intact à fort confinement ; l'imbrication des grains polygonaux crée à elle seule une résistance confinée. IJRMMS 71:117, DOI 10.1016/j.ijrmms.2014.07.005 **[V-DOI, chiffres d'extrait]**.
- Wang & Cai 2019 : l'imbrication des grains rigides/élastiques **surestime** le résiduel, d'où le besoin d'un φ_res de contact quasi nul ; les grains cassables (joints intra-grain) rendent le post-pic sans cet artifice **[V]**.
- Conséquence pour un GBM à 3 mm sur 20 mm (7 grains en largeur, `README.md:36-39`) : l'imbrication est surreprésentée (peu de grains, pas de fragmentation intra-grain si les joints intra restent forts).

### 3.3 La chute post-pic de 26 % à 50 MPa est-elle reproduite ?

- Expérimental Bohus : chute 30 / 26 / 25 / 7 % à 20 / 50 / 75 / 100 MPa (§0) ; le résiduel à 50 MPa (≈ 345 MPa de déviateur) reste **58 %** du pic. À titre de comparaison, sur les granitoïdes de Bátaapáti (38 essais MFS), q_res = 11,38 σ₃^0,693 (Sci. Rep. 2025, DOI 10.1038/s41598-025-14419-9 **[V-URL]**), soit ≈ 171 MPa à 50 MPa **[calcul personnel]** — le Bohus est nettement plus résistant en résiduel.
- **Aucune étude FDEM 2D vérifiée ne montre la courbe triaxiale post-pic d'un granite à 50 MPa reproduite quantitativement** ; Lisjak, Grasselli, Vietor 2014 (IJRMMS 65:96–115, DOI 10.1016/j.ijrmms.2013.10.006 **[V]**) quantifient l'effet du confinement en triaxial, mais sur l'argile à Opalinus. Les seuls modèles discrets qui revendiquent post-pic + résiduel + dilatance sur granite sont des GBM à blocs inélastiques (Sinha 2020) ou à grains cassables (Wang & Cai 2019) **[V]**.
- Localement, bulk élastique + joints : chute 52–64 % (`README.md:108-117`), aucun mécanisme de résistance résiduelle progressive. Les leviers documentés : μ_res (`jointResidualMu`, pic → résiduel par f(D), `DOCUMENTATION_rockim.md:177`), φ_joint (hétérophase 0,81 chez Grasselli), imbrication (nombre de grains), plasticité de bulk bornée (`law = mc`, `:644`), et vitesse de chargement (l'inertie lisse le post-pic, mémoire).

---

## 4. (iv) Séquence de calibration et paramètres typiques

### 4.1 Séquences publiées

- **Tatone & Grasselli 2015** (IJRMMS 75:56–72, DOI 10.1016/j.ijrmms.2015.01.011 **[V-DOI]**) : procédure de référence ; quatre groupes de paramètres (bulk, amortissement visqueux, rupture, interaction élastique), calés sur UCS, brésilien (BDS) et **biaxial** (Geomechanica, page « FEMDEM model calibration », lue). Le matériau calibré, les vitesses et les valeurs de pénalité du papier : **[NON VERIFIE]** (texte inaccessible ; un extrait de recherche l'associe au granite de Stanstead, mais probablement par confusion avec Zhao 2015).
- **Oliveira, Pinto, Mazzinghy 2021** (Irazu, REM 74, DOI 10.1590/0370-44672020740067 **[V]**), en suivant Tatone : (1) taille d'élément, (2) vitesse de platine, (3) amortissement visqueux, (4) pénalités de contact (≈ **100 E**), (5) énergies de rupture ajustées simultanément sur BTS et UCS. Vitesse 0,03 m/s, élément 0,75 mm.
- **Min, Fukuda et al. 2020** (Appl. Sci., PDF lu) : E, ν, ρ pris de l'expérience ; P_n = 10E, P_open = P_tan = 100E, P_overlap = 1000E ; calibration UCS + BTS ; G_fI 30–60 N/m, G_fII 90–120 N/m ; le frottement de contact est ensuite le paramètre clé du résiduel en cisaillement.
- Autres procédures dédiées : Yan & Tong 2020, calibration des pénalités (Int. J. Geomech. 20(7), DOI 10.1061/(ASCE)GM.1943-5622.0001686 **[V]**) ; Deng, Liu, Lu 2022, calibration des paramètres de joint (EFM 276:108924, DOI 10.1016/j.engfracmech.2022.108924 **[V-DOI]**) ; un extrait de recherche indique qu'une valeur fiable de G_I s'obtient par un brésilien à chemin de rupture pré-tracé **[NON VERIFIE : source exacte]**.
- **GBM** : Wang & Cai 2019 proposent une procédure par type de grain **[V]** ; Ghazvinian, Diederichs, Quey 2014 en 3DEC (JRMGE 6:506, DOI 10.1016/j.jrmge.2014.09.001 **[V]**).

### 4.2 Paramètres finaux publiés pour des granites (bulk élastique)

| Étude / phase | E [GPa] | ν | ft [MPa] | c [MPa] | φ ou μ | G_fI [J/m²] | G_fII [J/m²] | Pénalités |
|---|---|---|---|---|---|---|---|---|
| Zhao 2015, quartz | 83,1 | 0,26 | 11,4 | 24,2 | 51,8° | 907 | 1814 | non extrait |
| Zhao 2015, biotite | 17,2 | 0,30 | 4,2 | 24,2 | 51,8° | 599 | 1198 | |
| Zhao 2015, feldspath | 56,4 | 0,45 [valeur extraite automatiquement, douteuse] | 5,5 | 24,2 | 51,8° | 310 | 620 | |
| Aboayanah 2024, feldspath | 60 | 0,32 | 8 | 33 | μ 1,27 | 390 | 690 | rupture 564 GPa ; normal 112 GPa/m |
| Aboayanah 2024, biotite | 29,3 | 0,36 | 5,5 | 30 | μ 1,27 | 599 | 1198 | 172 ; 344 |
| Aboayanah 2024, quartz | 83,1 | 0,17 | 14 | 35 | μ 1,27 | 907 | 1810 | 832 ; 166 |
| Aboayanah 2024, joints homophases (Bt-Bt / Fsp-Fsp / Qz-Qz) | — | — | 4,95 / 7,2 / 12,6 | 27 / 31 / 31,5 | μ 1,14 | 449 / 1,98 / 680 | 899 / 465 / 1361 | |
| Aboayanah 2024, joints hétérophases (Bt-Fsp / Bt-Qz / Qz-Fsp) | — | — | 8,5 / 8,5 / 2,0 | 32,5 / 31 / 32 | μ 0,81 | 1,81 / 1,07 / 300 | 382 / 382 / 1450 | |
| Zhang 2021, quartz intra / Qz-Qz / Qz-Fsp | 94,5 | 0,08 | 20 / 15 / 10 | 25 / 20 / 20 | μ 1,2 / 1,1 / 0,9 | 900 / 700 / 50 | 1800 / 1400 / 500 | |
| Euser 2019, Hwangdeung (homogène, HOSS) | 55 | 0,15 | 9,2 | 55,4 | 35° | non extrait | | élément 1,0 mm |
| Lisjak 2013, matériau modèle | 3 | 0,29 | 3 | 15 | 35° (int.) / 35° (fracture) | 2 | 10 | rupture 15 GPa ; contact 30 GPa·m |

Constantes de pratique : G_fII/G_fI ≈ **2** dans le groupe Grasselli (toutes les lignes ci-dessus) contre 10–20 par défaut dans rockim (`DOCUMENTATION_rockim.md:133`, `README.md:21`) ; pénalité 100–1000 E en 2D (Min 2020 ; extraits citant Tatone) contre 20 E/h dans rockim (`:171`, complaisance 4–5 % sur E). Les valeurs quartz 83,1/14/35 et biotite 29,3/5,5/30 du deck `q3_gbm_P050.cfg` (`README.md:32-35`) sont celles d'Aboayanah 2024 ; le feldspath y est recalé pour l'équivalence des moyennes.

### 4.3 Séquence recommandée pour Bohus (synthèse des sources + contraintes locales)

1. **E, ν par phase** (littérature ci-dessus), puis correction déformation plane E/(1−ν²) : viser ≈ 72 GPa dans le deck pour 77 apparents (`README.md:87-88`).
2. **ft et G_fI** sur la traction directe (le brésilien de rockim n'est pas encore valide, mémoire `project-rockim-gbm.md`) en gardant **ℓ_cz = E G_f / ft² constant** (`README.md:98-104`) ; G_fII/G_fI ≈ 2 comme repère littérature.
3. **c, φ des joints** sur l'UCS (UCS/σ_t réel ≈ 6,9 ; la littérature GBM met φ_joint à 13–52° selon la représentation).
4. **Triaxial 20 et 50 MPa** : μ_res (résiduel), φ_joint (pente), puis seulement si nécessaire un bulk MC borné (§8), avec la règle d'ordre des enveloppes (joints nettement sous le bulk ; mémoire).
5. **Prédire 75 et 100 MPa** sans retouche ; la chute de 7 % à 100 MPa est le test le plus discriminant du résiduel.

---

## 5. (v) Vitesse de chargement, pas de temps, quasi-staticité

| Étude | Vitesse de platine | ε̇ équivalent | Justification donnée |
|---|---|---|---|
| Lisjak et al. 2013 (Y-Geo) | 0,05 m/s | — | « constant strain rate » **[V]** |
| Zhao et al. 2015 (Stanstead, 0–30 MPa) | **0,25 m/s** | **≈ 2,31 s⁻¹** | — **[V]** |
| Aboayanah et al. 2024 (Irazu GBM) | 0,125 m/s | 1,8 s⁻¹ | « high compared to laboratory testing, however … reasonable run time without compromising the results » ; dt 10⁻⁸–10⁻⁷ s **[V]** |
| Zhang et al. 2021 (Irazu GBM) | 0,05 m/s | — | **[V]** |
| Oliveira et al. 2021 (Irazu) | 0,03 m/s | — | fixée lors de la calibration UCS **[V]** |
| Euser et al. 2019 (HOSS) | 0,1 m/s | — | « constant strain rate such that quasi-static behavior is maintained » **[V]** |
| Min, Fukuda et al. 2020 (Y-HFDEM 3D) | 0,01 m/s + amortissement critique | — | « to approximately achieve a quasi-static loading condition » **[V]** ; 0,05 m/s dans d'autres travaux du groupe Fukuda **[NON VERIFIE : quel article]** |
| Yang et al. 2022 (Abaqus 3D-FDEM) | 0,001 m/s | — | **[V]** |
| Liu & Deng 2019 (FDEM, étude systématique) | **≤ 0,5 m/s** en laboratoire ; **≥ 27–28 éléments** dans le diamètre | — | EFM 211:442–462, DOI 10.1016/j.engfracmech.2019.02.007 **[V]** |
| rockim (sonde 4) | 0,25 m/s sur 40 mm | ≈ 5–6 s⁻¹ (`README.md:40`) | post-pic « lissé par construction », dampingLocal 0,7 amplifie (mémoire) |

Lecture : la vitesse de rockim est dans la fourchette publiée (Zhao 2015 est à 0,25 m/s, mais sur une éprouvette plus haute : ε̇ 2,3 s⁻¹ contre 5–6), et sous le plafond de Liu & Deng ; c'est le **rapport vitesse / hauteur** et l'amortissement qui posent problème, pas la vitesse absolue. Diviser pullV par 5–10 (ε̇ ≈ 0,5–1 s⁻¹) ramène dans le régime des références Irazu.

- **Régime dynamique matériau** : la revue Zhang & Zhao 2014 (RMRE 47:1411, DOI 10.1007/s00603-013-0463-y **[V-DOI]**) situe l'amplification dynamique de la résistance en compression des roches à des ε̇ nettement supérieurs à 10 s⁻¹ **[NON VERIFIE : seuil exact]** ; à 5 s⁻¹ l'effet est donc inertiel/numérique (lissage du post-pic), pas constitutif.
- **Mass scaling** : aucune des études FDEM de laboratoire ci-dessus n'en utilise ; elles jouent sur l'amortissement (critique, local). Étude systématique du mass/time scaling en mécanique des roches : Heinze et al. 2016, Tectonophysics 684:4–11, DOI 10.1016/j.tecto.2015.10.013 **[V-DOI, contenu NON VERIFIE]**. Règle de thèse déjà en place : pas de mass scaling global (`CLAUDE.md`, règle 4).
- **Vérification de quasi-staticité** : en DEM granulaire, énergie cinétique / énergie potentielle < 10⁻⁵ et nombre d'inertie < 10⁻⁴ (arXiv 1506.00439, extrait **[NON VERIFIE]**) ; en FDEM, aucun seuil chiffré publié n'a été trouvé — la pratique est la comparaison de deux vitesses (Aboayanah) ou la borne de Liu & Deng. Critères locaux disponibles : bilan de platines (|F_top|−|F_bot|)/moyenne, dampWork ≤ 0, énergie cinétique négligeable devant le travail (`DOCUMENTATION_rockim.md:1047-1048`, mémoire). Proposition opérationnelle : KE/W_ext < 1 % avant le pic, bilan de platines < 2–5 %, et invariance du pic à 3 % près quand on divise pullV par 2.

---

## 6. (vi) GBM : rapports inter/intra-granulaires, taille et dispersion des grains

### 6.1 Rapports utilisés

| Étude | Représentation | ft_inter / ft_intra | c_inter / c_intra | μ_inter / μ_intra | G_fI inter / intra |
|---|---|---|---|---|---|
| Aboayanah 2024 (Irazu) | FDEM 3 phases | homophase **0,9** (7,2/8 ; 12,6/14 ; 4,95/5,5) ; Bt-Fsp 8,5 (≈ 1) ; **Qz-Fsp 2,0 → 0,14–0,25** | ≈ 0,9–1,0 | 1,14/1,27 = 0,90 ; 0,81/1,27 = **0,64** | Fsp-Fsp 1,98/390 = **0,005** ; Qz-Qz 0,75 ; Qz-Fsp 300/907 = 0,33 |
| Zhang 2021 (Irazu) | FDEM 4 phases | homophase 0,75 ; hétérophase **0,5** | 0,8 | 0,92 ; 0,75 | 0,78 ; **0,056** |
| Hofmann, Babadagli, Yoon, Zang, Zimmermann 2015 (PFC2D, Aue) | grains multi-particules, contacts intra vs inter distincts | **[NON VERIFIE : valeurs]** | | | EFM 147:261–275, DOI 10.1016/j.engfracmech.2015.09.008 **[V-DOI]** |
| Peng, Wong, Teh, Li 2018 (PFC2D, Bukit Timah) | smooth joints aux joints de grain | **[NON VERIFIE : valeurs]** ; erreurs < ±6 % sur les macro-propriétés | | | RMRE 51:135, DOI 10.1007/s00603-017-1316-x **[V-DOI]** |
| Li, Zhang, Li, Zhao 2018 (GB-DEM) | PFC, calibration quasi-statique puis SHPB | **[NON VERIFIE]** | | | RMRE 51:3785, DOI 10.1007/s00603-018-1566-2 **[V-DOI]** |
| Saadat & Taheri 2019 (PFC, Aue) | contacts intra-grain cohésifs adoucissants | **[NON VERIFIE]** | | | Comput. Geotech. 111:89, DOI 10.1016/j.compgeo.2019.03.009 **[V]** |
| Wang & Cai 2018 (3DEC) | inter + intra explicites | étude paramétrique : l'hétérogénéité inter/intra pilote pic **et** trajet des fissures **[V]** ; ratios **[NON VERIFIE]** | | | |
| Bahrani 2014 (PFC) | cohésion inter = 0 | 0 | 0 | — | UCS < 50 % intact, 80 % confiné |
| Wang, Zhou et al. 2021 (GBM, étude paramétrique de la résistance des joints minéraux) | — | **[NON VERIFIE]** | | | EFM 241:107388, DOI 10.1016/j.engfracmech.2020.107388 **[V-DOI]** |
| rockim `q3` | `gbAlpha*` = 1 (`README.md:28-31`) ; ancienne calibration α = 0,5 (mémoire, invalidée) | 1 | 1 | 1 | 1 |

Retenir : la littérature FDEM-GBM affaiblit **peu** les joints homophases (0,75–0,9) mais fortement les joints **hétérophases** (0,14–0,5 en ft, 0,05–0,3 en G_fI, 0,64–0,75 en frottement). Dans rockim cela correspond à `gbHeteroFactor` ≪ `gbAlphaTen` (`DOCUMENTATION_rockim.md:145-147`).

### 6.2 Taille et dispersion des grains

- **Taille** : Aboayanah 2024, UCS / CI / CD = A·log(d) + B, décroissants (chiffres §2.2) **[V]** ; Peng, Wong, Teh 2021 (JRMGE 13:755, DOI 10.1016/j.jrmge.2021.01.011 **[V]**) : le DEM-GBM donne au contraire une résistance **croissante** avec d, et il faut **dégrader les paramètres de joint de grain avec d** pour retrouver la tendance expérimentale ; Saadat & Taheri 2019 : plus gros grains → plus forte résistance des éprouvettes pré-fissurées (Aue) **[V]**. Le sens de l'effet de taille est donc un **artefact possible du modèle** : à surveiller lors du passage 2 → 3 mm.
- **Dispersion de taille (polydispersité)** : Peng, Wong, Teh 2017 (JGR 122:1054, DOI 10.1002/2016JB013469 **[V]**) : le pic et E **augmentent** quand le modèle passe d'hétérogène à homogène en taille ; les fissures de traction aux joints diminuent et les fissures intra-grain augmentent ; la fissuration se développe préférentiellement dans le quartz ; en traction directe, la dispersion ne change guère la courbe mais déplace la fracture. Nicksiar & Martin 2014 : la distribution de taille agit sur CI **[V-DOI]**. Liu, Cai, Huang 2018 : indice d'hétérogénéité géométrique → CI, CD, pic **[V]**. Peng, Wong, Teh 2017 IJRMMS 100:207 (DOI 10.1016/j.ijrmms.2017.10.004 **[V-DOI]**) : rapport taille de grain / taille de particule.
- **Discrétisation et nombre de grains** : Aboayanah : maillage 0,3 mm pour d ≤ 1,7 mm, 0,5 mm au-delà (h/d ≈ 0,1–0,2), éprouvette 61 × 137 mm, 3135 grains (granite gris) ; Zhang 2021 : 0,5 mm pour d 1,5–3 mm ; Liu & Deng : ≥ 27–28 éléments dans le diamètre. rockim : 0,18 d (mémoire), 113 grains de 3 mm sur 20 × 40 mm — la réserve « c'est le NOMBRE de grains » (`README.md:36-39`) est confirmée par la pratique publiée (≈ 20–30 grains en largeur).

---

## 7. (vii) Travaux FDEM/DEM sur le Bohus

- **Continuum FEM (Abaqus), école KTH–Grenoble–ENS** : Saadati, Forquin, Weddfelt, Larsson, Hild 2014 (KST-DFH plasticité + endommagement anisotrope, percussion, DOI 10.1002/nag.2235 **[V]**) ; Saadati et al. 2015 (fissures préexistantes, IJNAMG, DOI 10.1002/nag.2331 **[V-DOI]**) ; Saadati et al. 2016 (traction dynamique, Adv. Mater. Sci. Eng., DOI 10.1155/2016/6279571 **[V-URL]**) ; Shariati et al. 2018 (DP à dilatance variable, quasi-œdométrique + indentation, DOI 10.1007/s00603-018-1646-3 **[V]**) ; Shariati et al. 2022 (endommagement anisotrope + Weibull en indentation, DOI 10.1002/nag.3303 **[V]**) ; Shariati, Bouterf, Saadati, Larsson, Hild 2022 (indentation in situ + DVC : « compressible elastoplasticity should be accounted for », DOI 10.1007/s00603-022-02991-9 **[V-DOI]**) ; thèses KTH de Saadati (2015) et Shariati (2019) sur DiVA **[NON VERIFIE : identifiants]**.
- **DEM / FDEM / GBM sur Bohus** : trois recherches ciblées (« Bohus granite » × PFC/UDEC/grain-based/DEM/FDEM/Y-Geo/Irazu/bonded particle/Voronoi) n'ont retourné **aucun** travail ; un extrait mentionne seulement une caractérisation expérimentale des microfissures du Bohus sous chargement cyclique **[NON VERIFIE]**. Conclusion : à ma connaissance, la calibration FDEM-GBM du Red Bohus serait une première ; les seules données constitutives publiées (DP à dilatance variable, surface de charge jusqu'à 750 MPa, ft Weibull) sont celles de l'école KTH — cohérentes avec la carte DP-DFH de la thèse (`DOCUMENTATION_rockim.md:647`).

---

## 8. VERDICT — quelle représentation du bulk pour un granite à 20–50 MPa ?

**Représentation de base recommandée : bulk élastique par phase (GBM) avec joints intra / homophases / hétérophases distincts.** C'est le standard vérifié de tous les FDEM-GBM de granite (Y-Geo, Irazu, Y-HFDEM, HOSS) et il suffit, dans la littérature, à reproduire E, UCS, BTS, **σ_ci ≈ 0,45 et σ_cd ≈ 0,85–0,88 UCS** et la transition axiale/cisaillement avec la taille de grain (Aboayanah 2024 **[V]**), parce que la non-linéarité pré-pic d'un granite est physiquement due au contraste de raideur biotite/quartz-feldspath, aux joints hétérophases très faibles (ft 0,14–0,5, G_fI 0,05–0,3 de l'intra-grain) et aux microfissures préexistantes — pas à une plasticité de grain. À 20–50 MPa, le Bohus reste fragile (chute 30 et 26 %, `targets`), ce qui est le domaine de validité de cette représentation.

**Là où le bulk élastique atteint sa limite, et ce que dit la littérature** : (1) le **résiduel** — en FDEM 2D à joints Coulomb, il n'est porté que par μ_res et l'imbrication ; Wang & Cai 2019 **[V]** montrent que des grains élastiques imposent un frottement résiduel de contact quasi nul pour rendre le post-pic et « ne peuvent pas capturer le comportement fragile », et Bahrani 2014 que l'imbrication seule donne 80 % de la résistance intacte sous confinement ; (2) la **concavité de l'enveloppe** (pentes locales 13,9 puis 6,5 pour Bohus, mémoire) — un Coulomb de joint linéaire ne la produit pas ; (3) l'**ε_pic ≈ 1 %** — le bulk élastique casse à 0,65–0,87 % (`README.md:71-73`). Les seuls modèles discrets qui revendiquent post-pic, résiduel et dilatance ensemble sont à blocs inélastiques (Sinha 2020 **[V]**) ou à grains cassables ; pour le Bohus, l'école KTH n'a jamais pu se passer d'une plasticité compressible à dilatance variable, mais pour des pressions bien au-delà de 50 MPa (Shariati 2018/2022 **[V]**).

**Décision proposée** : mener les trois calibrations (homogène, Weibull, GBM) à **bulk élastique**, en calibrant le résiduel par `jointResidualMu` et le frottement hétérophase (0,64–0,75 du pic selon Grasselli), et n'introduire la **plasticité MC non associée du bulk** (`law = mc`, Ye 2025 **[V]**, ψ petit, enveloppe de bulk nettement au-dessus de celle des joints) que comme **correction bornée pour 50 MPa et la prédiction 75/100 MPa** — jamais comme substitut au mécanisme de fissuration, et jamais en attendant d'elle l'arc σ_ci → σ_cd, qui doit venir de l'hétérogénéité (GBM, joints hétérophases faibles, `jointPrebrokenFrac`), la dispersion Weibull des seuls ft/c étant démontrée inefficace localement. Priorités de protocole avant tout jugement sur le bulk : ε̇ ÷ 5–10 (les références Irazu sont à 0,03–0,125 m/s sur des éprouvettes 2–3 fois plus hautes ; Liu & Deng ≤ 0,5 m/s), G_fII/G_fI ≈ 2, pénalité 100 E, et ≥ 20 grains en largeur — les quatre écarts qui, d'après cette revue, expliquent avant les lois de bulk les 52–64 % de chute et les pics décalés observés.

---

## 9. Références vérifiées (DOI)

Mahabadi et al. 2012 Int. J. Geomech. 10.1061/(ASCE)GM.1943-5622.0000216 · Mahabadi, Randall, Zong, Grasselli 2012 GRL 10.1029/2011GL050411 · Mahabadi, Tatone, Grasselli 2014 JGR 10.1002/2014JB011064 · Lisjak et al. 2013 GJI 10.1093/gji/ggt221 · Lisjak & Grasselli 2014 JRMGE 10.1016/j.jrmge.2013.12.007 · Lisjak, Grasselli, Vietor 2014 IJRMMS 10.1016/j.ijrmms.2013.10.006 · Tatone & Grasselli 2015 IJRMMS 10.1016/j.ijrmms.2015.01.011 · Zhao et al. 2015 GJI 10.1093/gji/ggv370 · Abdelaziz, Zhao, Grasselli 2018 10.1016/j.compgeo.2018.07.003 · Aboayanah et al. 2024 RMRE 10.1007/s00603-024-03789-7 · Lisjak et al. 2024 ARMA 10.56952/arma-2024-0620 · Zhang et al. 2021 Materials 10.3390/ma14143969 · Oliveira et al. 2021 REM 10.1590/0370-44672020740067 · Fukuda et al. 2019 IJNAMG 10.1002/nag.2934 · Fukuda et al. 2020 RMRE 10.1007/s00603-019-01960-z · Liu & Deng 2019 EFM 10.1016/j.engfracmech.2019.02.007 · Deng, Liu, Lu 2022 EFM 10.1016/j.engfracmech.2022.108924 · Deng, Liu, Lu 2023 CG 10.1016/j.compgeo.2022.105138 · Deng et al. 2024 EABE 10.1016/j.enganabound.2024.105924 · Yan & Tong 2020 IJG 10.1061/(ASCE)GM.1943-5622.0001686 · Ye et al. 2025 IJRMMS 10.1016/j.ijrmms.2025.106233 · Knight et al. 2020 CPM 10.1007/s40571-020-00349-y · Rougier et al. 2020 IJNME 10.1002/nme.6255 · Euser et al. 2019 RMRE 10.1007/s00603-019-01773-0 (arXiv 1805.06032) · Lei, Latham, Xiang 2016 RMRE 10.1007/s00603-016-1064-3 · Yang et al. 2022 Front. Earth Sci. 10.3389/feart.2022.998521 · Lan, Martin, Hu 2010 JGR 10.1029/2009JB006496 · Nicksiar & Martin 2013 Eng. Geol. 10.1016/j.enggeo.2012.12.007 · Nicksiar & Martin 2014 RMRE 10.1007/s00603-013-0451-2 · Ghazvinian, Diederichs, Quey 2014 JRMGE 10.1016/j.jrmge.2014.09.001 · Bahrani, Kaiser, Valley 2014 IJRMMS 10.1016/j.ijrmms.2014.07.005 · Hofmann et al. 2015 EFM 10.1016/j.engfracmech.2015.09.008 · Chen & Konietzky 2014 Tectonophysics 10.1016/j.tecto.2014.06.033 · Peng, Wong, Teh 2017 JGR 10.1002/2016JB013469 · Peng, Wong, Teh 2017 IJRMMS 10.1016/j.ijrmms.2017.10.004 · Peng et al. 2018 RMRE 10.1007/s00603-017-1316-x · Peng, Wong, Teh 2021 JRMGE 10.1016/j.jrmge.2021.01.011 · Li, Zhang, Li, Zhao 2018 RMRE 10.1007/s00603-018-1566-2 · Wang & Cai 2018 CG 10.1016/j.compgeo.2018.04.016 · Wang & Cai 2019 IJRMMS 10.1016/j.ijrmms.2019.01.008 · Liu, Cai, Huang 2018 CG 10.1016/j.compgeo.2017.11.013 · Saadat & Taheri 2019 CG 10.1016/j.compgeo.2019.03.009 · Sinha, Shirole, Walton 2020 JGR 10.1029/2019JB018844 · Wang, Zhou et al. 2021 EFM 10.1016/j.engfracmech.2020.107388 · Song 2025 IJNAMG 10.1002/nag.70077 · Zhang & Zhao 2014 RMRE 10.1007/s00603-013-0463-y · Heinze et al. 2016 Tectonophysics 10.1016/j.tecto.2015.10.013 · Li, Zhang, Jiang 2024 TAFM 10.1016/j.tafmec.2024.104408 · Saadati et al. 2014 IJNAMG 10.1002/nag.2235 · Saadati et al. 2015 IJNAMG 10.1002/nag.2331 · Shariati et al. 2018 RMRE 10.1007/s00603-018-1646-3 · Shariati et al. 2022 IJNAMG 10.1002/nag.3303 · Shariati et al. 2022 RMRE 10.1007/s00603-022-02991-9 · Sci. Rep. 2024 granite strain-softening 10.1038/s41598-024-72834-w · Sci. Rep. 2025 résiduel granitoïdes 10.1038/s41598-025-14419-9.

Non vérifiés ou inaccessibles : contenu de Tatone & Grasselli 2015 (matériau, vitesses, pénalités) ; valeurs numériques de Hofmann 2015, Peng 2018, Li 2018, Wang & Cai 2018 (ratios) ; AbuAisha et al. 2015 éq. 7.5 ; chiffres Solidity de Yang (0,18 / 1,85) contre l'article ; DOI de Min, Fukuda et al. 2020 (Appl. Sci.) ; seuil dynamique de Zhang & Zhao 2014 ; ν = 0,45 du feldspath dans la table de Zhao 2015 (extraction automatique).

Fichiers locaux cités : `C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_f2\calib_quick\README.md` ; `C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_f2\DOCUMENTATION_rockim.md` ; `C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_f1\calib_triax3d\targets_triax_bohus.json` ; `C:\Users\fuzquianoalricabi\simulations\CONTINUUM\calib_bohus_triax\exp_qc\seuils_sbm_bohus.json` ; `C:\Users\fuzquianoalricabi\simulations\CONTINUUM\calib_bohus_triax\exp_qc\experimental_data_red_bohus_clean.json`.