# TÂCHE A1 — Pétrographie réelle du Red Bohus (granite de Bohus, Bohuslän, Suède)

*Rapport du 2026-09-02. Convention : [VÉRIFIÉ] = document lu intégralement ou page/PDF consulté ; [NON VÉRIFIÉ] = non confirmé sur une source primaire ; [ESTIMATION] = extrapolation de l'agent, à ne pas citer comme donnée.*

---

## 0. Réponse courte

| Question | Réponse | Confiance |
|---|---|---|
| D'où viennent les **62/31/7** ? | Du docstring de `rockim_f1/calib_triax3d/gen_gbm.py:8-11` : « Composition Red Bohus **PROVISOIRE** (littérature granite de Bohus, À CONFIRMER par la pétrographie ORCHYD) : feldspath 62 %, quartz 31 %, biotite 7 % ; grains ~2 mm », repris tel quel dans `CAMPAGNE.md:360` et `calib_quick/README.md:14`. **Aucune source primaire** n'est citée. Coïncidence troublante : dans Dumoulin et al. 2024 §2, c'est le **Kuru Grey** qui est décrit à « feldspars 62 wt% … quartz 31.7 wt% » — les 62/31 semblent avoir été empruntés au granite finlandais [hypothèse, NON VÉRIFIÉ]. | — |
| D'où vient le **« 60/35/5 (poids) »** ? | Verbatim de Dumoulin et al. 2024 §2 : « feldspar … *of the order of* 60 wt%, quartz *of the order of* 35 wt%, biotite *of the order of* 5 wt% » (ordre de grandeur, pas un comptage de points) → `calibration_redbohus/README.md:12` et `targets_redbohus.json:6`. | [VÉRIFIÉ] |
| D'où vient le **« 1-3 mm »** ? | `calibration_redbohus/README.md:13` et `phd/HETEROGENEITE.md:19` seulement. **Absent de Dumoulin 2024**, qui classe Red Bohus « medium-grained » avec la définition « medium (1–5 mm) » et note des feldspaths « larger than several millimetres ». Le 1-3 mm est une paraphrase non sourcée. | [NON VÉRIFIÉ] |
| Les 62/31/7 sont-ils faux ? | **Non, quasi justes** : quatre comptages de points EN 12407 sur les carrières rouges de Brastad (Skarstad, Broberg, « Brastad ») donnent en moyenne **feldspath 62 %, quartz 32 %, biotite 4 %, opaques 1–2 %**. Seule la biotite (7 %) est surestimée : 4–6 % en réalité (6 % au maximum si l'on y agrège opaques + chlorite). | [VÉRIFIÉ] |
| Tailles de grain par minéral ? | **Aucune mesure publiée par minéral** n'a été trouvée pour le Red Bohus. Seules des classes qualitatives convergentes : « medium-grained » (1–5 mm, ORCHYD), « Kornstruktur : Medel » sur une échelle fin < 3 mm / grossier > 10 mm (Stenkartoteket), « medium to coarse grained, some porphyric parts » (SGU), feldspaths pluri-millimétriques, tendance porphyrique du feldspath K. Le tableau §6 est donc une **estimation raisonnée** à confirmer par une mesure d'interceptes sur une photo de carotte (1 h de travail, voir §7). | [ESTIMATION] |

---

## 1. Identification de la roche testée

- **Roche ORCHYD** : « Red Bohus granite (medium grained) is extracted near the cities of **Brastad and Lysekil** in Sweden. The batholith dates from 922 million years, is 100 km long and 22 km wide » — Dumoulin S., Thenevin I., Kane A., Rouabhi A., Latham J.-P., Jahangir E., Sellami H. (2024), *A complete experimental study on hard granites: microstructural characterization, mechanical response, and failure criterion*, Geomech. Energy Environ. 40, 100592, **DOI 10.1016/j.gete.2024.100592** (CC-BY, texte intégral lu) [VÉRIFIÉ]. Données : Zenodo 10.5281/zenodo.10617548 (description DataCite lue ; le zip lui-même n'a pas été ouvert).
- Les carrières de Brastad appartiennent à **Hallindens Granit AB** (Stenbrottet Skarstad 501, 454 92 Brastad) : variétés **Skarstad Red Bohus** (58.474°N 11.507°E), **Broberg Red Bohus** (58.473°N 11.440°E), et la fiche Stenkartoteket « Röd Bohus Brastad » (carrière de Prästtorp, Lysekils kommun, 58°23'16"N 11°30'04"E) [VÉRIFIÉ, fiches lues]. Ce sont donc bien les fiches de **ce** granite. Laquelle des carrières a fourni le bloc ORCHYD n'est pas dit dans l'article [NON VÉRIFIÉ].
- **Contexte géologique** : granite post-cinématique sveconorvégien, 920 ± 5 Ma (Eliasson & Schöberg 1991, Precambrian Res. 51, DOI 10.1016/0301-9268(91)90107-L), massif composite de plusieurs intrusions « with rather homogeneous mineralogy », monzogranite ; « colour ranges from red to reddish grey … grain size ranges from medium grained to coarse grained and even with some porphyric parts » — Schouenborg B. & Eliasson T. (SGU), EGU 2015, abstract EGU2015-14883 [VÉRIFIÉ, PDF lu]. Prolongement norvégien = granite d'Iddefjord.
- **Attention à l'homonymie** : le « Bohus granite » de la lignée KTH/Atlas Copco (Saadati, Shariati, Weddfelt) a la composition **quartz 33 / plagioclase 33 / feldspath K 29 / biotite 6 vol%** « tested by SP » (Saadati 2015, thèse KTH n° 89, p. 4 [VÉRIFIÉ] ; Shariati et al. 2019 RMRE 52:645, DOI 10.1007/s00603-018-1646-3 §1 [VÉRIFIÉ]). C'est, au pour-cent près, la fiche Stenkartoteket du **Grå Bohus Evja** (Kungshamn, testé par SP : plag 33 / qz 33 / Kfs 29 / bt 4) [VÉRIFIÉ]. La roche KTH est donc vraisemblablement un Bohus **gris** de la région Ävja/Kungshamn, pas le rouge de Brastad [hypothèse, NON VÉRIFIÉ]. Les paramètres DFH « historiques » du dépôt (m = 24, E = 52 GPa, ν = 0,15-0,25) viennent de cette lignée.

---

## 2. Composition modale

### 2.1 Comptages de points EN 12407 sur les carrières rouges de Brastad (les meilleures données)

| Variété (labo, année) | Quartz | Feldspath K | Plagioclase | Biotite | Opaques | Autres | Source |
|---|---|---|---|---|---|---|---|
| Röd Bohus **Brastad** (CBI/SP, fiche 2011) | 30 | 32 | 30 | 6 | 2 | — | sten.se/wp-content/kartotek/pdf/Bohus_Brastad.pdf [VÉRIFIÉ] |
| Röd Bohus **Broberg** (CBI/SP, fiche 2011) | 35 | 30 (microcline) | 29 | 4 | 1 | — | sten.se/…/Bohus_Broberg.pdf [VÉRIFIÉ] |
| **Skarstad** Red Bohus (RISE, 2024-05-14) | 32 | 36 | 27 | 3 | 0,8 | chlorite 0,2 | hallindensgranit.se/…/SKARSTAD-Red-Bohus-Tekniskt-datablad-eng-2024-05-14.pdf [VÉRIFIÉ] |
| **Broberg** Red Bohus (RISE, 2024-05-14) | 30 | 35 | 30 | 4 | 1,4 | 0,4 | hallindensgranit.se/english/broberg-red-bohus/ [VÉRIFIÉ, page lue] |
| **Moyenne des 4** | **31,8** | **33,3** | **29,0** | **4,3** | **1,3** | ~0,3 | — |
| Grå Bohus Evja (SP, 2011) — pour comparaison | 33 | 29 | 33 | 4 | — | — | sten.se/…/BKS_Gra_BohEvjaPol.pdf [VÉRIFIÉ] |

Remarques : (i) EN 12407 = comptage sur lame mince, fractions **surfaciques** — exactement ce que rockim attend ; (ii) feldspath total = **62 %** (Kfs 33 + plag 29), quartz **32 %**, biotite **4 %**, opaques+chlorite **1,5 %** ; (iii) la dispersion inter-carrières est de ±3 % par phase — c'est l'incertitude à retenir.

### 2.2 Autres sources

- Dumoulin 2024 §2 : feldspath ~60 wt%, quartz ~35 wt%, biotite ~5 wt% (« of the order of ») [VÉRIFIÉ]. Converti en volume avec ρ = 2,57/2,65/3,0 : feldspath 61 %, quartz 35 %, biotite 4,4 % — cohérent avec §2.1. Conclusion de l'article : « roughly the same mineral content (quartz 30 %, feldspar 60 %, mica 10 %) » pour les trois granites.
- Saadati 2015 / Shariati 2019 (Bohus KTH, vraisemblablement gris) : qz 33 / plag 33 / Kfs 29 / bt 6 vol% [VÉRIFIÉ].
- Accessoires (SGU/Wikipédia d'après Eliasson & Schöberg 1991) : magnétite, apatite, zircon, titanite, grenat, monazite ; localement prehnite, calcite, chlorite (altération) [VÉRIFIÉ sur la page Wikipédia, non sur l'original]. Muscovite magmatique locale (granites à deux micas), biotite = « det helt dominerende mørke mineral », pas de hornblende (rapakivi.dk, d'après Asklund 1947) [VÉRIFIÉ sur le site, source secondaire].
- Fourchette « quartz 25–35 vol%, feldspath alcalin 40–50 %, plagioclase 10–20 %, biotite 5–10 % » renvoyée par le moteur de recherche (page Grokipedia, générée par IA) : **incompatible** avec les comptages (Kfs 30–36 %, plag 27–33 %) → **à rejeter** [NON VÉRIFIÉ].
- Vérification de cohérence densité : 0,32×2,65 + 0,33×2,56 + 0,29×2,65 + 0,045×3,0 + 0,015×5 ≈ **2,67** g/cm³ vs 2,63–2,66 mesurés → composition plausible.

---

## 3. Texture et tailles de grain

### 3.1 Ce que disent les sources (toutes qualitatives)

| Source | Énoncé | Statut |
|---|---|---|
| Dumoulin 2024 §2 | « For granite rocks, the average grain size is used to define three rock categories: fine (<1 mm), **medium (1–5 mm)** and coarse (>5 mm) … Kuru Grey, Red Bohus and Sidobre are fine-, medium- and coarse-grained, respectively » ; « Red Bohus granite is a light reddish rock, with interlocking crystal grains … **Feldspars have sizes that can be larger than several millimetres** (Fig. 3b) ; iron oxides coloured feldspar minerals in orange. **No anisotropy is observable.** » ; Fig. 3(d) « focus on a **large quartz crystal** of Red Bohus granite showing truncated failures » ; « grain boundaries … more faceted ones for Sidobre and Red Bohus » (vs imbriquées pour Kuru). Repères des voisins : Kuru ≈ 1 mm toutes phases ; Sidobre : quartz en amas de 5 mm, microcline ~1 cm, plagioclase 2–3 mm, biotite ~3 mm. | [VÉRIFIÉ] |
| Dumoulin 2024 §5.2, Table 6 | « The Red Bohus specimens had a **moderately visible foliation structure** », angle de foliation θ = 31–61° mesuré sur les 12 éprouvettes triaxiales (2_2 … 2_13). ⚠ Contredit le « no anisotropy » du §2 : il existe une orientation préférentielle faible (fluidalité magmatique des feldspaths K, cf. rapakivi.dk : « tydeligt flydemønster med ensretning af de lidt større kalifeldspatter »). | [VÉRIFIÉ] |
| Stenkartoteket 2011 (Brastad, Broberg, Evja) | « Kornstruktur » noté **3 = Medel** sur l'échelle 1 (Fin < 3 mm) … 5 (Grov > 10 mm) ; Brastad : « **Kristallina anhopningar förekommer** » (des amas cristallins se rencontrent) ; « Ytliga korngränser » 2–3 (joints de grains peu à moyennement apparents). | [VÉRIFIÉ] |
| RISE 2024 (Skarstad / Broberg) | Skarstad : « **Medium- to coarse-grained texture** » ; Broberg : « Mediumgrained ». | [VÉRIFIÉ] |
| Schouenborg & Eliasson (SGU) 2015 | « medium grained to coarse grained and even with some **porphyric** parts ». | [VÉRIFIÉ] |
| rapakivi.dk (d'après Asklund 1947, SGU ser. C, 52 p. sur le Bohus) | « De fleste er nogenlunde **enskornede**, men der er også mange **porfyriske** typer » ; « Korngrænserne er skarpe ». | [VÉRIFIÉ, secondaire] |
| Lecomte & Schuman 2026, Constr. Build. Mater. 511, 145169, DOI 10.1016/j.conbuildmat.2026.145169 | « typical **coarse-grained**, polymineralic igneous rock, with **millimetric grains** of quartz, feldspar (pink) and mica » ; cartes EBSD (la zone cartographiée contient ~50 % feldspath, 4 % quartz — zone locale, non représentative). | [VÉRIFIÉ, abstract + extraits ; corps payant] |
| Iddefjord (même batholite, Norvège) | « medium coarse with grain size of 1–5 mm » (page earthcache, sans source) | [NON VÉRIFIÉ] |
| Petersson & Eliasson 1997 (Lithos 42:123) ; Petersson et al. 2014 (Lithos 196-197:99) | Episyénites « typically medium-grained, somewhat **K-feldspar porphyritic** … granite texture partially preserved by the feldspar framework of the protolith » (donc le protolithe est un granite à tendance porphyrique en Kfs). | [VÉRIFIÉ, abstracts + extraits] |
| Åkesson, Hansson & Stigh 2004, Eng. Geol. 72:131, DOI 10.1016/j.enggeo.2003.07.001 | Bohus = « **isotropic granite** » avec « a well-developed existing microcrack pattern » corrélé au joint set horizontal ; nouvelles fissures intragranulaires surtout dans les **feldspaths** ; micas et opaques favorisent la propagation (contraste de E). Pas de tailles de grains dans l'abstract. | [VÉRIFIÉ, abstract] |
| Saadati 2015 thèse, Fig. 4 (CT 20×20 mm) | Sur l'image, taches claires (biotite/opaques) ≈ 1–2 mm, grains gris (feldspath/quartz) ≈ 3–6 mm. | [ESTIMATION visuelle] |

### 3.2 Synthèse texturale (à citer avec les réserves ci-dessus)

Granite **équigranulaire à faiblement porphyrique**, grain **moyen (classe 1–5 mm), localement grossier**, feldspath K (microcline perthitique) = phase la plus grosse (pluri-mm, phénocristaux occasionnels ≈ 1 cm), quartz en grains/amas de quelques mm, plagioclase un peu plus fin que le Kfs, biotite en paillettes ≈ 1 mm ; joints de grains **nets et facettés** (pas d'imbrication type Kuru) ; **fissuration intragranulaire préexistante** (quartz surtout, arrêtée par les feldspaths — Dumoulin Fig. 3d ; Saadati « many pre-existing cracks »), faible fluidalité des Kfs (θ 31–61°).

**Aucune distribution de tailles par minéral n'est publiée** : ni l'article ORCHYD, ni les fiches EN 12407 (qui ne publient que les fractions), ni Åkesson 2004 (abstract) n'en donnent. La colonne « taille » du tableau final est donc une estimation.

---

## 4. Propriétés physiques

| Grandeur | Valeur | Source |
|---|---|---|
| Densité (UCS, Ø50×100) | 2626–2640, **moy. 2634 kg/m³** | Dumoulin 2024 Table 1 [VÉRIFIÉ] |
| Densité (triax, Ø40×100) | 2645–2650 kg/m³ | Dumoulin 2024 Table 6 [VÉRIFIÉ] |
| Densité apparente EN 1936 | 2640 (Brastad 2011), 2620 (Broberg 2011), **2647** (Skarstad 2024), 2662 (Broberg 2024) | fiches [VÉRIFIÉ] |
| Densité (KTH, spalling) | 2660–2671 kg/m³ | Saadati et al. 2016, DOI 10.1155/2016/6279571 [VÉRIFIÉ] |
| ⚠ Dépôt | `targets_redbohus.json:269` = **2620** ; `DOCUMENTATION_rockim.md:647` (carte DP-DFH) = 2620 | → passer à **2640 ± 15** |
| Porosité | « below 0.5 % » (Dumoulin) ; **0,3 % ouverte** EN 1936 (RISE 2024) ; ~0,2 % (Saadati 2015 p. 4) ; absorption d'eau 0,1–0,13 % pds | [VÉRIFIÉ] |
| Vitesse P | **5515 m/s** (moy. UCS), 5331–5497 (triax) ; « no open porosity » d'après la corrélation de del Río | Dumoulin 2024 [VÉRIFIÉ] |
| Vitesse de barre C₀ | 4050 m/s (éprouvette de spalling) | Saadati 2016 [VÉRIFIÉ] |

---

## 5. Données mécaniques publiées

### 5.1 ORCHYD (Dumoulin et al. 2024) — la source des essais de Fernando

| Essai | Résultats bruts | Moyenne |
|---|---|---|
| UCS Ø50×100 mm, ε̇ = 10⁻⁶ s⁻¹ (Mines Paris) | 112,1 / 158,2 / 117,9 / 114,8 MPa | **125,8 MPa** (dispersion « probably due to some defects in the block … or mineral defects ») |
| E local (jauges ~1 cm) | 64,25 / 60,0 / 50,46 / 54,46 GPa | **57,3 GPa** |
| ν local | 0,11 / 0,15 / 0,20 / 0,22 | 0,17 |
| E global (LVDT platines) | 17,3–22,4 GPa (souplesse machine — inutilisable) | — |
| BTS Ø50×25, 200 N/s | 10,2 / 10,6 / 11,4 / 9,0 MPa | **10,3 MPa** (= 10,27 recalculé par Fernando) |
| Triaxial Ø40×100 (SINTEF/NTNU, GCTS RTR-4000), σ₁ pic | σ₃ 20 : 422,4/427,9/424,3 ; 50 : 652,0/647,9/647,1 ; 75 : 777,3/772,3/787,7 ; 100 : 894,7/900,0/903,1 MPa | q = 404,9 / 599,0 / 704,1 / 799,3 MPa (= `targets_redbohus.json:201-206`) |
| Angle de rupture β | 19–30° (≈25° à σ₃ ≥ 50) | — |
| Critère proposé (éq. 5) | σ_C^ref = 424,8 MPa à σ₃ = 20 MPa, α = 0,54 (pression), β = 0,019 (vitesse) | — |
| Classe post-pic | **Class II** à tous confinements ; « influence of confining pressure on axial strength is highest for Red Bohus … confining pressure closes transgranular cracks » | — |

Toutes ces lignes : [VÉRIFIÉ] (texte intégral lu ; Tables 1, 2, 6, 9).

### 5.2 Lignée KTH / Atlas Copco (Bohus « gris » probable, voir §1)

| Grandeur | Valeur | Source |
|---|---|---|
| E (traction et compression directes) | **52 GPa** ; ν = **0,15** | Shariati 2019 §2.2 ; Saadati 2016 Table 1 [VÉRIFIÉ] |
| Résistance en traction directe quasi-statique | ≈ **8 MPa** (éprouvette de spalling) | Saadati 2016 [VÉRIFIÉ] |
| Traction dynamique (spalling, 70 s⁻¹) | **18,9 MPa** (modèle DFH : 19,5) ; ratio dyn/QS ≈ 2,5 | Saadati 2016 [VÉRIFIÉ] |
| Flexion 3 points 40×40×150 | σ_w = **18,7 MPa**, σ_sd = 1,0, V_eff = 189–195 mm³, **m = 23** | Shariati et al. 2022 (arXiv 2201.01870, p. 8-9) ; Saadati 2016 Table 1 [VÉRIFIÉ] |
| m par indentation (12 essais) | m = 24 (3 points faibles écartés ; m = 6 si on les garde — « two populations of defects ») | Shariati 2022 arXiv p. 10 [VÉRIFIÉ] — confirme `DOCUMENTATION_rockim.md:193` |
| Traction uniaxiale vs taille | indépendante du volume ; point-load décroît avec le volume | Wijk, Rehbinder & Lögdström 1978, Rock Mech. 10:201, DOI 10.1007/BF01891959 [VÉRIFIÉ, résumé] |
| Plasticité confinée | surface DP linéaire jusqu'à p ≈ 750 MPa, dilatance ψ(p) décroissante | Shariati 2019 [VÉRIFIÉ] |
| Rigidité de spécimens pré-fissurés | 30–35 GPa au lieu de 52 | Saadati 2016 [VÉRIFIÉ] |

### 5.3 Fiches industrielles (EN 1926 sur cubes, 6 éprouvettes typ.)

| Variété | UCS EN 1926 | Flexion EN 12372 | Source |
|---|---|---|---|
| Röd Bohus Brastad (2011) | **243 MPa** | 17,0 MPa | fiche [VÉRIFIÉ] |
| Röd Bohus Broberg (2011) | 244 MPa | 18,4 MPa | fiche [VÉRIFIÉ] |
| Skarstad Red Bohus (RISE 2024) | 218 ± 22 MPa (borne basse 167) | 14,8 ± 0,9 | fiche [VÉRIFIÉ] |
| Broberg Red Bohus (RISE 2024) | 199 MPa (borne basse **91** !) | 18,2 | fiche [VÉRIFIÉ] |
| Grå Bohus Evja (SP 2011) | 219 MPa | 16,7 | fiche [VÉRIFIÉ] |

Lecture : l'UCS sur cubes (200–245 MPa) est ~1,7× l'UCS ORCHYD sur cylindres 2:1 (126 MPa) — effet de forme + défauts du bloc ORCHYD (Dumoulin : « probably due to some defects in the block »). La borne basse de 91 MPa chez Broberg confirme une population à défauts (cohérent avec le m ≈ 6 « deuxième population » de Shariati et avec les UCS 112–158 de Fernando).

### 5.4 Ce qui n'existe PAS dans la littérature (ou n'a pas été trouvé)

- **K_Ic macroscopique du Bohus** : aucune valeur trouvée. Lecomte & Schuman 2026 mesurent une ténacité **par phase** par nano-indentation cube-corner (corps de l'article payant, valeurs non lues) [NON VÉRIFIÉ]. Pour ordre de grandeur (granite générique, nano-indentation) : minéraux 3,1–6,2 MPa·m^½, interfaces 0,7–4,3 (Yao, Liu, Chen 2026, Eng 7:130, DOI 10.3390/eng7030130 [VÉRIFIÉ, résumé]) — ce ne sont **pas** des K_Ic macro.
- **σ_ci / σ_cd publiés** : aucun pour le Bohus. Les seuls existants sont ceux de Fernando (`CONTINUUM/calib_bohus_triax/exp_qc/seuils_sbm_bohus.json`) : CI ≈ 0,52–0,60 du pic, CD ≈ 0,59–0,75 selon l'essai [VÉRIFIÉ, fichier lu].
- **Modules par minéral mesurés sur le Bohus** : existent (Lecomte & Schuman 2026, EBSD + nano-indentation Berkovich, feldspath avec ISE classique, quartz/mica sans RISE intrinsèque) mais derrière paywall [NON VÉRIFIÉ]. Les valeurs du dépôt (quartz 83,1 / biotite 29,3 GPa, Table 2 d'Aboayanah 2024) restent le choix par défaut.
- **Distribution de tailles de grain par minéral** : aucune (voir §3).

---

## 6. Tableau PRÊT À L'EMPLOI pour rockim (GBM, 2D déformation plane)

### 6.1 Option A — 3 phases (structure actuelle des decks `q3_gbm_*.cfg`)

| Phase rockim | Fraction d'aire | Confiance | Taille moyenne d_eq [mm] | Confiance | Justification |
|---|---|---|---|---|---|
| `feldspar` (Kfs + plagioclase) | **0,62** | **haute** (4 comptages EN 12407 : 59–63 %) | **3,5** (fourchette 2,5–4,5) | **basse** | phase la plus grosse ; « larger than several millimetres » ; classe 1–5 mm ; Kfs porphyrique occasionnel ≈ 10 mm à ne pas modéliser |
| `quartz` | **0,32** | **haute** (30–35 %) | **2,5** (1,5–3,5) | **basse** | « large quartz crystal » avec fissures internes ; grains/amas de quelques mm |
| `biotite` (+ opaques + chlorite) | **0,06** | **haute** (bt 3–6 % + opaques 1–2 %) | **1,0** (0,5–1,5) | **basse** | paillettes mm ; taches CT 1–2 mm ; ⚠ hmin et dt suivent ce grain (`calib_quick/README.md:141-142` : dt ÷ 3,6 à σ = 0,8) — 1,2 mm est un compromis coût acceptable |
| **Global** `grainSize` (moyenne pondérée par l'aire) | — | — | **3,0** | moyenne | 0,62×3,5 + 0,32×2,5 + 0,06×1,0 = 3,0 mm → le « grains 3 mm » de `q3_gbm_P050.cfg` est **physiquement défendable**, pas seulement un choix de coût |
| `grainSizeSpread` (sd de ln d, global) | **0,5** (balayer 0,3 / 0,5 / 0,7) | basse | — | [ESTIMATION] : sd intra-phase ≈ 0,4 + variance inter-phases (moyennes 3,5/2,5/1,0 aux poids 0,62/0,32/0,06) → sd globale ≈ 0,51 ; 0,7–0,8 pour les faciès porphyriques signalés par SGU/Asklund ; c'est le cas c du §5.16 de `DOCUMENTATION_rockim.md` (réalisé 0,491) |

À remplacer dans le dépôt : `phase.biotite.fraction = 0.07` → **0.06**, `phase.quartz.fraction = 0.31` → **0.32** (`gen_gbm.py:49-56`) ; `rho = 2620` → **2640** ; commentaire « grains ~2 mm » → « d_eq ≈ 3 mm (classe 1–5 mm) ».

### 6.2 Option B — 4 phases (si l'on veut séparer les deux feldspaths, comme certains GBM)

| Phase | Fraction d'aire | d_eq [mm] |
|---|---|---|
| Feldspath K (microcline perthitique) | 0,33 | 4,0 |
| Plagioclase (oligoclase–andésine) | 0,29 | 2,5 |
| Quartz | 0,32 | 2,5 |
| Biotite + opaques | 0,06 | 1,0 |

Fractions [VÉRIFIÉ ±3 %] ; tailles [ESTIMATION]. Intérêt limité tant que E_Kfs ≈ E_plag dans la table de phases ; utile surtout si l'on veut affaiblir spécifiquement les joints Qz–Kfs (astuce d'Aboayanah).

### 6.3 Paramètres macroscopiques associés (rappel sourcé)

| Paramètre | Valeur recommandée | Source |
|---|---|---|
| ρ | 2640 kg/m³ | §4 |
| Porosité | 0,3 % (négligeable, pas de « pores » à modéliser) | §4 |
| E apparent cible | 57 GPa (jauges UCS, Dumoulin) vs 77,7 GPa (fit triaxial LVDT de Fernando) vs 52 GPa (KTH) — **le module croît avec la fermeture des fissures préexistantes** ; garder E, ν comme sorties (règle déjà posée dans `targets_redbohus.json:267`) | §5 |
| BTS | 10,3 MPa ; traction directe QS ≈ 8 MPa (ratio DTS/BTS ≈ 0,55–0,8) | §5 |
| Weibull m (traction, méso) | 23–24 ; population secondaire de défauts (m apparent 6) = argument pour `jointPrebrokenFrac` | §5.2 |
| Foliation | faible mais réelle (θ 31–61°) — à ignorer en première calibration, mais à garder en tête si le faciès de rupture est systématiquement incliné | §3.1 |

---

## 7. Ce qui manque et comment l'obtenir à bas coût

1. **Mesure des tailles par minéral (priorité 1, ~1 h)** : photographier une face de carotte ORCHYD (ou utiliser les photos des éprouvettes rompues, Dumoulin Fig. 16) avec une règle, puis comptage d'interceptes linéaires par phase (ImageJ) : d_eq par phase, sd de ln d, fraction d'aire de contrôle. Ceci remplace toute la colonne « [ESTIMATION] » du §6.
2. **Rapport pétrographique EN 12407 complet** : la norme impose une description des tailles de grains ; RISE l'a produit pour Hallindens le 2024-05-14 (fiches Skarstad/Broberg). Demande à faire à Hallindens Granit AB (info@hallindensgranit.se) ou à RISE. [non fait]
3. **Livrable ORCHYD D2.1** (« Project Specifications », Cazenave, Gerbaud, Sellami, Thenevin, Velmurugan 2021, https://www.orchyd.eu/repository/) — cité comme réf. 5 par Dumoulin ; peut contenir la pétrographie détaillée (lames minces Mines Paris). [non consulté]
4. **Zenodo 10.5281/zenodo.10617548** : le zip (Excel + scripts Python) reproduit les figures, dont Fig. 3 (lames minces) — vérifier si des images à l'échelle y sont incluses. [site inaccessible pendant la tâche : 403/504]
5. Lecomte & Schuman 2026 (accès institutionnel PSL) : E, H, K_Ic **par phase du Bohus** — remplacerait la Table 2 d'Aboayanah par des valeurs propres à la roche.
6. Asklund B. (1947) *Svenska stenindustriområden I-II — Gatsten och kantsten*, SGU ser. C : 52 pages sur le Bohus, avec descriptions pétrographiques par carrière ; non numérisé en ligne (à demander à la bibliothèque SGU).

---

## 8. Sources consultées (statut)

- Dumoulin et al. 2024, GETE 40:100592, DOI 10.1016/j.gete.2024.100592 — texte intégral lu via ScienceDirect (CC-BY) [VÉRIFIÉ]
- Zenodo 10.5281/zenodo.10617548 — métadonnées DataCite lues ; fichiers non ouverts
- Hallindens Granit AB, fiches techniques RISE 2024-05-14 Skarstad / Broberg Red Bohus (EN 1926/1936/12372/12407/13755) [VÉRIFIÉ]
- Sveriges Stenindustriförbund, Stenkartoteket 2011 : Röd Bohus Brastad, Röd Bohus Broberg, Grå Bohus Evja (CBI/SP) [VÉRIFIÉ]
- Saadati M. 2015, thèse KTH n° 89 (kappa, 44 p.) [VÉRIFIÉ] ; Saadati et al. 2014 IJNAMG 38:828 DOI 10.1002/nag.2235 [résumé] ; Saadati et al. 2016 AMSE 6279571 DOI 10.1155/2016/6279571 [VÉRIFIÉ] ; Saadati et al. 2018 J. Test. Eval. 46(1) DOI 10.1520/JTE20160072 [non lu]
- Shariati et al. 2019 RMRE 52:645 DOI 10.1007/s00603-018-1646-3 [VÉRIFIÉ, OA] ; Shariati et al. 2022 IJNAMG 46:374 (arXiv 2201.01870) [VÉRIFIÉ] ; Shariati et al. 2022 RMRE 55:7369 DOI 10.1007/s00603-022-02991-9 [résumé + annexes] ; Shariati 2019 lic. thesis KTH [VÉRIFIÉ, kappa]
- Åkesson, Hansson & Stigh 2004 Eng. Geol. 72:131 DOI 10.1016/j.enggeo.2003.07.001 [résumé + extraits] ; Lindqvist, Åkesson & Malaga 2007 Mater. Charact. 58:1183 DOI 10.1016/j.matchar.2007.04.012 [résumé — pas de données Bohus]
- Wijk, Rehbinder & Lögdström 1978 Rock Mech. 10:201 DOI 10.1007/BF01891959 [résumé]
- Schouenborg & Eliasson 2015, EGU2015-14883 [VÉRIFIÉ] ; Eliasson & Schöberg 1991 Precambrian Res. [cité] ; Petersson & Eliasson 1997 Lithos 42:123 DOI 10.1016/S0024-4937(97)00040-6 [résumé] ; Petersson et al. 2014 Lithos 196-197:99 DOI 10.1016/j.lithos.2014.02.025 [résumé + extraits]
- Lecomte & Schuman 2026 CBM 511:145169 DOI 10.1016/j.conbuildmat.2026.145169 [résumé + extraits] ; Yao, Liu & Chen 2026 Eng 7:130 DOI 10.3390/eng7030130 [résumé] ; Peng, Wong & Teh 2017 JGR 122:1054 DOI 10.1002/2016JB013469 [résumé]
- SGU « Landskapsstenar i Götaland » ; sv.wikipedia « Bohusgranit » (« vanligen jämnt småkornig » — contredit par SGU/ORCHYD, faible poids) ; rapakivi.dk/bohusgranit [secondaires]
- Fichiers du dépôt : `rockim_f1/calib_triax3d/gen_gbm.py:8-11,49-62` ; `rockim_f1/calib_triax3d/CAMPAGNE.md:360` ; `rockim_f2/calib_quick/README.md:14,138-142` ; `rockim/rockim_p1/calibration_redbohus/README.md:10-13` ; `…/targets/targets_redbohus.json:2-7,269` ; `…/tools/extract_targets.py:68-73` ; `phd/HETEROGENEITE.md:19` ; `rockim_f2/DOCUMENTATION_rockim.md:193,647,1254-1315` ; `CONTINUUM/calib_bohus_triax/exp_qc/seuils_sbm_bohus.json` ; bibs Overleaf (`refs_dfh.bib` : saadati2014/2015/2018, shariati2022, forquin2010 ; `refs_forage.bib` : saadati2020) — **manquent** : dumoulin2024 (dans `refs_complement.bib:99` sous « Dumoulin & Kane », à vérifier), akesson2004, wijk1978, shariati2019, saadati2016, schouenborg2015, lecomte2026.

Échecs d'accès à signaler : Zenodo (403/504), HAL hal-04700616 (Anubis anti-bot), archive NVA Oslo (PDF 25 Mo réservé), ResearchGate, StoneContact (403), SGU K70 (pas de Bohus), Grokipedia (403).
