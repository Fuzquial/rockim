# Nouveaux bancs d'essai pour rockim (FDEM 2D et 3D, granite sous percussion)

Préparé le 2026-10-06. Ce document propose des bancs qui ne sont pas encore traités dans la
campagne V&V (`docs/VV_campagne.md`, `vv/`) ni dans le rapport-guide (`docs/rapport_guide/sections/`).
Il contient l'inventaire de la recherche, avec les candidats écartés et la raison de chaque rejet,
puis sept fiches classées, puis les decks prêts avec leur commande et leur durée estimée. Les
références exactes sont dans `nouveaux_bancs.bib` (même dossier) ; `rockim.bib` n'est pas modifié.
Aucun calcul complet n'a été lancé : seuls des essais de démarrage de quelques pas ont tourné
(`OMP_NUM_THREADS = 2`, binaire Linux du 2026-10-06 compilé depuis `src/` courant).

Conventions. Toute valeur chiffrée porte sa source : article (tableau, figure, page), deck existant
de rockim, ou calcul fait ici à partir de valeurs sourcées (formule donnée). Une valeur « lue sur
figure » porte son incertitude de lecture. Les critères d'acceptation sont écrits ici, avant tout
calcul ; ils ne seront pas modifiés après coup.

## 1. Méthode de recherche

1. Existant. Lecture de `docs/VV_campagne.md` (P1 à P5, B1 à B9), des `README.md` de `vv/*/`
   (tableau d'inventaire de `vv/B_impact_coupe/README.md` compris), des titres de sections de
   `docs/rapport_guide/sections/*.tex`, de `DOCUMENTATION_rockim.md` (§1, §5.1, §5.4 bis, §5.7, §5.8,
   §5.21) et de `GUIDE_rockim.md`.
2. Bibliographie du doctorant. Première et deuxième pages des 209 PDF de
   `phd_geothermie/bibliographie/` extraites par `pdftotext` (index dans le scratchpad de la
   session, `biblio_txt/`), puis lecture complète des articles candidats.
3. Web. Recherches ciblées (Pekeris, Petersson, Zhou-Molinari-Ramesh, Saadati). Le proxy de la
   session bloque arxiv.org, academic.oup.com, pnas.org et plusieurs sites de documentation : les
   valeurs qui n'ont pu être relues que dans un extrait de moteur de recherche sont signalées comme
   telles, et la bibliographie locale est préférée chaque fois qu'elle contient la source.

Critères de choix, dans l'ordre de la demande : (a) référence analytique ou essai plutôt que code
contre code ; (b) indépendance vis-à-vis du calage de rockim ; (c) jeu de données complet publié ;
(d) faisable sans développement, ou développement chiffré ; (e) coût, cas COURT (moins de 30 min à
4 fils, en cloud) ou LONG (poste local, 14 fils).

## 2. Ce qui existe déjà (à ne pas refaire)

| Domaine | Bancs existants | Où |
|---|---|---|
| Élément, ondes | onde 1D dans une barre confinée (d'Alembert), balayage de pénalité ; patch test et Kirsch 3D ; Kirsch, Lamé, Sneddon en 2D | `vv/V1_onde`, `vv/P1_elastique`, `tunnel_edz/`, `bench_abuaisha/` |
| Contact | bloc glissant (Xiang, Fukuda) ; sphère sur plan (Hertz) ; collisions ; Signorini, chaîne de masses, raclage | `vv/P2_bloc`, `vv/P2_sphere`, r03 |
| Joint, rupture | joint isolé (modes I, II, cycle fermé) ; fissure pressurisée (Guo §2.4) ; fissure en aile (Lisjak, MLER) | r03, `vv/P3_fissure_guo`, `vv/P3_aile_lisjak` |
| Essais de laboratoire | flexion trois points non entaillée et brésilien de Guo ; UCS, triaxial, brésilien, SHPB de Yan 2023 (contre code) ; calibration Red Bohus | `vv/P4_essais_guo`, r05a, r07a |
| Articles | AbuAisha 2017 ; tunnel de Wang 2024 ; St Anne et Kuru (Yang 2025, 2026) ; Saksala 2011 contre Abaqus ; coupe PDC (Heilman 2024, LeBaron 2023) | `vv/B_articles`, `vv/B_impact_coupe` |

Lacunes que ces bancs laissent, d'après `docs/VV_campagne.md` §4 : aucune onde de surface ni bord
absorbant vérifié ; aucun essai dynamique de laboratoire comparé à rockim (« Aucun essai SHPB n'est
comparé à rockim ») ; aucune validation sur des données qui n'ont pas servi au calage ; énergie de
fissuration jamais comparée à Gf multiplié par l'aire rompue ; objectivité de la fragmentation non
quantifiée.

## 3. Inventaire des candidats

### 3.1 Retenus (classés)

| Rang | Banc | Référence | Type de référence | Brique testée | Développement | Coût |
|---|---|---|---|---|---|---|
| 1 | B10 écaillage du granite de Bohus | Saadati et al. 2016 | essai (vélocimétrie laser, jauges) | traction dynamique : DIF, loi cohésive, ondes | aucun | COURT |
| 2 | B11 brésilien dynamique du granite de Kuru | Saksala et al. 2013 ; Wessling et Kajberg 2022 ; Padmanabha et al. 2023 | essais à quatre vitesses + un second laboratoire | effet de vitesse en traction indirecte, montage SHPB 2D | aucun | COURT (par run) |
| 3 | P1.4 problème de Lamb | Lamb 1904, Pekeris 1955, Boussinesq | analytique | ondes P, S, Rayleigh en 3D, bords absorbants | aucun | COURT (fem3d), LONG (fdem3d) |
| 4 | P3.4 anneau en expansion | Levy et Molinari 2010 (loi de Zhou, Molinari et Ramesh) ; Grady 1982 | loi maîtresse numérique et analytique | fragmentation dynamique, convergence cohésive (Camacho-Ortiz) | vitesse initiale radiale, 0,5 à 1 jour | COURT |
| 5 | P3.5 poutre entaillée de Petersson (et extension à l'effet d'échelle) | Petersson 1981 ; Wendner et al. 2015 | essai + identité énergétique exacte | énergie de rupture Gf, objectivité | maillage entaillé (script) | LONG |
| 6 | B12 taillant à trois boutons sur granite de Kuru | Saksala et al. 2014 | essai (force-pénétration, cratères) | percussion multi-boutons, train de tiges | taillant maillé (préparation, pas de code) | LONG |
| 7 | B13 coalescence de deux fissures (gypse) | Bobet et Einstein 1998 | essai | amorçage et coalescence en compression | aucun (2D, `preBrokenJoints`) | LONG |

### 3.2 Écartés

| Candidat | Raison principale |
|---|---|
| Essais de Mazars, de Brokenshire (torsion d'une poutre entaillée), de Nooru-Mohamed (traction-cisaillement), panneau en L de Winkler | béton, quasi statique 3D long en explicite ; le pilotage biaxial de Nooru-Mohamed n'existe pas dans rockim ; rien de propre au granite que la poutre de Petersson (rang 5) ne teste déjà |
| Kalthoff et Winkler (plaque d'acier entaillée impactée) | acier maraging, référence d'essai réduite à un angle de branchement ; maillage 3D fin coûteux ; aucune valeur vérifiable dans la bibliographie locale (Kalthoff 2000 absent) |
| Cônes de Hertz dans le verre (Chaudhri, Lawn, loi d'Auerbach) | verre ; la loi d'Auerbach traduit une statistique de défauts de surface que rockim ne porte pas ; cône à résoudre sous le rayon de contact (maillage très fin) ; sources absentes du dépôt |
| Plaques à trous multiples (Haeri) | comparaison de trajets de fissures sur photographies, sans grandeur chiffrée |
| Klepaczko et Brara 2001 ; Erzar et Forquin 2010 (écaillage du béton) | même essai que B10, sur béton ; Saadati et al. 2016 fait le même essai, dans le même laboratoire (3SR), sur granite de Bohus |
| Saadati et al. 2014, impact sur la tranche (EOI) | comparaison qualitative des densités de fissures ; fissures structurelles de fabrication non quantifiées (texte de l'article) |
| Mahabadi et al. 2012 (Y-Geo, brésilien à microstructure réelle) | exige la carte minéralogique de l'éprouvette ; comparaison code contre essai déjà couverte par Yan 2023 en GBM |
| Bancs publiés de HOSS et d'Irazu | Knight et al. 2020 constatent eux-mêmes l'absence de bancs standard ; l'élément de Lei et Rougier 2016 est déjà repris en P1 |
| Camacho et Ortiz 1996 seul | fondu dans P3.4 : l'anneau est le test de convergence cohésive de référence (Molinari et al. 2007) |
| Wu et al. 2021 (indentation dynamique confinée du granite) | confinement biaxial vrai ; fdem3d n'a qu'une pression scalaire (`confiningPressure`), limite déjà notée pour Heilman (B9) |
| Jiang et al. 2020 (FDEM d'indentation) | code contre code (joints bilinéaires d'Abaqus) |
| Zhang et al. 2022 (impact de cylindres de granite, PLOS ONE) | seuls UCS, E, ν et ρ publiés (§2) : ft et Gf à caler sur la même campagne ; références surtout qualitatives (angle du plan de rupture, distributions de fragments) ; à garder en réserve, données ouvertes (figshare) |
| Wendner et al. 2015 seul (béton, quatre tailles) | conservé comme extension du rang 5 (effet d'échelle), pas comme banc séparé |
| Liao et al. 2021 (dynamique d'un indenteur vibrant) | dynamique non linéaire du système de forage, pas de fissuration |
| Lundberg 1973 (transfert d'énergie en percussion) | théorie du train de tiges ; le contact d'outil et la chaîne de masses sont déjà vérifiés (r03) |
| Saksala 2011, Aising et Yang, Heilman et LeBaron, Guo, Lisjak, AbuAisha, Wang, Yan | déjà traités (section 2) |

## 4. Fiches

### Banc 1 (B10) : écaillage d'un barreau de granite de Bohus

Référence. M. Saadati, P. Forquin, K. Weddfelt, P.-L. Larsson, On the Tensile Strength of Granite
at High Strain Rates considering the Influence from Preexisting Cracks, Advances in Materials
Science and Engineering 2016, 6279571, doi:10.1155/2016/6279571 (PDF local). Paramètres
complémentaires : Saadati et al. 2014, IJNAMG (preprint local `IJNAMG2014-ccsd.pdf`, Table 1).

Objet. Un projectile d'aluminium frappe une barre de Hopkinson collée au barreau de granite ;
l'onde de compression se réfléchit en traction sur la face arrière libre et écaille le barreau.
Un vélocimètre laser mesure la vitesse de la face arrière ; la chute après le premier maximum
(vitesse de « pullback » ΔV_pb) donne la résistance dynamique par la formule de Novikov,
σ_dyn = ½ ρ C₀ ΔV_pb (éq. 1 de l'article).

Ce que le banc teste dans rockim. La traction dynamique pure à ε̇ ≈ 70 s⁻¹, le régime de la
percussion selon l'article (§1, §2.1) : insertion adaptative en traction, loi cohésive, facteur
dynamique (`strainRateDIF`), propagation d'onde dans un barreau à faces libres. C'est le premier
banc de rockim qui compare une résistance dynamique à une mesure.

Indépendance. Les paramètres de rupture du deck viennent de l'article (ft quasi statique) et d'un
jeu de rockim (Gf) qui n'a pas été calé sur un essai dynamique. Réserve : le calage Red Bohus de
rockim vise σ_t = 18,3 MPa (« Saadati/Shariati QS », en-tête de
`configs/calib/R_tens.cfg`), voisin des 18,9 MPa de l'écaillage ; le deck B10 n'utilise pas ce
jeu (ft = 8 MPa de l'article), mais il faudra dire lequel des deux jeux sert à la percussion.

Géométrie. Barreau à section carrée de 39,9 × 39,9 mm (aire d'un disque de Ø45 mm), longueur
140 mm. Hypothèse : l'article ne cote pas l'éprouvette (renvoi à Saadati 2015, non disponible) ;
140 mm est compatible avec la jauge G1 placée à 116 mm de la face arrière (§2.2), et Ø45 mm est
l'ordre de grandeur des éprouvettes d'écaillage du même laboratoire (Erzar et Forquin 2010), à
confirmer. La section carrée permet le générateur `box3dbc` existant ; à aire égale, la célérité
de barre et l'impédance sont les mêmes.

Maillage. Tétraèdres non structurés, h = 2,5 mm, 62 217 éléments
(`tools/make_unstructured_mesh.py box3dbc 0.0399 0.0399 0.14 0.0025 …`). Longueur de zone
cohésive ℓ_cz = E G_f / f_t² = 43,6 × 10⁹ × 70 / (8 × 10⁶)² = 48 mm, donc h < ℓ_cz/2 avec une
marge de 10. Pas de mass scaling, `dampingLocal = 0`.

Matériau et provenance.

| Grandeur | Valeur | Source |
|---|---|---|
| ρ | 2 660 kg/m³ | Saadati 2016, texte sous l'éq. 1 |
| C₀ (célérité de barre) | 4 050 m/s | idem |
| E | 43,63 GPa = ρ C₀² | calculé, pour que le barreau propage à C₀ |
| ν | 0,15 | Saadati 2014, Table 1 |
| f_t quasi statique | 8 MPa | Saadati 2016, §2.2 (« same size », renvoi [15]) |
| G_f | 70 J/m² | non publié ; jeu Red Bohus de rockim (`configs/calib/R_tens.cfg`) ; sensibilité 35 et 140 J/m² |
| c, φ | 25 MPa, 40° | défauts de rockim, sans effet attendu en traction uniaxiale |
| DIF | `yang-fig2`, figé à l'insertion | DOCUMENTATION §5.4 bis ; variante sans DIF |

Chargement et conditions aux limites. Traction morte sur la face d'entrée égale à la contrainte
incidente de la jauge G1 (fig. 3 de l'article, courbe pleine, numérisée à ±1,5 MPa et ±1 µs) : pied
à t = 0, pic de 41,5 MPa à 35 µs, retour à zéro à 55 µs (table complète dans le deck). Faces
latérales et face arrière libres. Contrôle de cohérence de la lecture : la vitesse élastique de la
face arrière vaut 2σ/(ρC₀) = 2 × 41,5 × 10⁶ / (2 660 × 4 050) = 7,70 m/s, contre un pic mesuré de
7,75 m/s (fig. 2a, lu à ±0,1 m/s), soit 0,6 % d'écart.

Solution de référence (essai, fig. 2a, lue à ±0,1 m/s et ±1 µs). Premier maximum 7,75 m/s à
79 µs (horloge de la figure) ; minimum suivant 4,25 m/s vers 97 µs ; rebond 4,95 m/s vers 107 µs ;
ΔV_pb = 3,5 m/s (texte) ; σ_dyn = 18,9 MPa à ε̇ ≈ 70 s⁻¹ (texte). Épaisseur de l'écaille déduite ici
de la période de rebond, h_s = C₀ Δt/2 avec Δt = 107 − 79 = 28 µs entre maxima : h_s = 57 mm
(±4 mm pour ±2 µs de lecture). Le modèle DFH de l'article prédit 19,5 MPa (code, non retenu comme
référence).

Grandeurs à extraire. Vitesse de la face arrière v(t) = dU_top_z/dt et dU_rear_z/dt (sommet le
plus proche du centre, à 0,18 mm) ; premier maximum, minimum, ΔV_pb, σ_Novikov rockim ; position
du plan d'écaillage (trames VTU, joints rompus) ; déformations aux jauges g1, g2, g3 (différences
de déplacements) contre la fig. 2b ; résidu B4 ; joints rompus en traction et en cisaillement.

Critères d'acceptation (variante `dif`, maillage de 2,5 mm).

1. Pic élastique de la vitesse arrière : 7,75 m/s ± 5 % (les deux variantes).
2. Temps de montée de la vitesse arrière, de 10 % au pic : 22 µs ± 3 µs (fig. 2a : 57 à 79 µs).
3. ΔV_pb = 3,5 m/s ± 20 %, soit σ_Novikov entre 15,1 et 22,7 MPa.
4. Plan d'écaillage à 57 ± 10 mm de la face arrière.
5. Résidu B4 sous 0,1 % de l'échelle.

La variante `nodif` n'a pas de verdict : elle mesure la part de l'inertie et de la loi cohésive
seules. Ordre de grandeur attendu, non critère : le DIF de Yang vaut au plus 1,85 (borne de
l'éq. `eq-dif` de r02) ; avec f_t = 8 MPa, le DIF seul ne porte la résistance qu'à 14,8 MPa.

Coût. 62 217 tétraèdres, dt = 6,61 × 10⁻⁹ s, 22 700 pas pour T = 150 µs. Essai de démarrage :
303 pas en 24,5 s à 2 fils (initialisation comprise), soit environ 25 min à 2 fils et 15 min à
4 fils par variante, avant le surcoût des joints insérés. COURT.

Réserves. Géométrie de l'éprouvette supposée ; G_f non publié ; la traction imposée suppose que
l'onde de retour n'atteint pas la face d'entrée avant l'écaillage (2L/C₀ = 69 µs après le passage
du pic, l'écaillage se produit au premier retour) ; le collage barre-éprouvette n'est pas
modélisé ; un seul essai publié (dispersion inconnue).

### Banc 2 (B11) : brésilien dynamique du granite de Kuru

Références. T. Saksala, M. Hokka, V.-T. Kuokkala, J. Mäkinen, Numerical modeling and
experimentation of dynamic Brazilian disc test on Kuru granite, IJRMMS 59 (2013) 128-138,
doi:10.1016/j.ijrmms.2012.12.018 (PDF local) : référence principale. A. Wessling, J. Kajberg,
Materials 15 (2022) 8264, doi:10.3390/ma15228264 (PDF local) : second laboratoire, même granite.
V. Padmanabha et al., RMRE 56 (2023) 109-128, doi:10.1007/s00603-022-03075-4 (PDF local) :
dix essais sur granite de Malsburg, référence secondaire de tendance.

Objet. Résistance en traction indirecte d'un disque comprimé entre deux barres de Hopkinson, à
quatre vitesses de percuteur (5, 10, 15, 20 m/s), et sa dépendance à la vitesse de chargement.

Ce que le banc teste dans rockim. Le scénario `shpb` du solveur fdem 2D (barres, disque, jauges)
et l'effet de vitesse en traction (`strainRateDIF`, inertie, loi cohésive), avec le jeu Kuru de
rockim. Premier banc qui compare rockim à un essai SHPB (lacune de VV_campagne §3.3).

Indépendance. Le jeu Kuru de rockim vient de Yang et al. 2026, Table 1 (repris de
`configs/yang2026_kuru_train1_v5.cfg`), calé par ces auteurs sur d'autres essais. Le banc compare
de plus le rapport σ_dyn/σ_QS, tous deux calculés par rockim, au rapport expérimental, ce qui
retire le biais de f_t. Réserve à lever avant calcul : vérifier que la figure 2b de Yang 2025, d'où
vient l'exposant `yang-fig2`, n'inclut pas les essais de Saksala 2013 ou de Wessling 2022 ; sinon
la variante `dif` n'est plus indépendante.

Géométrie (Saksala 2013, §2.3 et Table 1). Disque D = 40,8 mm, épaisseur L = 16 mm ; barres
d'acier Ø22 mm de 1,2 m, percuteur de 0,4 m, amortisseur de caoutchouc. Dans rockim (2D, déformation
plane) : barre incidente de 1,6 m (pour que l'impulsion la plus longue, 280 µs × 5 048 m/s =
1,41 m, ne se superpose pas à sa réflexion sur la jauge M1 placée à 0,8 m du disque), barre
transmise de 1,2 m.

Maillage. Disque `shpbDiscElemSize` = 0,5 mm, barres 2 mm (73 521 éléments au total, hmin 0,19 mm
au lissage de la jante). ℓ_cz = 60 × 10⁹ × 50 / (10,98 × 10⁶)² = 24,9 mm, h < ℓ_cz/2 avec une marge
de 25. Étape quasi statique : disque Gmsh Delaunay de 14 151 triangles, méplats de 2 × 5°
(`gen_mesh.py`, recette de `calib_quick/make_disc_mesh.py`).

Matériau. ρ = 2 626 kg/m³, E = 60 GPa, ν = 0,24, f_t = 10,98 MPa, c = 29,84 MPa,
tan φ = 1,85, G_I = 50 J/m², G_II = 1 000 J/m² (Yang 2026, Table 1). Barres : ρ = 7 850 kg/m³,
E = 200 GPa, ν = 0,29 (acier de Yang 2026 ; Saksala 2013 écrit seulement « high strength steel »).

Chargement. Vitesse imposée à l'extrémité de la barre incidente, trapèze de rockim
(`shpbPulse = trapezoid`) ajusté sur les approximations linéaires de la fig. 2a de Saksala 2013
(lecture ±10 µs, ±10 MPa), V₀ = σ_inc/(ρ_b c_b) avec ρ_b c_b = 39,62 MPa·s/m :

| Deck | Essai de Saksala | σ_inc (MPa) | Montée, palier, descente (µs) | V₀ (m/s) | τ (µs) | Palier |
|---|---|---|---|---|---|---|
| `shpb_v05_*` | 1 (5 m/s) | 115 | 120, 45, 115 | 2,902 | 280 | 0,161 |
| `shpb_v10_*` | 6 (10 m/s) | 225 | 55, 150, 35 | 5,678 | 240 | 0,625 |
| `shpb_v15_*` | 11 (15 m/s) | 285 | 45, 135, 40 | 7,193 | 220 | 0,614 |
| `shpb_v20_*` | 16 (20 m/s) | 400 | 35, 140, 30 | 10,095 | 205 | 0,683 |

Le trapèze de rockim est symétrique : montée et descente sont prises égales à leur moyenne. La
comparaison se fait à la vitesse de chargement mesurée dans le calcul ; une erreur sur l'impulsion
déplace le point le long de la droite de référence sans fausser le verdict.

Solution de référence. Saksala 2013, fig. 5 : ajustement linéaire imprimé σ_t = 0,017 Ṡ + 21
(σ_t en MPa, Ṡ en GPa/s), sur les moyennes lues à environ 27,7 ; 37,0 ; 44,8 et 51,5 MPa pour
Ṡ ≈ 420, 1 000, 1 370 et 1 850 GPa/s (lecture ±1 MPa, ±30 GPa/s) ; vitesses de déformation 6,7 ;
16,7 ; 23,3 et 30 s⁻¹ (§2.4) ; f_t0 quasi statique ≈ 13 MPa (§2.2). Dispersion des essais à
5 m/s : de 25,5 à 29 MPa environ (fig. 5, cercles ; l'étoile à 23,5 MPa est la simulation de
Saksala). Wessling 2022, Table 5 : quasi statique 16,9 MPa ;
dynamique 38,1 ± 0,9 MPa au pic, 26,6 ± 3,2 MPa à l'amorçage probable (surcharge de 30 %),
ε̇ = 115 ± 34 s⁻¹, Ṡ ≈ 1 450 GPa/s (texte, pente de l'onde transmise), disque Ø25 × 12,5 mm. À 1 450 GPa/s, l'ajustement de
Saksala donne 45,7 MPa : les deux laboratoires diffèrent de 17 % sur le même granite, ce qui fixe
l'ordre de l'incertitude inter-laboratoires. Padmanabha 2023, Table 2 (granite de Malsburg,
f_t QS = 8,38 ± 2,5 MPa, Table 1) : dix essais de 16,7 à 48,5 MPa entre 2,7 et 27,1 s⁻¹.

Grandeurs à extraire. σ_t(t) = 2 F₂/(π D) avec F₂ = E_b ε_M2 × 0,022 m (force transmise par
mètre d'épaisseur, colonne `epsM2`) et contrainte au centre `sxxC` ; pic σ_peak ; Ṡ = pente entre
20 et 80 % du pic ; instant d'amorçage (premier joint rompu, `nBroken`) ; faciès (fissure
diamétrale, écrasement aux appuis) ; σ_QS de l'étape 0 (`sigmaTpeak`).

Critères d'acceptation (variante `dif`).

1. Pour chacun des quatre decks : |σ_peak/(0,017 Ṡ + 21) − 1| ≤ 15 %, Ṡ mesuré dans le calcul.
2. Pente dσ_peak/dṠ des quatre points entre 0,0085 et 0,034 MPa/(GPa/s) (facteur 2 autour de
   0,017).
3. Rapport σ_peak(v20)/σ_QS,rockim entre 3,0 et 5,0 (essai : 51,5/13 = 3,96).
4. Fissure amorcée au centre du disque (premier joint rompu à moins de D/4 du centre).
5. Équilibre des forces avant le pic : |F₁ − F₂|/F₂ ≤ 10 % sur la montée (l'équilibre est la
   condition de validité du SHPB rappelée par Padmanabha 2023, §2.1 ; le seuil de 10 % est fixé ici).

Variante `nodif` : sans verdict, elle sépare l'effet structurel (inertie, surcharge) de l'effet
matériau. Wessling 2022 et Padmanabha 2023 : comparaison informative, pas de verdict (autres
éprouvettes, autre granite pour Padmanabha).

Coût. 73 521 éléments, dt = 3,0 × 10⁻⁹ s, 224 000 à 249 000 pas selon le deck (T de 672 à
747 µs). Essai de démarrage : 6 663 pas en 65 s à 2 fils, soit environ 40 min à 2 fils et 20 à
25 min à 4 fils par run. COURT par run, huit runs (quatre vitesses × deux variantes). Étape 0
quasi statique : 14 151 triangles, dt = 2,98 × 10⁻⁹ s, au plus 1,34 million de pas (arrêt après le
pic), 2,7 ms par pas à 2 fils : une heure au plus à 2 fils, à classer en cloud long ou poste local.

Réserves. Déformation plane contre un disque de L/D = 0,39 (Saksala 2013 travaille en contrainte
plane et valide ce choix en 3D, §4.3) ; propriétés des barres supposées ; le DIF de rockim est
plafonné à 1,85 alors que le rapport expérimental atteint 3,96 : le banc dira si l'inertie et la
surcharge comblent l'écart, ce qui est précisément ce qu'il doit établir ; frottement barre-disque
non publié (0,1 pris).

### Banc 3 (P1.4) : problème de Lamb, onde de surface et bords absorbants

Références. H. Lamb, Phil. Trans. R. Soc. A 203 (1904) 1-42 ; C. L. Pekeris, The seismic surface
pulse, PNAS 41 (1955) 469-480, doi:10.1073/pnas.41.7.469 (forme fermée pour ν = 1/4) ;
J. D. Achenbach, Wave Propagation in Elastic Solids (1973) pour l'équation de Rayleigh ;
K. L. Johnson, Contact Mechanics (1985) pour la solution de Boussinesq.

Objet. Force ponctuelle verticale P en échelon (montée en cosinus de 2 µs) au centre de la surface
d'un bloc élastique ; déplacement vertical de la surface à quatre distances.

Ce que le banc teste dans rockim. Les ondes P, S et de Rayleigh en 3D (l'onde de Rayleigh n'est
vérifiée nulle part), les charges ponctuelles par groupes (§5.21) et les bords absorbants de Lysmer
avec ressorts (`absorbing = all`, §5.8), qu'aucun banc ne vérifie (VV_campagne §4). En percussion,
l'onde de Rayleigh porte une part importante de l'énergie rayonnée par l'insert.

Géométrie et maillage. Bloc 240 × 240 × 120 mm, source au centre de la face supérieure ; stations
de surface à r = 20, 30, 40, 50 mm sur l'axe x et à r = 40 mm sur la diagonale (contrôle
d'isotropie). Tétraèdres non structurés gradués : 1,5 mm jusqu'à 58 mm de la source, 8 mm au-delà
de 110 mm, source et stations incluses comme sommets exacts (`lamb_pekeris/gen_mesh.py`, 686 392
tétraèdres, 1,5 min de maillage). La longueur d'onde de Rayleigh associée à la montée de 2 µs vaut
5,5 mm (2 766,9 m/s × 2 µs), soit 3,7 éléments.

Matériau. ρ = 2 650 kg/m³, E = 60 GPa, ν = 0,25 (matériau de `vv/V1_onde`). En fdem3d, f_t et c
portés à 10¹² Pa pour rester continu sous la singularité de la force ponctuelle (aucun joint
inséré, comme la variante adaptative de V1).

Solution de référence. Célérités, avec μ = E/(2(1+ν)) = 24 GPa : c_S = √(μ/ρ) = 3 009,4 m/s ;
c_P = √3 c_S = 5 212,5 m/s (ν = 1/4) ; c_R/c_S = √(2 − 2/√3) = 0,919 40 (racine de l'équation de
Rayleigh pour ν = 1/4, Achenbach 1973), c_R = 2 766,9 m/s. Arrivées aux stations : onde P à 3,84 ;
5,76 ; 7,67 ; 9,59 µs, onde de Rayleigh à 7,23 ; 10,84 ; 14,46 ; 18,07 µs pour r = 20, 30, 40,
50 mm (plus le demi-temps de montée). Déplacement statique de Boussinesq w = (1 − ν) P/(2π μ r)
(Johnson 1985, §3.2) : 248,7 ; 165,8 ; 124,3 ; 99,5 nm pour P = 1 kN. Pekeris 1955 donne la
réponse complète pour ν = 1/4 ; d'après la littérature secondaire consultée, le déplacement de
surface y rejoint la valeur statique après le passage de l'onde de Rayleigh. Cette propriété et la
forme fermée doivent être transcrites de l'article original avant le calcul (le PDF n'a pas pu être
téléchargé) ; le critère 3 en dépend. Fenêtre sans réflexion : la première onde réfléchie (P sur
un côté) atteint la station de 50 mm après 190 mm de trajet, soit 36,4 µs ; le fond, après 245 mm
(47 µs). Fenêtre de Lamb : 0 à 35 µs ; fenêtre des bords : 35 à 80 µs.

Grandeurs à extraire. U_r20_z … U_r50_z, U_d40_z (history.csv) ; instants de première arrivée
(seuil 2 % de w_Boussinesq) et du minimum de w ; moyenne de w sur [t_R + 5 µs, 34 µs] ; même
moyenne sur [40, 80] µs ; postes Lysmer et ressorts du bilan.

Critères d'acceptation (fem3d, maillage de 1,5 mm).

1. c_P mesurée (pente des arrivées sur les quatre stations de l'axe) à 2 % de 5 212,5 m/s.
2. c_R mesurée (pente des instants du minimum de w) à 2 % de 2 766,9 m/s.
3. Plateau après Rayleigh à 5 % de Boussinesq aux stations de 30, 40 et 50 mm (sous réserve de la
   transcription de Pekeris, ci-dessus).
4. Isotropie : écart L² relatif entre U_d40_z et U_r40_z sur 0-35 µs inférieur à 2 %.
5. Bords absorbants : sur 40-80 µs, |w − w_Boussinesq|/w_Boussinesq ≤ 10 % à r = 40 mm avec
   `absorbing = all` ; le témoin `lamb_fem3d_noabs` (côtés libres, fond encastré) est rapporté sans
   verdict.
6. fdem3d adaptatif sans joint inséré : mêmes signaux que fem3d à 10⁻³ près en relatif (V1 a
   mesuré 10⁻⁴), résidu B4 sous 10⁻² % ; fem3d : résidu sous 10⁻⁶ %.

Coût. fem3d : dt = 1,96 × 10⁻⁸ s, 4 090 pas, 0,13 s par pas à 2 fils mesuré sur 154 pas : 9 min à
2 fils, environ 5 min à 4 fils (COURT). fdem3d : 2,75 millions de nœuds, dt = 4,47 × 10⁻⁹ s,
17 900 pas, de l'ordre de 0,5 s par pas à 2 fils d'après les 68 pas de l'essai de démarrage
(initialisation comprise) : plusieurs heures à 2 fils, LONG (poste local). Le maillage (34 Mo) n'est
pas laissé dans le dépôt : il se régénère par `gen_mesh.py`.

Réserves. Le bloc n'est pas un demi-espace : seule la fenêtre de 0 à 35 µs est exacte ; la force
ponctuelle sur un sommet régularise la singularité à l'échelle de h ; les ressorts de Deeks et
Randolph (défaut `absorbSpringFactor = 1`) ne restituent pas exactement le champ statique de
Boussinesq, d'où la tolérance de 10 % du critère 5.

### Banc 4 (P3.4) : anneau en expansion, nombre de fragments

Références. S. Levy, J.-F. Molinari, Dynamic fragmentation of ceramics, signature of defects and
scaling of fragment sizes, JMPS 58 (2010) 12-26, doi:10.1016/j.jmps.2009.09.002 (PDF local) ;
F. Zhou, J.-F. Molinari, K. T. Ramesh, IJSS 42 (2005) 5181-5207 et IJF 139 (2006) 169-196 ;
D. E. Grady, J. Appl. Phys. 53 (1982) 322-325 ; G. T. Camacho, M. Ortiz, IJSS 33 (1996) 2899-2938 ;
J.-F. Molinari et al., IJNME (2007), doi:10.1002/nme.1777 (convergence en énergie).

Objet. Un anneau mince reçoit une vitesse radiale initiale v = ε̇ R ; il se fragmente en traction
circonférentielle uniforme. On compte les fragments en fonction de ε̇.

Ce que le banc teste dans rockim. L'insertion adaptative cohésive en régime de fragmentation
multiple, l'objectivité vis-à-vis du maillage et la conservation de l'énergie de rupture (énergie
dissipée = 2 A G_c par fragment créé, à la limite du maillage fin). C'est le test de convergence
de la méthode de Camacho et Ortiz que demande VV_campagne §4, lacune 5.

Développement nécessaire. Une clé de vitesse initiale radiale par groupe, par exemple
`groupVelRadial.<nom> = ε̇ cx cy`, posant v = ε̇ (x − cx, y − cy, 0) sur les nœuds du corps, sur le
modèle de `groupVel.<nom>` (`src/Fdem3dSolver.cpp`, boucle « vitesse initiale par groupe », une
vingtaine de lignes). Chiffrage : 0,5 à 1 jour avec un repère de suite (énergie cinétique initiale
exacte ½ ρ V ε̇² ⟨r²⟩ et bilan B4). Le mode fdem 2D n'a ni groupes ni scénario libre (§5.21, « Ce
qui n'est pas fait ») : le banc se fait en fdem3d, anneau mince.

Géométrie et maillage. Paramètres de Levy et Molinari (§2.1) : circonférence 50 mm (R = 7,96 mm),
section 0,3 × 0,3 mm (choix ici, anneau mince), E = 275 GPa, ν = 0,3, ρ = 2 750 kg/m³,
c = √(E/ρ) = 10 000 m/s, σ_c = 300 MPa, G_c = 100 N/m (§2.3). Échelles (éq. 12 à 14) :
t₀ = E G_c/(σ_c² c) = 30,6 ns ; s₀ = c t₀ = 0,306 mm ; ε̇₀ = σ_c/(E t₀) = 3,57 × 10⁴ s⁻¹.
Maillage : h ≤ s₀/2 = 0,15 mm le long de l'anneau (ℓ_cz = E G_c/σ_c² = s₀), soit environ 8 000
tétraèdres ; trois réalisations de maillage par vitesse ; un tirage de Weibull faible
(`fdem3d_tension_weibull.cfg` pour la syntaxe) en variante.

Chargement. Vitesse radiale initiale, aucune condition aux limites. ε̇/ε̇₀ = 0,5 ; 1 ; 2 ; 5 ; 10 ;
20 ; 50 ; 100.

Solution de référence. Taille moyenne normalisée s̄ = s/s₀ en fonction de ε̄̇ = ε̇/ε̇₀ :
loi de Zhou, Molinari et Ramesh, s̄ = 4,5/(1 + 4,5 ε̄̇^(2/3)) (Levy et Molinari 2010, éq. 18) ;
Grady, s̄ = (24/ε̄̇²)^(1/3) (éq. 16). Le nombre de fragments vaut N = L/(s̄ s₀).

Grandeurs à extraire. Nombre de fragments `nFrag` à convergence (plus aucun joint qui s'endommage),
énergie de rupture dissipée (`eGc`) rapportée à 2 × aire de section × G_c × (N − 1), énergie
cinétique initiale.

Critères d'acceptation. Pour ε̄̇ ≥ 1 : rapport s̄_rockim/s̄_ZMR entre 0,67 et 1,5 (moyenne des trois
maillages) ; pente de log s̄ contre log ε̄̇ entre −0,75 et −0,55 pour ε̄̇ ≥ 10 (asymptote −2/3 de
Zhou et de Grady) ; énergie dissipée à 10 % de 2 A G_c (N − 1) sur le maillage le plus fin ;
écart du nombre de fragments entre les deux maillages les plus fins inférieur à 10 %.

Coût. Environ 8 000 tétraèdres, dt de l'ordre de 10⁻⁹ s, quelques microsecondes simulées : quelques
minutes par run (estimation, pas d'essai de démarrage possible sans la clé). COURT.

Réserves. Matériau céramique fictif (choisi pour comparer à la loi publiée telle quelle) : le banc
vérifie la méthode, pas le granite ; une version granite se déduit par les mêmes échelles (avec le
jeu Kuru, s₀ = 25 mm, t₀ = 5,2 µs et ε̇₀ ≈ 35 s⁻¹, calcul fait ici avec E = 60 GPa, G = 50 J/m²,
f_t = 10,98 MPa, c = √(E/ρ) = 4 780 m/s). La loi de Zhou concerne le cas homogène ; la perturbation du maillage non structuré
suffit à briser la symétrie, à vérifier.

### Banc 5 (P3.5) : poutre entaillée de Petersson, énergie de rupture

Références. P.-E. Petersson, Crack Growth and Development of Fracture Zones in Plain Concrete and
Similar Materials, rapport TVBM-1006, Lund (1981) ; A. Hillerborg, M. Modéer, P.-E. Petersson,
CCR 6 (1976) 773-781. Extension : R. Wendner et al., Materials and Structures 48 (2015)
3603-3626, doi:10.1617/s11527-014-0426-0 (PDF local).

Objet. Flexion trois points d'une poutre de béton entaillée à mi-hauteur, courbe force-flèche
jusqu'à la séparation complète.

Ce que le banc teste dans rockim. L'énergie de rupture : le travail extérieur total à séparation
complète vaut G_f b (D − a), identité exacte pour une loi cohésive (lacune 3 de VV_campagne §4) ;
la position et le niveau du pic contre l'essai ; l'objectivité en maillage. P4.1 (Guo) teste une
poutre non entaillée, sans référence en énergie.

Géométrie et matériau (recherche web du 2026-10-06, à confirmer sur la source avant calcul). Poutre
2 000 × 200 × 50 mm, entaille de 100 mm ; E = 30 GPa, ν = 0,2, f_t = 3,33 MPa, G_f = 124 N/m. La
portée n'a pas pu être confirmée. Longueur caractéristique ℓ_ch = E G_f/f_t² = 335 mm : h = 5 mm
sur le ligament et 20 mm ailleurs respecte largement h < ℓ_cz/2.

Chargement. fdem3d `scenario = loads` : vitesse imposée sur une bande `box.` au-dessus du ligament,
appuis `fix.` sur des lignes de sommets (groupes physiques de courbe du maillage). Durée de
chargement à fixer pour rester quasi statique (au moins 20 allers-retours d'onde dans la poutre),
avec `stopPeakDrop` inopérant ici (scénario `loads`) : arrêt par T.

Solution de référence. Identité énergétique : G_f b (D − a) = 124 × 0,05 × 0,1 = 0,62 J. Courbe
d'essai de Petersson à numériser sur la source, avec son incertitude de lecture, avant le calcul.

Critères d'acceptation. Travail extérieur moins énergie élastique et cinétique résiduelles égal à
0,62 J à 5 % près ; énergie cohésive dissipée (`eGc`) à 5 % de la même valeur ; pic de force à
10 % de la valeur numérisée ; flèche au pic à 15 % ; écart de pic entre h = 5 et 2,5 mm sur le
ligament inférieur à 3 %.

Développement. Aucun dans le code ; un script de maillage Gmsh avec entaille et groupes
d'appuis (une demi-journée).

Coût. Quelques dizaines de milliers de tétraèdres, mais un chargement quasi statique en
explicite : plusieurs centaines de milliers de pas. LONG (poste local).

Réserves. Béton et non granite ; dimensions et paramètres vérifiés seulement par des extraits de
moteur de recherche. Extension proposée : les quatre tailles de Wendner et al. 2015 (D = 40 à
500 mm, Table 4 ; G_f par la méthode du travail de rupture, Table 3) pour un effet d'échelle sur un
même gâchage, données téléchargeables (URL donnée dans le résumé de l'article).

### Banc 6 (B12) : taillant à trois boutons sur granite de Kuru

Référence. T. Saksala, D. Gomon, M. Hokka, V.-T. Kuokkala, Numerical and experimental study of
percussive drilling with a triple-button bit on Kuru granite, IJIE 72 (2014) 56-66,
doi:10.1016/j.ijimpeng.2014.05.006 (PDF local).

Objet. Taillant à trois boutons monté au bout d'une barre incidente (Ø32 mm, 1 190 mm), frappé par
un percuteur (Ø22 mm, 200 mm) à 10, 16 et 22 m/s ; courbes force-pénétration (fig. 9), cratères et
écaillage entre boutons (fig. 5). Rayon de la surface de contact des boutons : 4,5 mm (texte
qui précède la fig. 4).

Ce que le banc teste dans rockim. La percussion à plusieurs inserts, l'interaction des cratères et
l'écaillage latéral, avec un train maillé (corps nommés et `groupVel`, comme
`configs/yang2026_kuru_train1_v5.cfg`). C'est la seule mesure de force en percussion indépendante
du calage de rockim (déjà relevé dans `vv/B_impact_coupe/README.md`).

Matériau. Jeu Kuru de rockim (Yang 2026, Table 1) ; Saksala 2014 donne E = 67 GPa, f_t0 =
11,4 MPa, f_c0 = 386,8 MPa, ν = 0,26, ρ = 2 630 kg/m³, G_I = 100 J/m² (Table 1) pour comparaison.

Critères d'acceptation (proposés). À 10 et 16 m/s : force maximale à 20 % de la fig. 9,
pénétration maximale à 25 % ; écaillage entre boutons absent à 10 m/s et présent à 22 m/s (fig. 5 ;
à 16 m/s il apparaît dans deux essais sur cinq, §4). À 22 m/s, rapport sans verdict (le modèle de
Saksala lui-même s'en écarte).

Développement. Aucun dans le code ; préparation du maillage du taillant (géométrie de la fig. 4b),
un à deux jours. Coût : LONG (train de plus d'un mètre, comparable aux runs Kuru de Yang, plusieurs
heures sur 14 fils).

Réserves. Valeurs à numériser sur figures ; contact initial des trois boutons sensible au
parallélisme (l'article le note, §4) ; G_I de Saksala (100 J/m²) double de celui du jeu rockim.

### Banc 7 (B13) : coalescence de deux fissures en compression

Référence. A. Bobet, H. H. Einstein, Fracture coalescence in rock-type materials under uniaxial
and biaxial compression, IJRMMS 35 (1998) 863-888 (absent de la bibliographie locale).

Objet. Éprouvettes de gypse à deux fissures préexistantes, contraintes d'amorçage et de
coalescence selon l'inclinaison, l'espacement et le recouvrement des fissures.

Ce que le banc teste dans rockim. Amorçage des fissures en aile et secondaires, et leur
coalescence, contre l'essai (P3.3 ne compare qu'à la MLER et au code de Lisjak). Mise en œuvre
directe par le deck de `vv/P3_aile_lisjak` (2D, `preBrokenJoints`, platines).

Référence chiffrée. À extraire de l'article avant toute préparation : la source n'a pu être lue
ni localement ni en ligne, et aucune valeur n'est reprise ici.

Coût. 2D, compression lente : de l'ordre d'une heure par géométrie à 2 fils par analogie avec P3.3.
LONG si plus de quatre géométries.

## 5. Decks prêts et lancement

Trois bancs ont leurs decks dans `configs/nouveaux_bancs/` (chemins relatifs à la racine du
dépôt). Chaque deck a démarré sans erreur ni clé non lue (sauf `seed`, lecture conditionnelle
légitime) lors d'un essai de quelques pas.

### 5.1 Préparation

```bash
mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j && cd ..
pip install gmsh          # sous Linux, gmsh demande aussi libglu1-mesa
# B10 : maillage déjà présent (meshes/nb_ecaillage_h2p5.msh, 3,2 Mo) ; pour le régénérer :
python3 tools/make_unstructured_mesh.py box3dbc 0.0399 0.0399 0.14 0.0025 \
    docs/rapport_guide/simulations_a_lancer/meshes/nb_ecaillage_h2p5.msh 1
# B11 : maillage du disque QS déjà présent (meshes/nb_disc41_h05_f5_s1.msh) ; decks SHPB régénérables :
python3 docs/rapport_guide/simulations_a_lancer/configs/nouveaux_bancs/bresilien_dyn_kuru/gen_mesh.py
python3 docs/rapport_guide/simulations_a_lancer/configs/nouveaux_bancs/bresilien_dyn_kuru/gen_decks.py
# P1.4 : maillage à générer (34 Mo, 1,5 min), non laissé dans le dépôt
python3 docs/rapport_guide/simulations_a_lancer/configs/nouveaux_bancs/lamb_pekeris/gen_mesh.py 1.5 1
```

### 5.2 Commandes

Le lanceur existant prend un sous-dossier comme groupe (un calcul à la fois, reprise des calculs
finis, sorties dans `sim_out/<groupe>/`) :

```bash
bash docs/rapport_guide/simulations_a_lancer/lancer.sh nouveaux_bancs/ecaillage_bohus 4
bash docs/rapport_guide/simulations_a_lancer/lancer.sh nouveaux_bancs/bresilien_dyn_kuru 4
bash docs/rapport_guide/simulations_a_lancer/lancer.sh nouveaux_bancs/lamb_pekeris 14
```

Le lanceur prend les decks par ordre alphabétique : dans `lamb_pekeris/`, `lamb_fdem3d.cfg` (LONG)
passe avant les deux decks fem3d. En cloud, lancer ces deux-là seuls :

```bash
for c in lamb_fem3d lamb_fem3d_noabs; do
  OMP_NUM_THREADS=4 ./build/rockim docs/rapport_guide/simulations_a_lancer/configs/nouveaux_bancs/lamb_pekeris/$c.cfg \
      sim_out/nouveaux_bancs/lamb_pekeris/$c > sim_out/nouveaux_bancs/lamb_pekeris/$c.log 2>&1
done
```

Un deck seul :

```bash
OMP_NUM_THREADS=4 ./build/rockim \
    docs/rapport_guide/simulations_a_lancer/configs/nouveaux_bancs/ecaillage_bohus/ecaillage_dif.cfg \
    sim_out/nouveaux_bancs/ecaillage_dif
```

### 5.3 Durées estimées et répartition

Estimations tirées des essais de démarrage à 2 fils sur la machine de la session (4 cœurs), puis
divisées par environ 1,7 pour 4 fils ; l'insertion de joints ajoute un surcoût non mesuré.

| Deck | Solveur | Taille | Pas | Durée à 2 fils | Durée à 4 fils | Où |
|---|---|---|---|---|---|---|
| `ecaillage_bohus/ecaillage_dif.cfg` | fdem3d | 62 217 tets | 22 700 | 25 min | 15 min | cloud court |
| `ecaillage_bohus/ecaillage_nodif.cfg` | fdem3d | 62 217 tets | 22 700 | 25 min | 15 min | cloud court |
| `bresilien_dyn_kuru/shpb_v05_dif.cfg`, `_nodif` | fdem 2D | 73 521 tri | 249 000 | 40 min | 25 min | cloud court |
| `bresilien_dyn_kuru/shpb_v10_dif.cfg`, `_nodif` | fdem 2D | 73 521 tri | 236 000 | 38 min | 22 min | cloud court |
| `bresilien_dyn_kuru/shpb_v15_dif.cfg`, `_nodif` | fdem 2D | 73 521 tri | 229 000 | 37 min | 22 min | cloud court |
| `bresilien_dyn_kuru/shpb_v20_dif.cfg`, `_nodif` | fdem 2D | 73 521 tri | 224 000 | 36 min | 21 min | cloud court |
| `bresilien_dyn_kuru/bd_qs.cfg` | fdem 2D | 14 151 tri | ≤ 1,34 M | ≤ 60 min | ≤ 35 min | cloud long ou poste local |
| `lamb_pekeris/lamb_fem3d.cfg` | fem3d | 686 392 tets | 4 090 | 9 min | 5 min | cloud court |
| `lamb_pekeris/lamb_fem3d_noabs.cfg` | fem3d | 686 392 tets | 4 090 | 9 min | 5 min | cloud court |
| `lamb_pekeris/lamb_fdem3d.cfg` | fdem3d | 686 392 tets | 17 900 | plusieurs heures | 1 à 2 h | poste local (14 fils) |

Total cloud court : deux écaillages, huit SHPB et deux Lamb fem3d, environ 4 h 30 à 4 fils en
série. Poste local : le brésilien quasi statique (s'il n'est pas lancé en cloud) et Lamb fdem3d.
Ordre conseillé : `lamb_fem3d` (vérification la moins chère, valide la chaîne des groupes ponctuels),
puis `bd_qs` (dénominateur de B11), puis les SHPB `dif`, puis l'écaillage.

## 6. Réserves générales

- Les impulsions et les courbes de référence lues sur figures l'ont été à l'œil sur le rendu des
  PDF ; une numérisation outillée (même procédé que `vv/P3_fissure_guo`) doit les remplacer avant
  le dépouillement, avec les incertitudes affichées ici comme bornes.
- B10 et B11 font porter le verdict sur le DIF de Yang ; le choix de l'exposant (`yang-fig2`) et la
  provenance de la figure 2b de Yang 2025 conditionnent leur indépendance (fiche 2).
- Les critères de P1.4 sur le plateau statique attendent la transcription de Pekeris 1955.
- Les bancs 4 à 7 n'ont pas de deck : le 4 attend la clé de vitesse radiale, les 5 et 6 un script
  de maillage, le 7 l'extraction des valeurs d'essai.
