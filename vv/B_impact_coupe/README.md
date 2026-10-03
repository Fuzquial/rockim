# B_impact_coupe : reproduction d'articles d'impact et de coupe (partie III)

Ce dossier porte deux bancs de la partie III de la campagne (`docs/VV_campagne.md`) : B8, impact
d'un bouton sous une impulsion de tige (Saksala 2011), et B9, coupe au cutter PDC 3D (Heilman et al.
2024, avec les essais de LeBaron et al. 2023 en référence secondaire). Les bancs d'insert unique de
Yang et al. 2025 et 2026 (B7), de Guo (B1) et de Lisjak (B2) sont traités ailleurs.

Fichiers : `b_impact_coupe.py` (actions `prepare`, `run`, `analyse`, `all`, fonction `JOBS`),
`reference_publiee.json` (valeurs publiées, page ou figure, incertitude de lecture, drapeau
`lu_sur_figure`), `resultats.json` et `fig_b8_*`, `fig_b9_*` après la campagne, `resultats_smoke.json`
et `fig_*_smoke.*` pour l'essai de fumée, `meshes/`, `out/` et `out_smoke/`.

## Inventaire et sélection (2026-10-04)

Les titres de `bibliographie/` ont été extraits par `pdftotext -l 1`. Les articles retenus comme
candidats sont les suivants.

| Article | Type | Données publiées | Faisable dans rockim | Indépendance |
|---|---|---|---|---|
| Saksala 2011, IJNAMG 35 | impact, code (FEM axisymétrique) | Table I complète, impulsion, géométrie (fig. 5), courbes F(u) (fig. 6i, 8f) | oui : fem3d, `law = saksala2011` (défauts = Table I), `quarterModel`, `toolShape = flat` et `sphere`, `toolPulse*` | code contre code, rien n'est calé sur ces calculs |
| Saksala et al. 2014, IJIE 72 | impact, essai (Kuru, taillant à trois boutons, barre de Hopkinson) | F(t) et contrainte dans la tige à 10, 16 et 22 m/s, cratères | partiel : un outil analytique n'a qu'un bouton ; il faudrait un taillant maillé en fdem3d (corps nommés), coûteux | paramètres de Saksala calés en quasi-statique, pas sur ces essais |
| Saksala 2012, IJNAMG | impact, code 3D | un bouton (1,13 mm, 19 kN), puis multi-boutons avec trou | un bouton oui ; multi-boutons non | la loi y est modifiée (énergies de rupture) : pas la loi portée |
| Aising et al. 2024 (IJRMMS, RMRE) | impact, essai (insert unique) | enfoncement, cratère, COR | oui, déjà monté | données de B7 (Yang), exclu ici |
| Saadati et al. 2014 ; Shariati et al. 2021 | indentation, Bohus | courbe F(h) expérimentale, faciès | oui en principe (indenteur lent à grande masse) | indépendant du calage Bohus de rockim, mais la géométrie de l'essai est dans une référence absente de la bibliographie |
| Jiang et al. 2020 et 2021 | impact d'indenteur, essai et FDEM | volumes et profondeurs de cratère | non évalué en détail : paramètres de joints bilinéaires d'Abaqus | à vérifier |
| Wu et al. 2021 | indentation dynamique confinée | endommagement par tomographie | confinement triaxial vrai absent en fdem3d | à vérifier |
| Heilman et al. 2024, ARMA 24-0238 | coupe, code (HOSS 3D) | Table 1, géométrie, F(course) pour 12 cas (fig. 2, 3), 3 cas confinés (fig. 4) | oui sans confinement : fdem3d, `toolShape = pdc`, deck `configs/cut3d_heilman.cfg` ; le confinement anisotrope (30,29 et 17,85 MPa) est hors de portée | paramètres calés sur UCS et brésilien 2D, pas sur la coupe |
| LeBaron et al. 2023, GRC 47 | coupe, essai (raclage, granitoïde FORGE) | force moyenne et MSE aux passes 0,508, 1,016, 1,524 mm, garde 20 degrés | oui, mêmes runs que Heilman | aucun calage de Heilman ni de rockim sur ces essais ; vitesse 400 fois plus faible |
| Liu et al. 2025, Petroleum 11 | coupe, code (FDEM à grains) | forces 3D et MSE en fonction de la passe | partiel : tranche à grains Voronoï ; paramètres calés sur leurs propres essais | moyenne |
| Labra, Rojek et Oñate 2017 | coupe au disque de tunnelier, DEM/FEM | essai LCM | non : pas de disque roulant | sans objet |
| Rojek 2014 | coupe thermomécanique DEM, pic de drague | forces et usure | non : pic et thermique absents | sans objet |
| Sun et al. 2025 ; Li et al. 2021 ; Dong et Chen 2018 | coupe percussive PDC | fractales, faciès, ROP | non : outil à chargement axial et en torsion combinés absent | sans objet |

Classement retenu.

| Rang | Impact | Coupe |
|---|---|---|
| 1 | Saksala 2011 : jeu de données complet, même loi que rockim, coût faible | Heilman 2024 : jeu complet, deck existant, 12 courbes chiffrées |
| 2 | Saksala et al. 2014 : seule mesure de force en percussion indépendante du calage, mais un taillant maillé est nécessaire | LeBaron 2023 : seule mesure de coupe sur la même roche, mêmes passes et même garde |

Les deux premiers rangs sont mis en œuvre. LeBaron est dépouillé sur les runs de Heilman sans calcul
supplémentaire.

## B8 : impact d'un bouton, Saksala 2011

Référence. T. Saksala, Numerical modelling of bit-rock fracture mechanisms in percussive drilling
with a continuum approach, IJNAMG 35 (2011) 1483-1505, section 4.1.

Cas. Un bouton rigide de diamètre 10 mm, plat (`sk11_cyl`) ou hémisphérique (`sk11_hemi_R10`), est
poussé par une onde de contrainte de 200 MPa dans une tige de 166,7 mm² (1 000 mm² répartis sur six
boutons), montée en 10 µs, palier jusqu'à 100 µs, descente en 10 µs. Le bouton porte une masse de
calcul de 0,1 kg et un amortisseur ρ_b c_b A qui représente la tige longue (éq. 11 de l'article).
Le massif est homogène, sans tirage de Weibull. Modèle rockim : fem3d, quart de bloc de 100 × 100 ×
150 mm (`quarterModel = true`, masse, force incidente et impédance divisées par quatre), maillage
gmsh raffiné autour du point d'impact (0,6 mm jusqu'à 8 mm du point, 6 mm à 60 mm, 51 984
tétraèdres), bords absorbants de Lysmer, contact de Signorini sans frottement, `dampingLocal = 0`.

Paramètres et provenance.

| Grandeur | Valeur | Provenance |
|---|---|---|
| E, ν, ρ, ft, c₀, φ | 60 GPa ; 0,2 ; 2 600 kg/m³ ; 13 MPa ; 37,5 MPa ; 30° | Table I, p. 1492 |
| paramètres de cap, viscosité, endommagement | défauts de `law = saksala2011` | Table I (port vérifié à 8e-14 contre la VUMAT) |
| σ_A, t_rise, t_dur, t_des | 200 MPa ; 10 µs ; 100 µs ; 10 µs | p. 1492 |
| A_rod, ρ_b, m_bit | 166,7 mm² ; 7 800 kg/m³ ; 0,1 kg | Table I |
| c_b | 5 189 m/s | non publiée : acier, E = 210 GPa |
| rayon du bouton hémisphérique | 10 mm | cote de la fig. 5d, ambiguë : `sk11_hemi_R5` teste 5 mm |
| domaine | rayon 100 mm, hauteur 150 mm | lus sur la fig. 5a, à ±10 mm |

Dans rockim, `toolPulseForce` est l'amplitude de l'onde incidente F_inc = A σ_A et la source vaut
F = 2 F_inc − Z v, ce qui est l'éq. 11 de Saksala. La lecture de t_dur comme début de la descente
donne une impulsion de 110 µs, proche de 2L/c = 116 µs pour le piston de 0,3 m de l'article.

Grandeurs comparées. Force maximale F_max (moyenne glissante sur 2 µs, quart multiplié par quatre),
enfoncement maximal u_max, enfoncement résiduel u_res (enfoncement quand la force retombe sous 2 % du
pic). Références : bouton plat 55 kN, 0,35 mm et 0,245 mm ; bouton hémisphérique 32,5 kN, 0,80 mm et
0,73 mm. F_max et u_max sont imprimés dans le texte ; u_res est lu sur les figures 6(i) et 8(f).

Critères, écrits avant le calcul : |écart| ≤ 15 % sur F_max, 20 % sur u_max, 25 % sur u_res, pour
`sk11_cyl` et `sk11_hemi_R10`. Le budget couvre la lecture, le passage de l'axisymétrique CST aux
tétraèdres 3D et la sensibilité au maillage que Saksala déclare lui-même (fig. 13). `sk11_hemi_R5`
est une sensibilité, sans verdict.

Coût. dt = 7,9e-9 s, 27 900 pas, environ 0,05 s par pas à un fil sur la machine chargée : 10 à
15 minutes par run à deux fils, 30 à 45 minutes pour les trois runs.

Commande.

```bash
cd ~/rockim
python3 vv/B_impact_coupe/b_impact_coupe.py all --only sk11_cyl sk11_hemi_R10 sk11_hemi_R5 --threads 2
```

## B9 : coupe au cutter PDC 3D, Heilman 2024 et LeBaron 2023

Référence principale. E. Heilman et al., FDEM modeling of PDC cutter-rock interactions for
geothermal drilling applications, ARMA 24-0238 (2024), code HOSS, figure 3B (garde 20 degrés,
10 m/s). Référence secondaire, expérimentale. A. LeBaron et al., Numerical simulation of rock
cutting for the Utah FORGE geothermal project, GRC Transactions 47 (2023), figures 3 et 4.

Cas. Éprouvette de granite FORGE de 40 × 30 × 20 mm, entaille de départ de 3 mm de long à la
profondeur de passe plus un jeu de 0,284 mm, fond encastré. Cutter rigide de 13 mm de diamètre et
2,5 mm d'épaisseur, garde de 20 degrés (`backRakeDeg = -20`, convention de §5.6 quater : le corps
monte en arrière), à 10 m/s, passes 0,508, 1,016 et 1,524 mm. Modèle rockim : fdem3d, `scenario =
shear`, `toolShape = pdc`, `toolContact = signorini`, insertion adaptative, contact potentiel,
maillage `tools/make_cut3d_mesh.py` à 0,5 mm dans le couloir de coupe (21 000 à 21 700 tétraèdres).
Le deck reprend `configs/cut3d_heilman.cfg` avec deux différences voulues : `dampingLocal = 0`, règle
des bancs de force (le Cundall freine l'impulsion de contact, §5.6 ter, banc T0b), et
`jointSecantRatchet = on`, parce que `jointShearUnload = origin` seul crée de l'énergie (avertissement
du solveur, §5.4).

Paramètres et provenance.

| Grandeur | Valeur | Provenance |
|---|---|---|
| E, ν, ρ | 48,26 GPa ; 0,23 ; 2 610 kg/m³ | Heilman, Table 1 |
| ft, c, φ | 10,62 MPa ; 44,815 MPa ; atan 0,8 = 38,66° | Heilman, Table 1 |
| GfI, GfII/GfI | 152,9 J/m² ; 645,3/152,9 = 4,220 | Heilman, Table 1 |
| volume | élastique (`crushCap = 1e12`) | Heilman, section 2 |
| frottement cutter-roche | 0,8 | Heilman, Table 1 |
| chanfrein | 0 | présent chez Heilman mais non coté |
| maille fine | 0,5 mm | Heilman : 0,25 mm (trois à huit fois plus cher) |

Grandeurs comparées. La course est comptée depuis le premier contact : la face de coupe penche vers
l'avant et touche l'arête de l'entaille avant que l'arête du cutter n'atteigne la face verticale
(0,37 mm plus tôt à la passe de 1,016 mm). La force est la dérivée du travail de l'outil par rapport
à la course sur une fenêtre de 0,05 mm, ce qui écarte le pic d'impulsion du premier toucher. Contre
Heilman : pic F_pic et sa position x_pic, chute de la force après le pic, ordre des pics selon la
passe. Pics lus sur la figure 3B : 14,8 kN à 1,35 mm, 21,5 kN à 1,65 mm et 30,8 kN à 1,80 mm (±0,5 kN,
±0,05 mm). Contre LeBaron : force moyenne sur la course de 0,2 à 2,0 mm et MSE ; essais : 934 ± 245,
1 668 ± 483 et 2 157 ± 560 N (moyenne ± un écart-type, figure 3, lbf convertis).

Critères, écrits avant le calcul. Heilman : |écart| ≤ 25 % sur F_pic, |x_pic − référence| ≤ 0,5 mm,
force retombée sous 20 % du pic avant 3,0 mm de course, ordre F(0,508) < F(1,016) < F(1,524).
LeBaron : force moyenne dans la bande de ±1 écart-type publiée. Ce second critère est noté sans
ajustement s'il échoue : l'essai est mené à 25 mm/s, le calcul à 10 m/s, et les pics de Heilman sont
déjà dix fois plus grands que les moyennes de LeBaron.

Coût. dt = 1,58e-9 s, 221 000 pas pour 3,5e-4 s, de 7 à 25 ms par pas à un fil selon le nombre de
joints insérés : une à deux heures et demie par run à deux fils, trois à sept heures pour les trois
passes. L'option `--etendu` ajoute la garde de 10 degrés (figure 2B), au même coût.

Commande.

```bash
cd ~/rockim
python3 vv/run_queue.py B_impact_coupe --slots 3 --threads 2
python3 vv/B_impact_coupe/b_impact_coupe.py analyse
```

## Essai de fumée (2026-10-04)

`python3 vv/B_impact_coupe/b_impact_coupe.py all --smoke --threads 1`, `OMP_NUM_THREADS = 1`, 56 s de
calcul au total. Maillages grossiers (2,5 mm pour B8, 1,2 mm pour B9) ; ces chiffres valident la
chaîne, ils ne sont pas des résultats.

| Run | Maillage | Durée simulée | Mesure | Temps |
|---|---|---|---|---|
| sk11_cyl | 3 120 tétraèdres | 220 µs | F_max 46,9 kN (−14,7 %), u_max 0,420 mm (+20,0 %), u_res 0,272 mm | 20 s |
| sk11_hemi_R10 | 3 120 tétraèdres | 220 µs | F_max 24,2 kN (−25,5 %), u_max 0,764 mm (−4,5 %), u_res 0,710 mm | 20 s |
| he24_br20_d1p016 | 2 901 tétraèdres | 30 µs (0,17 mm de course) | 6,75 kN à 0,09 mm, 274 joints rompus, résidu B4 1,0e-3 % | 16 s |

La courbe du bouton plat a déjà la forme de la figure 6(i) de Saksala. Celle du bouton hémisphérique
reste à 2 kN jusqu'à 0,22 mm : un seul nœud porte le contact sur ce maillage, défaut que le maillage
de 0,6 mm supprime. La coupe s'engage sans pompe visible (résidu B4 de 1e-3 %).

## Réserves et capacités manquantes

- Confinement de Heilman (30,29 MPa vertical, 17,85 MPa horizontal, figure 4) : fdem3d n'a qu'une
  pression scalaire `confiningPressure` ; le cas confiné n'est pas reproductible sans un tenseur de
  contraintes in situ (§5.6 quater, phase 2).
- Taillant multi-boutons (Saksala 2014, Saksala 2012) : l'outil analytique de fem3d et fdem3d n'a
  qu'une forme ; seul un taillant maillé en fdem3d avec corps nommés le permettrait.
- Indicateurs de pompe du contact outil (`injection outil`, `v nodale max`) : imprimés en 2D
  seulement ; en 3D la pompe se lit sur le résidu B4 et sur les postes de `history.csv`.
- Chanfrein du cutter de Heilman non coté : arête vive ici.
- B8 compare deux codes qui partagent la même loi : il vérifie la chaîne loi, contact et source à
  impédance, pas la physique du granite.
