# P3 aile : amorçage de fissures en aile (Lisjak 2013)

## Objectif

Vérifier que rockim en mode `fdem` (FDEM 2D de Munjiza, déformation plane) reproduit l'amorçage de fissures en aile depuis une fissure préexistante inclinée sous compression uniaxiale. Deux références :

- B2 (comparaison de codes) : le résultat Y-Geo de Lisjak, σ1c = 7,0 MPa et angle de branchement 64° ;
- P3.3 (validation physique) : la fourchette des modèles de mécanique linéaire élastique de la rupture (MLER) de la Table 3.3, de 5,1 à 14,6 MPa, et l'angle de contrainte circonférentielle maximale, 70,5°.

Source : A. Lisjak, thèse de doctorat, Université de Toronto, 2013, chapitre 3, §3.4.3, pages imprimées 36 à 39 (page PDF = page imprimée + 20). Toutes les valeurs extraites sont dans `lisjak_reference.json`.

## Montage de référence

- Éprouvette homogène 50 × 100 mm, fissure rectiligne de 2c = 5 mm inclinée de γ = 45° sur σ1, au centre.
- Maillage de Delaunay, arête moyenne h = 0,70 mm, 21 600 triangles ; pas de temps 5·10⁻⁹ s.
- Deux platines rigides se déplaçant en sens opposés à 0,05 m/s ; frottement platine-échantillon 0,1.
- Table 3.1 : ρ = 2 300 kg/m³, E = 3 GPa, ν = 0,29, φi = 35°, c = 15 MPa, φf = 35°, ft = 3 MPa, GIc = 2,0 J/m², GIIc = 10 J/m², μ = 7,4·10³ kg/(m·s), pn = 30 GPa·m, pt = 3 GPa/m, pf = 15 GPa.
- Amortissement : μ = 2 μc avec μc = 2h√(ρE) (éq. 3.21) ; recalculé μc = 3 677 kg/(m·s), publié 3,7·10³.
- Résultat Y-Geo : ailes aux deux pointes (événements EA 1 et 2, t = 1,86 à 1,90 ms), σ1c = 7,0 MPa lu sur les réactions des platines, soit 21 % de la résistance ultime (33 MPa dans le texte, 32 MPa dans la légende de la fig. 3.4c) ; angle de branchement 64°.

## Modèles MLER

Notations (compression positive, σ2 = 0) : contrainte normale et cisaillement sur la fissure

σn = (σ1 + σ2)/2 − (σ1 − σ2)/2 · cos 2γ ; τ = (σ1 − σ2)/2 · sin 2γ ; τ* = τ − μ σn ; F = 2c τ*

contrainte normale distante sur le plan de l'aile : σN = σ1 sin²(θ − γ) + σ2 cos²(θ − γ).

Paramètres de la Table 3.2 : 2c = 5 mm, γ = 45°, μ = 0,7, θ = 64°, l = 0,7 mm, σ2 = 0. Amorçage quand KI = KIc (éq. 3.22), avec KIc = √(E GIc) (éq. 3.23) ; Lisjak publie 0,077 MPa√m, le calcul donne 0,0775 (et 0,0809 avec E' = E/(1 − ν²)).

| modèle | formule de KI | publié (MPa) | recalculé, KIc = 0,077 (MPa) | statut |
|---|---|---|---|---|
| Cotterell et Rice 1980 | ¾ [sin(θ/2) + sin(3θ/2)] τ* √(πc) | 5,1 | 5,07 | reproduit |
| Zaitsev 1985 | π τ* c sin θ / (2√(πl)) | 6,9 | 6,82 | candidat non vérifié |
| Horii et Nemat-Nasser 1986 | F sin θ / √(π(l + l*)) − σN √(πl), l* = 0,27c | 14,6 | 14,57 | reproduit |
| Lehner et Kachanov 1996 | F / √(πl) − σN √(π(l + 0,27c)) | 8,5 | 8,53 | candidat non vérifié |

Lisjak ne donne pas les formules, et les articles originaux de Zaitsev et de Lehner et Kachanov ne sont pas dans le dépôt. Les deux formules marquées « candidat » sont des reconstructions qui retombent sur les valeurs publiées à 1 % près ; elles ne sont pas démontrées et ne servent qu'à l'ordre de grandeur. Le critère P3.3 utilise les valeurs publiées. Avec E' au lieu de E, toutes les valeurs montent de 4,5 %.

Angle de contrainte circonférentielle maximale en mode II pur (Erdogan et Sih) : 3 cos θ = 1, soit θ = 70,53°.

Le calcul est dans `lefm_aile.py` (`python3 lefm_aile.py` régénère `lisjak_reference.json`).

## Mise en œuvre rockim

- Maillage : gmsh, rectangle 50 × 100 mm avec la fissure en ligne incluse (`embed`), algorithme de Delaunay (`Mesh.Algorithm = 5`), MSH 2.2, lu par `mesh = file`. La taille gmsh vaut 1,09 h pour que l'arête moyenne soit h : 23 670 triangles à h = 0,7 mm (Lisjak : 21 600).
- Fissure préexistante : `preBrokenJoints` (§5.11) sur le segment, `preBrokenTol = 1e-6`, `preBrokenAngleDeg = 5` ; seules les arêtes de la ligne incluse naissent rompues (D = 1). Frottement de la fissure `jointResidualMu = 0.70` (tan 35°). `gcSurfaceRefresh = eager`, recommandé par le solveur avec cette clé.
- Chargement : `scenario = tension`, `loading = platens`, `pullV = -0.1` (fermeture totale, soit deux platines à 0,05 m/s), `pullRamp = 1e-4`. La vitesse de fermeture de 0,1 m/s est cohérente avec Lisjak : E' ε = 6,1 MPa à 1,87 ms.
- Matériau : `rho = 2300`, `E = 3e9`, `nu = 0.29`, `ft = 3e6`, `cohesion = 15e6`, `frictionDeg = 35`, `Gf = 2.0`, `gfShearFactor = 5` ; adoucissement `jointSoftening = yan` (courbe de Munjiza, éq. 3.4, constantes 0,63, 1,8 et 6,0).
- Frottement platine : `contactMu = 0.1`.
- Amortissement : `bulkViscosity = 3.7e3`, soit µ/2 : Lisjak pose C = µI avec µ = 7,4e3 kg/(m·s) (convention de Munjiza, contrainte visqueuse µD) et rockim écrit 2µD, comme établi pour Solidity au §5.4 quinquies de la documentation ; `dampingLocal = 0`, `meanTensionCapFactor = 0` (pas de garde-fou absent de Y-Geo), `jointXi = 0.01` (règle maison en quasi-statique).
- Pénalité des joints : `jointPenaltyFactor = 20` (défaut rockim) ; variantes à 2,5 et 5, qui encadrent l'équivalent de pf = 15 GPa = 5E de Y-Geo selon que sa raideur vaut pf/(2h) ou pf/h (voir les doutes).
- Sorties utilisées : `history.csv` (colonne `sigma` = réaction de la platine supérieure / W), `fdem_final_joints.csv` (`tBreak`, `breakMode`), `fdem_seismic.csv` (`tYield`, avec `microseismic = true`).

## Métriques

- σ1c : contrainte axiale à l'instant `tBreak` du premier joint non préexistant rompu dont le milieu est à moins de 1,5 h d'une pointe. Indicateur secondaire : contrainte à l'entrée en endommagement (`tYield`) du premier joint proche d'une pointe, borne basse.
- Angle de branchement θ : angle signé (sens trigonométrique) entre la tangente sortante de la fissure et le vecteur pointe → extrémité lointaine du premier joint rompu. Pour cette fissure (bas-gauche vers haut-droite), le côté en traction est le sens trigonométrique aux deux pointes ; une valeur négative signale le mauvais côté.
- Trajet : amas connexe de joints rompus issu de chaque pointe en fin de calcul ; longueur de corde, angle de la corde par rapport à σ1 et angle signé.

## Variantes

| nom | h (mm) | fermeture (m/s) | T (ms) | rôle |
|---|---|---|---|---|
| h0p7 | 0,7 | 0,1 | 3,5 | référence Lisjak |
| h1p0 | 1,0 | 0,1 | 3,5 | maillage plus grossier |
| h0p5 | 0,5 | 0,1 | 3,5 | maillage plus fin |
| h1p0_v0p05 | 1,0 | 0,05 | 7,0 | vitesse divisée par deux |
| h1p0_pf2p5 | 1,0 | 0,1 | 3,5 | `jointPenaltyFactor = 2.5` |
| smoke_h2_v0p5 | 2,0 | 0,5 | 1,2 | essai de chaîne seulement |

T = 3,5 ms couvre environ 1,8 fois l'amorçage Y-Geo et reste loin de la rupture (9 ms).

## Critères (fixés avant tout calcul, le 2026-10-04)

Sur la variante h0p7 :

- P3.3 : 5,1 MPa ≤ σ1c ≤ 14,6 MPa ;
- B2 : |σ1c / 7,0 − 1| ≤ 15 %, soit σ1c entre 5,95 et 8,05 MPa ;
- angle : θ à 10° près de 64° ou de 70,5°, aux deux pointes, du bon côté (θ > 0) ;
- les deux pointes amorcent.

Sensibilité : |σ1c(h) / σ1c(0,7) − 1| ≤ 15 % pour h = 1,0 et 0,5 mm ; |σ1c(v/2) / σ1c(v) − 1| ≤ 10 % à h = 1,0 mm (quasi-statisme).

## Lancer

Depuis la racine rockim :

```
python3 vv/P3_aile_lisjak/p3_aile_lisjak.py prepare            # maillages et decks
python3 vv/run_queue.py P3_aile_lisjak --slots 4 --threads 2   # campagne
python3 vv/P3_aile_lisjak/p3_aile_lisjak.py analyse            # resultats.json, figures
python3 vv/P3_aile_lisjak/p3_aile_lisjak.py all --smoke --threads 1   # essai de chaîne
```

Sorties : `out/<variante>/`, `resultats.json`, `fig_sigma_amorcage.pdf`, `fig_ailes.pdf` (et leurs PNG).

## Essai de chaîne (2026-10-04, un fil)

`smoke_h2_v0p5` : 3 180 triangles, 3 joints préexistants sélectionnés, dt = 8,5·10⁻⁹ s, 141 418 pas, 33 s. Ailes aux deux pointes au même instant (1,14 ms), σ1c = 16,7 MPa, entrée en endommagement à 13,8 MPa, θ = +34° du bon côté, cordes de 2,7 mm alignées sur σ1 à 1° près. Bilan d'énergie au résidu 10⁻¹³ %, platines équilibrées à 2·10⁻⁴ %. Ce maillage (3 arêtes sur la fissure) et cette vitesse (5 fois Lisjak) ne valent rien pour les critères : ils montrent seulement que la chaîne fonctionne.

Le pas de temps est commandé par la borne diffusive de la viscosité ρh²/(4μ), pas par les ressorts : il varie comme h². Coût estimé de la campagne avec 2 fils par run : environ 1 h pour h0p7, 3 à 4 h pour h0p5, 15 à 30 min pour chaque variante à h = 1 mm.

## Doutes et limites

- Convention de la viscosité : rockim applique 2μD. La documentation (§5.4 bis) assimile le μ de Munjiza au μ de rockim (ξ = 2 vaut le critique 2h√(Eρ)). Si Y-Geo applique μD, la valeur équivalente serait 3,7·10³ et le pas doublerait.
- Pénalité des joints : dans Y-Geo, pf = 15 GPa = 5E, et l'ouverture au pic op = 2h ft / pf donnerait une raideur pf/(2h), soit un facteur rockim de 2,5. Cette lecture n'est pas vérifiée, d'où la variante `h1p0_pf2p5`.
- `contactMu` règle à la fois le frottement des platines (0,1) et le contact général des fissures neuves rompues ; Lisjak prend 35° pour celles-ci. Effet attendu faible avant l'amorçage, réel sur le trajet tardif.
- Les pénalités de contact pn et pt de Lisjak ne sont pas transposées (défauts rockim).
- La détection par `tBreak` (rupture complète, D = 1) correspond aux événements EA de Lisjak ; son σ1c est lu au temps de pic d'énergie cinétique, entre la plastification et la rupture de l'élément. Écart attendu de l'ordre de la largeur d'un événement (environ 0,01 ms, soit 0,06 MPa).
- Le maillage est plus dense que celui de Lisjak (23 670 contre 21 600 triangles) parce que l'arête moyenne est calée sur 0,70 mm.
- Formules de Zaitsev et de Lehner et Kachanov non vérifiées (voir plus haut).
