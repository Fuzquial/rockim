# Ajouts rockim pour la reproduction de Yan, Zheng & Wang (2023)

Base : `rockim_gbm.tar.gz` (version allégée d'origine). Toutes les additions sont
pilotées par une clé de config dont le **défaut reproduit le comportement
d'origine**. Non-régression vérifiée : `verify_fdem_tension` 9,84423 MPa,
`verify_fdem_voronoi_tension` 10,7066 MPa, `verify_fem_bar` c = 4382,34 m/s,
`verify_dem_tension` 9,98908 MPa, `verify_fdem3d_tension` 9,68194 MPa — identiques
au dépôt d'origine.

## Clés ajoutées

| clé | défaut | effet |
|---|---|---|
| `jointSoftening` | `linear` | `yan` = loi cohésive exponentielle f(D), éq. 11 (a=0,63 b=1,8 c=6,0) |
| `jointFrictionScaled` | `false` | multiplie le terme de Coulomb par f(D) (forme littérale de l'éq. 10) — **non testé sur cas complet** |
| `gravity` | `0` | force de volume nodale ρgA/3 |
| `brazilianStopAfterPeak` / `brazilianStopDelay` | `false` | arrêt de l'essai après la chute post-pic |
| `ucsStopAfterPeak` / `ucsStopDelay` | `false` | idem en compression |
| `gaugeLoFrac` / `gaugeHiFrac` | `0.3` / `0.8` | bande de lecture de la jauge élastique du brésilien |
| `pullDelay` | `0` | délai avant fermeture des plateaux (équilibrage de σ3 en triaxial) |
| `writeJointMode` | `false` | ajoute `failMode` aux VTU de joints |
| `shpbPulse` / `shpbPulseV0` / `shpbPulseTau` | — | impulsion imposée (`halfsine` \| `trapezoid`) |
| `shpbNoDisc` | `false` | barre seule, sans disque (vérification de propagation) |
| `absorbFactor` | `1` | facteur du dashpot de Lysmer (`2` = éq. 21 littérale, **réfléchit 33,3 %**) |

Nouveautés non pilotées par clé : `scenario = shpb` + `geometry = shpb` (montage
barre–disque–barre), champs `sigmaXY` et `epsXX` dans les VTU d'éléments, colonne
`breakMode` dans les VTU de joints, colonnes `peakLocked`, `sxxC`, `syyC` dans
`history.csv`.

## Contrôles passés sur CETTE archive

- `gravity = 9.8`, bande pesante : uy en tête = −6,00866e−4 m contre 6,00866e−4 m
  analytique (éq. 20).
- `jointSoftening = yan` : ∫f(D)dD = 0,386307, pic 9,82577 MPa (contre 9,84423 en
  linéaire).
- SHPB barre seule, T = 4e−4 s : pic ε = −9,31277e−4 à t = 289,2 µs contre
  −9,31415e−4 = −V₀/c analytique, soit **+0,015 %** ; arrivée au seuil 1 % à
  176,7 µs pour 179,12 µs théoriques.

## Fusion

Les six patchs (3 de phase 1 + `mode_rupture`, `ucs_yan`, `patch_shpb`) touchent
tous `FdemSolver.cpp`. Les trois derniers ayant été écrits en parallèle sur la même
base, ils entrent en conflit à l'application séquentielle (7 hunks rejetés). La
fusion livrée ici est une **fusion à trois voies** des fichiers complets, avec 9
conflits résolus à la main — dont deux implémentations indépendantes du mode de
rupture, conservées toutes les deux : `failMode` (continu dans [0,1]) et
`breakMode` (1 = traction, 2 = cisaillement).

## Configs

`configs_yan/` — 42 fichiers : familles en pénalité du brésilien (`bd_*`),
sensibilité au maillage (`mesh_*`), compression uniaxiale (`ucs_*`), triaxiale
(`tx_*`), bande pesante (`strip_*`), SHPB (`shpb_*`).
