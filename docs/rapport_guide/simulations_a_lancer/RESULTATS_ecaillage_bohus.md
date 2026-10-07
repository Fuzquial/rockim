# Résultats — banc B10, écaillage d'un barreau de granite de Bohus (Saadati et al. 2016)

Calculs du 2026-10-06 (UTC) : `ecaillage_dif` de 23 h 20 à 23 h 35, `ecaillage_nodif` de 23 h 42 à
23 h 57 (relancé après un redémarrage du conteneur qui avait interrompu le premier essai au bout de
3 min). Session cloud Linux, binaire `f638b0f` (code = `e9ba9a6`), **4 fils OpenMP**, 14 min chacun.
Sorties brutes dans `sim_out/nouveaux_bancs/ecaillage_bohus/` (hors git).

Maillage : `nb_ecaillage_h2p5.msh` régénéré ici par la commande de NOUVEAUX_BANCS.md §5.1 avant que
la session de rédaction ne le versionne (62 217 tétraèdres, même nombre ; h inscrit
0,497/1,216/1,853 mm).

Extraction (script de la session) : vitesse = différence centrée de `U_rear_z` (sommet à 0,18 mm du
centre de la face arrière) et de `U_top_z` (moyenne sur la face arrière, 377 sommets) ; premier
maximum, minimum dans les 40 µs suivantes, rebond dans les 30 µs suivantes ;
σ_Novikov = ½ ρ C₀ ΔV_pb avec ρ = 2 660 kg/m³ et C₀ = 4 050 m/s ; plan d'écaillage = position des
joints rompus (`breakMode` > 0) dans la dernière trame de joints, distance à la face arrière (z = 0,14 m).

## Verdict (variante `dif`, critères de NOUVEAUX_BANCS.md, fiche B10)

| critère | rockim (point arrière) | rockim (moyenne de face) | essai | cible | verdict |
|---|---|---|---|---|---|
| 1. pic élastique de vitesse arrière | 7,579 m/s (−2,2 %) | 7,502 m/s (−3,2 %) | 7,75 m/s | ± 5 % | **PASSE** |
| 2. montée 10 % → pic | 25,2 µs | 24,4 µs | 22 µs | 22 ± 3 µs | **ÉCHOUE de 0,2 µs** au point, passe en moyenne de face |
| 3. ΔV_pb, σ_Novikov | 4,285 m/s → 23,08 MPa | 3,769 m/s → 20,30 MPa | 3,5 m/s → 18,9 MPa | 15,1-22,7 MPa | **ÉCHOUE de 0,4 MPa** au point (+22 %), passe en moyenne de face (+7 %) |
| 4. plan d'écaillage | médiane 55,2 mm (quartiles 50,0-62,8) | — | 57 mm | 57 ± 10 mm | **PASSE** |
| 5. résidu B4 | −5,9e-14 J (8,7e-13 %) | — | — | < 0,1 % | **PASSE** |

Trois critères passent nettement. Les deux autres sont à la limite et dépendent du point de lecture :
au sommet arrière seul, ils sortent de la bande de quelques pour cent ; en moyenne sur la face, ils y
rentrent. L'essai mesure par vélocimètre laser en un point, ce qui plaide pour la lecture ponctuelle ;
le sommet est à 0,18 mm du centre et un seul sommet porte un bruit de maillage. **À trancher par la
session de rédaction**, avec une lecture moyennée sur un petit disque (Ø 2 à 5 mm, taille du spot).

## Les deux variantes

| grandeur | `dif` | `nodif` | essai |
|---|---|---|---|
| pic de vitesse arrière (point) | 7,579 m/s à 68,6 µs | 7,579 m/s à 68,6 µs | 7,75 m/s |
| minimum suivant | 3,295 m/s à 103,1 µs | 4,959 m/s à 102,3 µs | 4,25 m/s vers 97 µs |
| rebond | 4,317 m/s à 114,6 µs | 5,338 m/s à 116,4 µs | 4,95 m/s vers 107 µs |
| ΔV_pb (point / face) | 4,285 / 3,769 m/s | 2,620 / 2,385 m/s | 3,5 m/s |
| σ_Novikov (point / face) | 23,08 / 20,30 MPa | 14,11 / 12,85 MPa | 18,9 MPa |
| plan d'écaillage, médiane (quartiles) | 55,2 mm (50,0-62,8) | 52,7 mm (46,4-71,2) | 57 mm |
| joints insérés / rompus | 15 854 / 1 594 | 18 989 / 1 819 | — |
| mode de rupture | 1 575 traction, 19 cisaillement | 1 786 traction, 33 cisaillement | — |
| fragments, volume détaché | 2 ; 3,8e-8 m³ | **5 ; 8,13e-5 m³** | écaille détachée |
| DIF de traction à l'insertion | médiane 1,833 (1,53-1,85), ε̇ médian 89 s⁻¹ | — | ε̇ ≈ 70 s⁻¹ |
| résidu B4 | 8,7e-13 % | 5,7e-13 % | — |

Lecture :

- Sans DIF, la résistance apparente vaut 14,1 MPa pour f_t = 8 MPa : l'inertie et la loi cohésive
  portent seules un facteur 1,76. Avec DIF (médiane 1,83, 35 % des insertions au plafond de 1,85),
  23,1 MPa. L'essai (18,9 MPa) tombe entre les deux, plus près de la variante `dif`.
- Le DIF de Yang est ici appliqué à un ε̇ médian de 89 s⁻¹, au-dessus des 70 s⁻¹ de l'article, et
  pour plus d'un tiers des facettes au plateau (> 100 s⁻¹) : le facteur ne discrimine plus.
- **Sans DIF, l'écaille se détache** : 8,13e-5 m³, soit un bloc de 39,9 × 39,9 × 51 mm (cohérent avec
  la médiane de 52,7 mm). Avec DIF, rien ne se détache dans les 150 µs (3,8e-8 m³), la fissure est
  amorcée mais pas traversante. L'essai montre une écaille. À vérifier sur les trames : faut-il un
  T plus long pour la variante `dif` ?
- L'épaisseur tirée de la période de rebond (C₀Δt/2, 87 à 135 mm selon la lecture) ne concorde pas
  avec la position des joints rompus (≈ 55 mm) : le rebond simulé n'est pas une réverbération propre
  dans l'écaille, notamment avec DIF où celle-ci ne se détache pas. La position des joints rompus
  est la mesure retenue pour le critère 4.
- 50,6 % des facettes insérées en `dif` restent dormantes (jamais ouvertes au-delà de 0,05 dnF) ;
  le solveur signale un seuil d'insertion ou un `insertionHoldSteps` trop bas.

## Complément du 2026-10-07 : lecture au spot et T doublé

Deux calculs supplémentaires, avec un groupe `box.spot` (et `force.spot = 0 0 0`, pour que
`U_spot_*` soit écrit dans `history.csv`) centré sur la face arrière. **Le spot de Ø 3 mm demandé
n'est pas réalisable sur ce maillage** : `box.` retient les faces extérieures dont les trois
sommets sont dans la boîte, et une boîte de 3 × 3 mm ne contient aucune face à h = 2,5 mm (le
solveur refuse le deck : « aucune face exterieure dans la boite »). Le plus petit spot qui marche
est une boîte de **6 × 6 mm** (aire d'un disque de Ø 6,8 mm). Un spot de Ø 3 mm demanderait un
raffinement local de la face arrière.

Le conteneur était environ 2,3 fois plus lent par pas après deux redémarrages (même calcul, même
nombre de pas : 690 s contre 311 s) : les deux calculs ont été coupés à 30 min.

| variante | lecture | pic (m/s) | montée 10 % → pic (µs) | ΔV_pb (m/s) | σ_Novikov (MPa) |
|---|---|---|---|---|---|
| `dif` (T = 300 µs, coupé à 129,8 µs) | point arrière | 7,592 | 25,7 | 4,283 | 23,07 |
| | **spot 6 × 6 mm** | **7,593** (−2,0 %) | **26,2** | **4,260** | **22,94** |
| | moyenne de face | 7,568 | 25,9 | 3,822 | 20,59 |
| `nodif` (T = 150 µs, coupé à 144,6 µs) | point arrière | 7,579 | 25,2 | 2,620 | 14,11 |
| | **spot 6 × 6 mm** | **7,516** (−3,0 %) | **24,1** | **2,537** | **13,66** |
| | moyenne de face | 7,502 | 24,4 | 2,385 | 12,85 |

Critères de la fiche B10 (variante `dif`) relus au spot, sans les changer :

| critère | spot 6 × 6 mm | cible | verdict au spot |
|---|---|---|---|
| 1. pic | 7,593 m/s (−2,0 %) | 7,75 ± 5 % | passe |
| 2. montée | 26,2 µs | 22 ± 3 µs | échoue (+1,2 µs au-delà de la bande) |
| 3. σ_Novikov | 22,94 MPa | 15,1-22,7 MPa | échoue de 0,24 MPa (+21 % sur 18,9) |

Le spot se comporte comme le point arrière (écart < 1 %), pas comme la moyenne de face : la
moyenne de face, qui passait, adoucit le signal par les bords du barreau. La lecture fidèle à un
vélocimètre laser donne donc **trois critères sur cinq**, les deux autres hors bande de peu. Le
calcul `dif` à T doublé reproduit le premier à mieux que 0,1 % sur la partie commune (23,07 contre
23,08 MPa au point), ce qui confirme la répétabilité à 4 fils.

**Détachement de l'écaille : non tranché.** À 129,8 µs (`dif`, T doublé coupé), 2 fragments et
3,8e-8 m³ détachés, comme à 150 µs dans le premier calcul. Mais `nodif`, qui se détachait à 150 µs
(8,1e-5 m³), ne s'est pas encore détaché à 144,6 µs (5,7e-9 m³) : le détachement survient tard,
entre 145 et 150 µs. Les 130 µs atteints en `dif` ne disent donc rien. À relancer sur le poste
local avec T = 300 µs (environ 30 min à 4 fils sur une machine normale).

## Non fait

Déformations aux jauges g1, g2, g3 contre la fig. 2b (les déplacements `U_g*` sont dans
`history.csv`) ; sensibilité à G_f (35 et 140 J/m²) ; lecture moyennée sur un spot.
