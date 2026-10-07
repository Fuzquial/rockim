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

## Non fait

Déformations aux jauges g1, g2, g3 contre la fig. 2b (les déplacements `U_g*` sont dans
`history.csv`) ; sensibilité à G_f (35 et 140 J/m²) ; lecture moyennée sur un spot.
