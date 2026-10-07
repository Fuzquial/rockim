# Résultats — groupe `yan` (rejeu de Yan et al. 2023 avec le binaire courant)

Calculs du 2026-10-06, 22 h 05 à 22 h 19 (UTC), session cloud Linux (4 cœurs).
Binaire : `f638b0f` (code identique à `f0209ef`/`e9ba9a6`, dernier changement de `src/` le 03/10),
compilé en Release (CMake, g++, Eigen 3.4.0). Lancement :
`bash docs/rapport_guide/simulations_a_lancer/lancer.sh yan 4`, soit **4 fils OpenMP** : les sorties
ne sont donc pas garanties identiques au bit au calcul série (`lancer.sh` : 1 fil = bit-identique).
Sorties brutes dans `sim_out/yan/` (hors git) ; chiffres extraits des `.log` et des `history.csv`.

## Synthèse

| essai | grandeur | rejeu 06/10 | référence | écart |
|---|---|---|---|---|
| bande pesante, adaptatif | tassement du sommet | −6,00866e-4 m | −6,0087e-4 m (éq. 20) | −0,0007 % (Yan : 0,23 %) |
| bande pesante, intrinsèque 100E | tassement du sommet | −6,18412e-4 m | −6,0087e-4 m | +2,92 % |
| bande pesante, intrinsèque 1E | tassement du sommet | −2,28024e-3 m | −6,0087e-4 m | +279,5 % (Yan : 208 %) |
| brésilien adaptatif | σt = 2P/(πDt) | 5,3031 MPa | août : 5,30 MPa | +0,06 % |
| compression uniaxiale adaptative | UCS | 51,114 MPa | août : 51,07 ; macOS 04/10 : 51,43 ; Yan ≈ 51 | +0,09 % / −0,61 % |
| compression uniaxiale adaptative | E apparent (plateaux) | 15,54 GPa | août : 15,5 GPa | +0,3 % |
| compression uniaxiale adaptative | joints rompus / insérés | 333 / 1 351 | août : 327 / 1 288 | +1,8 % / +4,9 % |
| triaxial adaptatif | φ (ajustement 0-40 MPa) | 22,92° | août : 22,8° ; Yan : 22,87° | +0,12° / +0,05° |
| triaxial adaptatif | c (même ajustement) | 16,74 MPa | entrée : 16,4 MPa | +2,1 % |

Lecture : **la reproduction d'août tient avec le binaire courant.** Sur ce conteneur, l'UCS
retombe à 51,11 MPa (+0,09 % sur août, dans la tolérance de 0,29 %), et les joints rompus à 333
(+1,8 %). Les +34 % de joints rompus et les 51,43 MPa du rejeu macOS du 04/10 ne viennent donc pas
de l'évolution du code. L'hypothèse restante est celle du chapitre : le tirage de Voronoï, qui
dépend de la bibliothèque standard (libstdc++ ici, libc++ sous macOS ; la chaîne Windows d'août
n'a pas été vérifiée).

## Bande pesante (`strip_adap`, `strip_i1`, `strip_i100`)

Déplacement vertical moyen des 32 nœuds du sommet (ligne `[FDEM] gravity … mean top displacement`).
Durées : 4,3 s, 6,1 s et 23,4 s.

| deck | uy sommet (m) | écart à l'éq. 20 |
|---|---|---|
| `strip_adap` | −6,00866e-4 | −0,0007 % |
| `strip_i100` | −6,18412e-4 | +2,92 % |
| `strip_i1` | −2,28024e-3 | +279,5 % |

L'adaptatif retrouve la solution exacte. L'intrinsèque à 1E est plus souple que chez Yan (279 % contre
208 %), un écart de même sens et de même ordre ; à vérifier si la définition de la pénalité (E/h,
épaisseur de référence) diffère de l'article.

## Brésilien adaptatif (`bd_adap`)

Durée : 1 min 51 s. D = 50,8 mm, t = 1 m.

- pic P = 438 793 N à t = 0,9815 ms, aucun joint rompu au pic (pic verrouillé à la chute) ;
- σt = 2P/(πDt) = **5,3031 MPa**, soit 4,08 × ft (1,3 MPa) ; août : 5,30 MPa ;
- faciès : 88 des 140 joints rompus à moins de 0,15 R de l'axe de chargement (62,9 % diamétral),
  |x−xc|/R moyen 0,232 ; 8 fragments ; 117 intra + 23 homo rompus.
- jauge élastique sur la bande σt ∈ [0,3 ; 0,8] ft : σxx centre / σt = 0,992 (PASS) ;
  au dernier pas avant la première rupture : 0,700 (FAIL, bande 0,85-1,25), avec 4,7 % des joints
  au-dessus de D = 0,01 (rochet sous-critique de la pénalité intrinsèque, signalé par le solveur).
  Même comportement attendu qu'en août, à documenter dans r05a si ce n'est pas déjà fait.

## Compression uniaxiale adaptative (`ucs_adap`)

Durée : 43 s.

- UCS = **51,114 MPa** (39,3 × ft) ; le maximum glissant et le pic de la courbe coïncident ici ;
- E apparent (régression σ-ε entre 20 et 60 % du pic, avant toute rupture, 126 points) :
  plateaux 15,54 GPa, faces 15,94 GPa, extensomètre 15,89 GPa ;
- 1 351 / 3 525 joints insérés (38,3 %), 333 rompus (165 traction, 168 cisaillement, 50,5 % en
  cisaillement) ; 297 intra + 36 homo.

## Triaxial adaptatif (`tx_adap_0/10/20/40`)

Durées : 53 s, 1 min 51 s, 3 min 00 s, 4 min 01 s. Le pic est pris **sur la courbe** : maximum de
`sigma` avant la première chute post-rupture (σ < 0,85 × max courant avec des joints rompus), et non
le maximum glissant du solveur, qui attrape la sonnerie après rupture.

| σ3 (MPa) | σ1 au pic, courbe (MPa) | ε plateaux au pic | max glissant du solveur (MPa) | août |
|---|---|---|---|---|
| 0 | 51,06 | 0,00325 | 51,10 | 51,1 |
| 10 | 73,57 | 0,00451 | **131,83** (sonnerie) | 73 (112 enregistrés) |
| 20 | 94,51 | 0,00568 | 101,77 | 95 ou 100,5 selon la planche |
| 40 | 142,25 | 0,00760 | 157,53 | — |

Ajustement linéaire σ1 = 50,51 + 2,276 σ3 → **φ = 22,92°, c = 16,74 MPa** (résidus +0,54, +0,30,
−1,53, +0,69 MPa). Entrée : φ = 23°, c = 16,4 MPa ; août : φ = 22,8° ; Yan : 22,87°.
À 20 MPa, le pic de courbe (94,5 MPa) tranche entre les deux planches d'août en faveur de 95 MPa.
Confinement à 10 MPa restitué à −0,045 % en fin de rampe.

Le maximum glissant reste faux à 10 MPa (131,8 MPa contre 73,6 réels, +79 %) : le correctif
« pic verrouillé à la chute » ne suffit pas en triaxial. À signaler, puisque les chiffres du `.log`
ne sont pas utilisables tels quels pour les essais confinés.

## Non fait

- Compression à maille 0,75 mm (maillage de l'article) : poste local (~1 h, STRATEGIE.md).
- Profil diamétral du brésilien avec la bonne référence σyy = −(2P/πDt)[4/(1−r²) − 1] : non retracé
  dans cette session (les champs sont dans `sim_out/yan/bd_adap/`).
