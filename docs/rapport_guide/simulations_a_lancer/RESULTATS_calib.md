# Résultats — groupe `calib` (jeu Red Bohus de juillet avec le binaire courant)

Calculs du 2026-10-06, 22 h 36 à 22 h 43 (UTC), session cloud Linux, binaire `f638b0f`
(code = `e9ba9a6`), **4 fils OpenMP**. Sorties brutes dans `sim_out/calib/` (hors git).
Les quatre decks reprennent le jeu de `configs/bohus_gbm_calibrated.cfg` (vérifié : mêmes clés,
seuls `pullV`, `pullRamp`, `T` et le confinement changent, conformément à la recette de son en-tête).

## Synthèse

| essai | grandeur | rejeu 06/10 | archive de juillet | essai (Dumoulin / BTS) | écart à juillet | écart à l'essai |
|---|---|---|---|---|---|---|
| `R_tens` | σt (traction directe) | 16,67 MPa | 17,07 MPa (cible 18,3) | 10,27 MPa (BTS) | −2,3 % | +62 % |
| `R_ucs` | UCS | **33,03 MPa** | 139,1 / 141,5 MPa | 126,6 MPa | **−76 %** | −74 % |
| `R_tx20` | σ1 au pic, brut | 69,22 MPa | 186 MPa | 424,8 MPa | −63 % | −84 % |
| `R_tx50` | σ1 au pic, brut | 122,01 MPa | 247 MPa | 649,0 MPa | −51 % | −81 % |

**La traction tient, la compression s'effondre.** Avec le binaire courant, le jeu de juillet ne donne
plus qu'un quart de son UCS. Le chiffre n'est donc pas « à situer » : il n'est plus reproductible.

## Le décrochage de la compression précède l'historique git

Le même `R_ucs.cfg`, à 1 fil, avec deux binaires plus anciens compilés pour l'occasion :

| binaire | date | UCS |
|---|---|---|
| `f638b0f` (courant, 4 fils) | 06/10 | 33,03 MPa |
| `139df12` (`origin/main`) | 25/08 | 33,04 MPa |
| `5f8f71b` (premier commit de l'historique) | 19/08 | 33,04 MPa |

Le résultat est identique depuis le **premier commit sous git (19/08)**. Le changement qui fait passer
l'UCS de 139 à 33 MPa est donc antérieur à l'historique : c'est l'un des correctifs ou changements de
défauts d'août importés avec l'arbre, entre la calibration du 31/07 et le 19/08. On ne peut pas le
bissecter depuis ce dépôt ; il faudrait les binaires `rockim_f1`/`rockim_f2` du poste local.

Forme de la courbe `R_ucs` (ce n'est pas une rupture fragile) : montée linéaire jusqu'à
31,7 MPa à t = 0,35 ms (fin de rampe 0,30 ms), **pic de 33,02 MPa à t = 0,376 ms avec 0 joint
rompu**, plateau à 32-33 MPa jusqu'à 0,54 ms avec 2 à 9 joints rompus, puis décroissance lente vers
17-22 MPa avec 72 joints rompus sur 2 054 en fin de calcul. Le plafond arrive avant toute rupture :
c'est une saturation (endommagement diffus des joints sans rupture, ou plasticité/écrasement du
continu), pas une fissure. Le rapport 0,97 × ft (34 MPa) imprimé par le solveur est suggestif : la
compression plafonne au niveau de la résistance en traction des joints. Piste à vérifier en premier :
un défaut de `jointShearUnload`, de la loi d'adoucissement ou du critère de cisaillement des joints
(c = 13,6 MPa, φ = 13,4°, `gfShearFactor` = 10) changé en août.

## Biais ν/(1−ν)σ3 des mors bloqués

Mesuré sur `history.csv` à la fin de la rampe de confinement (t = 0,10 ms), avant la charge axiale :

| σ3 | σyy lu à 0,10 ms | ν/(1−ν)σ3 (ν = 0,25) | écart |
|---|---|---|---|
| 20 MPa | 6,56 MPa | 6,67 MPa | −1,6 % |
| 50 MPa | 16,44 MPa | 16,67 MPa | −1,4 % |

Le biais est confirmé : les mors bloqués latéralement pendant la mise en confinement ne portent que
ν/(1−ν)σ3 au lieu de σ3. Le déviateur au pic vaut donc 69,22 − 6,56 = 62,66 MPa (20 MPa) et
122,01 − 16,44 = 105,57 MPa (50 MPa), soit σ1 corrigé = σ3 + q = **82,7 MPa** et **155,6 MPa**.
Même corrigés, les pics restent très loin de juillet et de l'essai. Confinement atteint au cœur :
−20,006 MPa (−0,03 %) et −50,004 MPa (−0,009 %).

Avec les pics corrigés et l'UCS (σ3 = 0, 20, 50) : σ1 ≈ 33 + 2,45 σ3 au mieux, soit une pente
proche de juillet (2,39) mais une ordonnée divisée par quatre.

## Détails

| deck | durée | pic | joints rompus | remarque |
|---|---|---|---|---|
| `R_tens` | 46 s | 16,668 MPa (0,49 ft, référence gbAlphaTen = 0,5) | 36 / 2 054 | chute nette à 0,35 ms |
| `R_ucs` | 88 s | 33,026 MPa | 72 / 2 054 | plateau sans rupture, voir plus haut |
| `R_tx20` | 148 s | 69,218 MPa | 50 / 2 054 | 1 joint rompu au pic ; plateau 55-58 MPa ensuite |
| `R_tx50` | 142 s | 122,007 MPa | 33 / 2 054 | 0 joint rompu au pic ; décroissance lente |

Le maximum glissant du solveur coïncide ici avec le pic de la courbe (pas de sonnerie).
