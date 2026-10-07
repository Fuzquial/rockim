# Bancs d'insertion des joints, rockim 2D : adaptative contre intrinsèque — résultats

Campagne du 2026-10-07. Binaire `/home/user/rockim/build_contact/rockim` (commit b869dae),
`mode = fdem`, `scenario = tension`, `mesh = file`. Critères écrits avant calcul :
[CRITERES.md](CRITERES.md) (cinq addenda datés, chacun écrit avant les calculs qu'il
concerne ; l'écart de procédure du banc 1 y est déclaré). Figure de synthèse :
[fig_synthese_insertion.pdf](fig_synthese_insertion.pdf). Aucun fichier source du dépôt n'a
été modifié ; rien n'a été commité.

Exécution : `OMP_NUM_THREADS=2 timeout 1200 nice -n 10 build_contact/rockim <deck>.cfg <deck>_out`
(script `scripts/run_queue.sh`, un ou deux calculs à la fois à côté du rejeu tunnel). La durée
de chaque calcul est dans `<deck>.dur` (rc = 124 : arrêt à 20 min), les chiffres dans
`<banc>_resultats.csv`. 120 calculs au total (dont la première série `b1`) ; trois ont atteint la limite de 20 min et un a été arrêté à la main (déclarés plus bas). Les
`.vtu` ont été supprimés après dépouillement ; les segments des joints rompus et les courbes
sont conservés dans les `.npz`.

Matériau commun : E = 50 GPa (200 GPa dans `b1s`), f_t = 10 MPa, c = 25 MPa, joints `yan`
(z-curve, éq. 11 de Yan et al. 2023) sauf mention, pénalité p = facteur·E/h_loc
(h_loc = moyenne des diamètres inscrits des deux triangles). Contrainte macroscopique
σ = |F_mors|/(W·t).

## Résumé

1. Loi de joint (banc 1) : adaptatif et intrinsèque donnent la même loi, à 0,7 % près, en
   traction pure, en cisaillement pur et en mode mixte. L'énergie dissipée vaut G_I·L et G_II·L
   à 0,3 % près (adaptatif 1,0000 ; intrinsèque 10 E 1,0027), l'ouverture de rupture vaut
   δ_c à 0,6 % près (linéaire) et 1,7 % près (`yan`). Seul défaut : en mode mixte avec
   l'adoucissement `linear`, le pic dépasse f_t de 6 à 10 % et l'énergie varie de 55 % entre
   10 E et 100 E, dans les deux schémas. C'est donc la loi qui est en cause, pas l'insertion ;
   avec `yan` le défaut disparaît.
2. Barre en traction (banc 2) : l'adaptatif finit sur une seule fissure traversante dès
   h ≤ ℓ_cz/8, mais avant le pic il insère des joints partout. À h = ℓ_cz/16, 393 joints sur
   1 478 (27 %) sont insérés avant le pic, et 81 % des joints insérés sont à plus de 2h de la
   fissure finale. Le pic monte avec le raffinement (1,08 → 1,15 → 1,24 f_t) : il ne converge
   pas, alors que la première insertion a lieu à 1,00 f_t à toutes les tailles. L'énergie
   converge (1,47 → 1,50 G_I·b, 2,5 %). L'intrinsèque converge sur le pic (1,07 → 1,09 f_t,
   1,2 %) ; son énergie n'est pas jugeable à la maille fine (barre non séparée à 20 min).
3. Contact corrigé (banc 2) : aucun effet sur le pic ni sur le motif en adaptatif (énergie
   ±0,2 %). Le travail du contact reste ≤ 0 partout. Le correctif coûte ×1,3 à plus de ×10 de
   temps de calcul.
4. Plaque entaillée (banc 3) : les deux schémas donnent la même trajectoire, droite depuis la
   pointe (même chemin de 19 joints à h = 2,5 mm). La longueur vaut 1,06 à 1,08 fois le
   ligament et la part de propagation est de 75 à 79 % ; le seuil de 80 % est manqué de peu,
   par les deux schémas. Le pic adaptatif dépasse celui de l'intrinsèque de 12 à 14 %.
   L'adaptatif insère dès 0,25 f_t (singularité) et laisse 65 % de ses joints insérés loin de
   la fissure.
5. Sensibilité (banc 4) : à h = 2,5 mm, aucune des cinq variantes ne change le résultat
   (pic ±0,5 %, énergie ±1,4 %, 9 rompus, une fissure). `max` insère 15 % de joints de plus.
   À h = 5 mm, `elliptic`, `volume` et `max` font chuter l'énergie de 32 à 36 % : la barre
   s'y casse en une fissure au lieu de deux. C'est une bifurcation du maillage grossier, pas
   un effet systématique du critère.

## Banc 1 — joint isolé (mini-maillage de 8 CST)

Montage : `scripts/mesh_b1.py` (bande de 10 mm, joint horizontal ou à 45°, six arêtes
parasites presque verticales), ν = 0, mors latéralement libres ; φ = 0, G_I = 70 J/m²,
G_II = 700 J/m², L = 10 mm (mode I) ou 14,14 mm (45°). Courbes lues sur 400 trames
(`scripts/analyse_b1.py b1s`). Durées : 0,25 à 2,6 s par calcul.

La version initiale `b1` (E = 50 GPa, bande de 80 mm) est instable après le pic : la bande
stocke au pic f_t²H/(2E) = 80 J/m², plus que G_I = 70 J/m², et le joint s'ouvre de 0 à 18 µm
en 66 µs. Ses énergies sont justes (W/G_I·L = 0,997 à 1,009 ; W/G_II·L = 0,9995 à 1,0004),
ses ouvertures de rupture ne sont pas mesurables. Les chiffres retenus viennent de la variante
stable `b1s` (E = 200 GPa, addendum 1) :

| trajet | schéma | pénalité / loi | pic | W [J/m] | W/(G·L) | δ de rupture [µm] (loi) | durée [s] |
|---|---|---|---|---|---|---|---|
| I | adaptatif | 100 / linear | 1,0172 f_t | 0,70000 | 1,0000 | 13,91 (14,00) | 1,54 |
| I | adaptatif | 100 / yan | 1,0172 f_t | 0,70000 | 1,0000 | 17,82 (18,12) | 1,55 |
| I | adaptatif | 10 / linear | 1,0172 f_t | 0,69997 | 1,0000 | 13,92 | 0,58 |
| I | adaptatif | 10 / yan | 1,0172 f_t | 0,69997 | 1,0000 | 17,81 | 0,57 |
| I | intrinsèque | 100 / linear | 1,0155 f_t | 0,70018 | 1,0003 | 13,91 | 1,48 |
| I | intrinsèque | 100 / yan | 1,0155 f_t | 0,70018 | 1,0003 | 17,82 | 1,61 |
| I | intrinsèque | 10 / linear | 1,0039 f_t | 0,70192 | 1,0027 | 13,99 | 0,52 |
| I | intrinsèque | 10 / yan | 1,0039 f_t | 0,70192 | 1,0027 | 17,88 | 0,71 |
| II | adaptatif | 100 / linear | 1,0095 c | 9,8998 | 1,0000 | 55,72 (56,00) | 2,58 |
| II | adaptatif | 100 / yan | 1,0095 c | 9,8997 | 1,0000 | 71,91 (72,48) | 2,13 |
| II | adaptatif | 10 / linear | 1,0081 c | 9,8998 | 1,0000 | 55,42 | 0,81 |
| II | adaptatif | 10 / yan | 1,0095 c | 9,8998 | 1,0000 | 71,69 | 0,91 |
| II | intrinsèque | 100 / linear | 1,0093 c | 9,9002 | 1,0001 | 55,67 | 2,43 |
| II | intrinsèque | 100 / yan | 1,0094 c | 9,8999 | 1,0000 | 71,59 | 2,33 |
| II | intrinsèque | 10 / linear | 1,0097 c | 9,8998 | 1,0000 | 55,80 | 0,73 |
| II | intrinsèque | 10 / yan | 1,0096 c | 9,8999 | 1,0000 | 71,90 | 0,71 |
| mixte, mors bloqués | adaptatif | 100 / linear | 1,1019 f_t | 1,5893 | 1,61 G_I·L | — | 2,04 |
| mixte, mors bloqués | adaptatif | 100 / yan | 1,0248 f_t | 1,6235 | 1,64 G_I·L | — | 1,97 |
| mixte, mors bloqués | adaptatif | 10 / linear | 1,0566 f_t | 2,4662 | 2,49 G_I·L | — | 0,77 |
| mixte, mors bloqués | adaptatif | 10 / yan | 1,0247 f_t | 1,6232 | 1,64 G_I·L | — | 0,72 |
| mixte, mors bloqués | intrinsèque | 100 / linear | 1,1009 f_t | 1,5871 | 1,60 G_I·L | — | 2,06 |
| mixte, mors bloqués | intrinsèque | 100 / yan | 1,0249 f_t | 1,6250 | 1,64 G_I·L | — | 2,16 |
| mixte, mors bloqués | intrinsèque | 10 / linear | 1,0590 f_t | 2,4831 | 2,51 G_I·L | — | 0,83 |
| mixte, mors bloqués | intrinsèque | 10 / yan | 1,0248 f_t | 1,6207 | 1,64 G_I·L | — | 0,82 |

Valeurs attendues de δ de rupture : 2G/f pour la loi linéaire, G/(f·0,386307) pour `yan`.
En mode mixte, le pic est σ_n/f_t au joint (τ = σ_n, soit τ/c = 0,41).

Le pic est lu sur la réaction des mors ; celle-ci porte un bruit inertiel de ±2 %
(oscillation de 9,53 à 10,02 MPa entre deux trames avant le pic). Les 1,7 % au-dessus de f_t
en mode I sont ce bruit : l'insertion adaptative se déclenche à σ = f_t (banc 2 : première
insertion à 1,00 f_t à toutes les tailles).

| critère | verdict |
|---|---|
| B1-a pic mode I (1 ± 0,02) | conforme : 1,0039 à 1,0172 |
| B1-b pic mode II (1 ± 0,02) | conforme : 1,0081 à 1,0097 |
| B1-c pic mixte (1 ± 0,03) | conforme avec `yan` (1,0247-1,0249) ; non conforme avec `linear` (1,057-1,102), dans les deux schémas |
| B1-d énergie mode I (1 ± 0,02) | conforme : 1,0000 (adaptatif), 1,0003 à 1,0027 (intrinsèque) |
| B1-e énergie mode II (1 ± 0,02) | conforme : 1,0000 à 1,0001 |
| B1-f énergie mixte | dans [G_I·L ; G_II·L] partout ; adaptatif/intrinsèque ≤ 0,7 % : conforme. 10 E/100 E : conforme avec `yan` (0,02 %), non conforme avec `linear` (+55 %) |
| B1-g ouverture de rupture (δ_c ± 10 %) | conforme sur `b1s` : 13,91 contre 14,00 µm (linear), 17,82 contre 18,12 µm (`yan`), 55,7 contre 56,0 µm (II) ; non mesurable sur `b1` (instable) |
| B1-h pénalité 10/100 E (≤ 2 %) | conforme en modes I et II (écart maximal 1,2 % sur le pic intrinsèque, dans le bruit) |

Diagnostics complémentaires :
- Avec les mors latéralement libres (`b1s_mix`), le trajet mixte n'est pas piloté : le bloc
  supérieur part en ouverture pure sous `yan` (δ_s = 0,3 µm, W = 1,00 G_I·L) et en glissement
  pur sous `linear`. D'où la variante à mors bloqués (addendum 2). Même à mors bloqués, les
  CST se déforment et le trajet n'est pas proportionnel : sous `yan`, l'ouverture précède le
  glissement.
- Avec `linear` en mode mixte, la traction de cisaillement chute presque à zéro dès
  l'insertion : glissement brusque de 13 µm en 9 µs, puis rechargement jusqu'à 5,5 MPa avant la
  rupture. La dissipation dépend alors de la pénalité. Le même comportement apparaît en
  intrinsèque : c'est la branche linéaire de la loi en mode mixte qui en est la cause.

## Banc 2 — barre 20 × 30 mm en traction, maillage non structuré

ℓ_cz = E·G_I/f_t² = 20 mm (G_I = 40 J/m²) ; h = 10 / 5 / 2,5 / 1,25 mm
(ℓ_cz/2, /4, /8, /16 ; maillages gmsh `meshes/bar_h*.msh`, 32 à 1 012 triangles). Pénalité
100 E/h pour les deux schémas. « brut » : `contact = potential` sans correctif ; « corr » :
`contactCandidates = vertex`, `gcBirth = offset`, `potForceExact = true`. La série `b2`
s'arrête à T = 1,25 ms. La série `b2L` (adaptatif seulement, T = 3 ms, addendum 3) est la
référence de l'adaptatif, car les barres adaptatives h = 10 et 5 mm ne sont pas séparées à
1,25 ms. Séparée : σ_fin < 0,05 f_t. « hors fiss. » : part des joints insérés (adaptatif) ou
endommagés D > 0 (intrinsèque) dont le milieu est à plus de 2h de la fissure principale.

| calcul | σ_pic/f_t | σ_fin/f_t | séparée | insérés / D>0 | rompus | N_fiss | étendue x/W | hors fiss. | W/(G_I·b) | propag. | durée [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| b2L_h10_adapt_brut | 1,223 | 0,001 | oui | 8 | 6 | 1 | 1,00 | 0,00 | 1,741 | 0,80 | 3,5 |
| b2L_h10_adapt_corr | 1,223 | 0,001 | oui | 8 | 6 | 1 | 1,00 | 0,00 | 1,741 | 0,80 | 6,5 |
| b2L_h5_adapt_brut | 1,084 | 0,000 | oui | 26 | 10 | 2 | 1,00 | 0,27 | 2,018 | 0,78 | 7,6 |
| b2L_h5_adapt_corr | 1,084 | 0,000 | oui | 26 | 10 | 2 | 1,00 | 0,27 | 2,018 | 0,78 | 37,7 |
| b2L_h2.5_adapt_brut | 1,148 | 0,000 | oui | 87 | 9 | 1 | 1,00 | 0,60 | 1,466 | 0,75 | 95,5 |
| b2L_h2.5_adapt_corr | 1,148 | 0,000 | oui | 87 | 9 | 1 | 1,00 | 0,60 | 1,464 | 0,75 | 904,1 |
| b2L_h1.25_adapt_brut | 1,240 | 0,001 | oui | 405 | 17 | 1 | 1,00 | 0,81 | 1,502 | 0,75 | 621,1 |
| b2_h10_adapt_brut | 1,223 | 0,698 | non | 8 | 4 | 2 | 0,67 | 0,00 | 1,056 | 0,67 | 1,0 |
| b2_h10_adapt_corr | 1,223 | 0,697 | non | 8 | 4 | 2 | 0,67 | 0,00 | 1,056 | 0,67 | 1,4 |
| b2_h5_adapt_brut | 1,085 | 0,179 | non | 22 | 7 | 2 | 0,81 | 0,36 | 1,228 | 0,83 | 3,3 |
| b2_h5_adapt_corr | 1,085 | 0,181 | non | 22 | 7 | 2 | 0,81 | 0,36 | 1,228 | 0,83 | 10,2 |
| b2_h2.5_adapt_brut | 1,150 | 0,000 | oui | 87 | 9 | 1 | 1,00 | 0,60 | 1,466 | 0,75 | 19,7 |
| b2_h2.5_adapt_corr | 1,150 | 0,000 | oui | 87 | 9 | 1 | 1,00 | 0,60 | 1,464 | 0,75 | 183,3 |
| b2_h1.25_adapt_brut | 1,242 | 0,002 | oui | 405 | 17 | 1 | 1,00 | 0,81 | 1,502 | 0,75 | 123,9 |
| b2_h1.25_adapt_corr | 1,240 | 0,005 | oui | 405 | 17 | 1 | 1,00 | 0,81 | — (arrêté) | 0,75 | 1 200 (rc 124) |
| b2_h10_intri_brut | 1,192 | 0,019 | oui | 10 | 4 | 2 | 1,00 | 0,00 | 1,416 | 0,67 | 1,1 |
| b2_h10_intri_corr | 1,193 | 0,018 | oui | 10 | 4 | 2 | 1,00 | 0,00 | 1,418 | 0,67 | 1,6 |
| b2_h5_intri_brut | 1,054 | 0,008 | oui | 24 | 5 | 1 | 1,00 | 0,33 | 1,178 | 1,00 | 3,6 |
| b2_h5_intri_corr | 1,054 | 0,006 | oui | 24 | 5 | 1 | 1,00 | 0,33 | 1,177 | 1,00 | 9,0 |
| b2_h2.5_intri_brut | 1,073 | 0,011 | oui | 122 | 11 | 2 | 1,00 | 0,64 | 1,570 | 0,90 | 30,6 |
| b2_h2.5_intri_corr | 1,073 | 0,006 | oui | 124 | 12 | 2 | 1,00 | 0,64 | 1,628 | 0,91 | 124,2 |
| b2_h1.25_intri_brut | 1,086 | 0,059 | non | 371 | 17 | 2 | 1,00 | 0,75 | 1,345 | 0,81 | 130,6 |
| b2_h1.25_intri_corr | 1,086 | 0,060 | non | 371 | 17 | 2 | 1,00 | 0,75 | 1,345 | 0,81 | 749,9 |

Chronologie de l'adaptatif (`history.csv`) : la première insertion a lieu à 1,002 à
1,015 f_t, quelle que soit la maille. Les insertions s'accumulent ensuite sans rupture
jusqu'au pic : 4, 15, 81 et 393 joints insérés au pic pour h = 10, 5, 2,5 et 1,25 mm. Le
premier joint ne rompt qu'après le pic, quand σ est retombé à 0,14 à 0,34 f_t. Le pas de
temps est le même dans les deux schémas (2,17 ns à h = 1,25 mm), car la pénalité
d'insertion vaut aussi 100 E.

| critère | adaptatif (b2L) | intrinsèque (b2) |
|---|---|---|
| B2-a pic ∈ [0,95 ; 1,10] | conforme à h = 5 mm seulement (1,084) ; non conforme à 10 ; 2,5 ; 1,25 mm (1,22 ; 1,15 ; 1,24) | conforme à 5 ; 2,5 ; 1,25 mm (1,05 ; 1,07 ; 1,09) ; non conforme à 10 mm (1,19) |
| B2-b convergence du pic (≤ 3 %) | non conforme : +8,0 % entre 2,5 et 1,25 mm | conforme : +1,2 % |
| B2-c énergie ∈ [1,0 ; 1,5] | conforme à 2,5 mm (1,466) ; non conforme à 10 ; 5 ; 1,25 mm (1,74 ; 2,02 ; 1,502, ce dernier à la limite) | conforme à 10 et 5 mm (1,42 ; 1,18) ; non conforme à 2,5 mm (1,57 / 1,63) ; non jugeable à 1,25 mm (non séparée) |
| B2-d convergence de l'énergie (≤ 10 %) | conforme : 1,466 → 1,502, +2,5 % | non jugeable (maille fine non séparée) ; 5 → 2,5 mm : +33 % |
| B2-e une fissure traversante, ≤ 10 % des rompus hors d'elle | conforme à 10 ; 2,5 ; 1,25 mm ; non conforme à 5 mm (deux fissures, 60 / 40 %) | conforme à 5 et 1,25 mm (1 rompu sur 17 hors fissure) ; non conforme à 10 et 2,5 mm |
| B2-f localisation des insertions / endommagements | 10 mm : localise (0,00) ; 5 mm : mixte (0,27) ; 2,5 mm : nuclée partout (0,60) ; 1,25 mm : nuclée partout (0,81) | 10 mm : localise ; 5 mm : mixte (0,33) ; 2,5 mm : nuclée partout (0,64) ; 1,25 mm : nuclée partout (0,75) |
| B2-g contact corrigé (pic et W à ≤ 1 %, travail du contact ≤ 0) | conforme : pic identique, W ±0,15 %, travail −4e-9 à −4e-11 J/m | conforme à 10 ; 5 ; 1,25 mm ; non conforme à 2,5 mm (W +3,7 %, 11 → 12 rompus) ; travail ≤ 0 partout |

Coût du correctif : ×1,3 à ×5 sur les maillages grossiers, ×9 en adaptatif à 2,5 mm
(95 → 904 s), plus de ×10 en adaptatif à 1,25 mm (124 s sans correctif, arrêté à 1 200 s avec
le correctif à 89 % de T, après la rupture), ×5,7 en intrinsèque à 1,25 mm (131 → 750 s).

## Banc 3 — plaque à entaille latérale (SENT), 40 × 80 mm, a = 10 mm

Maillages `meshes/sent_h*.msh` (`scripts/mesh_sent.py`, 1 621 et 5 593 triangles).
h = 2,5 mm : contact corrigé, T = 2 ms ; les deux calculs atteignent 20 min (rc 124, arrêt à
91 % et environ 70 % de T, après la rupture : dépouillement sur la dernière trame, à 1,6 et
1,2 ms ; énergie non disponible, faute de résumé). h = 1,25 mm : contact sans correctif,
T = 1 ms (addenda 4 et 5). σ est rapportée à la section brute (W = 40 mm). « N_fiss » :
composantes connexes brutes / après fusion des composantes distantes de moins de 1,5 h (la
fissure adaptative h = 1,25 mm saute un triangle à x = 27 mm). W est rapporté à G_I × ligament
(30 mm).

| calcul | σ_pic/f_t | σ_fin/f_t | insérés / D>0 | rompus | N_fiss | 1er rompu – pointe [mm] | x_max [mm] | écart moyen / max à y_entaille [mm] | L/ligament | propag. (rupture) | propag. (insertion) | W/(G_I·lig.) | durée [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| b3_h2.5_adapt | 0,694 | 0,000 | 127 | 19 | 1 | 0,52 | 40,0 | 2,23 / 4,42 | 1,074 | 0,78 | 0,63 | — | 1 200 (rc 124) |
| b3_h2.5_intri | 0,609 | 0,014 | 145 | 19 | 1 | 1,38 | 40,0 | 2,23 / 4,42 | 1,074 | 0,78 | — | — | 1 200 (rc 124) |
| b3_h1.25_adapt | 0,740 | 0,015 | 526 | 29 | 2 / 1 | 1,40 | 40,0 | 0,64 / 2,42 | 1,082 | 0,75 | 0,58 | 1,786 | 473 |
| b3_h1.25_intri | 0,661 | 0,062 | 393 | 30 | 1 / 1 | 0,75 | 40,0 | 1,06 / 2,61 | 1,061 | 0,79 | — | 1,652 | 900 |

Chronologie : en adaptatif, la première insertion a lieu dès σ = 0,24 à 0,25 f_t (brut), à
cause de la concentration en pointe. 57 et 161 joints sont insérés au pic, et aucun n'est
rompu avant le pic dans les deux schémas.

| critère | adaptatif | intrinsèque |
|---|---|---|
| B3-a amorçage à < 2h de la pointe | conforme (0,52 et 1,40 mm) | conforme (1,38 et 0,75 mm) |
| B3-b traversée de la pointe au bord droit | conforme à 2,5 mm ; à 1,25 mm, conforme après fusion de deux composantes séparées d'un triangle (non conforme en topologie stricte) | conforme |
| B3-c écart moyen ≤ 4 mm, maximal ≤ 8 mm | conforme (2,23 / 4,42 ; 0,64 / 2,42) | conforme (2,23 / 4,42 ; 1,06 / 2,61) |
| B3-d longueur / ligament ∈ [1,0 ; 1,3] | conforme (1,074 ; 1,082) | conforme (1,074 ; 1,061) |
| B3-e propagation ≥ 80 %, hors fissure ≤ 10 % | non conforme de peu (0,78 ; 0,75) ; hors fissure 0 % | non conforme de peu (0,78 ; 0,79) ; hors fissure 0 % |

À h = 2,5 mm, la fissure principale est le même chemin de 19 joints dans les deux schémas (même
écart, même longueur). Seuls le pic (+14 %) et les instants diffèrent.

## Banc 4 — sensibilité de l'adaptatif au critère (barre du banc 2, contact corrigé, T = 3 ms)

Le calcul de base reproduit `b2L_h*_adapt_corr` exactement (1,084 / 26 / 10 / 2,018 à 5 mm ;
1,148 / 87 / 9 / 1,464 à 2,5 mm). Les journaux confirment que chaque clé est active (ligne
`[FDEM]` d'annonce).

| calcul | σ_pic/f_t | insérés | rompus | N_fiss | part de la fissure princ. | hors fiss. | propag. (rupture) | propag. (insertion) | W/(G_I·b) | durée [s] |
|---|---|---|---|---|---|---|---|---|---|---|
| b4_h5_base | 1,084 | 26 | 10 | 2 | 0,60 | 0,27 | 0,78 | 0,72 | 2,018 | 30 |
| b4_h5_elliptic | 1,087 | 24 | 8 | 2 | 0,50 | 0,38 | 0,86 | 0,70 | 1,381 | 36 |
| b4_h5_volume | 1,083 | 22 | 7 | 2 | 0,71 | 0,36 | 0,83 | 0,67 | 1,301 | 36 |
| b4_h5_max | 1,080 | 24 | 6 | 2 | 0,83 | 0,38 | 0,80 | 0,65 | 1,309 | 46 |
| b4_h5_hold5 | 1,083 | 26 | 10 | 2 | 0,60 | 0,27 | 0,78 | 0,72 | 2,018 | 28 |
| b4_h5_tip1.6 | 1,084 | 27 | 10 | 2 | 0,60 | 0,22 | 0,78 | 0,73 | 2,052 | 28 |
| b4_h2.5_base | 1,148 | 87 | 9 | 1 | 1,00 | 0,60 | 0,75 | 0,53 | 1,464 | 530 |
| b4_h2.5_elliptic | 1,145 | 89 | 9 | 1 | 1,00 | 0,60 | 0,75 | 0,55 | 1,470 | 568 |
| b4_h2.5_volume | 1,149 | 87 | 9 | 1 | 1,00 | 0,60 | 0,75 | 0,53 | 1,467 | 518 |
| b4_h2.5_max | 1,142 | 100 | 9 | 1 | 1,00 | 0,61 | 0,75 | 0,58 | 1,485 | 572 |
| b4_h2.5_hold5 | 1,148 | 86 | 9 | 1 | 1,00 | 0,60 | 0,75 | 0,53 | 1,463 | 497 |
| b4_h2.5_tip1.6 | 1,148 | 88 | 9 | 1 | 1,00 | 0,59 | 0,75 | 0,54 | 1,464 | 510 |

| variante | prédiction écrite | h = 5 mm | h = 2,5 mm | verdict « sensible » (pic > 2 % ou rompus / N_fiss / W > 10 %) |
|---|---|---|---|---|
| elliptic | pic ≤ base, ≥ 0,9 base ; plus d'insérés | pic +0,3 %, 24 insérés : prédiction fausse sur les deux points | pic −0,3 % ; 89 insérés (+2) | sensible à 5 mm (W −32 %, rompus −20 %) ; insensible à 2,5 mm |
| volume | neutre, pic ± 1 %, rompus ± 20 % | pic −0,1 % ; rompus −30 % | identique (W +0,2 %) | sensible à 5 mm (W −36 %) ; insensible à 2,5 mm |
| max | pic ≤ base ; plus d'insérés | pic −0,4 % ; 24 insérés (−2) | pic −0,5 % ; 100 insérés (+15 %) | sensible à 5 mm (W −35 %, rompus −40 %) ; insensible à 2,5 mm (W +1,4 %) |
| hold 5 | pic ± 1 %, motif inchangé | conforme (identique) | conforme (identique) | insensible |
| tip 1,6 | N_fiss ≤ base, propagation ≥ base, pic ± 2 % | conforme (N_fiss 2 = 2, propagation égale, pic identique) | conforme (identique) | insensible |

À 5 mm, les trois variantes « sensibles » font toutes passer la barre de deux fissures
(W = 2,02 G_I·b) à une fissure dominante (W = 1,30-1,38) : c'est un basculement du nombre de
fissures sur un maillage à quatre éléments dans la largeur. À 2,5 mm, plus aucune variante ne
change ni le pic, ni le nombre de rompus, ni la fissure.

## Réponses aux questions

L'adaptatif localise-t-il sur une fissure ou nuclée-t-il partout ? Les deux, à deux
instants différents. La rupture finale est localisée : une fissure traversante dès
h ≤ ℓ_cz/8, sans aucun joint rompu hors d'elle, plus nettement que l'intrinsèque (deux
composantes à 2,5 et 1,25 mm). En revanche, l'insertion est diffuse avant le pic, et de plus
en plus diffuse quand la maille s'affine : 27 % des joints de la barre sont insérés au pic à
ℓ_cz/16, 81 % des insérés sont à plus de 2h de la fissure, et la part de propagation à
l'insertion n'est que de 53 à 61 %. Dans un champ uniforme, toutes les facettes proches de
l'horizontale franchissent f_t ensemble et l'adaptatif les insère toutes. Sur 405 joints
insérés, 17 rompent ; les autres restent endommagés (D > 0) loin de la fissure. On retrouve
le « tapis » d'insertions décrit dans BILAN_insertion_adaptative.md §1, cette fois sur un
banc où la réponse attendue est connue.

Converge-t-il en maille ? Pour l'énergie, oui (1,466 → 1,502 G_I·b entre ℓ_cz/8 et ℓ_cz/16,
+2,5 %). Pour le pic, non : il remonte de 1,08 à 1,15 puis 1,24 f_t de ℓ_cz/4 à ℓ_cz/16, alors
que la première insertion reste à 1,00 f_t. Le pic mesure l'instant où le tapis de joints
insérés et endommagés finit par localiser ; plus il y a de facettes, plus ce tapis porte de
charge. L'intrinsèque converge sur le pic (+1,2 %) ; sur l'énergie, il n'est pas jugeable à
la maille fine. Sur la plaque entaillée, le pic adaptatif monte de 6,6 % entre 2,5 et
1,25 mm, l'intrinsèque de 8,5 % : ni l'un ni l'autre n'a convergé à ces tailles. La
trajectoire, elle, est stable et identique dans les deux schémas.

Le contact corrigé change-t-il quelque chose ici ? Pas sur ces bancs quasi statiques : le
travail du contact reste de l'ordre de 1e-10 J/m, l'adaptatif garde exactement le même pic et le même motif (énergie ±0,15 %), et un seul calcul intrinsèque (2,5 mm) bascule d'un joint rompu. Son coût
est en revanche élevé : de ×5 à plus de ×10 sur les maillages fins.

Le critère compte-t-il ? À 2,5 mm, non : `elliptic`, `volume`, `max`, `insertionHoldSteps = 5`
et `insertionTipFactor = 1,6` donnent tous le même pic (±0,5 %), les mêmes 9 rompus et la même
énergie (±1,4 %). Les choix de critère ne corrigent pas le tapis d'insertions (86 à 100 insérés
pour 9 rompus). Ce point rejoint le bilan d'août : le tapis vient du champ, pas du critère.

## Réserves

- Ordre de procédure : le banc 1 a été lancé une fois avant l'écriture de CRITERES.md (relance
  identique au bit près ; écart déclaré dans CRITERES.md).
- Quatre modifications de protocole ont été faites en cours de campagne, chacune déclarée avant
  les calculs concernés : variante stable `b1s` et trajet mixte à mors bloqués, barre ramenée
  de 40 à 30 mm (erreur d'un facteur 2 dans le calcul d'énergie élastique de CRITERES.md),
  prolongation à T = 3 ms des calculs adaptatifs, contact sans correctif et T = 1 ms pour la
  plaque à h = 1,25 mm. Aucun seuil n'a été modifié.
- Trois calculs ont atteint la limite de 20 min et un a été arrêté à la main : `b2_h1.25_adapt_corr` (89 % de T),
  `b3_h2.5_adapt` (91 %), `b3_h2.5_intri` (environ 70 %) et `b3_h1.25_adapt` avec correctif
  (arrêté à 24 %, remplacé par le calcul sans correctif). Les trois premiers ont dépassé la
  rupture, mais leur énergie n'est pas disponible. `b2L_h1.25_adapt_corr` n'a pas été lancé.
- La barre intrinsèque h = 1,25 mm n'est pas séparée à T = 1,25 ms (σ_fin = 0,06 f_t) ; un
  calcul plus long aurait dépassé 20 min (750 s pour 1,25 ms avec le correctif).
- Le pic est lu sur la réaction des mors, avec un bruit inertiel de ±2 % (Cundall 0,7,
  pullV = 0,02 m/s). Les écarts de pic inférieurs à 2 % ne sont pas significatifs.
- N_fiss repose sur le partage d'un sommet initial. Une fissure qui traverse un triangle sans
  joint rompu (b3 h = 1,25 mm adaptatif) compte pour deux composantes ; la fusion à 1,5 h n'a
  été appliquée qu'aux calculs b3 h = 1,25 mm et b4 h = 2,5 mm, car les `.vtu` des autres
  étaient déjà supprimés quand elle a été ajoutée.
- « Endommagé D > 0 » (intrinsèque) et « inséré » (adaptatif) sont comparés comme deux mesures
  du même seuil d'enveloppe : en adaptatif, nInserted = nDamaging à chaque instant jusqu'au pic. Ce n'est
  pas strictement le même objet.
- Un seul tirage de maillage par taille (graine gmsh 1) : la dispersion d'un tirage à l'autre
  n'est pas mesurée. La bifurcation une/deux fissures à 5 mm montre qu'elle n'est pas
  négligeable sur les maillages grossiers.
- Les calculs ont tourné à côté du rejeu tunnel (jusqu'à trois processus à deux fils sur quatre
  cœurs) : les durées sont des durées murales, majorées par ce partage.

## Fichiers

- `CRITERES.md` ; `RESULTATS.md` ; `fig_synthese_insertion.pdf` (et `.png`).
- `scripts/` : `mesh_b1.py`, `mesh_sent.py`, `make_decks_b1.py`, `make_decks_b234.py`,
  `run_queue.sh`, `analyse_b1.py`, `analyse_frac.py`, `fig_synthese.py`.
- `meshes/` : `b1_t0.msh`, `b1_t45.msh`, `bar_h*.msh`, `sent_h*.msh`.
- `b1/`, `b1s/`, `b2/`, `b2L/`, `b3/`, `b4/` : decks `.cfg`, journaux `.log`, durées `.dur`,
  `<banc>_resultats.csv`, courbes (`b1*_courbes.npz`) et segments de fissures
  (`<banc>_fissures.npz`, `*_depouille.npz`).
- La première série `b1` (avant CRITERES.md) est dans `b1_premier/` (journaux seulement).
- Dans cette copie, `meshFile` des decks pointe vers `docs/rapport_guide/bancs_insertion/meshes/` ;
  les journaux `.log` gardent les chemins du dossier de travail d'origine (scratchpad).
