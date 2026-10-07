# Surface initiale et coupe a -0,5 mm : comparaison a 165 us

Run `out_stanne2025_rock137`, trame 11, temps exact 164,995 us. Post-traitement seul ; aucune modification du run.

La surface initiale est le plan z = 0 dans le repere de l'impact (z = 0,15 m dans les VTU). Pour chaque facette reellement rompue (`tBreak >= 0`, hors interfaces `bonded`), on retient son arete dont les deux sommets appartiennent a ce plan dans la geometrie initiale. Les doublons geometriques sont elimines. Une rencontre en un sommet seul ne donne pas une ligne. Aucun triangle coplanaire n'a ete trouve.

Le filtre precedent exigeait des sommets strictement de part et d'autre du plan et donnait zero segment a z = 0. La nouvelle extraction donne 301 aretes uniques. Elle est stable pour des tolerances de 1e-7, 1e-6 et 1e-5 mm. Un essai synthetique verifie une rencontre sur arete, sur sommet seul, une arete dupliquee et la detection d'un triangle coplanaire ; le filtre strict sert de controle negatif.

| Operation | Segments | Somme des longueurs de traces | Rayon maximal des traces |
|---|---:|---:|---:|
| Aretes de fissures sur la surface initiale | 301 | 311,62 mm | 14,84 mm |
| Intersection avec z = -0,5 mm, geometrie initiale | 700 | 402,06 mm | 14,56 mm |
| Memes aretes de surface aux positions courantes exportees | 301 | 311,61 mm (3D) | 14,85 mm (XY) |

Le nombre de segments ne compte pas les fissures : le decoupage differera entre une coupe interieure et les aretes de bord. La somme de leurs longueurs n'est pas la longueur d'une radiale. Le rayon maximal des traces n'est ni un rayon de cratere ni la longueur curviligne d'une fissure. Il differe de la distance centre-pointe de 16,26 mm du reseau 3D, car cette pointe peut se trouver sous la surface.

**Observation visuelle :** les principales radiales debouchent effectivement sur la surface initiale. Elles persistent sur la coupe a -0,5 mm avec des decalages et des ramifications locales. La coupe sous la surface est plus dense dans le noyau. La geometrie courante projetee en XY est voisine ; cela ne signifie pas que le deplacement vertical soit nul (deplacement 3D maximal des sommets suivis : 0,828 mm).

Les arêtes du panneau de surface sont des traces de joints rompus, et non le filaire de tout le maillage. Le reseau triangule central est donc une observation de nombreuses ruptures a la surface dans ce modele. Aucune ouverture minimale n'a ete imposee : une interface rompue peut etre refermee ou tres peu ouverte, et sa trace n'est pas necessairement visible sur une photographie experimentale.

Le panneau courant suit les aretes qui appartenaient a la surface initiale. Ce n'est pas une reconstruction des nouvelles surfaces visibles apres ejection, ni une evaluation des occultations par l'insert ou par les fragments. Les couleurs des modes sont indicatives (`jointBreakModeRef = slipF` dans le run).

Figure PNG et PDF : `output/pdf/stanne_surface_165us/comparaison_surface_coupe`. Mesures et hashes des deux VTU : `output/pdf/stanne_surface_165us/mesures.json`. Script : `tools/fig_surface_compare.py`, parametres `--run`, `--frame`, `--out`.
