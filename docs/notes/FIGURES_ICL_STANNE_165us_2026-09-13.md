# St Anne : fissures radiales et rendu comparable aux figures ICL

Analyse du 13 septembre 2026, run `out_stanne2025_rock137`, g1y19. Instant fige : trame 11, **164,995 microsecondes**, et non 170 microsecondes. Les trames 6, 8, 9, 10 et 11 servent a la chronologie. Aucun changement du solveur, des decks ou des sorties du run actif.

## Ce que les donnees permettent de dire

**Des branches de rupture radiales connectees existent dans ce run.** Le filtre utilise `tBreak >= 0` et exclut les interfaces `bonded`. A 165 us : 5 023 facettes rompues, 59 composantes connexes ; la composante principale en contient 4 928. L'adjacence est definie par une arete partagee dans la geometrie initiale, apres reunification des sommets dedoubles. La partition a ete verifiee par deux algorithmes independants (graphe creux et union-find de `crack_paths.py`).

Pour un noyau de rayon 6 mm, trois grandes branches satisfont les criteres geometriques de radialite du script (`moyenne |n.e_r| < 0,35` et `moyenne |n_z| < 0,6`, ponderees par les aires) :

| Branche | Facettes | Distance centre-pointe | Azimut |
|---|---:|---:|---:|
| A | 99 | 16,26 mm | 173 degres |
| B | 51 | 14,84 mm | 94 degres |
| C | 35 | 13,06 mm | -107 degres |

Ce sont des distances centre-pointe, pas des longueurs curvilignes ni des longueurs de fissures exclusivement en surface. Les nombres de facettes portent sur des surfaces ramifiees, pas une chaine de 99 facettes consecutives. Avec un rayon de noyau de 8 puis 10 mm, la branche A reste radiale et atteint toujours 16,26 mm (51 puis 24 facettes hors noyau). La branche B persiste egalement. Les classements des petites branches sont plus sensibles au seuil.

Cela depasse l'indice visuel de quelques triangles sur le maillage grossier : **la localisation peripherique est maintenant visible et connectee a cette resolution**. Il reste a verifier sa convergence et les criteres quantitatifs de Yang. Le noyau conserve un reseau dense : le rendu ne fait pas disparaitre cette observation.

## Pourquoi les anciens plots donnent une impression differente

1. **Aretes de toutes les facettes** : le wireframe dessine les limites numeriques des triangles, y compris au milieu d'une meme surface de fissure. Une projection melange aussi les facettes profondes et superficielles. Des milliers de traits donnent une impression de grillage.
2. **Surfaces opaques** : la rangee superieure de la fig. 14 ICL utilise visuellement des facettes remplies et ombrees. Dans les nouveaux rendus, les aretes sont desactivees, les surfaces sont opaques et les facettes cachees sont occultees par le rendu 3D. Aucun lissage ni suppression de fissures.
3. **Coupe versus tranche** : `|y| < 5 mm` est une tranche de 10 mm d'epaisseur projetee, pas une coupe. Une intersection triangle-plan produit des lignes fines et revele les branches. Les nouveaux plans sont explicites : `z = -0,5 mm`, `y = 0`, `x = 0`, `z = -3 mm`.
4. **References du mauvais essai** : les anciens cercles "Yang crater 7" et "Yang radiales 10" proviennent du Kuru a 9 m/s. Ils ne conviennent pas au St Anne a 10,66 m/s et ont ete retires.
5. **Geometrie initiale versus courante** : les figures livrees utilisent la geometrie initiale afin de ne pas confondre propagation de fissure et deplacement des fragments. Une coupe a profondeur fixe dans la geometrie deformee peut avoir une cavite centrale differente. Les nouvelles coupes ne sont pas une extraction de la surface libre, ni une mesure du cratere.
6. **Orientation versus physique** : qualifier toutes les facettes subverticales de "fissures" et les autres de "broye" n'est pas un diagnostic mecanique. Des fissures laterales sont plutot horizontales ; ce run a `bulkDamage = off`. Le gris de la figure des branches signifie uniquement "reste du reseau", pas pulverisation.

## Figures livrees

Dans `output/pdf/stanne_ICL_165us/`, chaque figure existe en PNG haute resolution et PDF vectoriel :

- `01_evolution_ICL` : chronologie 90, 120, 135, 150, 165 us ; surfaces vues de dessus, coupe horizontale puis coupe verticale. Cadrage constant dans chaque rangee ; barres d'echelle sur les coupes. Les vues 3D ont leur propre cadrage.
- `02_effet_du_rendu` : memes facettes en aretes projetees, surfaces opaques et coupes exactes.
- `03_branches_connectees` : vue 3D complete et identification de trois branches, avec distances centre-pointe.
- `mesures_et_provenance.json` : mesures, sensibilite au rayon du noyau et SHA-256 des trames utilisees.

Les couleurs rouge/jaune suivent `breakMode` exporte. **Le deck effectif porte encore `jointBreakModeRef = slipF`** : la partition traction/cisaillement est indicative tant que la normalisation a la rupture n'est pas corrigee. La connectivite et les distances ne dependent pas de cette etiquette. Les couleurs A/B/C de la troisieme figure indiquent des branches, pas des modes de rupture.

Script reproductible : `tools/fig_icl_stanne_review.py`. Depuis la racine du depot :

```powershell
python tools/fig_icl_stanne_review.py --run out_stanne2025_rock137 --frame 11 --out output/pdf/stanne_ICL_165us
```

## Sources visuelles consultees

- `bibliographie/main.pdf`, Yang et al. 2025, page 11, figures 14-16 : **reference materiau et vitesse pertinente**, St Anne, piston 10,66 m/s.
- `bibliographie/ICL.pdf`, Yang et al. 2026, page 10, figures 14-15 : reference de presentation complementaire, mais Kuru, autres vitesses et pulverisation active.

Les nouvelles figures sont inspirees de leur presentation ; la coupe horizontale a -0,5 mm est notre choix explicite de visualisation. Ce n'est pas une affirmation que leur rangee de surface emploie exactement ce filtre. Aucune donnees de Yang n'a ete fabriquee ou superposee.
