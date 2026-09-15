# Configuration et maillage du run St Anne 2025 — 15/09/2026

Les deux planches qui décrivent le calcul **avant** ses résultats. Elles se refont depuis le seul
fichier de maillage, sans les trames et sans relancer :

```
python tools/fig_montage_impact.py meshes/impact_yang_train1_rock137_hxt.msh \
       --out results/fig/config_stanne/montage.png --v0 10.66 --cas stanne2025 \
       --dt 1.70663e-9 --spp 0.588
python tools/fig_mesh3d.py meshes/impact_yang_train1_rock137_hxt.msh \
       --out results/fig/config_stanne/maillage.png --impact 0 0 0 --zoom 0.03
```

Le maillage n'est pas dans le dépôt, il se régénère à l'identique, graine fixée :
`python tools/make_impact_mesh.py meshes/impact_yang_train1_rock137_hxt.msh 1.0 2e-5 1.37 gap=2e-5 quality=hxt train=fixed`

## `montage.png` — l'assemblage et sa chronologie

Trois panneaux. **(a)** la chaîne à six corps, chacun d'une couleur : piston, bit, insert carbure,
circlip, plaque, roche. **(b)** le zoom sur l'insert hémisphérique de rayon 8,51 mm posé sur la
roche, avec le jeu initial de 0,02 mm. **(c)** la chronologie, avec le coût en pas de temps.

| Événement | Instant |
|---|---:|
| Le piston part, jeu 0,02 mm | 0 µs |
| Il touche le bit | 1,9 µs |
| L'onde traverse le bit, 242 mm à 5 048 m/s | 47,9 µs |
| **L'onde charge la roche** | **50 µs** |
| Fin de la phase de charge, Yang §5.1 | 254 µs |
| Arrêt des médianes, des radiales, des latérales | 291, 388, 482 µs |

Les cinquante premières microsecondes sont un préambule où **rien n'atteint la roche**. C'est ce que
le panneau (c) rend visible, et c'est la raison pour laquelle un run court ne montre jamais rien.

⚠️ Le paramètre `--cas` choisit l'article qui fournit les repères. `stanne2025` pour ce run.
Le défaut, `kuru2026`, pose les repères du granite et mettrait une **fausse citation** sous la
figure : c'est le défaut corrigé le 15/09.

## `maillage.png` — la gradation, et son contrôle

Trois panneaux. **(a)** une coupe **exacte** par le plan médian : chaque tétraèdre traversé est
réellement découpé, pas approximé par une tranche d'épaisseur fixe qui mélangerait les tailles.
**(b)** le zoom sur les 30 mm autour de l'impact. **(c)** la gradation **réalisée** en fonction de la
distance au point d'impact, médiane et bande 5-95 %. Ce troisième panneau est le contrôle falsifiant :
il montre la taille obtenue, pas celle demandée.

| | Tout le maillage | Cœur, r < 12,5 mm |
|---|---:|---:|
| Tétraèdres, roche | 90 878 | **14 030** |
| Arête médiane, roche | 5,01 mm | **1,42 mm** |
| Diamètre inscrit médian, roche | 1,74 mm | 0,49 mm |

Total 108 667 tétraèdres et 22 607 nœuds sur les six corps.

**Qualité** : médiane 0,80 sur l'échelle où 1 est le tétraèdre régulier, minimum 0,176, et **aucun
élément sous 0,10**. Le pas de temps est commandé par le diamètre inscrit minimal, 0,2528 mm, qui se
trouve dans l'**insert** en carbure et non dans la roche. Trente des trente-cinq éléments sous 0,3 mm
sont dans l'insert : c'est lui qui paie le pas de temps de 1,71 ns, pas le maillage de roche.

## Coût

0,588 s par pas mesuré sur 14 fils, soit 29 h 21 pour les 300 µs du run, 179 641 pas.
Atteindre les 482 µs où Yang voit s'arrêter les latérales coûterait environ 47 h de plus.
