# Données du run St Anne 2025 sur maillage rock137 — archive du 14/09/2026

Les **données texte** du run `out_stanne2025_rock137`, terminé le 14/09 à 12 h 05 après 29 h 21 de
calcul sur 14 fils. 4,2 Mo, contre 2,0 Go pour le dossier de sortie complet.

Ce sont **ces fichiers qui portent tous les chiffres** cités dans
`docs/COMPARAISON_yang2025_stanne_2026-09-14.md` et `docs/ECARTS_guo2014_rockim_2026-09-13.md` §10.
Toutes les courbes se refont à partir d'eux seuls, sans les trames.

| Fichier | Taille | Contenu |
|---|---:|---|
| `history.csv` | 0,43 Mo | 2 022 pas enregistrés, 38 colonnes : temps, forces de contact par paire, positions et vitesses de chaque corps, contrainte à la jauge, compteurs de rupture, postes d'énergie |
| `frames.csv` | 1 ko | l'instant de chacune des 22 trames |
| `fdem3d_final_elements.csv` | 3,73 Mo | l'état élément par élément à 300 µs |
| `config_effective.cfg` | 0,02 Mo | **le deck tel que le solveur l'a lu**, defaults résolus compris |
| `journal_solveur.log` | 0,01 Mo | sortie du solveur, dont le bilan d'énergie final et les avertissements |

`config_effective.cfg` est le fichier à citer dans un rapport, pas le deck source : il dit ce qui a
réellement tourné, y compris les valeurs par défaut que le deck ne mentionne pas.

## Refaire une courbe

```
python tools/fig_fp_direct.py    out_stanne2025_rock137 --stem results/fig/x_fp
python tools/fig_kinetics.py     out_stanne2025_rock137 --stem results/fig/x_cin
python tools/fig_retournement.py out_stanne2025_rock137 --stem results/fig/x_ret
```

Ces trois-là ne lisent que `history.csv`. Il suffit donc de recréer un dossier
`out_stanne2025_rock137/` contenant les fichiers de ce répertoire.

Les figures de **géométrie** (joints, coupes, cratère, branches, contrainte, revues ICL et surface)
lisent les trames VTU, qui ne sont pas ici. Elles se rejouent en relançant le calcul,
`docs/REPRODUIRE_stanne_radiales_2026-09-14.md`.

## Le bilan d'énergie, pour mémoire

Énergie cinétique 60,12 → 3,28 J, résidu du budget −9,8e-8 J, soit 1,5e-7 % de l'échelle.
Postes : rupture 34,5 J, frottement 27,8 J, joints 14,5 J, élastique stocké 18,1 J.
12 090 facettes rompues sur 209 380, 4 089 fragments, volume détaché 2 153 mm³.
