# Où sont les figures, et ce que chacune montre

Écrit le 13/09, **mis à jour le 14/09 à midi** quand le run St Anne s'est terminé.

435 figures sont versionnées, en PNG de relecture et en **PDF vectoriel Computer Modern**
directement inclusible en LaTeX. Elles vivent aux deux mêmes endroits dans les deux dépôts :

| | Dépôt `Fuzquial/rockim`, branche `g1` | Archive `phd_geothermie` |
|---|---|---|
| chemin | `results/fig/`, `output/pdf/` | `FDEM/rockim_g1/results/fig/`, `.../output/pdf/` |

Prendre toujours le `.pdf`, jamais le `.png`, pour un document LaTeX.

## 1. L'ÉTAT FINAL du run St Anne 2025 : `results/fig/stanne300/`

**Mis à jour le 14/09 à midi, quand le run s'est terminé.** Trame 21, t = 300 µs, maillage rock137
(108 667 tétraèdres, cœur à 1,42 mm). **C'est le jeu à citer par défaut** ; les dossiers `stanne192`
à `stanne295` sont des états intermédiaires, gardés pour tracer une évolution et pour rien d'autre.

| Fichier | Ce qu'il montre |
|---|---|
| `evolution.pdf` | **les mesures du cratère contre le temps, sur quatorze trames** : rayon de peau et extension totale, profondeur, volume détaché, portée des radiales, avec le retournement en trait vertical. La planche de synthèse |
| `retournement.pdf` | la datation du retournement, 284,3 µs, avec la force de contact en second axe et la pause de 260 µs marquée d'une croix rouge. **Lire la légende** : un croisement de vitesse sous charge maintenue n'est pas un retournement |
| `joints.pdf` | les 12 090 facettes rompues seules, six vues, couleur par mode de rupture |
| `coupes.pdf` | coupes exactes en y = 0, 5, 10 mm et z = −0,5, −3, −8 mm ; la trace est l'intersection facette/plan, pas une projection |
| `sig_coupe.pdf` | σ₁ en coupe verticale, lissée, roche seule — la chute de contrainte sur les lèvres des radiales s'y lit |
| `sig_dessus.pdf` | σ₁ en vue de dessus à z = −2 mm, non lissée |
| `cratere.png` | vue de dessus sectorisée en 16, avec la portée des radiales secteur par secteur |
| `branches.pdf` | branches connectées, classées radiale / conique / horizontale ; chiffres dans les deux `.csv` à côté |
| `fp.pdf` | force-pénétration, courbe mesurée au contact **et** courbe estimée par le train rigide |
| `cinetique.pdf` | vitesses d'indentation et de rebond à la manière de Yang (pente des portions linéaires), plus la jauge |

Deux compagnons dans `output/pdf/stanne_300us/` : `ICL/` (revue en trois planches, avec
`mesures_et_provenance.json` qui porte les empreintes SHA-256 des trames lues) et `surface/`
(comparaison des traces en surface et sous la surface, avec ses mesures).

Les chiffres qui vont avec ces planches sont dans
`docs/COMPARAISON_yang2025_stanne_2026-09-14.md`, et les données brutes dans
`results/data/stanne2025_rock137/`.

## 2. Les états antérieurs du même run

Dix jeux intermédiaires, `stanne192` `stanne210` `stanne225` `stanne240` `stanne255` `stanne266`
`stanne270` `stanne277` `stanne285` `stanne295`, et dans `output/pdf/` les séries `stanne_ICL_*`,
`stanne_surface_*` et `stanne_210us` à `stanne_285us`. **Le fichier `evolution.pdf` du dossier final
résume tout cela en une planche** : préférez-le, et ne descendez dans les jeux intermédiaires que
pour montrer une image à un instant précis.

Ancienne formulation, conservée parce qu'elle reste vraie : les radiales
n'occupent pas encore tous les secteurs à 192 µs, elles le font à 210 µs.

`output/fragments_stanne/` porte la planche de la dernière trame, le gif du mouvement des fragments
et la vérification à 120 µs.

## 3. Le reste

Dans `results/fig/` à plat :

- `coupes_s1_v3_*.pdf`, `kinetics_s1_v3*.pdf` — le run s = 1 de la nuit du 12 au 13 ;
- `yang_s*_intrinseque_*.pdf`, `yang_s*_adaptif_moyenne_*.pdf`, `yang_s*_adaptif_max_*.pdf` — la
  comparaison des trois critères d'insertion, celle qui a établi que `facetAverage = max` est le
  critère des impacts ;
- `stress_stanne_{roche,lisse,dessus}.pdf` — les premières cartes de contrainte, avant la
  systématisation en jeux datés.

## 4. Refaire une figure

Les commandes exactes sont dans `REPRODUIRE_stanne_radiales_2026-09-14.md` §4. Toutes prennent le
dossier de sortie en premier argument et un `--stem` pour le nom de sortie ; par défaut elles lisent
la **dernière trame complète**. Aucune ne demande de recompiler.

Trois d'entre elles ne lisent que `history.csv`, donc elles se refont depuis les seules données
versionnées, sans les trames :

```
python tools/fig_fp_direct.py    <dossier> --stem results/fig/x_fp
python tools/fig_kinetics.py     <dossier> --stem results/fig/x_cin
python tools/fig_retournement.py <dossier> --stem results/fig/x_ret
```

Il suffit de recréer un dossier contenant les fichiers de `results/data/stanne2025_rock137/`.
Les autres, celles de **géométrie**, ont besoin des trames VTU, qui pèsent 2 Go, ne sont pas dans
les dépôts et restent sur le poste de calcul. Il faut rejouer le run pour les refaire.

## 5. Les chiffres qui vont avec

Les mesures ne sont pas dans les figures, elles sont dans
`COMPARAISON_yang2025_stanne_2026-09-14.md`, qui confronte le run aux sept critères publiés, et dans
`ECARTS_guo2014_rockim_2026-09-13.md` §10. Les données brutes sont dans
`results/data/stanne2025_rock137/`, et chaque planche a ses `.csv` et `.json` posés à côté.
Citer ceux-là plutôt que de relire une valeur sur une image.
