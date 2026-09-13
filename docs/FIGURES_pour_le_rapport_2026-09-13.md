# Où sont les figures, et ce que chacune montre — 13/09/2026

257 figures sont versionnées, en PNG de relecture et en **PDF vectoriel Computer Modern**
directement inclusible en LaTeX. Elles vivent aux deux mêmes endroits dans les deux dépôts :

| | Dépôt `Fuzquial/rockim`, branche `g1` | Archive `phd_geothermie` |
|---|---|---|
| chemin | `results/fig/`, `output/pdf/` | `FDEM/rockim_g1/results/fig/`, `.../output/pdf/` |

Prendre toujours le `.pdf`, jamais le `.png`, pour un document LaTeX.

## 1. L'état le plus avancé du run St Anne 2025 : `results/fig/stanne210/`

Trame 14, t = 210 µs, maillage rock137 (108 667 tétraèdres, cœur à 1,42 mm).
C'est le jeu à citer par défaut.

| Fichier | Ce qu'il montre |
|---|---|
| `joints.pdf` | les 8 777 facettes rompues seules, six vues, couleur par mode de rupture |
| `coupes.pdf` | coupes exactes en y = 0, 5, 10 mm et z = −0,5, −3, −8 mm ; la trace est l'intersection facette/plan, pas une projection |
| `sig_coupe.pdf` | σ₁ en coupe verticale, lissée, roche seule — la chute de contrainte sur les lèvres des radiales s'y lit |
| `sig_dessus.pdf` | σ₁ en vue de dessus à z = −2 mm, non lissée |
| `cratere.png` | vue de dessus sectorisée en 16, avec la portée des radiales secteur par secteur |
| `branches.pdf` | branches connectées, classées radiale / conique / horizontale ; chiffres dans les deux `.csv` à côté |
| `fp.pdf` | force-pénétration, courbe mesurée au contact **et** courbe estimée par le train rigide |
| `cinetique.pdf` | vitesses d'indentation et de rebond à la manière de Yang (pente des portions linéaires), plus la jauge |

Deux compagnons dans `output/pdf/stanne_210us/` : `ICL/` (revue en trois planches, avec
`mesures_et_provenance.json` qui porte les empreintes SHA-256 des trames lues) et `surface/`
(comparaison des traces en surface et sous la surface, avec ses mesures).

## 2. Les états antérieurs du même run

`results/fig/stanne192/` et, dans `output/pdf/`, les séries `stanne_ICL_165us`, `stanne_ICL_180us`,
`stanne_ICL_192us` et leurs `stanne_surface_*`. Utiles pour une planche d'évolution : les radiales
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
  systématisation en `stanne210/`.

## 4. Refaire une figure

Les commandes exactes sont dans `REPRODUIRE_stanne_radiales_2026-09-14.md` §4. Toutes prennent le
dossier de sortie en premier argument et un `--stem` pour le nom de sortie ; par défaut elles lisent
la **dernière trame complète**. Elles n'ont besoin ni de recompiler ni de relancer le calcul, mais
elles ont besoin des trames VTU, qui ne sont **pas** dans les dépôts et restent sur le poste de
calcul.

## 5. Les chiffres qui vont avec

Les mesures ne sont pas dans les figures, elles sont dans
`ECARTS_guo2014_rockim_2026-09-13.md` §10 pour St Anne, et dans les `.csv` et `.json` posés à côté de
chaque planche. Citer ceux-là plutôt que de relire une valeur sur une image.
