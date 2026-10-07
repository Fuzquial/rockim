# bench_phases : banc court du matériau par phase en éléments finis 3D

## 1. Objet

Banc de vérification de la capacité « propriétés par phase minérale » du solveur `mode = fem3d` (feldspath, quartz, biotite du Red Bohus). Huit configurations tournent en moins d'une seconde chacune sur quelques centaines à quelques milliers de tétraèdres.

Question : l'indice de phase atteint-il réellement la loi de comportement, la masse, le pas de temps stable et les bilans, et toute clé de phase mal posée est-elle refusée avec un message juste ?

## 2. Statut

Référence, écrite le 6 septembre 2026 avec la capacité elle-même ; le banc doit être vert avant toute campagne à phases. Les calculs sont élastiques, sans joints ni contact entre fragments : la correction du contact du 7 octobre 2026 ne les concerne pas, et la calibration de juillet non plus. Le binaire de l'époque (`rockim_g0`, Windows) n'a pas été rejoué sous Linux pour ce banc.

## 3. Résultat principal

Les neuf contrôles passent (`VERDICT_depouille.txt`) :

- F2 et F5 : des phases identiques à la fiche globale rendent `history.csv` identique octet pour octet, sur le chemin fichier et sur le chemin Voronoï ; la variante à contraste diffère bien.
- F3 : barre de deux couches en série, module apparent 44,158 GPa contre 43,678 GPa par la borne de Reuss (1,10 %) ; rapport série / couche dure 0,5253, qui vaudrait 1 si la loi n'était pas indexée par la phase.
- F3d : couches permutées, −1,48 % (asymétrie des mors). F4 : pas de temps divisé par 0,8497, valeur attendue √(60/83,1).
- Les 16 cas fautifs sont refusés avec le bon message (`VERDICT_erreurs.txt`) ; 8 configurations de référence restent identiques au bit près (`VERDICT_bitid.txt`).

## 4. Figures

Aucune figure : les sorties du banc sont les verdicts texte ci-dessus. La vérification correspondante est résumée dans le rapport-guide (chapitre de calibration, sous-section « Modèle à grains »).

## 5. Contenu du dossier

| chemin | contenu |
|---|---|
| `LISEZ_MOI.md` | description détaillée des contrôles, de ce qu'ils falsifient et de ce que le banc ne prouve pas |
| `make_bar2.py` | barre à deux couches nommées (`bar2.msh`) ou sans noms (`bar1.msh`) |
| `f2_*.cfg`, `f3a_*.cfg` à `f3d_*.cfg`, `f5a_*.cfg` à `f5c_*.cfg` | les huit configurations (neutralité, Reuss, commutativité, Voronoï) |
| `log_*.txt` | journaux des huit calculs |
| `depouille.py`, `VERDICT_depouille.txt` | neuf contrôles et leur verdict |
| `erreurs.py`, `VERDICT_erreurs.txt` | seize cas qui doivent échouer, et leur verdict |
| `bitid_phases_fem3d.json`, `VERDICT_bitid.txt` | non-régression bit à bit de huit configurations de référence |
| `suite_avant_HEAD.txt`, `suite_apres_phases.txt` | suite de vérification avant et après la capacité |

## 6. Rejouer

Depuis la racine du dépôt :

```sh
export OMP_NUM_THREADS=1
python3 bench_phases/make_bar2.py bench_phases/bar2.msh
python3 bench_phases/make_bar2.py bench_phases/bar1.msh --sans-noms
for d in f3a_homogene_dur f3b_homogene_mou f2_deux_phases_identiques f3c_serie_reuss \
         f3d_serie_permutee f5a_voronoi_1phase f5b_voronoi_3phases_egales f5c_voronoi_3phases_contrastees; do
  build/rockim bench_phases/$d.cfg bench_phases/out_$d > bench_phases/log_$d.txt 2>&1
done
python3 bench_phases/depouille.py bench_phases   # 9 contrôles
python3 bench_phases/erreurs.py build/rockim      # 16 cas fautifs (défaut : build/rockim.exe)
```

La version PowerShell d'origine est dans `LISEZ_MOI.md`. Chapitre du rapport-guide : `docs/rapport_guide/sections/r07a_calib.tex` (sous-section « Vérifications » du modèle à grains), dans [rapport_guide_rockim.pdf](../docs/rapport_guide/rapport_guide_rockim.pdf).

Limites rappelées par `LISEZ_MOI.md` : ce n'est pas un modèle à grains cohésif (nœuds partagés, aucune frontière de grain ne s'ouvre) ; les tétraèdres linéaires sans B-bar sur-raidissent avec ν, ce qui sous-estime le contraste ; F3 est élastique et ne dit rien de l'objectivité de maillage.
