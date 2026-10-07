# gbm_uniaxial : compression élastique d'une microstructure à grains contre un témoin homogène

## 1. Objet

Essai uniaxial élastique en éléments finis 3D (`mode = fem3d`, aucun joint) sur une éprouvette de 16 × 16 × 32 mm découpée en 1 001 grains de Voronoï, 40 724 tétraèdres. Même maillage, même chargement : seule change la répartition du matériau, trois phases minérales du Red Bohus (feldspath 62 %, quartz 31 %, biotite 7 %) contre un matériau unique aux moyennes de Voigt.

Question : le seul contraste de raideur entre grains concentre-t-il les contraintes, et les maxima se placent-ils près des frontières de grain ?

## 2. Statut

Étude terminée le 6 septembre 2026, première application de la capacité « matériau par phase » vérifiée par `bench_phases/`. Calcul élastique, sans joints ni contact entre fragments : la correction du contact du 7 octobre 2026 et la non-reproductibilité de la calibration de juillet ne le concernent pas. Deux réserves écrites dans `UNIAXIAL.md` restent ouvertes : l'éprouvette est un parallélépipède et non un cylindre, et le maillage intragranulaire est l'éventail par défaut de la tessellation interne ; `gbm_neper_like.py` a été écrit pour les lever, sans nouveau calcul archivé.

## 3. Résultat principal

- Module apparent : 69,24 GPa pour le modèle à grains, 71,32 GPa pour le témoin (Voigt 71,21, Reuss 66,77) ; le modèle à grains est 2,91 % plus souple (`UNIAXIAL.md` §1).
- Concentration : le quantile 99 % de la contrainte équivalente vaut 1,33 fois la médiane contre 1,04 pour le témoin ; le maximum est 1,53 fois celui du témoin (`UNIAXIAL.md` §2).
- Par phase, rapport à la contrainte du témoin : quartz 1,16, feldspath 0,97, biotite 0,50, dans l'ordre des modules (`UNIAXIAL.md` §3).
- Les maxima ne se rapprochent presque pas des frontières de grain : −7 % de distance au joint pour le 0,1 % le plus chargé, ±1 % au-delà (`UNIAXIAL.md` §4). Les écarts sont des bornes basses : les tétraèdres linéaires sur-raidissent davantage les phases à ν élevé (0,17 à 0,36).

## 4. Figures

Aucune figure versionnée : `depouille_uniaxial.py`, `fig_maillage.py` et `fig_mesh.py` écrivent dans un dossier `figures/` absent du dépôt, à partir des sorties VTU (`fem3d_0021.vtu`), non versionnées.

## 5. Contenu du dossier

| chemin | contenu |
|---|---|
| `UNIAXIAL.md` | note de résultats : module apparent, concentration, répartition par phase, distance des maxima aux joints de grain |
| `U_gbm.cfg`, `U_homo.cfg` | les deux essais (trois phases ; témoin aux moyennes de Voigt) |
| `log_gbm.txt`, `log_homo.txt` | journaux des deux calculs |
| `depouille_uniaxial.py` | dépouillement élément par élément des deux calculs, écrit `UNIAXIAL.md` et la planche de synthèse |
| `fig_maillage.py`, `fig_mesh.py` | figures du maillage à grains |
| `gbm_neper_like.py` | générateur de microstructure cylindrique (volumes physiques nommés par phase) qui remplace la tessellation interne |

## 6. Rejouer

Depuis la racine du dépôt :

```sh
build/rockim gbm_uniaxial/U_gbm.cfg out_U_gbm
build/rockim gbm_uniaxial/U_homo.cfg out_U_homo
python3 gbm_uniaxial/depouille_uniaxial.py
```

Le script lit `out_U_gbm/` et `out_U_homo/` dans le répertoire courant et y écrit `UNIAXIAL.md` et `figures/uniaxial_synthese.*` : lancé depuis la racine, il ne touche pas la note versionnée de ce dossier. Chapitre du rapport-guide : `docs/rapport_guide/sections/r07a_calib.tex`, sous-section « Vérifications » du modèle à grains, dans [rapport_guide_rockim.pdf](../docs/rapport_guide/rapport_guide_rockim.pdf).
