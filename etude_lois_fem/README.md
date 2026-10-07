# etude_lois_fem : briques constitutives en éléments finis continus (percussion et coupe sous confinement)

## 1. Objet

Étude des lois de comportement continues du solveur `mode = fem3d`, brique par brique : méridien courbe (loi puissance), coupure de traction de Rankine, cap de compaction, endommagement compressif, suppression d'éléments. Les briques sont comparées une clé à la fois dans le noyau `dpr` / `saksala` de rockim, à carte Red Bohus égale et sans calibration, sur deux essais : la percussion 3D d'un insert (R = 7,94 mm, 16 J) et la coupe en tranche plane par une lame rigide, sous 0, 50 et 100 MPa de confinement.

Question : quelle brique change les observables convergentes d'un impact et d'une coupe, et à partir de quelle pression chacune devient inerte (moins de 10 %) ou indispensable (plus de 20 %) ?

## 2. Statut

| date | état |
|---|---|
| 3 septembre 2026 | proposition (`PROPOSITION_etude_lois_fem.md`), phase A (9 calculs), codage des briques par addition |
| 4-5 septembre 2026 | matrice de percussion v1 puis v2 (27 calculs) ; dépouillement et verdicts |
| 5-6 septembre 2026 | coupe en file, non jugée ; ajout du matériau par phase (voir `bench_phases/`) |

Étude en cours, interrompue après la percussion ; la coupe n'a pas été dépouillée. Calculs en éléments finis continus, sans joints ni contact entre fragments : la correction du contact du 7 octobre 2026 et la non-reproductibilité de la calibration de juillet ne les concernent pas. Les sorties, les verdicts détaillés (`matrice_v2/VERDICTS_v2.md`, `phaseA/RESULTATS_phaseA.md`) et l'enquête préalable (`enquete/`) sont restés sur le poste du doctorant ; ce dossier n'en garde que le journal, la proposition et les scripts.

## 3. Résultat principal

Percussion, matrice v2 (`JOURNAL.md`, entrées du 5 septembre 18:25 et 19:00) :

- À P = 0, seule la coupure de traction compte : passer à l'apex du cône (23 MPa) donne +22 % de force de pic, un coefficient de restitution doublé et un volume très endommagé divisé par 2,7 ; méridien, cap et endommagement compressif sont inertes à 3 %.
- Dès 50 MPa, la traction s'éteint (volume à D > 0,9 de 11 863 à 140 mm³) et le méridien devient la seule brique active : le méridien linéaire donne +13 % de pic à 50 MPa, +15 % à 100 MPa, −26 % d'énergie dissipée.
- Le cap de compaction est de second ordre (−5 % de pic, +8 % d'énergie à 100 MPa) ; l'endommagement compressif est inerte en percussion (d_c ≤ 0,29), ce qui falsifie la prédiction « dominant à 100 MPa ».
- Deux défauts du noyau d'origine ont été trouvés et corrigés par des clés opt-in : retour au-delà de l'apex du cône (`dpApex`) et coupure de Rankine pilotée par la déformation (`rankineDrive`) (`doc_section_5_18.md`).

## 4. Figures

Aucune figure versionnée : `fig_apercu.py`, `fig_matrice.py` et `fig_FP.py` écrivent dans un dossier de sortie du poste du doctorant (`fig1` à `fig7`), à partir de sorties non versionnées.

## 5. Contenu du dossier

| chemin | contenu |
|---|---|
| `PROPOSITION_etude_lois_fem.md` | proposition d'étude du 3 septembre : thèse, briques, essais, coût, décisions, table triaxiale analytique |
| `JOURNAL.md` | journal du chantier : binaires, chronologie, phases A et B, matrices v1 et v2, verdicts, chantier du matériau par phase |
| `doc_section_5_18.md` | section 5.18 de la documentation : clés des briques et du solveur fem3d ajoutées |
| `make_cfgs.py` | génération des configurations (liste blanche de clés) |
| `extract_fem3d.py` | dépouillement des calculs fem3d |
| `fig_apercu.py`, `fig_matrice.py`, `fig_FP.py` | figures de la phase A, de la matrice et de la force-pénétration |
| `meshes/T1_c05_clean.msh` | maillage de percussion sans nœud orphelin |
| `meshes/check_orphans.py`, `drop_orphans.py` | détection et retrait des nœuds orphelins (défaut de la « broche fantôme » du 5 septembre) |

## 6. Rejouer

Le contrôle analytique du méridien tourne sans maillage :

```sh
build/rockim selftest-triax triax.csv
```

Les configurations de la matrice sont produites par `make_cfgs.py` ; le générateur de maillages et `matrice/make_matrice.py` cités dans le journal ne sont pas dans le dépôt. Les clés sont documentées dans `doc_section_5_18.md` et dans `DOCUMENTATION_rockim.md` §5.18. L'étude n'a pas de chapitre propre dans le rapport-guide ; la comparaison des lois continues sous Abaqus figure dans `docs/rapport_guide/sections/r07a_calib.tex` (« Septembre et calibrations continues »).
