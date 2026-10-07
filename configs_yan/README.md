# configs_yan/ — reproduction de Yan et al. (2023) et calibration 2D Red Bohus

Decks 2D (`mode = fdem`) de la campagne qui reproduit Yan, Zheng et Wang (IJRMMS 169, 2023) et
Yan, Wang et Jiao (IJRMMS 169, 2023, 105439), puis s'en sert pour calibrer le granite de Red Bohus.
95 configurations et 4 générateurs Python. Les ajouts au code faits pour cette campagne sont
décrits dans [`docs/notes/CHANGES_YAN.md`](../docs/notes/CHANGES_YAN.md) ; les résultats dans la
section « Insertion adaptative : reproduction de Yan et al. 2023 » du
[rapport-guide](../docs/rapport_guide/rapport_guide_rockim.pdf).

Différence avec les deux autres dossiers : [`configs/`](../configs/README.md) rassemble tous les
decks du projet ; [`configs_bench/`](../configs_bench/README.md) archive un seul banc de pénalité.
Ici, chaque essai est décliné en séries qui ne changent qu'une clé.

## Convention des suffixes

| Suffixe | Sens |
|---|---|
| `_adap` | insertion adaptative des joints (schéma de Yan) |
| `_pN`, `_iN` | insertion intrinsèque, facteur de pénalité `jointPenaltyFactor = N` |
| `mesh_aXXX`, `mesh_iXXX` | brésilien, taille de maille XXX µm, adaptatif (`a`) ou intrinsèque (`i`) |
| `_0`, `_10`, `_20`, `_40` | confinement en MPa |
| `_s4211`, `_s7`, `_s99` | graine du tirage du maillage |

## Inventaire

| Série | Fichiers | Nombre | Objet | Statut |
|---|---|---:|---|---|
| Préréglage | `article_exact_base.cfg` | 1 | toutes les formes littérales de Yan 2023, à inclure en fin de deck | référence |
| Brésilien | `bd_yan_calibre`, `bd_adap`, `bd_p1` à `bd_p500`, `bd_mesh_*`, `mesh_*` | 19 | section 3.2, fig. 11-14 : balayage de la pénalité et de la taille de maille | étude |
| Compression simple | `ucs_adap`, `ucs_p1` à `ucs_p500` | 10 | fig. 17, éprouvette 36 × 72 mm | `ucs_adap` dans la suite de tests ; étude |
| Balayage de la loi cohésive | `ucs_sw_*` (généré par `gen_ucs_sweep.py`) | 8 | sensibilité de l'UCS homogène à ft, c, GfI, GfII | étude |
| Triaxial | `tx_adap_*`, `tx_p50_*`, `tx_p100_*`, `tx_gbm3_*` (généré par `gen_tx_gbm3.py`) | 15 | confinements 0 à 40 MPa, milieu homogène ou GBM à trois phases | étude |
| Calibration Red Bohus | `cal20_*` (généré par `gen_calib20.py`) | 8 | triaxial à σ3 = 20 MPa, lot de dépistage | étude |
| Effet de taille | `size_L0` à `size_L4` (généré par `gen_size_study.py`) | 11 | peut-on calibrer sur une petite éprouvette ? trois graines | étude |
| Bande pesante | `strip_adap`, `strip_i1` à `strip_i500` | 10 | section 3.1 : vérification en milieu continu | étude |
| SHPB | `shpb_mini`, `shpb_barre`, `shpb_barre_fac2`, `shpb_complet_adaptatif`, `shpb_p1` à `shpb_p500` | 13 | section 3.4, fig. 22-25 : barres de Hopkinson | `shpb_mini` dans la suite de tests ; étude |

Les deux decks lus par `tools/verify_suite.py` (niveau `full`) sont `ucs_adap.cfg` et
`shpb_mini.cfg`. Les fichiers `results/*.err` portent le nom des decks `cal20_*`, `tx_gbm3_*` et
`ucs_sw_*` : ce sont les sorties d'erreur, vides, de leurs runs.
