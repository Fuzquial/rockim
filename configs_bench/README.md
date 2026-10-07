# configs_bench/ — banc de pénalité du 11/09/2026

Archive d'une mesure ponctuelle : de combien la pénalité intrinsèque des joints assouplit-elle le
module apparent ? Barreau de 20 × 20 × 40 mm, 5 534 tétraèdres, joints rendus incassables
(`ft = c = 1e12`) pour isoler la seule complaisance élastique des joints, traction lente
(`mode = fdem3d`, `scenario = tension`). Les cinq decks ne diffèrent que par
`jointPenaltyFactor`.

Différence avec les deux autres dossiers : [`configs/`](../configs/README.md) rassemble tous les
decks du projet, [`configs_yan/`](../configs_yan/README.md) la campagne Yan 2023 et la calibration
2D. Ce dossier-ci contient aussi les journaux des runs, ce qui en fait une archive de résultat.

| Deck | Facteur de pénalité | Journal | Contrainte macroscopique maximale | Temps de calcul |
|---|---:|---|---:|---:|
| `pen_f5.cfg` | 5 | `log_f5.txt` | 7,49 MPa | 89 s |
| `pen_f10.cfg` | 10 | `log_f10.txt` | 8,22 MPa | 123 s |
| `pen_f20.cfg` | 20 | `log_f20.txt` | 8,60 MPa | 232 s |
| `pen_f40.cfg` | 40 | `log_f40.txt` | 8,81 MPa | 233 s |
| `pen_f80.cfg` | 80 | `log_f80.txt` | 8,92 MPa | 339 s |

À déplacement imposé égal, la contrainte atteinte croît avec la pénalité et sature : l'écart entre
f = 5 et f = 80 est de 16 %, entre f = 40 et f = 80 de 1,2 %. La mention `[FAIL]` des journaux
est attendue : le contrôle automatique compare le pic à `ft = 1e12`, que des joints incassables
n'atteignent jamais.

Le maillage `meshes/bench_penalite.msh` n'est pas versionné (règle `*.msh` du `.gitignore`).
Aucune figure n'a été tirée de ce banc. L'étude de pénalité du rapport-guide
(`fig_penalite`, section sensibilités) porte sur un autre banc, la célérité d'onde dans une barre.
