# Fumee `masses` — exe `rockim_g1y16.exe`, T = 2e-06 s, frames = 1, OMP_NUM_THREADS = 4

Cout estime du deck complet = pas(T_deck/dt) x ms/pas ; 77.73 ms/pas partage (results/yang_bench_s25_v3P.log : 5 498 s / 70 733 pas, 14 fils) ; 17.0 ms/pas seul (estimation du cadrage, non mesuree).

| deck | maillage | demarre | code | dt (ns) | pas (2 us) | T deck (us) | pas total | cout partage | cout seul | mur fumee | avert. | residu B4 (J) | masses piston / bit / insert (kg) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `t3_masses_temoin.cfg` | `impact_yang_train1_rock25_hxt.msh` | oui | 0 | 2.983 | 671 | 200 | 67043 | 5211 s (1.4 h) | 1140 s (19 min) | 75 s | 1 | -4.07059e-12 | 1.058 / 1.288 / 0.06458 |
| `t3_masses_rho.cfg` | `impact_yang_train1_rock25_hxt.msh` | oui | 0 | 2.983 | 671 | 200 | 67043 | 5211 s (1.4 h) | 1140 s (19 min) | 72 s | 1 | -3.71664e-12 | 1.173 / 1.43 / 0.06458 |
| `t3_masses_rhoE.cfg` | `impact_yang_train1_rock25_hxt.msh` | oui | 0 | 2.983 | 671 | 200 | 67043 | 5211 s (1.4 h) | 1140 s (19 min) | 61 s | 1 | -3.81611e-12 | 1.173 / 1.43 / 0.06458 |
