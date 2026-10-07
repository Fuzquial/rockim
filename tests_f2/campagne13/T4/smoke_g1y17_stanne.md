# Fumee `g1y17_stanne` — exe `rockim_g1y17.exe`, T = 3e-06 s, frames = 1, OMP_NUM_THREADS = 4

Cout estime du deck complet = pas(T_deck/dt) x ms/pas ; 77.73 ms/pas partage (results/yang_bench_s25_v3P.log : 5 498 s / 70 733 pas, 14 fils) ; 17.0 ms/pas seul (estimation du cadrage, non mesuree).

| deck | maillage | demarre | code | dt (ns) | pas (3 us) | T deck (us) | pas total | cout partage | cout seul | mur fumee | avert. | residu B4 (J) | masses piston / bit / insert (kg) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `stanne2025_bench_s25_visc0.cfg` | `impact_yang_train1_rock25_hxt.msh` | oui | 0 | 4.027 | 745 | 450 | 111735 | 8685 s (2.4 h) | 1899 s (32 min) | 93 s | 1 | -4.24238e-12 | 1.058 / 1.288 / 0.06458 |
| `stanne2025_bench_s25_visc.cfg` | `impact_yang_train1_rock25_hxt.msh` | oui | 0 | 4.027 | 745 | 450 | 111735 | 8685 s (2.4 h) | 1899 s (32 min) | 100 s | 2 | -3.96394e-12 | 1.058 / 1.288 / 0.06458 |
