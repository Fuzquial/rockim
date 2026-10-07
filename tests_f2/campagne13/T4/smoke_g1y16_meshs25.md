# Fumee `g1y16_meshs25` — exe `rockim_g1y16.exe`, T = 2e-06 s, frames = 1, OMP_NUM_THREADS = 4, meshFile surcharge -> `meshes/impact_yang_s2.5_pose.msh`

Cout estime du deck complet = pas(T_deck/dt) x ms/pas ; 77.73 ms/pas partage (results/yang_bench_s25_v3P.log : 5 498 s / 70 733 pas, 14 fils) ; 17.0 ms/pas seul (estimation du cadrage, non mesuree).

| deck | maillage | demarre | code | dt (ns) | pas (2 us) | T deck (us) | pas total | cout partage | cout seul | mur fumee | avert. | residu B4 (J) | masses piston / bit / insert (kg) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `stanne2025_bench_s25_visc0.cfg` | `impact_yang_s2.5_pose.msh` | oui | 0 | 4.165 | 481 | 450 | 108052 | 8399 s (2.3 h) | 1837 s (31 min) | 12 s | 1 | -7.97383e-14 | 0.7767 / 1.08 / 0.06313 |
| `stanne2025_bench_s25_visc.cfg` | `impact_yang_s2.5_pose.msh` | oui | 0 | 4.165 | 481 | 450 | 108052 | 8399 s (2.3 h) | 1837 s (31 min) | 10 s | 2 | 2.03237e-13 | 0.7767 / 1.08 / 0.06313 |
