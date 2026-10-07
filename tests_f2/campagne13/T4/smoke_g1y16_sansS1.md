# Fumee `g1y16_sansS1` — exe `rockim_g1y16.exe`, T = 2e-06 s, frames = 1, OMP_NUM_THREADS = 4, cle retiree `jointBreakModeRef`

Cout estime du deck complet = pas(T_deck/dt) x ms/pas ; 77.73 ms/pas partage (results/yang_bench_s25_v3P.log : 5 498 s / 70 733 pas, 14 fils) ; 17.0 ms/pas seul (estimation du cadrage, non mesuree).

| deck | maillage | demarre | code | dt (ns) | pas (2 us) | T deck (us) | pas total | cout partage | cout seul | mur fumee | avert. | residu B4 (J) | masses piston / bit / insert (kg) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `yang2026_bench_s25_v4_A.cfg` | `impact_yang_s2.5_pose.msh` | oui | 0 | 2.828 | 708 | 300 | 106098 | 8247 s (2.3 h) | 1804 s (30 min) | 16 s | 1 | 5.7454e-15 | 0.7767 / 1.08 / 0.06313 |
| `yang2026_bench_s25_v4_B.cfg` | `impact_yang_s2.5_pose.msh` | oui | 0 | 4.079 | 491 | 300 | 73556 | 5718 s (1.6 h) | 1250 s (21 min) | 9 s | 1 | 5.76865e-14 | 0.7767 / 1.08 / 0.06313 |
| `yang2026_bench_s25_v4_B1.cfg` | `impact_yang_s2.5_pose.msh` | oui | 0 | 2.828 | 708 | 300 | 106098 | 8247 s (2.3 h) | 1804 s (30 min) | 15 s | 1 | 6.93889e-15 | 0.7767 / 1.08 / 0.06313 |
| `yang2026_bench_s25_v4_B2.cfg` | `impact_yang_s2.5_pose.msh` | oui | 0 | 4.079 | 491 | 300 | 73556 | 5718 s (1.6 h) | 1250 s (21 min) | 34 s | 1 | 3.91354e-15 | 0.7767 / 1.08 / 0.06313 |
| `yang2026_bench_s25_v4_C.cfg` | `impact_yang_s2.5_pose.msh` | oui | 0 | 4.079 | 491 | 300 | 73556 | 5718 s (1.6 h) | 1250 s (21 min) | 8 s | 2 | -3.93158e-14 | 0.7767 / 1.08 / 0.06313 |
| `yang2026_bench_s25_v4_D.cfg` | `impact_yang_s2.5_pose.msh` | oui | 0 | 2.828 | 708 | 300 | 106098 | 8247 s (2.3 h) | 1804 s (30 min) | 14 s | 1 | -6.70852e-14 | 0.7767 / 1.08 / 0.06313 |
