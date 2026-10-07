# Fumee `g1y16_refus` — exe `rockim_g1y16.exe`, T = 2e-06 s, frames = 1, OMP_NUM_THREADS = 4

Cout estime du deck complet = pas(T_deck/dt) x ms/pas ; 77.73 ms/pas partage (results/yang_bench_s25_v3P.log : 5 498 s / 70 733 pas, 14 fils) ; 17.0 ms/pas seul (estimation du cadrage, non mesuree).

| deck | maillage | demarre | code | dt (ns) | pas (2 us) | T deck (us) | pas total | cout partage | cout seul | mur fumee | avert. | residu B4 (J) | masses piston / bit / insert (kg) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `yang2026_bench_s25_v4_A.cfg` | `impact_yang_s2.5_pose.msh` | NON | 1 | 2.828 | 708 | 300 | 106098 | 8247 s (2.3 h) | 1804 s (30 min) | 0 s | 1 | - | - / - / - |
| `yang2026_bench_s25_v4_B.cfg` | `impact_yang_s2.5_pose.msh` | NON | 1 | 4.079 | 491 | 300 | 73556 | 5718 s (1.6 h) | 1250 s (21 min) | 0 s | 1 | - | - / - / - |
| `yang2026_bench_s25_v4_B1.cfg` | `impact_yang_s2.5_pose.msh` | NON | 1 | 2.828 | 708 | 300 | 106098 | 8247 s (2.3 h) | 1804 s (30 min) | 0 s | 1 | - | - / - / - |
| `yang2026_bench_s25_v4_B2.cfg` | `impact_yang_s2.5_pose.msh` | NON | 1 | 4.079 | 491 | 300 | 73556 | 5718 s (1.6 h) | 1250 s (21 min) | 0 s | 1 | - | - / - / - |
| `yang2026_bench_s25_v4_C.cfg` | `impact_yang_s2.5_pose.msh` | NON | 1 | 4.079 | 491 | 300 | 73556 | 5718 s (1.6 h) | 1250 s (21 min) | 0 s | 2 | - | - / - / - |
| `yang2026_bench_s25_v4_D.cfg` | `impact_yang_s2.5_pose.msh` | NON | 1 | 2.828 | 708 | 300 | 106098 | 8247 s (2.3 h) | 1804 s (30 min) | 0 s | 1 | - | - / - / - |

Refus / erreurs de `yang2026_bench_s25_v4_A.cfg` :

    [rockim] error: 1 cle du deck 'C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1\tests_f2\campagne13\T4\yang2026_bench_s25_v4_A_g1y16_refus\yang2026_bench_s25_v4_A_smoke.cfg' refusee (unknownKeys = error) ; poser unknownKeys = warn pour continuer avec un avertissement :
    - cle 'jointBreakModeRef' (ligne 352 du deck) inconnue de rockim (aucun getter du code ne la lit, aucune cle connue a moins de 2 caracteres)

Refus / erreurs de `yang2026_bench_s25_v4_B.cfg` :

    [rockim] error: 1 cle du deck 'C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1\tests_f2\campagne13\T4\yang2026_bench_s25_v4_B_g1y16_refus\yang2026_bench_s25_v4_B_smoke.cfg' refusee (unknownKeys = error) ; poser unknownKeys = warn pour continuer avec un avertissement :
    - cle 'jointBreakModeRef' (ligne 381 du deck) inconnue de rockim (aucun getter du code ne la lit, aucune cle connue a moins de 2 caracteres)

Refus / erreurs de `yang2026_bench_s25_v4_B1.cfg` :

    [rockim] error: 1 cle du deck 'C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1\tests_f2\campagne13\T4\yang2026_bench_s25_v4_B1_g1y16_refus\yang2026_bench_s25_v4_B1_smoke.cfg' refusee (unknownKeys = error) ; poser unknownKeys = warn pour continuer avec un avertissement :
    - cle 'jointBreakModeRef' (ligne 353 du deck) inconnue de rockim (aucun getter du code ne la lit, aucune cle connue a moins de 2 caracteres)

Refus / erreurs de `yang2026_bench_s25_v4_B2.cfg` :

    [rockim] error: 1 cle du deck 'C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1\tests_f2\campagne13\T4\yang2026_bench_s25_v4_B2_g1y16_refus\yang2026_bench_s25_v4_B2_smoke.cfg' refusee (unknownKeys = error) ; poser unknownKeys = warn pour continuer avec un avertissement :
    - cle 'jointBreakModeRef' (ligne 352 du deck) inconnue de rockim (aucun getter du code ne la lit, aucune cle connue a moins de 2 caracteres)

Refus / erreurs de `yang2026_bench_s25_v4_C.cfg` :

    [rockim] error: 1 cle du deck 'C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1\tests_f2\campagne13\T4\yang2026_bench_s25_v4_C_g1y16_refus\yang2026_bench_s25_v4_C_smoke.cfg' refusee (unknownKeys = error) ; poser unknownKeys = warn pour continuer avec un avertissement :
    - cle 'jointBreakModeRef' (ligne 373 du deck) inconnue de rockim (aucun getter du code ne la lit, aucune cle connue a moins de 2 caracteres)

Refus / erreurs de `yang2026_bench_s25_v4_D.cfg` :

    [rockim] error: 1 cle du deck 'C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1\tests_f2\campagne13\T4\yang2026_bench_s25_v4_D_g1y16_refus\yang2026_bench_s25_v4_D_smoke.cfg' refusee (unknownKeys = error) ; poser unknownKeys = warn pour continuer avec un avertissement :
    - cle 'jointBreakModeRef' (ligne 359 du deck) inconnue de rockim (aucun getter du code ne la lit, aucune cle connue a moins de 2 caracteres)
