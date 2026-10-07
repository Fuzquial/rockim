# Balayage statique des decks — regle C1 des cles (w21, lecteurs w22) — 2026-09-05 20:50:38

Script : `tools/scan_decks.py` (analyse STATIQUE : registre `tools/keys_by_mode.json` du 2026-09-05 20:44:31, 389 cles connues avec leurs lecteurs, 150 propres a un seul mode ; obsoletes `tools/obsolete_keys.json` : `fragBrushV` -> `fragBrushV0`). Racines : `C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_f2`, `C:/Users/fuzquianoalricabi/simulations/CONTINUUM/calib_bohus_triax/cdp_rockim`. Dossiers `out_*`, `orig/` (copies avant edition) et `selftest_*` (bancs falsifiants, fautes voulues) ignores. Aucun run de rockim : les lectures conditionnelles ne sont pas jugeables ici (regle w22 : cle du registre dont le mode courant ou le code partage est lecteur = legitime ; cle qu'aucun lecteur du mode ne lit = « sans effet en mode »). Les gardes PRE-INIT de `main.cpp` (thermal/bedding/law/phases/mesh/hydro hors mode) sont modelisees ; les gardes de MAILLAGE (C3 : noeud orphelin, tet plat) ne le sont pas (il faudrait lire le .msh). 76 fichier(s) sans cle de solveur (point materiel `rockim matpoint`, cartes materiau, decks temporaires de selftest) comptes a part, non juges.

## Bilan

| decks de solveur lus | sans faute | avec faute(s) | cles fautives | obsoletes | autre mode (lecteurs) | gardes pre-init main.cpp | fautes de frappe | inconnues | cles dynamiques (non jugees) | fichiers sans cle de solveur (non juges) |
|---|---|---|---|---|---|---|---|---|---|---|
| 881 | 878 | 3 | 9 | 0 | 0 | 1 | 0 | 8 | 400 | 76 |

Lecture : un deck « avec faute » serait REFUSE par `rockim_f2w22.exe` (`unknownKeys = error`, defaut) — sauf si la cle fautive est en realite lue par un getter que le registre ne voit pas (cle construite) ; le verdict definitif est celui du binaire. Une cle dynamique (`phase.<nom>.E`, `groupBond.<A>.<B>`, `gauge.<groupe>`...) n'est jugee qu'au run.

Fichiers sans cle de solveur (matpoint / cartes / temporaires, non juges) : `CDP/audit_crush/biax.cfg`, `CDP/audit_crush/biax0_30.cfg`, `CDP/audit_crush/biax0_300.cfg`, `CDP/audit_crush/biax0_5000.cfg`, `CDP/audit_crush/biax0_50000.cfg`, `CDP/audit_crush/oedo.cfg`, `CDP/audit_crush/oedo30.cfg`, `CDP/audit_crush/oedo_30.cfg`, `CDP/audit_crush/oedo_300.cfg`, `CDP/audit_crush/oedo_5000.cfg`, `CDP/audit_crush/oedo_50000.cfg`, `CDP/audit_crush/tens50.cfg`, `CDP/audit_crush/tens50_e6_gfi.cfg`, `CDP/audit_crush/tens50_e6_tab001.cfg`, `CDP/audit_crush/tens50_e6_tab01.cfg`, `CDP/audit_crush/tens50_e6_tab05.cfg`, `CDP/audit_crush/tens50_res_gfi.cfg`, `CDP/audit_crush/tens50_res_tab001.cfg`, `CDP/audit_crush/tens50_res_tab01.cfg`, `CDP/audit_crush/tens50_res_tab05.cfg`, `CDP/audit_crush/triax.cfg`, `CDP/audit_crush/triax300_30.cfg`, `CDP/audit_crush/triax300_300.cfg`, `CDP/audit_crush/triax300_5000.cfg`, `CDP/audit_crush/triax300_50000.cfg`, `CDP/audit_crush/triax300_600_e30.cfg`, `CDP/audit_crush/triax50.cfg`, `CDP/audit_crush/uni_30.cfg`, `CDP/audit_crush/uni_300.cfg`, `CDP/audit_crush/uni_5000.cfg`, `CDP/audit_crush/uni_50000.cfg`, `CDP/audit_crush/uniax.cfg`, `CDP/carte_confines_Gc.cfg`, `CDP/carte_confines_Gc_matpoint.cfg`, `CDP/carte_confines_Wexc.cfg`, `CDP/carte_confines_Wexc_matpoint.cfg`, `CDP/carte_toutes_Gc.cfg`, `CDP/carte_toutes_Gc_matpoint.cfg`, `CDP/carte_toutes_Wexc.cfg`, `CDP/carte_toutes_Wexc_matpoint.cfg`, `CDP/matpoint_cdp_hist.cfg`, `CDP/matpoint_cdp_hist_compband_lc1.cfg`, `CDP/matpoint_cdp_hist_compband_lc4.cfg`, `CDP/matpoint_cdp_hydro_cap.cfg`, `CDP/matpoint_cdp_inverse.cfg`, `CDP/revue/bigsteps_ten.cfg`, `CDP/revue/psi55_biax.cfg`, `CDP/revue/tiny_lc.cfg`, `CDP/se_cdp/mp_SE_cdp_K0.cfg`, `CDP/se_cdp/mp_SE_cdp_P050.cfg`, `CDP/se_cdp/mp_SE_cdp_P300.cfg`, `CDP/selftests/cdp_w11.csv.tmp.cfg`, `CDP/selftests/cdp_w12.csv.tmp.cfg`, `CDP/selftests/selftest_cdp_w12.csv.tmp.cfg`, `CDP/selftests/selftest_triax_w12.csv.tmp.cfg`, `CDP/selftests/triax_w11.csv.tmp.cfg`, `CDP/selftests/triax_w12.csv.tmp.cfg`, `CDP/tmp/calA_confines.cfg`, `CDP/tmp/calA_inverse.cfg`, `CDP/tmp/calA_toutes.cfg`, `CDP/tmp/calB_confines_Gc.cfg`, `CDP/tmp/calB_confines_Wexc.cfg`, `CDP/tmp/calB_toutes_Gc.cfg`, `CDP/tmp/calB_toutes_Wexc.cfg`, `CDP/tmp/cmp_inverse.cfg`, `CDP/tmp/cmp_pente.cfg`, `etude_lois_fem/selftests/triax_w1a.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w1b.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w1c.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w1d.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w2.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w3.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w5.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w6.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w7.csv.tmp.cfg`, `etude_lois_fem/selftests/triax_w9.csv.tmp.cfg`

## Fautes par cle (toutes ; a decider par Fernando, sauf les obsoletes deja corrigees)

| type | cle | diagnostic | decks | ou (deck:ligne) |
|---|---|---|---|---|
| garde | `thermal` | 'thermal' n'est implemente que pour mode = fdem (garde main.cpp) | 1 | `tests_f2/tg_3d.cfg:29` |
| unknown | `dfhKPsi` | inconnue de rockim | 2 | `bench_impact/configs/impact3d_dpdfh.cfg:39`, `bench_impact/configs/impact3d_dpdfh_gros.cfg:39` |
| unknown | `dfhPsi0` | inconnue de rockim | 2 | `bench_impact/configs/impact3d_dpdfh.cfg:38`, `bench_impact/configs/impact3d_dpdfh_gros.cfg:38` |
| unknown | `dfhPsiMax` | inconnue de rockim | 2 | `bench_impact/configs/impact3d_dpdfh.cfg:40`, `bench_impact/configs/impact3d_dpdfh_gros.cfg:40` |
| unknown | `dfhPsiVar` | inconnue de rockim | 2 | `bench_impact/configs/impact3d_dpdfh.cfg:37`, `bench_impact/configs/impact3d_dpdfh_gros.cfg:37` |

## Par deck (decks avec au moins une faute : 3)

| deck | mode | fautes |
|---|---|---|
| `bench_impact/configs/impact3d_dpdfh.cfg` | fdem3d | `dfhPsiVar` (l.37, unknown : inconnue de rockim) ; `dfhPsi0` (l.38, unknown : inconnue de rockim) ; `dfhKPsi` (l.39, unknown : inconnue de rockim) ; `dfhPsiMax` (l.40, unknown : inconnue de rockim) |
| `bench_impact/configs/impact3d_dpdfh_gros.cfg` | fdem3d | `dfhPsiVar` (l.37, unknown : inconnue de rockim) ; `dfhPsi0` (l.38, unknown : inconnue de rockim) ; `dfhKPsi` (l.39, unknown : inconnue de rockim) ; `dfhPsiMax` (l.40, unknown : inconnue de rockim) |
| `tests_f2/tg_3d.cfg` | fdem3d | `thermal` (l.29, garde : 'thermal' n'est implemente que pour mode = fdem (garde main.cpp)) |

