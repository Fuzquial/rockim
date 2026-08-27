# bench_polyaxial — le repere COMPRESSED-SHEAR-TO-FAILURE (Guo 2014, §3.5)

Le test qui manquait a la suite : un joint COMPRIME mene a la rupture en
cisaillement. C'est l'absence de ce repere qui a laisse la plage de mode II
cohesion-seule (l'erreur de transcription corrigee par `jointShearRange =
coulomb`, commits cf562d1/48fd6ba) passer 40 tests verts pendant des mois —
aucun repere de la suite ne chargeait ce regime.

Source : Guo 2014 (these Imperial), §3.5 « Polyaxial compression tests »,
Table 3.3 et fig. 3.15-3.18. Transpose en TRIAXIAL equivalent (confinement
lateral -3,75 MPa + compression axiale) : Mohr-Coulomb ignore sigma2, raison
meme pour laquelle Guo a choisi sigma2 = -9,375 intermediaire — le cercle de
Mohr est identique.

## Lancer

    rockim bench_polyaxial/polyaxial_guo_t3.cfg out_poly_t3
    rockim bench_polyaxial/polyaxial_guo_t1.cfg out_poly_t1

Cube 500 mm, ~68 000 joints, ~2 min sur 2 threads (2026-08-27, conteneur
4 coeurs). Materiau Table 3.3 : E 30 GPa, nu 0,27, ft 1,0 MPa, Gf 60 J/m2
(unique, modes I et II), mu 0,6. t3 : c 3,8 MPa / phi 12,4 deg ; t1 :
c 0,5 MPa / phi 33,4 deg. Seuil theorique (au confinement effectif -4,0
MPa) : sigma1f = N phi * 4,0 + UCS ~ 15,6-15,7 MPa pour les deux jeux.

## Valeurs mesurees le 2026-08-27 (binaire 48fd6ba, 18x18x18 jitter 0,3)

| deck | loi | rompus / 68 040 | mode | pic macro |
|---|---|---|---|---|
| t3 | coulomb (deck) | 7 824 | 100 % cisaillement | 16,7 MPa |
| t3 | cohesion (cle retiree) | 5 104 | 100 % cisaillement | 17,3 MPa |
| t1 | coulomb (deck) | 11 893 | 100 % cisaillement | 23,3 MPa |
| t1 | cohesion (cle retiree) | **0** | — | **33,1 MPa = 2,1 x le seuil** |

Le verdict du repere est la ligne t1/cohesion : avec la plage 3 GfII/c
(cohesion seule), delta_c = 3*60/0,5e6 = 360 um et le cube ENCAISSE plus de
deux fois sa contrainte de rupture de Mohr-Coulomb sans UN SEUL joint rompu
— « plus la cohesion est faible, plus la roche est incassable ». C'est le
mecanisme exact qui interdisait le noyau broye de l'impact a insert (BILAN
du 2026-08-27). La loi publiee (plage divisee par fs = c + tan(phi)
|sigma_n|, these eq. 2.24 + 2.30) casse 11 893 joints, tous en cisaillement.

Critere de PASS propose pour la suite : t1_coulomb rompt > 5 000 joints en
mode cisaillement dominant ET t1_cohesion (controle) n'en rompt aucun ; t3 :
pic macro dans [15 ; 18] MPa pour les deux lois. A affiner en repere
chiffre (err_pct) lors de l'integration au tier fast.
