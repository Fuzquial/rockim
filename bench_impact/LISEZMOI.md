# bench_impact — impacts a insert unique (spec 005)

Reference : Yang, Xiang, Naderi, Wang, Aising, Ugarte, Latham — IJRMMS 191
(2025) 106125 (StAnne/Rhune, 7 criteres) et IJRMMS 206 (2026) 106660 (granite
Kuru, pulverisation). Donnees d'essai : Mines Paris-PSL.

## Chaine complete

    python tools/make_impact_mesh.py meshes/impact_s15.msh 1.5
    rockim_i3.exe bench_impact/configs/impact_stanne_s15.cfg out_imp_stanne
    python bench_impact/tools/fig_impact.py out_imp_stanne --stem bench_impact/fig_stanne_s15
    python bench_impact/tools/gif_impact.py out_imp_stanne --out bench_impact/gif_stanne_s15.gif

Quatre corps (volumes physiques gmsh) : rock / insert carbure / bit / piston
acier. Insert brase au bit (`groupBond.bit.insert = joints`), piston lance par
`groupVel`, roche fixee en base, insertion ADAPTATIVE + DIF `yang-fig2`.
L'echelle 1,5 est la variante econome (~40 k tets) ; l'echelle 1 reproduit le
maillage de leur fig. 6 (~230 k tets).

Les 7 criteres de leur Table 3 sont imprimes par fig_impact.py ; les
predictions sont inscrites dans le bandeau du deck AVANT lancement.

## Reserves connues

- s = 1,5 elargit la barre maillage de leur Table 3 (~10 % sur la
  fissuration) ; la replique fidele demande s = 1 et une nuit dediee ;
- plaque de charge et circlip omis (leur poids sur l'outil n'est pas
  applique dans leurs simulations non plus) ;
- memoire phd/FDEM.md : instabilite latente 3D en phase DEBRIS (gelee le
  2026-08-07) — verifier le bilan d'energie de fin de run avant de citer
  la phase posterieure a ~500 us.
