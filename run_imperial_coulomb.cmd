@echo off
rem ---------------------------------------------------------------------------
rem LE RUN DECISIF — replique Imperial College AVEC la plage de mode II publiee
rem (2026-08-27). Deck = copie exacte du run du 26/08 + les deux seules cles
rem  jointShearRange = coulomb  et  jointFrictionScaled = 1  (diff verifie).
rem
rem Binaire DEDIE rockim_cl.exe (commits 4cbd74f..3ab1748) : l ancien exe du
rem run du 26/08 reste intact pour la comparaison A/B.
rem
rem Verdict A/B a historique egal : t_sim ~ 106 us (le run out_imperial s est
rem arrete a 1,0597e-4 s). Si le cisaillement y est toujours nul, c est le
rem correctif qui est refute, pas le run.
rem
rem Le log porte la banniere `[FDEM3D] jointShearRange = coulomb` — une cle mal
rem orthographiee serait IGNOREE EN SILENCE par le lecteur de config.
rem ---------------------------------------------------------------------------
cd /d %~dp0
rockim_cl.exe bench_impact\configs\impact_imperial_coulomb.cfg out_imperial_coulomb
