# bench_impact/donnees — sorties de run COMPRESSEES (.npz)

Les runs d'impact tournent sur la machine du doctorant ; le depouillement se
fait ailleurs. `tools/pack_run.py` comprime un dossier de sortie complet
(history.csv + frames VTU, ~355 Mo) en UN .npz de quelques Mo, transportable
par git.

    python tools/pack_run.py out_pulv_coulomb bench_impact/donnees/P1.npz

Ce qui est garde : history.csv en entier (float32) ; par frame, les elements
de ROCHE proches de l'impact (positions, bulkD, vonMises) ; les joints de la
DERNIERE frame seulement — tBreak etant cumulatif, elle porte toute la
chronologie de fissuration ; l'enveloppe verticale de l'outil.
Ce qui est jete : les tets d'outil, les champs redondants, les frames de
joints intermediaires. Mesure sur le jumeau du 29/08 : 224 Mo -> 4,5 Mo.

NB : le dossier s'appelle `donnees` et non `runs`, `runs/` etant ignore par
`.gitignore` (les sorties BRUTES ne doivent jamais etre committees — seule
la forme compressee l'est).

Ces .npz sont des artefacts de TRAVAIL, lies a leur deck. A la cloture d'une
etude, l'archive definitive part dans la base de these (`phd_geothermie`,
convention des `*_data.npz`) ; les VTU bruts ne quittent jamais la machine
de calcul.
