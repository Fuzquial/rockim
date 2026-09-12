# Campagne de correction du 13/09/2026 (après-midi) — cadrage pour les agents

Demande de Fernando : « prends tout cela en considération et lance une grande campagne de correction à
partir de ces documents ». Sources à lire AVANT de travailler (dans `docs/`) :
`DIAGNOSTIC_ICL_independant_2026-09-12.md`, `COMPLEMENT_YANG_sources_et_corrections_2026-09-12.md`,
`ECARTS_guo2014_rockim_2026-09-13.md` (§5 = relecture des deux premiers), `AUDIT_2026-09-13.md` §6-§10.
Dépôt : `C:\Users\fuzquianoalricabi\simulations\FDEM\rockim_g1` (git, branche `g1`). Binaire courant
`rockim_g1y16.exe` ; le prochain sera **`rockim_g1y17.exe`**.

## Règles NON NÉGOCIABLES (toutes tâches)

1. **Croissance par addition** : toute capacité nouvelle est une clé opt-in ; le chemin par défaut reste
   textuellement intact et bit-identique (ancre `tools/bitid.py`, 8 decks + `--refs tools/bitid_refs_jointlaw.json
   --only fdem3d_yang_v2_court`). Ne rien enlever, ne rien renommer.
2. **C++, LaTeX, chemins Windows : uniquement via les outils Edit/Write**, jamais par heredoc Bash ni `python -c`
   (Bash mange les anti-slash). Après toute écriture : aucun caractère de code < 32 hors tab/LF/CR.
3. **Miroir 2D** : une clé de loi de joint ou de sortie ajoutée en 3D (`src/Fdem3dSolver.cpp`) est ajoutée en 2D
   (`src/FdemSolver.cpp`) quand la notion existe ; sinon dire pourquoi dans le commentaire.
4. **Registre des clés** : après avoir ajouté une clé, lancer `python tools/gen_keys_by_mode.py` (régénère
   `include/rockim/KeysByMode.hpp` et `tools/keys_by_mode.json`).
5. **Documentation** : une ligne par clé dans `DOCUMENTATION_rockim.md` (tableau des clés, même style que
   `jointPenaltyLength`), et des puces dans `CHANGELOG.md` sous le titre
   `### Campagne de correction du 13/09 (après-midi) — \`rockim_g1y17.exe\`` (Edit, jamais Write du fichier entier).
6. **Ne JAMAIS** : lancer un run long (> 5 min), tuer un processus, toucher aux dossiers `out_yang_bench_*`,
   `out_yang2026_v3`, aux journaux `results/queue_*.log`, ni au fichier en cours `results/yang_bench_*.log`
   (des bancs tournent sur 14 fils jusqu'à ~17 h). Pas de `git commit` ni `git add` (l'orchestrateur committe).
7. **Build** : `powershell -ExecutionPolicy Bypass -File tools/build.ps1 -Jobs 4` (incrémental, depuis la
   racine) pour vérifier que ça compile ; le build propre (`-Clean`) et l'ancre sont faits UNE fois à la fin
   par l'agent « Build ». Ne pas copier `build/rockim.exe` sur un nom `rockim_g1y*.exe` (réservé à Build).
8. **Mini-tests** : chaque changement de code est accompagné d'un test court (< 2 min, `OMP_NUM_THREADS=4`,
   sorties dans `C:\Users\FUZQUI~1\AppData\Local\Temp\claude\...\scratchpad` ou `tests_f2/campagne13/`)
   avec un critère FALSIFIANT (une variante qui DOIT échouer ou une égalité qui DOIT tenir à 1e-9).
9. Le run s = 1 de référence (dernière trame 18, `out_yang2026_v3/fdem3d_joints_0018.vtu`, filtre
   `tBreak >= 0`) et le témoin s = 2,5 (`out_yang_bench_s25_v3P/`, 200 µs) servent de données d'essai.
10. Rapport de fin de tâche : fichiers touchés, clés ajoutées (nom, défaut, modes), tests exécutés avec
    leurs chiffres, ce qui n'a PAS été fait et pourquoi. Pas d'affirmation non mesurée.

## Tâches

### S1 — Instrumentation de rupture et de contrainte (solveur, 3D + 2D)
Motif : DIAGNOSTIC §5 et §6.1 ; ECARTS §5. (a) Champs VTU des joints : `dead` (0/1) et `openMax`
(ouverture géométrique normale maximale sur les points, m) en plus de `damage`, `tBreak`, `breakMode`.
(b) Clé opt-in `jointBreakModeRef = slipF | slipRef` (défaut `slipF`, historique) : sous `slipRef`, la
partition rn/rs de `failMode`/`bmode` à la rupture est normalisée par la même plage que le moteur
(`slipRef`, pression courante) et non par `J.slipF`. (c) Champ VTU des éléments `pMean` (pression moyenne,
traction > 0) et, dans le journal à chaque trame, le nombre d'éléments écrêtés par `meanTensionCapFactor`
et l'excès max (compteur, aucune force changée). (d) Mini-test : sur la trame 18 du s = 1 rejouée 1 µs
(`T` court sur le deck `configs/yang2026_impact_v3_plastic.cfg` ? trop long : utiliser le deck s = 2,5
`configs/yang2026_bench_s25_v3P.cfg` avec `T = 5e-6`, `frames = 1`) vérifier que `dead` ⊂ `tBreak >= 0`,
que le compteur d'écrêtage est imprimé, et que les sorties sont bit-identiques quand les clés sont absentes.

### S2 — Force de contact entre corps nommés (solveur 3D)
Motif : DIAGNOSTIC §6.1 (« exporter les forces de contact piston/bit et insert/roche »). Clé opt-in
`contactForcePairs = insert:rock bit:piston ...` : à chaque ligne de `history.csv`, la somme des forces de
contact (normale + tangentielle, composantes x y z) exercées par le premier corps sur le second, colonnes
`Fc_<a>_<b>_x/y/z`. Somme faite dans la boucle de contact par potentiel (paires élément-élément → groupes
via `elemGroup_`), réduction OpenMP sûre (accumulateurs par fil). Mini-test : deck s = 2,5 court
(`T = 3e-5`, l'onde n'a pas atteint la roche : `Fc_insert_rock` doit rester nulle, `Fc_bit_piston_z`
non nulle après le contact piston/bit à ~6,7 µs ; et la réaction de la roche à 200 µs sur le témoin
`out_yang_bench_s25_v3P` si le deck peut être rejoué 50 µs : comparer à l'estimateur quantité de mouvement
de `tools/fig_fp.py`).

### S3 — Banc de joint cinématique (solveur 3D, scénario nouveau)
Motif : DIAGNOSTIC §6.2 (« valider une loi de joint unique sur un petit banc 3D »). `scenario = jointbench`
(fdem3d) : deux tétraèdres partageant une facette (un seul joint), tétra A fixe, tétra B en translation
rigide prescrite le long d'un chemin ; sortie `jointbench.csv` (t, dn, ds, sigma, tau, D par point, rompu).
Clés : `jbMode = tension | shear | mixed`, `jbAmp` (m), `jbNormal` (m, offset normal constant pour les
essais de cisaillement sous pression), `jbUnloadAt` (fraction de jbAmp où l'on décharge jusqu'à 0 puis
recharge), `jbRate` (m/s). Pas de dynamique parasite : masses grandes ou `dampingLocal` fort, ou mieux,
chemin quasi statique. Critères FALSIFIANTS à imprimer : (i) sous `jointShearUnload = solidity` la
décharge avant rupture retrace la charge (|σ_recharge − σ_charge| < 1e-6 ft) ; (ii) sous `plastic` un
glissement résiduel subsiste après décharge en cisaillement ; (iii) l'aire sous σ(dn) jusqu'à rupture vaut
Gf (± 2 %) sous `jointDeltaC = exact`/`plastic` et ≈ 1,159 Gf sous `solidity` (ot = 3 Gf/ft, ∫z = 0,3863) ;
(iv) la règle deux points sur trois : avec une rotation de B (option `jbTilt`), les points cèdent à des
instants différents et le joint ne meurt qu'au deuxième. Comparer `plastic`, `origin`, `solidity` sur le
même maillage et la même pénalité (`jointPenaltyLength = edge`, facteur 25).

### S4 — Pulvérisation : mesures alternatives et essais élémentaires (solveur 3D)
Motif : DIAGNOSTIC §3 ; COMPLEMENT §7. Clés opt-in : `bulkDamageLength = inscribed | edge` (défaut
`inscribed` = h_e actuel ; `edge` = arête moyenne du tétra) et `bulkDamageStrain = deviatoric | principal |
total` (défaut `deviatoric` = ε_vm actuel ; `principal` = plus grande déformation principale en valeur
absolue ; `total` = norme de Frobenius √(2/3)‖ε‖ trace comprise). Essais élémentaires (deck fdem3d ou
fem3d à quelques tétras, scénario existant le plus proche : chercher `scenario = triax`/`tension` dans
`configs/` et `tests_f2/`) : compression isotrope (D doit rester 0 sous `deviatoric`, s'armer sous
`total`), compression uniaxiale, cisaillement ; imprimer δ_m, D, et le seuil en déformation pour chaque
mesure (attendu : 14 µm / 0,51 mm = 2,7 % sous `inscribed`, 14 µm / 1,37 mm = 1,0 % sous `edge`).

### T1 — Fissures connectées et cratère (outil Python)
Motif : DIAGNOSTIC §5 (rayon max des centroïdes ≠ longueur de fissure). `tools/crack_paths.py` : à partir
du VTU des joints (filtre `tBreak >= 0`, `bench_impact/tools/imp_lib.py`), composantes connexes des
facettes rompues (adjacence par arête partagée), séparation noyau (composante contenant l'axe, r < 6 mm)
/ fissures périphériques ; pour chaque composante périphérique : longueur radiale (r max − r min),
orientation moyenne (radiale si |n·e_r| < 0,35), profondeur, nombre de facettes ; « longueur de fissure
radiale » au sens de Yang = distance du centre à la pointe de la plus longue composante radiale connectée
au noyau ; rayon de cratère = frontière extérieure des facettes de surface (z > −1 mm) du noyau. Sortie :
tableau + figure (vue de dessus, une couleur par composante, PDF vectoriel + PNG). Appliquer à la trame 18
du s = 1 et à la dernière trame de `out_yang_bench_s25_v3P`.

### T2 — Estimateurs de Yang dans les figures (outil Python)
Motif : COMPLEMENT §5. Dans `tools/fig_kinetics.py` et `tools/yang_report.py` : vitesse d'indentation =
pente de la portion linéaire (10-90 % de l'enfoncement) du déplacement de l'insert ; vitesse de rebond =
pente de la portion linéaire remontante après le retournement ; contrainte de référence = pic de la
PREMIÈRE onde à mi-bit ; imprimer les deux estimateurs (pente et max instantané) côte à côte avec les
fenêtres. Vérifier sur `out_yang2026_v3` (attendu ≈ 6,86 m/s pente, 7,37 max).

### T3 — Série de maillages à train figé et masses (mailleur)
Motif : DIAGNOSTIC §4 et §6.5 ; ECARTS lignes 23-24. Avec `tools/make_impact_mesh.py` (argument SR =
échelle de la roche, s = échelle du train ; preset `quality=hxt`) : générer `meshes/impact_yang_train1_rock25_hxt.msh`,
`_rock137` (SR = 1, l'actuel), `_rock073` (SR = 0,73 : arête médiane visée 1,0 mm dans la boule), tous avec
le MÊME train à s = 1 ; mesurer avec `tools/mesh_quality.py` (h min, arête médiane dans la boule R 12,5,
nombre de tétras par corps) et imprimer les masses par corps (volume × ρ : 2626 roche, 7850 acier, 15250
carbure) ; comparer aux volumes analytiques des cylindres et à Yang (piston 1,173 kg, bit 1,509 kg) ;
expliquer l'écart (facettisation, longueurs, évidement) et proposer la correction (densité par corps via
deux phases acier, ou géométrie). Rapport `docs/MAILLAGE_serie_2026-09-13.md`. Ne pas lancer de run.

### T4 — Decks de conformité et cas St Anne (decks, sans lancement)
Motif : COMPLEMENT §2 (« commencer par St Anne »), DIAGNOSTIC §4 (`meanTensionCapFactor = 0`). Écrire :
(a) `configs/stanne2025_bench_s25.cfg` — impact à insert unique sur St Anne (Yang 2025 Table 4 : ρ 2731,
E 57 GPa, ν 0,31, ft 7,0, c 18,8, tanφ 1,0 → frictionDeg 45, GI 12, GII 800, glissement 0,6 ; ARMA 2024 :
p0 3 000 GPa → `jointPenaltyLength = edge`, facteur p0/(2E) = 26,3 ; `bulkDamage = off` ; DIF 2025 ;
amortissement INCONNU : deux variantes `_visc0` et `_visc` avec `bulkViscosity = 2000` marquée
exploratoire) sur le maillage s = 2,5 à train figé (T3) ; (b) `configs/yang2026_bench_s25_v4_*.cfg` :
la série Kuru à partir du témoin A avec `meanTensionCapFactor = 0` et `jointBreakModeRef = slipRef` ; (c)
pour chaque deck, un en-tête avec la liste EXACTE des clés qui diffèrent du témoin, le coût estimé (à
partir de `results/yang_bench_s25_v3P.log` : 5 498 s pour 200 µs à 14 fils avec la machine partagée ; ~20
min seul) et le critère de succès ; (d) un tableau récapitulatif dans `docs/DECKS_conformite_2026-09-13.md`
pour la validation de Fernando. Vérifier chaque deck par un lancement de 2 µs (`T = 2e-6`, `frames = 1`)
qui doit démarrer sans erreur de clé.

### B — Build propre, ancre, tests (après S1-S4)
`tools/build.ps1 -Clean -Jobs 8`, `sha256`, copie en `rockim_g1y17.exe`, `python tools/bitid.py --exe
rockim_g1y17.exe --json results/bitid_g1y17.json > results/bitid_g1y17.log`, puis `--refs
tools/bitid_refs_jointlaw.json --only fdem3d_yang_v2_court`. Exiger 8/8 + 1/1 IDENTIQUE ; sinon dire
quel deck diffère et ne PAS copier l'exe. Rejouer les mini-tests de S1-S4 avec `rockim_g1y17.exe`.

### V — Vérification adversariale (après B)
Pour chaque tâche S1-S4 : deux relecteurs indépendants (lentille « bit-identité des défauts » : lire le
diff `git diff` des sources et prouver que sans la clé rien ne change ; lentille « conformité au cadrage » :
le comportement annoncé est-il celui codé, les tests sont-ils falsifiants, le 2D est-il en miroir ?).
Verdict REFUTÉ par défaut si un doute subsiste.

### Z — Synthèse
Critique de complétude contre les deux documents du 12/09 (chaque recommandation : faite / préparée /
non faite et pourquoi) ; rapport `docs/CAMPAGNE_correction_2026-09-13_RAPPORT.md` ; liste des decks à
valider par Fernando avec clés et coût.
