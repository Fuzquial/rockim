# Rejeu du tunnel EDZ (Wang et al. 2024) avec le contact corrigé

Date : 2026-10-07. Binaire : `/home/user/rockim/build_contact/rockim`, compilé depuis b869dae
(HEAD du dépôt ; aucune source plus récente que le binaire ; les clés `contactCandidates`,
`gcBirth = offset`, `potForceExact` sont reconnues). Maillage : `meshes/tunnel_hs_red.msh`
(22 730 triangles, identique octet pour octet à celui du 6 octobre). Decks : ceux de
`simulations_a_lancer/configs/tunnel/`, copiés sans changement (`red_*_A.cfg`) ; la variante B
ajoute en fin de deck `contactCandidates = vertex`, `gcBirth = offset`, `potForceExact = true`.
Lancés depuis `/home/user/rockim`, `nice -n 10`, 40 min au plus. Sorties, journaux et JSON
dans ce dossier ; chaîne : `chain.txt`, `run_chain*.sh`.

Extraction (`extract.py`, outils du dépôt sans modification) : `analyse()` de
`figures/tunnel/fig_tunnel_rejeu.py` (mêmes définitions que `rejeu_metrics.json`), puis
`tunnel_edz/tools/block_sizes.py --t 0.25 --r 25` (0,38 pour l'intrinsèque) et
`tunnel_edz/tools/nucleation_vs_propagation.py`, enfin le bilan du journal. Le script a été
vérifié sur les sorties du 6 octobre : il redonne les chiffres de D3 §9
(`metrics_ref6oct.json`).

## Ce qui a tourné

| deck | A (clés au défaut) | B (contact corrigé) |
|---|---|---|
| red_adapt | 4 fils, 8 min 43 s ; et 1 fil, 8 min 49 s | 4 fils, 24 min 07 s |
| red_tip16 | 4 fils, 6 min 35 s | 4 fils, 23 min 11 s |
| red_d07 | 4 fils, 17 min 15 s (machine partagée) | 4 fils, 10 min 13 s |
| red_intr | 4 fils : arrêté à t = 0,060 s après 13 min ; repris à 1 fil, 24 min 48 s | non achevé (voir plus bas) |

Intrinsèque B : à 4 fils, t = 0,024 s en 16 min (relâchement à 0,15 s), arrêté. Diagnostic
court à 1 fil, T = 0,004 s (`diag/`, avant toute excavation) : A 12 s, `potForceExact` seul
11 s, `gcBirth = offset` seul 11 s, `contactCandidates = vertex` seul 67 s, B complet 131 s,
soit ×10,5. Cause mesurée : pendant le tassement in situ initial, 2 264 joints sur 33 984
passent à D > 0 (D max 0,095) dès t = 0,004 s ; en intrinsèque ce seuil déclenche l'anneau de
sommets (Fukuda 2021), et 13,0 millions d'éléments-pas sont ajoutés en 1 316 pas, environ
9 900 éléments par pas, 44 % du maillage. Projection pour T = 0,38 s : de l'ordre de 3 h à 1 fil.
Hors budget.

Machine : 4 cœurs, partagée à partir de 9 h environ avec un autre calcul rockim (environ 2 fils,
`bancs_insertion`). À 4 fils, l'intrinsèque A allait environ 3 fois moins vite qu'à 1 fil (fils en
attente aux barrières sur une machine surchargée) ; l'adaptatif A met le même temps à 4 fils et
à 1 fil. Les durées sont donc indicatives ; seuls les rapports B/A pris au même moment
(adapt : ×2,8 ; tip16 : ×3,5) ont un sens.

## Contrôle A contre le 6 octobre

- red_adapt A à 1 fil et red_intr A à 1 fil redonnent le 6 octobre au bit près : mêmes
  rompus (6 397 ; 13 674), mêmes lignes de résumé, mêmes blocs, même propagation. Le binaire
  corrigé, clés au défaut, reproduit bien l'ancien comportement.
- À 4 fils, l'ordre de sommation change et le calcul diverge : adapt 6 107 rompus (−4,5 %),
  R_EDZ p95 −2,3 % ; tip16 −2,1 % de rompus, p95 +0,9 % ; d07 −0,1 %. C'est l'ordre de
  grandeur de la variabilité multifil, à garder en tête pour lire les écarts B/A.
- Seul le résidu B4 change : 13,5 %, 5,2 %, 41,1 % et 1,1 % le 6 octobre, entre 10⁻¹³ et
  10⁻¹⁰ % maintenant, à physique identique (adapt A 1 fil). C'est la correction du bilan des
  rouleaux de b869dae (contrainte appliquée avant la correction leapfrog), pas un effet du
  contact. Le résidu B4 ne discrimine donc plus les variantes ; le poste à lire est le contact.

## Tableau complet

Toutes les grandeurs à t − t_s = 0,23 s (T = 0,25 s ; 0,38 s pour l'intrinsèque). « 6 oct. » :
binaire f0209ef, 1 fil (D3 §9). A : b869dae clés au défaut, 4 fils sauf mention. B : b869dae
avec les trois clés, 4 fils. Rompus : D ≥ 0,999 ; modes = `breakMode` du solveur. Blocs :
disque r ≤ 25 m ; le premier des « 5 plus gros blocs » est le massif intact (sauf en
intrinsèque, où aucun bloc ne domine).

| calcul | rompus | traction / cisaill. | R_EDZ p95 / max [m] | demi-axes h / v [m] | longueur [m] | U paroi moy. / p95 [m] | blocs (r ≤ 25 m) | mono-élém. | 5 plus gros blocs [m²] | propagation |
|---|---|---|---|---|---|---|---|---|---|---|
| adaptatif, 6 oct. | 6397 | 2385 / 4012 | 15,02 / 17,41 | 12,69 / 13,93 | 2425 | 0,094 / 0,128 | 2483 | 74,7 % | 1448,0 ; 16,0 ; 11,5 ; 10,9 ; 10,5 | 29,2 % |
| adaptatif, A (4 fils) | 6107 | 2268 / 3838 | 14,67 / 17,16 | 12,99 / 13,25 | 2323 | 0,091 / 0,124 | 2341 | 73,0 % | 1457,7 ; 12,7 ; 10,9 ; 10,8 ; 10,0 | 27,2 % |
| adaptatif, A (1 fil) | 6397 | 2385 / 4012 | 15,02 / 17,41 | 12,69 / 13,93 | 2425 | 0,094 / 0,128 | 2483 | 74,7 % | 1448,0 ; 16,0 ; 11,5 ; 10,9 ; 10,5 | 29,2 % |
| adaptatif, B | 4654 | 1696 / 2958 | 13,73 / 15,91 | 11,83 / 12,17 | 1769 | 0,083 / 0,127 | 1607 | 67,1 % | 1548,8 ; 14,0 ; 11,7 ; 8,1 ; 7,5 | 24,2 % |
| pointe k = 1,6, 6 oct. | 11846 | 3553 / 8293 | 19,71 / 27,66 | 16,28 / 17,01 | 4770 | 0,135 / 0,186 | 5765 | 85,1 % | 930,3 ; 49,6 ; 36,3 ; 26,0 ; 24,2 | 42,3 % |
| pointe k = 1,6, A | 11594 | 3441 / 8153 | 19,88 / 26,68 | 16,85 / 17,53 | 4753 | 0,142 / 0,216 | 5590 | 83,4 % | 636,0 ; 102,1 ; 89,9 ; 63,1 ; 47,4 | 42,3 % |
| pointe k = 1,6, B | 9870 | 2718 / 7151 | 19,28 / 24,78 | 17,08 / 16,55 | 3875 | 0,122 / 0,172 | 4236 | 76,0 % | 1042,4 ; 53,1 ; 27,8 ; 25,1 ; 16,0 | 35,1 % |
| amortissement 0,7, 6 oct. | 1663 | 596 / 1067 | 9,63 / 10,52 | 8,71 / 8,84 | 637 | 0,025 / 0,035 | 473 | 58,4 % | 1754,9 ; 10,0 ; 5,4 ; 5,2 ; 4,7 | 32,8 % |
| amortissement 0,7, A | 1662 | 617 / 1045 | 9,47 / 10,53 | 8,73 / 8,28 | 638 | 0,026 / 0,041 | 451 | 57,2 % | 1761,5 ; 7,5 ; 7,3 ; 7,1 ; 5,0 | 33,8 % |
| amortissement 0,7, B | 1575 | 623 / 952 | 9,49 / 10,58 | 8,51 / 8,85 | 602 | 0,026 / 0,036 | 389 | 52,2 % | 1760,4 ; 9,0 ; 6,3 ; 5,7 ; 4,2 | 32,8 % |
| intrinsèque, 6 oct. | 13674 | 879 / 12795 | 24,17 / 54,15 | 20,28 / 18,75 | 6979 | 0,186 / 0,248 | 6473 | 79,5 % | 121,1 ; 92,3 ; 63,9 ; 49,7 ; 41,0 | 52,5 % |
| intrinsèque, A (1 fil) | 13674 | 879 / 12795 | 24,17 / 54,15 | 20,28 / 18,75 | 6979 | 0,186 / 0,248 | 6473 | 79,5 % | 121,1 ; 92,3 ; 63,9 ; 49,7 ; 41,0 | 52,5 % |

Bilans (fin de run, J/m). W = « net work injected by general contact » (positif = injecté).
F = frottement dissipé. W + F = part normale (positive = injection par le potentiel normal).
Énergie de naissance = énergie de recouvrement des paires nées en recouvrement (nouvelle ligne
de b869dae ; « materialisable » sous la rampe au défaut, « neutralisee » sous `gcBirth = offset`).

| calcul | travail net du contact W [J/m] | dont frottement dissipé F [J/m] | part normale W + F [J/m] | énergie de naissance [J/m] (paires) | énergie cohésive [J/m] | Cundall [J/m] | résidu B4 | durée |
|---|---|---|---|---|---|---|---|---|
| adaptatif, 6 oct. | 2,69·10⁴ | 1,08·10⁶ | 1,11·10⁶ | non imprimée | 4,29·10⁵ | 3,19·10⁶ | 13,5 % | 7 min 38 s |
| adaptatif, A (4 fils) | −1,95·10⁵ | 8,64·10⁵ | 6,68·10⁵ | 8,42·10⁵ (9940, materialisable) | 4,09·10⁵ | 2,73·10⁶ | 2,35·10⁻¹² % | 8 min 43 s |
| adaptatif, A (1 fil) | 2,69·10⁴ | 1,08·10⁶ | 1,11·10⁶ | 1,30·10⁶ (9928, materialisable) | 4,29·10⁵ | 3,19·10⁶ | 1,58·10⁻¹² % | 8 min 49 s |
| adaptatif, B | −7,44·10⁵ | 6,26·10⁵ | −1,18·10⁵ | 9,33·10¹ (17790, neutralisee) | 2,97·10⁵ | 1,64·10⁶ | 1,47·10⁻¹¹ % | 24 min 07 s |
| pointe k = 1,6, 6 oct. | 1,94·10⁶ | 2,48·10⁶ | 4,42·10⁶ | non imprimée | 9,57·10⁵ | 8,17·10⁶ | 5,2 % | 8 min 42 s |
| pointe k = 1,6, A | 1,53·10⁶ | 2,39·10⁶ | 3,92·10⁶ | 4,35·10⁶ (11543, materialisable) | 9,34·10⁵ | 7,98·10⁶ | 2,80·10⁻¹² % | 6 min 35 s |
| pointe k = 1,6, B | −1,51·10⁶ | 1,17·10⁶ | −3,39·10⁵ | 5,40·10¹ (32914, neutralisee) | 6,71·10⁵ | 3,47·10⁶ | 3,90·10⁻¹¹ % | 23 min 11 s |
| amortissement 0,7, 6 oct. | −2,51·10⁴ | 9,24·10⁴ | 6,74·10⁴ | non imprimée | 9,35·10⁴ | 9,38·10⁵ | 41,1 % | 4 min 34 s |
| amortissement 0,7, A | −4,09·10⁴ | 1,05·10⁵ | 6,41·10⁴ | 8,28·10⁴ (1995, materialisable) | 9,17·10⁴ | 9,04·10⁵ | 1,28·10⁻¹² % | 17 min 15 s |
| amortissement 0,7, B | −1,23·10⁵ | 1,05·10⁵ | −1,79·10⁴ | 1,39·10¹ (5566, neutralisee) | 8,39·10⁴ | 7,57·10⁵ | 2,43·10⁻¹³ % | 10 min 13 s |
| intrinsèque, 6 oct. | 1,89·10⁷ | 8,94·10⁶ | 2,79·10⁷ | non imprimée | 2,66·10⁶ | 3,13·10⁷ | 1,1 % | 22 min 02 s |
| intrinsèque, A (1 fil) | 1,89·10⁷ | 8,94·10⁶ | 2,79·10⁷ | 2,91·10⁷ (7145, materialisable) | 2,66·10⁶ | 3,13·10⁷ | 2,03·10⁻¹⁰ % | 24 min 48 s |

Sous B, les renaissances (paires retrouvées après un pas sans recouvrement) neutralisent en
plus 41 J/m (adapt), 26 J/m (tip16) et 0,5 J/m (d07). Candidats ajoutés par `vertex` :
9,6·10⁷, 2,0·10⁸ et 4,2·10⁷ éléments-pas.

Lecture énergétique. Au défaut, l'énergie de naissance matérialisable dépasse l'énergie
cohésive : 0,84 à 1,30 MJ/m contre 0,41 à 0,43 (adapt), 4,35 contre 0,93 (tip16), 29,1 contre
2,66 (intrinsèque). Sous B, elle tombe à 14–93 J/m et la part normale du contact devient
négative (−0,02 à −0,34 MJ/m) ; d'après la documentation, le biais du compteur est négatif,
donc ces valeurs ne signalent plus d'injection. Le contact n'injecte plus d'énergie dans
aucune des trois variantes B qui ont tourné.

## Réponses aux cinq questions

1. Adaptatif sans facteur de pointe, contact corrigé : l'anneau reste granulé. Blocs
   mono-élément 67,1 % (1 079 blocs sur 1 607), contre 73,0 % en A (4 fils) et 74,7 % le
   6 octobre. Aire moyenne hors massif 0,21 m², deuxième bloc 14,0 m². La correction retire
   environ 6 points de mono-éléments, 24 % des rompus (4 654 contre 6 107), 24 % de longueur
   fissurée, et 6 à 7 % du rayon d'EDZ (p95 13,73 contre 14,67 m). Ces écarts dépassent la
   variabilité multifil (environ 5 %), mais le motif ne change pas : on n'obtient pas de blocs
   métriques. D'un fil à l'autre, la part de mono-éléments varie de 1,7 point (adapt A
   contre 6 octobre : 73,0 contre 74,7 % ; tip16 : 83,4 contre 85,1 %) ; la baisse de 6 points
   dépasse ce bruit, sur un seul tirage.

2. Facteur de pointe k = 1,6, contact corrigé : son effet persiste et grandit. Par rapport à
   adapt B, tip16 B a ×2,12 rompus (9 870), ×2,19 de longueur, R_EDZ p95 19,28 contre 13,73 m
   (+40 %), U paroi moyen ×1,47, une propagation de 35,1 % contre 24,2 %, 76,0 % de blocs
   mono-élément contre 67,1 %, et des blocs secondaires plus gros (53, 28 et 25 m² contre 14, 12
   et 8). En A, le rapport des rompus était ×1,90. Hypothèse du doctorant (« codé pour corriger
   une surfracturation ») : aucun des huit calculs ne la soutient. Le facteur ne réduit la
   fissuration dans aucun cas ; il l'augmente, avant comme après la correction. La correction
   du contact, elle, réduit la fissuration (adapt −24 %, tip16 −15 %). Le commentaire du deck
   donne un autre motif : relever la part de propagation (43,7 % contre 56,8 % pour
   l'intrinsèque en production). Sur ce critère, k = 1,6 fait encore passer la propagation de
   24 à 35 %, mais au prix d'une fissuration doublée et d'un anneau plus granulé. Dire si c'est
   un apport ou un biais demande un critère extérieur (convergence en maillage, mesure de
   terrain), que ce rejeu ne fournit pas.

3. Intrinsèque : non tranché avec le contact corrigé, le calcul B n'a pas tenu dans le budget
   (voir plus haut). Avec l'ancien contact (A, au bit près le 6 octobre), les fissures atteignent
   le bord (R_max 54,15 m pour un demi-côté de 50 m). Mais ce calcul matérialise 29,1 MJ/m
   d'énergie de naissance, 11 fois l'énergie cohésive (2,66 MJ/m), et sa part normale de contact
   vaut +27,9 MJ/m. La fissuration jusqu'au bord vient donc d'un calcul qui injecte de l'énergie ;
   on ne sait pas si elle survit à la correction. À noter aussi : sous B, l'intrinsèque déclenche
   l'anneau de sommets dans 44 % du maillage avant l'excavation, à cause des D > 0 du tassement
   in situ.

4. Rayon d'EDZ face aux 19 m de Wang : adapt B donne R_max 15,91 m et p95 13,73 m. La correction
   l'éloigne un peu de 19 m (A : 17,16 / 14,67 m ; 6 octobre : 17,41 / 15,02 m). tip16 B donne
   p95 19,28 m et R_max 24,78 m : il encadre 19 m, mais avec deux fois plus de rompus et un
   max qui dépasse la référence de 30 %. d07 B : 10,58 / 9,49 m. Cette comparaison ne peut pas
   départager les variantes : maillage trop grossier, T non convergé, définition du rayon
   différente de celle de Wang (max et p95 des milieux de joints rompus).

5. Amortissement 0,7 : l'écart reste important mais se réduit à environ 3. Rapports adapt/d07
   en B : rompus 2,95, longueur 2,94, U paroi moyen 3,25, R_EDZ p95 1,45. En A : 3,67, 3,64,
   3,48 et 1,55 ; le 6 octobre : 3,85, 3,80, 3,78. d07 bouge à peine avec la correction
   (−5 % de rompus, dans la variabilité) alors que l'adaptatif à 0,15 perd 24 %. Une partie
   (environ 20 %) de l'ancien « facteur 4 » venait donc de l'injection du contact dans le calcul
   peu amorti. Le reste tient à l'amortissement lui-même.

## Réserves

- Maillage réduit (22 730 triangles, h_fine = 0,40 m > ℓ_cz/2 = 0,278 m) : il viole la règle
  d'objectivité dx < ℓ_cz/2. La part de blocs mono-élément dépend directement de la taille
  d'élément ; ces chiffres ne valent que comme tendance.
- T = 0,25 s (0,23 s après le début du relâchement) : la convergence de paroi n'est pas
  atteinte (+20 % de U sur le dernier cinquième du run au 6 octobre), donc EDZ et U sont des
  valeurs à un instant donné, pas des états stabilisés.
- 4 fils : les calculs ne sont pas reproductibles au bit près ; un seul tirage par variante.
  Les écarts B/A sous environ 5 % (tip16 p95 −3 %, d07 tout entier) ne sont pas significatifs.
- Machine partagée pendant une partie des calculs : les durées ne sont pas comparables entre
  elles au-delà des rapports B/A donnés.
- Intrinsèque B absent : les questions 3 et la comparaison intrinsèque/adaptatif sous contact
  corrigé restent ouvertes. Un calcul B de l'intrinsèque coûte environ 3 h à 1 fil sur ce
  maillage.
- `jointShearUnload = origin` reste signalé non conservatif par le solveur (toutes variantes).
- Les blocs (block_sizes) et la propagation sont pris à la dernière trame et entre trames
  successives (12 ou 18 trames) : la part de propagation dépend de la cadence des trames.

## Fichiers

- `metrics_fix.json` (A, A 1 fil, B), `metrics_ref6oct.json` (6 octobre repassé par le même
  script), `m_red_*.json`, `tables.md`, `table.py`, `extract.py`.
- `diag/` : diagnostic de coût de l'intrinsèque (T = 0,004 s, 1 fil, une clé à la fois).
- `out_red_intr_*_arrete/` et journaux `*_arrete.log` : calculs interrompus, non exploités.
- Figure : `/home/user/rockim/docs/rapport_guide/figures/tunnel/fig_tunnel_fix.pdf`, script
  `fig_tunnel_fix.py` à côté (aperçu PNG ici). `apercu.py` et `apercu_adapt_AB.png` ont été
  écrits dans ce dossier à 7 h 54 par un autre processus ; ils ne viennent pas de ce rejeu.
