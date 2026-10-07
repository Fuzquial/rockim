# Revue adversariale : correction du contact par potentiel (rockim, diff non commité du 2026-10-07)

Lecture seule. Aucun fichier du dépôt modifié. Un seul essai hors dépôt (`scratchpad/revue/t.cpp`, équilibre
des forces de `exactGradientForces` à |X| ~ 100, voir §Vérifié).

Décompte : 0 bloquant, 6 à corriger, 8 mineurs.

## À corriger

### C1. `gcBirth = offset` : e0 périmé quand le recouvrement disparaît sans passer par la branche E ≤ e0
- `src/FdemSolver.cpp:7360` (`if (!pot::pairForce(...)) continue;` avant tout accès à `potFt_`),
  `:7384` (`if (isNewH) H.e0 = ePair;`), `include/rockim/PotentialContact.hpp:374`. `potFt_` n'est jamais
  effacé (seule occurrence : `try_emplace`, l. 7368).
- Le cliquet ne s'applique que si la paire est évaluée AVEC recouvrement et E ≤ e0. Dans trois cas e0 reste
  à sa valeur de naissance (grande) :
  1. le recouvrement passe de E > e0 à vide en un pas (fragment éjecté, choc) ;
  2. la paire sort du jeu de candidats (`act_` recomposé, élément hors boîte `einb`) ;
  3. un élément est inversé ou dégénéré (`pairForce` rend false).
- Scénario : paire née d'un joint mort avec S0 profond (e0 = 1e3 J/m), séparée brutalement, puis nouvelle
  approche. Comme `isNewH` est false, la paire garde e0, U = 0 et la force est nulle jusqu'à E > e0. La
  pénétration est libre jusqu'au niveau d'énergie de la première naissance : c'est un contact manqué, pas
  une création d'énergie. Même dans le cas propre (séparation progressive), e0 reste égal au dernier E non
  nul et n'est jamais remis à 0.
- Correction : sous `birthOffset_` seulement (le défaut reste intact), dater la dernière évaluation de
  chaque paire (`H.lastEval`). Si elle n'a pas été évaluée au pas précédent, la traiter comme une
  renaissance (e0 = E courant). Ou bien : quand `pairForce` rend false pour une paire énumérée, poser
  `e0 = 0` dans l'entrée si elle existe.

### C2. Affirmation fausse : « e0 = 0 (paire née par approche) redonne Munjiza à l'identique »
- `include/rockim/PotentialContact.hpp:372`, `DOCUMENTATION_rockim.md:907`, CHANGELOG, IMPLEMENTATION.md.
- La naissance a lieu à la première évaluation où `pairForce` a rendu true, donc avec une aire > 1e-300 :
  e0 = E(S_1) > 0 pour TOUTE paire, y compris une paire née par approche. Aucune paire n'est donc en Munjiza
  pur sous `offset`. Une paire en contact persistant garde un facteur 1 − (e0/E)² < 1 et une pénétration
  décalée d'environ la pénétration du premier pas (v_rel·dt) tant qu'elle ne s'est pas séparée.
- Gravité : comportement voisin de la rampe, qui démarre aussi à sc = 0, mais la doc décrit une propriété
  que le code n'a pas.
- Correction : corriger le texte. Si l'on veut vraiment e0 = 0 pour les paires nées par approche, il faut
  un critère explicite, par exemple naissance d'un joint mort, ou E(S0) au-dessus d'un seuil.

### C3. « Défaut bit-identique » : faux sous `budgetAbortPct > 0` avec `lateralRollers = true`
- Le correctif des rouleaux (`src/FdemSolver.cpp:9283-9308`, `:9409-9426`) modifie `biasW_` et `lysWork_`.
  Ces deux postes entrent dans `resid` de `checkEnergyAbort()` (`:8270-8282`), qui peut arrêter le run.
- Sur un deck à rouleaux avec moniteur actif, l'instant d'arrêt change, donc la trajectoire écrite aussi
  (avant : résidu parasite de −136 005 J/m sur le smoke). C'est voulu (l'ancien résidu était faux), mais la
  promesse « trajectoire bit-identique, bilan seulement » est conditionnelle.
- Aucun deck existant ne combine les deux : les tunnels posent `budgetAbortPct = 0`, et les decks b7 et s25
  sont en fdem3d.
- Correction : le dire dans le CHANGELOG et la DOCUMENTATION.

### C4. Registre des clés périmé : `contactCandidates` y est déclarée propre à `fdem`
- `include/rockim/KeysByMode.hpp:30` et `:772`, `tools/keys_by_mode.json`.
- `Fdem3dSolver.cpp` lit pourtant `gets("contactCandidates")`. J'ai relancé `scan()` du générateur à blanc :
  lecteurs = {fdem, fdem3d}, donc la clé est commune. Le registre a été généré avant le portage 3D.
- Effet aujourd'hui : aucun à l'exécution. `keysbymode::check()` n'est plus appelé, et KeyGuard ne juge que
  les clés non consommées, alors que le 3D la consomme toujours. En revanche :
  - `gen_keys_by_mode.py --scan-decks` déclarerait « refusé » tout deck fdem3d portant la clé ;
  - les compteurs `n_common` et `n_mode_specific` sont faux.
- Correction : relancer `tools/gen_keys_by_mode.py`, puis remettre les séparateurs comme fait précédemment.

### C5. Tests : l'intégration au solveur des trois clés n'est jamais testée en tier fast
- P3 (`src/FdemSolver.cpp:11330-11440`) réimplémente lui-même la naissance : son propre `born`, son propre
  `e0` et `pot::birthOffsetScale`. Il ne passe pas par `potentialContact()`, `potFt_`, `isNewH`, `cplDf`, le
  frottement ou l'ordre `exactGradientForces` puis facteur. Exemples :
  - une inversion de la ligne 7384 (`H.e0` laissé à 0) resterait verte en fast ;
  - un oubli de l'appel `exactGradientForces` dans le solveur aussi.
- Les deux repères full sont faibles :
  - `contactfix_tension_2d` verrouille une valeur, il ne contrôle pas une propriété ;
  - `zeroload_contactfix_2d` (gcwork 0 ± 1e-12) ne distingue pas les clés du défaut, puisqu'à charge nulle
    le défaut ne crée rien non plus.
- Il manque :
  - un test de complétude de `vertex` : éventail avec un joint qui meurt, où la paire à sommet commun naît
    profonde sous `active` et pas sous `vertex` ;
  - un repère pour le résidu B4 des rouleaux (reconnu dans IMPLEMENTATION) ;
  - des tests de refus des clés : valeurs invalides, `offset` + `gcBirthTau`, `vertex` sous penalty.
- La regex `birthe` (`tools/verify_suite.py:51`) est définie mais aucun repère ne l'utilise.
- Correction : ajouter un mini-deck 2D rapide (deux blocs, un joint qui meurt sous compression), avec trois
  critères :
  - `birthe` grand sous `ramp`, petit sous `vertex` ;
  - `gcwork ≤ 0` sous les trois clés ;
  - un résidu B4 à rouleaux avec `lateralRollers`.

### C6. `vertex` sous insertion intrinsèque : l'anneau couvre tout le maillage dès t = 0
- `src/FdemSolver.cpp:7243-7250` et le miroir 3D (`Fdem3dSolver.cpp`, bloc 1b, boucle sur `jt_`).
- En intrinsèque, tous les joints sont `!bonded` dès l'initialisation (`bonded` n'est mis à true qu'en
  `adaptive` ou `none`, l. 1103-1104). Le critère « joint INSÉRÉ » rend donc tous les éléments candidats à
  chaque pas. Les deux références citées déclenchent plus tard : Guo active l'anneau à la rupture du joint,
  Fukuda au début de l'adoucissement.
- Coût mesuré par l'auteur : ×2 à ×4 en 2D, ×13 en 3D à charge nulle. De plus, les paires à sommet commun
  du continuum intact sont évaluées en permanence (recouvrements réels mais minuscules, −1,5e-23 J/m).
- Correction : en intrinsèque, déclencher sur un joint qui s'adoucit ou qui est mort (D > 0 ou `dead`), et
  garder « inséré » pour l'adaptatif. À défaut, documenter que `intrinsic + vertex` équivaut à « tous les
  éléments ».

## Mineurs

- M1. `src/FdemSolver.cpp:10179` : le bilan dit « matérialisée sans travail par la rampe… à hauteur de
  1 − relax ». C'est trompeur : à géométrie figée la rampe matérialise TOUT E(S0) quel que soit τ ;
  1 − relax n'est que la part du premier pas. De même, l'avertissement n'apparaît que si τ < dt (l. 1455+),
  ce qui laisse croire que la rampe est neutre au-delà. Reformuler.
- M2. `DOCUMENTATION_rockim.md:481` (ligne `contact`, modifiée par le diff) garde « Relève de naissance par
  aire… ne matérialise pas son énergie potentielle — signe absorbant garanti ». Cela contredit la nouvelle
  ligne `gcBirth = offset` et l'enquête.
- M3. `offset` sans `potForceExact` est accepté sans un mot. Or P3 mode 2 mesure 0,12 de dérive dans ce cas,
  et IMPLEMENTATION dit que « les deux correctifs vont ensemble ». Ajouter un avertissement à l'init
  (`FdemSolver.cpp:~1441`).
- M4. Cohérence 2D/3D des messages :
  - en 3D, `gcBirth = offset` est refusé avec « ramp | penalty | relay » (`Fdem3dSolver.cpp:1849`), sans
    dire que la valeur existe en 2D seulement ;
  - le message d'erreur 3D de `contactCandidates` ne détaille pas les valeurs, contrairement au 2D ;
  - pas de ligne « naissances » ni d'avertissement τ < dt en 3D.
- M5. `potForceExact` est lu par `getb` (`src/Config.cpp:126`), qui n'accepte que `1/true/yes/on` en
  minuscules : `potForceExact = True` ou `= ture` passe silencieusement à false. Convention préexistante,
  mais une clé opt-in nouvelle en hérite.
- M6. `overlapIntegrals`, `pairEnergy` et `exactGradientForces` (`PotentialContact.hpp:275-358`) ne testent
  pas l'orientation ni la dégénérescence des triangles (`AffBary::set` divise par den sans garde). Ils ne
  sont sûrs que parce que tous les appelants passent d'abord par `pairForce`, qui refuse den ≤ 1e-300. Ces
  fonctions étant publiques, ajouter la même garde.
- M7. P3 : le commentaire « A monte, B descend : on referme » (l. 11360) et la doc (« puis pressés ») ne
  décrivent pas la trace (`potx.csv`, mode 1). E décroît d'abord de 0,544 à 0, les blocs se détendent et le
  cliquet joue, puis le contact se réengage (Uc max 1,95e-3, soit environ 0,5 H0). Le test est
  discriminant ; c'est seulement la description qui est inexacte.
- M8. Coût sur le chemin par défaut : `pairEnergy` (9 clips) est désormais calculé à chaque naissance de
  paire dans tous les modes, pour le compteur. Aucun effet sur les forces ni sur l'ordre flottant, coût
  négligeable. À signaler seulement parce que la doc parle d'un « compteur pur ».

## Vérifié sans défaut

1. **Chemin par défaut** (clés au défaut) : forces inchangées.
   - `ePair` ne sert qu'au compteur, et la branche `else if (birthPenalty_)` / `else` rampe est identique.
   - `if (act_.empty() && !candVertex_)` équivaut à l'ancien test.
   - La boucle des paires 2D est séquentielle (pas d'OpenMP) et triée ; aucun conteneur non ordonné n'est
     itéré ; `potFt_` est seulement agrandi (champ `e0`).
   - Les réductions OpenMP `bias` et `lw` d'`integrate()` gardent le même ordonnancement statique ; seule
     leur valeur sur ROLLERX change, et elle n'agit que via C3.
2. **Signe et forme du terme −p·I_A·∇λ_a** :
   - P1 compare à des différences finies de `pairEnergy`, qui est indépendant de `pairForce`
     (2,2e-10 contre 0,157) ;
   - l'intégration par régions argmin est exacte (φ linéaire, centroïde) ;
   - repère local pour I_A et I_B ;
   - un recouvrement vide rend 0 sans correction ;
   - un élément inversé est écarté en amont.
   - Essai complémentaire : à |X| ~ 100 et h = 1e-4 ou 1e-6, la somme des forces nodales après correction
     reste au niveau de Munjiza (5e-13 relatif). L'absence de repère local dans `AffBary` de
     `exactGradientForces` ne dégrade donc pas la 3e loi.
3. **Offset** :
   - U(e0) = 0 et U′(e0) = 0, donc force continue ;
   - le cliquet est monotone et n'agit que sur la branche U = 0 ;
   - la règle de dérivation −U′(E)·∇E est cohérente avec les forces exactes ;
   - pas de course de données (2D séquentiel, 3D non porté).
4. **`vertex`** :
   - table `vElems_` construite sur `vOf_`, donc couvrant les copies dédoublées ;
   - tous les nœuds ont un `elemOf_` (seul site de remplissage, l. 2521-2523) ;
   - joints insérés en adaptatif pris en compte au pas suivant (balayage de `jt_`) ;
   - exclusion des paires à joint vivant inchangée, et `jt_` complet dès l'init (toutes les arêtes
     intérieures) ;
   - 3D : 4 sommets par tet, 3 nœuds par face.
5. **Rouleaux** :
   - les deux chemins 2D (groupé, par nœud) sont traités ;
   - ROLLERX est le seul drapeau qui annule une vitesse APRÈS le terme de biais : FIXED, DRIVEX et
     PRESCRIBED font `continue` avant ;
   - ressort `kAbsX` et Cundall : contributions nulles puisque v_x = 0 ;
   - en 3D, `uMask_` exclut déjà les axes imposés (l. 6989-6994, 7102).
6. **Validation des clés 2D** :
   - valeurs invalides refusées avec un message explicite ;
   - `offset`, `vertex` et `potForceExact` refusés hors `contact = potential` ;
   - `offset` + `gcBirthTau` refusé ;
   - `potForceExact` dans un deck fdem3d refusé par KeyGuard (lecteur fdem seul).
