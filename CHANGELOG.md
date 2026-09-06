# Changelog de rockim (arbre g0)

Format : [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/). Une section par tag (mesure A3/A5 du
plan de robustesse du 2026-09-05). L'arbre `g0` est sous git depuis le tag `g0-0.1.0` (mesure A1) ; la section `[Non publié]`
reçoit les lignes exigées par les règles déjà en vigueur — dont **toute ancre de bit-identité changée**
(`tools/bitid_refs.json`, règle de `tools/BITID.md`).

## [Non publié]


### Ajouté — `law = dfhplus` : le périmètre de DP-DFH sur une énergie libre postulée (2026-09-06)

Étape 1 du chantier **DFH+**, suite du banc `thermobench`. **`DpDfhLaw` n'est pas touchée** — ni son
code, ni ses clés, ni ses résultats (`selftest-dpdfh` rejoué OK, `thermobench dpdfh` reproduit
7 172 / 59 483, pire −15 240 J/m³, 1 124 significatives, 0 / 3 574 en symétrie) : `dfhplus` est une
**loi séparée**, dans un fichier séparé. Compte rendu chiffré : **`docs/DFHPLUS_etape1.md`**.

- **Nouveaux fichiers** `include/rockim/MatLawDfhPlus.hpp`, `src/MatLawDfhPlus.cpp` (la loi ;
  l'en-tête porte l'énergie libre, ses trois dérivées et le choix de la contrainte effective) et
  `src/DfhPlusBench.cpp` (`rockim selftest-dfhplus`, sept bancs). Une branche
  `else if (kind == "dfhplus")` dans `MatLaw::make`, trois entrées dans `CMakeLists.txt`, une
  commande dans `main.cpp`.
- **L'énergie libre**, avec **décomposition spectrale** de ε^e (Miehe et al. 2010) — seule la partie
  **positive** est dégradée, l'unilatéralité est naturelle :
  `ρψ = (λ/2)[g_v <tr ε^e>₊² + <tr ε^e>₋²] + G A:(ε^e⁺)² + G ‖ε^e⁻‖²`, avec `A = Σ_i (1−D_i) n_i⊗n_i`
  (repère figé) et `g_v` = moyenne **harmonique** des trois intégrités (couplage en série : g_v → 0
  dès qu'un D_i → 1). À D = 0 c'est **exactement** l'élasticité linéaire.
- **σ = ∂ρψ/∂ε^e** analytique (projection spectrale de Daleckii–Krein), vérifiée contre les
  différences finies sur ρψ : **2,0e−9**. **Y_i = −∂ρψ/∂D_i** analytique, **≥ 0 inconditionnellement**
  (deux carrés × λ > 0, G > 0) et **borné** à saturation (g_v ≤ 3(1−D_i) ⇒ le préfacteur ≤ 9).
  L'obscuration est branchée sur `Y_i` par la contrainte équivalente d'énergie
  `σ^eq = √(2 E Y_i / c_ν)`, calibrée pour valoir exactement σ̄_ii en traction uniaxiale à D = 0.
- **Plasticité sur `∂ρψ/∂ε|_{D=0} = C : ε^e`** — la seule mesure indépendante de D, donc la seule qui
  laisse la surface de charge en place quand la roche se fissure ; linéaire, donc le critère, le
  retour radial et l'apex de `dpdfh` s'y transposent sans modification. **Mesuré : compression
  uniaxiale et triaxiaux 20/50/100 MPa identiques à `dpdfh` à 2,4–2,6e−10.**
- **Clés.** Les neuf clés matériau de `dpdfh` sont reprises **à l'identique** (`dfhBetaDeg`,
  `dfhDCoh`, `dfhPsiDeg`, `dfhWeibullM`, `dfhSigW`, `dfhZeff`, `dfhK`, `dfhS`, `dfhDeld`, plus
  `dfhPsiVar`/`dfhPsi0`/`dfhKPsi`/`dfhPsiMax`) : **les cartes sont interchangeables**. Deux clés
  propres : `dfhpPsiClamp` (défaut `true`, écrêtage d'admissibilité de la dilatance avec
  prédicteur-correcteur) et `dfhpVolInteg` (`harmonic` | `min` | `none`). Registre des clés régénéré
  (`tools/gen_keys_by_mode.py`).
- **RÉSULTAT CENTRAL** — `rockim thermobench dfhplus` → **PASS, 0 violation** sur les cinq tests, sur
  les **mêmes** 12 000 tirages où `dpdfh` échoue : symétrie **0 / 2 078**, dissipation **0 / 120 000**
  (pire déficit **−4,7e−3 J/m³ = 3,8e−6 % de Gf/ℓc**, contre −15 240 J/m³ = 17,6 % pour `dpdfh`),
  réduction 0 / 96 000 (3,3e−15), objectivité 0 / 480 000 (2,3e−14), continuité 0 / 120 000.
  **Avec le MÊME instrument** que `dpdfh` (`--probe`, sonde de décharge) : **1 / 79 826** contre
  7 172 / 59 483, et **0 violation significative** contre 1 124.
- **Exposant de vitesse 3/(m+3) CONSERVÉ** : m = 6 / 12 / 24 → **0,3333 / 0,1989 / 0,1107** contre
  les cibles 0,3333 / 0,2000 / 0,1111 (−0,02 / −0,55 / −0,41 %), aussi près que `dpdfh` sur le même
  banc. La régression n'est faite que dans le **régime de fragmentation multiple** (λ_frag V_el > 10) :
  en dessous l'élément casse sur un défaut et le pic est rate-indépendant.
- **Écarts mesurés avec `dpdfh`** (banc 6) : pic de traction **identique à 5e−13** (il est atteint à
  l'amorçage, où les deux lois sont la même élasticité) ; aire sous la courbe de traction **+0,23 %** ;
  branche adoucissante **plus raide de +17 % au maximum** (vers D ≈ 0,65 : le terme λ ne peut pas être
  dégradé direction par direction sans casser la symétrie majeure) puis **plus relâchée** en fin de
  course (σ à D = 0,98 : 0,027 contre 1,094 MPa) ; cisaillement dégradé par la moyenne **arithmétique**
  des intégrités (conséquence du potentiel) au lieu du `min(f_i, f_j)` **postulé** par la VUMAT ;
  chemin **non coaxial** : +12,5 % de travail total, et **0 violation** là où `dpdfh` concentrait ses
  cinq pires déficits.
- **Contrôle qui DOIT échouer** : `dfhpPsiClamp = false`. Sur la carte de la thèse (Ψ = 15°)
  l'écrêtage ne mord jamais et la dissipation plastique est positive sur les 371 incréments **avec ou
  sans** la clé — c'est un résultat, pas un échec (sur la surface q̄/p̄ > tan β = 1,266 ≫ tan Ψ = 0,268).
  En écoulement **associé** (Ψ = β = 51,7°) la clé sert et le facteur est **21** sur le pire déficit
  (−8 127 → −386 J/m³). Elle ne suffit pas : les 252 incréments résiduels ont `q_proj < 0`, c'est le
  terme **déviatorique** qui est négatif — point ouvert, `docs/DFHPLUS_etape1.md` §9.

### Ajouté — `thermobench` : énergie libre exposée, et deux discriminants (2026-09-06)

- **Trois virtuelles AJOUTÉES à `MatLaw`** — `hasFreeEnergy()`, `freeEnergy(eps, state)`,
  `damageForces(eps, state, Y[])` — avec un défaut « non exposée » : **zéro effet sur les lois
  existantes** (croissance par addition, principe VIII ; `thermobench dpdfh` et
  `thermobench elastic` rendent exactement les mêmes chiffres qu'avant). Quand une loi les expose,
  le test 2 bascule de l'estimateur par sonde à la **valeur exacte** : les limites L1/L2/L3
  disparaissent et **plus aucun incrément n'est exclu** (0 contaminé, 0 non relâché sur `dfhplus`,
  contre 49 650 + 10 867 sur `dpdfh`). Nouvelle option **`--probe`** : force la sonde même sur une
  loi qui expose son énergie libre, pour comparer deux lois **avec le même instrument**.
- **Discriminant de QUADRATURE (test 2).** Un incrément flagué à 8 sous-pas est rejoué à **32** :
  l'erreur du trapèze au *coin* d'une réponse C¹ par morceaux est divisée par 16, une vraie violation
  ne bouge pas. Sur `dfhplus` : **58 flags à 8 sous-pas, tous effacés à 32** (rapport grossier/fin
  65 ≈ 8², la signature exacte de la quadrature).
- **Discriminant TEMPS / DÉFORMATION (test 5), et une lecture corrigée.** `r = ‖dσ‖/(M‖dε‖)` n'est un
  critère de continuité de σ(ε) qu'**à temps gelé** : un mécanisme piloté par le temps (obscuration,
  dx ~ dt) fait tomber σ sans que ε bouge, et le raffinement **ne peut pas** l'effacer puisqu'il
  divise dt en même temps que dε. Le verdict porte désormais sur la valeur à temps gelé (dt × 1e−6),
  la valeur à temps courant reste publiée. **Conséquence : les 84 flags de `dpdfh` et les 87 de
  `dfhplus` disparaissent tous**, et le pire r tombe à **0,99999 ≤ 1** — exactement la borne
  élastique — pour les deux lois. La conclusion du compte rendu précédent (« discontinuité de la loi
  à saturation ») était **fausse** : les deux réponses σ(ε) sont continues.

### Ajouté — `rockim thermobench` : banc thermodynamique générique des lois (2026-09-06)

Étape 1 du chantier **DFH+**. La relecture adverse de
`CONTINUUM/loi_dfh_plus/loi_DFH_plus.pdf` (66 p., 146 équations) a établi que le cadre de DP-DFH est
**sur-déterminé** — contrainte effective à la Lemaitre *plus* équivalence en énergie *plus* énergie
libre écrite indépendamment : trois énoncés pour deux libertés — et que plusieurs dissipations peuvent
y devenir négatives. Décision de Fernando : reconstruire la loi sur une **énergie libre postulée
unique** (σ = ∂ρψ/∂ε^e, Y_k = −∂ρψ/∂D_k, décomposition spectrale de la déformation). Un tel cadre se
**prouve** au point matériel — donc le banc s'écrit **avant** la loi, comme l'exige la règle « banc
court avant run long ».

- **Nouveaux fichiers** `include/rockim/ThermoBench.hpp` et `src/ThermoBench.cpp`, câblés dans
  `main.cpp` à côté des `selftest-*` et ajoutés à `CMakeLists.txt`. **Aucune loi n'est touchée** : le
  banc n'utilise que `MatLaw::stress(eps, MatState&, dt, lc)`. `dpdfh` est inchangé — code, clés et
  résultats.
- **Cinq tests** sur des états et chemins **tirés au hasard** (graine fixe, 12 000 tirages par défaut,
  six familles : traction, compression, **triaxial 0–300 MPa**, cisaillement, chemins **non
  coaxiaux**, états **fortement endommagés** jusqu'à D = 0,9999 ; dt 1e-9 – 1e-5 s, ℓ_c 0,5 – 2 mm,
  repère et position aléatoires) : (1) **symétrie majeure** de la tangente numérique centrée, avec
  séparation des incréments élastiques (verdict) et inélastiques (information — un écoulement non
  associé est légitimement asymétrique) ; (2) **positivité de la dissipation**, ρψ estimée par sonde
  de **décharge élastique à temps gelé**, avec détection et exclusion des incréments contaminés ou non
  relâchés ; (3) **réduction** — `--ref <loi>` compare deux lois sur les mêmes chemins à 1e-12, sans
  `--ref` c'est la réduction élastique ; (4) **objectivité** sous rotation rigide superposée ;
  (5) **continuité** — ‖dσ‖ / (M‖dε‖) sur l'incrément puis découpé en 8. Détail, tolérances et
  **limites explicites de l'estimateur d'énergie libre** : en-tête de `src/ThermoBench.cpp` et
  DOCUMENTATION §3.4.
- **Déterminisme** : un générateur par tirage → CSV **identique à 1, 4 et 8 fils** (vérifié).
  0,7 s pour 12 000 tirages à `OMP_NUM_THREADS = 4`.
- **Les deux contrôles du banc.** `thermobench elastic` → **PASS**, 0 violation, tout au bruit machine
  (réduction 0,0 exactement, continuité r ≤ 0,999998). `thermobench dpdfh` → **ÉCHEC**, code 1 :
  **7 172 / 59 483** violations de dissipation (12,1 %), dont **1 124 dépassent 1 % de G_f/ℓ_c**, pire
  cas **−15,2 kJ/m³ = 17,6 % de G_f/ℓ_c** sur un chemin **non coaxial** à D = 0,9999 ; et **84 /
  120 000** violations de continuité, pire cas un saut de **109,3 MPa** en un incrément, **inchangé au
  raffinement ×8**. Objectivité et réduction élastique **exactes** (5e-14, 8e-16) ; tangente
  **élastique symétrique** (2,9e-11) — le défaut est dans le couplage des mécanismes, pas dans
  l'élasticité endommagée. Contrôle du test 3 : `--ref dpdfh` sur `dpdfh` → 0 / 480 000, écart 0,0
  exactement ; `--ref elastic` → 125 567 / 480 000. Relevé des sept lois du dépôt : DOCUMENTATION §8,
  point 11.
- **Bit-identité** : `python tools/bitid.py --exe build/rockim.exe --threads 4` → **8/8 IDENTIQUE**
  (le banc n'ajoute qu'une commande ; ancre inchangée).

### Ajouté — port de la branche orpheline `insertion-pointe` (décision 3)

Quatre commits de `insertion-pointe` (2ead636), jamais fusionnés, dont **aucune des six clés
n'existait dans `f2`** : contenu réellement orphelin. Deux decks de `f2`
(`bench_impact/configs/impact3d_dpdfh.cfg` et `..._gros.cfg`) les posaient déjà et étaient donc
**refusés** par `rockim_f2w21` (garde C1). Tout arrive en **clé opt-in à défaut bit-identique**.
Détail complet, mesures et bancs : [`docs/PORT_INSERTION_POINTE.md`](docs/PORT_INSERTION_POINTE.md).

- **`insertionTipFactor` (défaut 1) et `insertionTipDamage` (défaut 0,5)** — insertion préférentielle
  en POINTE de fissure, 2D et 3D. Une facette dont un sommet porte déjà un joint inséré et endommagé
  (`D >= insertionTipDamage`) voit son enveloppe DIVISÉE par `insertionTipFactor` ; une facette en
  terrain vierge garde l'enveloppe nominale, donc **l'amorçage — et la résistance macroscopique
  mesurée — restent intacts** (mesuré : pic 51,0807 → 51,0807 MPa). Motivation : l'adaptatif ne
  propage que 43,7 % de ses ruptures contre 58,9 % pour l'intrinsèque (chiffre CORRIGÉ le 2026-09-06 : « 56,8 % » était une coquille de recopie, la mesure d'origine est 58,9 % — `BILAN_insertion_adaptative.md` l. 102, 127 et 154, et `PHD.md`), la moyenne sur deux CST
  écrasant la singularité de pointe. À `insertionTipFactor = 1` aucun test supplémentaire n'est
  évalué : **chemin d'origine au bit près**. Gardes : facteur < 1 refusé, `insertionTipDamage` hors
  [0,1] refusé, facteur > 1 sans `insertion = adaptive` refusé. Compte
  `propagations / nucléations` imprimé en fin de run.
- **`insertion = none`** — le CONTINUUM PUR (2D et 3D) : aucun joint n'existe ni ne peut naître, les
  copies de nœuds restent liées pour toujours (= éléments finis à nœuds partagés). Remplace le
  bricolage `insertion = adaptive` + `ft = 1e12`, qui a fait passer un impact 3D DP-DFH de 53 J à
  **−89 GJ en 10 µs** le 2026-08-25 (89 424 joints « inatteignables » activés par un élément
  distordu, portant `dnE = ft/pj = 8 cm`). Trois conséquences câblées : joints liés comme en
  adaptatif, groupes liés intégrant comme UN nœud, et **ressort de pénalité des joints sorti du
  budget de pas de temps** (mesuré : dt × 1,715 sur le banc 2D).
- **`dfhPsiVar` (défaut 0), `dfhPsi0`, `dfhKPsi`, `dfhPsiMax`** — dilatance VARIABLE ψ(p̄) du DP-DFH,
  la forme de `vumat_hole.f` l. 336-339 : `ψ = clamp(dfhPsi0 − dfhKPsi·p̄[MPa], 0, dfhPsiMax)`. La
  constante #5 de la carte (`dfhPsiDeg`, 15 deg) est alors morte. Clé absente ⇒ `psiDeg` fixe,
  trajectoires inchangées au bit près.
- **`dfhD` et `dfhTini` dans les `.vtu`** sous `law = dpdfh`, **en 2D seulement** : DMAX = max(D1,D2,D3) des
  trois endommagements directionnels du repère figé (SDV 4-6 de la VUMAT, ce que lisent les
  extracteurs du banc 6) et l'instant du premier amorçage. **Ajout de sortie pur.**
  RÉSERVE (relecture adverse du 2026-09-06) : la branche d'origine écrivait ces deux champs
  dans les DEUX solveurs (`git show insertion-pointe:src/Fdem3dSolver.cpp` l. 3421) ; le port
  ne les a mis qu'en 2D, le `writeFrame()` de `Fdem3dSolver` reste muet sur l'endommagement
  DP-DFH. À porter en miroir (une quinzaine de lignes) avec une passe `bitid` de confirmation —
  aucun deck de l'ancre n'étant en `dpdfh`, l'ajout y serait neutre.
- **Decks** `tunnel_edz/configs/tunnel_tip13.cfg`, `tunnel_tip16.cfg`, `tunnel_tip20.cfg` (calibration
  du facteur 1,3 / 1,6 / 2,0 sur le tunnel EDZ), repris tels quels.
- **Bancs courts** `tests_f2/insertion_pointe/` (4 decks + `check_tip.py`, ~10 min à OMP 2) et
  `tests_f2/psivar/` (2 decks point-matériel + `check_psivar.py`, ~2 s). Tous deux impriment un
  verdict OK/ÉCHEC par critère. `check_psivar.py --falsify` (écrite à la relecture adverse du
  2026-09-06 : l'option était annoncée dans l'en-tête du script et n'existait pas) rejoue le
  deck d'essai avec `dfhPsiVar = 0` : écart max **+0,000e+00**, trace **bit-identique** au
  témoin — le critère B échoue comme il doit, la clé est bien le seul pilote de ψ(p).
- **`etude_lois_fem/meshes/drop_orphans.py` et `check_orphans.py`** : ces outils sont NOMMÉS par le
  message d'erreur du lecteur de maillage mais n'avaient pas été repris à la naissance de `g0`
  (seuls les `*.py` de la racine de `etude_lois_fem/` l'avaient été). Repris de `rockim_f2`.

### Corrigé
- `bench_impact/configs/impact3d_dpdfh.cfg` et `impact3d_dpdfh_gros.cfg` **chargent** (ils étaient
  refusés clé par clé). Trois réparations : les six clés existent maintenant ; la **commande exacte
  de régénération** du maillage est dans l'en-tête (les `.msh` n'existaient nulle part, ni dans `g0`
  ni dans `rockim_f2`, et les comptes annoncés — 98 858 et 143 451 tétraèdres — venaient de
  paramètres jamais consignés : les comptes réellement obtenus, 105 498 et 150 535, sont écrits) ;
  les decks pointent sur `*_clean.msh` car gmsh écrit le point du champ de taille (= le point
  d'impact) en nœud 0-D orphelin que le lecteur refuse.
- Message d'annonce du bloc de liaison des joints : il disait `adaptive insertion: N bonded edges`
  même sous `insertion = none`. Il dit maintenant `insertion = none: N bonded edges` dans ce cas.
  **Sortie seule**, aucun flottant touché.

### Ajouté — port de la branche orpheline `dif-intrinseque` (décision 3)
- `bench_impact/tools/fig_bilan.py` (159 lignes) : la planche du run de référence — partition de
  l'énergie contre les chiffres publiés (ARMA 24-0952 : fissuration 2,6 %, frottement 64,9 % de
  49,3 J) et les sept critères de la Table 3 de Yang et al. rapportés à leur fourchette. C'est,
  avec la coquille de `fig_fp.py`, **tout** ce que `dif-intrinseque` apportait encore : son C++ est
  un sous-ensemble strict de `joint-handoff` (`git diff --stat joint-handoff...dif-intrinseque --
  src include` est vide).
- **Corrigé** `bench_impact/tools/fig_fp.py` : coquille de chaîne non brute `label="moyenne (30 $\mu$s)"` →
  `label=r"..."` (portée de `dif-intrinseque`).

### Non porté — branche `joint-handoff` (décision 3), lecture ligne par ligne
- **Aucune ligne de `src/` ni de `include/` n'est reprise.** Les trois « pertes probables » de la
  revue sont **toutes les trois déjà dans `g0`**, et dans les trois cas ÉTENDUES par `f2` : le poste
  d'énergie séparé `brushWork_` (que `f2` peut en plus faire entrer dans le bilan par
  `energyBodyForces`) ; le court-circuit `if (muCRes_ < 0.0) return muC_;` de `contactResidualMu`
  (devenu `return mu` — le rendre à `muC_` **annulerait** le frottement par phase de `f2`, ce serait
  une régression) ; le budget de pas de temps Signorini A1 qui sort `kp_` du CFL (présent dans les
  deux solveurs, `f2` ajoutant les platines, le contact général SHPB, le potentiel,
  `dtBudgetTangential` et deux bornes visqueuses). Sur les 106 lignes de `joint-handoff` absentes de
  `g0`, 98 sont des lignes de `main` que `f2` a réécrites pour son propre compte et 8 seulement ont
  été écrites par la branche — toutes du refactor. Verdict ligne par ligne :
  [`docs/PORT_JOINT_HANDOFF.md`](docs/PORT_JOINT_HANDOFF.md).

### Bit-identité
- `python tools/bitid.py --exe build/rockim.exe --threads 4` contre l'ancre de naissance
  `tools/bitid_refs.json` : **8/8 IDENTIQUE** après le port (`results/bitid_apres_port_insertion_pointe.json`).
  L'ancre est **inchangée** (aucun `--update`).
- **Trou de l'ancre comblé à la main** : aucun de ses 8 decks ne pose `law = dpdfh`, elle ne pouvait
  donc rien dire du port de ψ(p). Preuve dédiée : sur le point matériel DP-DFH
  `tests_f2/psivar/mp_dpdfh_psivar_off.cfg`, `build/rockim.exe` et le binaire d'avant le portage
  `rockim_f2w21.exe` rendent la **même trace au SHA-256 près**
  (`341fac848845623c2db4a702897f4efab6a797842d843964b6cf83274c997580`).

### Relecture adverse (2026-09-06)

Contre-expertise indépendante de la naissance et des trois ports. **Verdict : conforme.** Ce qui a
été vérifié à la main, et ce qui a été corrigé.

Vérifié :
- `rockim_f2` **intact** : aucun fichier de `src/`, `include/` ni aucun `*.exe` postérieur à
  20 h 54 (dernière écriture de Fernando, `src/main.cpp` ; dernier exe `rockim_f2w22.exe` à 20 h 51),
  soit plus d'une heure avant la création du worktree `g0` (22 h 01). Les seuls fichiers de
  `rockim_f2` touchés depuis sont des SORTIES des campagnes en cours
  (`etude_lois_fem/heterogeneite/figures/`, `etude_lois_fem/bitid_w21/`). Les six worktrees
  (`rockim_p1`, `rockim_f2_wt`, `rockim_p2`…`p4`, `studio_wt`) sont dans l'état qu'ils avaient déjà.
- `rockim_g0` est bien un worktree de `FDEM/rockim` (gitdir `rockim/rockim_p1/.git/worktrees/rockim_g0`)
  sur la branche `g0`, `git status --porcelain` **vide** après les trois commits de port.
  Aucun `*.exe`, `*.obj`, `*.vtu`, `*.pdb` n'est suivi ; les `*.msh` suivis sont 14 fichiers,
  11,4 Mo au total, le plus gros 4,67 Mo (les 4 exemptions déclarées, plus 10 fichiers déjà suivis
  par `f2-2026-09-02` — conservés au titre de la décision 2). Rien n'est poussé sur `origin`.
- **Rebuild complet** `tools/build.ps1 -Clean` (reconfiguration CMake + 13 objets + link) puis
  `python tools/bitid.py --exe build/rockim.exe --threads 4` : **8/8 IDENTIQUE** contre l'ancre
  héritée, inchangée (rapport `results/bitid_relecture_adverse_2026-09-06.json`). L'exe reconstruit
  fait 2 068 480 octets comme le précédent mais un SHA-256 différent : MSVC horodate le PE, c'est
  l'algèbre qui est reproductible, pas l'octet du binaire.
- **`diff -r` de `src/` et `include/` entre `rockim_f2` et `rockim_g0`** : exactement 6 fichiers
  diffèrent, et ce sont les 6 du port — `FdemSolver.cpp/.hpp`, `Fdem3dSolver.cpp/.hpp`, `MatLaw.cpp`,
  `KeysByMode.hpp` (régénéré). Aucune autre divergence.
- **Opt-in relu dans le code, pas dans le compte rendu** : `tipFactor_ = 1.0` et `noJoints_ = false`
  par défaut ; sous `tipFactor_ == 1` le drapeau `tipBias` est faux, `vertTip_` n'est jamais rempli
  et le facteur vaut `1.0` exactement (multiplication neutre au bit près) ; `psiVar` défaut `false`.
  Le registre `KeysByMode.hpp` se régénère à l'identique (hors horodatage) par
  `tools/gen_keys_by_mode.py`.
- **`docs/PORT_JOINT_HANDOFF.md` recompté indépendamment** et dix de ses verdicts recontrôlés dans
  le code : voir le paragraphe « 4bis » ajouté à ce document. La conclusion « rien à porter » tient.
- **Banc ψ(p) rejoué** avec l'exe reconstruit : les cinq lignes du tableau de
  `docs/PORT_INSERTION_POINTE.md` §4.2 sont reproduites au chiffre près.

### Corrigé — relecture adverse (2026-09-06)
- `tests_f2/psivar/check_psivar.py` : l'option **`--falsify`** était annoncée dans l'en-tête du
  script ET dans `mp_dpdfh_psivar_on.cfg` mais **n'existait pas** (`unrecognized arguments`). Elle est
  écrite : elle rejoue le deck d'essai avec `dfhPsiVar = 0` et exige que la trace redevienne
  bit-identique au témoin. Mesuré : écart max `+0,000e+00`, `sha256` identique — le critère B échoue
  comme il doit. Un banc dont la variante négative n'existe pas ne prouve rien : c'est exactement la
  règle « un réglage que l'on peut écrire, que le solveur accepte, et qui ne fait rien ».
- `docs/PORT_INSERTION_POINTE.md` §4.1 citait un rapport `results/bitid_g0-0.2.0_insertion_pointe.json`
  **qui n'existe pas** (le fichier s'appelle `bitid_apres_port_insertion_pointe.json`, et il n'y a pas
  de tag `g0-0.2.0`). Référence corrigée.
- `dfhD` / `dfhTini` étaient annoncés « 2D et 3D » ici et dans `docs/PORT_INSERTION_POINTE.md` §2.4 ;
  ils ne sont **qu'en 2D**. La branche d'origine les écrivait dans les deux solveurs : le port est
  incomplet sur ce point. Titre et texte corrigés, réserve écrite aux deux endroits.
- `docs/PORT_JOINT_HANDOFF.md` disait que « la routine `nodeSig` a été sortie vers
  `ToolSignorini.hpp` » : seule l'ALGÈBRE impulsion / saut de vitesse (`toolsig::impulse`) l'a été,
  la lambda `nodeSig` et toute la géométrie restent dans chaque solveur, à dessein. Rédaction
  corrigée aux §3 et §4 ; le verdict REFACTOR, lui, est confirmé.
- `CHANGELOG.md` : l'historique hérité de `rockim_f2` suivait la section `[g0-0.1.0]` **sans titre de
  version** et se lisait donc comme en faisant partie. Un titre
  `## [hérité de rockim_f2]` l'en sépare. Préambule remis à jour (`g0` est sous git depuis A1).
## [g0-0.1.0] — 2026-09-05 — naissance de g0

### Naissance
- **`g0` naît de `rockim_f2`.** Branche `g0` créée dans le dépôt
  `FDEM/rockim` à partir de `f2-2026-09-02` (commit `4d65395`, la branche la plus proche de l'arbre
  vivant `rockim_f2`), worktree `FDEM/rockim_g0`. L'**état complet** de `rockim_f2` au 2026-09-05 22 h
  (sources `src/` + `include/`, `tools/`, `tests_f2/`, `configs/`, `configs_yan/`, les documents de
  travail dont `PLAN_ROBUSTESSE_2026-09-05.md` et `ETAT_DES_LIEUX_2026-09-05.md`, les `*.py`/`*.md` de
  `etude_lois_fem/`) y est recopié. `diff -r` de `src/` et de `include/` entre `rockim_f2` et
  `rockim_g0` : **0 différence**.
- **Rien n'est supprimé** (décision 2 du plan) : `rockim_f2` n'est ni modifié ni recompilé, les
  worktrees `rockim_p1`…`rockim_p4` restent en place. Ne sont pas repris dans g0 les produits
  (`*.exe`, `*.obj`, `out_*/`, `*.vtu`, gros CSV, figures) ni les dossiers de campagne
  (`etude_lois_fem/matrice*`, `heterogeneite`, `bitid_w*`, `phaseA/B`, `bench_impact/out*`,
  `calib_quick`, `panel_contact*`, `tunnel_schisto`, `g2d`/`g3b`/`gif*`) : ils restent dans
  `rockim_f2`, gelé mais conservé. `.gitignore` réécrit en conséquence.
- **Maillages de la bit-identité versionnés** : les quatre `.msh` dont dépendent les decks
  `tests_f2/bitid/*.cfg` (`meshes/cut3d_h50.msh`, `meshes/impact_kuru_s15.msh`,
  `meshes/t1_toolcontact.msh`, `etude_lois_fem/meshes/T1_c05_clean.msh` — 7,5 Mo au total) sont
  exemptés de la règle `*.msh` du `.gitignore`, sans quoi `tools/bitid.py` ne serait pas rejouable
  depuis un clone. Le cinquième (`Q1_c05_clean.msh`, deck `fem3d_cdp_PQ_court`) reste **hors dépôt**,
  atteint par le chemin relatif inchangé `../../CONTINUUM/calib_bohus_triax/cdp_rockim/perc3d/` :
  `rockim_g0` étant, comme `rockim_f2`, un dossier de `FDEM/`, le chemin des decks n'a pas eu à être
  touché (aucun deck de bit-identité n'a été modifié).

### Ajouté
- **CMake + Ninja** (décision 9, mesure D1a) : `CMakeLists.txt` réécrit pour émettre **exactement** les
  drapeaux de `build_f2.cmd` — `FLAGS = /nologo /EHsc /O2 -std:c++17 -MT -openmp` plus
  `-D_USE_MATH_DEFINES -DNOMINMAX -I include -I ../rockim/eigen-3.4.0` — c'est-à-dire en neutralisant
  les ajouts par défaut de CMake sous MSVC : `/DWIN32 /D_WINDOWS /W3 /GR`, `/Ob2 /DNDEBUG` (qui
  désarmerait les `assert()`) et surtout `/MD` (CRT dynamique) là où `cl` sans drapeau prend `/MT`
  (vérifié : `rockim_f2w21.exe` n'importe que `KERNEL32.dll` et `VCOMP140.DLL`, l'exe CMake aussi).
  `/W4 /WX` viendront au lot D1b, séparément. `tools/build.ps1` enchaîne `vcvars64` → `cmake -G Ninja`
  → `cmake --build`, avec les `cmake.exe`/`ninja.exe` livrés par Visual Studio 2022.
- `build_f2.cmd` est **conservé** un temps comme secours (décision F du plan), il n'est plus la chaîne
  de référence.

### Bit-identité
- Ancre **inchangée** (`tools/bitid_refs.json`, héritée de `rockim_f2`, prise sur `rockim_f2w18.exe`,
  `OMP_NUM_THREADS = 4`, 8 decks). `python tools/bitid.py --exe build/rockim.exe --threads 4` sur
  l'exe produit par CMake + Ninja : **8/8 IDENTIQUE** (2026-09-05 22 h 43 ; exe CMake sha256 `92ef42b6e8168322...` ; rapport
  `results/bitid_g0-0.1.0.json`, journal `results/bitid_g0-0.1.0.log.txt`). La chaîne CMake reproduit donc l'algèbre de la
  chaîne `build_f2.cmd` : g0 naît sur l'ancre de f2, et cette ancre est **la nouvelle ancre de
  naissance** de l'arbre g0 (tag `g0-0.1.0`).


## [hérité de rockim_f2] — historique antérieur à la naissance de g0

Sections reprises telles quelles du `CHANGELOG.md` de `rockim_f2` : elles décrivent l'arbre
`rockim_f2`, gelé mais conservé, et non des changements faits dans `g0`. Sans ce titre elles
se lisaient comme faisant partie de `[g0-0.1.0]` (relecture adverse du 2026-09-06).

### Ajouté
- 2026-09-05 — `tools/bitid.py` (bit-identité générique, mesure B6), `tests_f2/bitid/` (8 decks
  courts : 3 fem3d, 3 fdem3d, 2 fdem), `tools/exe_manifest.py` + `tools/exe_manifest.json`
  (empreinte SHA-256 des 46 `rockim_f2*.exe` existants), `tools/BITID.md`.
- **ancre bitid changée : première ancre** — prise sur `rockim_f2w18.exe` existant (05/09 13:45,
  `build_f2.cmd`, SHA-256 `655a87329c68b0d88e499a6eebfa73543d5869d9373685ea6430d7c5f1e8c008`),
  `OMP_NUM_THREADS = 4`, 8 decks ; remplace les scripts `etude_lois_fem/bitid_wNN.sh` (ancre changée
  deux fois, hachages non comparables entre nombres de fils). Détail dans `tools/BITID.md`.
  Le dépôt git de `rockim_f2` n'existe pas encore (A1) : le json sera committé avec l'état du 05/09.
- 2026-09-05 (w21, C1/C7) — `include/rockim/KeyGuard.hpp` (audit des clés consommées après `init()` : obsolète,
  autre mode, famille dynamique, suggestion de Levenshtein ≤ 2, inconnue ; politique `unknownKeys = error | warn` ;
  écriture de `config_effective.cfg`) ; `Config::unusedKeys/lineOf/consumedKeys/effective/path/size` ;
  `tools/obsolete_keys.json` (table des clés obsolètes → nom courant, une entrée par renommage, jamais retirée) ;
  registre `tools/gen_keys_by_mode.py` étendu (`kKnown`, `kDynamicPrefix`, `kCommonPrefix`, `kObsolete` dans
  `KeysByMode.hpp` ; `prefixes_dynamic`, `obsolete` dans le json) ; `tools/scan_decks.py` (balayage statique des
  decks selon la règle, `--fix-obsolete`) + rapport `tools/scan_decks_2026-09-05.md` ; banc falsifiant
  `etude_lois_fem/bitid_w21/selftest_cles/` (17 decks, `run_cles.py`, `SELFTEST_cles.md`) ; DOC §8.11.
- **C7** : chaque run écrit `<outputDir>/config_effective.cfg` (clés consommées, valeur deck ou `(defaut)`, tri
  stable) et imprime `[rockim] cles : N consommees (n du deck, m au defaut), k du deck non lues`.
- 2026-09-05 (w22) — `KeysByMode.hpp` : table `kReaders` (lecteurs de chaque clé : solveurs / `shared`), champ
  `readers` du json enrichi (`shared` pour `ALWAYS_COMMON`) ; `Config::deck()`, `Config::seal()` ; banc
  `selftest_cles` : decks `c12`–`c15` (clés à plusieurs lecteurs hors mode, `warn`) et partie C (rejeu de
  `config_effective.cfg`, `history.csv` bit-identique sur les six modes) ; `tools/scan_decks.py` modélise les gardes
  pré-init de `main.cpp` et met à part les fichiers sans clé de solveur (matpoint, cartes matériau).

### Corrigé
- 2026-09-05 (C0) — `fragBrushV` → `fragBrushV0` dans 7 decks de `bench_impact/configs/` (`impact_kuru9*.cfg` ×4,
  `impact_pulv_a/b.cfg`, `impact_stanne_fidele_s15.cfg`) et les 7 mêmes de `FDEM/rockim/rockim_f2_wt/` (originaux
  dans `etude_lois_fem/bitid_w21/orig/`). **Les runs Kuru / pulvérisation / St Anne antérieurs tournaient SANS
  brossage des fragments** (la clé était inerte depuis l'origine) : leurs résultats changent si on les rejoue.
  Restent à corriger (worktrees d'archive, décision Fernando) : `rockim_p1` et `rockim_p2` `bench_impact/configs/impact_pulv_a/b.cfg`, `impact_stanne_fidele_s15.cfg`.

### Modifié
- 2026-09-05 (w21) — **clé inconnue = erreur (décision Fernando 2026-09-05, 20:00)** : toute clé du deck non
  consommée pendant l'initialisation et hors registre (ni commune, ni du mode courant) rend le code 1 avec un message
  qui explique (obsolète → nouveau nom ; autre mode → le mode ; faute de frappe → suggestion ; sinon « inconnue de
  rockim »), toutes les clés fautives listées d'un coup ; `unknownKeys = warn` : mêmes messages en WARNING, le run
  continue. Seul changement de défaut du chantier C. Le contrôle par mode de w20 (`keysbymode::check` avant init) est
  **fondu** dans cet audit (même sous-chaîne de message ; `warn` s'y applique) ; `keysWithPrefix("hydro")` n'est plus
  appelé en mode fdem (il marquerait consommées les clés `hydro*`). Bit-identité : ancre **inchangée** (w18, 4 fils),
  `tools/bitid.py --exe rockim_f2w21.exe --threads 4` : voir `etude_lois_fem/bitid_w21/bitid_w21.log`.
- 2026-09-05 (w22, relecture adverse de w21 — D1, D2) — **règle des lecteurs** : une clé du deck non consommée est
  légitime ssi le mode courant ou le code partagé la lit (`kReaders`) ; sinon `cle 'X' sans effet en mode Y : cle du
  mode Z seulement` / `cles des modes Z1, Z2`. w21 tenait pour « commune = légitime » toute clé lue par ≥ 2 solveurs :
  `jointSoftening`, `insertion`, `bulkDamage`, `fragBrushV0`, `gravity`… passaient en fem3d sans message (D1). 0 deck
  existant touché (957 balayés). **`config_effective.cfg` rejouable** : défauts du code en lignes commentées
  `# cle = valeur (defaut)`, clés du deck non consommées mais légitimes gardées actives ; w21 activait les défauts et le
  rejeu échouait en fem3d (`tensionShearRetention`) et fdem (`dampingLocalAfter`) sur les gardes « satellite orpheline »
  (D2). **Coût** : `Config::seal()` après l'audit — plus de verrou ni de chaîne de défaut dans les getters appelés à
  chaque pas (`confineGaugeTime`). `src/main.cpp` repassé en CRLF (original). Bit-identité : ancre **inchangée** (w18,
  4 fils), `tools/bitid.py --exe <abs>/rockim_f2w22.exe --threads 4` : **8/8 IDENTIQUE** (21:15,
  `etude_lois_fem/bitid_w22/bitid_w22.log`) ; banc `selftest_cles` 47/47 (A 21, C 6, B 20).
- 2026-09-05 (relecture de B5/B6, ancre **inchangée** : mêmes 8 × hachages) — `tools/bitid.py` :
  sorties par défaut dans un dossier temporaire (`tempfile.mkdtemp`, convention de `verify_suite.py`)
  au lieu de `tests_f2/bitid/out_*` ; `.vtu` supprimés aussi sur échec/timeout ; détection Apex One
  par motif étroit (l'ancien « Acc » aurait pris `fragBrushAccel` pour un refus) ; chemins de decks en
  séparateur posix. `tools/bitid_refs.json` : métadonnées seulement (clé orpheline `partial` retirée,
  `cfg` en posix). `tests_f2/bitid/fdem3d_kuru9_court.cfg` : en-tête complété (quatre clés du deck
  source absentes, non documentées avant). `tests_f2/bitid/.gitignore` (`out_*/`).
- 2026-09-05 19:55 — ancre bitid partielle : `fem3d_cdp_PQ_court` rebasé sur `Q1_c05_clean.msh` (w20 refuse le maillage à nœud orphelin) ; history/frames inchangés (e43ca95062fbc58b / 6670a60ed53648e0), vtu e2e42e15c8fbdf2d ; ancre toujours rockim_f2w18.exe, threads 4. Passe complète `tools/bitid.py --exe rockim_f2w20.exe` : 8/8 IDENTIQUE (w20 = gardes C3/C4/C6, sans effet sur les decks valides).
