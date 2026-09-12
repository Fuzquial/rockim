# Changelog de rockim (arbre g1)

Format : [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/). Une section par tag (mesure A3/A5 du
plan de robustesse du 2026-09-05). L'arbre `g0` est sous git depuis le tag `g0-0.1.0` (mesure A1) ; la section `[Non publié]`
reçoit les lignes exigées par les règles déjà en vigueur — dont **toute ancre de bit-identité changée**
(`tools/bitid_refs.json`, règle de `tools/BITID.md`).

## [Non publié] — arbre g1, 2026-09-11 soir : le point sur les impacts Yang/Solidity

Point demandé par Fernando (« à chaque fois une nouvelle erreur, c'est infini »). Document :
`docs/ETAT_yang2026_2026-09-11.md`. Binaires `rockim_g1y.exe` puis **`rockim_g1y2.exe`** (même
physique, plus l'impression du nœud qui commande le pas). Bit-identité : `results/bitid_g1y2.log`.

### Ajouté (opt-in, défaut bit-identique)

- **`groupContinuum.<corps> = true`** (fdem3d, mesh = file) : corps sans joints — facettes
  intérieures `Joint::perm`, jamais évaluées ni insérées, hors budget CFL ; liaison de nœuds par
  groupes armée en schéma intrinsèque dès qu'il en existe (`integrate()` bascule sur la branche par
  groupes si `nPerm_ > 0`). Sur le deck Yang s = 1 : 39 803 facettes sur 232 408 retirées des
  ressorts, acier et carbure continuums EF exacts. Gain mesuré sur dt : ×1,06 seulement — le pas
  est repris aussitôt par les slivers de la **roche** (diamètre inscrit 0,226 mm, K = 98,6 % de
  joints), voir ETAT §5.1.
- `tools/make_impact_mesh.py` : arguments nommés `gap=`, `gapr=` (jeux piston/bit et
  insert/roche), `opt=N` (passes Netgen + Relocate3D), `algo3d=` (10 = HXT). Défauts inchangés.
  Mesuré : `opt=2` sans effet sur le sliver (0,227 mm), HXT **pire** (0,183 mm, dt 0,80 ns).
- `tools/yang_report.py` : les sept critères de Yang et la santé du run depuis `history.csv`.
- Decks : `configs/yang2026_impact.cfg` **v2** (onze écarts de la v1 corrigés, ETAT §3),
  `configs/yang2026_bench_s25.cfg` et `_court.cfg` (banc court s = 2,5, 10 563 tets, toute la
  chronologie en minutes). Maillages `meshes/impact_yang_s1_pose.msh` (120 185 tets, jeux 0,02 mm)
  et `impact_yang_s2.5_pose.msh`.
- Registre `tools/keys_by_mode.json` régénéré (préfixe `groupContinuum.`) ; ligne dans
  `DOCUMENTATION_rockim.md`.

### Avertissements nouveaux (impression seule)

- Outil analytique **fantôme** : `groupVel.<corps>` posé et `toolShape ≠ none`.
- **Résistances héritées** : `phase.<x>.E` posé sans `.ft` ni `.cohesion` quand E > 1,5 E_global
  (un carbure dont les joints cassent à 11 MPa).
- Translation du maillage et **étendue en z de chaque corps** dans le repère solveur ; jauge hors de
  la moitié centrale d'un corps.
- **Nœud qui commande le pas de temps** : corps, phase, h inscrit, part élément / joints / contact.

### Modifié — contact par potentiel 3D PARALLÈLE (`rockim_g1y3.exe`, commit séparé)

`Fdem3dSolver::potentialContact()` était « série et déterministe » : sur le deck Yang s = 1 la
grille + la boucle des paires pesaient 53 % du pas (14 fils inutilisés), et le banc s = 2,5 a
ralenti **×8** dès l'amorçage de la fracture (6 → 0,73 µs simulées par minute). Découpage en trois
phases, **bit-identique par construction** : A (série) exclusion des paires portées par un joint
vivant et entrée du cache par paire (`potFt_`, insertion non thread-safe, pointeurs de nœuds
stables) ; B (parallèle, `schedule(static)`) axe séparateur + clip du polyèdre — géométrie pure,
n'écrit que le `sepAxis` de sa paire, scratch `static thread_local` dans `PotentialContact.hpp` ;
C (série, ordre canonique) compteurs, relève de naissance, frottement incrémental, assemblage —
le corps historique inchangé. AABB en parallèle ; génération des paires par fil (listes
concaténées puis triées : l'ordre canonique vient du tri, pas du découpage).

Contrôles : (1) même deck (`configs/yang2026_bench_s25_court.cfg`), même 14 fils, `history.csv`
de la version série (`rockim_g1y`) et de la parallèle (`rockim_g1y3`) comparés au caractère près
sur leurs **59 instants communs jusqu'à 74,6 µs**, fracture comprise (2 000+ joints rompus) :
**0 cellule différente** (`tools/ab_history.py`) ; (2) ancre de bit-identité 8 decks à 4 fils : **8/8 IDENTIQUE**
(`results/bitid_g1y3.log`). Vitesse : 75 µs du banc s = 2,5 en **438 s** (26 525 pas, 16,5 ms/pas
dont contact 12,3 ms, `ROCKIM_PROF=1`), contre ~28 ms/pas avant fracture et ~230 ms/pas après
pour la version série sur le même deck (mesures contendues par un autre run).

### Corrigé — insertion adaptative + `jointElastic = parabolic` (`rockim_g1y4.exe`, nuit)

À l'insertion adaptative le décalage de continuité `dn0` était calculé sur la branche linéaire
(`dn0 = s/pj`) quelle que soit la branche élastique : sous la parabole de Guo (éq. 2.31,
t = ft(2r − r²)) le joint naissait avec une traction 2s − s²/ft > s, soit +25 % de ft pour une
insertion pilotée par le cisaillement à s = ft/2 (exact seulement au pic, s = ft). Désormais
r = 1 − √(1 − s/ft), dn0 = r·dnE. Linéaire (défaut) : inchangé au caractère près ; la combinaison
n'existait dans aucun deck d'ancre. Trouvé en préparant la variante adaptative du deck Yang
(`docs/ADAPTATIF_impact_2026-09-11.md` §2.3).

### Ajouté — variantes de deck pour le point sur l'adaptatif (nuit)

`configs/yang2026_impact_adaptive.cfg`, `yang2026_bench_s25_adaptive.cfg` (quatre lignes changent :
`insertion = adaptive`, `insertionPenaltyFactor = 4`, `jointElastic = parabolic` et `gcBirth =
penalty` retirés — le second pour le triplet interdit du 02/09) ; `yang2026_bench_s25_sep.cfg`,
`_adaptive_sep.cfg` (`jointDeath = separation` + `jointResidualMu = 0.18` : un joint rompu comprimé
reste porté par la loi de joint au lieu de devenir une paire de contact). dt s = 1 en adaptatif :
**2,43 ns** contre 1,00 ns en intrinsèque. `tools/bench_compare.py` : tableau des bancs côte à côte.

### Corrigé — borne KE du garde-fou portée au 2D ; `jointShearUnload = origin` crée de l'énergie en 2D (`rockim_g1y10.exe`, 2026-09-12)

`FdemSolver::checkEnergyAbort` : même borne physique qu'en 3D (KE ≤ KE₀ + travail des sources,
hydro comprise), même clé `budgetAbortPct` ; ancre **8/8 IDENTIQUE** (`results/bitid_g1y10.log`).
Trouvé en lançant le 2D intrinsèque « comme Yang » : création d'énergie dès le contact (94 % des
joints rompus à 155 µs, résidu B4 muet). Vingt variantes de bissection
(`configs/impact2d_bis_V*.cfg`, ETAT §12) : joints incassables et continuum pur sains, une seule
clé de loi en cause — **`jointShearUnload = origin`** (éq. 18 de Yan, secante à l'origine dont
le seuil s_E suit |σ_n| courant : ressort paramétrique non conservatif) ; `plastic` est sain. Le
deck 2D du 03/09 relancé tel quel : Cundall + dashpot absorbent 1,6× le travail de l'outil ;
sans amortissement il abort à 58 µs — la campagne 2D du 03/09 est à relire. `jointShearRange =
coulomb` est attaché à `origin` par la garde 2D, donc indisponible en 2D tant qu'`origin` n'est
pas repris. Le 3D (`Fdem3dSolver`, huit bancs en `origin`, borne jamais déclenchée) est à
instruire séparément. Decks 2D propres : `configs/impact2d_kuru_{intrinseque,adaptatif}_plastic.cfg`
(600 µs : e = 0,760 et 0,784 pour 0,83 chez Yang, 19 et 5 joints rompus). `tools/impact2d_report.py`,
`meshes/impact2d_grad_clean.msh`.

### Ajouté — `facetAverage = nodal` (`rockim_g1y9.exe`, 2026-09-12) — le critère de Camacho-Ortiz, et il insère MOINS

Traction transmise par la facette liée, par partition des forces nodales (Camacho & Ortiz 1996,
Pandolfi & Ortiz 2002) : forces internes des copies du côté A du plan, **moins l'inertie**
(accélération du pas précédent, `accN_`, écrite dans `integrate()` sous la clé), forme symétrique
½[Σ_A − Σ_B] pour partager la charge extérieure d'un sommet de surface, attribuée à la facette au
prorata de son aire tributaire, t = −F/A_f, pré-filtre `max` ≥ 50 % du seuil, exige `insertion =
adaptive`. **Validé** sur un champ uniforme (`configs/tension3d_adaptive_*.cfg` : pic 9,72 MPa
pour ft = 10, même instant de rupture que `arith` et `max`). **Sur l'impact il insère le moins de
tous** : 33 facettes, 0 rompue (g1y8 sans inertie : 79 / 0 ; `arith` 728 / 27 ; `max` 2 781 /
412) — la force transmise est l'équilibre de tout le patch nodal, la mesure la plus lisse ; le
joint intrinsèque lit l'élément seul, dont l'analogue est `max`. Conservé comme option de
référence pour les champs réguliers ; `max` reste le critère des impacts
(`docs/ADAPTATIF_impact_2026-09-11.md` §8). Défaut inchangé ; ancre **8/8 IDENTIQUE** (`results/bitid_g1y9.log`).

### Ajouté — `facetAverage = max` (`rockim_g1y7.exe`, nuit) — LA correction de l'adaptatif en impact

Le critère d'insertion évaluait la traction sur la MOYENNE des contraintes des deux tétraèdres
(§2.1 éq. 10 de la note, `arith` ou `volume`). Sous l'insert, l'élément qui porte l'anneau de
traction hertzien a pour voisin un élément comprimé : la moyenne divise sa traction par deux et
rien ne s'insère (premier joint à 90 µs contre 64 en intrinsèque, 27 joints rompus contre 6 090
à 100 µs, banc s = 2,5). `max` retient le plus chargé des deux (rapport max de σ_n/ft_dyn et
|τ|/f_s, même DIF, même enveloppe) ; σ_n, τ et f_s passés à l'insertion sont ceux de cet élément.
Mesuré : premier joint à **79 µs**, **412** rompus à 100 µs (×15), 2 781 facettes insérées (×4),
cône à 8 mm de profondeur au lieu d'une peau de 2 mm, bilan fermé (KE 31,5 → 30,1 J). Ni la loi
de volume `dpr` ni le cap de compaction n'avaient bougé ces chiffres. Opt-in ; `facetStress()`
servie ailleurs reste la moyenne. Ancre : **8/8 IDENTIQUE** (`results/bitid_g1y7.log`). 2D : à porter.

### Ajouté — `lawPhase = <phase>` (`rockim_g1y6.exe`, nuit)

La loi de volume (`law = dpr | cdp | …`) ne porte que sur les éléments de la phase désignée, les
autres restent élastiques. Jusqu'ici `law` + plusieurs `phases` était refusé, ce qui interdisait
toute loi de volume sur le montage d'impact à six corps (roche, acier, carbure). Sans la clé :
garde historique inchangée. Motif : le banc adaptatif laisse les tétraèdres isolés sous l'insert
porter 1,4 GPa sans céder (`docs/ADAPTATIF_impact_2026-09-11.md` §5) — la loi de volume est la
correction, et elle ne concerne que la roche.

### Corrigé — garde-fou d'énergie aveugle à une création (`rockim_g1y5.exe`, nuit)

`budgetAbortPct` testait le résidu B4, qui inclut la correction leapfrog `biasW_` : sur le banc
s = 2,5 `jointDeath = separation` + `gcBirth = penalty`, l'énergie cinétique est passée de 31 à
**255 J** pour 42,8 J incidents (contact « −203 J », 7 000 joints rompus jusqu'à 67 mm de l'axe)
et le résidu affichait **[OK] à 8·10⁻¹¹ %**, le poste « intégration » (+68,9 J) absorbant la
création. Ajout, sous la même clé et la même tolérance, d'une **borne physique** : KE ≤ KE₀ +
travail des sources extérieures (outil, liaisons, pesanteur et tri, confinement) ; au-delà, arrêt
avec le hotspot. Sans `budgetAbortPct` : rien ne change (bit-identique, ancre 8/8
`results/bitid_g1y5.log`). La variante `separation + gcBirth = penalty` est disqualifiée en l'état
(cause non élucidée : la naissance calée d'un joint mort par ouverture ? à instruire).

### Mesuré — l'adaptatif en impact (nuit, `docs/ADAPTATIF_impact_2026-09-11.md`)

Banc s = 2,5 à 110 µs, physique identique : l'adaptatif tourne (dt ×2,2, calcul ×5-15 moins
cher, bilan fermé) mais **27 joints rompus contre 6 090** en intrinsèque à 100 µs, une peau de
3 mm au lieu d'un cône broyé de 24 mm ; les tétraèdres isolés sous l'insert portent 1,36 GPa,
facettes toutes rompues autour : le critère d'insertion fait son travail, c'est le continuum
élastique qui ne s'écrase pas. `jointFrictionMobilised = damage` ne change rien (607 facettes).
Correction = loi de volume sur la roche (`lawPhase`, ci-dessus) ; premier essai `law = dpr` au §7
du document.

### Ajouté — la loi de joint de Solidity mot à mot et `jointPenaltyLength` (`rockim_g1y16.exe`, 13/09, 13 h - 14 h)

Motif : le run s = 1 du 13/09 (deck v3-P, `plastic`) confronté à yang2026 relu (`docs/AUDIT_2026-09-13.md` §10) :
lit broyé **lâche** (2 892 facettes mortes à 173 µs, réaction 20-30 kN, bit 7,37 m/s, nPulv = 0) là où
Yang a un noyau **solide** à joints réversibles (~57 kN, 5,62 m/s, ~360 éléments pulvérisés). Demande de
Fernando (« Fais 1 à 3 ») : coder leur loi telle quelle, prendre leur pénalité, tester au banc s = 2,5.

- **`jointShearUnload = solidity`** (fdem, fdem3d) : transcription littérale de `Sigma_tau` (`Y3Dfd.c`
  l. 1078-1300, code public d'Imperial) — z **sans mémoire** (le joint guérit), tractions fonctions de
  l'ouverture et du glissement **courants**, op/ot/sp/st recalculés à chaque pas, f_s = c en traction et
  c + tanφ·|2 pj dn| en compression, z-curve 0,63/1,8/6 en dur (`kSolidityZ`), dpefm = 0 (aucun
  frottement de joint), point rompu = traction nulle, joint rompu à **deux points sur trois au même pas**
  puis mort immédiate. Gardes : parabolic + munjiza + majority + jointXi 0 + jointEtaN/S 0, pas de camacho.
  Rend inertes ratchet, normalProxy, shearRange, deltaC, residualMu, frictionMobilised, frictionScaled,
  jointDeath. Bloc placé en tête de la boucle des points, `continue` avant les branches historiques :
  le chemin par défaut est textuellement intact.
- **`jointPenaltyLength = min | local | edge`** (fdem3d) : la longueur h de pj = pf·E/h. `edge` =
  moyenne des trois arêtes de la facette = le `el` de Solidity (`S_N_direction`) et le h de Guo éq. 2.25,
  d'où pf = p0/(2E) sous parabolic (25 pour 3 000 GPa à 60 GPa). `min` = comportement historique.
- Bancs courts (scratchpad, deck `fdem3d_visc_yan_3d` + parabolic/munjiza/midedge/majority, 4 fils,
  155 s) : `solidity` tourne, résidu B4 −3e-15 J, pic macro 9,95 MPa = ft, 52 joints rompus à 140 µs
  (plastic : 60 ; poste joints 0,137 J contre 0,158 J) ; la variante `jointFailRule = any` est
  **refusée** par la garde.
- Decks du banc s = 2,5 à 300 µs : `configs/yang2026_bench_s25_v3P_300.cfg` (A, témoin),
  `yang2026_bench_s25_solidity.cfg` (B : solidity + edge + pf 25),
  `yang2026_bench_s25_solidity_visc.cfg` (C : B + `bulkViscosity = 4000`). File d'attente
  `tools/queue_benches_s25.sh` derrière le run s = 1 (un gros job à la fois) et derrière l'ancre.
- Registre régénéré (429 clés). Ancre 8/8 + 9ᵉ deck **1/1 IDENTIQUE** : `results/bitid_g1y16.log`
  (13 h 15).

### Ajouté — contact de Solidity dans les bancs, banc D, liste des écarts Guo/Solidity/rockim, chasse aux slivers (13/09, 13 h - 14 h 30)

Demande de Fernando à 13 h : « Kill s1 » (fait, 182,9 µs, 19 trames, 2 éléments pulvérisés apparus à
183 µs), « le même contact que Solidity pour comparer », « la liste des écarts avec la thèse de Guo, juge
si c'est handicapant, puis s'il faut copier Solidity », « optimise le mailleur pour éviter les slivers ».

- Bancs B et C réécrits avec **le contact de Solidity** : `potPenaltyFactor = 0.25` (= p0/200 = 15 GPa,
  `Y3Did.c` `penalty = MINIM(d1pepe)/200`) et `gcBirth = penalty` (recalage de naissance [0,01 ; 3]) ;
  banc **D** = témoin A + ce seul contact (`configs/yang2026_bench_s25_v3P_contact.cfg`,
  `tools/queue_bench_D.sh` derrière la file A/B/C). Banc C corrigé : η = 4 000 de Guo éq. 2.6 (T = … + η D)
  vaut `bulkViscosity = 2000` dans la convention 2 μ D de rockim.
- **`docs/ECARTS_guo2014_rockim_2026-09-13.md`** : 26 lignes Guo (chapitre 2 lu dans le texte) /
  Solidity (code) / rockim, avec le jugement « handicapant ? » et la décision. Six écarts comptent : la
  mémoire du joint (banc B), la longueur de pénalité et le pas de temps (`edge`, mailleur), l'amortissement
  (banc C), la **résolution** (notre s = 1 rend une arête médiane de 1,37 mm dans la boule, 14 722 tétras
  au lieu de ~35 000 : un s ≈ 1,37 par rapport à Yang), la pulvérisation (δ_m et Cd inconnus). Faut-il
  copier Solidity : la physique qui se copie l'est (opt-in) ; le code public n'a rien de plus.
- **`tools/mesh_quality.py`** : diamètre inscrit h = 6 V/Σ aires par tétra et par corps, pires éléments.
  Mesure : les slivers de la roche (h 0,226-0,30 mm) ont deux nœuds sur z = 0 et deux à 0,7-1,1 mm sous
  la surface ; Netgen/Relocate3D ne bougent pas les nœuds de surface. `tools/make_impact_mesh.py` :
  options `algo2d=`, `optthr=`, `smooth=`, `optgmsh=`, `opt2d=`, `lap2d=` et le **preset `quality=hxt`**
  (HXT + seuil 0,5) : roche 0,2265 → 0,270 mm (+19 %), slivers < 0,3 mm 25 → 5, insert 0,205 → 0,250 ;
  les autres combinaisons font moins bien (tableau dans ECARTS §4). Maillage `meshes/impact_yang_s1_pose_hxt05.msh`
  généré pour le run suivant. MMG3D absent de gmsh 4.15.2.

### Corrigé — filtre des facettes rompues, bancs d'attribution B1/B2 (13/09, 15 h, relecture des diagnostics indépendants du 12/09)

- **`bench_impact/tools/imp_lib.py` `broken()`** : lisait `damage ≥ 0,999`, qui vaut 1 dès qu'UN point
  d'intégration cède alors que la facette ne rompt qu'à deux points sur trois (`jointFailRule = majority`) ;
  lit désormais `tBreak ≥ 0` (repli sur l'ancien filtre si le champ manque). Trame 18 du s = 1 : 3 357
  rompues et non 3 737 (380 faux positifs, 10 %). Toutes les figures « joints » et les comptes antérieurs
  du 13/09 sont ~10 % trop hauts ; conclusions inchangées.
- Bancs **B1** `configs/yang2026_bench_s25_law.cfg` (témoin + loi solidity seule) et **B2**
  `yang2026_bench_s25_pen.cfg` (témoin + `edge` + facteur 25 seuls), file `tools/queue_bench_E.sh` derrière D.
- `docs/ECARTS_guo2014_rockim_2026-09-13.md` §5 : relecture point par point des deux diagnostics (filtre,
  étiquette de mode, cisaillement pré-pic linéaire sous `plastic`, h_e = diamètre inscrit dans δ_m, 3 000 GPa
  = St Anne, cap caché `meanTensionCapFactor = 3`, estimateur de vitesse d'indentation 6,86 m/s, St Anne
  d'abord).

### Ajouté — le lot de l'audit (`rockim_g1y14.exe`, 13/09, 2 h 30 - 3 h) : proxy σ_n, ratchet armé au cap, garde de phase, raideur de contact par phase, joints ×6

Fernando : « refais un audit pour vérifier s'il n'y a pas d'autres erreurs ou optimisations ». Quatre
agents en parallèle (loi de joint, contact/dt/bilan, performance, deck contre Yang 2026), rapports
et décisions dans `docs/AUDIT_2026-09-13.md`. Codé, tout opt-in ou bit-neutre par construction :
- **`jointNormalProxy = penalty | law`** (2D + 3D ; audit A #1, BLOQUANT) : sous `jointElastic =
  parabolic` la loi transmet 2·pj·dn en compression mais s_E, la plage coulomb et l'amorçage du DIF
  lisaient pj·dn — σ_n divisé par deux dans tout ce qui fixe une résistance de mode II. `law` pose
  le facteur 2 (Solidity lit σ_tmp = pe·o/el, la vraie contrainte). Le run s = 1 de 2 h 20 (g1y13)
  portait ce défaut : relancé.
- **ratchet de mode II armé au cap** (audit A #4) : `Joint::ksr` n'est mémorisé qu'une fois
  τ_env < pj·s_max — un joint inséré sous traction ne peut plus naître verrouillé sur un τ_lim
  transitoire. Change le comportement de `jointSecantRatchet = on` seulement (clé du 12/09).
- **`bulkDamagePhase = <phase>`** (audit D) : les éq. 3-4 de Yang n'avaient aucune garde de phase —
  acier et carbure endommageables (δ0 = 14 µm contre h·ε ≈ 9-13 µm à 175 MPa). Avertissement si
  absent avec plusieurs phases.
- **`potStiffnessByPhase = max | min`** (audit B #10) : potP_ et potKt_ étaient bâtis sur
  phases_.maxE() = 600 GPa pour toutes les paires — roche/roche 10× trop raide. `min` : pénalité
  de la paire = potPenaltyFactor × min(E_A, E_B), k_t au prorata ; le budget de dt garde potKt_.
- **Fusion parallèle des forces de joint** (audit C #1) : l'ancienne fusion fTL_/touchedTL_ était
  série (30,7 ms sur 32,4 ms de joints à s = 1, 95 %). Fusion par nœud, même ordre t = 0..nT−1 :
  bit-identique par construction. Remise à zéro de f_ en parallèle (C #7). Ping-pong cur/nxt des
  fragments dans `pot3::pairForce` (C #3, même geste que sur Poly3).
- Registre régénéré (428 clés). Ancre 8/8 + 9ᵉ deck : voir results/bitid_g1y14*.log.
Non codé, consigné dans l'audit : dtContactAudit, gcBirthWork, difLengths, relay figé/mesuré, CSR de
grpsOfVert_, AVX2, remaillage 0,7/1,0 mm, masses du train de frappe (−6 à −9 %), forme de D (éq. 4
vs Fig. 4a), T et fragBrush pour les critères 3 et 5.

### Ajouté — `facetAverage = max` en 2D (`rockim_g1y13.exe`, 12/09, 2 h) — le port du critère des impacts

Item en attente depuis le 11/09 (ADAPTATIF §8). `FdemSolver::insertionSweep()`, branches OpenMP et
série : sous `max`, le critère lit le plus chargé des deux triangles (rapport max(σ/ft_dyn, |τ|/f_s),
même DIF, même enveloppe) et REMPLACE (σ, τ, f_s) de la moyenne ; l'offset de continuité dn0 porte
alors la traction du triangle qui déclenche, comme en 3D (g1y7). Chemin mort sous `arith`/`volume`.
`rockim_g1y12.exe` (relay seul) a été compilé à cheval sur l'édition de `FdemSolver.hpp` — jeté,
jamais utilisé (piège COMDAT) ; `g1y13` = rebuild complet propre avec relay + max.

### Ajouté — `gcBirth = relay` (`rockim_g1y12.exe`, 12/09, 2 h) — le relais joint mort → contact à force continue, sans l'injection du premier contact

D4/N7 du conseil. `penalty` (Solidity, Y3Did.c l. 915-964) cale la pénalité de la paire naissante
sur fn_joint/fn_contact pour TOUTES les paires — y compris celles qui n'ont pas de joint mort
(premier contact piston/bit), où fn_joint = 0 donne un facteur 1 sur un recouvrement déjà formé :
½kδ₀² injectés, 1 J, abort à 2,9 µs (`_certif`). `ramp` (défaut) part de zéro sur gcBirthTau pour
toutes les paires — sous l'insert, où 90 % des joints meurent en compression (2 à 3,7 MN lâchés
cumulés sur le banc s = 2,5), c'est une coupure du chemin d'effort qui nourrit le terme leapfrog.
`relay` = `penalty` pour les paires nées d'un joint mort (`c.jI >= 0`, continuité de force, raideur
tangentielle re-échelonnée sur les seules paires calées `H.penScale >= 0`), `ramp` pour les autres.
`gcBirthTau` reste actif (rampe des paires sans joint). Défauts bit-identiques (`ramp` et `penalty`
inchangés au test près). Corrigé au passage : `gcBirthPenMin/Max` étaient imprimés avant d'être lus
(la bannière disait toujours 0,01 / 3,0). Decks : `configs/yang2026_bench_s25_{ratchet,plastic_coulomb}_relay.cfg`.

### Ajouté — la loi de joint passée au conseil : `jointSecantRatchet`, coulomb sur `plastic`, avertissement `origin` (`rockim_g1y11.exe`, 12/09, 1 h 30)

Fernando : « tout doit être mathématique et physique ; corrige les problèmes mathématiques, écris les
équations correctement, regarde les vraies équations, puis l'implémentation numérique ». Plan écrit
(`docs/PLAN_loi_joint_2026-09-12.md` : rockim telle que codée, Solidity `Sigma_tau` transcrite
verbatim depuis leur source GitHub, Yang 2026 §2, six constats M1-M6), soumis à un conseil de quatre
critiques + synthèse séparée (verdict REVISE, liste N0-N7 ; sorties dans le plan et ETAT §13).

Ce que le conseil a établi et que le code porte maintenant :
- **M1 — `jointShearUnload = origin` est non conservatif** : τ = τ_lim(σ_n)/s_max · s dès que le cap
  est actif, ∂τ/∂δ_n ≠ ∂σ/∂δ_s, un cycle à glissement fixé crée ½(k₂−k₁)s². Mesuré : 16 J en 81 µs
  (référence), 0,1 J sans les conventions Solidity (B4), sain sous `plastic` (B1). **AVERTISSEMENT**
  imprimé au démarrage (2D et 3D), formule et chiffre ; rien n'est refusé (principe VIII, 40 decks
  archivés restent rejouables ; `budgetAbortPct` coupe par la mesure).
- **M7 (sceptique) — le DIF continu est une raideur, pas seulement un critère** : sous
  `strainRateDIFArm = continuous`, `refreshDif` réécrit ft et dnE à chaque pas, donc l'enveloppe
  parabolique de mode I et la sécante de l'éq. 17 suivent ft(t) non monotone — une seconde pompe, en
  mode I, quelle que soit la branche de cisaillement. D5 du plan (« rien à corriger ») était faux.
- **`jointSecantRatchet = on | off`** (nouvelle clé commune, défaut off bit-identique) : les sécantes
  de décharge des éq. 17 (mode I, toutes branches) et 18 (mode II, branche `origin`) deviennent
  NON CROISSANTES dans le temps (`Joint::knr`, `Joint::ksr` par point d'intégration). C'est la
  condition Φ ≥ 0 d'une loi d'endommagement : ni la remontée du DIF ni celle de τ_lim avec la
  compression ne peuvent restituer plus que la charge n'a stocké ; le DIF et la pression continuent
  d'élever l'ENVELOPPE atteignable en charge croissante. `origin` + ratchet = la forme dissipative de
  l'éq. 18 (D3 du plan), `plastic` + ratchet = plastic avec le mode I protégé du DIF.
- **M4 — garde `jointShearRange = coulomb` ⇒ `origin` levée** (2D et 3D) : la normalisation
  3 G_II/f_s(σ_n) de Solidity est une longueur de référence du moteur, pas une raideur ; sur `plastic`
  le moteur devient |s_p| / max(2 s_E, s_F·c/f_s) (dissipatif, D ratchet, pj constant). `slipRef`
  vaut `J.slipF` au bit près quand la clé est absente ou sous `origin`.
- Registre `tools/keys_by_mode.json` + `KeysByMode.hpp` régénérés (425 clés, `jointSecantRatchet`
  commune). Build complet dans `obj_g1y11` (jamais de lien partiel après un .hpp).
- **Ancre 8/8 IDENTIQUE** (`results/bitid_g1y11.log`, `bitid_g1y11.json`). **Ancre latérale créée** :
  `tools/bitid_refs_jointlaw.json` (deck `fdem3d_yang_v2_court`, 339,6 s à 4 fils, prise sur
  `rockim_g1y10.exe`, sha 08e46e22…) — raison : aucun des 8 decks principaux n'exerce `parabolic` +
  `yang` + `coulomb` + `midedge` + `majority` + DIF continu, le chemin même que ce lot touche
  (remarque des critiques pré-mortem et ingénieur). `bitid.py` : champ `side=True`, joué seulement
  sur `--only` ou sur un `--refs` autre que l'ancre principale ; la passe 8/8 par défaut est inchangée.
- **N0 du conseil** (DIF gelé / désarmé sur B4, 90 µs) : joints −0,10 / −0,15 / −0,18 J, leapfrog
  +1,73 / +1,66 / +1,23 J, 1 375 / 1 396 / 1 550 rompus — le DIF n'est pas une pompe mesurable sur une
  fracture à 95 % en cisaillement comprimé (M7 reste vrai en traction ouverte, le ratchet le couvre).

Décisions du conseil non codées (volontairement) : refus de `origin` sous `budgetAbortPct`
(cassait 40 decks et l'ancre `cut3d_heilman_court`), `JointLaw.hpp` (le « miroir 2D/3D exact » est
faux : etaN/etaS, jsOn_, midedge 0,5/0,5 vs (1−tq, tq), majority 3 vs 2 points — à consolider à froid),
`gcBirth = relay` (seulement si un candidat retenu tue > 500 joints en compression), le pilote
`joint_cycle` (point-matériel, aveugle à `midedge` — à porter au niveau facette et fondre dans une
instrumentation ∮τ·dδ_s par facette).

Decks du conseil : `configs/yang2026_bis3d_B4_difins.cfg` / `_difoff.cfg` (N0 : DIF gelé à
l'enveloppe / désarmé sur B4), `_B4bd.cfg` (chaînon manquant : sans midedge ET majority),
`configs/yang2026_bench_s25_plastic_coulomb.cfg` (N3, la loi D1 enfin mesurée, 200 µs),
`configs/yang2026_bench_s25_ratchet.cfg` (N4 : origin + coulomb + ratchet, 200 µs),
`tests_f2/bitid/fdem3d_yang_v2_court.cfg` (N1 : 9ᵉ deck d'ancre portant les cinq clés du deck v2,
bit lancé à 9 m/s, 20 µs ; référence latérale `tools/bitid_refs_jointlaw.json`).

### Ajouté — figures « rien que les joints » (12/09, 0 h - 0 h 30)

Demande de Fernando : valider les fissures à l'œil, joints seuls, arêtes seules, puis des coupes.
- `tools/fig_joints_only.py out_dir [--frame k] [--stem] [--title] [--lim 30] [--depth 25] [--faces]` :
  six vues des facettes rompues — rangée 1 par MODE (rouge traction, jaune cisaillement : 3D, dessus,
  coupe |y| < 5 mm), rangée 2 par ORIENTATION (sub-verticales |n_z| ≤ 0,6 = fissures, sub-horizontales
  = broyé). Rendu par arêtes par défaut (faces transparentes), `--faces` pour le plein. Cercles de
  référence : insert 8,51 mm, Yang 9 m/s cratère 7 mm et radiales 10 mm. Imprime rompus, % traction,
  % sub-verticales, extension radiale, rayon de cratère, profondeur.
- `tools/fig_joints_cuts.py out_dir [--ys 0,5,10] [--zs -1,-3,-6,-10] [--lw 2]` : COUPES EXACTES —
  chaque facette rompue est intersectée avec le plan, la trace est un segment (pas un triangle
  projeté) ; verticales x-z à y = 0, 5, 10 mm et y-z à x = 0, horizontales à z = −1, −3, −6, −10 mm ;
  nombre de traces et longueur cumulée par coupe.
- Les deux réutilisent `bench_impact/tools/imp_lib.py` (read_vtu, broken, frames). Figures dans
  `results/fig/joints_*`, `joints_aretes_*`, `coupes_*` pour les trois bancs s = 2,5 (intrinsèque
  80 µs, adaptatif + max 110 µs, adaptatif moyenne 110 µs). Lecture (ETAT §13) : l'intrinsèque est
  une cuvette de cisaillement 3 à 5× trop large (l'image de la pompe d'énergie), l'adaptatif + max
  un cône broyé au rayon de Yang (7 mm) avec anneau de traction à 12-15 mm ; aucune radiale nulle
  part — à 2,5 mm de maille une fissure de 10 mm fait quatre arêtes.

### Mesuré (ETAT §5-7)

- Coût par pas v1 (14 fils, 120 k tets) : éléments 11 ms, joints 32 ms, contact 50 ms — joints et
  contact **ne scalent pas** avec les fils (×1,05 et ×1,04 de 7 à 14).
- dt v1 0,942 ns (joints du carbure) → v2 1,00 ns (slivers de la roche) ; CFL seule 4,67 ns ;
  pf 10 : 1,37 ns.
- Banc court s = 2,5 : voir ETAT §7.

## [Non publié] — arbre g1, portage de la loi de la note 2026

Ouverture de l'arbre `g1` (copie conforme de `g0` le 2026-09-11) pour porter la loi de la note de
travail *« Lois constitutives proposées pour un FDEM hybride à insertion adaptative »* (septembre
2026) : matrice viscoplastique de Mohr-Hoek parfaite + endommagement de **compression seul**, joints
cohésifs **extrinsèques** initialement rigides, mode mixte de Benzeggagh-Kenane, dépendants de la
vitesse. Contrat d'implémentation : `docs/SPEC_loi_note_2026.md`. Audit de départ :
`../rockim_g0/AUDIT_loi_adaptative_2026-09-11.md`.

Règle du lot : **croissance par addition** (principe VIII). Toute capacité nouvelle arrive par une clé
opt-in dont le défaut reproduit `g0` bit pour bit ; aucun défaut existant n'est modifié — `capP0`,
`jointXi`, `dampingLocal` et `erodeD` s'éteignent **dans le deck**, pas dans le code.

### Ajouté

- `include/rockim/JointTsl.hpp` — noyau **partagé** de la loi cohésive extrinsèque, fonctions pures,
  appelées à l'identique par le solveur 2D et le solveur 3D de sorte que la loi ne puisse pas diverger
  entre les deux par recopie : critère elliptique Φ_F (eq. 12), DIF (eq. 13), séparation et traction
  effectives de Camacho-Ortiz (eq. 16), tampon d'insertion avec G_C de Benzeggagh-Kenane **gelé au
  ratio de mode à l'insertion** (eq. 17-18), partition (eq. 19), cap de cisaillement
  `(1−D)·cohésion + D·μ⟨−t_n⟩`, `h_e = (12 V_e/√2)^(1/3)` et borne de snap-back (eq. 6-7).
- `MatState::lcComp` — longueur de bande de la **compression**, seul point de contact inter-fichiers
  du lot. `≤ 0` = comportement de `g0` (β_c prend `lc`), bit-identique.

### Ajouté — les 21 clés de la loi (lots A/B/C du 2026-09-11)

Matrice (`MatLaw`) : `bulkTensionDamage`, `mhForm`, `dpFlowForm`. Solveurs FDEM 2D et 3D, mêmes
noms : `facetAverage`, `facetRate`, `insertionCriterion`, `insertionHoldSteps`, `difExpT`, `difExpS`,
`jointTSL`, `jointMixLaw`, `jointBKEta`, `jointFrictionMobilised`, `jointEtaN`, `jointEtaS`,
`jointViscousInCriterion`, `gbCombine`, `jointWeibullXu`, `jointWeibullScale`, `compBandLength`,
`energyBreakdown`. Toutes opt-in, défaut = `g0`. Registre régénéré (`tools/gen_keys_by_mode.py`,
422 clés). **Bit-identité : 8/8 IDENTIQUE** (`tools/bitid.py --exe rockim_g1.exe`, ancre `g1ref`),
plus un neuvième contrôle hors ancre : `configs/fdem3d_percussion.cfg` (outil libre) rend
24 754,2 N / 1 259 joints / 1,90983e6 J/m³ avec les deux binaires. Trois refus nouveaux, tous sur
des clés jusqu'ici inertes : `erodeD/erodeEpv/erodeDc` en FDEM ; `jointTSL = camacho` sans
`insertion = adaptive` ; `strainRateDIF` armé avec `jointEtaN/S > 0`.
Deck de la loi complète : `configs/loi_note_2026.cfg` (0 clé non lue, 21 482 tets, 408 grains).
Trois points où les lots ont eu raison contre le contrat : `difExpS` (voir SPEC §2) ; le Newton du
méridien puissance en σ₃ = 0 (tangente verticale, dépassement de 2,4 % corrigé par un Newton borné) ;
l'affichage de s_MH (facteur 5,6 d'unités attrapé au banc avant livraison).

### Corrigé — la loi initialement rigide était instable (seconde passe, 2026-09-11)

**Ce qui a été observé.** `configs/loi_note_2026.cfg` explose à ~91 µs : bloc de 40 mm, nœuds à 18 m
(`out_note2026`). Le résidu d'énergie affichait `[OK] 9,2e-6 %` — c'est une identité comptable, le
poste « éléments » absorbait l'énergie créée. Deux contre-essais : `jointFrictionMobilised = off`
explose plus tôt (`_fricoff`, trame 8) ; `dtFactor = 0,043` (le pas historique, ÷3,49) explose à
la **même trame**, moins violemment (1 m, `_dtsafe`). Une instabilité CFL ordinaire est binaire :
ce n'en était pas une. Mesure à la trame 9 de `_dtsafe`, joints insérés : raideur sécante
t_ins/(D·δ_f) p99,9 = **210 × pj**, max **549 × pj** ; 651 joints à 0 < D < 1e-3 ; `eGc` de −45 J
à −1 707 J entre 84 et 100 µs pour 3,2 J injectés par l'outil.

**Diagnostic.** Une loi initialement rigide a une raideur de charge/décharge **non bornée**
(t_ins/δ_max → ∞ pour un joint qui s'est ouvert d'un rien puis recharge), hors de tout budget CFL
quel que soit dt ; et à l'insertion `jtsl::split()` rendait 0 à δ_m = 0 là où la facette liée
transmettait t_ins — un Dirac de −t_ins par insertion, ×30 000. Pathologie classique des lois
extrinsèques en explicite (Papoulia, Sam & Vavasis 2003), anticipée par la note au §2.4.

**Correctif** (`include/rockim/JointTsl.hpp`, `YangDif.hpp`, puis câblage 2D/3D) : branche
ascendante courte `jointTSLRise` (défaut 1e-3 sous `camacho`, la valeur de la note) — décalage
δ_0 = rise·δ_f posé dans la direction de la traction tamponnée, raideur bornée k₀ = t_ins²/(2·rise·G_C),
traction continue à l'insertion, ∫t dδ inchangé = G_C (banc `tests_f2/check_jointtsl_header.cpp`,
14/15 puis 15/15 après correction d'un seuil du banc). Budget CFL : **toutes** les facettes, liées
comprises, pèsent max(k₀, kPara·pj)·A0/3.

**Retrait.** Le « gain ×3,49 sur dt, saturé » annoncé par le lot B est retiré : il venait d'exclure
les facettes liées du budget, ce que le lot C avait refusé en 2D avec la bonne raison. Le ×1,96 de
l'audit comparait à `insertion = none` (aucun joint) et n'était pas une comparaison valide.

**Recette de la seconde passe** : build 2 352 128 o, registre 423 clés (`jointTSLRise` → fdem, fdem3d),
`tools/bitid.py --exe rockim_g1.exe` → **8/8 IDENTIQUE** (`results/bitid_g1_rise.json`).

**Troisième passe (même jour) — le vrai mécanisme de l'explosion, et son correctif.** La branche
ascendante était nécessaire (raideur bornée, plus de Dirac) mais **pas suffisante** : avec elle, le
deck explose encore (0,66 m). Cinq isolations sur `configs/loi_note_2026.cfg` ont désigné la source
(détail : `docs/RETOUR_v3_2026-09-11.md` §1.4 à 1.4 quater) : (b) matrice de la note + joints
**pénalisés** de g0 → 100 µs sains, 3,3 % d'insertions ; (a) volume **purement élastique** + joints
Camacho → **explose**, 67 GPa, 97 % d'insertions ; (c) cône DP de g0 → même trajectoire ; (d) η ÷ 1000
et (e) φ_j = 0 en cours. La matrice — viscosité et retour principal compris — est innocentée ; la
source est la branche Camacho : **le frottement de Coulomb y était appliqué comme une traction de
magnitude constante [D]·μ⟨−t_n⟩ dirigée par le déplacement de glissement** (3D : `tau = (lim/dsEffN)·
dsEffV` avec `lim = tauCoh + D·μ⟨−t_n⟩` ; 2D : `tau ± fr` selon le signe de `dtg`) — un ressort sec à
force constante, discontinu à l'origine, sur chaque joint inséré ou rompu comprimé. La v2 le
suggérait (*« s'ajoute »*), les deux lots l'ont implémenté littéralement ; la v3 §3.4 écrit la bonne
forme (opposé à la **vitesse**, régularisé). **Correctif** (opt-in, branche pénalisée intacte au
caractère près) : la part frottante devient un **cap sur une traction d'essai de collage avec retour
de glissement**, le mécanisme de g0 — `rockim::camachoFrictionSlider` (3D, `Fdem3dSolver.hpp:72-136`,
banc `tests_f2/check_camacho_friction3d.cpp` 13/13 : dissipation par cycle = 2·f_cap·(P − 2f_cap/pj),
jamais négative, continue en δ_s = 0) et `rockim::jfric::capReturn` (2D, `FdemSolver.hpp:93-155`, banc
`tests_f2/check_camacho_friction.cpp`). Build `rockim_g1fix.exe` 2 353 664 o ; **bit-identité 8/8 IDENTIQUE**
(`results/bitid_g1fix.json`). Validation du correctif sur les mini-bancs 4×4×4 (384 tétraèdres,
`configs/mini/`) : `elas_cam_T80` passe de 8,15 mm / 6 812 MPa / 663 rompues / **+402 J créés** à
**6,008 mm / 79,7 MPa / 50 rompues / 0,21 J** (bloc initial 6,000 mm) ; `note_camnobk`, qui explosait
à ~9 µs, termine. Contrôle 2D (`configs/_ucs2d_camacho_check.cfg`, UCS adaptatif + camacho corrigé) :
pic 41,7 MPa — la valeur de l'ancre — résidu d'énergie **2,5·10⁻¹³ %**.

**Lecture finale du mécanisme (six isolations).** La divergence demande DEUX ingrédients :
(1) la **source**, le terme frottant non conservatif — le retirer (`frictionDeg = 0`, isolation e)
suffit à stabiliser ; (2) le **gain**, l'incapacité de la matrice à céder — l'amarre vaut
`D·μ ×` pression de contact, donc sa puissance créée croît comme le carré de cette pression : matrice
qui plafonne ⇒ création sous la dissipation (isolation d, `saksalaEta` ÷ 1000 : **saine**, joints qui
dissipent −13,3 J, 15,9 % d'insertions) ; matrice élastique (a), cône inatteignable (c) ou
sur-contrainte `s_MH·κ̇` non saturée au-delà de 10³ /s ⇒ emballement exponentiel. Les isolations (a)
et (c) ne dissociaient donc pas « viscosité » de « matrice qui ne cède pas » : elles les confondaient
(un solide élastique est le cas limite d'une viscosité infinie).

**Résorbé au passage** : `difCompressionYang(edot, n)` dans `YangDif.hpp` (les deux transcriptions
locales des lots sont supprimées — une seule dans le dépôt, comme l'en-tête l'exige) ;
`jtsl::shearCap` porte l'enveloppe de Yang (les deux duplications manuelles du facteur D sont
supprimées).

### Modifié — ancre de bit-identité

**Ancre reprise sur `rockim_g1ref.exe`** (sha256 `76a3dbf1a283729d…`, 2 263 040 octets, build de
référence de l'arbre `g1` à sources inchangées), 8/8 decks à OMP 4. Raison : l'ancre héritée datait de
`rockim_f2w18.exe` (2026-09-05) et ne valait déjà plus pour `g0` — ce CHANGELOG le consignait lui-même
(« rend 5/8, les trois cas `fem3d` ayant changé par les seules colonnes ajoutées, forces de pic
inchangées : 7916,75 / 6775,14 / 30390,2 N ; les cinq cas `fdem` restés identiques »). Mesure refaite
le 2026-09-11 : **identique au mot près**, ce qui confirme que la divergence est bien de la métrologie
et non de la physique, et que le **chemin FDEM — celui que ce lot modifie — n'avait pas bougé**.
L'ancre héritée est conservée sous `tools/bitid_refs_w18_herite.json` pour la traçabilité.
Le maillage `etude_lois_fem/meshes/T1_c05_clean.msh`, absent de la copie initiale, a été rapatrié
(sans lui le deck `fem3d_dpr_T1_court` échouait au lieu de se comparer).

## [Non publié]

### Corrigé — lot A/B du guide de correction du 2026-09-06

Trois verrous étaient annoncés ; **deux sont levés**, le troisième (`CDP-06`,
tessellation portable) ne l'est pas — son test d'acceptation demande une seconde
plateforme. Le **ré-ancrage unique** des empreintes attend donc `CDP-06` : en
l'état, `python tools/bitid.py --exe rockim_fix.exe` rend **5/8**, les trois cas
`fem3d` ayant changé par les seules colonnes ajoutées (forces de pic inchangées
à la précision affichée : 7916,75 / 6775,14 / 30390,2 N), les cinq cas `fdem`
étant restés **identiques**.

#### `HET-03` — le ressort absorbant avalait le confinement (A1)

`Fem3dSolver` mémorise l'état de référence `uConf_` et le ressort de
Deeks–Randolph travaille sur `u - uConf` au lieu de `u`.

**Le correctif du guide, pris à la lettre, ne passe pas son propre test.** Figer
`uConf_` à `confineGaugeTime` en laissant le ressort actif *avant* ne change
rien à la jauge (mesuré : toujours −75,31 MPa) : le ressort a déjà mangé sa part
pendant la rampe, et l'état figé est l'état dégradé. Il faut que le ressort soit
**inerte pendant la mise en pression** — `uConf_` suit `u_` tant que
`confineGaugeTime` n'est pas franchi, puis se fige. C'est aussi ce que dit la
physique : le massif environnant subit la même compression statique que le bloc.

- **Test d'acceptation : PASSE.** Banc T3 (`B1_T3_springs_DOIT_ECHOUER.cfg`,
  bloc 48 × 48 × 32, P = 100 MPa, sans outil) avec `absorbSpringFactor = 1` :
  jauge **−99,983 MPa** contre une consigne de −100, identique au témoin sans
  ressorts (−99,983). Avant : −75,31 MPa.
- **Non-régression à P = 0 : PASSE.** Percussion 30 µs, `absorbing = all`,
  `absorbSpringFactor = 1`, `confiningPressure = 0` : 625 lignes d'historique,
  **zéro écart** sur les 27 colonnes préexistantes, trois VTU **bit-identiques**
  au binaire de référence. La garde est explicite : sans confinement,
  `uConfSet_` passe à vrai au premier pas avec `uConf_ = u_ = 0`.
- L'avertissement « confinement + ressorts : poser absorbSpringFactor = 0 » est
  levé dans le solveur, et le refus correspondant dans
  `etude_lois_fem/make_cfgs.py` aussi ; les deux `absorbSpringFactor = 0` codés
  en dur y sont retirés (le défaut du code, 1, redevient le défaut des decks).

#### `HET-15` — le bilan d'énergie du bloc (A2)

Les deux compteurs de `Fdem3dSolver` sont portés en `fem3d` — énergie élastique
**stockée** dans les ressorts (sur `u - uConf`, cohérent avec A1) et travail
**cumulé** des amortisseurs de Lysmer — **plus un troisième que le guide n'avait
pas identifié** : le travail de l'amortissement local `dampingLocal`.

Sur le run de référence (`matrice_v2/C_T1_R_P000`, 300 µs, 98 342 tétraèdres),
rejoué ici et reproduisant le symptôme au millième près :

| postes comptés | résidu | % du travail de l'outil |
|---|---|---|
| dissipations de loi seules (état antérieur) | 3,069 J | **21,29 %** |
| + ressorts (0 J) et Lysmer (2,271 J) | 0,798 J | 5,54 % |
| + amortissement local (0,168 J) | 0,630 J | **4,37 %** |

**Le critère d'acceptation du guide (résidu < 1 %) n'est pas atteint.** Il reste
0,63 J sans compteur, et l'explication la plus plausible est l'énergie élastique
encore **stockée** dans le bloc à 300 µs — un quatrième poste, absent des deux
solveurs. À traiter séparément.

Trois colonnes sont ajoutées à `history.csv`, **en toute fin de ligne et sous
condition** : `uSpring`, `wLysmer` si `absorbing` est posé, `wDampLocal` si
`dampingLocal > 0`. Les decks de `bench_phases` n'emploient ni l'un ni l'autre :
leur `history.csv` reste **bit-identique**, contrairement à ce que le guide
annonçait.

#### Lot B — les six items triviaux

- **B1 (`CDP-04`)** — `Material.hpp` : les contrôles de signe passent **avant**
  le contrôle de finitude de `sqrt(E/rho)`. *Réserve mesurée* : le symptôme
  décrit (message accusant `rho` alors que `E` est négatif) **n'est pas
  reproductible** au HEAD — les deux binaires rendent
  `Material (global): E must be > 0`, et `bench_phases/erreurs.py` donne
  **16/16 avant comme après**. La correction reste juste sur le fond (les signes
  avant les grandeurs dérivées) ; elle ne corrige pas un symptôme observable.
- **B2 (`CDP-03`)** — `ThermoBench.cpp` : la déformation du pas est conservée
  (`epsSeq`) et versée dans les colonnes `e11..e13` du test 4, à la place de
  `sgA[k]` qui est une **contrainte**. *Réserve* : aucune loi ne produit
  aujourd'hui de violation du test 4, donc aucune ligne n'exerce ce chemin.
- **B3 (`ORCH-02`, `HET-16`)** — `--help` et `-h` sont reconnus (ils rendaient
  `Config: cannot open '--help'`), et l'usage liste les **onze** sous-commandes
  `selftest-*` au lieu de quatre.
- **B4 (`CDP-15`)** — sans objet : `LISEZ_MOI.md`, `CHANGELOG.md`,
  `DOCUMENTATION_rockim.md` et `bench_phases/LISEZ_MOI.md` disent déjà **seize**.
  Seul le message du commit `10344ec` dit quatorze, et il est immuable.
- **B5 (`CDP-16`)** — les cinq `*.w0` sont retirés du dépôt (`git rm --cached`,
  fichiers conservés sur disque) et `*.w0` est ajouté au `.gitignore`.
- **B6** — `FemSolver::Elem::B` reçoit un initialiseur par défaut (les quatre
  `missing field 'B' initializer` de clang), `oxm` et `kbn_` inutilisés sont
  retirés. *Non fait* : `Tessellation.cpp:40 polyArea`, dont la suppression n'est
  pas sûre — une lambda homonyme la masque plus bas et trois sites l'appellent
  hors de la portée de la lambda. À reprendre sous clang, seul compilateur qui
  lève ces avertissements (MSVC n'en émet aucun).

#### Lot C — `dfhplus` : la prémisse ne tient plus

Le guide écarte `dfhplus` parce que `thermobench dfhplus` rendrait **ÉCHEC** avec
dix violations (pire cas 8,88·10⁻²). **Non reproductible au HEAD** : la commande
rend **PASS, zéro violation**, y compris avec `--draws 400 --seed 7`. Les huit
commits qui séparent le tag `g0-0.1.0` du HEAD contiennent la relecture adverse
`3277033` (« un défaut MAJEUR, quatre mineurs »), qui l'a vraisemblablement
corrigé. La variante V3 n'a donc pas lieu d'être codée sur ce motif ; la
saturation de l'endommagement de traction sous compression triaxiale, elle, n'a
pas été recontrôlée.


### Ajouté — matériau PAR PHASE en éléments finis (`mode = fem3d`, 2026-09-06)

Capacité **opt-in** : sans la clé `phases`, le comportement est **strictement inchangé** (le banc
`bench_phases` le prouve octet à octet, voir plus bas). Elle rend possible l'expérience qui manquait :
`matWeibullM` fait varier la **résistance** sur un bloc de raideur **uniforme**, donc sans contraste
élastique aucune contrainte ne pouvait se concentrer et le calcul ne **pouvait pas** localiser. Un
contraste de raideur (quartz 83,1 / feldspath 70 / biotite 29,3 GPa) le peut.

- **Une instance de loi par phase**, même clé `law` pour toutes (`laws_[e.phase]->stress(...)`).
  Imposé par le code : `MatLaw` fige `lam_`, `G_`, `K_`, `adp_`, `kdp_` dans son **constructeur** et
  `stress()` ne reçoit pas de matériau — un simple champ `phase` sur l'élément ne changerait **rien** à
  l'élasticité. Les **clés d'option de loi** (`erodeD`, `dfh*`, `cdp*`, `meridian`…) restent **globales**
  (lues une fois dans le Config) ; les écrire par phase est refusé avec ce message. `lcMax_` reste
  **global** pour toutes les phases : la garde de bande de fissuration est alors **conservative** (elle
  peut refuser trop, jamais laisser passer un snap-back structurel).
- **Deux sources de phase.** (a) `mesh = voronoi` en fem3d : **réutilise la tessellation du FEMDEM**
  (`Tessellation3`, mêmes clés, même graine) ; la seule différence est en aval — le FEMDEM **dédouble**
  les nœuds pour poser ses joints, le continuum prend les sommets **tels quels**, donc maillage à
  **nœuds partagés** et frontière de grain = simple **saut de propriétés**. (b) `mesh = file` +
  `$PhysicalNames` dim 3 : le lecteur Gmsh de fem3d **lisait les ntags puis les jetait** — il garde
  désormais le tag physique ; phase homonyme ou `groupPhase.<groupe> = <phase>`.
  Divergence **assumée** avec fdem3d : un groupe sans phase est une **ERREUR** en fem3d (fdem3d retombe
  sur la phase 0 avec un WARNING) — en continuum la phase EST le matériau, un nom mal tapé effacerait
  silencieusement le contraste élastique, c'est-à-dire l'effet même que l'on mesure.
- **Tout ce qui devient par phase**, sous peine de chiffres faux en silence : masse nodale condensée
  (les **deux** sites), pas de temps critique (`hmin / max_p c_P` — sur `mat_.cP()` le dt aurait été
  29 % trop grand avec du quartz au-dessus d'une fiche à 50 GPa, et le schéma ne serait pas devenu
  instable mais **bruyant**, bruit qui se lit comme de la fissuration), pénalité de contact outil
  (`max_p E`), **impédances de Lysmer** prises sur la phase de l'élément **propriétaire** de la face
  (les deux boucles, boîte et cylindre ; `G` descendu dans les boucles), viscosité de volume
  (`rho_p`, `rho_p c_d,p`, tables précalculées).
- **Gardes ajoutées** (une clé lue et sans effet est le motif interdit n° 1) : `phases` sur
  `mesh = grid` refusé ; `phases` sur `mesh = file` sans `$PhysicalNames` refusé ; clés de **joint**
  (`gb.<a>.<b>.*`, `gbAlpha*`, `gbHeteroFactor`, `groupBond.<A>.<B>`, `contactMu.<phase>`, `groupVel.`,
  `gauge.`) refusées en fem3d avec la raison (nœuds partagés = aucun joint) ; audit complet de la
  famille `phase.<nom>.<propriété>` (quatre fautes distinguées, quatre messages : nom de phase
  inconnu, clé de **loi** écrite par phase, `phase.<nom>.law`, et fiche de phase posée **sans** clé
  `phases` — sans elle `PhaseSet` nomme sa fiche unique « rock » et `phase.rock.E` serait acceptée,
  consommée et parfaitement inerte) ; `groupPhase.<groupe>` hors `mesh = file` refusée ;
  `phaseWeibull` posée sans objet (pas de `matWeibullM`, ou une seule phase) refusée ; **audit de masse** `sum(m)` contre `sum_p rho_p V_p` recalculé indépendamment
  (seuil 1e-9, lève) ; en `mesh = voronoi`, aire des faces extérieures = aire de la boîte, refus des
  faces vues plus de deux fois (elles étaient **avalées en silence**), refus si aucune face n'est
  partagée entre grains distincts ; en `mesh = file` multi-groupes, refus si aucune face conforme
  entre groupes (deux corps disjoints se traverseraient — fem3d n'a **pas** de contact entre corps) ;
  compactage **déterministe** des sommets orphelins de la tessellation (jamais d'épinglage FIXED
  silencieux, qui donnerait des broches fantômes encastrées dans le bloc).
- **`phaseWeibull` (false, nouvelle clé)** : la composition `matWeibullM` × phases est **refusée par
  défaut**. Les deux hétérogénéités se multiplient et un run qui localise ne serait plus attribuable —
  or c'est précisément le constat qui motive tout le chantier.
- **Sorties.** Champs cellulaires `.vtu` `phase`, `grain`, `matE`, `matRho`, `matFt` (propriétés
  **réellement vues** par l'élément, `ftScale` compris : la seule façon de vérifier que le contraste
  est arrivé jusqu'à la **loi**, pas seulement jusqu'à la couleur), activés automatiquement dès qu'il
  y a une microstructure — condition **fausse** pour tous les decks existants, `.vtu` inchangés.
  Console : table des phases, fractions **réalisées vs visées**, bilan érosion / endommagement /
  dissipation **par phase**, qualité des tets. Verdict PASS/FAIL du scénario `tension` **supprimé** en
  multiphase (la résistance du bloc n'est celle d'aucune phase ; un chiffre sous l'intitulé
  « vérification » y serait pire qu'une absence de chiffre). **Aucune colonne** ajoutée à
  `history.csv`.
- **Limites dites d'avance, pas des bugs** : ce n'est **pas** un GBM cohésif (aucun joint, donc aucune
  frontière intrinsèquement faible — la frontière concentre, elle ne s'ouvre pas) ; **verrouillage
  volumique** des tets linéaires sans B-bar, corrélé à la phase si `nu` varie (avertissement imprimé,
  il **sous-estime** le contraste) ; `lc = V0^(1/3)` sur les tets en cône de la tessellation, plats
  surtout aux frontières de grain (`lc` médian et rapport diamètre inscrit / lc imprimés).
- **Fichiers** : `include/rockim/Fem3dSolver.hpp`, `src/Fem3dSolver.cpp`, `src/main.cpp` (deux
  barrières élargies à fem3d), `include/rockim/KeysByMode.hpp` **régénéré**
  (`grainSize`/`grainJitter`/`grainSeeding`/`lloydIters`/`refineLevels`/`vertexMergeFrac` étendus à
  fem3d, `phaseWeibull` ajoutée). **Aucune** modification de `MatLaw`, `MatState`, `Material`,
  `Tessellation3`, `Fdem3dSolver` — c'est ce qui borne le risque, et `dpdfh` n'est pas effleurée.
- **Banc `bench_phases/`** (quelques secondes, aucun run long) : `depouille.py` **9/9** — neutralité
  **octet à octet** (2 puis 3 phases identiques = deck sans `phases`, sur les chemins fichier **et**
  voronoï) avec sa variante qui **doit** échouer ; **borne de Reuss** en forme fermée sur une barre à
  deux couches nommées (E_app 44,158 GPa contre 43,678 attendu, **1,10 %**), le seul contrôle qui
  distingue une capacité qui **marche** d'une capacité qui **affiche** ; commutativité par
  `groupPhase` (−1,48 %) ; CFL sur `c_P` **max** (dt ×0,8497 = √(60/83,1), exact).
  `erreurs.py` **16/16** cas fautifs refusés **avec le bon message** (le test vérifie le **message**,
  pas seulement le code de retour : un refus qui envoie chercher au mauvais endroit coûte une
  demi-journée).
- **Preuves de non-régression.** `tools/bitid.py` **8/8 IDENTIQUE** contre l'ancre w18
  (`OMP_NUM_THREADS = 4`, exe issu d'une **reconstruction complète** objets supprimés, `.hpp`
  modifié) : pics retrouvés exactement 7 916,75 / 6 775,14 / 30 390,2 / 19 541,3 / 205 788 N — les
  trois decks `fem3d` couvrent masse, pénalité, Lysmer cylindre, viscosité de volume et cinématique
  hencky, et `fdem3d_kuru9` (6 corps, 3 phases, `groupBond`) atteste que le FEMDEM n'a pas bougé.
  `tools/verify_suite.py` tier fast (`OMP_NUM_THREADS = 1`) : **48/48 avant** (build propre de `HEAD`)
  et **48/48 après**, test par test et **valeur mesurée par valeur mesurée** — aucun échec de
  plateforme préexistant, aucun ajouté. Rapports dans `bench_phases/`.


### Corrigé — relecture adverse de l'étape 1 DFH+ (2026-09-06)

Relecture indépendante (autre graine, 10⁵ tirages, implémentation Python écrite depuis les formules
annoncées, exe du commit précédent recompilé). **Documentation seulement — aucun changement de
comportement** ; `bitid` 8/8 IDENTIQUE avant et après. Compte rendu : `docs/DFHPLUS_etape1.md` §11.

- **⛔ DÉFAUT MAJEUR — `dfhplus` NE REPREND PAS le périmètre de DP-DFH au-delà de σ₃ ≈ 160 MPa**
  (§11.11), c'est-à-dire qu'elle **échoue au critère d'acceptation de l'étape 1**. Le compte rendu
  concluait « compression et triaxiaux IDENTIQUES » sur σ₃ = 0 / 20 / 50 / 100 MPa : toute la
  campagne est **sous** le seuil où le mécanisme s'allume. Balayage `matpoint triax` jusqu'à
  600 MPa : identique à **+0,00 %** jusqu'à σ₃ = 150 MPa, puis `dfhplus` **sature l'endommagement de
  traction (d_t = 0,9999) sous compression triaxiale** là où `dpdfh` reste **exactement à 0**, et la
  résistance confinée tombe de **−1,8 % (175 MPa) à −20,7 % (600 MPa)**. À σ₃ = 600 MPa l'amorçage a
  lieu sous **σ_lat = −600 MPa**. **Cause établie analytiquement** : `Y_i = G‖ε⁺n_i‖²` est piloté par
  la partie positive de la **déformation** élastique, et `ε^e_lat > 0 ⟺ |σ_ax|/σ₃ > (1−ν)/ν = 2,45`,
  condition vraie sur **toute** la surface DP ; avec `σ^eq = 0,98030·E·ε^e_lat` le seuil de Weibull
  (112,529 MPa) est franchi à **σ₃ = 175,8 MPa prédit**, mesuré entre 150 et 175 MPa (accord 0,2 %).
  C'est le défaut classique des découpages spectraux **en déformation** sous confinement. **La bande
  0–100 MPa de la campagne de confinement est sûre, mais la zone broyée sous percussion travaille à
  460–750 MPa : le régime est atteint dans l'application visée.** L'admissibilité thermodynamique
  n'est PAS en cause (0 dissipation négative, ρψ convexe, σ = ∂ρψ/∂ε) : c'est le **moteur
  d'endommagement à l'intérieur du cadre** qui est à revoir. **`dfhplus` ne doit pas être utilisée
  au-delà de σ₃ ≈ 150 MPa, ni sur un maillage de percussion, avant traitement.**

- **Réfuté et corrigé : la sonde de décharge du `thermobench` A des faux positifs.** La limite (L1)
  affirmait qu'une dissipation négative mesurée « prouve » une violation. `dfhplus` la réfute : la
  sonde flague **11** incréments (pire **−2 033 J/m³**) là où la **même loi** mesurée par son énergie
  libre exposée en a **0** (pire −0,0245 J/m³), tous à D = 0,9999 et 8/11 non coaxiaux. Cause : la
  décharge est de direction **isotrope figée** face à une compliance anisotrope — biais de signe non
  contrôlé, non couvert par l'argument « énergie stockée croissante ». Corrigé dans l'en-tête de
  `src/ThermoBench.cpp` (L1) et `DOCUMENTATION_rockim.md` §3.4. **Les comptes de la sonde sont des
  bornes supérieures** ; les violations de `dpdfh` au-dessus du plafond d'artefact (≈ 2 kJ/m³, soit
  5,4 % de ses lignes retenues, pire −35,5 kJ/m³ = 57 % de G_f/ℓc) restent avérées.
- **Portée du PASS de `dfhplus` précisée** : il vaut **à ψ = 15°**. En écoulement **associé**
  (ψ = β = 51,7°, régime atteint par `dfhPsiVar`), `dfhplus` **échoue au test 5** (6 incréments,
  pire r = 3,17 à temps gelé) — défaut **hérité** du retour DP partagé, `dpdfh` échoue sur les
  **mêmes 6 incréments** (r = 2,73). En contrepartie, l'écrêtage `dfhpPsiClamp` y est enfin
  **falsifié et validé** : désarmé 1 dissipation négative significative (−2 280 J/m³), armé 0.
- **Point ouvert n° 5 clos : ρψ est CONVEXE** (hessienne de Voigt-Mandel SDP sur 3 000 états
  aléatoires + 7 familles dégénérées ; plus petite valeur propre 6,0·10⁶ Pa à saturation totale).
  L'opérateur tangent reste elliptique.
- **Contrainte latérale résiduelle sous-déclarée** : σ₂₂/(Eε) va de −5,6 % à D = 0,5 (seul chiffre
  publié) à **−22,5 % à saturation**, soit −25,3 MPa de confinement parasite au pic de traction.
- **`dfhpVolInteg = min` ne passe pas le banc** (1 / 10⁶, −0,1996 J/m³, 0 significative) : `min(1−D_i)`
  n'est pas dérivable et ρψ n'est que C⁰ en D au changement d'argmin. `harmonic` (défaut) et `none`
  passent. À réserver aux ablations.
- `dpdfh` **confirmée intacte** : source de `DpDfhLaw` identique octet pour octet, `selftest-dpdfh`
  et 52 004 lignes de `matpoint` (4 chemins × 4 confinements) **identiques au bit** entre l'exe
  recompilé à `f631676` et l'exe livré. `rockim_f2` sans aucune trace du chantier.

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
