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

### Campagne de correction du 13/09 (après-midi) — `rockim_g1y17.exe`

Cadrage : `docs/CAMPAGNE_correction_2026-09-13.md` (tâches S1-S4 solveur, T1-T4 outils/decks, B build et
ancre, V vérification adversariale, Z synthèse), issu des deux diagnostics indépendants du 12/09 et de
`docs/ECARTS_guo2014_rockim_2026-09-13.md` §5. Les agents ajoutent leurs puces ci-dessous.

- **T2 (deuxième passe) — contrôle croisé des estimateurs de Yang : la pente est-elle une mesure ou un
  artefact de fenêtre ?** (`tools/yang_estimators.py`, `tools/yang_report.py`, `tools/test_yang_estimators.py` ;
  aucune source C++, aucune clé solveur, registre inchangé). La première passe de T2 étant déjà livrée
  (pentes 10-90 %, rebond, pic de la 1re onde), cette passe la **vérifie par des estimateurs indépendants**.
  Trois fonctions ajoutées hors du chemin par défaut : `slope_theilsen` (médiane des pentes de toutes les
  paires de points sur la MÊME fenêtre), `window_sensitivity` (fenêtres 5-95, 10-90, 20-80, 30-70, 40-60 %),
  `wave_transit` / `transit_from_log` (transits 1D du train, **géométrie lue dans le journal du run** :
  `bit : z = [...]`, `piston : z = [...]`, `gauge bit : ... z = [...]`, `groupVel.piston`). Option `--robust`
  de `yang_estimators.py` et de `yang_report.py` ; **sans l'option, sortie inchangée** — le rapport
  `yang_report.py out_yang2026_v3 results/yang2026_v3.log` est octet pour octet celui de la première passe
  (`tests_f2/campagne13/T2/yang_report_v3_T2.log`) et le tableau de `fig_kinetics.py` est identique ligne
  pour ligne. **Mesuré** (`tests_f2/campagne13/T2/test_yang_estimators_passe2.log`, **45 PASS / 0 échec**,
  rc = 0, `OMP_NUM_THREADS = 4`, ~25 s, aucun solveur lancé) : (1) **la jauge est bien à mi-bit** — bande
  `[0,28 ; 0,31]` m (centre 0,295) contre bit `z = [0,17322 ; 0,41502]` m (milieu 0,29412 ; L = 241,8 mm),
  écart 0,9 mm = 0,36 % de L, et le garde-fou du solveur est muet dans les deux journaux ; (2) **Theil-Sen
  égale les moindres carrés** à 1,07 / 2,24 % (s = 1 insert / bit) et 1,07 / 1,20 % (s = 2,5) ; (3)
  **dispersion des fenêtres larges** 2,36 % (insert) et 5,55 % (bit) sur le s = 1, 3,56 % et 4,77 % sur le
  témoin — mais +10,55 % en fenêtre étroite 40-60 % : le chiffre ne se cite pas sans sa fenêtre ; (4) les
  estimateurs naïfs DOIVENT différer et diffèrent : p_max/t_pmax −22,7 à −30,9 %, max instantané +11,9 à
  +46,9 % ; (5) **l'écart à Yang survit au choix de fenêtre** — s = 1 : +22,0 % (insert) et +16,8 % (bit) sur
  10-90 %, encore +12,2 % sur la fenêtre la plus favorable (6,31 m/s contre 5,62) ; le témoin s = 2,5 tombe à
  −0,3 % (5,60 m/s) ; (6) **première onde** : premier signal physiquement possible à la bande 20,40 µs (bord
  haut, c_P = 5 778 m/s), arrivée au centre 26,00 µs (c_barre = 5 048 m/s), retour de la réflexion du bas du
  bit 74,25 µs — le pic retenu (34,27 µs en s = 1 ; 67,99 µs en s = 2,5) est **avant** ce retour, donc la
  contrainte de référence n'est pas une superposition (marge de 6,3 µs seulement pour le témoin), et il vaut
  l'impact 1D ρcv/2 = 178,30 MPa à −1,3 % (175,95) et −3,0 % (173,03) ; la fenêtre détectée (66,0 / 71,7 µs)
  est plus COURTE que le créneau de piston (103,0 µs) : le critère de retombée ne délimite pas le créneau
  complet ; (7) **contre-exemples G** : un déplacement bilinéaire 3 puis 9 m/s viole le critère de dispersion
  (21,74 %) et une jauge muette avant 80 µs viole le critère de première onde (pic à 160 µs) — les critères
  peuvent donc tomber. **Fait nouveau sur le rebond** : le balayage des 44 runs d'impact du dépôt montre
  qu'**aucun n'enregistre le retournement** ; le seul qui le franchit, `out_yang_bench_s25_plastic` (300 µs,
  p_max 0,904 mm à 264,2 µs), n'a récupéré que 2,43 % de p_max — non mesurable au sens de la première passe ;
  forcé (`min_amp = 0,02`) il donne 0,915 m/s et la vitesse instantanée finale est +1,341 m/s, contre 4,65 m/s
  chez Yang. Une mesure conforme (pente après 450 µs) demande 450 µs, soit 2,25 × le témoin (≈ 12 400 s à
  14 fils partagés, ~45 min seul) — **non lancé** (règle 6). NON fait : aucune figure nouvelle (la planche
  `fig_kinetics` de la première passe reste la livraison graphique) ; les vitesses de transit supposent l'acier
  du deck (200 GPa, 7850, 0,29) et une propagation 1D dans un bit cylindrique — les élargissements de section
  (circlip, plaque, épaulement) ne sont pas modélisés dans `wave_transit`.

- **B — build propre, ancre de bit-identité et rejeu des mini-tests S1-S4 avec `rockim_g1y17.exe`**
  (`tools/build.ps1 -Clean -Jobs 8` le 12/09 à 16 h 40, 17 unités recompilées, lien OK ; `build/rockim.exe`
  copié en **`rockim_g1y17.exe`**, sha256 `c5bfcb9b1890c3a661bdfd54a937482ca41e459ff79be7a3377229278fa50caf`,
  2 492 416 octets — même taille que le dernier incrémental S4 `c099f1a2…`, hachage différent = horodatage PE
  du build propre). **Ancre `python tools/bitid.py --exe rockim_g1y17.exe` : 8/8 IDENTIQUE** contre
  `rockim_g1ref.exe` (`results/bitid_g1y17.log`, `.json` ; 4 fils, machine partagée par un banc 14 fils :
  cdp_PQ 105 s pic 7 916,75 N, dpr_T1 255 s pic 6 775,14 N, sk2011 110 s pic 30 390,2 N, kuru9 233 s,
  heilman 86 s pic 19 541,3 N, visc_yan 258 s, toolcontact 179 s pic 205 788 N, ucs_yan 181 s) **+ 1/1
  IDENTIQUE** sur le 9e deck `fdem3d_yang_v2_court` (`--refs tools/bitid_refs_jointlaw.json`, ancre
  `rockim_g1y10.exe`, 376 s, `results/bitid_jointlaw_g1y17.json`). Les clés S1-S4 absentes ne changent donc
  aucun octet des 9 chemins ancrés. Rejeu des mini-tests des quatre tâches avec `rockim_g1y17.exe`
  (`OMP_NUM_THREADS = 4`, scripts, journaux et sorties des vérificateurs dans `tests_f2/campagne13/B/`) :
  **S1** 23 runs rc = 0, `check_s1.py` identique ligne pour ligne à la première passe (hors chemins) — dead ⊂
  tBreak ≥ 0 partout, 149 mortes / 355 rompues (UCS 2D 6 ms), openMax 330/336 ouvertes + 6 fermées
  (heilman), slipRef vs slipF 2/392 étiquettes en 2D et 0/336, 0/200 en 3D, cap 3 ft = 0 écrêté partout,
  cap 1e-9 ft 3D = 148 el. (max 211), excès 22,88 MPa ; les 19 sorties (history.csv + dernier VTU joints)
  sont **octet pour octet** celles de l'exe S1 `cc80c1ef…`. **S2** `check_s2.py` TOUS OK (contact
  piston/bit à 2,234 µs, +7 369 N, max 2,721e5 N à 8,907 µs, antisymétrie 0 exact, Σ Fc = grpFz à 0,97 N,
  impulsion 2,20393 vs 2,20543 N·s soit 0,08 %, témoin identique hors colonnes Fc_*), deux refus rc = 1 ;
  sorties identiques à l'exe S2 `c9f5953c…` ; les deux runs de 30 µs ont pris 430 s et 374 s sous la charge
  partagée (308 s et 295 s à la première passe). **S3** 23/23 OK (`W_I` = 1,15987 Gf solidity, 1,00095 Gf
  plastic/origin, glissement résiduel 11,98 µm sous plastic, tilt : joint mort au deuxième point), 9
  `jointbench.csv` identiques à l'exe S3 `8ef6e476…`. **S4** 43/43 OK (seuils 1,68995 % / 0,55552 %,
  isotrope : deviatoric D = 0, total/principal 48/48 à D = 0,9 ; témoins `*_ref` identiques à
  `rockim_g1y16.exe`), 23 sorties identiques à l'exe S4 `c099f1a2…`. Aucun exe copié sous un autre nom,
  aucun banc arrêté, aucun `git add`.
- **S4 — pulvérisation : mesures alternatives de δm et essais élémentaires** (`src/Fdem3dSolver.cpp`,
  `include/rockim/Fdem3dSolver.hpp`, miroir 2D `src/FdemSolver.cpp`, `include/rockim/FdemSolver.hpp` ; motif
  DIAGNOSTIC §3 et COMPLEMENT §7 : la longueur de référence et la mesure de déformation de l'éq. 3-4 de Yang 2026
  ne sont pas publiées ; rockim mesurait δm = h_inscrit·εvm, soit un seuil de 2,7 % sur le s = 1 et **rien** en
  compression isotrope). Trois clés opt-in (fdem, fdem3d), toutes **refusées sans `bulkDamage = yang`** (une clé
  inerte et muette est le piège n° 1) : **`bulkDamageLength = inscribed | edge`** (défaut `inscribed` = hEl_ = 6V/A ;
  `edge` = moyenne des six arêtes du tétraèdre, des trois arêtes en 2D, vecteur `hEdge_` rempli sous l'option seule),
  **`bulkDamageStrain = deviatoric | principal | total`** (défaut `deviatoric` = √(2/3)‖dev ε‖ ; `principal` =
  max |εi| via `maxAbsEigSym3`, 2D sur ε1, ε2, εzz = 0 ; `total` = √(2/3)‖ε‖ trace comprise), **`bulkDamageProbe`**
  (bool, sortie seule : seuils en déformation δ0/h et δF/h pour les deux longueurs au démarrage, max δm / max D /
  éléments armés à chaque trame et au résumé, champ VTU `bulkDm`). Dans `elementForces` la ligne historique
  `dm = hEl_[eI] * std::sqrt(2.0 / 3.0) * ed.norm()` est conservée mot pour mot sous `bdLen_ == 0 && bdStrain_ == 0`
  (2D idem) ; les variantes sont dans la branche `else`. Registre des clés régénéré (3 clés `fdem fdem3d`).
  **Mesuré** (build incrémental 4 fils, exe `tests_f2/campagne13/S4/rockim_s4test.exe` sha256 `c099f1a2…`, decks
  `make_decks.py` : boîte fdem3d 4 mm, grid 2×2×2 = 48 tétraèdres de Kuhn d'arête 2 mm, joints incassables, crushCap
  et meanTensionCap neutralisés, Kuru δ0 = 14 µm ; `check_s4.py` → `RESULTATS_check_s4.txt`) : (A) seuils imprimés
  **1,68995 % (inscrit 0,8284 mm) / 0,55552 % (arête 2,5202 mm)** = δ0/h exacts du tétraèdre de Kuhn à 3e-7 près
  (rapport 3,04, l'analogue des 2,7 % / 1,0 % du diagnostic) ; (B) **compression isotrope 2,5 GPa** (pression suiveuse
  sur les six faces, base bloquée en z) : sous `deviatoric` et sous `edge` seul, max δm = 0,81 / 2,47 µm = **εm
  résiduel 9,8e-4** (transitoire de rampe), **D = 0, nPulv = 0, bdWork = 0** ; sous `total`, `principal`, `edge`+`total`
  : 48/48 éléments à D = 0,9, bdWork −9,5 / −11,5 / −1,0 J — la mesure à trace s'arme, la déviatorique jamais ;
  confinement latéral seul et uniaxial (−0,6 m/s, mors libres) : toutes les mesures s'arment ; (C) **bit-identité** :
  les quatre decks témoins sans clé S4 (iso, biax, uni, uni2d) donnent `history.csv` et dernier VTU **IDENTIQUES**
  octet pour octet entre l'exe S4 et `rockim_g1y16.exe` ; avec `bulkDamageProbe` seul, `history.csv` identique et
  VTU = témoin + le seul champ `bulkDm` ; ancre `fdem3d_kuru9_court` (bulkDamage = yang, refs par défaut, ancre g1ref) **IDENTIQUE** 166 s
  (`RESULTATS_bitid.txt` ; les sept autres decks et `fdem3d_yang_v2_court` sont laissés à l'agent Build) ; (D) **recomposition indépendante** de δm par élément depuis la géométrie des VTU
  (F = dx·dX⁻¹, polaire exacte par SVD, ε = U − I, les six combinaisons h × εm) contre le champ `bulkDm` à la
  dernière trame avant armement : la mesure active du deck coïncide à **< 6e-4** (uniaxial 5e-5, biaxial 3e-4 à
  6e-4, isotrope 2e-4 à 6e-4) et **chaque mesure inactive s'en écarte** (3,9 % minimum : `total` contre
  `deviatoric` en uniaxial à ν = 0,25 ; 20 % `principal` ; ×3,04 `edge`) — une mesure mal codée serait vue ;
  (E) refus : `bulkDamageStrain` sans `bulkDamage`, `bulkDamageStrain = vonmises`, `bulkDamageLength = min` → code 1,
  message attendu. Sous pression imposée la chute de raideur fait **sauter** la boîte (iso `principal` : 13,3 µm à la
  trame 6, 259 µm à la trame 8, courant final 139 µm) : `bulkDm` est un max historique, la comparaison se fait
  avant l'armement. **NON fait** : un essai de **cisaillement simple** (aucun scénario fdem3d ne prescrit un
  déplacement tangentiel de mors ; le confinement latéral seul, σxx = σyy = −p, σzz = 0, tient lieu de troisième
  état déviatorique) ; les essais sont sur maillage grid de Kuhn (48 tétraèdres congruents), pas sur un maillage gmsh ;
  aucune des trois mesures n'est la définition des auteurs (COMPLEMENT §7) — ce sont des hypothèses instrumentées.

- **S3 — banc de joint cinématique, solveur 3D, scénario nouveau** (`src/Fdem3dSolver.cpp`,
  `include/rockim/Fdem3dSolver.hpp` ; motif DIAGNOSTIC §6.2 « valider une loi de joint unique sur un petit banc 3D »).
  **`scenario = jointbench`** (fdem3d) : deux tétraèdres réguliers d'arête `jbEdge` partageant une facette (un joint,
  trois points), tétraèdre A fixe, tétraèdre B **déplacé rigidement** en tête de pas le long d'un chemin prescrit —
  **aucun nœud libre**, donc ni masse ni amortissement ni dynamique parasite ; le pas de temps est un pas
  d'échantillonnage (jbAmp/(jbRate·jbSteps)). Clés (toutes fdem3d, refusées hors du scénario) : **`jbMode`** (tension \|
  shear \| mixed), **`jbAmp`** (2e-5 m), **`jbNormal`** (0 m, offset normal posé d'abord par une rampe), **`jbUnloadAt`**
  (0 ; fraction de jbAmp où l'on décharge à 0 puis recharge, alignée sur les pas), **`jbRate`** (1 m/s), **`jbTilt`**
  (0° ; rotation de B proportionnelle à s, axe du plan de la facette à bras distincts), **`jbEdge`** (1e-3 m),
  **`jbSteps`** (4000). Sortie `jointbench.csv` (t, s, phase, dn/ds/σ/τ/D par point, nFail, broken, dead, Fn, Ft) et
  quatre critères falsifiants imprimés au résumé `[JOINTBENCH]` : (i) retrace de la recharge sur la charge
  (échantillons appariés au même dn/ds, sans interpolation), (ii) résiduel à traction nulle sur la décharge, (iii) aire
  sous σ(dn) [τ(ds)] jusqu'à tBreak contre Gf et contre la loi codée, (iv) instants de rupture des points et du joint.
  Refusés : `insertion = adaptive | none`, `jointTSL = camacho`, `mesh`/`W`/`D`/`H`/`nx`/`ny`/`nz`. Deux lignes de relevé
  `if (jbOn_)` dans `processJoint` (branches plastic/origin et solidity), un `jbDrive()` en tête de `step()`, un
  `jbRecord()` avant `integrate()`, un `jbReport()` en tête de `finalize()`, trois retours anticipés (`buildMesh` →
  `jbBuildMesh`, `placeTool` → `toolNone_`, `setupBoundaries`), `jbSetupPath()` après `computeStableDt()` et l'exclusion
  du WARNING « mesh = grid + scénario de fissuration » — tous sous `jbOn_` ; hors du scénario `jbOn_ = false` et rien
  d'autre ne s'exécute (les huit `cfg_.has("jb*")` de refus ne lisent que le deck). Registre des clés régénéré (8 clés `fdem3d`). **Pas de miroir 2D** (scénario, pas clé de loi ;
  le banc vise la loi 3D de Solidity et sa règle nfail > 1 à trois points). **Mesuré** (build incrémental, 4 fils, decks
  `tests_f2/campagne13/S3/` générés par `make_decks.py`, joint Kuru ft 10,98 / c 29,84 MPa / Gf 50 / GfII 1 000, pj =
  25 E/edge = 1,5e15 Pa/m, conventions de Solidity, 0,2-0,3 s par deck ; `check_s3.py` → `RESULTATS_check_s3.txt`) :
  traction 20 µm, décharge à 10 µm — `solidity` : W_I = **1,15987 Gf** (attendu par la loi codée 1,15999, écart
  0,011 %), retrace **5e-10 ft** sur 6 000 paires [PASS] ; `plastic` et `origin` : W_I = **1,00095 Gf** (attendu
  1,00107, 0,012 %), retrace **FAIL** (0,9999 ft : sécante de décharge) — la variante qui doit échouer ; cisaillement
  30 µm sous σ_n = −60 MPa (jbNormal −20 nm), décharge à 12 µm — `plastic` : glissement résiduel **11,98 µm** aux trois
  points (|s_p| interne 18,4 µm en fin de run après la charge inverse), retrace FAIL (3,7 c), W_II = 1,70 GfII
  (frottement compris) ; `solidity` : retrace 8e-11 c sur 4 800 paires [PASS], résiduel ≤ 1e-18 m, W_II = 1,1677 GfII
  (1,1593 + 2/3 f_s²/pj à f_s = 141 MPa) ; `origin` : résiduel ≤ 1e-18 m, retrace FAIL (4,7 c), W_II = 1,0066 GfII ;
  basculement 0,5° (bras aux milieux d'arêtes 0,28 / 0,63 / 0,76 mm) — `solidity` + majority : points rompus à 10,26 /
  10,715 µs / jamais, joint mort à **10,715 µs = le deuxième** [PASS] ; `plastic` + majority : 8,855 / 9,245 / jamais,
  joint à 9,245 µs [PASS] ; `plastic` + **any** : joint à **8,855 µs = le premier** (conforme à `any`), et W_I = 0,935 Gf
  (−6,6 % : le D partagé casse la facette avant que les deux autres points aient dissipé) ; sans basculement les trois
  instants sont confondus (non discriminant, comme prévu). Refus : `jbAmp` en percussion et `mesh = grid` sous
  jointbench → code 1 avec le message attendu. Le recalcul indépendant de W_I depuis `jointbench.csv` égale le chiffre
  du solveur. Ancres bitid avec l'exe final (sha256 `8ef6e476…`, `RESULTATS_bitid.txt`) : `fdem3d_kuru9_court`
  (refs par défaut, ancre g1ref) **IDENTIQUE** 156 s ; `fdem3d_yang_v2_court` (refs jointlaw, ancre g1y10) **IDENTIQUE**
  369 s — les six autres decks de l'ancre sont laissés à l'agent Build (règle 7). Le rejeu des mêmes decks sans le
  scénario ne passe par aucune ligne nouvelle (`jbOn_ = false`).
  NON fait : le trajet à pression **variable** et le relais joint/contact du diagnostic (le banc n'a pas de nœud
  libre) ; la branche élastique (7 nm à pj = 25 E/mm) n'est pas résolue à 4 000 pas (1,5 pas ; l'adoucissement l'est à
  0,01 %) ; aucun banc 2D.
- **S1 — instrumentation de rupture et de contrainte, 3D + 2D** (`src/Fdem3dSolver.cpp`, `src/FdemSolver.cpp`,
  headers ; motif DIAGNOSTIC §5-§6.1, ECARTS §5). Trois ajouts, tous **sortie seule** (aucune force ne les lit) :
  (a) clé **`writeRuptureFields`** (false) → VTU des joints + `dead` (0/1) et `openMax` (m, ouverture normale
  géométrique max sur les points et le temps, sans dn0), VTU des éléments + `pMean` (Pa, tr(σ)/3 assemblée, traction > 0) ;
  (b) clé **`jointBreakModeRef`** = `slipF` (défaut, historique) \| `slipRef` : sous `slipRef` la partition rn/rs de
  l'étiquette de rupture (`failMode`/`breakMode`/`rnB`/`rsB`, branche `plastic`) normalise le glissement plastique par
  la plage COURANTE du moteur (`slipRef`, pression du point sous `jointShearRange = coulomb`), relevée point par point
  après le retour radial, au lieu de `J.slipF` ; inerte sous `origin`/`solidity`/`camacho` ; (c) **compteur du cap**
  `meanTensionCapFactor` : quand il est actif (défaut 3 !), le démarrage l'annonce et chaque trame imprime le nombre
  d'éléments écrêtés au dernier pas, le max sur un pas depuis la trame précédente et l'excès max (MPa) — même condition,
  même valeur écrêtée qu'avant (compté par fil, réduit après la boucle). Registre des clés régénéré (`fdem fdem3d`).
  **Mesuré** (build incrémental sha256 `cc80c1ef…`, 4 fils, machine partagée ; decks et vérificateur dans
  `tests_f2/campagne13/S1/`, résultats `RESULTATS_*.txt`) : `tools/bitid.py` **4/4 IDENTIQUE** à l'ancre g1ref
  (`fdem_ucs_yan_adaptive_court` 56 s, `fdem3d_kuru9_court` 172 s, `fdem3d_cut3d_heilman_court` 67 s pic 19 541 N,
  `fdem3d_visc_yan_3d` 182 s). Clés armées sans `coulomb` : `history.csv`/`frames.csv` byte-identiques et 0 valeur de
  `damage`/`tBreak`/`breakMode`/`vonMises` changée sur UCS 2D 6 ms (355 rompues), kuru9 3D, heilman 3D (origin, 2
  rompues), visc_yan 3D (200 rompues) → **0 étiquette changée quand slipRef ≡ slipF, comme prévu**. `dead ⊂ tBreak ≥ 0` :
  0 violation partout (UCS 2D : 149 mortes / 355 rompues ; heilman/visc : 0 morte, `jointDeath = separation`).
  `openMax` : UCS 2D 355/355 rompues ouvertes (max 23 µm) ; **heilman plastic + coulomb : 330/336 ouvertes, 6 rompues
  mais FERMÉES** — la distinction demandée par le diagnostic ; visc 200/200. `slipRef` vs `slipF` sous `coulomb` :
  UCS 2D (392 rompues) **2 étiquettes changées** (traction 194 → 192, cisaillement 198 → 200), `tBreak`/`damage`/
  `vonMises`/`frames.csv` identiques, `history.csv` ne diffère QUE par les colonnes de recensement `nBrokTen`/`nBrokShear`
  (133 lignes) ; heilman 3D (336 rompues, 335 cisaillement) et visc 3D (200 rompues, toutes traction) : 0 changée,
  `history.csv` byte-identique. Compteur du cap : à 3 ft, 0 élément écrêté sur tous ces decks ; falsifiant 2D
  `meanTensionCapFactor = 0.01` → max 2 éléments sur un pas, excès 0,41 MPa ; 3D kuru9 `0.01` → 0 (pm max 377 Pa dans
  la roche à 10 µs), `1e-9` → 148 el. au dernier pas de la trame 1 (max 211), 113 à la trame 2, excès max 22,9 MPa.
  NON fait : le test sur le deck s = 2,5 `T = 5e-6` du cadrage (aucune rupture à 5 µs : l'onde n'atteint pas la roche ;
  remplacé par les quatre decks bitid ci-dessus qui rompent) ; aucune relecture de la trame 18 du s = 1 (les nouveaux
  champs exigent un nouveau run). ⚠️ un exe copié dans `%TEMP%` est refusé par Apex One (WinError 5, 11 essais) : les
  mini-tests ont tourné sur une copie de `build/rockim.exe` placée dans le dépôt, supprimée ensuite.
- **S2 — force de contact entre corps nommés, solveur 3D** (`src/Fdem3dSolver.cpp`, `include/rockim/Fdem3dSolver.hpp` ;
  motif DIAGNOSTIC §6.1 « exporter les forces de contact piston/bit et insert/roche »). Clé opt-in
  **`contactForcePairs`** (—, fdem3d) = `a:b c:d ...` : trois colonnes `Fc_<a>_<b>_x/y/z` par paire dans `history.csv`
  (après les jauges `szz_*`, avant `eEl`), somme **au pas courant** des forces de contact général — normale + tangentielle,
  potentiel de Munjiza (phase C) ET pénalité nœud-face — exercées par a sur b, en N. Sommée dans les boucles d'assemblage,
  qui sont **série par construction** (l'annonce du cadrage « réduction OpenMP par fil » ne s'applique pas : la phase
  parallèle du potentiel ne calcule que la géométrie, l'assemblage reste dans l'ordre canonique des paires — donc le
  même chiffre quel que soit le nombre de fils, sans atomique). Table `fcIdx_` (nGroups²) → indice de paire, `fcAccum(eLo,
  eHi, F)` inline, remise à zéro en tête de `generalContact` comme `grpF_`. Refusés : `a:a` (somme nulle par
  construction), corps inconnu, forme sans `:`, maillage sans corps nommés ; exclus : outil analytique, joints vivants
  (`groupBond`). **Pas de miroir 2D** : `FdemSolver` n'a pas de corps nommés (`elemGroup_`) — le registre des clés
  (régénéré, `fdem3d`) refuse la clé en 2D. **Mesuré** (build incrémental sha256 `c9f5953c…`, 4 fils, machine partagée ;
  decks, script et vérificateur dans `tests_f2/campagne13/S2/`, résultats `RESULTATS_*`) : deck s = 2,5 `T = 3e-5`,
  `frames = 1`, 8 paires + `trackGroup = bit` (10 610 pas, 308 s ; témoin sans la clé 295 s ; débit 15 ms/pas mesuré
  d'abord sur 1 µs). (C1) `Fc_insert_rock` et `Fc_rock_insert` = **0 exact** sur les 2 123 lignes (l'onde n'a pas atteint
  la roche) ; (C2) premier contact piston/bit à **t = 2,234 µs** (0,02 mm / 9 m/s = 2,2 µs, et non « ~6,7 µs » écrit au
  cadrage), `Fc_bit_piston_z` = +7 369 N puis max **272 kN à 8,9 µs**, jamais négative (1 965 lignes > 0, 0 < 0) ; (C3)
  `Fc_piston_bit + Fc_bit_piston` = **0 exact** (x, y, z) ; (C4) `grpFz` (V2/B2) = Σ `Fc_<X>_bit_z` sur les cinq partenaires
  du bit à l'arrondi 6 chiffres de `history.csv` près (écart max 0,97 N pour une tolérance d'arrondi de 1 N ; un terme
  manquant vaudrait ~1e5 N) — la plaque porte sur le bit (1,5–3 kN) ; (C5) estimateur quantité de mouvement de
  `tools/fig_fp.py` appliqué au piston : ∫`Fc_bit_piston_z` dt = 2,2039 N·s contre m_p Δv_p = 2,2054 N·s, **écart 0,08 %**
  (m_p = 0,7767 kg lue au journal) ; (C6) sans la clé : `history.csv` **identique** ligne à ligne une fois les colonnes
  `Fc_*` retirées (2 124 lignes), 4 VTU sha256 identiques, journaux identiques hors annonce, temps mur et comptage des
  clés (compteurs `potential stats` égaux : 146 351 300 paires). Variantes qui DOIVENT échouer : `rock:rock` → code 1
  « nulle par construction », `insert:foo` → code 1 « corps 'foo' inconnu ». Ancre `tools/bitid.py --refs
  tools/bitid_refs_jointlaw.json --only fdem3d_yang_v2_court` : **1/1 IDENTIQUE** (324 s) contre `rockim_g1y10.exe`.
  NON fait : la réaction de la roche à 200 µs sur le témoin `out_yang_bench_s25_v3P` — aucun redémarrage n'existe dans
  le solveur (grep `restart`/`checkpoint` : rien), un rejeu de 50 µs en régime fracturé dépasserait de loin les 5 min ;
  la comparaison à l'estimateur de `fig_fp.py` est donc faite sur le piston à 30 µs. Les deux runs de 30 µs ont pris
  ~5 min chacun (machine partagée par les bancs), à la limite de la règle 8.
- **S2bis — le canal insert/roche est prouvé vivant, et le biais de l'estimateur `fig_fp.py` est chiffré**
  (aucune ligne de C++ changée : `contactForcePairs` était déjà implémentée et ancrée ; ajouts = deck
  `tests_f2/campagne13/S2/s25_fcrock.cfg`, `run_s2bis.sh`, `check_s2bis.py`, résultats `RESULTATS_check_s2bis.txt`,
  `RESULTATS_s25_fcrock.log`, `RESULTATS_runs_s2bis.log` ; exe `rockim_g1y17.exe`, `OMP_NUM_THREADS = 4`). Motif : C1
  ci-dessus mesure `Fc_insert_rock` = 0 **exact** — une colonne toujours nulle ne prouve pas qu'elle est branchée, et la
  réaction de la roche restait sans chiffre. Banc d'**instrumentation** (et non de physique : la vitesse n'est pas celle
  de l'article) : même maillage s = 2,5 et même loi, train lancé à −20 m/s (`groupVel.bit/insert/circlip`), `T = 4e-6`,
  `frames = 1`, `trackGroup = insert`, 6 paires autour de l'insert — 79 s, 1 416 lignes. (C7a) `Fc_insert_rock_z`
  devient non nulle à **t = 1,009 µs** et vaut **2 150 N au maximum (t = 4,0 µs)** ; **1 059 lignes non nulles, 1 059
  négatives** (l'insert pousse la roche vers le bas : le signe attendu, partout) ; (C7b) `Fc_insert_rock +
  Fc_rock_insert` = **0 exact** sur x, y, z ; (C7c) Σ `Fc_<X>_insert_z` sur les cinq partenaires possibles de l'insert
  = `grpFz` (`trackGroup = insert`) à **0 N d'écart, exactement**, sur les 1 416 lignes, la somme n'étant pas triviale
  (max |grpFz| = 2 150 N). (C8) sur le banc 30 µs de S2, les colonnes `Fc` permettent enfin de confronter l'estimateur
  « quantité de mouvement » de `tools/fig_fp.py` aux forces **extérieures mesurées** du système {piston, bit, insert,
  circlip} : la seule force extérieure y est le contact **plaque/bit, jusqu'à 3 214 N** (roche 0 N, poids −19 N) ;
  l'estimateur **exact** (Σ m_g v_z,g par corps, circlip non suivi) referme le bilan à **1 990 N max / 1 052 N RMS,
  soit 0,73 % / 0,39 %** de l'échelle (272 kN) ; l'estimateur **publié** de `fig_fp.py` (train supposé rigide,
  m_train dv_bit/dt) s'en écarte de **13,0 kN max / 6,3 kN RMS = 4,8 % / 2,3 %** : c'est le prix mesuré de l'hypothèse
  de corps rigide. Contrôle FALSIFIANT : système faux (piston exclu, sa poussée devient extérieure) → résidu 96,8 kN,
  **49 ×** celui du système correct. (C9) sur le témoin 200 µs `out_yang_bench_s25_v3P` (antérieur à la clé, donc sans
  colonnes `Fc`), les deux estimateurs donnent **65,1 kN à 193,2 µs** (`fig_fp.py`) contre **43,9 kN à 197,8 µs**
  (exact, 3 corps) : **48 % d'écart** entre deux lectures de « la réaction de la roche », à surveiller avant toute
  comparaison à Yang. NON fait : la confrontation directe des colonnes `Fc_insert_rock` à cet estimateur sur 200 µs —
  pas de redémarrage dans le solveur et un rejeu ≥ 50 µs dépasse le budget de 5 min (règle 6). Build incrémental
  `tools/build.ps1 -Jobs 4` : « ninja: no work to do », `build/rockim.exe` sha256 `f6b07967…`.
- **T3 — série de maillages à train figé et masses du train** (`docs/MAILLAGE_serie_2026-09-13.md`).
  Mesuré d'abord : le train DÉPENDAIT de SR (`Mesh.MeshSizeExtendFromBoundary = 0` → les faces latérales
  et l'intérieur de l'insert, du bit et du piston suivent le champ de la roche) : piston 1 033 tétras /
  1,057 kg à SR = 1 contre 209 / 0,777 kg à SR = 2,5 pour le même s — c'est le « piston 26,5 % plus léger »
  du diagnostic du 12/09 §4. `tools/make_impact_mesh.py` : arguments nommés opt-in **`train=fixed`** (train
  maillé seul, roche cachée, sous son champ à l'échelle s ; roche dans un second modèle à l'échelle SR ;
  fusion v2.2 écrite par le script), `srfar=` (échelle du seul champ lointain), `verbose=` ; défaut
  bit-identique (sha256 `a0047b76…` du s = 2,5 régénéré). Deux schémas en un modèle réfutés en chemin
  (`Restrict` : le train varie avec SR par l'ordre des tirages aléatoires ; passes par visibilité :
  `generate(d)` efface les faces du groupe caché). Vérifié : train identique (nœuds, surface, tétras, masses)
  à SR = 2,5 / 1 / 5 sous Delaunay et HXT ; lu par `rockim_g1y16.exe` (2 µs, code 0, masses = outil).
  `tools/mesh_quality.py` : `--ball R` (reproduit 14 722 tétras / 1,372 mm), `--masses` (reproduit les six
  masses du résumé du solveur sur trois maillages), `--yang`, `--rho`. Série générée (s = 1, hxt) :
  `meshes/impact_yang_train1_rock{25,137,073}_hxt.msh` — roche 6 221 / 90 878 / 233 671 tétras, arête
  médiane dans la boule 3,59 / 1,42 / 1,03 mm, h min 0,765 / 0,270 / 0,200 mm ; train commun 17 789 tétras ;
  la roche du `rock137` est identique à celle de `impact_yang_s1_pose_hxt05.msh`. Masses : la
  facettisation coûte −6,0 % au piston et −4,0 % au bit, et les cylindres exacts du dessin pèsent déjà moins
  que Yang (piston 1,126 contre 1,173 kg, bit 1,342 — 1,421 avec insert et circlip — contre 1,509) ;
  correction proposée = phases `steelPiston` (ρ 8 703) / `steelBit` (8 715 ou 9 195 selon le périmètre du
  1,509 kg), de préférence (ρ, E) scalés ensemble. Aucun run lancé.

- **T4 — decks de conformité et cas St Anne** (aucun run > 2 µs, aucune clé de solveur ajoutée, registre
  inchangé ; motif COMPLEMENT §2, DIAGNOSTIC §4) : `configs/stanne2025_bench_s25_visc0.cfg` (impact insert
  unique sur le calcaire St Anne, Yang 2025 Table 4 : ρ 2731, E 57 GPa, ν 0,31, ft 7,0, c 18,8, φ 45°, GI 12,
  GII 800, glissement 0,6 ; pénalité ARMA 24-0952 p0 = 3 000 GPa → `jointPenaltyLength = edge` + facteur
  26,316 = p0/(2E) ; `bulkDamage = off` (donc `contactDamageCoupling` retirée, le code l'exige avec
  `bulkDamage = yang`) ; DIF 2025 ; `meanTensionCapFactor = 0` ; 10,66 m/s ; 450 µs ; maillage T3 à train
  figé) et `_visc.cfg` (+ `bulkViscosity = 2000`, exploratoire, amortissement ARMA sans unité) ; série
  `configs/yang2026_bench_s25_v4_{A,B,B1,B2,C,D}.cfg` = les bancs v3 + `meanTensionCapFactor = 0` +
  `jointBreakModeRef = slipRef` (S1, **exige g1y17**). Outils : `tools/make_conformity_decks.py` (génération
  idempotente, liste des clés différant du témoin calculée par diff), `tools/deck_smoke.py` (fumée
  T/frames/meshFile surchargés, dépouillement, coût). Mesuré (g1y16, 4 fils, `tests_f2/campagne13/T4/`) :
  St Anne 2/2 démarrent sur le maillage s = 2,5 historique (T3 absent au moment du test), dt 4,165 ns →
  108 052 pas pour 450 µs ≈ 2,3 h partagé / 31 min seul ; v4 tels qu'écrits 6/6 **refusés** sur la seule clé
  `jointBreakModeRef` (garde `unknownKeys = error`, voulu), 6/6 démarrent une fois la ligne retirée ; dt 2,828 ns
  (A, B1, D : 106 098 pas / 300 µs ≈ 2,3 h) et 4,079 ns (B, B2, C : 73 556 pas ≈ 1,6 h). Tableau pour
  validation : `docs/DECKS_conformite_2026-09-13.md`. La ligne `meanTensionCapFactor` de la doc (défaut réel 3,
  compteur d'écrêtage) a été corrigée par S1 pendant T4 ; T4 n'y ajoute que le §9 (outils).

- **T2 — estimateurs de Yang dans les figures** (`tools/yang_estimators.py` nouveau, `tools/fig_kinetics.py`,
  `tools/yang_report.py`, test `tools/test_yang_estimators.py` ; motif COMPLEMENT_YANG §5) : vitesse
  d'indentation = pente de p(t) entre 10 % et 90 % de l'enfoncement maximal ; vitesse de rebond = pente
  de la portion remontante (10-90 % de l'amplitude récupérée) après le retournement, **non mesurable**
  déclaré si le retournement n'est pas dans l'enregistrement ou < 5 % de p_max récupéré ; contrainte de
  référence = pic de la PREMIÈRE onde à mi-bit (départ > 5 % du max global, fin sous 10 % du max courant).
  Les max instantanés historiques restent imprimés à côté, avec les fenêtres ; fenêtres hachurées et
  droites ajustées sur la figure (`--body insert|bit`). Test falsifiant : pentes synthétiques 6 et 4 m/s
  retrouvées à 1e-12, max instantané 7,4 ≠ pente 6,0 sur le même signal, variante naïve « max global »
  se trompe d'onde (200 ≠ 160 MPa) comme attendu, troncatures à 210/216 µs → non mesurable, à 300 µs →
  4 m/s. Mesuré : s = 1 `out_yang2026_v3` v_ind pente **6,86 m/s** (insert, 64,8-169,4 µs, r² 0,998),
  6,57 m/s (bit) contre 7,37 max (bit) ; jauge 175,95 MPa à 34,3 µs = pic de la 1re onde (21-87 µs) ;
  rebond non mesurable (bit −5,73 m/s à 183 µs). Témoin s = 2,5 : 6,01 / 5,60 m/s en pente, 7,61 / 6,27
  en max, 173,0 MPa à 68 µs. `yang_report.py` n'imprime plus une vitesse négative comme « rebond ».
  Sorties : `tests_f2/campagne13/T2/`. Aucune clé solveur (registre inchangé), aucune source C++ touchée.

- **T1 — `tools/crack_paths.py` : fissures connectées et cratère** (outil Python, aucune clé, aucun
  effet solveur ; DIAGNOSTIC §5). Composantes connexes des facettes rompues (`tBreak ≥ 0`) par arête
  partagée après unification des nœuds dupliqués (coordonnées de la trame 0) ; noyau = plus grande
  composante touchant r < 6 mm, bras du noyau au-delà, orientation par |n·e_r| < 0,35 **et** |n_z| < 0,6
  (radiale) / |n_z| ≥ 0,6 (horizontale) / conique ; longueur de fissure radiale de Yang = pointe du plus
  long bras radial connecté au noyau ; rayon de cratère = frontière des facettes de surface (< 1 mm) du
  noyau, r max et moyenne par 12 secteurs ; géométrie de référence par défaut. Figure PDF + PNG
  (vue de dessus une couleur par composante, coupe par l'axe et la pointe radiale), CSV. Cadrage
  **ajusté aux données** (`--lim 0` / `--depth 0` par défaut, 15 % de marge, plancher 10 mm : les
  ±20 mm fixes laissaient 60 % de la figure vide) et titre portant les **deux** rayons de cratère
  (r max et moyenne par secteur) : les CSV sont inchangés, seule la figure l'est (vérifié).
  `imp_lib.broken_mask()` ajouté (même règle que `broken()`, qui l'appelle : sortie inchangée) et
  avertissement dans `imp_lib.metrics()` — `radial` / `crater` sont des étendues spatiales, pas les
  métriques de Yang (9,16 mm contre 8,26 mm de radiale connectée à la trame 18) : valeurs rendues
  inchangées (mesuré : `radial` 9,1627, `crater` 9,1627, `depth` 10,1449 mm, n = 3 357), seul le
  commentaire est nouveau ; le second défaut de lecture du DIAGNOSTIC §5 est ainsi signalé sur place.
  **Mesuré** — s = 1 trame 18 (180 µs) : 3 357 facettes, **22 composantes**, noyau 3 328 (99,1 %),
  17 centrales secondaires + 4 périphériques de 1 à 5 facettes, 16 bras dont 4 radiaux ; **Yang =
  8,26 mm** (bras de 3 facettes, azimut 20°) contre 9,16 mm par l'ancienne métrique (r max des
  centroïdes) ; le plus long bras (84 facettes, pointe 9,66 mm) a |n_z| = 0,60 : pas une radiale ;
  cratère r max 8,71 mm, moyenne par secteur 7,23 mm (min 5,26, couverture 100 %), profondeur 10,74 mm.
  Témoin s = 2,5 à 200 µs : 294 facettes, 3 composantes, noyau 292, 6 bras **tous coniques** (|n·e_r|
  0,38-0,80) : Yang = 0 ; cratère r max 14,01 mm, moyenne par secteur 7,89 mm (ancienne métrique
  11,97 mm). Ce r max est un **point aberrant du maillage grossier**, mesuré : 2 facettes seulement
  sur les 105 de surface l'atteignent à 0,5 mm près, elles partagent le même sommet extérieur et
  leurs aires (4,52 et 4,74 mm²) valent 1,5 × la médiane ; le rapport moyenne par secteur / r max
  vaut **0,563** contre **0,830** au s = 1 (12 facettes sur 885 dans les 0,5 mm du max, sans sommet
  commun) — le r max seul n'est pas une mesure de cratère, la moyenne par secteur l'est.
  Tests : `--selftest` (population synthétique, variante « raccord retiré » qui fait tomber Yang de
  16,5 à 0 mm et le noyau de 16 à 4 facettes, égalités à 1e-9, **17/17 OK**) ; `--check` (partition
  scipy = union-find indépendant) OK sur les deux runs ; `tests_f2/campagne13/T1/verif_T1.py`
  (**10/10 OK**, 20 s) : `broken_mask()` est un refactor pur — masque identique au code en ligne
  d'avant et `broken()` rend `pts[con[masque]]` à 0 exactement sur les deux VTU —, variante
  falsifiante du filtre `damage ≥ 0,999` du 12/09 (3 737 contre 3 357 au s = 1, +11,3 % ; 346 contre
  294 au s = 2,5, +17,7 %) avec inclusion `tBreak ⊂ damage` (0 facette hors), et les deux critères
  de cratère ci-dessus. Rejeu bit-identique des CSV (`tests_f2/campagne13/T1/verif_s1*`,
  `verif_s25*` contre `results/fig/`, 0 ligne différente) ; 26 s et 6 s machine partagée par les
  bancs (10,4 s et 3,0 s machine libre). Réserve mesurée : `CMU Serif` n'est pas installée sur cette
  machine, le PDF sort en **STIXGeneral** (mathtext `cm`) — repli commun à toutes les figures du
  dépôt, convention `font.serif` inchangée.
  Sorties `results/fig/crack_paths_out_yang2026_v3*` et `crack_paths_out_yang_bench_s25_v3P*`.

### Corrigé et ajouté — S1/S2 fermés, `jbMode = cycle` (`rockim_g1y18.exe`, 13/09 nuit)

Suite de la campagne, après la critique indépendante du 13/09 (points repris dans
`docs/ECARTS_guo2014_rockim_2026-09-13.md` §7 et §8).

- **S1 clos** (le correcteur avait appliqué ses quatre points avant l'arrêt de la campagne) : le journal du
  compteur d'écrêtage (bannière + ligne par trame) et le champ `e.pm` sont désormais **sous
  `writeRuptureFields`** — sans la clé, journal et calcul textuellement inchangés ; les deux affirmations de
  documentation contredites par la mesure sont corrigées (`jointBreakModeRef` change bien les colonnes
  `nBrokTen`/`nBrokShear` de `history.csv` ; `pMean` en 2D est un hybride, σzz n'étant ni écrêtée ni
  multipliée par la pulvérisation).
- **S2 clos** : `contactForcePairs` n'était analysée que dans `buildMeshFile()` et restait **silencieusement
  inerte** en `mesh = grid`/`voronoi` et sous `scenario = jointbench` (le piège des clés inertes du dépôt).
  Refus explicite ajouté dans `buildMesh()`, avant tout maillage.
- **`jbMode = cycle` + `jbNormal2` + `jbCycles`** (fdem3d) : cycle **fermé** compression-glissement à pression
  variable, le test que la critique réclamait et que `jbMode = shear` ne fait pas. Critère (v) de `jbReport` :
  travail net par cycle, cumul, rapport W/amplitude, avec la convergence en Δt par `jbSteps`. Neuf decks
  `tests_f2/campagne13/S3bis/cyc*_{plastic,origin,solidity}[_dt2,_dt4].cfg`.
  **Verdict mesuré** : `solidity` crée +8,92e-9 J/cycle (élastique, 3,3 % de l'amplitude) et +1,37e-5 J/cycle
  (endommagé, 21 %), **convergé en Δt** ; `plastic` et `origin`+ratchet donnent un résidu exactement ∝ Δt qui
  s'annule. C'est aussi la première validation sur cycle fermé du correctif `jointSecantRatchet` du 12/09.
- Ancre de bit-identité de `rockim_g1y18.exe` : `results/bitid_g1y18.log`.

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
