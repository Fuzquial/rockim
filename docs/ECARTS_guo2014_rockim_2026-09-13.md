# Les écarts entre rockim et le modèle de Guo 2014 / Solidity — liste, jugement, décision (13/09/2026, 14 h)

Demande de Fernando : « va dans la thèse de Guo, fais la liste des écarts, juge si c'est handicapant ou non
de garder l'écart, et ensuite s'il faut copier Solidity ». Sources lues dans le texte : Guo 2014 chapitre 2
entier (`bibliographie/Guo-L-2014-PhD-Thesis_split/texte/02_…txt`, éq. 2.1-2.61, §2.3.1-2.5), §1.4 ; le
code public de Solidity (`Y3Dfd.c` Sigma_tau l. 1078-1300, `S_N_direction` ; `Y3Did.c` contact) ; rockim
tel que configuré par le deck v3-P (`configs/yang2026_bench_s25_v3P.cfg`) et par les bancs B/C/D du 13/09.
« Handicapant » = empêche de reproduire l'impact de Yang 2026 (réaction ~57 kN, bit 5,62 m/s, ~360
éléments pulvérisés, radiales à 10 mm).

## 1 — Le tableau

| # | Élément | Guo 2014 (thèse) | Solidity (code public) | rockim (v3-P → bancs) | Écart | Handicapant ? | Décision |
|---|---|---|---|---|---|---|---|
| 1 | Éléments finis | tétraèdres linéaires 4 nœuds, 1 point d'intégration, déformation constante (§2.3.1) | idem | idem | aucun | — | rien |
| 2 | Loi de volume | néo-hookéen T = μ/J (B − I) + λ/J ln J I + η D, grandes déformations (éq. 2.2-2.8) | CauchyTet4 idem | co-rotationnel petites déformations (défaut) ; `bulkModel = neohookean` = éq. 2.6 disponible (exclusif de `law`) | forme de la loi | **non** pour la roche intacte (ε < 2 %) ; **discutable** pour les fragments écrasés (det F 0,5-0,7 sous l'insert) et pour la mesure δ_m de la pulvérisation, qui lit « leur » ε | option : banc E = C + `bulkModel = neohookean` |
| 3 | Viscosité η D (éq. 2.6) | η « paramètre visqueux », valeur non donnée | négligeable en dur dans le code public ; ARMA 2024 : « mass damping 4 000 (calcaire) / 5 000 (grès) » | `bulkViscosity` = σ_v = 2 μ D co-roté (μ = η/2) ; v3-P : 0 | valeur du granite inconnue ; nature (η de l'éq. 2.6 ou c massique en 1/s) incertaine | **oui si c'est une dissipation majeure** (η = 4 000 = 12 % du critique de Munjiza à 1,4 mm : sensible sur le rebond) | banc C : η = 4 000 → `bulkViscosity = 2000` ; demander à Xiang la nature du « mass damping » |
| 4 | Joints intrinsèques partout, nœuds dédoublés (§2.3.1) | idem | roche : idem ; acier et carbure : `groupContinuum` (continuum EF exact, sans joints) | acier/carbure sans complaisance de joint | **non** (à p0 = 50 E leur complaisance vaut ~2 % ; le nôtre est plus proche de l'acier réel, et c'est ×1,06 sur le pas de temps) | garder |
| 5 | Repère du joint : plan moyen N16-N25-N34, ouverture = projection sur z', glissement dans le plan (éq. 2.11-2.22) | idem | normale de facette courante, ouverture et glissement par paire de nœuds, points milieux d'arête (`jointQuadrature = midedge`) | détail géométrique, identique aux petites rotations | **non** | garder |
| 6 | Intégration du joint : 3 points milieux d'arête, poids 1/3, f16 = S/6 (σA + σC) (Table 2.2, éq. 2.34-2.36) | idem | `midedge` : At = S/3 par point, traction répartie pour moitié sur chaque paire = exactement éq. 2.34 | aucun | — | rien |
| 7 | Enveloppe de cisaillement f_s = c − σ_n tanφ, cut-off c − ft tanφ en traction (éq. 2.24) | fs = c en traction (clamp), c + tanφ\|σ\| en compression | `jointShearEnvelope = yang` = la thèse ; sous `solidity` = le code | thèse ≠ code en traction | **non** (traction + cisaillement mixte rare ; les radiales sont en mode I pur) | garder les deux |
| 8 | Longueur h de la pénalité = **moyenne des arêtes** du joint (éq. 2.25-2.26) | el = moyenne des trois arêtes de la facette (`S_N_direction`) | **avant** : hmin_ = le sliver (0,226 mm) pour TOUS les joints → pj uniforme, ~240 E équivalent ; **depuis g1y16** : `jointPenaltyLength = edge` | convention | **oui** : ×2,5 sur le pas de temps, raideur de joint dépendante du sliver et non de l'élément | `edge` dans B, C ; à généraliser au deck s = 1 suivant |
| 9 | Valeur de p0 : E ≤ p0 ≤ 10 E (éq. 2.28) | Yang : 3 000 GPa = 50 E (hors plage de Guo, Turon α ≈ 50) | v3-P : pf 20 sur hmin ≈ 240 E ; B/C : pf 25 sous `edge` = 50 E = Yang | Guo vs Yang | **non** en soi (c'est le choix publié de Yang) ; le coût est le pas de temps | garder 50 E pour comparer ; tester 10 E ensuite (dt ×2,2) |
| 10 | δ_c = 3 Gf/f depuis zéro, Gf ≈ ⅓ f δ_c (éq. 2.29-2.30) | ot = max(2 op, 3 Gf/ft) depuis op | `jointDeltaC = solidity` (code) ; `guo` (thèse) disponible | op ≪ ot (0,1 % pour Kuru) | **non** | rien |
| 11 | z-curve a = 0,63, b = 1,8, c = 6 (éq. 2.32) | idem, en dur | `jointSoftening = munjiza` (mêmes constantes) | aucun | — | rien |
| 12 | **D (éq. 2.33) : fonction des δ COURANTS, « 0 otherwise »** — sans mémoire, le joint redevient neuf en se refermant | idem : z recalculé à chaque pas, nfail compté au pas courant | v3-P : D = **max historique** + glissement plastique (règle de Fukuda / Yan 2023) ; g1y16 : `jointShearUnload = solidity` = littéral | **le gros écart** : chez nous un joint endommagé reste endommagé, chez eux il guérit | **oui** : c'est l'hypothèse du lit lâche (2 892 facettes mortes à 173 µs, réaction 20-30 kN) contre le noyau solide (~57 kN) | **banc B tranche** |
| 13 | Compression : σ = (p0/h) δn, linéaire, sans plafond (éq. 2.31 l. 1) | idem | `jointElastic = parabolic` : 2 pj dn | aucun | — | rien |
| 14 | Cisaillement : τ = z f_s avec f_s(σ_n) ; à z = 0 le joint ne transmet plus rien, ni cohésion ni frottement | dpefm = 0 en dur | v3-P : `jointResidualMu = 0` + `jointDeath = damage` (équivalent) ; `solidity` : littéral | aucun (équivalent) | — | rien |
| 15 | Rupture « dès que l'état de contrainte atteint l'enveloppe » (p. 70) | nfail > 1 : deux points sur trois au même pas (l. 1175) | `jointFailRule = majority` = le code | thèse ≠ code | **non** | garder le code |
| 16 | Après rupture : contact par potentiel, éq. 2.49 (Munjiza), détection NBS, Coulomb (§2.3.4) | pénalité = min(pe)/200, kt = 2/7 pe, μ = μd·min(1 − D_i, 1 − D_j), recalage de naissance [0,01 ; 3], aucun amortissement | potentiel de Munjiza (même intégrale), détection propre, `potStiffnessByPhase = min`, `potTangentFactor 1,4286`, `contactDamageCoupling = solidity`, `gcBirth = penalty` ; raideur : 1 E (v3-P) → **0,25 E** (bancs B, C, D) | raideur ×4 avant le 13/09 | **moyennement** : leur contact est plus mou que le nôtre, donc ce n'est pas la cause du déficit de réaction ; mais il fixe la dissipation et la mobilité des fragments | bancs B, C, D à 0,25 E ; D isole cet effet seul |
| 17 | Intégration temporelle : Euler avant sur v puis x avec v_{t+1} (éq. 2.53-2.55) = schéma symplectique | idem | leapfrog / Verlet vitesse : même famille, même ordre | aucun | — | rien |
| 18 | Pas de temps : 0,1 h √(ρ/E), h = plus petite arête (éq. 2.56) ≈ 21 ns à 1 mm ; DEM π/5 √(m/k), k = E h (éq. 2.57-2.60) — **sans les ressorts de joint** | Yang : 2,5 ns (empirique, validé 2025) | `dtFactor = 0,15` × période du ressort nodal le plus raide, **joints compris** : 1,0 ns à s = 1, commandé par le sliver de la roche | notre règle compte les joints et le sliver | **oui pour le coût** (×2,5 sur le pas), non pour la physique | `edge` (fait), mailleur (§3), puis tester dtFactor 0,25 sous le garde-fou d'énergie |
| 19 | Masse nodale concentrée | idem | idem | — | — | rien |
| 20 | Aucun amortissement numérique hors η | idem | `dampingLocal = 0`, `jointXi = 0` | aucun | — | rien |
| 21 | DIF (absent de Guo) | Yang éq. 1-2, réévalué à chaque pas (dpeftdif) | `strainRateDIF = yang`, `strainRateDIFArm = continuous` | aucun | — | rien (coquille de l'exposant 0,07 notée) |
| 22 | Pulvérisation (absente de Guo) | Yang éq. 3-4 dans la version interne : Cd, ε_d, définition de δ_m non publiés | `bulkDamage = yang` avec Cd = 1 et δ_m = h ε_vm déviatorique ; **inerte** dans nos runs (nPulv 0 jusqu'à 180 µs, 2 à 183 µs) | définition de δ_m, Cd, ε_d | **oui potentiellement** (leur noyau à D = 0,9 = ~360 éléments = 40 mm³) mais **en aval** de la réaction : sous 24 kN le continuum reste sous le seuil | après B : si la réaction monte, regarder si nPulv suit ; sinon demander Cd et δ_m |
| 23 | **Maillage** : 1 mm dans la demi-boule R 12,5, 2 mm à R 25, 10 mm au bord, insert 0,7 ; 230 788 tétras ; dt 2,5 ns (Yang) | — | même champ de tailles, mais gmsh rend une **arête médiane de 1,37 mm** dans la boule (14 722 tétras au lieu de ~35 000 pour un vrai 1 mm) ; 120 185 tétras au total | notre « s = 1 » est **1,37× plus grossier** que le leur dans la zone broyée | **oui** : fissures, pulvérisation (nPulv compte des éléments), cratère — tout se lit à l'échelle de l'élément | SR ≈ 0,73 dans `make_impact_mesh.py` (×2,5 tétras dans la boule, dt ×0,73) — à décider avec le coût |
| 24 | Masses des corps : piston 1,173 kg, bit 1,509 kg (Yang 2026) | — | 1,057 et 1,367 kg (−10 % chacun) | géométrie du maillage (facettisation, longueurs) | **moyennement** (−10 % d'énergie et de quantité de mouvement ; le rapport piston/bit est le bon) | caler ρ de l'acier par corps (deux phases acier) ou corriger la géométrie |
| 25 | Sensibilité au maillage (§2.4) : la zone plastique doit couvrir ≥ 3 éléments, sinon sa longueur est celle de l'élément | — | Kuru : l_pz ~ 0,4 E Gf/ft² ≈ 10 mm → 7 éléments à 1,37 mm, 10 à 1 mm | — | **non** | rien |
| 26 | Groupes de contact après rupture (éq. 2.39-2.46) : les tétras des six nœuds du joint | i1elbe marqués à la rupture | `gcActivation = adaptive` (paires nées des joints morts + voisinage) | implémentation | **non** | rien |

## 2 — Ce que dit le tableau

Sur 26 lignes : 13 sans écart ou sans effet (1, 5, 6, 10, 11, 13, 14, 17, 19, 20, 21, 25, 26), 5 écarts non
handicapants à garder (4, 7, 9, 15 et le détail 2 pour la roche intacte), 2 écarts moyens (16 contact,
24 masses), et **6 écarts qui comptent** :

- **12, la mémoire du joint** — le seul écart de physique des joints qui reste : c'est l'hypothèse du lit
  lâche, testée par le banc B ;
- **8 et 18, la longueur de pénalité et le pas de temps** — corrigés par `jointPenaltyLength = edge`, le
  reste est le mailleur (§3) ;
- **3, l'amortissement** — valeur et nature inconnues pour le granite, banc C ;
- **23, la résolution** — notre s = 1 est un s ≈ 1,37 ; rien ne le compensera ;
- **22, la pulvérisation** — nos paramètres sont leurs Table 1, mais δ_m et Cd sont à eux.

## 3 — Faut-il copier Solidity ?

Copier veut dire deux choses.

**Copier leur physique** : les trois écarts qui comptent et qui se copient (12 la loi réversible, 16 le
contact à p0/200, 3 l'amortissement η D) sont copiés depuis le 13/09, en clés opt-in bit-identiques, et
les bancs B, C, D les mesurent un par un contre le témoin A. Les autres écarts sont nuls, non handicapants,
ou inconnus (22). Il n'y a **plus rien à copier** dans le code public : ce qu'il contient est transcrit,
et ce qui manque (DIF interne, pulvérisation, deck du granite, amortissement) n'y est pas.

**Copier leur code** (compiler Solidity et lui donner notre maillage) : possible, LGPL, format Y3D ; mais
sans DIF, sans pulvérisation et sans deck, il reproduirait au mieux le banc B. Ça ne vaut que comme
contrôle croisé de rockim, pas comme outil.

**Décision proposée** : ne pas copier davantage avant les bancs. Si B ramène la réaction vers 45-57 kN, la
loi réversible était la cause, et on garde `solidity` comme loi de comparaison (pas comme loi de rockim :
elle n'a pas de potentiel en mode mixte et dissipe Gf d'un coup à la rupture, le budget B4 le mesure). Si B
ne bouge pas, l'écart n'est pas dans la loi de joint : on passe à 23 (résolution, SR 0,73), à 24 (masses),
puis à 22 (pulvérisation, en demandant δ_m et Cd).

## 4 — Le mailleur (chasse aux slivers, mesure du 13/09)

`tools/mesh_quality.py` (diamètre inscrit h = 6 V / Σ aires, par corps) sur `impact_yang_s1_pose.msh` :
les 25 pires tétras de la roche (h 0,226-0,30 mm, arêtes 0,8-1,4 mm, h/arête 0,15-0,21 contre 0,41 pour un
tétra régulier) ont **deux nœuds sur la face z = 0 et deux à 0,7-1,1 mm sous elle** — la première couche
sous la surface ; les tétras intérieurs ont h ≥ 0,323 mm. Netgen et Relocate3D ne déplacent pas les nœuds
de surface : `opt = 2` était sans effet (0,2267). MMG3D n'est pas dans ce gmsh (4.15.2 : « MISSING DATA »).

| variante (s = 1) | roche h min [mm] | < 0,3 mm | insert h min | commentaire |
|---|---|---|---|---|
| B (Delaunay + Netgen, actuel) | 0,2265 | 25 | 0,205 | référence |
| H10 : HXT + `Mesh.OptimizeThreshold` 0,5 | **0,270** | **5** | 0,250 | +19 %, meilleure |
| H10b : HXT + seuil 0,7 + 2 passes gmsh + Netgen | 0,247 | 14 | 0,214 | l'optimiseur ré-abîme |
| T7 : Delaunay + seuil 0,7 + 3 passes + Netgen | 0,258 | 16 | 0,222 | +14 % |
| F4R : frontal + Relocate2D + Netgen | (en cours) | | | |

Le gain accessible avec gmsh est de l'ordre de +20 % sur h_min (dt +9 % sous `edge`, +20 % sous `min`) :
les slivers de surface ne partent pas sans remailleur de qualité. Le vrai levier du pas de temps est 8
(`edge`, ×2,2) ; et la résolution (23) coûtera dt ×0,73 quoi qu'il arrive.

## 5 — Relecture des deux diagnostics indépendants du 12/09 (lus le 13/09 à 15 h)

`docs/DIAGNOSTIC_ICL_independant_2026-09-12.md` et `docs/COMPLEMENT_YANG_sources_et_corrections_2026-09-12.md`
(auteur non identifié, non commités). Ce qu'ils apportent, vérifié ici :

- **Filtre des facettes rompues (§5 du diagnostic) : vrai, corrigé.** `imp_lib.broken()` prenait `damage ≥ 0,999`,
  or sous `jointFailRule = majority` `damage` est le MAX des trois points et vaut 1 dès qu'UN point cède ; la
  facette ne rompt qu'à deux points sur trois (`tBreak ≥ 0`). Vérifié sur la trame 18 du s = 1 : 3 737 contre
  **3 357** rompues, 380 faux positifs (10 %). `imp_lib.broken()` lit désormais `tBreak` (repli sur l'ancien
  filtre si le champ manque). Les figures et comptes de ce document et du rapport du matin sont donc ~10 % trop
  hauts ; les conclusions (cône compact, rien au-delà de 9 mm) ne changent pas.
- **Étiquette traction/cisaillement à la rupture : vrai.** Le moteur `plastic` normalise par `slipRef` (pression
  courante) et l'étiquette exportée par `J.slipF` : les proportions rouge/jaune sont indicatives, pas à comparer
  quantitativement à leur fig. 14. À corriger dans l'instrumentation (sortie seule).
- **Branche élastique de cisaillement sous `plastic` : vrai.** Elle est linéaire de pente pj (`tauTr = pj (s − s_p)`)
  là où Guo éq. 2.31 transposée au cisaillement donne la parabole de pente initiale 2 pj ; le glissement au pic
  s_p = f_s/pj est le même. `jointElastic = parabolic` ne touchait que le mode I. La branche `solidity`
  transcrit la parabole. Ligne 13 du tableau à lire ainsi : compression identique, cisaillement pré-pic différent.
- **Pulvérisation : vrai et précisé.** rockim mesure δ_m = h_e ε_vm avec h_e = **diamètre inscrit** (médiane
  0,51 mm dans la boule), pas l'arête (1,37 mm) : le seuil δ_0 = 14 µm vaut donc 2,7 % de déformation
  déviatorique, et D = 0,9 est atteint à κ = 106 µm (forme rationnelle de Camanho), soit 21 %, et non à 400 µm.
  En compression isotrope dev(ε) = 0 : le modèle ne s'arme jamais. Ligne 22 complétée.
- **3 000 GPa n'est pas la pénalité du Kuru : vrai.** C'est la valeur ARMA 2024 du calcaire St Anne (E = 57 GPa,
  soit 52,6 E) ; le granite n'est pas publié. Les bancs B/B2 l'emploient comme proxy et le disent.
- **« Mass Damping Coefficient » 4 000 : unités et opérateur inconnus, vrai.** Le banc C (η D, 2 000 dans la
  convention 2 μ D) est exploratoire, pas une transposition établie. À demander avec le deck.
- **det F ≈ 0,96 dans la zone centrale à 180 µs : accepté.** Ma mention « det F 0,5-0,7 sous l'insert » (ligne 2)
  venait d'un commentaire du code, pas d'une mesure de ce run : retirée.
- **`meanTensionCapFactor = 3` (cap caché à 32,9 MPa sur la pression moyenne en traction) : vrai, actif dans
  tous les runs.** Sous joints intrinsèques à ft = 11 MPa (DIF ≤ 1,85) un élément porte au plus ~20 MPa de
  traction principale : le cap ne devrait pas s'armer, mais les VTU n'exportent pas la contrainte des éléments
  (champ `velocity` seul) : invérifiable a posteriori. À poser à 0 dans tout deck de conformité et à instrumenter
  (compteur d'écrêtage dans le journal).
- **Vitesse d'indentation : estimateur différent, vrai.** Yang 2025 prend la pente de la portion linéaire du
  déplacement ; avec cet estimateur (10-90 % de l'enfoncement, 65-169 µs) le s = 1 donne **6,86 m/s** (insert)
  au lieu du maximum instantané 7,37. L'écart à 5,62 reste (+22 %), moins qu'annoncé (+31 %).
- **Jumeau s = 2,5 : piston 26 % plus léger, vrai.** La comparaison s = 1 / s = 2,5 n'est pas une étude de
  résolution ; le mailleur sait figer le train (argument SR distinct de s) : à faire pour la prochaine série.
- **Banc B = un paquet (loi + longueur + facteur + contact + naissance) : vrai.** D isole le contact ;
  **B1** (loi seule, `yang2026_bench_s25_law.cfg`) et **B2** (pénalité seule, `_pen.cfg`) ajoutés à la file
  (`tools/queue_bench_E.sh`, après D).
- **St Anne d'abord (impact 2025 sans pulvérisation, paramètres ARMA publiés) : d'accord.** C'est le cas
  discriminant du moteur joint/contact, sans les inconnues de l'éq. 3-4. À monter après la série s = 2,5.
- **Guo 2020 p. 40 : les fractures formées restent permanentes** : cohérent avec la transcription (guérison
  pré-rupture seulement, mort à la rupture).

Ce qui ne change pas : la chaîne lit lâche → réaction faible → bit rapide → ni radiales ni pulvérisation reste
la lecture des données ; les deux documents ne la contestent pas, ils demandent qu'elle soit attribuée par des
sensibilités séparées (fait : A, B, B1, B2, C, D) et mesurée avec les bons filtres (fait).

## 6 — Résultats des bancs s = 2,5 à 300 µs (13/09, 14 h - 17 h 15, rockim_g1y16, 14 fils)

Maillage `impact_yang_s2.5_pose.msh` (piston 0,777 kg = 31,5 J, train 1,158 kg), témoin = deck v3-P.
Réaction = dérivée de la quantité de mouvement piston + train (`tools/fig_fp.py`) ; v_ind = pente 10-90 % du
déplacement de l'insert (estimateur de Yang 2025, `tools/fig_kinetics.py` T2) ; rompus = `tBreak ≥ 0`.

| Banc | Ce qui change par rapport à A | t_fin | p [mm] | F_roche max [kN] | v_ind pente [m/s] | rompus | eJnt [J] | verdict |
|---|---|---|---|---|---|---|---|---|
| A témoin | — | 300 µs | 1,18 | 53,5 | 5,20 | 375 | −2,2 | pas de retournement à 300 µs (bit 2,45 m/s) |
| D | contact de Solidity seul (0,25 E, `gcBirth = penalty`) | 300 µs | 1,22 | 50,7 | 5,54 | 424 | −3,0 | **sans effet** : le contact n'est pas la cause |
| B2 | pénalité seule (`edge`, facteur 25 = 50 E) | 300 µs | 1,12 | 60,7 | 5,06 | 424 | −3,2 | léger raidissement (+13 % de F), bit 1,9 m/s à 300 µs |
| B1 | **loi `solidity` seule** | **75 µs** | 0,14 | — | — | 2 028 | **+43,5** | **ENERGY ABORT** : 18,4 J créés (KE 49,9 J > 31,5 initiale) |
| B | loi + pénalité + contact | **83 µs** | — | — | — | 2 020 | **+47,4** | **ENERGY ABORT** : 22,2 J créés (70 % de l'énergie du piston) |
| C | B + viscosité η D (`bulkViscosity = 2000`) | **92 µs** | — | — | — | 1 560 | **+37,6** | **ENERGY ABORT** : 4,4 J créés malgré la viscosité |

Lecture :
- **La loi de joint de Solidity, transcrite mot à mot, crée de l'énergie** dès la première vague de ruptures
  (B1 : le poste joints passe de −0,03 J à +43 J entre 74 et 82 µs pendant que 2 → 2 028 joints cèdent, soit
  ~20 mJ par joint). B1 l'attribue à la loi seule : ni le contact (D sain), ni la pénalité (B2 sain), ni le
  relais `penalty` (D sain) n'en sont responsables. Mécanisme : τ = z·f_s(σ_n) sans mémoire — pas de potentiel
  en mode mixte ; sous l'insert f_s atteint ~1 GPa (c + 1,85 × 500 MPa) sur un glissement au pic de ~1 µm, et
  un joint comprimé qui glisse sous une contrainte normale oscillante rend au retour plus qu'il n'a reçu à
  l'aller (la pathologie mesurée sur `origin` le 12/09 : 16 J en 81 µs). Ce n'est pas une erreur de
  transcription : c'est Sigma_tau tel qu'écrit dans Y3Dfd.c. Le code public n'a pas de bilan d'énergie.
- **Conséquence pour la comparaison** : la loi réversible ne peut pas être « la cause du noyau solide » au sens
  d'une physique meilleure ; si leur code fait la même chose, une part de leur réaction plus haute (~57 kN) peut
  être de l'énergie créée. Hypothèse à tester : bancs B', B1', C' **sans garde-fou** (`budgetAbortPct = 0`,
  file `tools/queue_bench_F.sh`, en cours) — mesurer la réaction, le bit et les fissures à côté de l'énergie
  créée, comme leur code le ferait.
### 6 bis — la variante sans garde-fou tranche l'hypothèse (13/09, 17 h 15 - 22 h, banc B')

`configs/yang2026_bench_s25_solidity_noabort.cfg` = banc B avec `budgetAbortPct = 0` : la loi fait ce
qu'elle fait chez eux, dont le code n'a aucun bilan d'énergie. Arrêté à la main à 251 µs sur 300 (11 trames
conservées) parce que le verdict est acquis et monotone :

| B' sans garde-fou | 80 µs | 100 µs | 150 µs | 200 µs | 250 µs |
|---|---|---|---|---|---|
| Réaction de la roche [kN] | 21,9 | 11,4 | 5,6 | −3,0 | 1,4 |
| Joints rompus | 842 | 6 246 | 10 145 | 14 133 | 14 240 |
| Énergie créée (poste joints) [J] | 16 | 369 | 1 320 | 2 732 | 2 745 |

**2 745 J créés pour 31,5 J apportés par le piston (×87), la quasi-totalité des joints de la roche rompus,
et la réaction qui s'effondre de 22 kN à ~1 kN.** Conclusion : l'hypothèse « une part de leurs ~57 kN
pourrait être de l'énergie créée » est FAUSSE dans cette forme — ici l'énergie créée détruit la roche et
tue la réaction au lieu de la gonfler. La loi de Solidity transcrite mot à mot n'est pas utilisable comme
loi de rockim et n'explique pas leur noyau solide. Elle reste utile comme point de comparaison contrôlé
(avec garde-fou) et comme preuve que leur formulation, telle que publiée, n'a pas de potentiel en mode mixte.
B1' (loi seule) et C' (avec viscosité) NON lancés : même divergence, déjà mesurée par leurs versions avec
garde-fou (B1 +18 J à 75 µs, C +4,4 J à 92 µs) ; décision de Fernando le 13/09 à 22 h pour libérer la machine.

- **Le témoin ne se retourne pas à 300 µs** (Yang 255 µs, 1,0 mm) ; B2 s'en approche (bit 0,94 m/s à 300 µs,
  estimateur T2). Aucun banc ne montre de radiales (fissures ≤ 16-17 mm de rayon de centroïdes pour des
  facettes de 3,4 mm : mailles trop grosses pour en juger, voir `tools/crack_paths.py` de T1).
