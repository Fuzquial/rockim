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
| 9 | Valeur de p0 : E ≤ p0 ≤ 10 E (éq. 2.28) | Yang : 3 000 GPa = 50 E (hors plage de Guo, Turon α ≈ 50) | v3-P : pf 20 sur hmin ≈ 240 E ; B/C : pf 25 sous `edge` = 50 E = Yang | Guo vs Yang | **non** en soi ; mais **la pénalité du Kuru n'est PAS publiée** : 3 000 GPa est la valeur ARMA 2024 du calcaire St Anne (52,6 E). 50 E est donc un **choix de comparaison** de notre part, pas « la pénalité de Yang » (correction de la relecture V) | garder 50 E pour comparer ; tester 10 E ensuite (dt ×2,2) |
| 10 | δ_c = 3 Gf/f depuis zéro, Gf ≈ ⅓ f δ_c (éq. 2.29-2.30) | ot = max(2 op, 3 Gf/ft) depuis op | `jointDeltaC = solidity` (code) ; `guo` (thèse) disponible | op ≪ ot (0,1 % pour Kuru) | **non** | rien |
| 11 | z-curve a = 0,63, b = 1,8, c = 6 (éq. 2.32) | idem, en dur | `jointSoftening = munjiza` (mêmes constantes) | aucun | — | rien |
| 12 | **D (éq. 2.33) : fonction des δ COURANTS, « 0 otherwise »** — sans mémoire, le joint redevient neuf en se refermant | idem : z recalculé à chaque pas, nfail compté au pas courant | v3-P : D = **max historique** + glissement plastique (règle de Fukuda / Yan 2023) ; g1y16 : `jointShearUnload = solidity` = littéral | **le plus visible**, mais pas le seul : la réponse tangentielle AVANT le pic diffère aussi (branche linéaire de pente pj sous `plastic` contre la parabole de pente initiale 2 pj, voir §5) | **oui** : c'est l'hypothèse du lit lâche (2 892 facettes mortes à 173 µs, réaction 20-30 kN) contre le noyau solide (~57 kN) | **banc B tranche** |
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
| 23 | **Maillage** : 1 mm dans la demi-boule R 12,5, 2 mm à R 25, 10 mm au bord, insert 0,7 ; 230 788 tétras ; dt 2,5 ns (Yang) | — | même champ de tailles, mais gmsh rend une **arête médiane de 1,37 mm** dans la boule (14 722 tétras au lieu de ~35 000 pour un vrai 1 mm) ; 120 185 tétras au total | notre « s = 1 » est **1,37× plus grossier** que le leur dans la zone broyée | **à vérifier** (et non « handicapant » tranché) : tout se lit à l'échelle de l'élément — fissures, pulvérisation (nPulv compte des éléments), cratère. Mais une taille *annoncée* de 1 mm chez Yang n'établit pas une arête médiane de 1 mm, la série T3 garde un insert à 1,43 mm d'arête pour 0,7 annoncés, et « sept éléments dans la zone de processus » ne remplace pas une convergence de la fissuration (relecture V) | SR ≈ 0,73 dans `make_impact_mesh.py` (×2,5 tétras dans la boule, dt ×0,73) — à décider avec le coût |
| 24 | Masses des corps : piston 1,173 kg, bit 1,509 kg (Yang 2026) | — | 1,057 et 1,367 kg (−10 % chacun) | géométrie du maillage (facettisation, longueurs) | **moyennement** (−10 % d'énergie et de quantité de mouvement ; le rapport piston/bit est le bon) | caler ρ de l'acier par corps (deux phases acier) ou corriger la géométrie |
| 25 | Sensibilité au maillage (§2.4) : la zone plastique doit couvrir ≥ 3 éléments, sinon sa longueur est celle de l'élément | — | Kuru : l_pz ~ 0,4 E Gf/ft² ≈ 10 mm → 7 éléments à 1,37 mm, 10 à 1 mm | — | **à vérifier**, comme la ligne 23 : la longueur de zone de processus est une estimation analytique, pas une convergence mesurée (relecture V) | convergence de la fissuration à train figé, série T3 |
| 26 | Groupes de contact après rupture (éq. 2.39-2.46) : les tétras des six nœuds du joint | i1elbe marqués à la rupture | `gcActivation = adaptive` (paires nées des joints morts + voisinage) | implémentation | **non** | rien |

## 2 — Ce que dit le tableau (réécrit le 14/09 avec les mesures des §6 et §7)

Sur 26 lignes : 13 sans écart ou sans effet (1, 5, 6, 10, 11, 13, 14, 17, 19, 20, 21, 26), 4 écarts non
handicapants à garder (4, 7, 9, 15), 2 écarts moyens (16 contact, 24 masses), et **6 écarts qui comptent** —
avec, pour chacun, ce que les bancs du 13 et 14/09 en ont fait :

- **12, la mémoire du joint** — c'était « le gros écart » et l'hypothèse du lit lâche. **Tranché, mais pas
  dans le sens attendu** : la loi réversible de Solidity, transcrite mot à mot, **crée de l'énergie** sur un
  cycle fermé à pression variable (§7, +1,9e-8 J/cycle à D = 0, convergé en Δt) et désintègre la roche quand
  on la laisse faire (§6 bis). Ce n'est donc pas la correction à apporter. Écart voisin, reconnu au §5 : la
  **réponse tangentielle avant le pic** diffère aussi (linéaire de pente pj contre parabole de pente 2 pj).
- **8 et 18, la longueur de pénalité et le pas de temps** — corrigés par `jointPenaltyLength = edge`
  (pas de temps ×1,67 sur le maillage de la nuit) ; le banc B2 mesure que ce seul changement ne crée pas
  d'énergie et raidit la réaction de +13 %.
- **3, l'amortissement** — valeur et nature toujours inconnues pour le granite ; le banc C est exploratoire
  et ne conclut pas (sa loi de joint créait de l'énergie par ailleurs).
- **23 et 25, la résolution** — **à vérifier**, et non tranché : notre s = 1 a une arête médiane de 1,37 mm,
  mais une taille *annoncée* de 1 mm chez Yang n'établit pas la leur, et aucune convergence de la fissuration
  n'a été mesurée. C'est l'objet de la série à train figé (T3).
- **22, la pulvérisation** — nos paramètres sont leur Table 1, mais δ_m, Cd et ε_d sont à eux. Son inactivité
  peut venir de la **mesure de δ_m** (diamètre inscrit 0,51 mm, seuil = 2,7 % de déformation déviatorique)
  **indépendamment** de la réaction trop faible : les deux explications ne sont pas départagées.

**Statut de la chaîne causale.** « Lit trop lâche → réaction faible → ni pulvérisation ni radiales » reste une
**hypothèse cohérente avec les données**, pas une causalité démontrée. Les mesures qui la départageront sont
la force insert/roche directe (`contactForcePairs`), son impulsion, la pénétration et le retournement.

## 3 — Faut-il copier Solidity ? (réécrit le 14/09)

**Copier leur physique** : les trois écarts qui se copiaient (12 la loi réversible, 16 le contact à p0/200,
3 l'amortissement η D) l'ont été, en clés opt-in bit-identiques, et les bancs A, B, B1, B2, C, D les ont
mesurés un par un. **Verdict : copier la branche réversible n'est pas la solution** — elle crée de l'énergie
(§7), et ni le contact (D) ni la pénalité (B2) ne comblent le déficit de réaction.

Ce qu'il **reste** à copier : rien dans le code public, dont tout ce qui concerne le joint et le contact est
transcrit. Mais l'**équivalence complète n'est pas démontrée** pour autant, et ne peut pas l'être depuis ce
code : le DIF interne, la pulvérisation (éq. 3-4), le deck du granite et l'amortissement n'y sont pas. Dire
« il n'y a plus rien à copier » était donc trop fort ; la formulation juste est : **le code public ne contient
plus rien que nous n'ayons transcrit, et ce qui manque n'y est pas**.

**Copier leur code** (compiler Solidity et lui donner notre maillage) : possible, LGPL, format Y3D. Sans DIF,
sans pulvérisation et sans deck, il reproduirait au mieux le banc B — donc une divergence énergétique que leur
code ne mesure pas. Utile comme contrôle croisé, pas comme outil.

**Décision arrêtée le 14/09.** (1) `plastic` reste le **témoin de travail** ; (2) `origin` + ratchet reste une
**comparaison distincte** ; (3) `solidity` est **écartée des campagnes longues** et conservée comme point de
comparaison sous garde-fou et comme **test de régression** (les neuf decks du cycle) ; (4) la suite est la
**force insert/roche directe** pour localiser le déficit mécanique, puis **St Anne sans pulvérisation**
(tous ses paramètres sont publiés), puis la **convergence de maillage à train identique**.

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
| D | contact de Solidity seul (0,25 E, `gcBirth = penalty`) | 300 µs | 1,22 | 50,7 | 5,54 | 424 | −3,0 | **ne resout pas** : F max −5 %, v_ind +6,5 %, ruptures +13 % par rapport a A. Ce changement de contact ne corrige pas le deficit sur ce maillage avec `plastic` ; il n exclut ni une interaction avec une autre loi ni un defaut de transmission apres rupture (correction de la relecture V) |
| B2 | pénalité seule (`edge`, facteur 25 = 50 E) | 300 µs | 1,12 | 60,7 | 5,06 | 424 | −3,2 | léger raidissement (+13 % de F), bit 1,9 m/s à 300 µs |
| B1 | **loi `solidity` seule** | **75 µs** | 0,14 | — | — | 2 028 | **+43,5** | **ENERGY ABORT** : 18,4 J créés (KE 49,9 J > 31,5 initiale) |
| B | loi + pénalité + contact | **83 µs** | — | — | — | 2 020 | **+47,4** | **ENERGY ABORT** : 22,2 J créés (70 % de l'énergie du piston) |
| C | B + viscosité η D (`bulkViscosity = 2000`) | **92 µs** | — | — | — | 1 560 | **+37,6** | **ENERGY ABORT** : 4,4 J créés malgré la viscosité |

**Attention a la base de comparaison** (relecture V) : les ~57 kN de Yang sont une estimation de FREINAGE MOYEN deduite de la cinematique, alors que la colonne `F_roche max` donne des MAXIMA d un estimateur `m dv/dt` qui suppose le train RIGIDE — un maximum de 53,5 ou 60,7 kN ne prouve donc pas que la reaction est reproduite. Biais mesure par S2 : la plaque porte jusqu a 3 214 N sur le bit, soit 13,0 kN d erreur maximale (4,8 % de l echelle). La sortie `contactForcePairs` (force insert/roche mesuree dans le contact) devient la mesure principale a partir du run Kuru de la nuit.

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
qu'ecrit le CODE PUBLIC, qui n'a aucun bilan d'energie (leur version interne n'est pas connue : le §7 dit pourquoi on ne peut pas l'affirmer pour Yang). Arrêté à la main à **261,06 µs** (dernière ligne du CSV ; « 251 » dans la première rédaction était la ligne lue en cours de run) sur 300 (11 trames
conservées) parce que le verdict est acquis et monotone :

| B' sans garde-fou | 80 µs | 100 µs | 150 µs | 200 µs | 250 µs |
|---|---|---|---|---|---|
| Réaction de la roche [kN] | 21,9 | 11,4 | 5,6 | −3,0 | 1,4 |
| Joints rompus | 842 | 6 246 | 10 145 | 14 133 | 14 240 |
| Énergie créée (poste joints) [J] | 16 | 369 | 1 320 | 2 732 | 2 745 |

**2 745 J de travail CUMULE des tractions de joint (poste `eJnt` du bilan B4 — ni l énergie cinétique finale, ni le défaut global du bilan, qui valait 22,15 J pour B et 18,40 J pour B1 à l instant du garde-fou) pour 31,5 J apportés par le piston (×87), la quasi-totalité des joints de la roche rompus,
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

## 7 — Le banc de joint tranche la question de la « pompe » (S3 bis, 13-14/09, `rockim_g1y19`)

La critique indépendante refusait la conclusion « ce n'est pas une erreur de transcription » : le banc B1 isole
bien la branche, mais celle-ci change à la fois la mémoire, la réponse tangentielle avant le pic et les
conventions d'adoucissement, et elle interagit avec le DIF, la rupture et le contact. Le test demandé était
précis : **un cycle fermé compression-glissement à pression variable, sans rupture, sans contact, DIF
désactivé, dans la même routine de joint, puis une convergence en Δt**. `jbMode = shear` ne le fait pas : il
glisse à pression CONSTANTE.

**`jbMode = cycle`** (clés `jbNormal2`, `jbCycles`) parcourt le rectangle fermé du plan (s, δn) : glissement
**dans le plan** 0 → jbAmp à δn = jbNormal ; compression jbNormal → jbNormal2 à s = jbAmp ; retour du
glissement à δn = jbNormal2 ; décompression à s = 0. L'état final EST l'état initial, aux mêmes échantillons
(le chemin est paramétré par le numéro de pas). Deux tétraèdres prescrits : aucun contact, aucun amortisseur,
aucun DIF, aucune rupture. Le travail NET des tractions du joint sur un cycle doit donc être ≤ 0.

**Défaut trouvé et corrigé en cours de route (14/09)** : dans la première version, la direction de glissement
du mode `cycle` tombait dans le `else` de `mixed` (45°), si bien que la branche « glissement » ouvrait le joint
en traction jusqu'à dnE et l'endommageait (D = 1,1e-4). Le tableau ci-dessous est celui de la version corrigée,
glissement **dans le plan**, où **D vaut exactement 0** — c'est la condition que la critique posait pour que
l'argument soit net.

Cycle **strictement non endommagé** (jbAmp = 1,5e-8 m < s_E = 2,7e-8 m ; **D = 0 mesuré**) —
`tests_f2/campagne13/S3bis/cyce_*.cfg` :

| loi | W net du dernier cycle [J] | à Δt/2 | à Δt/4 | W/amplitude | limite Δt → 0 |
|---|---|---|---|---|---|
| `plastic` + ratchet | +2,4844e-10 | +1,2422e-10 | +6,2110e-11 | 7,5e-4 | **0** (W ∝ Δt au bit près) |
| `origin` + ratchet | +2,4844e-10 | +1,2422e-10 | +6,2110e-11 | 7,5e-4 | **0** (identique : sous le cap les deux branches sont la même élastique) |
| **`solidity`** | **+1,9286e-8** | **+1,9138e-8** | **+1,9064e-8** | **4,4e-2** | **+1,899e-8 ≠ 0** (Richardson : les deux extrapolations donnent 1,8990e-8) |

Cycle **endommagé** (jbAmp = 2e-6 m ; D 0,47 à 0,88) — `cyc_*.cfg` :

| loi | W net du dernier cycle [J] | à Δt/2 | à Δt/4 | W/amplitude | limite Δt → 0 |
|---|---|---|---|---|---|
| `plastic` + ratchet | −8,110e-5 | −8,123e-5 | −8,129e-5 | −0,210 | **−8,13e-5 : dissipatif, convergé** |
| `origin` + ratchet | +2,504e-7 | +1,252e-7 | +6,261e-8 | 1,8e-3 | **0** (W ∝ Δt ; cumul −6,0e-5 J, dissipatif) |
| **`solidity`** | **+5,2731e-4** | **+5,2714e-4** | **+5,2705e-4** | **+0,306** | **+5,270e-4 ≠ 0 : création convergée** |

**Ce que cela établit.** (1) La création d'énergie sous `solidity` **n'est pas un artefact de discrétisation** :
elle converge vers une valeur non nulle quand Δt → 0, dans les deux régimes, et vaut 4,4 % (cycle à D = 0) à
31 % (cycle endommagé) du travail échangé. Le couplage pression-cisaillement en est la seule source possible
dans ce montage. (2) `plastic` **passe ces essais** : résidu exactement proportionnel à Δt, donc nul à la
limite, et cycle endommagé franchement dissipatif. (3) `origin` **avec le ratchet passe ces essais aussi** —
c'est la première validation sur cycle fermé du correctif `jointSecantRatchet` du 12/09.

**Ce que cela n'établit pas.** (a) Que `plastic` ou `origin` soient « sains » en général : ces cycles valident
un comportement précis sur un trajet précis. Ils ne disent rien sous DIF variable, sous rotation, en trajet
mixte, ni à la transition joint → contact. (b) Que ce mécanisme explique **quantitativement** les milliers de
joules du banc B' : un joint de 1 mm sur un cycle de 15 nm ne s'extrapole pas à 20 000 joints sous un insert.
(c) Que la version **interne** de Yang présente la même divergence : leur code public n'a pas de bilan
d'énergie, donc personne ne l'aurait vue, mais ce n'est pas une preuve qu'elle y est.

**Décision qui en découle** : `solidity` est **écartée des campagnes longues de reproduction** et conservée
comme (i) point de comparaison contrôlé sous garde-fou et (ii) **test de régression** — les neuf decks du cycle
sont l'essai falsifiant à rejouer après toute modification de la loi de joint. `plastic` reste le témoin de
travail, `origin` + ratchet la comparaison distincte.

## 8 — Corrections apportées à ce document après la critique indépendante du 13/09

- §6 bis : le banc B' s'est arrêté à **261,06 µs** (dernière ligne du CSV), non 251 ; et les 2 745 J sont le
  **travail cumulé des tractions de joint** (poste `eJnt` du bilan B4), pas l'énergie cinétique finale ni le
  défaut global du bilan — le défaut mesuré au garde-fou valait 22,15 J pour B et 18,40 J pour B1.
- §6 : « le contact n'est pas la cause » était trop fort. D donne, par rapport à A : force maximale −5 %,
  vitesse d'indentation +6,5 %, ruptures +13 %. La conclusion tenable est : **ce changement de contact ne
  résout pas le problème sur ce maillage avec la loi `plastic`** ; il n'exclut ni une interaction avec une autre
  loi, ni un défaut de transmission des efforts après rupture.
- §6 : les ~57 kN de Yang sont une estimation de **freinage moyen** déduite de la cinématique, alors que le
  tableau donne des **maxima** ; 53,5 ou 60,7 kN de maximum ne prouvent donc pas que la réaction est reproduite.
  De plus `tools/fig_fp.py` assimile la vitesse de tout le train à celle du bit, ce qui est faux pendant la
  propagation des ondes (biais mesuré par S2 : la plaque porte jusqu'à 3 214 N sur le bit, soit 13,0 kN
  d'erreur maximale sur l'estimateur). **La sortie `contactForcePairs` devient la mesure principale**, contrôlée
  par l'impulsion et les quantités de mouvement de chaque corps.
- §1 ligne 23 : la résolution reste « **à vérifier** », pas « non handicapante ». Une taille annoncée de 1 mm
  chez Yang n'établit pas une arête médiane de 1 mm, la série T3 garde un insert à 1,43 mm d'arête médiane pour
  0,7 mm annoncés, et « sept éléments dans la zone de processus » ne remplace pas une convergence de la
  fissuration.

## 9 — Second tour de la relecture indépendante (14/09) : ce qui a été repris

La relecture lève sa réserve sur l'existence d'une création d'énergie constitutive sous `solidity`, et pose
trois réserves nouvelles. Toutes les trois sont traitées :

1. **« `plastic` est sain » → « passe ces essais ».** Formulation corrigée au §7 : ces cycles valident un
   comportement sur un trajet précis, pas la loi entière (DIF variable, rotation, trajet mixte, transition
   joint → contact restent non couverts). La demande « s'appuyer sur le cas strictement non endommagé et
   vérifier que son dommage reste nul » a révélé un **défaut de mon banc** : la direction de glissement du
   mode `cycle` tombait dans le `else` de `mixed` (45°), donc le « glissement » ouvrait le joint en traction
   jusqu'à dnE (D = 1,1e-4). Corrigé (glissement dans le plan, `rockim_g1y19`) ; les tableaux du §7 sont
   ceux de la version corrigée, avec **D = 0 mesuré**, et la conclusion en sort renforcée : +1,899e-8 J par
   cycle sous `solidity` contre un résidu ∝ Δt sous `plastic` et `origin`.
2. **Contradictions entre le début et la fin du document.** Les cinq passages cités sont réécrits sur place :
   ligne 9 (la pénalité du Kuru n'est pas publiée, 50 E est notre choix de comparaison), ligne 12 (« le plus
   visible, pas le seul » ; la réponse tangentielle pré-pic diffère aussi), ligne 25 (« à vérifier », cohérent
   avec la ligne 23), §2 et §3 entièrement réécrits, §6 bis (« ce que le code public écrit », et non « ce
   qu'elle fait chez eux »).
3. **La cause de l'absence de localisation reste à départager.** Dit explicitement au §2 : la chaîne « lit
   lâche → réaction faible → ni pulvérisation ni radiales » est une hypothèse cohérente, pas une causalité
   démontrée ; la mesure de δ_m peut empêcher la pulvérisation indépendamment, et la résolution peut agir sur
   la localisation.

Les quatre décisions techniques de la relecture sont adoptées telles quelles et inscrites au §3 : témoin
`plastic`, comparaison `origin` + ratchet, `solidity` hors campagnes longues mais gardée en régression,
exploitation de la force insert/roche directe, puis St Anne, puis convergence de maillage à train identique.
