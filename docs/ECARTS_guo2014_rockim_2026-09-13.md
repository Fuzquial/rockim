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
