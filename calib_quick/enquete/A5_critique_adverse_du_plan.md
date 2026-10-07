# A5 — Critique adverse du plan de calibration triaxiale Red Bohus (rockim_f2, 2026-09-02)

Périmètre lu : `calib_quick/README.md`, `DOCUMENTATION_rockim.md` (§5.2–5.4, 5.7, 5.16, 6.1, 8), `rockim_f1/calib_triax3d/PLAN_GBM.md` et `CAMPAGNE.md` (§0–2, 3ter), `rockim_p1/calibration_redbohus/README.md`, les decks `q1/q2/q3/q1t*`, les sorties `out_q1_homog_P050`, `out_q1t050_P050`, `out_q2_weibull_P050`, `out_q3_gbm_P050` (history.csv, fdem_final_joints.csv, `_log_q1t050.txt`), le code `src/FdemSolver.cpp` (mors, σ, amortissement, colonnes), les données `experimental_data_red_bohus_clean.json`, `seuils_sbm_bohus.json`, `targets_triax_bohus.json`. Aucun run lancé (règle : validation avant tout lancement) ; tous les chiffres « modèle » ci-dessous sont recalculés à partir des sorties existantes.

---

## 0. Les cinq constats qui invalident le plan en l'état

| # | Constat | Preuve | Gravité |
|---|---|---|---|
| A | **Les pics du README sont surestimés de 33 MPa** (σ₃ = 50) : `q = sigma − offset` avec offset = 16,67 MPa = ν/(1−ν)·σ₃ (mors bloqués pendant `pullDelay`), alors que le déviateur expérimental est σ₁ − σ₃ = sigma − 50. Cas 1 vaut **676 MPa (+13 %)**, non 709 (+18 %) ; q1t050 vaut **514 (−14 %)** ; GBM **480 (−20 %)**. Le facteur de résistance « ≈ 0,7 » devient ≈ 0,73–0,76 et ε_pic recule de 0,04 %. | `out_q1_homog_P050/history.csv` : sigma = 16,67 MPa pour t ∈ [1e-4 ; 3e-4] s ; `FdemSolver.cpp:7256-7262` (vg = 0 tant que t ≤ pullDelay → u_y = 0), `:4156` (sigma = \|gripFy\|/(W·thk)) ; `CAMPAGNE.md` §2 assume explicitement cet offset « pour rendre q comparable » — c'est faux : l'état avant déviateur n'est pas isotrope (σ_yy = 16,7 ≠ 50). | 3 |
| B | **Le post-pic n'est un observable ni côté essai ni côté modèle.** Essai : « chute 26 % » = moyenne de **0,427 / 0,330 / 0,021** (trois fins d'essai arbitraires) ; à σ₃ = 50, `eps_end_common` = 9 410 µε < `eps_peak` = 9 489 µε → n_rep = 2 au pic. Modèle : dans `_log_q1t050.txt`, **Cundall dissipe 2 757 J/m contre 1 368 J/m pour la cohésion** (contact 2 723 dont frottement 912) : le post-pic est à 80 % numérique (amortissement + restitution de contact). | `experimental_data_red_bohus_clean.json` (`chute_fraction` 2_5/2_6/2_7) ; `targets_triax_bohus.json` conf 50 ; `_log_q1t050.txt` lignes « energy budget ». | 3 |
| C | **Le pic est une propriété du maillage, pas des joints.** Première insertion = exactement l'enveloppe MC de la sonde sur l'arête la plus critique (prédit σ₁ = 50·6,79 + 2·25·2,6 = 470 MPa ; mesuré 468 ; pour ×0,5 : prédit 404, mesuré 403). Le pic est **1,45–1,6× au-dessus** de cet amorçage : c'est le blocage cinématique des triangles (les joints insérés ne peuvent glisser sans que leurs voisins cèdent). Le maillage frontal (R6 = 0,34) ne cassait pas du tout ; aucun contrôle en h n'a été fait. La sensibilité aux résistances est faible : **exposant mesuré 0,37–0,40** (676 → 514 pour ×0,5), et non « ≈ 0,7 » comme écrit au README ligne 115. | history.csv (`nInserted`, première insertion) ; README:77-80 (algo 6) ; README:115 | 3 |
| D | **Bruit de réalisation jamais mesuré.** Trois maillages `box20x40_h08_s1/s2/s3.msh` existent (`make_box_mesh.py`) et **aucun n'a été lancé** ; 0 réplicat GBM, 0 réplicat Weibull. La bande expérimentale (±2,6 MPa, 0,4 %) servira de poids dans la vraisemblance alors que le bruit de graine du modèle est vraisemblablement 5–20× plus grand → postérieur faussement étroit. | `calib_quick/meshes/` ; liste `out_*` (aucun `_s1/_s2/_s3`) ; `calibration_redbohus/README.md:22-23` (±2,8 / ±2,6 MPa) | 3 |
| E | **Non-identifiabilité structurelle.** 6 paramètres de joint (+ ℓ_I, ℓ_II, μ_res, α, tailles de grains) contre ≈ 4 observables robustes (2 pics, 2 ε_pic). UCS et BTS, cibles du protocole d'août, sont abandonnés. Weibull m = 8 : **709,0 vs 709,4** (aucun effet). Enveloppe : pente exp 13,9 → 6,5 → 4,2 → 3,8 ; joint MC φ = 48° donne Nφ − 1 = 5,8 avant amplification cinématique → la prédiction à 75/100 MPa surestimera (ordre +10–15 %) sauf bascule de mode traction → cisaillement, qui n'a pas eu lieu à 50 MPa (**77 % des joints rompus en traction** : 332/430). | README:71-73 ; `calibration_redbohus/README.md:18-36` ; `fdem_final_joints.csv` (breakMode) ; `CAMPAGNE.md` §3ter réserve 4 et tornado (φ 53° → pente 7,7 vs 6,5) | 3 |

---

## 1. Risques détaillés (gravité 1–3, test de contrôle, coût, remède)

Coûts estimés à partir des murs mesurés (README:59-66, 110-113) : cas 1 = 77 s, ×3 à ×7 quand la casse est abondante ; 14 threads.

### R1 — Déformation plane vs essai axisymétrique (gravité 2)

- **Module apparent.** Vérifié : le fit tangent donne 85,3 GPa = 80/(1−0,25²) (README:87-88, mon refit 85,3). ν du deck = 0,25 alors que ν_exp = 0,29–0,35 (`seuils_sbm_bohus.json`) : (1−ν²) = 0,916 à ν = 0,29 (et non 0,9375).
- **Conversion correcte.** Pendant la phase déviatorique σ_xx est constant, donc Δε_yy = Δσ_yy·(1−ν²)/E exactement (bulk homogène élastique). Il faut **garder E physique (77,7 GPa) et ν = 0,29 dans le deck** et comparer à la cible **ε_2D = ε_exp·(1−ν²)** — et non baisser E à 72 GPa comme suggéré au README:88, ce qui corromprait E pour tout transfert au forage (vitesses d'onde, raideur de contact).
- **Contrainte intermédiaire.** Au pic, σ_zz = ν(σ_xx + σ_yy) ≈ −0,25·(50 + 726) ≈ −194 MPa contre σ₂ = σ₃ = 50 dans l'essai. Les joints ne voient que les tractions dans le plan → le critère de joint est **insensible à σ_zz** : le modèle est un matériau « σ₁–σ₃ » et le pic n'est pas biaisé directement. Mais (a) la roche réelle est sensible à σ₂ (Haimson & Chang 2000, Westerly, DOI 10.1016/S1365-1609(99)00106-9) : le jeu calibré hérite de cette insensibilité, à assumer au transfert 3D ; (b) en GBM les ν par phase (0,17 / 0,25 / 0,36) donnent des σ_zz différents par phase → concentrations parasites aux joints de grains [NON VERIFIE, second ordre] ; (c) `crushCap`/`bulkDamage` intègrent σ_zz dans le déviateur (DOC:197) : un cap réactivé s'atteint plus tôt en déformation plane qu'en axisymétrie.
- **Consolidation non isotrope** : voir R2.
- **Test** : aucun run pour (a) ; pour (c) voir R12. **Remède** : conversion de la cible (0 run), ν = 0,29, documenter la limite σ₂ dans le manuscrit.

### R2 — Offset de consolidation (gravité 3, correction gratuite)

- Pendant `pullDelay`, les nœuds PRESCRIBED ont vg = 0 (`FdemSolver.cpp:7258`) : ε_yy = 0 sous σ_xx = −50 → σ_yy = ν/(1−ν)·σ_xx = −16,67 MPa (mesuré exactement). L'essai, lui, part d'un état isotrope à 50 MPa.
- **Conséquence** : q_README = q_vrai + (σ₃ − 16,67) = q_vrai + 33,3 MPa à 50, + 13,3 à 20. Les 12 lignes de résultats du README (tableaux lignes 69-73, 92-96, 108-113) sont à décaler ; la lecture « les trois représentations encadrent la cible (+18 / −14 %) » devient « +13 / −20 % ».
- **Test** : post-traitement seul — dans `plot_quick.py` remplacer `sig0 = mean(sigma avant pullDelay)` par `sig0 = confiningPressure` et prendre ε = 0 à l'instant où sigma = σ₃ (Δε ≈ 33,3 MPa / 85,3 GPa = 0,039 %). **Remède** durable : une consolidation à σ_yy = σ₃ avant le déviateur (mors pilotés en force pendant `pullDelay`, clé à créer — [NON VERIFIE] qu'elle existe) ; sinon la correction de post-traitement est exacte pour un bulk élastique et suffit.

### R3 — Taille d'éprouvette 20 × 40 mm (gravité 3)

- Dimensions réelles : carottes Φ50 d'après les disques BTS Φ49,4 mm (`extract_targets.py:120`) ; Φ50 × 100 pour le triaxial **[NON VERIFIE — article Dumoulin 2024, DOI 10.1016/j.gete.2024.100592, inaccessible (403)]**. `CAMPAGNE.md` §2 signale un maillage Φ50×100 3D prêt (264 603 tets, ×15).
- **Homogène** : ℓ_cz,I = E·G_f/f_t² = 80e9·40/(12e6)² = **22 mm ≈ W** ; ℓ_cz,II = E·20·G_f/c² = **102 mm ≫ W**. L'éprouvette est dans la zone de transition de la loi d'échelle (Bažant 1984, J. Eng. Mech. 110(4):518, [DOI NON VERIFIE]) : la résistance nominale dépend de W/ℓ_cz. Le balayage `q1s*` le prouve à l'envers (ℓ_cz ×4–16 à W fixe → pic +6–9 %, 3 joints rompus, README:92-100). Passer W de 20 à 50 mm change W/ℓ_I de 0,9 à 2,3 : le pic à paramètres égaux bougera [amplitude NON VERIFIEE, c'est le test].
- **GBM** : 113 grains de 3 mm → W/d = 6,7 < 10, sous le minimum ISRM (Bieniawski & Bernede 1979 : diamètre ≥ 10× le plus gros grain). `CAMPAGNE.md` §3ter réserve 2 jugeait déjà 255 grains insuffisants ; on est à 113. La biotite (7 %) est portée par ~8 grains (1 seul à σ = 0,8, DOC §5.16) : le contraste de phase « −14 % » est un tirage, pas une propriété.
- **Test** : (i) cas 1 sur 40 × 80 mm (`make_box_mesh.py --W 0.04 --H 0.08`, T = 4,1 ms pour la même ε à 0,25 m/s) : ×4 éléments × ×2 pas → ~10 min (peu de casse) à ~1 h ; (ii) 50 × 100 : ×6,25 × ×2,3 → ~18 min à ~2 h ; (iii) GBM 40 × 80 à grains 3 mm (~450 grains) → ~40 min [estimations NON VERIFIEES]. Critère : \|Δq_pic\| < 2σ_seed (R4).
- **Remède** : calibrer à 40 × 80 mm minimum (recommandation §3ter (b)) ; sinon mesurer et déclarer la loi d'échelle du modèle.

### R4 — Bruit de réalisation (gravité 3 tant qu'il n'est pas mesuré)

- Trois sources indépendantes : graine du maillage (`make_box_mesh.py`, champ de taille bruité ±10 %), `fieldSeed` du Weibull (DOC:190), `seed` du Voronoï/phases.
- **Mesure proposée** : rapport signal/bruit S = (∂q/∂ln f_t · Δln f_t) / σ_seed. Sensibilité mesurée : exposant 0,38 → un pas de 10 % sur (f_t, c) déplace le pic de 3,8 % ≈ 23 MPa. Si σ_seed(q_pic) ≳ 20 MPa (3 %), un point de calibration ne distingue plus deux jeux à 10 % près sans moyenner ≥ 3 graines (coût ×3). Mesurer aussi σ_seed(ε_pic) et σ_seed(E).
- **Test** : cas 1 sur `s1/s2/s3.msh` (3 × 77 s = 4 min, decks déjà à un `meshFile` près) ; cas 2 avec `fieldSeed` = 3 valeurs (4 min) ; cas 3 avec `seed` = 5 valeurs (5 × 5 min).
- **Remède** : bande de vraisemblance = max(σ_exp, σ_seed/√n) ; graine FIGÉE par point de plan d'expériences (nombres aléatoires communs) pour les dérivées/tornado ; les 3 graines du protocole d'août (`calibration_redbohus/README.md:65-69`) étaient la bonne pratique, le plan rapide l'a perdue.

### R5 — Vitesse de chargement 0,25 m/s, ε̇ = 6,25 s⁻¹ (gravité 2 pic, 3 post-pic)

- Pré-pic : 1,5 ms de montée pour 7 µs de transit d'onde (40 mm / 5 500 m/s) → quasi-statique en niveau ; la loi de joint est rate-indépendante (aucune clé DIF/viscosité dans les decks). Post-pic : la chute de 55 % se produit en 0,4 ms sous déplacement imposé ; la réserve 1 de `CAMPAGNE.md` §3ter (six ordres au-dessus du labo, 3× Aboayanah 2024) reste ouverte.
- **Test** : cas 1 avec `pullV = -0.05`, `T = 1.0e-2` (pas ×4,5 → ~6 min) ; comparer q_pic, ε_pic, et la chute à +0,1 % après pic. **Remède** : si Δq_pic < σ_seed, garder 0,25 m/s pour le criblage et ÷5 pour les runs de validation ; ne jamais lire le post-pic à 0,25 m/s.

### R6 — `dampingLocal = 0,7` (gravité 3 pour tout observable post-pic, 1 pour le pic)

- Cundall = force α·\|f\|·sign(v) (`FdemSolver.cpp:138`) : nulle en régime stationnaire, maximale dès qu'il y a déséquilibre, c'est-à-dire à la rupture. Budget `_log_q1t050.txt` : Cundall 2 757 J/m, contact 2 723 (dont frottement 912, le reste = restitution normale `gcRestitution = 0,2` défaut), cohésion 1 368, mors 2 746 → **la dissipation cohésive est 20 % du total**. La chute « ~60 % » (README:116) est celle de l'amortissement.
- **Test** : cas 1 avec `dampingLocal` 0,7 / 0,3 / 0,1 + `budgetAbortPct = 5` + `energyBodyForces = on` (3 runs, 80–240 s). Attendu : pic dans le bruit, post-pic différent ; un abort à 0,1 signifierait que 0,7 masque une instabilité.
- **Remède** : critère de validité par run « Cundall < 30 % de la dissipation cohésive au pic » (à sortir dans history, R15) ; le post-pic n'entre dans l'objectif qu'après ce contrôle et avec la vitesse de R5.

### R7 — Mors libres latéralement (gravité 1–2)

- `gripLateralFree` : vx du groupe intégré de la force nette (`:7264-7268`) = platine lubrifiée idéale. Les essais réels ont des platines acier frottantes ; l'effet de frettage à H/D = 2 est de quelques % [NON VERIFIE ; ordre de grandeur littérature Labuz & Bridell 1993].
- **Test** : `loading = platens` + `contactMu = 0.2` vs mors libres (2 runs, ~3 min). **Remède** : garder les mors libres (plus proches d'un essai ISRM lubrifié) et le dire.

### R8 — Rampe de confinement (gravité 1)

- 1e-4 s = 14 transits d'onde, jauge à 3e-4 : −50 MPa atteint à 5e-10 % (log), 0 joint rompu à t = pullDelay (history). Le vrai défaut est R2 (état non isotrope), pas la rampe. À vérifier à chaque run : ligne « confinement: target/achieved » et nBroken(pullDelay) = 0.

### R9 — Insertion adaptative (gravité 2)

- (a) Complaisance nulle : E apparent = E/(1−ν²) exact (mesuré). (b) Les joints **insérés** ont une pénalité `insertionPenaltyFactor = 4·E/h` (log : « activation penalty 4 E/h ») : au pic, 1 462 joints insérés (29 %, `nInserted`) portent une complaisance normale h/(4E) chacun → adoucissement pré-pic partiellement numérique [amplitude NON VERIFIEE]. (c) Le critère d'activation lit la contrainte d'éléments CST à 0,8 mm : la contrainte de pointe est lissée → insertion retardée, dépendante de h (lié à R3/C2). (d) DOC:181 : validé sur un UCS (47,8 MPa, E 99,1 %), jamais sur un triaxial confiné.
- **Test** : `insertion = intrinsic` (dt ÷2 → ~160 s) et `insertionPenaltyFactor = 20` (1 run) sur le cas 1. **Remède** : schéma figé (déjà la règle §5 d'août) mais l'écart intrinsèque/adaptatif à 50 MPa doit être connu et écrit.

### R10 — Pénalité 20 (gravité 1)

Sans effet en adaptatif (aucun joint à t = 0, DOC:171) ; seule la pénalité d'insertion (R9) compte. Ne pas la calibrer.

### R11 — Règle « ℓ_cz constant » (gravité 2 : juste mais incomplète)

- (i) Elle ne couvre que le mode I. ℓ_II = E·`gfShearFactor`·G_f/c² : dès que c et f_t varient indépendamment à `gfShearFactor` fixe, ℓ_II ∝ (f_t/c)² varie (±30 % sur c/f_t → ±70 % sur ℓ_II). Les runs `q1t*` gardaient c/f_t fixe, donc ne l'ont pas vu.
- (ii) ℓ_I n'est pas libre : G_f = K_Ic²(1−ν²)/E, K_Ic granite ≈ 1–2 MPa√m → G_f ≈ 12–48 J/m² [ordre de grandeur, NON VERIFIE pour Bohus] ; G_f = 40 est en haut de fourchette. Fixer ℓ_I par K_Ic, ne pas le calibrer.
- (iii) Régime : W (20 mm) < ℓ_II (102 mm) → régime « contrôlé par la résistance » : G_II est quasi inidentifiable sur le pic ; le fixer (ratio de littérature, Tatone & Grasselli 2015, DOI 10.1016/j.ijrmms.2015.01.011, [valeur NON VERIFIEE]).
- Surprise à intégrer : même à 50 MPa, **77 % des ruptures complètes sont en mode I** (332/430 ; 706/925 pour ×0,5 ; 220/315 en GBM). L'amorçage est Coulomb (c, φ), la propagation est traction (f_t, ℓ_I) : les deux longueurs comptent.
- **Test** : cas 1 avec c × 0,7 et `gfShearFactor` (a) fixe (b) recalculé (c/f_t)²·ℓ_II/ℓ_I (2 runs). **Remède** : paramétrer (f_t, c, φ, ℓ_I, ℓ_II) avec ℓ_I, ℓ_II figés ; dériver G_f = f_t²ℓ_I/E et `gfShearFactor` = (c/f_t)²·ℓ_II/ℓ_I dans le générateur de decks.

### R12 — `crushCap = 1e12`, bulk purement élastique (gravité 2)

- Toute la non-linéarité pré-pic doit venir des joints. Essai : sécante/tangente au pic = 62/77 = 0,81 (19 % de déformation non linéaire, mes fits sur 2_5–2_7) ; modèle : 81/85 = 0,95 (5 %). D'où le **déficit structurel de ε_pic ≈ 15 %** (0,83 % vrai vs 0,95), déjà constaté (« ~20 % de TOUTES les configs homogènes », `CAMPAGNE.md` §3ter). Un objectif « courbe entière » fera baisser E ou les résistances pour rattraper ε_pic → biais sur (c, φ).
- Biotite à E = 29 GPa sans plafond : von Mises > 500 MPa sans plasticité — non physique [seuil réel NON VERIFIE]. Le défaut 8·c = 200 MPa aurait tout plafonné ; la campagne 3D avait 400 MPa « découplé » (`CAMPAGNE.md` §2).
- **Test** : cas 1 avec `crushCap` 700 / 500 MPa (2 × 80 s) ; si le pic bouge, c'est un paramètre caché (Ye et al. 2025 calibrent un bulk MC, `calibration_redbohus/README.md:40-43`). **Remède** : soit cap figé et documenté, soit paramètre calibré ; dans les deux cas déclarer que la courbure pré-pic n'est pas une cible pour le modèle homogène (ou lui donner la cible `nDamaging`, R14).

### R13 — « Chute 26 % » mesurée en fin d'essai (gravité 3 comme cible)

- Réplicats à 50 : 0,427 / 0,330 / 0,021 ; à 20 : chute à +1 000 µε après pic = 0,17 / 0,26 / 0,40 (mes calculs). Dispersion ±0,1–0,2 sur une cible de 0,26.
- Côté modèle : μ_res = 0,25 → φ_r = 14°, N_φr − 1 = 0,64 → résiduel d'une bande traversante = 50·0,64 ≈ **32 MPa (5 % du pic)**. Le « résiduel » expérimental (443 MPa, 74 %) n'en est pas un : c'est un adoucissement interrompu. La cible compare un essai arrêté à une simulation pilotée par l'amortissement (R6). μ_res est donc **inidentifiable** ici.
- **Remède** : (i) cible = chute à Δε fixe après pic là où n_rep = 3 (à 20 MPa : +1 000 µε, bande ±0,1) ; (ii) traiter chaque fin d'essai comme donnée censurée (chute vraie ≥ mesurée) ; (iii) fixer μ_res à une valeur de littérature (0,18 Kuru, DOC:177) et le sortir de la calibration ; (iv) `stopPeakDrop` n'existe qu'en fdem3d (DOC:722, aucune occurrence dans `FdemSolver.cpp`) → le porter en 2D pour ne pas payer 0,4 ms de débris par run (le mur passe de 77 à 530 s avec la casse, README:123).

### R14 — σ_ci / σ_cd : mesurables ou non (gravité 2)

- Essai : CI = 55–57 %, CD = 62–73 % du pic (`seuils_sbm_bohus.json`), obtenus par ε_vol (méthode SBM).
- Modèle : `history.csv` de ce binaire porte `nInserted, nDamaging` (le §6.1 de la DOC est en retard) → proxy CI = première insertion = **62–69 %** (homogène), **30 %** (GBM) sur base déviateur vrai. Mais en homogène cette première insertion est **exactement l'enveloppe Coulomb** de l'arête la plus critique (470 → 468) : elle mesure (c, φ), pas une physique d'amorçage. CD (croissance instable) ≈ premier joint entièrement rompu : nBroken(pic) = 6 → CD ≈ 100 % du pic contre 72 % en essai : **le modèle homogène n'a pas de phase CD → pic**. Aucune ε_lat dans history → pas d'ε_vol → pas de CD au sens SBM.
- **Remède (code, 0 run)** : colonne `epsLat` (déplacement moyen des faces latérales / W) → ε_vol = ε_ax + ε_lat en déformation plane ; passer les courbes modèle dans `sbm_seuils.py` pour des seuils comparables. C'est l'observable qui rend Weibull et GBM identifiables (PLAN_GBM §3.3).

### R15 — Ce qui manque dans `history.csv` (gravité 2)

Présent (grips) : `t, gripFy, sigma, sigmaPeak, nBroken, nInserted, nDamaging`. Manquent : (1) ε_lat/ε_vol ; (2) réaction du mors bas (contrôle d'équilibre \|F_top\| − \|F_bot\|, existant en platens) ; (3) postes d'énergie par ligne (`eEl, eJnt, eGc, eFric, eCund` existent en percussion, `FdemSolver.cpp:7622` vs `:7639`) → indispensable pour R6 ; (4) `nBrokTen/nBrokShear` par ligne (platens seulement, `:7635`) ; (5) déplacement réel des mors (plot_quick suppose la rampe analytique) ; (6) D max / D moyen ; (7) `confAchieved` par ligne. `sigmaPeak` est un max glissant (DOC §6.1) : ne jamais le lire. `fdem_final_joints.csv` (tBreak, breakMode) permet le catalogue AE des ruptures complètes ; les temps d'insertion ne s'y trouvent pas [`tInsert` de rockim_f1, présence en f2j NON VERIFIEE].

### R16 — Étape Weibull non identifiable sur q(ε) (gravité 2)

m = 8 indépendant : 709,0 vs 709,4, 431 vs 430 rompus, première insertion 409 vs 468 MPa (seul effet visible : l'amorçage). Un système en compression est « parallèle » (redistribution), Weibull agit sur les systèmes « série » (BTS, traction) et sur CI. `CAMPAGNE.md` §3ter (« théorème du plancher ») l'avait établi. **Remède** : calibrer (m, `strengthCorrLength` ≥ taille de grain, `fieldSeed` × 3) sur la dispersion du BTS et sur le proxy CI (R14), pas sur le pic ; sans ces observables, l'étape 2 est 12 runs pour rien.

### R17 — Spécificités GBM (gravité 2)

- Grains 3 mm vs Bohus 1–3 mm (`calibration_redbohus/README.md:11-12`) ; W/d = 6,7 (R3).
- Cohésions par phase inversées par la contrainte « moyennes égales » : biotite c = 30 MPa > feldspath 19,44 (deck q3) — le minéral le plus faible est le second plus cohésif.
- Joints de grains non affaiblis (α = 1) : fraction intergranulaire des ruptures = 32/315 = 10 % ≈ leur proportion (985/8 609 = 11 %) → aucune préférence de trajet ; or le granite rompt aux joints de grains (Aboayanah et al. 2024, DOI 10.1007/s00603-024-03789-7). α_GB sera le levier fort, mais à 113 grains le réseau intergranulaire n'a pas de statistique de percolation.
- Polydispersité : dt ÷3,6 à σ = 0,8 (DOC §5.16) et biotite à 1–8 grains → bruit pur.
- Explosion paramétrique : 3 phases × 6 propriétés + α + tailles → inidentifiable sur 2 courbes. **Remède** : figer les propriétés de phases sur les tables de littérature (PLAN_GBM §1, Mahabadi/Villeneuve) et ne calibrer que α_GB (1–2 nombres) + éventuellement `gbHeteroFactor` ; 5 graines avant tout balayage.

### R18 — Enveloppe MC linéaire vs concave (gravité 2, prédiction 75/100)

Pentes locales exp : 13,9 / 6,5 / 4,2 / 3,8. Extrapolation linéaire de (20, 50) à 100 : 924 MPa vs 798 (+16 %). Le modèle ne peut courber l'enveloppe que par bascule de mode (traction → cisaillement) ; à 50 MPa il est encore à 77 % traction. **Test le moins cher du plan** : sonde 4 à `confiningPressure` = 20e6 et 100e6 (2 runs, 80–300 s) → pente propre du modèle avant de calibrer. **Remède** : si la pente modèle > 5 entre 50 et 100, annoncer la surestimation comme résultat attendu, ou envisager une enveloppe de joint non linéaire [inexistante dans rockim, NON VERIFIE].

### R19 — Coût et arrêt (gravité 1)

Le budget « < 5 min » n'est tenu que sans casse (77 → 530 s). Porter `stopPeakDrop` en 2D ou réduire T une fois ε_pic borné (le pic tombe avant 1,8 ms sur tous les runs).

### R20 — Reproductibilité (gravité 1)

Comparaisons fines à nombre de threads égal (DOC §8.3) : tous les runs de calibration à 14 threads, sans exception, et binaire figé (`rockim_f2l.exe` pour tout, puisque bit-identique clés absentes).

---

## 2. Stratégie d'identification recommandée (à la place de « tout sur q(ε) »)

1. **Figer** : ℓ_I (K_Ic), ℓ_II (ratio littérature), μ_res (0,18–0,25), pénalités, schéma adaptatif, h, taille d'éprouvette, ν = 0,29, E = 77,7 GPa.
2. **Observables** : q_pic(20), q_pic(50) corrigés (R2) ; ε_pic convertis (R1) ; proxy CI (`nInserted`, R14) ; BTS 10,3 MPa (réintroduit : seul observable qui fixe f_t) ; chute à Δε fixe censurée (R13) ; bande = max(σ_exp, σ_seed/√n).
3. **Paramètres calibrés** : φ (pente d'enveloppe — levier dominant du tornado d'août), c (niveau), f_t (BTS) — 3 paramètres pour 5–6 observables ; puis α_GB seul en GBM.
4. Résultat attendu à énoncer d'avance : surestimation à 100 MPa (R18) et déficit ε_pic en homogène (R12) — ce sont des verdicts sur la représentation, pas des échecs de calibration.

---

## 3. LISTE ORDONNÉE des contrôles à faire AVANT toute heure de calibration

Chaque ligne = une variable, deck dérivé du cas 1 (`q1_homog_P050.cfg`) sauf mention ; coûts à 14 threads ; à valider par Fernando avant lancement.

| ordre | contrôle | deck (diff au cas 1) | coût | décision attendue |
|---|---|---|---|---|
| **C0** | Corriger le dépouillement : q = sigma − σ₃, ε = 0 à sigma = σ₃, cible ε×(1−ν²) ; retabuler le README | `plot_quick.py` | 0 run, 1 h | nouveau point de départ (facteur ≈ 0,75) |
| **C1** | Bruit de graine du maillage | `meshFile = …_s1/_s2/_s3.msh` | 3 × 77 s | σ_seed(q_pic, ε_pic, E) → bande de vraisemblance, nombre de graines par point |
| **C2** | Objectivité en h | `make_box_mesh.py --h 0.5e-3` (et 1,2e-3) | 2 runs, ~5 min + 1 min | si Δq > 2σ_seed : calibrer au h du maillage de forage, pas à 0,8 mm |
| **C3** | Pente propre du modèle | `confiningPressure = 20e6` / `100e6` | 2 runs, 2–5 min | faisabilité de la prédiction 75/100 ; φ de départ |
| **C4** | Amortissement | `dampingLocal = 0.3` / `0.1` + `budgetAbortPct = 5` + `energyBodyForces = on` | 3 runs, ~8 min | domaine de validité du post-pic ; critère Cundall/cohésif |
| **C5** | Taille d'éprouvette | 40 × 80 mm, `T = 4.1e-3` (puis 50 × 100 si C5a bouge) | 10–60 min (+20–120) | taille de calibration ; loi d'échelle du modèle |
| **C6** | Vitesse | `pullV = -0.05`, `T = 1.0e-2` | ~6 min | rate-indépendance du pic ; vitesse des runs de validation |
| **C7** | Cap du bulk | `crushCap = 700e6` / `500e6` | 2 × 80 s | paramètre caché ou non |
| **C8** | Schéma d'insertion | `insertion = intrinsic` ; `insertionPenaltyFactor = 20` | 2 runs, ~5 min | écart à déclarer, schéma figé |
| **C9** | Mode II | c × 0,7 avec `gfShearFactor` fixe vs recalculé | 2 runs, ~5 min | règle (ℓ_I, ℓ_II) dans le générateur |
| **C10** | Code : `epsLat`, réaction bas, postes d'énergie, `nBrokTen/Shear` dans history grips ; `stopPeakDrop` 2D | `FdemSolver.cpp` (ajouts opt-in) | 0 run, ½ journée | rend CI/CD, R6 et le coût contrôlables |
| **C11** | Bruit de graine GBM | `q3` avec `seed` × 5 | 5 × 5 min | σ_seed GBM ; taille de grain/éprouvette pour la phase 3 |
| **C12** | Mors | `loading = platens`, `contactMu = 0.2` | 1 run, 2 min | frettage à déclarer |

Total : ~25 runs, 1,5–2 h de mur hors C5 grand format ; C0–C4 (≈ 25 min) suffisent pour décider si la calibration en 20 × 40 mm a un sens.

---

## 4. Sources

- Fichiers : `FDEM/rockim_f2/calib_quick/README.md` (l. 16-22, 59-73, 87-88, 92-124) ; `q1_homog_P050.cfg`, `q3_gbm_P050.cfg`, `q1t050_P050.cfg` ; `out_q1_homog_P050/history.csv` ; `calib_quick/_log_q1t050.txt` ; `out_*/fdem_final_joints.csv` ; `src/FdemSolver.cpp` (l. 112-140, 4150-4156, 6264-6270, 7256-7268, 7622-7646, 826) ; `DOCUMENTATION_rockim.md` (l. 126, 171-197, 677-726, 1246-1322, 1329-1345, 1448-1467) ; `rockim_f1/calib_triax3d/CAMPAGNE.md` (§0-2, §3ter l. 375-445) ; `rockim_f1/calib_triax3d/PLAN_GBM.md` ; `rockim_p1/calibration_redbohus/README.md` (l. 8-36, 61-75) ; `rockim_p1/calibration_redbohus/tools/extract_targets.py:120` ; `CONTINUUM/calib_bohus_triax/exp_qc/experimental_data_red_bohus_clean.json`, `seuils_sbm_bohus.json` ; `rockim_f1/calib_triax3d/targets_triax_bohus.json` ; `calib_quick/make_box_mesh.py`, `calib_quick/meshes/`.
- Littérature : Dumoulin et al. 2024, Geomech. Energy Environ. 40:100592, DOI 10.1016/j.gete.2024.100592 ; jeu de données Zenodo 10.5281/zenodo.10617548 ([page consultée](https://zenodo.org/records/10617548), dimensions d'éprouvettes absentes) ; Aboayanah et al. 2024, RMRE 57:4679–4706, DOI 10.1007/s00603-024-03789-7 ([lien](https://link.springer.com/article/10.1007/s00603-024-03789-7)) ; Haimson & Chang 2000, IJRMMS 37:285–296, DOI 10.1016/S1365-1609(99)00106-9 ([lien](https://www.sciencedirect.com/science/article/abs/pii/S1365160999001069)) ; Tatone & Grasselli 2015, IJRMMS 75:56–72, DOI 10.1016/j.ijrmms.2015.01.011 ([lien](https://www.sciencedirect.com/science/article/abs/pii/S1365160915000180)) ; Bieniawski & Bernede 1979, ISRM SM UCS (ratio diamètre/grain ≥ 10, [revue citée](https://www.sciencedirect.com/science/article/pii/S1674775524001914)) ; Bažant 1984 J. Eng. Mech. 110(4):518–535 [DOI NON VERIFIE] ; Ye 2025 / Bu 2026 / Jiang 2025 tels que cités dans `calibration_redbohus/README.md:38-51`.
- [NON VERIFIE] récapitulatif : dimensions Φ50×100 des triaxiaux Bohus ; amplitude de l'effet σ₂ pour Bohus ; K_Ic de Bohus ; valeur du ratio G_II/G_I chez Tatone & Grasselli ; seuil de plasticité de la biotite ; existence d'une clé de mors en force pendant `pullDelay` ; `tInsert` en f2j ; effet quantitatif de `insertionPenaltyFactor` sur ε_pic ; estimations de mur des runs 40×80 / 50×100 / GBM.