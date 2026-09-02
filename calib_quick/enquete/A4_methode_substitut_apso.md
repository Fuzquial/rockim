# TÂCHE A4 — Calibration par substitut + PSO adaptatif : littérature, dimensionnement et pipeline

*Rédigé le 2026-09-02 pour les trois calibrations Red Bohus (homogène → Weibull → GBM), σ₃ = 20/50 MPa au calage, 75/100 en prédiction. Toutes les valeurs chiffrées sont sourcées (DOI ou `fichier:ligne`) ; ce qui n'a pas pu être confirmé est marqué [NON VERIFIE].*

---

## 0. Recommandation en dix lignes

1. **Ne pas faire tourner un PSO (adaptatif ou non) directement sur rockim.** Un APSO à 20 particules × 30 générations = 600 jeux × 2 confinements × 2,4–8,8 min = 50–180 h, et chaque évaluation est **bruitée par la graine** (maillage, champ Weibull, Voronoï) : pBest/gBest chassent le bruit. C'est exactement le régime pour lequel la littérature « surrogate-assisted » existe (Jin 2011 ; Regis 2014 ; Sun et al. 2017).
2. **Architecture retenue : LHS → processus gaussien (Matérn 5/2 ARD + nugget) par observable → enrichissement séquentiel par lots → APSO comme optimiseur *sur le substitut* → postérieur MCMC sur le substitut → validation par runs directs.** C'est la continuité du pipeline d'août (`rockim_p1/calibration_redbohus/tools/analyze.py`, `enrich.py`) avec trois corrections : n₀ ≥ 10·d (août : 44 points en 6D → R² croisé 0,44 sur l'UCS, `tools/enrich.py:9-10`), traitement explicite du bruit de réalisation (réplicats + nugget), et objectif **sur la courbe** (6 descripteurs + RMSE pondérée).
3. **APSO (Zhan et al. 2009)** : à implémenter (≈ 70 lignes numpy, code §8.6) et à utiliser comme optimiseur interne du substitut et de l'acquisition. Son adaptativité (w, c₁, c₂, saut élitiste) est un vrai plus pour la robustesse aux optima locaux du substitut, mais **ce n'est pas d'elle que vient l'économie de calcul** — elle vient du substitut et du plan séquentiel.
4. **Forêt aléatoire : non comme substitut principal** (pas d'extrapolation, incertitude non calibrée, Bu et al. 2026 : meilleur R² direct mais mauvaise en inversion — `calibration_redbohus/METHODOLOGIE.md:123-125`). À garder comme **contrôle croisé** (importance par permutation, LOO) — scikit-learn 1.8.0 est présent.
5. **Tailles** (d = 6 homogène) : **n₀ = 60** LHS (`scipy.stats.qmc.LatinHypercube(optimization="random-cd")`) + **17 réplicats** pour le nugget + **6 lots × 8** d'enrichissement + **22 runs** de validation ≈ **270 runs ≈ 11–40 h** (146–530 s/run mesurés). GBM (d = 7) : ≈ 300 runs ≈ 25–45 h à grains 3 mm.
6. **Fonction objectif** : J = Σ_{σ₃} ½ Σ_k w_k ((y_k − y_k*)/τ_k)² avec 6 descripteurs (pic 0,30 ; RMSE-bande 0,20 ; ε_pic 0,15 ; σ_ci/pic 0,15 ; E 0,10 ; chute 0,10) et tolérances physiques τ_k (§6.2) → J ≈ 1 signifie « tout à la tolérance ».
7. **Identifiabilité** : deux pics confinés ne donnent que **deux nombres** (ordonnée et pente d'enveloppe MC) pour trois résistances (ft, c, φ) → vallée (ft, c). On la lève **sans nouvel essai coûteux** : σ_ci (plancher de frottement, `rockim_f1/calib_triax3d/CAMPAGNE.md:222-231`) fixe (c, φ) au second ordre ; **ft est fixé par le BTS (10,27 ± 0,98 MPa, `calibration_redbohus/README.md:21`) via un run de traction/brésilien de 30–45 s**, ce qui ramène d à 5 pour la partie chère.
8. **Reparamétrer** : (ln c, tan φ, ln(ft/c), ln ℓ_cz, ln G_II/G_I, μ_res) — ℓ_cz = E·G_f/ft² constant est la règle établie le 2026-09-02 (`calib_quick/README.md:98-104`).
9. **Validation 75/100** : attendre un biais **+8 % (75) et +16 % (100)** d'une enveloppe linéaire calée sur 20/50 (calcul §7.1 sur `targets_triax_bohus.json`) ; l'intervalle prédictif = tirages du postérieur × graines (22 runs). Un biais dans cette fourchette est une **erreur de forme de modèle**, pas de calibration — à documenter, non à corriger par les paramètres (`CAMPAGNE.md:412-415`).
10. **Avant tout lancement** : règle maison — liste exacte des clés + coût → validation Fernando (mémoire `feedback-config-validation-avant-run`). Le pipeline est un script Python lancé par Fernando dans son terminal (comme `campaign.py` en août).

---

## 1. Ce que l'existant impose (état des lieux chiffré)

| fait | valeur | source |
|---|---|---|
| Cibles courbes 20/50/75/100 : pic 404,0 / 599,2 / 703,0 / 798,2 MPa ; ε_pic 0,65 / 0,95 / 1,08 / 1,20 % ; chute 30 / 26 / 25 / 7 % ; grille commune 400 pts, bande ±1σ plancher 5 MPa | — | `rockim_f1/calib_triax3d/targets_triax_bohus.json` (`_meta.regle`), `CAMPAGNE.md:37-42` |
| Seuils SBM expérimentaux : CI ≈ 0,55–0,57 du pic (σ₃ = 20 et 50), CD ≈ 0,62–0,73 ; E 75–78 GPa ; ν 0,29–0,35 | — | `CONTINUUM/calib_bohus_triax/exp_qc/seuils_sbm_bohus.json` (essais 2_2…2_7) |
| UCS 126,6 ± 21,4 MPa (4 essais) ; BTS 10,27 ± 0,98 MPa ; q(20) 404,8 ± 2,8 ; q(50) 599,0 ± 2,6 | — | `rockim_p1/calibration_redbohus/README.md:18-25` ; ⚠️ `CAMPAGNE.md:413` parle d'un « UCS Bohus réel ~175 » — **contradiction non tranchée [NON VERIFIE]** |
| Coût réel d'un run cas 1 sur le maillage **isotrope** (algo 5, 3388 triangles, 541 571 pas) : **146 s** ; cas 2 : 148 s ; cas 3 GBM 3 mm : **295 s** ; cas 1 avec 1738 joints rompus : **530 s** | 14 threads | `rockim_f2/out_q1_homog_P050/run.log`, `out_q2_weibull_P050/run.log`, `out_q3_gbm_P050/run.log`, `calib_quick/_log_q1t025.txt`. ⚠️ Le « 77 s » de `calib_quick/README.md:63` est celui du maillage algo 6 **banni** (2930 triangles) |
| GBM à grains 2 mm (19 476 facettes, rockim_f1) : 400–1 800 s/run | 14 threads | `rockim_f1/calib_triax3d/gbm1_P050.log`, `gbm2_P020.log`, `out_scr_base_P020.log` |
| Weibull m = 8 sur les joints : **aucun effet** sur le pic (709,0 vs 709,4) ; m ∈ [3 ; 24] n'avance pas σ_ci (plancher de frottement) | — | `calib_quick/README.md:81-83` ; `CAMPAGNE.md:222-231` |
| Balayage (ft, c) à ℓ_cz constant : pic ∝ facteur^0,7 ; cible 599 → facteur ≈ 0,7 (ft ≈ 8,4, c ≈ 17,5, G_f ≈ 20) ; mais ε_pic 0,72–0,87 % vs 0,95 et chute 52–64 % vs 26 % | — | `calib_quick/README.md:106-124` |
| Tornado bi-confinement 2D (22 runs, régime pré-endommagé) : frictionDeg ≫ jointResidualMu > E > gfShearFactor > Gf, frac, ft ; +5° → +107 MPa ; ft/c +25 % → +8/+10 MPa | — | `CAMPAGNE.md:331-349, 426-429` |
| Sensibilité d(ln pic)/d(ln ft,c) ≈ 0,19 à frac 0,15 ; μ_res : +315 MPa par unité | — | `CAMPAGNE.md:302-330` |
| E_apparent = E/(1−ν²) en déformation plane (85,3 pour 80) → E_deck ≈ 72 GPa pour 77 apparent : **E se fixe analytiquement, pas par calibration** | — | `calib_quick/README.md:87-88` |
| Bruit de graine mesuré (BTS, 3 graines, GBM 5 mm) : 9,90 / 10,02 / 9,71 MPa → CoV ≈ 1,6 % ; (A3) 9,95 / 10,08 / 9,70 | — | `calibration_redbohus/phaseA_results.csv` |
| Extrapolation d'août (CALT calé sur UCS/BTS/tx20) : tx50 814 (+36 %), tx75 1 002 (+42 %), tx100 1 222 (+53 %) — l'enveloppe MC linéaire des joints **sur-prédit** les hauts confinements | — | `calibration_redbohus/points_results.csv` |
| Criblage d'août (GBM 5 mm) : BTS 28,3 → 32,3 pour ft 20 → 60 MPa (faible) mais 14,2 → 37,0 pour c 10 → 90 : le BTS de ce GBM était piloté par **c**, pas par ft | — | `calibration_redbohus/screen_results.csv` lignes `ft_lo/hi`, `cohesion_lo/hi` |
| PSO existant (Abaqus CDP) : « adaptatif » = rampes linéaires w 0,9 → 0,3, c₁ 2,5 → 0,5, c₂ 0,5 → 2,5 + bonus sur taux de succès ; v2 12 × 25, v3 8 × 20 | — | `python/PSO/pso_calibration.py:89-99, 144-156` ; `PSO_CALIBRATION_DOC.md:78-81, 108-112` |
| Outils Python présents : numpy 2.4.4, scipy 1.17.1 (`qmc.LatinHypercube`, `stats.sobol_indices`), scikit-learn 1.8.0 (GPR avec `alpha` par point, RF/ExtraTrees), pandas, matplotlib, joblib. **Absents** : pyDOE, SALib, pyswarms, scikit-optimize, GPy, pymoo, emcee | — | vérifié `python -c import…` le 2026-09-02 |
| Machine : 18 processeurs logiques ; août : 3 runs × 6 threads (`tools/campaign.py:23-24`) ; ici 1 × 14. Rendement 2 × 7 vs 1 × 14 sur `rockim_f2j.exe` : [NON VERIFIE] — à mesurer sur un lot de 2 runs avant le LHS | — | — |

**Lacune à combler avant de calibrer sur σ_ci** : `history.csv` (scénario tension/grips) ne contient que `t, gripFy, sigma, sigmaPeak, nBroken`, et `nBroken` ne compte que D ≥ 1 ; « ni le nombre de joints insérés ni l'endommagement maximal ne sont écrits — lacune connue » (`DOCUMENTATION_rockim.md` §6.1, lignes 1329-1343). Deux voies : (a) définir σ_ci **sur la courbe** (écart de 2 % à la tangente élastique — même opérateur appliqué à l'exp et à la sim, recommandé) ; (b) ajouter opt-in une colonne `nInserted` au scénario tension (existe déjà en shpb, ligne 1337) — développement mineur par addition.

---

## 2. APSO — Zhan, Zhang, Li, Chung (2009), IEEE TSMC-B 39(6) 1362-1381

DOI [10.1109/TSMCB.2009.2015956](https://doi.org/10.1109/TSMCB.2009.2015956) ; texte intégral lu sur l'eprint Glasgow ([eprints.gla.ac.uk/7645](https://eprints.gla.ac.uk/7645/1/7645.pdf)), notice [PubMed 19362911](https://pubmed.ncbi.nlm.nih.gov/19362911/). Numéros d'équations = ceux de l'article.

**2.1 PSO de base (éq. 1-2)** : v_i^d ← ω v_i^d + c₁ r₁ (pBest_i^d − x_i^d) + c₂ r₂ (gBest^d − x_i^d) ; x_i^d ← x_i^d + v_i^d ; |v^d| ≤ V_max^d = 20 % de l'étendue (p. 1363). Historique : Kennedy & Eberhart 1995 (10.1109/ICNN.1995.488968), inertie Shi & Eberhart 1998 (10.1109/ICEC.1998.699146, éq. 3 : ω = ω_max − (ω_max − ω_min) g/G, 0,9 → 0,4), constriction Clerc & Kennedy 2002 (10.1109/4235.985692, χ = 0,729, φ = 4,1).

**2.2 Estimation de l'état évolutif (ESE)** :
- Étape 1 (éq. 7) : distance moyenne de chaque particule aux autres, d_i = (1/(N−1)) Σ_{j≠i} ‖x_i − x_j‖ (euclidienne sur les D dimensions — **normaliser les dimensions à [0,1] avant**, sinon les paramètres à grande étendue dominent).
- Étape 2 (éq. 8) : facteur évolutif **f = (d_g − d_min)/(d_max − d_min) ∈ [0,1]**, d_g = distance de la particule gBest.
- Étape 3 (éq. 9a-d) : classification floue en quatre états, fonctions d'appartenance :

| état | μ(f) |
|---|---|
| S₁ exploration | 0 (f ≤ 0,4) ; 5f − 2 (0,4 < f ≤ 0,6) ; 1 (0,6 < f ≤ 0,7) ; −10f + 8 (0,7 < f ≤ 0,8) ; 0 (f > 0,8) |
| S₂ exploitation | 0 (f ≤ 0,2) ; 10f − 2 (0,2–0,3] ; 1 (0,3–0,4] ; −5f + 3 (0,4–0,6] ; 0 (f > 0,6) |
| S₃ convergence | 1 (f ≤ 0,1) ; −5f + 1,5 (0,1–0,3] ; 0 (f > 0,3) |
| S₄ saut (jumping-out) | 0 (f ≤ 0,7) ; 5f − 3,5 (0,7–0,9] ; 1 (f > 0,9) |

  Défuzzification « singleton » + base de règles suivant la séquence S₁ ⇒ S₂ ⇒ S₃ ⇒ S₄ ⇒ S₁ : en zone de recouvrement, on garde l'état précédent s'il est admissible (stabilité), sinon le successeur dans la séquence (p. 1367).

**2.3 Inertie adaptative (éq. 10)** : **ω(f) = 1/(1 + 1,5 e^{−2,6 f}) ∈ [0,4 ; 0,9]**, ω₀ = 0,9. ω suit f, pas le temps : grande en exploration/saut, petite en exploitation/convergence.

**2.4 Coefficients d'accélération (Table II, éq. 11-12)** : c₁ = c₂ = 2,0 au départ ; à chaque génération, selon l'état :

| état | c₁ | c₂ |
|---|---|---|
| S₁ exploration | + δ | − δ |
| S₂ exploitation | + 0,5 δ | − 0,5 δ |
| S₃ convergence | + 0,5 δ | + 0,5 δ |
| S₄ saut | − δ | + δ |

avec |c_i(g+1) − c_i(g)| ≤ δ, **δ ~ U[0,05 ; 0,1]** (éq. 11) ; c₁, c₂ bornés à **[1,5 ; 2,5]** ; somme bornée à [3,0 ; 4,0] et, si c₁ + c₂ > 4,0, renormalisation **c_i ← 4,0·c_i/(c₁ + c₂)** (éq. 12).

**2.5 Apprentissage élitiste (ELS, éq. 13-14, fig. 7)** — uniquement en état S₃ : choisir **une** dimension d au hasard de gBest, P^d ← P^d + (X_max^d − X_min^d)·N(0, σ²), **σ = σ_max − (σ_max − σ_min) g/G, σ_max = 1,0, σ_min = 0,1** ; borner au domaine ; évaluer ; si meilleur que gBest il le remplace, **sinon il remplace la pire particule** (diversité).

**2.6 Résultats revendiqués (Table III)** : sur 12 fonctions test en D = 30, N = 20 particules, l'adaptation seule réduit le nombre d'évaluations pour atteindre l'acceptation de ×1,5 à ×26 (GPSO sur la sphère : 106 367 → 6 578 FE).

**2.7 Ce que ça vaut pour nous** : sur le simulateur, chaque FE = 5–18 min et bruitée → inadapté. Sur le substitut (FE ≈ 1 ms), l'APSO est un excellent **optimiseur global d'acquisition** (Regis 2014 : PSO + substitut RBF, 10.1016/j.jocs.2013.07.004 ; Sun et al. 2017, SA-COSO, 10.1109/TEVC.2017.2675628 ; revue Jin 2011, 10.1016/j.swevo.2011.05.001). Réglage recommandé sur le substitut : N = 40, G = 300, V_max = 0,2, 20 redémarrages, coût total < 10 s.

---

## 3. Calibration assistée par substitut en DEM/FDEM — ce que la littérature a établi

| référence (DOI vérifié Crossref) | méthode | leçon pour nous |
|---|---|---|
| Yoon 2007, IJRMMS 44(6) 871-889, 10.1016/j.ijrmms.2007.01.004 | Plackett–Burman (criblage) + plan composite centré + surface de réponse quadratique + optimisation, PFC UCS | première formalisation DOE → RSM pour des micro-paramètres de roche ; la RSM quadratique est trop rigide pour nos réponses à seuils (runs « locked ») |
| Hanley et al. 2011, Powder Tech. 210 230-240, 10.1016/j.powtec.2011.03.023 | plans de Taguchi, DEM agglomérats liés | interactions entre paramètres non négligeables → un OAT ne suffit pas (confirmé par notre couplage arc↔pic) |
| Tatone & Grasselli 2015, IJRMMS 75 56-72, 10.1016/j.ijrmms.2015.01.011 | **procédure de calibration FDEM 2D** (Y-Geo) sur UCS/BTS : les paramètres cohésifs sont **attachés au maillage** | h figé pour toute la campagne (`CAMPAGNE.md:13-15`) ; les jeux calibrés portent l'étiquette (h, ε̇, 2D) |
| Benvenuti, Kloss, Pirker 2016, Powder Tech. 291 456-465, 10.1016/j.powtec.2016.01.003 | réseau de neurones comme substitut direct DEM | régime « grandes bases » ; pas le nôtre (≤ 200 points) |
| Rackl & Hanley 2017, Powder Tech. 307 73-83, 10.1016/j.powtec.2016.11.048 | **LHS + krigeage + optimisation** sur le métamodèle, DEM vrac | le schéma de base recommandé ici ; ils insistent sur la vérification finale par runs directs |
| Coetzee 2017, Powder Tech. 310 104-142, 10.1016/j.powtec.2017.01.015 | revue de la calibration DEM | non-unicité des jeux micro → exiger plusieurs essais/observables |
| Cheng, Shuku, Thoeni, Yamamoto 2018, Granular Matter 20:11, 10.1007/s10035-017-0781-y | **filtre quasi-Monte-Carlo séquentiel bayésien**, DEM triaxial | postérieur complet des micro-paramètres ; base de GrainLearning |
| Cheng et al. 2019, CMAME 350 268-294, 10.1016/j.cma.2019.01.027 ; GrainLearning : Cheng et al. 2024, JOSS 9(97) 6338, 10.21105/joss.06338 | **filtrage bayésien itératif** (SMC + mélange gaussien pour ré-échantillonner) | l'itération « échantillonner là où le postérieur a du poids » = notre enrichissement adaptatif |
| Do, Aragón, Schott 2018, Adv. Powder Tech. 29(6) 1393-1403, 10.1016/j.apt.2018.03.001 ; Mohajeri et al. 2020, 31(5) 1838-1850, 10.1016/j.apt.2020.02.019 | algorithme génétique **directement** sur le DEM | faisable seulement à ≤ 1 min/run ; sinon substitut |
| Richter et al. 2020, Powder Tech. 360 967-976, 10.1016/j.powtec.2019.10.052 | krigeage + optimisation pour la calibration standardisée | même schéma que Rackl–Hanley, avec protocole d'essais imposé |
| Qu et al. 2020, Powder Tech. 366 527-536, 10.1016/j.powtec.2020.02.077 | Adam « physics-informed » sur BPM | gradient approché par relations micro-macro : intéressant, mais nos réponses sont non lisses (rupture/locked) |
| Westbrink et al. 2021, Powder Tech. 379 602-616, 10.1016/j.powtec.2020.10.067 | apprentissage par renforcement multi-objectif | exotique ; pas de gain en petit budget |
| **Wang, Lu, Wan, Zhao 2021, Adv. Powder Tech. 32(2) 358-369, 10.1016/j.apt.2020.12.015** | **PSO amélioré** pour les micro-paramètres DEM roche ; conclusion (résumé vérifié) : « different sets of microparameters can be determined when few macroparameters are used » | la référence PSO-roche la plus proche de la demande de Fernando — et elle **démontre la non-unicité** avec peu de cibles macro |
| Fransen, Langelaar, Schott 2021, Powder Tech. 393 205-218, 10.1016/j.powtec.2021.07.048 | métamodèles DEM (type de métamodèle : [NON VERIFIE]) | — |
| Ji & Karlovšek 2022, IJMST 32(1) 121-136, 10.1016/j.ijmst.2021.11.003 ; Ji & Karlovšek 2022, Eng. Comput. 39 2001-2016, 10.1007/s00366-021-01564-8 | calibration + **analyse d'unicité** ; évolution différentielle optimisée | l'unicité doit être **mesurée** (postérieur/corrélations), pas supposée |
| Du, Liu, Lei, Liu 2023, RMRE 57 2195-2212, 10.1007/s00603-023-03680-x | calibration des micro-paramètres du **FDEM 3D** (méthode : [NON VERIFIE], résumé inaccessible) | à lire pour la phase 3D |
| Zhou, Xu, Gong, Ma 2025, Sci. Rep. 15, 10.1038/s41598-025-99480-0 (résumé vérifié) | cadre « analytique-optimal » BPM pour l'**enveloppe** : ft, c, φ locaux comme variables indépendantes ; estimation grossière par traction directe puis UCS, φ par m_i de Hoek–Brown, puis Adam | même hiérarchie que la nôtre : **ft par la traction, c et φ par l'enveloppe** |
| Ye et al. 2025, IJRMMS 194 106233, 10.1016/j.ijrmms.2025.106233 (titre vérifié : « Mohr-Coulomb strength and FDEM parameter determination of weathered granite via optimized neural network… ») | courbe entière + mode de rupture, NSGA-II (lecture d'août, `calibration_redbohus/README.md:40-43`) | objectifs sur la courbe, pas sur le pic |
| Bu et al. 2026, IJRMMS 199 106400, 10.1016/j.ijrmms.2026.106400 | 3 456 UDEC-BBM ; RF/SVR/GPR/DNN ; **RF gagne en R² direct, échoue en inversion** (`README.md:44-47`, `METHODOLOGIE.md:123-125`) | juger le substitut sur la **tâche inverse** |
| Jiang et al. 2025, Sci. Rep. 15:34923 (« 328 runs suffisent », `README.md:48-51`) | corrélation de Pearson → réduction, stacking | **[NON VERIFIE]** : DOI introuvable par Crossref le 2026-09-02 ; citée seulement d'après le README d'août |
| Aboayanah et al. 2024, RMRE 57 4679-4706, 10.1007/s00603-024-03789-7 | GB-FDEM + DIC ; seuils CI/CD comme cibles | σ_ci/σ_cd comme observables de calibration (déjà prévu, `PLAN_GBM.md`) |

**Fondements statistiques (DOI vérifiés)** : Sacks et al. 1989 (10.1214/ss/1177012413) ; Jones, Schonlau, Welch 1998, EGO/EI (10.1023/A:1008306431147) ; Kennedy & O'Hagan 2001, discrépance (10.1111/1467-9868.00294) ; Brynjarsdóttir & O'Hagan 2014 (10.1088/0266-5611/30/11/114007) ; Higdon et al. 2008, PCA + GP pour sorties de grande dimension = courbes (10.1198/016214507000000888) ; Loeppky, Sacks, Welch 2009, **n = 10·d** (10.1198/TECH.2009.08040) ; Gramacy & Lee 2012, nugget (10.1007/s11222-010-9224-x) ; Ankenman, Nelson, Staum 2010, krigeage stochastique/réplicats (10.1287/opre.1090.0754) ; Binois, Gramacy, Ludkovski 2018, hetGP (10.1080/10618600.2018.1458625) ; Picheny, Wagner, Ginsbourger 2013, critères d'enrichissement **sous bruit** (10.1007/s00158-013-0919-4) ; Ginsbourger, Le Riche, Carraro 2010, lots q-EI / « Kriging believer » (10.1007/978-3-642-10701-6_6) ; Srinivas et al. 2012, GP-UCB (10.1109/TIT.2011.2182033) ; Saltelli et al. 2010, indices de Sobol (10.1016/j.cpc.2009.09.018) ; Morris 1991, criblage (10.1080/00401706.1991.10484804) ; Rasmussen & Williams 2006 (10.7551/mitpress/3206.001.0001).

---

## 4. Choix du substitut et rôle de chaque brique

**4.1 Processus gaussien (krigeage) — substitut principal.** Régime d = 5–9, n = 100–250 : c'est *le* régime du GP. Il fournit la variance de prédiction, indispensable à l'enrichissement et au postérieur. Noyau Matérn 5/2 **ARD** (une longueur par paramètre = criblage gratuit, `METHODOLOGIE.md:276-294`) + `WhiteKernel` ou `alpha` par point (bruit de graine). Entrées normalisées [0,1] après transformation log des paramètres d'échelle. **Un GP par descripteur** (6 descripteurs × 2 confinements = 12 GP, ou 6 GP avec σ₃ en entrée supplémentaire — recommandé : σ₃ en entrée, cela lie les deux confinements et prépare la prédiction ; mais **ne jamais extrapoler le GP à 75/100** : la prédiction se fait par runs directs, §7).

**4.2 Forêt aléatoire / ExtraTrees — contrôle.** Sans hypothèse de lissité, robuste aux runs « locked », mais constante par morceaux, sans extrapolation, incertitude non calibrée. Usage : (i) comparaison LOO des erreurs de prédiction ; (ii) importance par permutation croisée avec les longueurs ARD ; (iii) si la RF bat nettement le GP en LOO sur un descripteur, c'est le signe d'une discontinuité (seuil de rupture) → ajouter un **classifieur** « rompt / ne rompt pas » (GP-classifier ou RF) et n'entraîner le GP de régression que sur les runs rompus.

**4.3 Alternatives écartées** : réseau de neurones (Benvenuti 2016, Bu 2026) — trop de points ; PCE — pas d'avantage à d ≤ 9 avec réponses non lisses ; RSM quadratique (Yoon 2007) — trop rigide.

**4.4 APSO** — optimiseur global de (a) la fonction d'acquisition à chaque lot, (b) J sur la moyenne du GP au final, (c) l'initialisation du MCMC. Alternative de secours dans scipy : `differential_evolution` (à faire tourner en parallèle pour confirmer le même optimum).

**4.5 Émulateur de courbe (option, étape 2)** : PCA sur les courbes interpolées sur la grille cible (4–5 modes ≈ 99 %) + GP par coefficient (Higdon 2008 ; déjà prévu `CAMPAGNE.md:134-138`). À faire **après** le pipeline à descripteurs, qui est plus robuste aux runs partiels.

---

## 5. Tailles, bruit, budgets

**5.1 Dimensions.**

| calibration | paramètres calibrés | d | figés |
|---|---|---|---|
| 1 homogène | ln c, tan φ, ln(ft/c) *(ou ft figé par BTS → d = 5)*, ln ℓ_cz, ln(G_II/G_I), μ_res | **6** | E = 77,7·(1−ν²) ≈ 72,8 GPa (analytique, `calib_quick/README.md:87-88`), ν 0,25, ρ, pénalité 20, ξ 0,01, h 0,8 mm, T 2,2 ms, pullV, damping |
| 2 Weibull | m ∈ [3 ; 12], ℓ_corr ∈ {0} ∪ [0,5 ; 3] mm (`jointWeibullM`, `strengthCorrLength`, `DOCUMENTATION_rockim.md` §5.4) | **2** (+ 6 hérités, gelés au MAP du cas 1) | cible spécifique : **largeur de bande inter-réplicats** + σ_ci ; 3 graines de champ (`fieldSeed`) par point |
| 3 GBM | s_ft, s_c (facteurs globaux sur les résistances de phase), tan φ, α_gb (= gbAlphaTen = gbAlphaCoh), gbHeteroFactor, ln ℓ_cz, μ_res | **7** | propriétés de phase (Table 2 d'Aboayanah, `PLAN_GBM.md`), fractions 62/31/7, grainSize 3 mm ; `grainSizeSpread` ∈ {0 ; 0,3 ; 0,6} traité comme **facteur** (LHS à 0,3 ; 0 et 0,6 rejoués au MAP avec 3 graines = 12 runs), tailles par phase selon la lame mince (`calib_quick/README.md:138-142`) |

**5.2 Plan initial.** Règle de Loeppky–Sacks–Welch : n₀ = 10·d → **60** (d = 6), **50** (d = 5), **70** (d = 7). LHS : `scipy.stats.qmc.LatinHypercube(d, scramble=True, optimization="random-cd", rng=…)` (critère de discrépance centrée ≈ maximin). Le plan d'août (44 points en 6D) était sous la règle et a donné R² croisé 0,44 sur l'UCS (`tools/enrich.py:9-10`) — d'où n₀ = 10d **minimum**. Coût par point = 2 runs (σ₃ = 20 et 50) ; le run à 20 MPa coûte plus cher (plus de joints rompus, `CAMPAGNE.md` : 1 817 s vs 1 529 s en f1) [NON VERIFIE pour f2].

**5.3 Bruit de réalisation — trois sources et comment les traiter.**
- *Sources* : (i) maillage (cas 1-2 : un `.msh` par graine — générer 5 maillages `box20x40_h08_algo5_s<k>.msh` avec le même champ de taille bruité ±10 % et des graines différentes) ; (ii) champ Weibull (`fieldSeed`, indépendant du maillage, `DOCUMENTATION_rockim.md:190`) ; (iii) Voronoï/phases (`seed`, ligne 124). Le GBM à **113 grains** est le plus bruité (réserve « nombre de grains », `CAMPAGNE.md:392-408`) — à mesurer, pas à supposer.
- *Mesure* : au **centre** du domaine, 5 graines ; en 4 points LHS, 3 graines → 17 points-réplicats × 2 confinements = **34 runs**. On en tire σ²_rep,k par descripteur (et sa dépendance au niveau via un modèle σ_rep ∝ y si visible).
- *Substitut* : krigeage **stochastique** (Ankenman 2010) — la valeur d'entraînement est la moyenne des réplicats, avec variance de bruit σ²_rep/n_rep pour les points répliqués et σ²_rep pour les autres ; en scikit-learn : `GaussianProcessRegressor(alpha=array_par_point)` (variance connue) **ou** `+ WhiteKernel()` (nugget homoscédastique appris ; Gramacy & Lee 2012 : garder un nugget même pour un code déterministe). Si σ_rep varie fortement dans le domaine (attendu au GBM près du seuil de rupture) → hetGP (Binois 2018) : approximation simple = second GP sur ln σ²_rep des réplicats [à coder à la main, pas de paquet].
- *Optimisation sous bruit* : viser le minimum de la **moyenne** du GP, pas de l'échantillon (Picheny 2013) ; le « meilleur run » brut n'est jamais le résultat.
- *Ordre de grandeur attendu* : CoV pic ≈ 1–3 % (BTS d'août : 1,6 %) pour les cas 1-2 ; 3–8 % pour le GBM à 113 grains [NON VERIFIE] ; le pic expérimental est à ±0,5 % → c'est le bruit numérique qui fixe la tolérance minimale sur le pic (τ_pic ≥ 2 σ_rep).

**5.4 Enrichissement séquentiel.** Critère : **LCB sur J** (Srinivas 2012), a(x) = Ĵ(x) − κ σ̂_J(x), κ = 2 (exploration) puis 1 (dernier lot) ; Ĵ et σ̂_J par propagation Monte-Carlo (200 tirages des GP) — l'EI de Jones 1998 est mal posé sous bruit (Picheny 2013), la LCB reste saine. Lots de **q = 8** points par la stratégie « Kriging believer » (Ginsbourger 2010 : on ajoute la prédiction comme pseudo-observation et on ré-optimise) + distance minimale 0,15 en coordonnées normalisées entre nouveaux points (comme `enrich.py`, DMIN 0,18). **R = 6 lots** (d = 6) → 48 points. Arrêt anticipé : amélioration de min Ĵ < 2 % sur deux lots **et** σ̂_J < 0,1 au minimum. À chaque lot : re-fit des GP, LOO, tableau ARD.

**5.5 Budgets (runs = points × 2 confinements ; coûts mesurés §1).**

| étape | cas 1 (homogène) | cas 2 (Weibull, au MAP du cas 1) | cas 3 (GBM 3 mm) |
|---|---|---|---|
| ancrage ft par traction/BTS | 3–5 runs (30–90 s) | — | 3 runs |
| réplicats (nugget) | 34 runs | inclus (3 graines/point) | 34 runs |
| LHS n₀ | 60 pts = 120 runs | 12 pts × 3 graines × 2 = 72 runs | 70 pts = 140 runs |
| enrichissement | 6 × 8 = 48 pts = 96 runs | 1 × 8 × 3 × 2 = 48 runs | 5 × 8 = 40 pts = 80 runs |
| validation (§7) | 22 runs | 12 runs | 22 + 12 (spread 0/0,6) = 34 runs |
| contrôle vitesse pullV ÷ 10 au MAP | 2 runs × 10 × coût | — | 2 runs × 10 × coût |
| **total runs** | **≈ 275** | **≈ 130** | **≈ 290** |
| **temps mur** (146–530 s ; 295–530 s GBM) | **11–40 h** (moyenne réaliste ≈ 18 h) | 5–19 h | **24–43 h** |

Soit **2–3 nuits** par calibration à 1 run × 14 threads ; si 2 × 7 threads donne > 1,3× de débit (à mesurer), diviser par ~1,3.

**5.6 Pourquoi pas plus, pourquoi pas moins.** Moins : sous 10·d le GP ne voit pas les interactions (août). Plus : au-delà de ~200 points en d = 6, la précision du GP est limitée par le nugget, pas par n — le gain marginal passe dans les réplicats et la validation.

---

## 6. Fonction objectif sur q(ε), normalisation, identifiabilité

**6.1 Descripteurs (même opérateur appliqué à la courbe simulée et à la courbe cible).**
- Prétraitement sim : q = sigma − offset de consolidation (moyenne juste avant `pullDelay`), ε = déplacement imposé (rampe cosinus analytique)/H (`calib_quick/plot_quick.py:56-70`) ; lissage Savitzky–Golay (fenêtre ≈ 2·10⁻⁴ de déformation, ordre 3) avant toute détection de pic (les oscillations à ε̇ ≈ 5 s⁻¹ sinon biaisent le max). ⚠️ ε est une déformation **globale** (mors) ; en août `epsGauge` (extensomètre intérieur) donnait 0,66 % là où `epsSpec` donnait 0,55 % (`tools/add_epspeak.py`, docstring). Si la cible expérimentale `eps_axial` est une jauge locale, il y a un biais métrologique de 10–20 % sur ε_pic **[NON VERIFIE — à trancher avant de pondérer ε_pic]** ; le « déficit systématique de ε_pic ~20 % » (`CAMPAGNE.md:430`) pourrait en être une part.
- k = 1 **q_pk** (pic lissé) ; k = 2 **ε_pk** ; k = 3 **E_app** = pente 20–50 % du pic (`simcurve.py:105-110`, `plot_quick.py:73-76`) ; k = 4 **r_ci = σ_ci/q_pk**, σ_ci = q au premier point où la courbe s'écarte de 2 % de la tangente élastique (opérateur identique exp/sim ; contrôle croisé avec CI_frac SBM 0,55–0,57) ; k = 5 **chute** = 1 − q(ε_end,commun)/q_pk (ε_end,commun du JSON : 0,764 % à 20, 0,941 % à 50 — la fenêtre T = 2,2 ms → 1,19 % les couvre) ; k = 6 **RMSE_bande** = √mean(((q_sim − q_exp)/σ_eff)²) sur la grille cible jusqu'à ε_end, avec **σ_eff = √(σ_exp² + (0,03·q_pk,exp)²)** — le plancher expérimental de 5 MPa (`targets…json` `_meta`) rend la RMSE de `simcurve.py:176-179` explosive pour tout modèle à 5 % près ; le terme 3 % est la discrépance de modèle (Kennedy & O'Hagan).
- Runs sans rupture (ex. `q1s035` : 3 rompus) : q_pk = q_end, chute = 0, drapeau `locked = nBroken_end − nBroken_0 < 20` ; on **garde** le point (frontière du domaine valide, leçon d'août `METHODOLOGIE.md:181-183`).

**6.2 Normalisation et pondérations.** Erreur normalisée e_k = (y_k − y_k*)/τ_k ; tolérances = ce que l'on accepterait physiquement (exp ⊕ bruit de graine ⊕ discrépance) :

| k | τ_k | justification |
|---|---|---|
| q_pk | 0,03·q_pk,exp (12 / 18 MPa) | exp ±0,5 % mais σ_rep 1–3 % et discrépance |
| ε_pk | 0,10·ε_pk,exp | métrologie incertaine (§6.1) |
| E_app | 0,05·E_exp ≈ 3,9 GPa | dispersion 75–78 GPa des essais |
| r_ci | 0,05 (absolu) | SBM : 0,547–0,569 |
| chute | 0,08 (absolu) | exp 0,30/0,26 ; post-pic lissé par ε̇ (réserve 1, `CAMPAGNE.md:384-391`) |
| RMSE_bande | 1 | déjà en unités de σ_eff |

**J(x) = Σ_{σ₃∈{20,50}} ½ Σ_k w_k e_k²**, w = (q_pk 0,30 ; RMSE 0,20 ; ε_pk 0,15 ; r_ci 0,15 ; E 0,10 ; chute 0,10), Σw = 1 → **J ≈ 1 = tout à la tolérance**, J < 0,5 = excellent. Raisons : le pic est la cible la plus reproductible et la plus « chère » physiquement ; la RMSE porte la forme (demande de Fernando : « la courbe entière ») mais recouvre E et pic, d'où 0,20 ; r_ci et ε_pk portent l'arc σ_cc→σ_ci→σ_cd que les pics ne voient pas ; E à 0,10 seulement parce qu'il est fixé analytiquement dans le cas 1 (il redevient informatif au GBM : complaisance des joints de grain). Pour le postérieur, la même somme sert de −2·log-vraisemblance (τ_k = écart-type total).

Variante multi-objectif (Ye 2025) : garder F₁ = (q_pk, E, RMSE), F₂ = (ε_pk, r_ci, chute) séparés et tracer le front sur le substitut (gratuit : tirage de 2·10⁵ points comme `analyze.py:127-140`) — utile pour **montrer** le compromis arc↔pic (`CAMPAGNE.md:422-425`) plutôt que de le cacher dans les poids.

**6.3 Identifiabilité : ce que deux confinements peuvent et ne peuvent pas dire.**
- *Comptage* : l'enveloppe des pics est quasi linéaire entre 20 et 50 : q_pk = q₀ + k·σ₃ → **deux nombres** (q₀ ≈ 274 MPa, k ≈ 6,5, `CAMPAGNE.md:413`, calcul §7.1). k ↔ tan²(45+φ/2) − 1 (k = 6,5 ↔ φ ≈ 50°) ; q₀ ↔ combinaison (c, ft) : à ℓ_cz constant, pic ∝ (ft·c)^0,7 (`calib_quick/README.md:115`). → **vallée (ft, c)** : toute paire à produit constant donne le même pic. Wang et al. 2021 (10.1016/j.apt.2020.12.015) le démontrent pour un PSO-DEM : « different sets… when few macroparameters are used ».
- *Comment la lever avec ce qu'on a* :
  1. **σ_ci** (onset) : l'insertion en cisaillement obéit à σ₁ ≥ σ₃·tan²(45+φ/2) + 2c·tan(45+φ/2) (`CAMPAGNE.md:222-231`, vérifié à 0,02 %) → l'onset à chaque σ₃ est une **relation supplémentaire sur (c, φ)** indépendante de ft. Avec q₀, k, onset₂₀, onset₅₀ : (c, φ) surdéterminés → identifiables ; ft reste faible.
  2. **ft par la traction** : un run `scenario = tension`, `pullV > 0` (`verifyFt`, `DOCUMENTATION_rockim.md` §5.7) ou un brésilien (`scenario = brazilian`, 26–45 s en août, `phaseA_results.csv`) relie ft au BTS mesuré 10,27 ± 0,98 MPa. Rapport BTS/ft du maillage isotrope h = 0,8 mm : [NON VERIFIE] — 3 runs (ft = 8, 12, 16 MPa) suffisent à l'établir. ⚠️ Au GBM d'août le BTS suivait **c** (14 → 37 MPa pour c 10 → 90) plus que ft (28 → 32) (`screen_results.csv`) : refaire ce contrôle sur chaque représentation.
  3. **UCS** : donne (c, φ) à σ₃ = 0 mais la valeur expérimentale est incertaine (126,6 ± 21,4 vs « ~175 ») et l'enveloppe est concave → l'UCS n'est **pas** compatible avec une MC linéaire calée sur 20/50 (ordonnée 274 vs UCS ≤ 175). À utiliser comme **contrôle de forme**, pas comme cible.
  4. ℓ_cz ↔ ε_pk et raideur de la chute ; G_II/G_I ↔ chute en cisaillement ; μ_res ↔ résiduel et chute en confiné (+315 MPa/unité de μ, `CAMPAGNE.md:317-330`) : ces trois-là sont identifiés par ε_pk, chute, RMSE — pas par les pics.
- *Diagnostic dans le pipeline* : (i) longueurs ARD ; (ii) indices de Sobol totaux sur le GP (`scipy.stats.sobol_indices`, Saltelli 2010) ; (iii) **matrice de corrélation du postérieur** et valeurs propres de sa covariance en coordonnées normalisées : un rapport λ_max/λ_min > 100 désigne une direction « molle » (attendue : ln c + 0,7·ln ft ≈ const si ft n'est pas ancré) ; (iv) profils de J le long de la vallée. C'est l'apport revendiqué en août (« postérieur + corrélations », `README.md:57-59`).

---

## 7. Validation à 75 et 100 MPa, intervalles

**7.1 Ce qu'il faut attendre.** Pente expérimentale 20→50 : (599,2 − 404,0)/30 = **6,51** ; 50→100 : **3,98** (`targets_triax_bohus.json`). Une enveloppe linéaire calée sur 20/50 prédit 761,9 MPa à 75 (**+8,4 %**) et 924,6 à 100 (**+15,8 %**). Avec le bulk élastique (`crushCap = 1e12`, seule option compatible avec `phases`, `calib_quick/README.md:19-20`), rien dans le modèle ne courbe l'enveloppe — sauf, au GBM, les joints hétérophases et la biotite [NON VERIFIE que cela suffise]. L'extrapolation d'août (+42/+53 %, `points_results.csv`) montre que l'écart peut être bien plus grand quand φ est trop haut. Verdict attendu : **biais positif de 8–16 % à 100 MPa = erreur de forme de modèle**, à assumer explicitement (`CAMPAGNE.md:412-415, 443-445`) ; ce n'est pas un échec de calibration. Levier de forme (hors calibration) : un plafond concave du bulk (`law = mc` en homogène, incompatible avec `phases`), ou une enveloppe de joint courbe — capacité à ajouter par addition si Fernando le décide.

**7.2 Protocole.** (a) MAP × 3 graines × {20, 50} = 6 runs → chaque descripteur dans τ_k ? (b) **Prédiction pure** : MAP × 3 graines × {75, 100} = 6 runs + 8 tirages du postérieur × 1 graine × {75, 100} = 16 runs → **22 runs**. (c) Intervalle prédictif à 90 % = quantiles 5/95 des 8 + 3 courbes par confinement (postérieur ⊕ graine), superposé à la bande exp ±1σ, descripteur par descripteur. Critères : pic exp dans l'intervalle **ou** |biais| ≤ 10 % (75) / ≤ 20 % (100), en cohérence avec le seuil « prédiction ±20 % » d'août (`METHODOLOGIE.md:105`). (d) Contrôle de vitesse au MAP : pullV ÷ 10 à σ₃ = 50 (2 runs, ~25–90 min chacun) → si le pic bouge de > 3 %, les paramètres portent l'étiquette ε̇ ≈ 5 s⁻¹.

**7.3 Ce que les intervalles disent.** Un intervalle étroit qui rate la cible = discrépance de forme (courbure) ; un intervalle large = paramètres mal contraints (vallée) — les deux sont des résultats de thèse, à condition d'être séparés, ce que le postérieur avec terme de discrépance (Brynjarsdóttir & O'Hagan 2014) permet.

---

## 8. Pipeline prêt à coder (Python 3.12, numpy/scipy/scikit-learn/joblib présents)

**8.1 Arborescence proposée** (`rockim_f2/calib_pipeline/`, un module par brique, tout redémarrable) :

```
calib_pipeline/
  space.py       # paramètres, bornes, transformations log, reparamétrisation (c, tanφ, ft/c, ℓcz, GII/GI, μres)
  decks.py       # écriture des .cfg depuis q1_homog_P050.cfg / q3_gbm_P050.cfg (points décimaux, clés exactes)
  runner.py      # file séquentielle 1 job × 14 threads (ou 2 × 7), reprise sur run.log 'wall time', OMP_NUM_THREADS
  features.py    # q(ε) depuis history.csv → 6 descripteurs + flag locked (opérateur commun exp/sim)
  objective.py   # J(x) pondéré, tolérances τ_k
  surrogate.py   # GP Matérn 5/2 ARD + bruit, LOO, RF de contrôle, ARD/Sobol
  apso.py        # APSO Zhan 2009 (numpy)
  enrich.py      # LCB + Kriging believer + distance minimale → lot de 8
  posterior.py   # Metropolis adaptatif sur le substitut, corrélations, valeurs propres
  validate.py    # MAP × graines, prédiction 75/100, intervalles, figures PDF Computer Modern
  base.csv       # une ligne par run : tag, seed, σ3, x (6-9 colonnes), descripteurs, locked, wall_s
```

**8.2 Espace et plan initial**
```python
import numpy as np
from scipy.stats import qmc
# bornes en coordonnées physiques ; log pour les échelles
SPACE = {"cohesion":(8e6, 40e6,"log"), "tanPhi":(0.8, 1.5,"lin"),      # φ 39–56°
         "ftOverC":(0.3, 0.8,"log"), "lcz":(0.05, 0.6,"log"),           # ℓcz [m] = E·Gf/ft²
         "gfShearFactor":(2, 30,"log"), "jointResidualMu":(0.2, 0.7,"lin")}
KEYS = list(SPACE)
def to_unit(x):   # physique -> [0,1]
    u=[]; 
    for k,v in zip(KEYS,x):
        lo,hi,s=SPACE[k]; u.append((np.log(v)-np.log(lo))/(np.log(hi)-np.log(lo)) if s=="log" else (v-lo)/(hi-lo))
    return np.array(u)
def to_phys(u):
    x=[]
    for k,ui in zip(KEYS,u):
        lo,hi,s=SPACE[k]; x.append(np.exp(np.log(lo)+ui*(np.log(hi)-np.log(lo))) if s=="log" else lo+ui*(hi-lo))
    return np.array(x)
def deck_values(x, E=72.8e9):     # -> clés rockim
    c,tanphi,r,lcz,gsf,mu = x
    ft = r*c; Gf = lcz*ft**2/E
    return {"cohesion":c,"frictionDeg":np.degrees(np.arctan(tanphi)),"ft":ft,"Gf":Gf,
            "gfShearFactor":gsf,"jointResidualMu":mu,"E":E}
U0 = qmc.LatinHypercube(d=len(KEYS), scramble=True, optimization="random-cd", rng=20260902).random(n=60)
```

**8.3 Descripteurs (opérateur commun exp/sim)**
```python
from scipy.signal import savgol_filter
def features(eps, q, tgt):            # eps [-], q [MPa] (offset déjà retiré), tgt = bloc JSON du confinement
    qs = savgol_filter(q, 21, 3)
    ipk = int(np.argmax(qs)); qpk = qs[ipk]; epk = eps[ipk]
    m = (qs[:ipk+1] > .2*qpk) & (qs[:ipk+1] < .5*qpk)
    E, b = np.polyfit(eps[:ipk+1][m], qs[:ipk+1][m], 1)            # MPa
    dev = (E*eps[:ipk+1]+b) - qs[:ipk+1]
    ici = np.argmax(dev > .02*qpk) if np.any(dev > .02*qpk) else ipk
    rci = qs[ici]/qpk
    eend = tgt["eps_end_common_microstrain"]*1e-6
    qend = np.interp(eend, eps, qs); drop = 1 - qend/qpk
    ge = np.array(tgt["eps_grid_microstrain"])*1e-6; gq = np.array(tgt["q_mean_MPa"]); gs = np.array(tgt["q_std_MPa"])
    k = ge <= min(eend, eps.max())
    seff = np.hypot(gs[k], .03*tgt["q_peak_mean_MPa"])
    rmse = np.sqrt(np.mean(((np.interp(ge[k], eps, qs) - gq[k])/seff)**2))
    return dict(qpk=qpk, epk=epk, E=E/1e3, rci=rci, drop=drop, rmse=rmse)   # E en GPa
```
Cible : appliquer `features` à (`eps_grid`, `q_mean_MPa`) du JSON → y* par confinement (ε_pk* ≈ 0,65/0,95 %, r_ci* à comparer aux 0,55–0,57 SBM).

**8.4 Objectif**
```python
W   = dict(qpk=.30, rmse=.20, epk=.15, rci=.15, E=.10, drop=.10)
def tol(ystar): return dict(qpk=.03*ystar["qpk"], epk=.10*ystar["epk"], E=.05*ystar["E"], rci=.05, drop=.08, rmse=1.0)
def J(y, ystar):                       # y, ystar : {σ3: {descripteur: valeur}}
    tot = 0.0
    for s3 in (20, 50):
        t = tol(ystar[s3])
        tot += .5*sum(W[k]*((y[s3][k]-ystar[s3][k])/t[k])**2 for k in W)
    return tot
```

**8.5 Substitut GP par descripteur, bruit connu, LOO, RF de contrôle**
```python
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C, WhiteKernel
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import cross_val_predict, KFold
def fit_gp(U, y, var_noise):           # U (n,d+1) avec σ3/100 en dernière colonne ; var_noise (n,) = σ²rep/n_rep
    ker = C(1.0,(1e-2,1e2))*Matern(length_scale=np.ones(U.shape[1]), length_scale_bounds=(5e-2,1e2), nu=2.5) \
          + WhiteKernel(1e-3,(1e-8,1.0))              # nugget résiduel EN PLUS du bruit connu
    gp = GaussianProcessRegressor(ker, alpha=var_noise, normalize_y=True, n_restarts_optimizer=8, random_state=0)
    return gp.fit(U, y)
def loo_score(model_fn, U, y):
    yp = cross_val_predict(model_fn(), U, y, cv=KFold(10, shuffle=True, random_state=0))
    return 1 - np.sum((y-yp)**2)/np.sum((y-y.mean())**2)
# descripteurs transformés : log(qpk), epk, E, rci, drop, log1p(rmse) ; un GP par descripteur, σ3 en entrée
# contrôle : ExtraTreesRegressor(500, min_samples_leaf=2) même LOO ; importance par permutation vs 1/longueur ARD
```
Criblage : après le LHS, longueur ARD ≫ 10 (en unités [0,1]) sur tous les descripteurs **et** indice de Sobol total < 0,02 (`scipy.stats.sobol_indices` sur la moyenne du GP) → paramètre figé, déclaré, et le GP ré-entraîné sans lui.

**8.6 APSO (Zhan et al. 2009) sur le substitut**
```python
def apso(fobj, d, lo, hi, N=40, G=300, rng=np.random.default_rng(0)):
    X = lo + rng.random((N,d))*(hi-lo); V = np.zeros((N,d)); vmax = .2*(hi-lo)
    F = np.array([fobj(x) for x in X]); P = X.copy(); Fp = F.copy(); g = np.argmin(F); gb, fg = P[g].copy(), Fp[g]
    w, c1, c2, state = .9, 2.0, 2.0, 1
    def memb(f):                                       # éq. 9a-d
        s1 = np.interp(f,[0,.4,.6,.7,.8,1],[0,0,1,1,0,0]); s2 = np.interp(f,[0,.2,.3,.4,.6,1],[0,0,1,1,0,0])
        s3 = np.interp(f,[0,.1,.3,1],[1,1,0,0]);         s4 = np.interp(f,[0,.7,.9,1],[0,0,1,1])
        return np.array([s1,s2,s3,s4])
    for gen in range(G):
        D = np.array([np.mean(np.linalg.norm((X[i]-X)/(hi-lo),axis=1)*d**-.5) for i in range(N)])*N/(N-1)
        f = (D[g]-D.min())/max(D.max()-D.min(),1e-12)          # éq. 8
        mu = memb(f); cand = np.flatnonzero(mu > 0) + 1
        if len(cand)==1: state = cand[0]
        elif state in cand: pass                                # stabilité
        else: state = next(s for s in [1,2,3,4,1] if s in cand and s == state%4+1) if any(s==state%4+1 for s in cand) else cand[np.argmax(mu[cand-1])]
        w = 1/(1+1.5*np.exp(-2.6*f))                           # éq. 10
        dl = rng.uniform(.05,.1)                               # éq. 11
        dc1, dc2 = {1:(dl,-dl), 2:(.5*dl,-.5*dl), 3:(.5*dl,.5*dl), 4:(-dl,dl)}[state]
        c1 = np.clip(c1+dc1,1.5,2.5); c2 = np.clip(c2+dc2,1.5,2.5)
        if c1+c2 > 4.0: c1, c2 = 4.0*c1/(c1+c2), 4.0*c2/(c1+c2)  # éq. 12
        if state == 3:                                          # ELS, éq. 13-14
            sig = 1.0 - .9*gen/G; k = rng.integers(d); Pn = gb.copy()
            Pn[k] = np.clip(Pn[k] + (hi[k]-lo[k])*rng.normal(0,sig), lo[k], hi[k]); fn = fobj(Pn)
            if fn < fg: gb, fg = Pn, fn
            else: worst = np.argmax(F); X[worst], F[worst] = Pn, fn
        r1, r2 = rng.random((N,d)), rng.random((N,d))
        V = np.clip(w*V + c1*r1*(P-X) + c2*r2*(gb-X), -vmax, vmax); X = np.clip(X+V, lo, hi)
        F = np.array([fobj(x) for x in X]); better = F < Fp; P[better], Fp[better] = X[better], F[better]
        g = np.argmin(Fp)
        if Fp[g] < fg: gb, fg = P[g].copy(), Fp[g]
    return gb, fg
```
(À tester d'abord sur Rastrigin/Rosenbrock en d = 6 contre `scipy.optimize.differential_evolution` ; 20 redémarrages ; sur le substitut, fobj = Ĵ ou la LCB.)

**8.7 Enrichissement par lot (LCB + Kriging believer)**
```python
def lcb_factory(gps, ystar, kappa=2.0, ns=200, rng=np.random.default_rng(1)):
    def a(u):
        ys = {s3:{} for s3 in (20,50)}; Js = []
        draws = {(s3,k): gps[k].sample_y(np.r_[u, s3/100][None,:], n_samples=ns, random_state=rng.integers(1e9)).ravel()
                 for s3 in (20,50) for k in W}
        for i in range(ns):
            y = {s3:{k: inv_transform(k, draws[(s3,k)][i]) for k in W} for s3 in (20,50)}
            Js.append(J(y, ystar))
        Js = np.array(Js); return Js.mean() - kappa*Js.std()
    return a
def next_batch(gps, U_train, q=8, dmin=0.15):
    batch = []
    for _ in range(q):
        u, _ = apso(lcb_factory(gps, YSTAR), d, np.zeros(d), np.ones(d))
        if all(np.linalg.norm(u-v) > dmin for v in list(U_train)+batch): batch.append(u)
        gps = believer_update(gps, u)          # ajoute la prédiction moyenne comme pseudo-observation (Ginsbourger 2010)
    return batch
```

**8.8 Postérieur (Metropolis adaptatif sur le substitut, comme `analyze.py:153-175` mais avec τ_k)**
```python
def logpost(u, gps, ystar):
    if np.any(u<0) or np.any(u>1): return -np.inf
    L = 0.0
    for s3 in (20,50):
        t = tol(ystar[s3])
        for k in W:
            m, sd = gps[k].predict(np.r_[u, s3/100][None,:], return_std=True)
            m = inv_transform(k, m[0]); sdk = sd_transform(k, m, sd[0])
            L += -.5*((m-ystar[s3][k])/np.hypot(t[k], sdk))**2       # tolérance ⊕ incertitude d'émulation
    return L
# chaîne : 40 000 pas, pas adapté à 25 % d'acceptation, burn-in 8 000, thin 5 ; 4 chaînes, R̂ < 1,05
# sorties : médiane/IC95 par paramètre, matrice de corrélation, valeurs propres de cov (identifiabilité), 8 tirages pour §7
```

**8.9 Ordre d'exécution et jalons (chaque lancement validé par Fernando : clés + coût)**
0. `runner.py bench` : 2 runs identiques 1×14 vs 2×7 → débit.
1. Ancrage ft : 3 runs traction/brésilien (ft 8/12/16 MPa) → relation BTS(ft) ; fixer ft ou borner ft/c.
2. Réplicats : centre × 5 graines + 4 points × 3 graines, 2 confinements (34 runs) → σ_rep par descripteur.
3. LHS 60 × 2 (120 runs, ~1 nuit). Fit, LOO (GP vs ExtraTrees), ARD, Sobol → figer ce qui est inerte.
4. Enrichissement 6 lots × 8 × 2 (96 runs, ~1 nuit) avec re-fit à chaque lot ; arrêt anticipé possible.
5. APSO sur Ĵ → MAP ; MCMC → postérieur, corrélations ; figures (PDF, Computer Modern).
6. Validation : 22 runs (§7.2) + contrôle pullV ÷ 10 (2 runs longs).
7. Cas 2 (Weibull) au MAP du cas 1 : LHS 12 × 3 graines × 2 (72 runs) sur (m, ℓ_corr), cible = bande + r_ci ; 1 lot d'enrichissement ; validation 12 runs.
8. Cas 3 (GBM) : mêmes étapes avec d = 7 (§5.1), départ centré sur le MAP homogène transposé (c, φ, ℓ_cz, μ_res), n₀ = 70, 5 lots, spread {0 ; 0,6} au MAP, contrôle à 2 mm/40 × 80 mm (3 graines, 2 confinements) pour la réserve « nombre de grains ».
9. Mise à jour PHD.md + phd/FDEM.md ; copie vers la base ; commit.

---

## 9. Points à trancher / [NON VERIFIE] récapitulatif

1. **UCS de référence** : 126,6 ± 21,4 (`calibration_redbohus/README.md:20`) ou « ~175 » (`CAMPAGNE.md:413`) — [NON VERIFIE].
2. **Nature de `eps_axial` expérimental** (globale vs jauge locale) → poids et tolérance de ε_pic — [NON VERIFIE].
3. **Rapport BTS/ft** sur le maillage isotrope h = 0,8 mm et sur le GBM 3 mm — [NON VERIFIE], 3 runs chacun.
4. **Bruit de graine** des trois représentations — à mesurer (17 points-réplicats) ; l'estimation « 3–8 % au GBM » est une attente, pas une mesure.
5. **Débit 2 × 7 vs 1 × 14** — [NON VERIFIE].
6. **Coût relatif de `grainSizeSpread = 0,3`** (dt suit hmin ; ÷ 3,6 mesuré à 0,8, `DOCUMENTATION_rockim.md` §5.16) — [NON VERIFIE] à 0,3.
7. `stopPeakDrop` n'est documenté qu'en fdem3d (`DOCUMENTATION_rockim.md:722`) — en 2D, garder T fixe (1,19 % couvre ε_end).
8. `jointPrebrokenFrac` (validé en 2D dans rockim_f1) est-il dans `rockim_f2j/f2l.exe` ? — [NON VERIFIE] ; s'il l'est, c'est un 7ᵉ paramètre candidat pour l'arc (frac ∈ [0,05 ; 0,20], `CAMPAGNE.md:281-285`), au prix d'un couplage arc↔pic déjà documenté.
9. Jiang et al. 2025 (Sci. Rep. 15:34923) — DOI non retrouvé ; à vérifier avant citation dans un livrable.
10. Méthode exacte de Du et al. 2023 (FDEM 3D) et de Fransen et al. 2021 — résumés inaccessibles ; DOI vérifiés seulement.

---

## 10. Références (DOI vérifiés Crossref le 2026-09-02 sauf mention)

- Zhan Z.-H., Zhang J., Li Y., Chung H.S.-H. (2009) Adaptive particle swarm optimization. *IEEE TSMC-B* 39(6) 1362-1381. 10.1109/TSMCB.2009.2015956 — [eprint Glasgow](https://eprints.gla.ac.uk/7645/1/7645.pdf), [PubMed](https://pubmed.ncbi.nlm.nih.gov/19362911/), [notice CityU](https://scholars.cityu.edu.hk/en/publications/adaptive-particle-swarm-optimization/).
- Kennedy J., Eberhart R. (1995) *Proc. ICNN'95* 4, 1942-1948. 10.1109/ICNN.1995.488968. — Shi Y., Eberhart R. (1998) *IEEE ICEC* 69-73. 10.1109/ICEC.1998.699146. — Clerc M., Kennedy J. (2002) *IEEE TEVC* 6(1) 58-73. 10.1109/4235.985692.
- Jin Y. (2011) *Swarm Evol. Comput.* 1(2) 61-70. 10.1016/j.swevo.2011.05.001. — Regis R.G. (2014) *J. Comput. Sci.* 5(1) 12-23. 10.1016/j.jocs.2013.07.004. — Sun C. et al. (2017) *IEEE TEVC* 21(4) 644-660. 10.1109/TEVC.2017.2675628.
- Sacks J. et al. (1989) *Stat. Sci.* 4(4). 10.1214/ss/1177012413. — Jones D.R. et al. (1998) *J. Global Optim.* 13(4) 455-492. 10.1023/A:1008306431147. — Kennedy M.C., O'Hagan A. (2001) *JRSS-B* 63(3) 425-464. 10.1111/1467-9868.00294. — Brynjarsdóttir J., O'Hagan A. (2014) *Inverse Problems* 30 114007. 10.1088/0266-5611/30/11/114007. — Higdon D. et al. (2008) *JASA* 103(482) 570-583. 10.1198/016214507000000888.
- Loeppky J.L., Sacks J., Welch W.J. (2009) *Technometrics* 51(4) 366-376. 10.1198/TECH.2009.08040. — Gramacy R.B., Lee H.K.H. (2012) *Stat. Comput.* 22(3) 713-722. 10.1007/s11222-010-9224-x. — Ankenman B., Nelson B.L., Staum J. (2010) *Oper. Res.* 58(2) 371-382. 10.1287/opre.1090.0754. — Binois M., Gramacy R.B., Ludkovski M. (2018) *JCGS* 27(4) 808-821. 10.1080/10618600.2018.1458625. — Picheny V., Wagner T., Ginsbourger D. (2013) *SMO* 48(3) 607-626. 10.1007/s00158-013-0919-4. — Ginsbourger D., Le Riche R., Carraro L. (2010) *Adaptation Learning and Optimization* 131-162. 10.1007/978-3-642-10701-6_6. — Srinivas N. et al. (2012) *IEEE Trans. Inf. Theory* 58(5) 3250-3265. 10.1109/TIT.2011.2182033. — Saltelli A. et al. (2010) *Comput. Phys. Commun.* 181(2) 259-270. 10.1016/j.cpc.2009.09.018. — Sobol' I.M. (2001) *Math. Comput. Simul.* 55 271-280. 10.1016/S0378-4754(00)00270-6. — Morris M.D. (1991) *Technometrics* 33(2) 161-174. 10.1080/00401706.1991.10484804. — Rasmussen C.E., Williams C.K.I. (2006) MIT Press. 10.7551/mitpress/3206.001.0001.
- Yoon J. (2007) *IJRMMS* 44(6) 871-889. 10.1016/j.ijrmms.2007.01.004. — Hanley K.J. et al. (2011) *Powder Technol.* 210(3) 230-240. 10.1016/j.powtec.2011.03.023. — Tatone B.S.A., Grasselli G. (2015) *IJRMMS* 75 56-72. 10.1016/j.ijrmms.2015.01.011. — Benvenuti L., Kloss C., Pirker S. (2016) *Powder Technol.* 291 456-465. 10.1016/j.powtec.2016.01.003. — Rackl M., Hanley K.J. (2017) *Powder Technol.* 307 73-83. 10.1016/j.powtec.2016.11.048. — Coetzee C.J. (2017) *Powder Technol.* 310 104-142. 10.1016/j.powtec.2017.01.015. — Cheng H. et al. (2018) *Granular Matter* 20:11. 10.1007/s10035-017-0781-y. — Cheng H. et al. (2019) *CMAME* 350 268-294. 10.1016/j.cma.2019.01.027. — Cheng H. et al. (2024) GrainLearning, *JOSS* 9(97) 6338. 10.21105/joss.06338. — Do H.Q., Aragón A.M., Schott D.L. (2018) *Adv. Powder Technol.* 29(6) 1393-1403. 10.1016/j.apt.2018.03.001. — Mohajeri M.J., Do H.Q., Schott D.L. (2020) *Adv. Powder Technol.* 31(5) 1838-1850. 10.1016/j.apt.2020.02.019. — Richter C. et al. (2020) *Powder Technol.* 360 967-976. 10.1016/j.powtec.2019.10.052. — Qu T. et al. (2020) *Powder Technol.* 366 527-536. 10.1016/j.powtec.2020.02.077. — Westbrink F. et al. (2021) *Powder Technol.* 379 602-616. 10.1016/j.powtec.2020.10.067. — Wang M., Lu Z., Wan W., Zhao Y. (2021) *Adv. Powder Technol.* 32(2) 358-369. 10.1016/j.apt.2020.12.015. — Fransen M.P., Langelaar M., Schott D.L. (2021) *Powder Technol.* 393 205-218. 10.1016/j.powtec.2021.07.048. — Ji S., Karlovšek J. (2022) *IJMST* 32(1) 121-136. 10.1016/j.ijmst.2021.11.003 ; *Eng. Comput.* 39 2001-2016. 10.1007/s00366-021-01564-8. — Du H., Liu Q., Lei D., Liu H. (2023) *RMRE* 57 2195-2212. 10.1007/s00603-023-03680-x. — Zhou Y., Xu ?, Gong ?, Ma ? (2025) *Sci. Rep.* 15. 10.1038/s41598-025-99480-0. — Wang M. et al. (2025) *Comput. Part. Mech.* 12 541-555. 10.1007/s40571-024-00820-0.
- Ye et al. (2025) *IJRMMS* 194 106233. 10.1016/j.ijrmms.2025.106233. — Bu et al. (2026) *IJRMMS* 199 106400. 10.1016/j.ijrmms.2026.106400. — Jiang et al. (2025) *Sci. Rep.* 15:34923 [NON VERIFIE]. — Dumoulin C., Thenevin I., Kane ?, Rouabhi A., Latham J.-P. (2024) *Geomech. Energy Environ.* 40 100592. 10.1016/j.gete.2024.100592 ; données Zenodo 10.5281/zenodo.10617548. — Aboayanah K.R. et al. (2024) *RMRE* 57 4679-4706. 10.1007/s00603-024-03789-7. — Mahabadi O.K. et al. (2012) *Int. J. Geomech.* 12(6) 676-688. 10.1061/(ASCE)GM.1943-5622.0000216. — Martin C.D., Chandler N.A. (1994) *IJRMMS* 31(6) 643-659. 10.1016/0148-9062(94)90005-1. — Diederichs M.S., Kaiser P.K., Eberhardt E. (2004) *IJRMMS* 41(5) 785-812. 10.1016/j.ijrmms.2004.02.003.

**Fichiers locaux cités** : `C:/Users/fuzquianoalricabi/simulations/FDEM/rockim_f2/calib_quick/README.md` ; `…/rockim_f2/DOCUMENTATION_rockim.md` ; `…/rockim_f2/out_q{1,2,3}_*/run.log` ; `…/rockim_f2/calib_quick/_log_q1t0{25,35,50}.txt` ; `…/rockim_f1/calib_triax3d/CAMPAGNE.md`, `PLAN_GBM.md`, `simcurve.py`, `targets_triax_bohus.json`, `gbm*.log`, `out_scr_base_P0*.log` ; `…/rockim/rockim_p1/calibration_redbohus/{README.md, METHODOLOGIE.md, lhs_results.csv, screen_results.csv, phaseA_results.csv, points_results.csv, tools/analyze.py, tools/enrich.py, tools/campaign.py, tools/add_epspeak.py}` ; `C:/Users/fuzquianoalricabi/simulations/CONTINUUM/calib_bohus_triax/exp_qc/seuils_sbm_bohus.json` ; `C:/Users/fuzquianoalricabi/python/PSO/{pso_calibration.py, PSO_CALIBRATION_DOC.md, objective.py}`.