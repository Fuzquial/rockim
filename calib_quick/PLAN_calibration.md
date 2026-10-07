# PLAN — trois calibrations Red Bohus sur les triaxiaux (2026-09-02, révisé 15:00)

*Décisions de Fernando (13:00) : (1) le bulk le plus pertinent physiquement ;
(2) ajuster sur σ₃ = 20 et 50 MPa, prédire 75 et 100 ; (3) tailles de grains
réelles du Red Bohus. Révision après l'enquête à cinq agents
(`enquete/A1…A5`) : la critique adverse A5 a invalidé deux conventions de
dépouillement (§2) et imposé cinq contrôles supplémentaires (§6). Rien n'est
lancé sans validation de la liste des clés et du coût.*

## 1. Le bulk : grains élastiques, toute l'inélasticité dans les joints

**Décision : bulk élastique par phase (E, ν, ρ du minéral), aucune
plasticité continue dans les grains ; `crushCap` inactif à 20–100 MPa.**

- Régime fragile : la transition fragile-ductile du granite est à
  σ₃ > 300 MPa à froid (Paterson & Wong 2005) ; ce que le Bohus montre
  (σ_ci → σ_cd → pic → chute 26–30 %) est de la microfissuration puis du
  glissement frictionnel — les joints cohésifs et leur frottement.
- C'est le standard vérifié de tous les FDEM-GBM de granite (A2 §1 : Y-Geo,
  Irazu, Y-HFDEM, HOSS ; Mahabadi 2012, Zhao 2015, Abdelaziz 2018,
  Aboayanah 2024, Fukuda 2020) ; il reproduit σ_ci ≈ 0,45 et σ_cd ≈ 0,85
  UCS par le contraste de raideur et les joints hétérophases faibles.
- Une plasticité MC dans le bulk ferait double emploi avec les joints et
  serait non identifiable à deux confinements (A4 §6.3) ; ceux qui l'ont
  ajoutée visaient les grandes déformations, les roches tendres/altérées
  (Ye 2025) ou des pressions ≫ 50 MPa (école KTH sur le Bohus, indentation).
- La preuve locale : `nInserted` (première insertion adaptative = amorçage
  Coulomb) donne CI/pic = 0,62 (homogène), **0,53 (Weibull m = 8, exp 0,55)**,
  0,30 (GBM α = 1) — l'arc pré-pic est déjà porté par les joints.
- Réversible : `law = mc` par phase est une addition de 60–100 lignes (A3
  §iii) si les joints ne décollent pas CD du pic ; à n'utiliser qu'en
  correction bornée pour 75/100 MPa, jamais comme substitut à la fissuration.

## 2. Deux conventions corrigées (A5, gravité 3)

1. **q = sigma − σ₃** (pas sigma − 16,7) : les mors bloqués pendant
   `pullDelay` laissent σ_yy = ν/(1−ν)σ₃ ; ε = 0 quand sigma atteint σ₃.
   Les anciens pics étaient surestimés de 33 MPa à 50 et 13 MPa à 20.
   Durable : consolidation isotrope (pression axiale pendant `pullDelay`,
   clé opt-in à créer, A3 item 7). `calib/extract.py` et `plot_quick.py`
   corrigés.
2. **E physique dans les decks** (77,7 GPa, ν 0,29) et **cibles
   ε × (1−ν²) = 0,916** (déformation plane à σ_xx constant :
   Δε_yy = Δσ_yy(1−ν²)/E). Garder E physique conditionne le transfert vers
   le forage (vitesses d'onde, raideurs de contact).

État corrigé (`q1v070` = ×0,7 des joints de la sonde 4, ℓ_cz 22 mm) :

| σ₃ | pic q (cible) | ε_pic (cible × 0,916) | CI/pic (exp) | CD/pic (exp) | chute |
|---|---|---|---|---|---|
| 50 | 567 (599, −5 %) | 0,83 % (0,87) | 0,67 (0,55) | 0,98 (0,72) | 0,63 (0,26) |
| 20 | 370 (405, −9 %) | 0,50 % (0,60) | 0,56 (0,57) | 0,98 (0,64) | 0,53 (0,30) |

Exposant de sensibilité du pic aux résistances : **0,40** (le pic est
1,45–1,6 × l'amorçage Coulomb : blocage cinématique des triangles, A5 C).

## 3. Observables et pourquoi

| observable | modèle | expérience | τ (tolérance) | poids |
|---|---|---|---|---|
| q_pic | max de q lissé | 404,8 / 599,0 ± 3 | 3 % | 1 |
| ε_pic | ε au pic | 0,65 / 0,95 % × 0,916 | 10 % | 0,5 |
| E_app | pente 20–50 % du pic | 77,7 → 84,8 en 2D | fixé (contrôle) | 0 |
| **CI** | première insertion (`nInserted`) | SBM 0,55–0,57 du pic | 0,05 abs. | 1 |
| CD | première rupture D ≥ 1 (proxy ; SBM vrai après `epsLat`, §7) | 0,62–0,73 | 0,10 | 0,25 |
| chute | à Δε = +0,1 % après pic, **à 20 MPa seulement** (n_rep = 3) ; censurée | 0,17–0,40 | 0,10 | 0,25 |
| RMSE | q(ε) sur la grille, σ_eff = √(σ_exp² + (0,03 q_pic)²) | bande ± 1σ | 1 | 1 |
| **BTS** | brésilien sur disque Gmsh Delaunay Ø40 à méplats 2 × 20° (`bts_v070b.cfg`, `mesh = file` accepté en disque depuis `rockim_f2n.exe`), plateaux 0,1 m/s ; BTS = k_band × σ_t nominal, k_band = rapport centre/nominal mesuré par la jauge élastique du solveur (0,89 ici = correction des méplats) ; fissure amorcée AU CENTRE (20 premiers joints rompus à r/R < 0,46 sur le diamètre chargé) ; ~50 s à 14 threads | 10,27 ± 0,98 MPa (disque plein) | 10 % | 1 |

Post-pic : à 0,25 m/s et dampingLocal 0,7, Cundall dissipe 2 × la cohésion
(A5 R6) → la chute ne pèse qu'après le contrôle C4 et la vitesse C6.

## 4. Les trois calibrations

### 4a. Homogène — 5 paramètres (μ_res sorti : non identifiable, A5 R13)

| paramètre | bornes | échelle | pourquoi |
|---|---|---|---|
| ft | 4 – 20 MPa | log | BTS 10,3 ; ancré par le brésilien |
| c | 8 – 40 MPa | log | niveau de l'enveloppe |
| φ | 35 – 55° | lin | pente de l'enveloppe (exp 6,5 ↔ φ ≈ 50°) |
| G_f | 10 – 60 J/m² | log | K_Ic 1–2 MPa√m → G_f = K²(1−ν²)/E ≈ 12–47 ; ℓ_cz = E G_f/ft² ≥ 3 h |
| G_II/G_I | 2 – 30 | log | Grasselli ≈ 2, rockim 10–20 |

Fixes : E 77,7, ν 0,29, ρ **2640** (A1 §4), μ_res **0,6** (φ_res 31°, frottement de
base du granite), pénalité d'insertion 4 (défaut ; `jointPenaltyFactor` est
INERTE en adaptatif, A3), ξ 0,01, Yan, insertion adaptative, h 0,8 mm,
éprouvette **selon C5**, pullV **selon C6**, dampingLocal **selon C4**.

Plan : LHS 50–60 points (10 d) + 17 réplicats de graine (nugget) +
enrichissement LCB par lots de 8 (Kriging believer) + MAP par APSO sur le
substitut + postérieur MCMC + validation 22 runs (A4 §5, §7). Décks par
`calib/make_decks.py` (G_f dérivé ou direct), file `calib/runner.py`,
observables `calib/extract.py`, substitut `calib/emulator.py`, `calib/apso.py`,
`calib/objective.py`.

### 4b. Weibull — (m ∈ [3 ; 12], ℓ_corr ∈ {0} ∪ [0,5 ; 3] mm), 3 `fieldSeed`
par point, au MAP du 4a ; cibles : **CI** et largeur de bande inter-réplicats
(le pic est insensible, démontré) ; `weibullScope = lcz` à ajouter (Gf ∝ stat²,
A3 item 4) pour disperser à ℓ_cz constant.

### 4c. GBM — **`grainMeshRandom = true` obligatoire** (le Delaunay intra-grain
par défaut est structuré : R6 0,55 ; corrigé en ajout dans `rockim_f2n.exe`,
DOC 5.16 bis — le cas 3 « −20 % » est à refaire : `q3r_gbm_P050.cfg`) ;
propriétés de phase figées (Aboayanah 2024, Table 2), fractions
**62 / 32 / 6** (A1 : quatre comptages EN 12407), ρ 2640 ; calibrés : s_ft, s_c
(facteurs globaux), φ, α_homo (gbAlphaTen = gbAlphaCoh), gbHeteroFactor
(hétérophase 0,14–0,5 en ft, 0,05–0,3 en G_f dans la littérature), G_f ; les
surcharges par paire `gb.<a>.<b>.*` existent (A3 §ii) pour la table à six
paires d'Aboayanah. Tailles (§5) avec `grainSizeSpread` 0,3/0,5/0,7 en
facteur, `phase.biotite.grainSize` posé. ⚠ W/d = 6,7 à 20 × 40 mm (ISRM ≥ 10,
littérature 20–30 grains en largeur) : le GBM se calibre à **40 × 80 mm
minimum** (≈ 450 grains, ~40 min/run) — à budgéter après 4a.

## 5. Pétrographie (A1, vérifiée / estimée)

| phase | fraction d'aire | confiance | d_eq [mm] | confiance |
|---|---|---|---|---|
| feldspath (Kfs 33 + plag 29) | **0,62** | haute (59–63 %) | **3,5** (2,5–4,5) | basse |
| quartz | **0,32** | haute (30–35 %) | **2,5** (1,5–3,5) | basse |
| biotite + opaques | **0,06** | haute (bt 3–6 + op. 1–2) | **1,0** (0,5–1,5) | basse |

Texture : équigranulaire à faiblement porphyrique, grain moyen (classe
1–5 mm), Kfs le plus gros (pluri-mm, phénocristaux ≈ 1 cm occasionnels),
joints de grains nets, microfissures préexistantes dans le quartz, faible
foliation (θ 31–61°, Dumoulin 2024 Table 6). **Aucune distribution de tailles
par minéral n'est publiée** ; le 62/31/7 du dépôt venait d'un docstring
« PROVISOIRE » (le 62/31 est celui du Kuru chez Dumoulin), le « 1–3 mm »
n'est pas sourcé. ρ = 2640 ± 15 (pas 2620). **À faire (1 h) : interceptes
linéaires par phase sur une photo de carotte ORCHYD** (ImageJ) → remplace
toute la colonne « estimé ».

## 6. Contrôles AVANT calibration (validation Fernando à chaque liste)

Faits (`jobs_controls.json`, 6 runs, 3 × 4 threads) :

| contrôle | résultat | décision |
|---|---|---|
| C1 bruit de réalisation (s1/s2/s3) | pic ± 1,1 %, CI ± 0,007, chute ± 0,025 | un maillage par point ; nugget |
| jeu ×0,7 à 50 / 20 (E 71) | −5 % / −9 % ; CI 0,67 / 0,56 | point de départ correct |
| C6 vitesse ÷ 4 (`q1u070slow`, 0,0625 m/s, ε̇ ≈ 1,3 s⁻¹) | pic 546 vs 567 (**−3,7 %**), ε_pic 0,80 vs 0,83 %, **CI identique** (379), chute 0,76 vs 0,63, rompus 890 vs 667 ; mur × 4 | criblage à 0,25 m/s avec biais +4 % sur le pic à déclarer ; validation à ÷ 4 ; CI est rate-indépendant ; la vitesse rapide LISSE le post-pic (elle ne le crée pas) |
| débit | 4 threads = 320 s, 14 threads = 146 s (A3 corrige le « 77 s » : maillage banni) | **1 × 14** |

Rendus (`jobs_controls2.json`, validés 15:20, `rockim_f2l.exe`, 1 × 14) :

| contrôle | résultat | décision |
|---|---|---|
| base E 77,7 à 50 / 20 | pic **591 (−1,4 %)** / 380 (−6 %) ; CI 0,64 / 0,54 (exp 0,55 / 0,57) ; CD 0,97 / 0,99 ; ε_pic 0,78 / 0,48 % (cibles 0,87 / 0,60) ; chute 0,58 / 0,47 | départ à 1 % du pic à 50 ; CD et ε_pic sont les cibles dures |
| C3 σ₃ = 100 | pic **881 vs 798 (+10 %)**, CI 0,75, **pas de rupture** (6 rompus à 1,05 %) ; pente 50→100 = 5,8 (exp 4,0) | la surestimation à 100 MPa est une erreur de forme (Coulomb linéaire), à déclarer ; prédire 75/100 sans retouche |
| C4 damping 0,3 / 0,1 | pic 568 / 560 (**−4 / −5 %** vs 0,7), CI identique (379), chute **0,63 / 0,77** (0,7 : 0,58), mur ≈ ×1,2 | l'amortissement n'engendre pas la chute (elle CROÎT quand il baisse) ; biais de +4–5 % sur le pic à 0,7, du même ordre que la vitesse ; garder 0,7 au criblage, déclarer |

| C8 pénalité d'insertion 20 | pic 594 (+0,5 %), ε_pic 0,73 % (−7 %), CI 0,64, chute 0,59, mur 426 s | la complaisance des joints insérés vaut ~7 % de ε_pic ; garder 4 (défaut), déclarer |
| C8 schéma intrinsèque | pic **548 (−7 %)**, CD 0,998, chute 0,47, pas de CI mesurable (pas d'insertion), mur 652 s | écart de schéma −7 % sur le pic à déclarer ; schéma adaptatif figé pour toute la calibration |
| C5 40 × 80 mm | **non fait** (Fernando 16:20 : arrêt après l'intrinsèque) — deck `q1v070L_P050.cfg` prêt, sortie interrompue en quarantaine `out_q1v070L_P050_INTERROMPU` | calibration à 20 × 40 mm ; l'effet de taille reste une réserve déclarée |
| cas 3 GBM refait | **abandonné** (Fernando 16:20 : « laisse tomber le GBM ») — deck `q3r_gbm_P050.cfg` conservé | la calibration 4c n'est pas engagée |

Liste initiale (`jobs_controls2.json`, 8 runs, base `q1v070_P050.cfg` : E 77,7,
ν 0,29, ft 8,4, c 17,5, G_f 20, φ 48, G_II/G_I 20, μ_res 0,25, `rockim_f2l.exe`, 1 × 14) :

| deck | variable | coût | décision attendue |
|---|---|---|---|
| `q1v070_P050` / `_P020` | base E physique aux deux confinements | 2 × ~3 min | recalage du départ |
| `q1v070_P100` (C3) | σ₃ = 100 | ~4 min | pente propre du modèle (exp 6,5 → 4,0) ; biais attendu +16 % |
| `q1v070d03` / `d01` (C4) | dampingLocal 0,3 / 0,1 | 2 × ~3 min | part numérique du post-pic ; domaine où la chute est un observable |
| `q1v070ip20` (C8) | insertionPenaltyFactor 20 | ~3 min | complaisance des joints insérés (29 % au pic) |
| `q1v070intr` (C8) | insertion = intrinsic | ~6 min (dt/2) | écart de schéma à déclarer |
| `q1v070L` (C5) | 40 × 80 mm, h 0,8 (13 548 tri., T 4,4 ms) | **30–40 min** | effet de taille (W/ℓ_cz 0,9 → 1,8) ; taille de calibration |

Total ≈ 1 h de mur. Puis (code, ½ journée, `rockim_f2m.exe`) : colonnes
`epsLat`/`epsVol` en grips (σ_cd par inversion de ε_v, opérateur SBM
identique exp/modèle), énergies par ligne, `nBrokTen/Shear` en grips,
`stopPeakDrop` en 2D (le mur passe de 146 à 530 s avec la casse : arrêter
après la chute), avertissement `jointPenaltyFactor` inerte, `weibullScope = lcz`.

## 7. Politique de lancement

Un run à la fois à 14 threads (`calib/runner.py … --parallel 1 --threads 14`),
listes validées par Fernando, binaire figé pour toute une calibration,
extraction automatique. Plan h1 (60 pts LHS, 120 decks, `runs_h1/`) à
**régénérer** sur la base `q1v070` (E 77,7) et l'espace §4a avant lancement.
