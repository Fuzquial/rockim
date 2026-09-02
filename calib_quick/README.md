# calib_quick — calibration rapide Red Bohus, triaxial σ₃ = 50 MPa

**2026-09-02.** Demande : trois simulations de moins de 5 minutes, une seule
variable entre elles, pour choisir la représentation avant de calibrer.
Binaire `rockim_f2j.exe`, 14 threads (`OMP_NUM_THREADS = 14` — sans la
variable rockim se bride à ~3 cœurs).

## Les trois cas

| cas | deck | représentation | ce qui change |
|---|---|---|---|
| 1 | `q1_homog_P050.cfg` | maillage **Gmsh Delaunay** 20 × 40 mm, h = 0,8 mm (2930 triangles), **aucun grain**, joints tous identiques, bulk élastique | — (référence) |
| 2 | `q2_weibull_P050.cfg` | cas 1 + `jointWeibullM = 8`, tirage indépendant par joint (`strengthCorrLength = 0`) | **la dispersion des joints** (ft et cohésion, moyenne 1) |
| 3 | `q3_gbm_P050.cfg` | Voronoï 3 phases (feldspath 62 / quartz 31 / biotite 7), grains 3 mm (113 grains), delaunay intra-grain 0,54 mm (5814 triangles) | **le contraste de phase**, aux moyennes pondérées égales au cas 1 (E 80, ft 12, c 25) |

Commun aux trois : chargement de la sonde 4 du 31 août (pullV −0,25 m/s,
rampe 0,2 ms, délai 0,3 ms, confinement 50 MPa rampé, mors libres
latéralement, dampingLocal 0,7, T = 2,2 ms → 1,19 % de déformation) ;
**bulk élastique** (`law = mc` est incompatible avec `phases`, donc le même
bulk partout, `crushCap = 1e12`) ; joints de la sonde 4, **non calibrés** :
ft 12 MPa, c 25, φ 48°, G_f 40, G_II/G_I 20, μ_res 0,25, pénalité 20,
adoucissement Yan, insertion adaptative.

## Choix à connaître

- **Pas de Weibull de volume en 2D** (`matWeibullM` n'existe qu'en fem3d) :
  le cas 2 ne disperse que les joints.
- **GBM à joints de grain non affaiblis** (`gbAlpha*` omis = 1,0,
  `gbHeteroFactor` omis = 1,0) : le cas 3 n'ajoute que le contraste de phase.
  L'affaiblissement des joints de grain est le levier suivant, pas une
  variable cachée de la comparaison.
- **Équivalence des moyennes** : feldspath ajusté à E = 84,2 GPa,
  ft = 11,73 MPa, c = 19,44 MPa pour que 0,62·feldspath + 0,31·quartz
  (83,1 / 14 / 35) + 0,07·biotite (29,3 / 5,5 / 30) redonne exactement
  80 / 12 / 25.
- **Grains de 3 mm et non 2** : pour tenir sous 5 min. 113 grains sur une
  éprouvette de 20 × 40 mm — la réserve 2 du §3ter de la campagne d'août
  (« c'est le NOMBRE de grains ») s'applique en plein ; c'est un premier tri,
  pas un résultat.
- **ε̇ ≈ 5 s⁻¹** (réserve 1 du §3ter) conservé : le post-pic est lissé par
  construction. À traiter à la reprise de la calibration (pullV ÷ 10).
- **Maillage : Delaunay pur (`Mesh.Algorithm = 5`) + champ de taille bruité ±10 %**,
  après un premier essai en frontal (algo 6) que Fernando a fait bannir à la
  vue du pavage quasi équilatéral — l'étude `mesh_algo_sweep.py` d'août
  l'avait mesuré, je ne l'avais pas appliqué. Contrôle chiffré sur
  20 × 40 mm, h = 0,8 mm (orientations d'arêtes en 36 classes ; R6 = ordre
  orientationnel 6-fold, 0 = isotrope, 1 = nid d'abeille) :

  | algo | triangles | pic/creux | angle min médian | R6 |
  |---|---|---|---|---|
  | 6 frontal (**banni**) | 3024 | **17,1** | 57,2° | **0,341** |
  | **5 Delaunay + bruit** | 3388 | **2,2** | 48,6° | **0,035** |

  Les sorties faites sur le maillage banni sont conservées sous
  `out_*_algo6_BANNI` (cas 1-2 : mêmes verdicts, aucune rupture) et ne
  doivent pas être citées. Le cas 3 (Voronoï + delaunay intra-grain) n'est
  pas concerné.

## Coût mesuré

| cas | éléments | dt | pas | mur |
|---|---|---|---|---|
| 1 | 2930 | 5,13 ns | 429 000 | **77 s** (14 threads) |
| 2 | 2930 | 5,13 ns | 429 000 | **80 s** |
| 3 | 5814 | 3,11 ns | 707 000 | |

### Réalisations du maillage (bruit de réalisation)

`tools/make_unstructured_mesh.py box2d … [seed]` promet de varier la
réalisation par `Mesh.RandomSeed` : **inopérant en 2D Delaunay** (seeds 1, 2, 3
→ trois fichiers strictement identiques, md5 `4628b4e3ab`). Le levier qui
marche est le **champ de taille bruité à phases seedées** :
`calib_quick/make_box_mesh.py --seeds 1 2 3` → `meshes/box20x40_h08_s{1,2,3}.msh`
(3413 / 3394 / 3384 triangles, pic/creux 1,9 / 1,8 / 1,5, R6 0,065 / 0,028 /
0,018, angle min médian 48–49°) — trois maillages isotropes et différents.
Decks témoins `q1_s{1,2,3}_P050.cfg` (cas 1 inchangé, maillage seul varie) :
c'est la mesure du bruit de l'objectif, à faire AVANT toute calibration.

**Mesure (2026-09-02, 4 réalisations : algo5 + s1/s2/s3, `rockim_f2l.exe`) :**

| observable | valeurs | écart-type | tolérance de calibration |
|---|---|---|---|
| pic q [MPa] | 709,4 / 717,8 / 704,2 / 699,1 | **8 (1,1 %)** | 3 % |
| ε_pic [%] | 0,871 / 0,885 / 0,852 / 0,852 | 0,016 (1,9 %) | 10 % |
| chute | 0,52 / 0,53 / 0,58 / 0,55 | 0,025 | 0,10 |
| CI/pic | 0,637 / 0,629 / 0,641 / 0,645 | 0,007 | 0,10 |
| CD/pic | 0,994 / 0,977 / 0,993 / 0,996 | 0,009 | 0,10 |
| joints rompus | 430 / 387 / 388 / 410 | 21 | — |

Le bruit de réalisation est **3 à 10 fois sous les tolérances** : un seul
maillage par point du plan suffit pour l'homogène (le nugget du GP l'absorbe) ;
les seeds ne sont nécessaires que pour Weibull et GBM (tirages internes).
Coût : 318–343 s à 4 threads contre 77 s à 14 threads pour le même run —
le scaling OpenMP est quasi linéaire à 3 000 éléments, donc **un run à
14 threads à la fois** rend plus que trois à 4 threads (110 s équivalents).

### ⚠️ Maillage intra-grain structuré (remarque de Fernando, 14:00)

Le Delaunay intra-grain du GBM place ses points intérieurs sur un réseau
triangulaire : sur le cas 3, **R6 = 0,548 et pic/creux 18,8 sur les arêtes
intra-grain** — plus structuré encore que le frontal banni. Corrigé en ajout
(`grainMeshRandom = true`, Poisson-disc, `rockim_f2n.exe`) : R6 0,007,
pic/creux 1,37, +15 % d'éléments, dt −15 %. Planche `fig_tess_grainmesh.png`.
**Le cas 3 ci-dessous (515 → 481 corrigé) est donc à refaire**
(`q3r_gbm_P050.cfg`) ; tous les decks GBM de calibration posent la clé.

### Contrôles C3–C8 (2026-09-02 15:20–16:25, `rockim_f2l.exe`, 1 × 14 threads)

Base `q1v070_P050.cfg` (×0,7 des joints de la sonde 4, E 77,7, ν 0,29,
ℓ_cz 22 mm), déviateur corrigé. Figure `fig_controls.png`
(`calib/plot_controls.py`). Tableau et décisions : `PLAN_calibration.md` §6.
Résumé : pic de base **591 à 50 (−1 %)**, 380 à 20 (−6 %), 881 à 100
(+10 %, Coulomb linéaire) ; bruit de maillage négligeable ; vitesse ÷ 4 →
−4 % ; amortissement 0,3/0,1 → −4/−5 % et chute PLUS forte ; pénalité
d'insertion 20 → pic +0,5 %, ε_pic −7 % ; schéma intrinsèque → pic −7 %.
**CD collé au pic dans tous les cas** (0,95–1,00 contre 0,72). Lot arrêté
par Fernando après l'intrinsèque : 40 × 80 mm non fait, GBM abandonné.

### Contrôle de vitesse (C6, 2026-09-02 15:10)

`q1u070slow_P050` = `q1u070_P050` à pullV −0,0625 (÷ 4, ε̇ ≈ 1,3 s⁻¹), T × 4 :
pic **546 vs 567 MPa (−3,7 %)**, ε_pic 0,80 vs 0,83 %, **CI identique** (379 MPa,
0,69 du pic), chute 0,76 vs 0,63 (890 vs 667 joints rompus), mur × 4.
Le pic porte un biais inertiel de ≈ +4 % à 0,25 m/s (≈ 1 tolérance) ; CI n'en
porte aucun ; la vitesse rapide lisse le post-pic. Décision : criblage à
0,25 m/s (biais déclaré), validation finale à ÷ 4.

## ⚠️ CORRECTION C0 (2026-09-02 14:30, critique adverse A5) — lire avant les tableaux

Les tableaux ci-dessous (« Résultats », « Balayage ») ont été calculés avec
q = sigma − (moyenne de sigma avant `pullDelay`). Or pendant `pullDelay` les
mors sont **bloqués** (ε_yy = 0) sous σ_xx = −σ₃, donc σ_yy = ν/(1−ν)·σ₃ =
16,7 MPa à 50 MPa (mesuré exactement) et non σ₃ : l'état de départ n'est pas
isotrope. Le déviateur comparable à l'essai est **q = sigma − σ₃**, avec
ε = 0 quand sigma atteint σ₃ (exact pour un bulk élastique). Les anciens q
sont donc **surestimés de 33 MPa à 50 MPa et de 13 MPa à 20 MPa**. Valeurs
corrigées (`calib/extract.py`, `plot_quick.py` corrigés ; les anciens
chiffres restent dans les tableaux comme trace) :

| run | q pic ancien | **q pic corrigé** | écart à 599 | CI/pic | CD/pic |
|---|---|---|---|---|---|
| 1 homogène | 709 | **676** | +13 % | 0,62 | 0,99 |
| 2 Weibull m = 8 | 709 | **676** | +13 % | **0,53** | 0,99 |
| 3 GBM α = 1 | 515 | **481** | −20 % | 0,30 | 0,95 |
| ×0,5 (ℓ_cz cst) | 548 | **514** | −14 % | 0,69 | 0,95 |
| ×0,7, E 71 (`q1u070`) à 50 | 596 | **567** | −5 % | 0,67 | 0,98 |
| ×0,7, E 71 à 20 (cible 405) | 381 | **370** | −9 % | 0,56 | 0,98 |

L'exposant de sensibilité du pic aux résistances est **0,40** (676 → 514
pour × 0,5), pas 0,7 comme écrit plus bas : le pic est 1,45–1,6 × l'amorçage
Coulomb de l'arête la plus critique (blocage cinématique des triangles), une
propriété du maillage autant que des joints (A5, constat C). Le durable :
une consolidation isotrope (pression axiale = σ₃ pendant `pullDelay`, clé à
créer en ajout) ; en attendant la correction de dépouillement est exacte.

Second changement de convention (A5, R1) : **E reste physique (77,7 GPa,
ν = 0,29)** dans les decks — nécessaire au transfert vers le forage (vitesses
d'onde, raideurs de contact) — et les **cibles en déformation sont
multipliées par (1−ν²) = 0,916** (déformation plane, σ_xx constant pendant le
déviateur : Δε_yy = Δσ_yy(1−ν²)/E exactement). Base des decks suivants :
`q1v070_P050.cfg` (E 77,7, ν 0,29, ft 8,4, c 17,5, G_f 20 → ℓ_cz 22 mm).

## Résultats sur le maillage isotrope (algo 5, `rockim_f2j.exe`)

| cas | pic q [MPa] | ε_pic | E [GPa] | chute | joints rompus | écart au pic exp (599) |
|---|---|---|---|---|---|---|
| 1 homogène | 709,4 | 0,87 % | 85,3 | 52 % | 430 | **+18 %** |
| 2 Weibull m = 8 | 709,0 | 0,87 % | 85,3 | 53 % | 431 | +18 % |
| 3 GBM α = 1 | 514,7 | 0,65 % | 83,7 | 58 % | 315 | **−14 %** |

Lecture :

- **Le maillage banni cachait la rupture** : sur le pavage quasi équilatéral
  les cas 1-2 ne cassaient pas (q 959, 0 rompu) ; sur le maillage isotrope le
  cas 1 casse à 709 MPa avec 430 joints rompus. Le chemin de fissure a besoin
  d'orientations d'arêtes variées.
- **Weibull m = 8 n'a aucun effet** (709,0 vs 709,4) : à ce niveau de
  dispersion la moyenne des joints pilote le pic. La dispersion ne devient un
  levier que via m ET les seuils σ_ci / σ_cd (initiation), pas via le pic.
- **Le contraste de phase seul fait −14 %** : la biotite (E 29 GPa) concentre
  les contraintes aux joints de grain et déclenche la rupture à 0,65 % ; les
  trois représentations encadrent la cible (+18 / −14 %).
- **E = 85 GPa n'est pas une erreur** : déformation plane, E/(1−ν²) = 80/0,9375
  = 85,3. Pour viser 77 apparent il faut E ≈ 72 GPa dans les decks.

### Balayage des résistances de joint à G_f fixe — ARTEFACT, à ne pas citer

| deck | ft, c | G_f | pic q | rompus | chute |
|---|---|---|---|---|---|
| `q1s050` | × 0,5 | 40 (fixe) | 752,7 (+26 %) | 84 | 24 % |
| `q1s035` | × 0,35 | 40 | 773,6 (+29 %) | 6 | 0 % |
| `q1s025` | × 0,25 | 40 | 773,5 (+29 %) | 3 | 0 % |

Des joints **plus faibles** donnent un pic **plus haut** : à G_f fixe,
ℓ_cz = E G_f / ft² est multiplié par 4, 8, 16 — les joints s'ouvrent tôt mais
adoucissent sur une longueur énorme et ne cassent plus (3 rompus), le bulk
élastique porte tout. Le balayage correct tient ℓ_cz constant : **G_f ∝ ft²**
(`q1t050/035/025_P050.cfg`, G_f = 10 / 4,9 / 2,5). Règle pour la calibration :
ft, c et G_f ne sont pas indépendants — paramétrer (ft, c, ℓ_cz) plutôt que
(ft, c, G_f).

### Balayage corrigé à ℓ_cz constant (G_f ∝ ft²) — résultat valide

| deck | ft, c | G_f | pic q | ε_pic | chute | rompus | mur |
|---|---|---|---|---|---|---|---|
| cas 1 | × 1 | 40 | 709,4 (+18 %) | 0,87 % | 52 % | 430 | 77 s |
| `q1t050` | × 0,5 | 10 | 547,5 (−9 %) | 0,72 % | 64 % | 925 | 235 s |
| `q1t035` | × 0,35 | 4,9 | 496,1 (−17 %) | 0,59 % | 66 % | 1274 | 367 s |
| `q1t025` | × 0,25 | 2,5 | 439,3 (−27 %) | 0,53 % | 62 % | 1738 | 530 s |

Monotone et sain : le pic suit les résistances (exposant ≈ 0,7 sur le
facteur), E inchangé (85,3), la chute ~ 60 % (bulk élastique + confinement :
aucun mécanisme de résistance résiduelle progressive). La cible 599 MPa
tombe entre × 1 et × 0,5 → **facteur ≈ 0,7 (ft ≈ 8,4, c ≈ 17,5, G_f ≈ 20)**
pour le cas homogène — mais ε_pic (0,72–0,87 % vs 0,95 %) et la chute
(52–64 % vs 26 %) montrent que le pic seul ne suffit pas : la non-linéarité
pré-pic (σ_ci → σ_cd) et la résistance résiduelle appellent un bulk
non élastique ou une résistance résiduelle de joint plus haute (μ_res, φ).
Le coût croît avec le nombre de joints rompus (77 → 530 s) : le budget de
5 min n'est tenu que pour les jeux qui cassent peu.

### Seuils CI / CD : les observables qui départagent les représentations

`history.csv` porte `nInserted` et `nDamaging` en insertion adaptative : le
**premier joint inséré** naît sur l'enveloppe de Mohr-Coulomb — c'est
l'amorçage de la microfissuration, le proxy de σ_ci (SBM : 0,55 du pic sur les
12 essais). La **première rupture complète** (D ≥ 1, `nBroken`) est le proxy
de σ_cd (exp 0,72). `calib/extract.py` les sort ; `calib/joint_frames.py`
donne l'histoire par frame depuis les VTU de joints.

| run | pic q | CI_frac (exp 0,55) | CD_frac (exp 0,72) | chute (exp 0,26) |
|---|---|---|---|---|
| 1 homogène | 709 | 0,64 | 0,99 | 0,52 |
| 2 Weibull m = 8 | 709 | **0,55** | 0,99 | 0,53 |
| 3 GBM α = 1 | 515 | 0,35 | 0,96 | 0,58 |
| ×0,5 (ℓ_cz cst) | 548 | 0,71 | 0,95 | 0,64 |
| ×0,35 | 496 | 0,74 | 0,92 | 0,66 |
| ×0,25 | 439 | 0,81 | 0,97 | 0,62 |

Lecture : le pic seul est dégénéré (Weibull = homogène à 0,1 %), mais **CI
sépare les trois** : la dispersion de Weibull abaisse CI sans toucher le pic
(0,64 → 0,55, exactement l'expérience), le contraste GBM amorce trop tôt
(0,35 avec α = 1 : il faudra des joints intra-granulaires plus forts que la
moyenne et des joints de grain plus faibles). **CD est collé au pic** dans
tous les cas (0,92–0,99 vs 0,72) : la croissance stable entre CD et le pic
n'existe pas dans le modèle actuel — c'est la cible dure de la calibration
(ℓ_cz, G_II/G_I, frottement résiduel).

## Polydispersité des grains (`rockim_f2l.exe`, 2026-09-02)

Demande de Fernando : « dire au Voronoï de faire varier les tailles tout en
gardant les proportions globales de chaque phase ». Clés `grainSizeSpread`
(sd de ln taille) et `phase.<nom>.grainSize` (taille cible par phase,
affinité) — méthode et validation dans `DOCUMENTATION_rockim.md` §5.16
(Laguerre à aires prescrites, Newton amorti ; deux impasses documentées).
Planche `fig_tess_polydisp.png` (`plot_tess.py`) : monodisperse (sd 0,07),
0,5 → 0,50 réalisé, 0,5 + tailles par phase (le quartz prend les gros grains,
la biotite les petits), 0,8 → 0,80. Decks d'essai `_poly_*.cfg`, sorties
`_poly_*_lag/`, bit-identité `_bit_j/` vs `_bit_l/`.

Pour la calibration GBM : `grainSize` = taille moyenne, `grainSizeSpread` ∈
{0 ; 0,3 ; 0,6} et `phase.quartz.grainSize` / `phase.biotite.grainSize` selon
la lame mince (Bohus : quartz 1–4 mm, feldspath 2–6 mm, biotite < 1,5 mm —
à confirmer sur la pétrographie de la campagne triaxiale). ⚠️ hmin suit le
plus petit grain : dt ÷ 3,6 à σ = 0,8.

## Dépouillement

```
python calib_quick/plot_quick.py
```

q = `sigma` de `history.csv` moins l'offset de consolidation (moyenne juste
avant `pullDelay`) ; ε = déplacement imposé des mors (rampe cosinus
analytique, `simcurve.py` de la campagne d'août) / H. Cibles :
`rockim_f1/calib_triax3d/targets_triax_bohus.json`, confinement 50 :
**pic 599,2 MPa à 0,95 %, chute 26 %, E ≈ 77 GPa**. Sortie :
`fig_quick_P050.png` + tableau pic / ε_pic / E / chute / joints rompus.
