# Documentation rockim — commandes, simulations, entrées, sorties

*Référence établie le 2026-08-11 par extraction exhaustive depuis le code source
(archive `rockim_yan3d_2026-08-11`). Chaque clé, défaut et sortie listés ici a été
vérifié dans le source — le code fait foi. Complète le `GUIDE_rockim.md` (prise en
main) et le `README.md` (théorie des modèles, EN) sans les remplacer.*

---

## 1. Vue d'ensemble

rockim est un exécutable unique (`rockim` / `rockim.exe`) qui lit un fichier de
configuration texte et écrit un dossier de résultats. Six solveurs derrière la même
interface, choisis par la clé `mode` :

| mode | modèle | fissuration | usage type |
|---|---|---|---|
| `fem` | FEM 2D déformation plane, CST | endommagement lissé + érosion | percussion/coupe rapides, onde de barre |
| `fem3d` | FEM 3D tets de Kuhn, co-rotationnel, **5 lois au choix** | endommagement + érosion | percussion 3D multi-lois (contrepartie Abaqus/VUMAT) |
| `dem` | BPM 2D (liaisons parallèles) | rupture de liaisons | comparaison discrète historique |
| `dem3d` | BPM 3D HCP | rupture de liaisons | idem 3D |
| `fdem` | **FDEM 2D Munjiza** + GBM | joints cohésifs explicites | l'outil de laboratoire complet (BD, UCS, triaxial, SHPB) |
| `fdem3d` | **FDEM 3D** tets + GBM 3D | joints cohésifs triangulaires | percussion/coupe 3D à fissures explicites |

Unités : **SI partout** (m, s, Pa, kg). Point décimal obligatoire (`0,5` est rejeté
avec un message nommant la clé). En 2D les forces sont par mètre d'épaisseur
(`thickness`, défaut 1 m).

## 2. Compiler

**Linux / macOS (CMake)** — Eigen est trouvé (paquet système) ou téléchargé :

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j            # produit ./rockim ; OpenMP détecté automatiquement
```

**Windows / MSVC** (depuis le dossier `rockim/` extrait, Eigen à côté) :

```bat
"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX ^
   /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim.exe
```

Seule dépendance : Eigen (headers seuls). OpenMP est optionnel (le code compile et
tourne en sériel sans).

## 3. Lancer

### 3.1 Ligne de commande

```bash
./rockim <config.cfg> [dossier_sortie]     # dossier par défaut : clé outputDir, sinon "out"
./rockim selftest-saksala2011 [out.csv]    # rejoue le harnais VUMAT (réf. Fortran, 8e-14)
./rockim selftest-dpdfh       [out.csv]    # idem DP-DFH (4.7e-12, tirages bit-identiques)
```

Le run affiche : bannière d'init (éléments, joints, dt, nombre de pas), progression
en %, puis le **résumé** (pics, casse, bilans d'énergie, verdicts PASS/FAIL).
Code retour 0 si le run va au bout, 1 sur erreur de config ou instabilité (NaN).

### 3.2 Variables d'environnement

| variable | effet |
|---|---|
| `OMP_NUM_THREADS` | nombre de threads. **1 = bit-identique au build sériel** ; N fixé = déterministe ; N différents = écarts d'associativité (±2 % sur les pics — comparer à threads égaux) |
| `ROCKIM_PROF=1` | profil par pas en fin de run (modes fdem et fdem3d) : elem / insert / joint / gcontact / tool en ms/pas |
| `RKM_NOTAU=1` | coupe la traction tangentielle des joints (interrupteur de bissection, fdem 2D) |
| `RKM_NOGC=1` | coupe le contact général (bissection, fdem 2D) |
| `ROCKIM_BRAZ_DEBUG=1` | profil du chemin de charge du brésilien au résumé |

### 3.3 Suite de non-régression

```bash
python3 tools/verify_suite.py --exe build/rockim              # tier fast (~3 min)
python3 tools/verify_suite.py --exe build/rockim --tier full  # + bit-repères longs, 3D (~30 min)
python3 tools/verify_suite.py --exe build/rockim --tier all   # + Voronoï 3D (~1 h)
python3 tools/verify_suite.py --exe build/rockim --only fdem3d --json rapport.json
```

Références en dur (baseline Linux 2026-08-11 — re-baseliner une fois sous MSVC),
tolérances par nature de test, contrôles à charge nulle et dampWork ≤ 0 inclus.
**À lancer après toute modification du code.**

### 3.4 Interface graphique

```bash
python tools/rockim_gui.py     # tkinter+matplotlib : édition de configs, lancement,
                               # suite de vérification en un clic, tracés intégrés
```

## 4. Le fichier de configuration

Format : `clé = valeur`, une par ligne ; `#` ouvre un commentaire ; clé répétée =
la dernière gagne (pratique pour surcharger une config de base par ajout en fin de
fichier). Les valeurs numériques sont parsées STRICTEMENT (virgule décimale ou
suffixe parasite → erreur nommant la clé). ⚠️ Les clés inconnues sont **ignorées
silencieusement** — relire l'orthographe en cas de comportement par défaut inattendu.
Les combinaisons incohérentes (ex. `phases` sans `mesh = voronoi`, `scenario =
brazilian` sans `geometry = disc`) arrêtent le run avec un message explicite.

## 5. Référence des clés

Colonne « portée » : modes qui lisent la clé. Défauts entre parenthèses.

### 5.1 Bloc commun

| clé (défaut) | rôle | portée |
|---|---|---|
| `mode` (fem) | fem \| fem3d \| dem \| dem3d \| fdem \| fdem3d | — |
| `scenario` (percussion) | percussion \| shear \| tension ; + `bar_wave` (fem) ; + `brazilian`, `shpb` (fdem) | tous |
| `geometry` (box ; disc si brazilian, shpb si shpb) | box \| disc (fdem) ; box \| cylinder (fem3d) | fdem, fem3d |
| `mesh` (grid) | grid \| voronoi (GBM) \| **file** (maillage non structuré importé, « à la Yan ») | fdem, fdem3d |
| `meshFile` (requis si mesh = file) | chemin d'un Gmsh MSH 2.2 ASCII (type 2 en 2D, type 4 en 3D) ; boîte translatée à l'origine, W/H/D relus de l'enveloppe ; générer via `tools/make_unstructured_mesh.py` | fdem, fdem3d |
| `T` (2.5e-4 fdem ; 2e-4 3D ; 2e-4 fem ; 2.5e-4 dem) | durée physique [s] | tous |
| `frames` (50) | nombre de frames VTU écrites | tous |
| `outputDir` (out) | dossier de sortie si absent de la CLI | tous |
| `W, H` (0.2×0.2 fdem ; 0.2×0.1 fem/dem) + `D` (3D) | dimensions du bloc [m] | tous |
| `thickness` (1.0) | épaisseur 2D [m] | fem, dem, fdem |
| `nx, ny` (64×64 fdem ; 96×48 fem) + `nz` (20×20×15 fdem3d ; 24×24×18 fem3d) | découpage grille | maillages grid |
| `seed` (12345 ; **42 en dem**) | graine du maillage (jitter, Voronoï, phases) | tous |
| `dtFactor` (0.2 fdem/dem ; 0.15 fdem3d ; 0.3 fem3d) / `cfl` (0.7, fem 2D) | fraction du pas critique | tous |
| `dampingLocal` (0.02 fdem dyn ; 0.7 fdem QS ; **0 en SHPB** ; 0.05 fdem3d/fem3d ; 0.02/0.7 dem) | amortissement de Cundall | tous sauf fem 2D |
| `gravity` (0) | force de volume ρg selon −y [m/s²], valeur positive | fdem |
| `extraContacts` (2 fdem ; 8 dem) | budget de contacts dans le dt stable | fdem*, dem* |

### 5.2 Matériau (partagé) et phases minérales (GBM)

Bloc global : `rho` (2650), `E` (50e9), `nu` (0.25), `ft` (10e6), `cohesion` (25e6),
`frictionDeg` (40), `Gf` (70), `gfShearFactor` (10 ; Gf_II = facteur × Gf_I).
Validation stricte : E, rho, ft, cohesion, Gf > 0 ; nu ∈ [0, 0.5) ; frictionDeg < 89.

Phases (fdem/fdem3d + `mesh = voronoi`) :

```
phases = quartz feldspar biotite          # déclare les noms
phase.quartz.fraction = 0.33              # fraction surfacique/volumique OBLIGATOIRE
phase.quartz.E = 94e9                     # toute propriété du bloc matériau est
phase.quartz.ft = 13e6                    # surchargeable par phase
```

Joints de grains = moyenne des deux phases × facteurs d'atténuation :
`gbAlphaTen`, `gbAlphaCoh`, `gbAlphaGf`, `gbAlphaE`, `gbAlphaFric` (tous 1.0) ;
frontières hétérophases : × `gbHeteroFactor` (1.0) en plus sur les résistances.
⚠️ un ft de joint nul rendrait le joint incassable — modéliser une frontière
pré-fissurée par un petit α (1e-3), jamais 0.

### 5.3 Maillage Voronoï / GBM (`mesh = voronoi`, fdem et fdem3d)

| clé (défaut) | rôle |
|---|---|
| `grainSize` (**requis**) | diamètre moyen de grain [m] |
| `grainSeeding` (hex) | hex (compact, ANISOTROPE) \| random (Poisson — recommandé dès que les trajets de fissures comptent) |
| `grainJitter` (0.5) | jitter du semis hex (fraction du pas) |
| `lloydIters` (2) | itérations de relaxation de Lloyd |
| `vertexMergeFrac` (0.12) | tolérance de contraction d'arêtes courtes × grainSize — **0.25 recommandé en 3D** (les éventails des faces minces contrôlent le dt) |
| `refineLevels` (0) | raffinement conforme intra-grain : 0..4 en 2D (×4/niveau), 0..2 en 3D (×8/niveau) |
| `grainMesh` (fan) | fan (éventail depuis le centroïde) \| delaunay (maillage non structuré intra-grain — la pratique Y-Geo/Irazu) — 2D |
| `grainElemSize` (0 = 0.18·grainSize) | taille cible du Delaunay intra-grain [m] ; **requis** pour `discMesh = native` |
| `meshJitter` (0 grid ; 0.25 disque natif) | désordre des nœuds (fraction de maille) |
| `meshMirror` (true) | fem3d : Kuhn miroité en damier (8 familles de diagonales) ; false = ancien maillage à l'identique |
| `discMesh` (cut) | fdem+disc : cut (découpe du box, jante en escalier) \| native (anneau exactement sur le cercle + remplissage hex + Delaunay) |

### 5.4 Joints cohésifs (fdem / fdem3d)

| clé (défaut) | rôle |
|---|---|
| `jointPenaltyFactor` (20) | pénalité intrinsèque p = facteur·E/h — complaisance ~4-5 % sur E, dt ∝ 1/√facteur |
| `jointXi` (0.05) | ratio d'amortissement du dashpot de joint (bilatéral sur joint intact, résultante écrêtée sur joint rompu, borne cd ≤ m/dt). **Règle maison : 0 pour les vérifications de loi, 0.01 en quasi-statique, 0.05 en impact** |
| `jointSoftening` (linear) | linear \| **yan** = f(D) exponentielle de Yan et al. 2023 (éq. 11), aires sous les branches = exactement GfI/GfII |
| `yanA`/`yanB`/`yanC` (0.63/1.8/6.0) | constantes de f(D) |
| `yanQuadN` (4096) | points de Simpson pour ∫f(D)dD (= 0.386307 aux défauts) |
| `jointFrictionScaled` (0) | 1 = le terme de Coulomb est aussi multiplié par f(D) (éq. 10 littérale — un joint broyé perd alors tout frottement résiduel) |
| **`insertion`** (intrinsic) | intrinsic \| **adaptive** = insertion dynamique extrinsèque (Yan et al. 2023) : aucun joint à t = 0, liaison cinématique exacte, activation quand σ_n ≥ ft ou \|τ\| ≥ c − σ_n·tanφ, continuité de contrainte à l'insertion. Gains mesurés : dt ×2, mur ×2.2–2.7, complaisance nulle |
| `insertionPenaltyFactor` (4) | pénalité des joints ACTIVÉS en mode adaptatif (décharge/contact) |
| `jointWeibullM` (0 = off) | m > 1 : ft et cohésion de chaque joint × facteur Weibull(m) de moyenne 1 (Gf non tiré) |
| `strengthCorrLength` (0) | 0 = tirages indépendants ; > 0 = champ gaussien corrélé (copule) de cette longueur [m] |
| `strengthCorrLengthB` (= A) + `strengthCorrAngleDeg` (0) | anisotropie orientée (foliation ; en 3D plan de texture incliné autour de y) |
| `fieldSeed` (seed+777) | graine du CHAMP, indépendante du maillage — deux maillages voient les mêmes zones faibles |
| `crushCap` (8·cohesion) | plafond élasto-plastique du déviateur du bulk (garde-fou, désactivé si une `law` est active) |

### 5.5 Lois de comportement (`law`, modes fem3d / fdem / fdem3d)

`law = elastic | dpr | saksala | saksala2011 | dpdfh` (défaut : dpr en fem3d ; absent
en fdem/fdem3d = bulk élastique + crushCap, bit-compatible avec l'historique).
En fdem 2D la loi 3D est utilisée en déformation plane exacte (ε_zz = 0). `law` est
incompatible avec `phases` (mono-matériau).

| loi | clés propres (défauts) |
|---|---|
| `dpr` | Drucker-Prager calé MC + Rankine crack-band : `erodeD` (0.98), `erodeEpv` (1.5) |
| `saksala` | + Perzyna et cap : `saksalaEta` (0.05e6 Pa·s), `capP0` (8·cohesion), `capH` (K) |
| `saksala2011` | portage VUMAT fidèle (vérifié 8e-14) : `skBetaDP` (0.0346), `skCres` (2.89e6), `skHdp` (−10e9), `skSdp`/`skSmr` (1e4), `skAt` (0.98), `skBetaT` (5000), `skPp0` (1040e6), `skPtr0` (377e6), `skDcap` (1e-9), `skWcap` (0.0433), `skNd` (7.5e-8) — **défauts = Table I du papier** ; poser E=60e9, nu=0.2, ft=13e6, cohesion=37.5e6, frictionDeg=30, rho=2600 |
| `dpdfh` | portage DP-DFH de la thèse (vérifié 4.7e-12) : `dfhBetaDeg` (51.7), `dfhDCoh` (153.3e6), `dfhPsiDeg` (15), `dfhWeibullM` (24), `dfhSigW` (120e6), `dfhZeff` (1e-9), `dfhK` (0.38), `dfhS` (4.18879), `dfhDeld` (1e9 = suppression OFF) — **défauts = carte Red Bohus** ; poser seulement E=52e9, nu=0.25, rho=2620 (ft/cohesion/frictionDeg ignorés par cette loi) |

Hétérogénéité des lois : `matWeibullM` (0 = off, fem3d) tire un facteur de résistance
par élément (i.i.d. ou champ corrélé via les mêmes clés strengthCorr*/fieldSeed).
⚠️ consommé par saksala2011 et dpdfh seulement ; et les tirages spatiaux de dpdfh
(hash du centroïde) ne sont renseignés qu'en fem3d à ce jour.

### 5.6 Outil et contact

| clé (défaut) | rôle | portée |
|---|---|---|
| `toolShape` | fem/dem : disc \| flat ; fdem : disc \| flat \| pdc (coupe) ; 3D : sphere \| flat (poinçon cylindrique) | percussion/shear |
| `toolRadius` (0.015 ; 0.01 fem) | rayon [m] | tous |
| `toolMass` (5.0 en 2D ; 0.5 en 3D) | masse [kg] (percussion = outil LIBRE) | tous |
| `toolWidth` (0.02) | largeur du flat 2D [m] | fem, dem |
| `impactSpeed` (8 ; 15 fem/dem 2D) | vitesse d'impact [m/s] | percussion |
| `toolGap` (1e-4) | jeu initial outil-surface [m] | tous |
| `toolX` (W/2 percussion ; −R−gap shear) + `toolY` (D/2, 3D) | position initiale | tous |
| `cutDepth` (0.004 ; 0.003 fem) / `cutSpeed` (10) | profondeur et vitesse de coupe | shear |
| `backRakeDeg` (20), `cutterLen` (0.013), `chamferLen` (0), `chamferDeg` (45) | couteau PDC 2D | fdem shear |
| `contactMu` (0.5 ; 0.3 fem) | frottement outil/platines | tous |
| `contactXi` (0.05 ; 0.1 dem) | amortissement du contact outil (fraction du critique) | tous |
| `contactVreg` (1e-3) | vitesse de régularisation du frottement tanh [m/s] | fem3d, fdem* |
| `kpFactor` (1.0) | pénalité outil = facteur·E·t | fem |

Contact général (débris, fdem/fdem3d) : `gcPenaltyFactor` (0.01 ×E·t — mou exprès),
`gcXi` (0.8), `gcRestitution` (0.2, quasi-plastique), `gcBirthTau` (1e-6 s, relaxation
de la pénétration de naissance) ; SHPB : `gcCell` (2·hDisc), `gcBoxMesh` (true),
`gcXwindow` (0.10 m autour du disque — le contact est le SEUL chemin de charge).

### 5.7 Essais quasi-statiques (fdem 2D sauf mention)

**Traction / compression (`scenario = tension`)** — `pullV` (0.05 ; **< 0 =
compression**), `pullRamp` (0 ; rampe cosinus [s] — sans rampe le transitoire casse
au mors), `pullDelay` (0 ; retarde l'axial, pour équilibrer un confinement),
`gripLateralFree` (false ; mors sans frottement latéral), `loading` (grips \|
platens — UCS/triaxial par PLATINES frottantes, la pratique Y-Geo), `verifyFt`
(true ; false = pas de PASS/FAIL contre ft). Aussi en fem3d/dem3d/fdem3d (grips).
Métrologie UCS (platens) : `gaugeLoFrac`/`gaugeHiFrac` (0.25/0.75, extensomètre
intérieur), `ucsStopAfterPeak` (false) + `ucsStopDelay` (5e-5 s).

**Brésilien (`scenario = brazilian`, exige `geometry = disc`)** — `discFlattenDeg`
(0 ; angle TOTAL 2α du disque aplati, correction k de Wang appliquée),
`brazilianLoading` (platens \| traction), `platenHalfWidth` (R, ou R·sinα si aplati),
`platenPenaltyFactor` (1.0 ×E·t), `platenTributary` (true ; poids par longueur de
jante tributaire), `pullV` = taux de FERMETURE total des deux platines,
`loadArcDeg` (7.5) + `loadRate` (20·ft/T) pour le mode traction,
`elasticGaugeLo`/`Hi` (0.3/0.8 ×ft ; bande de la jauge élastique du centre),
`brazilianStopAfterPeak` (false) + `brazilianStopDelay` (5e-5), `diametralBand`
(0.15 ×R ; verdict de diamétralité).

**Confinement (fdem, fdem3d)** — `confiningPressure` (0 = off, > 0 [Pa]),
`confiningRamp` (0 en 2D — **toujours en mettre une** ; 2e-4 en 3D), `confineFaces`
(sides \| all en 2D ; lateral \| all en 3D), `confineGaugeTime` (3×rampe ; instant de
la jauge σ_latéral atteint). Pression SUIVEUSE sur les faces extérieures d'origine
seulement (pas dans les fissures). tension + pullV < 0 + confinement = triaxial.

**SHPB (`scenario = shpb`, fdem 2D)** — `shpbIncidentLength` (2.0),
`shpbTransmitLength` (1.5), `shpbBarDiameter` (0.05), `shpbDiscDiameter` (0.05),
`shpbGap` (0), `shpbBarElemSize` (5e-3), `shpbDiscElemSize` (7.5e-4),
`shpbDiscSmooth` (8), `shpbPulse` (halfsine \| trapezoid), `shpbPulseV0` (5.2),
`shpbPulseTau` (2.2e-4), `shpbPulsePlateau` (0.5), `shpbMonitor1` (LIB−1),
`shpbMonitor2` (fin disque+1), `shpbGaugeHalfLength` (2·hBar), `shpbNoDisc` (false ;
barre seule = vérification d'onde), `absorbFactor` (1 = Lysmer classique ; 2 = éq. 21
de l'article, réfléchit 33 %), phases nommées par `shpbBarPhase`/`shpbRockPhase`
(bar/rock). ⚠️ après l'impulsion l'extrémité pilotée reste à v = 0 : ne pas prolonger
T au-delà de la rupture (rebroyage sans fin — chantier ouvert).

### 5.8 Frontières et conditions aux limites

| clé (défaut) | rôle | portée |
|---|---|---|
| `absorbing` (none) | none \| sides \| all — frontières de Lysmer + ressorts de Deeks-Randolph ; all remplace le fond encastré | tous |
| `absorbSpringFactor` (1.0) / `absorbSpringR` (W/2 latéral, H fond) | ressorts de rappel du champ lointain | tous |
| `absorbLayer` (2.2) | épaisseur de la couche absorbante DEM (×r) | dem* |
| `fixSides` (false) | flancs encastrés | fem |
| `bottomWall` (true hors tension) / `sideWalls` (false) | murs rigides DEM | dem* |
| `lateralRollers` (false) | rouleaux u_x = 0 sur les flancs (bande confinée de Yan §3.1) | fdem |
| `barV0` (1.0) / `barGaugeFrac` (0.8) | vérification onde de barre | fem, bar_wave |
| `damage` (true hors bar_wave) | endommagement on/off | fem |
| `erodeD` (0.98) / `strainCap` (0.15) | seuils d'érosion | fem |
| `stopOnRebound` (false) | (réservé — inopérant à ce jour) | fem |

### 5.9 DEM spécifiques

`packing` (hex 2D \| hcp 3D ; cubic/square pour les vérifications), `particleRadius`
(1.25e-3 2D ; 1.5e-3 3D), `bondTensile` (= ft), `bondCohesion` (= cohesion),
`bondFrictionDeg` (= frictionDeg), `bondRadiusFactor` (1.0 ; R_b = λ·r),
`bondStrengthScatter` (0 ; dispersion uniforme), `knFactor` (1.0), `ksRatio` (0.4).
⚠️ les propriétés macro d'un réseau régulier NE SONT PAS les propriétés de liaison :
calibration obligatoire avant tout usage quantitatif.

## 6. Sorties

Tous les fichiers vont dans le dossier de sortie. Fréquences : VTU toutes les
`T/frames`, une ligne d'historique tous les ~1/2000 du run.

### 6.1 `history.csv` — colonnes par scénario (mode fdem)

| scénario | colonnes |
|---|---|
| percussion / shear | `t, toolFx, toolFy, toolX, toolY, toolVx, toolVy, work, toolKE, nBroken, nFrag, detachedVol, specificEnergy` |
| tension (grips) | `t, gripFy, sigma, sigmaPeak, nBroken` |
| tension (platens) | + `epsPlaten, epsSpec, epsGauge, nBrokTen, nBrokShear, nFrag, confAchieved, peakLocked` |
| brazilian | `t, P, Pbot, drive, sigmaT, sigmaTpeak, nBroken, nFrag, sxxC, syyC, peakLocked` |
| shpb | `t, vDrive, epsM1, epsM2, sxxC, syyC, nBroken, nFrag, nInserted` |

fdem3d : `t, gripFz, sigma, sigmaPeak, nBroken` (tension) ; percussion/shear comme en
2D avec les trois composantes. fem/fem3d/dem/dem3d : variantes proches (force outil,
travail, casse). ⚠️ `sigmaPeak` est un max glissant qui attrape la sonnerie
post-rupture : recalculer les pics depuis les courbes échantillonnées, ou lire le pic
verrouillé (`peakLocked`).

### 6.2 Frames VTK (ParaView)

| fichier | contenu (champs par cellule sauf mention) |
|---|---|
| `fdem_XXXX.vtu` | maillage 2D : `vonMises, fragment, phase, grain, sigmaXX, sigmaYY, sigmaXY, epsXX` + `velocity` (nœuds) |
| `fdem_joints_XXXX.vtu` | joints (lignes) : `damage, tBreak, type` (0 intra/1 homo/2 hétéro), `ftScale, bonded, breakMode` (1 traction/2 cisaillement) + `failMode` si `writeJointMode = true` |
| `fdem3d_XXXX.vtu` | tets : `vonMises, fragment, phase, grain` + `velocity` |
| `fdem3d_joints_XXXX.vtu` | triangles : `damage, tBreak, type, ftScale, bonded, breakMode` (+ `failMode`) |
| `fem_XXXX.vtu` / `fem3d_XXXX.vtu` | `damage, vonMises, meanStress/pressure, kapDP, epvEq, ftScale, eroded` selon la loi |
| `dem*_particles/bonds_XXXX.vtu` | particules (Glyph→Sphere sur `radius`) et liaisons (`state`) |
| `frames.csv` | frame → temps et pose de l'outil (utilisé par make_gif) |

Astuces ParaView : seuiller `damage` = 1 sur les joints pour la surface de fissure ;
colorier par `grain`/`phase` pour le faciès GBM ; `bonded` = 1 montre ce que
l'insertion adaptative n'a pas encore activé.

### 6.3 Fichiers de fin de run

`fdem_final_elements.csv` (centroïdes, fragment, phase, grain),
`fdem_final_joints.csv` (`x1,y1,x2,y2,damage,type,breakMode,rn,rs,tBreak,bonded`),
`fdem_nodal_displacement.csv` (si gravité), `dem_fragments.csv`, `summary` sur stdout.

### 6.4 Le résumé stdout — les verdicts à lire

- bilan d'énergie : Ec du bloc, travail net du contact général (`gcWork`),
  **travail du dashpot de joint avec verdict** (`OK, dissipative` / `FAIL —
  INJECTED` : un chiffre positif invalide le run) ;
- adaptatif : nombre et fraction d'arêtes/faces insérées, casse ;
- GBM : fractions de phases atteintes, joints intra/homo/hétéro, **fraction
  intergranulaire de la casse** ;
- modes de rupture : « X tensile, Y shear (Z % shear) » ;
- confinement : σ_latéral visé/atteint (jauge au cœur, après équilibrage) ;
- brésilien : jauge élastique du centre (bande 0.85–1.25 = PASS), ratio de bande
  élastique, appui effectif (participation ratio), σ_t ISRM + verrouillage du pic,
  diamétralité de la fissure ;
- UCS platines : bilan de platines (|Ftop|−|Fbot|)/moyenne — quasi-statique si
  quelques % ;
- vérifications : lignes `[PASS]`/`[FAIL]` contre cibles analytiques.

## 7. Exemples types

**Percussion GBM 2D, insertion adaptative** (le premier essai de la séance) :

```
mode = fdem
scenario = percussion
mesh = voronoi
grainSize = 0.01
grainSeeding = random
refineLevels = 1
phases = quartz feldspar biotite
phase.quartz.fraction = 0.33      # + propriétés par phase, cf. §5.2
...
gbAlphaTen = 0.5
insertion = adaptive              # ← LA ligne. Tout le reste est inchangé.
toolRadius = 0.015
toolMass = 5.0
impactSpeed = 8.0
absorbing = all
```

**Traction de vérification** : `configs/verify_fdem_tension.cfg` (jointXi = 0,
cible ft exacte). **UCS par platines** : `scenario = tension`, `loading = platens`,
`pullV = -0.2`, `pullRamp` ≈ 10 transits d'onde, `ucsStopAfterPeak = true`.
**Triaxial** : + `confiningPressure = 10e6`, `confiningRamp`, `pullDelay` ≥ 3×rampe.
**Brésilien aplati** : `scenario = brazilian`, `discMesh = native`,
`discFlattenDeg = 20`, `grainMesh = delaunay`, `brazilianStopAfterPeak = true`.
**SHPB** : `configs_yan/shpb_complet_adaptatif.cfg`. Toutes les configs livrées dans
`configs/` (démos + verify_*) et `configs_yan/` (campagne article) sont commentées.

## 8. Pièges connus et règles maison

1. **Contrôle à charge nulle** après toute modification joint/contact/maillage :
   `pullV = 1e-12` → 0 joint cassé exigé, dampWork ≤ 0. Le test le plus discriminant.
2. `jointXi` : 0 en vérification de loi, 0.01 en quasi-statique, 0.05 en impact.
3. Comparaisons fines à **nombre de threads égal** ; certification à 1 thread.
4. Maillage structuré en FDEM = condition d'invalidité (trajets biaisés, divergence
   en phase débris) → Voronoï désordonné (`grainSeeding = random`).
5. En shear, l'outil doit démarrer HORS bloc (`toolX` négatif par défaut — ne pas
   copier un toolX de percussion).
6. Le brésilien sur disque plein s'amorce au contact (σ_t ≈ 0.5–0.7·ft, déficit
   structurel documenté) → disque aplati + vérifier `l_ch = E·Gf/ft²` vs R avant
   d'interpréter un BTS.
7. Un `.exe` fraîchement écrasé peut être verrouillé quelques secondes (antivirus).
8. Les paramètres des démos sont des ordres de grandeur NON calibrés ; la
   calibration (banc bayésien `tools/bayes_bench.py`, `tools/calibrate_bohus.py`)
   est un préalable à toute affirmation quantitative. Le jeu Red Bohus historique
   est invalidé (bug d'amortissement corrigé) — recalibration à refaire.
9. Reproductibilité : garantie par `seed` PAR binaire ; MSVC et libstdc++ tirent des
   nombres différents à graine égale (Voronoï, phases) — re-baseliner par plateforme.

## 9. Post-traitement fourni

`tools/plot_results.py` (champ + historiques sans ParaView), `tools/make_gif.py
<cfg> <run> [out.gif]` (animation avec outil dessiné), `tools/rockim_gui.py`
(tracés intégrés), `tools/export_abaqus.py <run> <out.inp>` (maillage frame 0 +
champ ftScale → deck Abaqus mm-t-s-MPa pour validation croisée iso-maillage),
`tools/yan_point.cpp` (∫f(D)dD en précision machine), `tools/verify_suite.py` (§3.3),
`tools/make_unstructured_mesh.py` (maillages simplexes non structurés uniformes via
Gmsh — `box3d W D H h out.msh [seed]` / `box2d W H h out.msh [seed]` — pour
`mesh = file` ; `pip install gmsh`).

**Règle maison maillage (2026-08-11)** : le maillage de BASE de toute étude est
DÉSORDONNÉ (`mesh = file` non structuré, ou `mesh = voronoi` si le sujet est le
GBM). Les grilles régulières sont réservées aux vérifications qui en ont besoin
par construction ; le solveur imprime un WARNING si un scénario de fissuration
part sur `mesh = grid` (mesuré : la grille de Kuhn 3D intrinsèque part en cascade
énergétique en phase débris là où le même cas sur maillage non structuré est sain).
