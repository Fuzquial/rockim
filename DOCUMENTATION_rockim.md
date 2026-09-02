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
| `meshFile` (requis si mesh = file) | chemin d'un Gmsh MSH 2.2 ASCII (type 2 en 2D, type 4 en 3D) ; boîte translatée à l'origine, W/H/D relus de l'enveloppe ; générer via `tools/make_unstructured_mesh.py` — variantes `box3d`, `box2d`, `bench1` (bloc + insert spherique), **`bench1g`** (idem GRADUE : `bench1g W D H R gap h hIns hFin rFin dFin out.msh [seed]`, champ de taille en rampe, fin dans un cylindre de rayon `rFin` et profondeur `dFin` sous l'axe d'impact — resserrer sur le rayon de contact de Hertz, pas sur l'etendue du champ visible), `tunnel` | fdem, fdem3d |
| **groupes physiques** (V1, 3D) | si le MSH porte des `$PhysicalNames` (dim 3), chaque volume physique devient un **corps** : AUCUN joint cohésif entre deux groupes (les faces deviennent extérieures, l'interaction passe par le contact général), matériau du corps = **phase homonyme** (ou `groupPhase.<nom> = <phase>` ; sans correspondance : phase 0 + WARNING — `phases` avec mesh = file **exige** des groupes nommés), `groupVel.<nom> = vx vy vz` (vitesse initiale du corps), `trackGroup = <nom>` (colonnes history : `grpZ`,`grpVz` — centroïde massique + vitesse moyenne — et V2/B2 : `grpFx,grpFy,grpFz` — force de contact NETTE sur le corps au pas courant, sommée dans les deux lois de contact : la F-δ se lit en direct — et `grpSzz` — jauge σzz moyenne volumique du corps), résumé par corps (KE, vz, masse) en fin de run. Avec `toolShape = none` (percussion), l'outil analytique est retiré : l'impacteur est un corps MAILLÉ du fichier — générateur fourni : `make_unstructured_mesh.py bench1 W D H R gap h hIns out.msh [seed]` (bloc + insert sphérique séparés de `gap`, groupes `rock`/`insert`) ; contrôle deux-corps au repos : `zeroload_bench1_3d` (0 casse, gcWork = 0 exact) | fdem3d |
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
| **`jointResidualMu`** (< 0 = non posée) | **Le coefficient de frottement RÉSIDUEL du joint rompu.** Le coefficient glisse du **pic** `tan(frictionDeg)` vers ce résiduel par la **même f(D)** que la cohésion : μ_eff = μ_res + (tanφ − μ_res)·f(D). rockim gardait jusqu'ici le frottement de **pic à vie**, ce qui verrouille une zone broyée sous forte compression. C'est la distinction que fait **Y-Geo** (AbuAisha et al. 2015, éq. 7.5 : un angle de frottement de **fracture** φ_f distinct de l'angle interne du pic) et que **Solidity** obtient autrement, en remettant le joint rompu au contact et à son glissement — **0,6** pour le calcaire, **0,18** pour le granite de Kuru, contre un frottement de pic de **1,85** : un facteur **10,3** entre pic et résiduel, et le papier granite dit explicitement que ce coefficient bas est ce qui permet aux fragments d'être éjectés et de cesser de porter le taillant. **`jointResidualMu` GÉNÉRALISE `jointFrictionScaled`** — μ_res = tan(frictionDeg) redonne le défaut, μ_res = 0 redonne `jointFrictionScaled = 1` — les deux clés sont donc **exclusives** (le run s'arrête si les deux sont posées). Les deux égalités sont **exactes**, vérifiées et verrouillées par `residualmu_equiv_defaut_2d` et `residualmu_equiv_scaled_2d` |
| **`jointShearUnload`** (plastic) | plastic \| **origin** = décharge ET recharge en cisaillement sur la **sécante à l'origine** passant par (s_max, τ_env(s_max)), éq. 18 de Yan et al. — symétrique exact de l'éq. 17 du mode I. `plastic` (défaut, inchangé) est une plasticité à retour radial : la décharge suit la sécante de pénalité et le glissement plastique est conservé. Les deux **coïncident en charge monotone** (le glissement au pic est le s_p = (c + tanφ·\|σ_n\|)/p de Munjiza, donc l'endommagement de mode II démarre au même instant) et ne diffèrent qu'à la décharge. ⚠️ l'éq. 18 place **tout** le cap dans la sécante, frottement de Coulomb compris : avec `jointFrictionScaled = 0` le glissement frottant devient réversible (aucune boucle d'hystérésis). **Forme littérale de l'article = `origin` + `jointFrictionScaled = 1`** |
| **`insertion`** (intrinsic) | intrinsic \| **adaptive** = insertion dynamique extrinsèque (Yan et al. 2023) : aucun joint à t = 0, liaison cinématique exacte, activation quand σ_n ≥ ft ou \|τ\| ≥ c − σ_n·tanφ, continuité de contrainte à l'insertion. Gains mesurés : dt ×2, mur ×2.2–2.7, complaisance nulle |
| **`gcActivation`** (full) | full \| **adaptive** = activation adaptative des faces de contact (Fukuda et al.) : `act_` ne contient que les faces qui **peuvent** toucher, au lieu de tout l'extérieur balayé à chaque pas. Trois règles, activation **monotone** : (C) peau endommagée — l'élément porte un joint cassé/mort, plus un anneau par sommet ; (A) autre corps à moins de `gcActMargin` cellules (composantes connexes par union-find sur les joints porteurs, recalculées quand nBroken change) — c'est ce qui arme le SHPB multi-corps dès t = 0 ; (B) voisinage d'une face ayant déjà **porté** une force (une face qui racle propage, une face inerte non). Balayage cadencé par v_max, borné par `gcActEvery`. Les faces libérées par joints morts entrent au **même pas** qu'en mode full (cache). **Mesuré : percussion 3D ×2,32 bit-identique** (1130 → 488 s, 4 % des faces activées), percussion 2D bit-identique, UCS −15 % aux mêmes chiffres, SHPB identique sur 83 % du run puis enveloppe chaotique (contrôle : full sous OMP=2 diverge 8× plus tôt). Approximation assumée (la même que Fukuda) : un continuum **intact** ne se replie pas sur lui-même |
| `gcActMargin` (2.0) | marge d'activation des règles A/B, en multiples de la cellule de détection |
| `gcActEvery` (64) | cadence maximale du balayage d'activation [pas] |
| **`contact`** (penalty) | penalty \| **potential** = contact général par **potentiel de Munjiza** (éq. 2-5 de Yan et al. 2023), **2D et 3D**. Paires d'**éléments** (et non nœud-face) : force normale distribuée F = p·∮(φ_A−φ_B)·n dΓ sur le bord du recouvrement — polygone triangle-triangle en 2D (φ = 3·min λ), **polyèdre tet-tet** en 3D (φ = 4·min λ, clip par les 4 demi-espaces + face de coupe reconstruite). Intégration **exacte** (subdivision aux plans de médiane : 6 en 2D, 12 en 3D), lumping nodal consistant, 3e loi de Newton **machine**, champ **conservatif** — collisions élastiques : ΔKE/KE₀ = 3,7e-12 (2D), 2,0e-8 (3D), transfert exact (selftest-potential2d/3d). Frottement tangentiel incrémental à ressort + cap de Coulomb (éq. 4-5, vectoriel en 3D), historique par paire. Détection O(N) type NBS (binning AABB), exclusion des paires liées par un joint **vivant**, compose avec `gcActivation`/`gcXwindow`. Relève de naissance par **aire/volume** de recouvrement (pen0_ du potentiel, τ = `gcBirthTau`) : une paire née en recouvrement (joint mort comprimé) ne matérialise pas son énergie potentielle — signe absorbant garanti (une rampe temporelle ferait l'inverse, mesuré +179 J/m). Gardes 3D : plancher de volume relatif (1e-12·min V) + contrôle de **fermeture** du polyèdre (les tets exactement tangents produisaient des slivers à faces non refermées — 5 joints cassés à charge nulle, attrapés par le contrôle zeroload). ⚠️ le gcWork peut porter un petit résidu positif (biais O(dt) du compteur + relève) — annoté dans le résumé, pas une pathologie en potentiel. L'outil analytique reste en pénalité (un outil MAILLÉ passe sous le potentiel via les groupes physiques + `toolShape = none`). SHPB : onde incidente identique au penalty à 3e-6 près ; zone broyée **conservative** → rebond plus élastique (percussion 2D : e 0,55 → 0,71, moins de casses) — écart de loi physique assumé. Coût par paire supérieur au nœud-face : combiner avec `gcActivation = adaptive`. **Perf (N1, 2026-08-14, tout bit-neutre)** : grille dense à seaux réutilisés + **ordre canonique des paires** (tri (eLo,eHi) — les sommes de forces ne dépendent plus de l'ordre de découverte ; réf shpb_mini_potential recalée 595 → 578, dernier changement d'ordre autorisé), **pré-filtre SAT complet 3D** (8 plans de faces + 36 axes d'arêtes croisées, cache du dernier axe séparateur par paire à la Baraff — jeu complet : s'il ne sépare pas, le recouvrement est réel), clip **sans copie** (ping-pong de pointeurs — l'ancien `P = Q` déplaçait ~9 Ko ×4 par clip). Compteurs `potential stats` au résumé (paires / joint-vivant / sep-hint / sep-face / sep-arête / clip-vide / clip-force, tGrid / tLoop) : sur percussion 3D T = 5e-5, 230 M de paires → 61 % réglées par le cache d'axe, 22 % de **clips VIDES** (~4 µs chacun — contacts rasants : recouvrement réel sous plancher). Mesures : 682 s (grille naïve) → 643 (seaux) → 477 (SAT faces) → 448 s (ping-pong) à T = 5e-5. Sur la **longue** T = 2e-4 (3 474 s vs 488 s pénalité, ~7×), les compteurs renversent le tableau : le poste dominant est l'**intégration exacte des 33 M de clips AVEC force** (paires de débris en contact permanent, ~70 µs pièce — la subdivision aux 12 plans — ≈ 2/3 du run), les clips vides ne pèsent que ~12 %, et le scan SAT complet n'y sépare plus rien (37 k sur 479 M — les axes d'arêtes sont dispensables en régime débris). Le critère N1 (≤ 1,3×) demande donc une refonte de l'INTÉGRATION en régime de contact persistant (quadrature moins chère = perte d'exactitude à arbitrer, warm-start du polyèdre, cadence des contacts stationnaires) — pistes au plan v2, décision à prendre |
| `potPenaltyFactor` (1.0) | pénalité normale du potentiel, en multiples de E·épaisseur (2D) / de E (3D) |
| `potTangentFactor` (1.0) | raideur tangentielle de l'éq. 4-5, en multiples de E·épaisseur (2D) / de E·hmin (3D) |
| `insertionPenaltyFactor` (4) | pénalité des joints ACTIVÉS en mode adaptatif (décharge/contact) |
| `jointWeibullM` (0 = off) | m > 1 : ft et cohésion de chaque joint × facteur Weibull(m) de moyenne 1 (Gf non tiré) |
| `strengthCorrLength` (0) | 0 = tirages indépendants ; > 0 = champ gaussien corrélé (copule) de cette longueur [m] |
| `strengthCorrLengthB` (= A) + `strengthCorrAngleDeg` (0) | anisotropie orientée (foliation ; en 3D plan de texture incliné autour de y) |
| `fieldSeed` (seed+777) | graine du CHAMP, indépendante du maillage — deux maillages voient les mêmes zones faibles |
| **`jointSizeEffect`** (0 = off) | **effet d'échelle statistique de Weibull** : `ft` et cohésion × (Zeff/V_J)^(1/m), **exactement la formule des VUMAT DP-DFH d'Abaqus** (`sig_k = sigw*(Zeff/V_el)**(1/m)`, `VUMATS/dfh/vumat_kstdfh_psivar.f`), le « point matériel » d'un joint étant la face entre deux éléments : V_J = moyenne des deux volumes adjacents (2D : moyenne des aires × `thickness`). **Obligatoire pour toute étude d'objectivité STRUCTURALE** — sans lui, raffiner désactive l'effet d'échelle (rapport DP-DFH §13.1, éq. 42 ; la variante à V_el figé de `vumat_psivar_rc99_veff1.f` est réservée aux contrôles au point matériel). Se compose multiplicativement avec la dispersion `jointWeibullM` : facteur total = Weibull(moyenne 1) × taille, replié dans `J.stat`, donc le champ `ftScale` des VTU montre le total. **Gf n'est PAS recalé** (énergie de fissuration = propriété du matériau) : monter `ft` à Gf fixé raccourcit la branche adoucissante, `dnF`/`slipF` sont recalculés. Vérifié : deux maillages dans le rapport 3,16 donnent des facteurs dans le rapport (3,16)^(1/24) à **0,03 %** près |
| `jointZeff` (1e-9 m³ = 1 mm³) | volume de RÉFÉRENCE auquel `ft`/`cohesion` de la config sont déclarés — même défaut que le `Zeff` des VUMAT (échelle de l'indentation). Doit être une constante **physique**, jamais déduite du maillage : un Zeff qui suivrait la moyenne du maillage ramènerait le facteur moyen à 1. ⚠️ en 2D, V_J dépend de `thickness`, souvent conventionnel (1 m) : déclarer Zeff en cohérence (le log imprime V_J) |
| `jointSizeEffectM` (= `jointWeibullM`) | exposant m du recalage, séparé pour permettre l'effet d'échelle **sans** dispersion, ou un m ≠ celui du tirage. Pour Bohus : **m ≈ 24** (§14.1 — la pente apparente m ≈ 6 vient d'une seconde population de défauts à écarter) |
| `historyFlush` (true) | vide `history.csv` après **chaque** ligne. Sans lui l'OS bufferise : le fichier reste vide jusqu'à la fin (impossible de suivre un run) et un run **tué** laisse une dernière ligne tronquée au milieu du tampon — constaté le 2026-08-14 sur `out_banc_mid` (26 colonnes au lieu de 28, terminées par `,-`). `histEvery` borne les lignes à ~2000 par run, le coût est négligeable. Vérifié : run normal 2127 lignes / 0 incomplète ; run tué par SIGKILL 496 lignes récupérées / **0 incomplète**. Purement I/O, bit-neutre |
| `jointSizeEffectClamp` (5) | bornage du facteur à [1/5, 5] avec avertissement et comptage — garde-fou contre un maillage très hétérogène, un Zeff mal choisi ou une épaisseur 2D non physique |
| `crushCap` (8·cohesion) | plafond élasto-plastique du déviateur du bulk (garde-fou, désactivé si une `law` est active) |
| `bulkDamage` (off) | **pulvérisation** (Yang et al. 2026, IJRMMS 206, éq. 3-4) : dégradation de raideur des tétraèdres, σ = Cd·(1−D)·σ̄, D linéaire (Camanho) en δm = h_e·ε_vm entre `bulkDamageDelta0` et `bulkDamageDeltaF` [m], irréversible, plafonné à `bulkDamageDmax` ; **2D et 3D** (déformation plane : ε_zz = 0 entre dans le déviateur), `law = elastic` seul. S'AJOUTE au crushCap (principe VIII) — le deck granite neutralise ce dernier (1e12). Colonnes `nPulv,bdWork` + champ VTU `bulkD` quand armé. Dissipation Y·dD ventilée dans le poste éléments |
| `bulkDamageDelta0` (1.4e-5) / `bulkDamageDeltaF` (4.0e-4) / `bulkDamageDmax` (0.9) / `bulkDamageCd` (1.0) | calibration Kuru Grey de l'article (leur 0,014/0,4 lus en mm, éléments de 1 mm) |
| `groupBond.<A>.<B>` (—) | **liaison entre corps** (3D, mesh = file) : l'interface conforme entre deux volumes physiques nommés reçoit des joints cohésifs (type GBM frontière, moyenne des phases × facteurs gb*) au lieu d'être remise au contact — le brasage insert/bit de la spec 005. Valeur : `joints`. L'insertion adaptative lie les nœuds de l'interface comme partout (rebindVertex) |
| `trackGroups` (—) | 3D, mesh = file : colonnes `z_<nom>,vz_<nom>` (centroïde massique, vitesse moyenne) par corps listé — vitesses d'indentation et de rebond du bit (spec 005). S'ajoute au `trackGroup` singulier existant |
| `gauge.<nom>` (—) | 3D : `"z0 z1"` — colonne `szz_<nom>`, σ_zz moyenné en volume dans la tranche [z0,z1] du corps (la jauge à mi-bit de leur fig. 8) ; tranche figée en configuration de référence |
| `jointSoftening = munjiza` | **alias** de `yan` : la f(D) de Yan et al. 2023 EST la z-curve de Munjiza 2004 (a = 0,63, b = 1,8, c = 6, ∫f dD = 0,386307), celle de Y-Geo et de Solidity (Yang et al.). Avec `jointShearUnload = origin`, le moteur √(rn²+rs²) de cette branche est l'ellipse mode I-II exacte de leur éq. 3 — le modèle cohésif de l'article est donc INTÉGRALEMENT disponible, insertion adaptative comprise |
| **`jointDeath`** (separation) | **QUAND le joint passe la main à l'algorithme de contact.** `separation` (défaut, historique) : le joint ne meurt qu'une fois franchement ouvert (`dnMax > 3·dnF`) ; un joint broyé qui glisse en compression reste **vivant** et sert de contact frottant de ses propres lèvres. `damage` : il meurt dès que **D ≥ 1**, quel que soit le signe de l'ouverture — la règle de **Guo (thèse Imperial 2014, §2.3.3)** : « the stress-displacement relation is not applied to this failed joint element anymore ; instead, the interaction between the fracture walls will be counted as contact forces that are calculated by the contact algorithm ». C'est ce relais qui, chez eux, achemine les 32 J de frottement entre fragments (65 % du budget d'impact, ARMA 2024). Sorties : ligne `relais joint->contact` au résumé — joints morts, part morte **en compression**, et charge normale lâchée au relais |

**Ce que le relais change, mesuré** (UCS `configs_yan/ucs_adap.cfg`, 2026-08-25,
`separation` → `damage`) :

| poste | separation | damage | |
|---|---|---|---|
| joints morts | 143 | 270 | |
| dont **en compression** | 11 (7,7 %) | 162 (60 %) | |
| charge lâchée au relais | 110 kN/m | 4 169 kN/m | ×38 |
| travail de contact | 0,765 J/m | **24,68 J/m** | **×32** |
| dont **frottement** | 0,0725 J/m | **2,160 J/m** | **×30** |
| UCS | 51,0395 MPa | 51,0395 MPa | inchangé |
| part de cisaillement | 48,6 % | 60,7 % | |
| résidu du bilan | −2,64e−12 J/m | −2,30e−12 J/m | OK |

Le relais achemine donc bel et bien la dissipation vers le contact — c'est le
mécanisme qui manquait — **sans dégrader le bilan d'énergie**. Et cet UCS tourne
à `contactMu = 0,1` seulement ; l'impact est à 0,6.

⚠️ **La crainte historique est levée, mais elle était fondée.** Le commentaire
du site de mort disait : *« killing it by slip hands interpenetrated faces to
the general contact, whose penalty then releases ½ k pen² of energy created from
nothing »*. C'était vrai avant la **relève de naissance `pen0_`** ajoutée au
chantier A3 ; le résidu mesuré ci-dessus montre qu'elle la neutralise.

⚠️ **Réserve ouverte.** Les 4 169 kN/m lâchés ne créent pas d'énergie mais
**disparaissent du chemin d'effort** le temps que `pen0_` décroisse
(`gcBirthTau`). Sans conséquence sur l'UCS, dont le pic précède le relais. À
mesurer sur l'impact, où le chemin d'effort sous l'insert est justement l'enjeu :
si le déficit s'y voit, il faudra une **continuité de traction** au relais,
miroir du `dn0` de l'insertion adaptative.

⚠️ **Constat sur le mode `separation` lui-même.** Il ne garantit PAS l'absence de
mort en compression : 11 joints sur 143 y meurent comprimés, parce que `dnMax`
est le **maximum sur les points d'intégration** — une interface en flexion,
béante d'un côté et comprimée de l'autre, franchit `dnMax > 3·dnF` avec une
résultante normale encore compressive.

### 5.4 bis Effets de vitesse : viscosité de volume et DIF

*Section ajoutée le 2026-08-25. Ces clés existaient depuis le 2026-08-18 et
n'avaient jamais été documentées — dette du principe VII soldée à l'occasion du
chantier « DIF intrinsèque ». Le code fait foi ; chaque ligne ci-dessous a été
relue dans `src/FdemSolver.cpp` et `src/Fdem3dSolver.cpp`.*

**Viscosité de volume** — une contrainte visqueuse newtonienne **2 μ D** (D = taux
de déformation co-rotée) est ajoutée au tenseur de Cauchy de chaque élément.
C'est le terme de l'éq. 6 de Yan et al. 2023, et c'est aussi le `η·D` de
l'éq. 2.6 de la thèse de Guo (Imperial College, 2014) dont le code Solidity de
Yang et al. est issu. ⚠️ **Attention à la convention du facteur 2** : rockim
applique `2 μ D` là où Guo écrit `η D`, donc **η = 2 μ**. Pour reproduire un η
publié, poser `bulkViscosity = η/2`.

| clé (défaut) | rôle | portée |
|---|---|---|
| `bulkViscosity` (0 = off) | μ **littéral** [Pa·s], le même pour tous les éléments. Exclusive avec `bulkViscosityXi` (le run s'arrête si les deux sont posées) | fdem, fdem3d |
| `bulkViscosityXi` (0 = off) | μ **calculé du maillage** : μ = ξ·h·√(E ρ) par élément. ξ = **2,0 vaut le critique de Munjiza** 2h√(Eρ) — c'est la valeur de la Table 1 de Yan et al. Le résumé imprime « soit 0,5·ξ × le critique » | fdem, fdem3d |
| `bulkViscosityGraded` (0) | 1 = μ **gradué** par élément (chaque tétra son h) ; 0 = μ **global**, pris à la médiane. Sur un maillage gradué, μ global fait payer le pas de temps du plus fin tétra partout — mais c'est la forme d'un η constant publié | fdem, fdem3d |
| `viscousInInsertion` (1) | 1 = le terme visqueux entre dans la contrainte d'essai du **critère d'insertion** ; 0 = le critère ne voit que la contrainte élastique. Argument du 0 : sinon le taux agit deux fois, comme contrainte d'essai ET comme seuil via le DIF. ⚠️ **clé 3D seulement** | fdem3d |

Le pas de temps porte une borne **diffusive** ρh²/4μ en plus de la borne
élastique : monter μ coûte du dt. Le travail visqueux est compté dans
`viscWork_`, **ventilé à l'intérieur du poste « éléments »** du bilan B4 (ce
n'est pas un poste de plus) et imprimé au résumé de fin de run avec son verdict
de signe. ⚠️ Il n'a **pas de colonne dans `history.csv`** : sur un run tué avant
la fin, la part visqueuse est irrécupérable.

**DIF (Dynamic Increase Factor)** — les résistances de joint sont multipliées par
un facteur fonction du taux de déformation, éq. 2 et 3 de Yang et al. 2025.
`DIF_traction` multiplie `ft` **et** `Gf` ; `DIF_compression` multiplie `cohesion`
**et** `GfII` — comme eux. Comme ft et Gf reçoivent le même facteur, la
**longueur de la branche adoucissante** kI·Gf/ft est invariante : seule la limite
élastique dnE = ft/pj bouge.

| clé (défaut) | rôle | portée |
|---|---|---|
| `strainRateDIF` (off) | off \| `yang` = leur éq. 3 **littérale**, exposant 0,07 \| `yang-fig2` = exposant **0,1707** déduit de leur figure 2b. ⚠️ L'exposant 0,07 imprimé ne raccorde pas la loi à ses bornes : elle saute de 1,516 à 1,85 en ε̇ = 10² /s, et en insertion extrinsèque ce saut est un **attracteur** (la population insérée s'empile juste sous 10² /s — mesuré : médiane 99,36 /s contre 40,22 avec `yang-fig2`). Trois repères de la suite verrouillent ce comportement | fdem, fdem3d |
| `strainRateTau` (1e-6 s) | constante de temps du **filtre exponentiel** du taux par élément (ε̇ = max des valeurs propres absolues de D co-rotée). Doit être > 0 si le DIF est actif | fdem, fdem3d |
| **`strainRateDIFArm`** (insertion) | **QUAND** le facteur est figé. `insertion` (défaut, comportement historique) : à l'instant de l'insertion — **exige `insertion = adaptive`**. `envelope` (2026-08-25) : au moment où le joint **quitte sa branche élastique**, c'est-à-dire là où il commence à s'endommager — **exige `insertion = intrinsic`**. Les deux sont l'analogue l'un de l'autre : en adaptatif le joint NAÎT au pic de l'enveloppe (continuité de contrainte, `dn0`), naissance et amorçage coïncident donc par construction ; en intrinsèque le joint est déjà là et seul l'amorçage subsiste. La table de validation refuse explicitement les deux croisements (`envelope` + adaptatif appliquerait le facteur deux fois ; `insertion` + intrinsèque est l'erreur historique, dont le message oriente désormais vers `envelope`) | fdem, fdem3d |

**Pourquoi l'armement intrinsèque ne peut PAS réutiliser le critère en
contrainte d'élément** (mesuré le 2026-08-25, gardé ici pour que le piège ne
soit pas retenté) : la première version armait sur le critère de
`insertionSweep()` — la contrainte moyenne des deux éléments contre l'enveloppe
de Mohr-Coulomb, exactement le critère de l'insertion adaptative. Elle ne
s'arme **jamais** : 0 joint gelé sur 6840, et 100 % des joints sollicités
s'endommagent sans DIF. La raison est structurelle et non un réglage : en
schéma intrinsèque le joint est le maillon faible et **écrête la contrainte que
ce critère surveille**, si bien que la moyenne des deux éléments n'atteint
jamais ft. Le critère partagé avec l'adaptatif est donc inutilisable en
intrinsèque, et l'armement porte sur la cinématique propre du joint.

Sorties : le résumé imprime `DIF intrinseque (armement a l enveloppe): N / M
joints geles ; K joints endommages SANS DIF`. **K est le contrôle falsifiable
de l'armement** — il vaut 0 par construction, et une valeur non nulle signale
que le critère d'armement a dérivé par rapport à la loi de joint. Repères
`dif_intrinseque_2d` (fast) et `dif_intrinseque_3d` (full), plus le contrôle à
charge nulle `zeroload_dif_intrinseque_2d` (aucun joint armé sous charge nulle).

**Enveloppe de cisaillement du joint**

| clé (défaut) | rôle | portée |
|---|---|---|
| `jointShearEnvelope` (yan) | `yan` = son éq. 8, le terme de frottement tombe à **zéro dès que la contrainte normale est en traction** ; `yang` = l'**éq. 1 de Yang et al.**, il décroît jusqu'au cut-off en ft : fs = c − tanφ·min(σn, ft). Les deux **coïncident exactement en compression** et ne diffèrent qu'en traction, où la forme de Yang AFFAIBLIT le cisaillement (−34 % au cut-off sur le banc de percussion). C'est ce qui gouverne le partage traction/cisaillement dans les zones tendues, donc le faciès radial. **La forme de l'article est `yang`** | fdem, fdem3d |
| `meanTensionCapFactor` (0 = off) | plafond sur la contrainte moyenne de l'élément, en multiples de `ft`. Garde-fou rockim, sans équivalent dans la littérature de référence : laisser éteint pour toute réplique | fdem, fdem3d |

### 5.4 quater Loi de joint : les deux dernières conventions de Guo

*Ajouté le 2026-08-25. Avec `jointShearEnvelope = yang`, `jointSoftening = yan`,
`jointShearUnload = origin` et une pénalité de 26,32 E/h, ces deux clés
achèvent le portage de la loi de joint de Solidity.*

| clé (défaut) | rôle | portée |
|---|---|---|
| **`jointElastic`** (linear) | `parabolic` = **Guo éq. 2.31** : la branche élastique vaut σ = ft·(2r − r²) avec r = δn/δnE, au lieu de la droite σ = pj·δn. Elle arrive au pic avec une **tangente nulle** — transition douce vers l'adoucissement, là où rockim a un coude — et part de l'origine avec la pente **2·pj**, des deux côtés de δn = 0 (la loi est C¹ à l'origine, la branche de compression devenant σ = 2·pj·δn, première ligne de son éq. 2.31). ⚠️ **Exige `jointSoftening = yan` ou `munjiza`** : la parabole va avec la z-curve et n'est implémentée que sur ce chemin. La combinaison est **refusée** plutôt que laissée sans effet | fdem, fdem3d |
| **`jointDeltaC`** (exact) | `guo` = **Guo éq. 2.30** : δc = 3·Gf/f mesuré **depuis zéro**, au lieu de δnE + Gf/(ft·∫f dD). Il approxime l'intégrale de la z-curve par 1/3 là où elle vaut **0,386307** : son modèle dissipe donc **1,159 fois son Gf nominal**. C'est SA convention, et ses Gf publiés ont été calibrés avec — il faut la reproduire pour retrouver ses chiffres | fdem, fdem3d |

**D'où vient le 26,32.** Leur « Penalty Number » de 3 000 GPa n'est qualifié
nulle part dans l'ARMA. Il est tranché par la citation que fait Guo à propos
de la pénalité : **Turon, Dávila, Camanho & Costa (2007)**, *Eng. Fract.
Mech.* 74:1665-1682 — le papier du rapport classique des modèles de zone
cohésive, **K = α·E/t avec α ≈ 50**. Or 3 000 / 57 = **52,6**. C'est donc la
pénalité des éléments **cohésifs**, rapportée au module de la **roche**, posée
par la règle de Turon — ni le contact, ni le carbure. Guo eq. 2.25 posant
δnp = 2·ft·h/p0, la raideur vaut p0/(2h) et l'équivalent rockim est
p0/(2E) = 26,32, le facteur 2 venant de **sa** convention.
⚠️ À noter : l'éq. 2.28 de Guo recommande E ≤ p0 ≤ 10E, ce qui **contredit**
le α ≈ 50 de Turon qu'il cite deux phrases plus haut. Les auteurs de l'article
ont suivi Turon, pas la thèse.

**Pourquoi la parabole rend la pénalité cohérente.** L'équivalence de pénalité
(§5.4, `jointPenaltyFactor` ≈ 26,32 pour leur p0 = 3 000 GPa) a été établie en
faisant coïncider **l'ouverture au pic** δnE = δnp. Avec la branche linéaire,
cela laisse la **raideur initiale** à la moitié de la leur. Avec la parabole, la
pente à l'origine vaut 2·ft/δnE : les deux quantités coïncident alors
**simultanément**. Les deux clés vont donc ensemble.

Mesures sur `verify_fdem_tension.cfg`, sous `jointSoftening = yan` :

| | err_pct | casses |
|---|---|---|
| yan seul | −1,70281 % | 24 |
| + `jointElastic = parabolic` | −2,49639 % | 24 |
| + `jointDeltaC = guo` | −2,33108 % | 24 |

Le **nombre de fissures ne bouge pas** : ces conventions déplacent la
complaisance et l'énergie dissipée par fissure, pas le compte.

⚠️ **Reste non porté : la quadrature.** Guo intègre le joint sur trois points
aux **milieux d'arêtes** (sa Table 2.2, poids 1/3) ; rockim intègre aux
**nœuds**. Les deux sont des règles à trois points de poids égaux, exactes pour
une variation linéaire ; elles ne diffèrent que sur la part non linéaire, donc
dans l'adoucissement. Non implémenté : cela demande de redéfinir les points
d'intégration, et l'état par point (`omax`, `smax`, `slip`) avec eux.

### 5.4 ter La loi de VOLUME : `bulkModel`

*Ajouté le 2026-08-25 — point 4 du tableau de comparaison à Yang et al.*

| clé (défaut) | rôle | portée |
|---|---|---|
| **`bulkModel`** (corotational) | `corotational` (défaut, historique) : décomposition polaire, déformation de **Biot** ε = sym(RᵀF) − I, σ = λ tr(ε) I + 2μ ε, assemblage P = R·σ. Exact en grandes **rotations**, valable en petites **déformations** seulement. `neohookean` : la loi de **Guo** (thèse Imperial 2014, **éq. 2.6**), celle du code **Solidity** de Yang et al. — `T = (μ/J)(B − I) + (λ/J)·ln(J)·I` avec B = FFᵀ et J = det F — assortie de l'assemblage **exact** P = J·T·F⁻ᵀ. Incompatible avec `law` (qui remplace déjà toute la loi de volume) | fdem, fdem3d |

**C'est un portage, pas une invention.** La formule est citée verbatim de la thèse
qui décrit leur code. La loi est **hyperélastique** — elle dérive de
W(F) = (μ/2)(tr B − 3) − μ·ln J + (λ/2)(ln J)², le néo-hookéen compressible de
Simo-Hughes — donc conservative, et la configuration initiale y est **naturelle**
(W(I) = 0, dW/dF(I) = 0).

**Elle redonne l'élasticité linéaire au premier ordre**, avec les *mêmes* λ et μ.
C'est un remplacement continu, pas un modèle concurrent. Écart vérifié
analytiquement hors solveur, en déformation uniaxiale :

| ε | −0,40 | −0,30 | −0,10 | +0,01 | +1e−4 |
|---|---|---|---|---|---|
| écart néo-hookéen / linéaire | **+59,8 %** | +37,6 % | +9,4 % | −0,82 % | −0,008 % |

Le signe compte : **en compression la loi se raidit**. Le terme (λ/J)·ln J diverge
quand J → 0, donc le matériau oppose une barrière infinie à l'écrasement et
l'élément ne peut plus s'inverser — ce que la loi linéaire ne fait pas, et c'est
la raison d'être du `crushCap`, garde-fou qui n'existe dans aucun code de
référence. Sous l'insert, det F tombe à **0,5–0,7** : c'est précisément là que
les deux lois cessent d'être interchangeables.

**L'assemblage vient avec, et c'est le point 5 du tableau.** La forme
co-rotationnelle assemble une contrainte de Cauchy sur une aire de **référence** :
il lui manque exactement le transport d'aire de Nanson, cof(U) = J·U⁻¹. Le
facteur d'écart est **J^(−2/3) en 3D** — soit +40,6 % sur la force interne à
det F = 0,6 — mais **J^(−1/2) en déformation plane**. ⚠️ Ne jamais écrire cet
exposant en dur : rockim passe par la forme générique `P = J·R·σ·U⁻¹`, correcte
dans les deux dimensions. Le signe de det F est conservé partout dans le chemin
des forces (en prendre la valeur absolue retournerait la force d'un élément
inversé et l'enfoncerait davantage) ; à det F ≤ 0 le solveur retombe sur
l'assemblage co-rotationnel.

En déformation plane, J = det(F₂ₓ₂) exactement et **T_zz = (λ/J)·ln J**, purement
volumique — et non la relation de Poisson ν(σ_xx + σ_yy) de la branche linéaire.

Repères : `bulkmodel_neohooke_2d` et `zeroload_neohooke_2d` (fast),
`bulkmodel_neohooke_3d` (full — indispensable, l'exposant de l'écart diffère
entre dimensions).

### 5.4 quinquies Les conventions lues dans le CODE de Solidity (2026-08-26)

> ⚠️ **AVERTISSEMENT DE LECTURE — porté le 2026-08-30, après contre-audit.**
> *À lire avant d'utiliser une seule des clés de cette section, et avant de citer
> une seule de ses lignes de code dans un manuscrit.*
>
> **1. La source est bien celle d'Imperial** — dépôt public
> `ImperialCollegeLondon/solidity-solver-open`, LGPL-3.0, lu le **2026-08-26**.
> Ce point a été contesté en interne pendant trois jours puis rétabli : voir
> [`chantier_imperial_2026-08-29/A03_resourcer_attributions.md`](etat_de_l_art/chantier/A03_resourcer_attributions.md)
> §2. Les noms de valeur `solidity` sont donc **exacts** et ne seront pas renommés.
>
> **2. Mais ce n'est PAS la version qui a produit l'article de 2026.** Le facteur
> d'endommagement d'élément y est câblé à zéro (`Y3Dfd.c` l. 749-751, `df = R0`)
> et le DIF y est neutre (`dpeftdif = R1`), alors que l'article publie les
> équations (3)-(4) d'un modèle d'endommagement. La lecture la plus simple :
> **version ouverte en retard sur la version interne** — banal pour un code de
> recherche. **Conséquence de méthode : lire une FORME ici et en conclure une
> implémentation de ce que décrit l'article de 2026 est une faute.** Ce n'est pas
> « le code de quelqu'un d'autre » — c'est bien le leur, même lignée, mêmes
> auteurs — c'est *une autre version*.
>
> **3. Trois statuts, et non deux.** Pour chaque convention ci-dessous, ne pas
> confondre : ce que disent les **articles publiés** ; ce que fait le **code
> public** ; ce que fait la **version interne** (inconnue, non consultable). Les
> relevés de cette section sont **tous du deuxième type**, sauf là où une source
> d'article est explicitement nommée (`jointFailRule`, `gcBirth` — voir ci-dessous).
>
> **4. Les citations `Y3D*.c l. NNNN` de cette section ne sont pas reproductibles
> telles quelles.** Le dépôt est activement maintenu (dernier push relevé le
> 2026-03-31) : **les numéros de ligne bougent.** Ils valent pour l'état lu le
> **2026-08-26** et n'ont pas été ré-ancrés sur un commit. Un rapporteur qui
> reclone aujourd'hui ne retrouvera pas nécessairement ces lignes. **Avant toute
> citation dans le manuscrit, ré-ancrer sur un commit** (action A3.1 de la fiche
> ci-dessus, non faite).
>
> **5. Ce que cet avertissement NE remet PAS en cause** : les mesures. Tous les
> repères de non-régression cités en fin de section ont été exécutés et leurs
> valeurs sont celles imprimées. La réserve porte sur l'**attribution** et sur la
> **portée** de ce qu'on peut en conclure, pas sur les nombres.

Le solveur d'Imperial College est **public** :
[`ImperialCollegeLondon/solidity-solver-open`](https://github.com/ImperialCollegeLondon/solidity-solver-open)
(LGPL-3.0, C, 17 000 lignes, format `.Y3D` — la lignée Munjiza de la thèse de
Guo et des articles de Yang *et al.*). Les trois clés ci-dessous ne sont plus
déduites d'un article : elles sont **relevées dans leur source**, fichier et
ligne cités. Toutes sont opt-in, tous les défauts restent bit-identiques.

| clé | valeurs | défaut |
|---|---|---|
| `jointDeltaC` | `exact` \| `guo` \| **`solidity`** | `exact` |
| `jointFailRule` | `any` \| **`majority`** | `any` |
| `strainRateDIFArm` | `insertion` \| `envelope` \| **`continuous`** | `insertion` |
| `gcBirth` | `ramp` \| **`penalty`** | `ramp` |
| `gcBirthPenMin` / `gcBirthPenMax` | bornes du facteur | 0.01 / 3.0 |
| `strainRateFilter` | `exponential` \| **`none`** | `exponential` |

**`jointDeltaC = solidity`** — leur code ne fait pas ce qu'écrit la thèse.
`Y3Dfd.c` l. 1098-1099 (mode I) et 1125-1126 (mode II) :

```c
op = R2*el*dpeft/dpepe;                 /* ouverture au pic   <-> dnE */
ot = MAXIM((R2*op),(R3*dpegfn/dpeft));  /* PLAGE d adoucissement      */
```

et la rupture est à `op + ot`. La convention `guo` (δ_c = 3G_f/f_t depuis zéro,
son éq. 2.30) oublie **et** l'offset `op` **et** le plancher `2·op`. Ce
plancher ne mord que si 3G_f/f_t < 2·dnE, c'est-à-dire en maillage **fin**
devant G_f/f_t — le régime d'un impact, pas celui d'un essai de traction. D'où
le repère `jointdeltac_solidity_2d` à −2,333 % contre −2,330 % pour `guo` :
sur un maillage grossier les deux sont indiscernables, et c'est normal.
Leur propre commentaire porte deux fois `/*need further investigation*/`.

**`jointFailRule = majority`** — `Y3Dfd.c` l. 1175 : `if((nfail>1)&&...)`. Un
seul point d'intégration au-delà de z ≥ 1 **ne tue pas** la facette ; il en
faut deux (sur trois en 3D, sur deux en 2D). La règle n'a de sens qu'avec un
endommagement **par point** : chez eux `z` est une variable locale de la boucle
d'intégration, donc un point rompu cesse de transmettre pendant que les autres
tiennent. rockim ne portait qu'un scalaire `J.D` par joint, déjà le *max* des
points — la clé arme `Joint::Dk[]` et rend chaque point autonome. `J.D` reste
tenu à jour comme le max, pour toutes les sorties.
Exige `jointQuadrature = midedge` (les points comptés doivent être les leurs).

> ✅ **Deuxième source, d'article celle-là** (ajoutée le 2026-08-30). Cette
> convention n'est pas seulement lue dans le code : le manuscrit UCL
> (`Manuscript_UCL_deposit.pdf`, **p. 14**) écrit qu'une facette est déclarée
> rompue quand « *at least two integration points have zero stress components* ».
> `nfail > 1` dans le code et « at least two » dans le texte concordent. C'est la
> seule clé de cette section qui possède **deux sources indépendantes** ;
> les autres n'ont que le code.

**`strainRateDIFArm = continuous`** — leur DIF n'est jamais gelé. `dpeftdif`
est une variable **locale de la boucle élément**, reprise à chaque pas
(l. 1448-1456), et le même facteur multiplie la résistance **et** son énergie
de rupture (f_t avec G_I, c avec G_II). Cela lève le seul point où la
réplication était réputée impossible : il n'y a pas d'instant de gel à deviner.
Conséquence d'implémentation : le facteur ne peut pas s'appliquer *en place*
sur f_t/coh/G_f/G_fII sans se composer indéfiniment — `snapBase()` sauvegarde
les valeurs de base une fois, et `refreshDif()` reconstruit à chaque pas. Le
contrôle falsifiant est le repère `dif_continuous_2d` : `difmed = 1,53036`,
sous le plafond 1,85 de `yang-fig2`. Une composition donnerait un nombre
astronomique.
⚠️ Sous cet armement, `edotIns` du résumé porte le taux **final**, pas celui
d'un instant de gel : il ne se lit pas comme celui des deux autres armements
(d'où 7,66 /s contre 57,2 pour `dif_intrinseque_2d`, même essai).

**Ce que la source CONFIRME** (rockim était déjà juste, rien à changer) : la
loi de volume `T = (μ/J)B + [(λ lnJ − μ)/J]I + η·D` (l. 716) ; la viscosité
`dpeks*D`, donc `bulkViscosity = η/2` puisque rockim écrit 2μD ; la z-curve
a = 0,63 b = 1,8 c = 6 (l. 1088-1090) ; le couplage elliptique
`SQRT(tmp1²+tmp2²)` (l. 1136) ; la branche élastique parabolique (2r−r²)
(l. 1274) et la raideur **double** en compression (l. 1265) ; les 3 points aux
milieux d'arêtes avec le poids A/6 sur chacun des deux nœuds de l'arête
(l. 900-917 et 940) ; et l'**absence totale d'amortissement dans
l'intégrateur** (`Y3Dsd.c`), la viscosité de volume étant leur seule
dissipation hors joints et frottement.

**Ce que la source montre DÉSACTIVÉ chez eux** : le DIF lui-même
(`dpeftdif = R1`, soit 1,0) et l'endommagement diffus du volume (`df = R0`,
avec le `(1−df)` qui ne s'applique **pas** au terme visqueux) — ce dernier
recoupant ce que dit leur article de 2026 sur le granite de Kuru à propos du
calcaire et du grès.

**Frottement de contact** — `Y3Did.c` l. 1017-1051 : leur loi est
*structurellement identique* à celle de rockim, ressort tangentiel à
glissement **mémorisé** puis retour radial de Coulomb, et non un amortisseur
en vitesse. Le rapport est `ktss = 2.0/(7.0)*penalty`, soit **k_t/k_n = 2/7**.
Dans rockim ce rapport vaut `potTangentFactor / potPenaltyFactor` : c'est une
valeur de config, pas une capacité — aucun code n'a été touché. Leur bloc de
frottement statique/dynamique à affaiblissement en vitesse existe mais est
**commenté** ; la loi active est `mu = mud*d_fact`, Coulomb constant.

**`gcBirth = penalty`** — la naissance d'un contact sur un joint rompu,
`Y3Did.c` l. 915-964. Les deux modes sont des philosophies opposées :

- `ramp` (défaut, historique) : on retranche un relevé de naissance (volume en
  3D, aire en 2D) qui décroît en `exp(-t/gcBirthTau)`, si bien que la force
  **part de zéro** et remonte. Sous un indenteur, où les joints meurent *en
  compression* sous forte charge, c'est une perte de portance à chaque rupture.
  **Le problème ET une rampe sont publiés** (relevé le 2026-08-30) :
  `Manuscript_UCL_deposit.pdf` **p. 17** décrit exactement la difficulté — « *when
  a shear fracture under normal compression is formed, the overlap between
  tetrahedral elements due to compression will generate an initial non-zero
  contact force f_contact^initial, which can cause instability problems* » — et
  publie son remède, **éq. (18)** : `f_contact = (n_c/n_total)·f_contact^initial`,
  avec « *n_total is the total time-steps for n_c (usually 10)* ». Leur rampe est
  donc **linéaire sur ~10 pas**, là où `gcBirthTau` pose une décroissance
  **exponentielle** : les deux ne se comparent pas terme à terme (voir la réserve
  d'échelle plus bas).
- `penalty` : au pas exact de la naissance ils lisent la force que le joint
  portait en mourant (`d1ejfc*`, ce que rockim enregistre déjà dans
  `Joint::fDeath`) et calent la pénalité de **cette paire** pour que la force
  du contact naissant l'égale — `d1pepe[icoup] = penalty·fn_joint/fn_contact`,
  bornée à [0,01 ; 3]. La force est **continue**, le facteur persiste ensuite
  pour la paire, et la raideur tangentielle le suit
  (`ktss = 2/7·d1pepe[icoup]`). Ils zèrent aussi l'effort tangentiel et le
  glissement mémorisé au pas de naissance — reproduit.

Exige `contact = potential` (le mécanisme y vit ; sous `contact = penalty`, le
défaut, la clé serait inerte et le solveur refuse). Exclusive de `gcBirthTau`.

> ⚠️ **Réserve d'échelle sur `gcBirthTau` (2026-08-30, contre-audit).** Comparer
> `gcBirthTau = 1e-6 s` aux « ~10 pas » de leur éq. (18) demande deux précautions.
> (a) **Leur rampe est linéaire, la nôtre exponentielle** : `relax_ = exp(-dt/τ)`
> n'a pas de « longueur », seulement une constante de temps ; le rapport n'a de
> sens qu'à un facteur près. (b) Le nombre de pas dépend du `dt` du run, et un
> chiffre de **~50×** a circulé en interne : il est **faux**, il importait le `dt`
> d'un run de gradient St Anne (1,93e-9 s) dans une ligne qui parle d'impact 3D.
> Sur les seuls `dt` 3D mesurés et publiés par le dépôt — **1,30e-8 s** (insertion
> intrinsèque) et **1,93e-8 s** (adaptative),
> [`chantier_imperial_2026-08-29/A11_dt_tangentiel.md`](etat_de_l_art/chantier/A11_dt_tangentiel.md)
> §5-6 — cela fait **77 et 52 pas**, soit un facteur **5 à 8**, pas 50. (c) Enfin
> `gcBirthTau` est **inerte** sous `gcBirth = penalty` : `relax_` n'est lu que dans
> la branche `else` de la naissance, et le solveur refuse même de poser les deux
> clés ensemble. Un balayage de τ ne renseigne donc **que** le mode `ramp`.

⚠️ **Le relevé de naissance n'était pas qu'une douceur.** L'en-tête de
`PotHist::aRef` documente sa vraie raison d'être : il empêche une **injection
d'énergie** sur une paire née en recouvrement (mesure historique : **+936 J/m
sans relevé**, +179 avec une rampe purement temporelle — d'où l'asymétrie
d'état retenue). `gcBirth = penalty` le supprime : le bilan d'énergie devient
donc le contrôle obligatoire, et le solveur l'écrit lui-même — *« en mode
penalty tout positif est une injection »*. Mesure du 2026-08-26 sur la
percussion 2D : résidu **−0,9174 J/m** en `ramp` contre **−0,994861** en
`penalty` — négatif donc dissipatif dans les deux cas, et même légèrement plus
dissipatif. Aucune injection sur cet essai, mais **c'est à revérifier sur tout
nouveau cas** : le repère `gcbirth_penalty_percussion_2d` verrouille ce résidu
précisément pour ça.

⚠️ **Piège de lecture du repère.** En traction pure les joints meurent sans
charge à relayer (`fDeath = 0`) et le facteur retombe à 1 pour toutes les
paires : `gcbirth_penalty_2d` ne mesure alors que la *suppression de la rampe*
(−1,67766 → −1,85267 %). Le rééchelonnement lui-même n'est exercé que sous
indenteur — d'où `gcbirth_penalty_percussion_2d` : 4 joints morts, **100 % en
compression**, 651 kN/m relayés, facteur moyen 1,031, travail de contact
0,213 → 0,274 J/m dont frottement 0,111 → 0,116. Le résumé imprime le facteur
moyen : **s'il colle à une borne, c'est le clamp — arbitraire chez eux — qui
décide à la place de la physique**, et il faut l'élargir pour le savoir.

**`strainRateFilter = none`** — le taux qui alimente le DIF. `Y3Dfd.c` l. 1448
prend le taux de l'élément **tel quel**, sans lissage. rockim filtrait par un
passe-bas de constante `strainRateTau`, et sa propre raison était explicite :
« le taux brut par élément est trop bruité pour **figer** un DIF dessus ».
Cette raison vise le *gel*. Sous `strainRateDIFArm = continuous` le facteur est
repris à chaque pas — un pic de bruit ne dure qu'un pas au lieu d'être gravé
dans le joint — et l'argument tombe. **Les deux clés se répondent : c'est
ensemble qu'elles font leur schéma.** Exclusive de `strainRateTau`, qui serait
sans effet (refusée plutôt qu'ignorée en silence).
Mesures : 2D `edotmed` 7,65775 (filtré) → 4,02148 (brut) ; 3D 0,094298 →
0,0746516. L'écart est exactement ce que le passe-bas retirait.

**`jointResidualMu` est confirmé absent chez eux** : leur `dpefm` vaut `0.0`
en dur (l. 1091). Un joint rompu ne porte aucun cisaillement ; tout le
frottement vient du contact. Pour une réplication fidèle, laisser la clé
désactivée.

Repères — tier fast : `jointdeltac_solidity_2d`, `jointfailrule_majority_2d`,
`dif_continuous_2d`, `gcbirth_ramp_2d`, `gcbirth_penalty_2d`,
`srfilter_none_2d` et leurs contrôles à charge nulle. Tier full :
`jointdeltac_solidity_3d`, `jointfailrule_majority_3d`, `dif_continuous_3d`,
`gcbirth_penalty_3d`, `srfilter_none_3d`, `zeroload_gcbirth_penalty_3d`, et
`gcbirth_penalty_percussion_2d` — le seul qui exerce vraiment le
rééchelonnement.

⚠️ `gcbirth_penalty_3d` ne verrouille **pas** un écart : en traction 3D les
deux modes donnent le même `err_pct` (−4,75889 %). Il verrouille le *fait* que
le mécanisme s'arme (118 paires calées). C'est délibéré et noté dans le repère.

### 5.5 Lois de comportement (`law`, modes fem3d / fdem / fdem3d)

`law = elastic | dpr | saksala | saksala2011 | dpdfh` (défaut : dpr en fem3d ; absent
en fdem/fdem3d = bulk élastique + crushCap, bit-compatible avec l'historique).
En fdem 2D la loi 3D est utilisée en déformation plane exacte (ε_zz = 0). `law` est
incompatible avec `phases` (mono-matériau).

| loi | clés propres (défauts) |
|---|---|
| `dpr` | Drucker-Prager calé MC + Rankine crack-band : `erodeD` (0.98), `erodeEpv` (1.5) |
| `mc` | **Mohr-Coulomb élasto-plastique de Ye et al. (IJRMMS 194, 2025) — la loi « MC-FDEM »** : vrai critère à arêtes en contraintes principales (pas l'approximation lisse de `dpr`), retour de Clausen à 4 régions (plan principal, arête de compression, arête d'extension, apex), écoulement **non associé** par la dilatance ψ. Toute la fissuration reste dans les joints, toute la dissipation plastique dans le bulk. Clés : `mcCohesion` (défaut = `cohesion`), `mcFrictionDeg` (défaut = `frictionDeg`), `mcDilationDeg` (défaut 0 ; ψ = φ = associé). Vérifiée analytiquement par `rockim selftest-mc` : σc = 2c·cos φ/(1−sin φ), σt = 2c·cos φ/(1+sin φ), σ₁ = N·σ₃ − σc — **écart 1e-9 % en compression et triaxial**, 0,13 % en traction (résolution du pas). Repère de suite `selftest_mc` |
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
(sides \| all \| bore en 2D ; lateral \| all en 3D), `confineGaugeTime` (3×rampe ; instant de
la jauge σ_latéral atteint). Pression SUIVEUSE sur les faces extérieures d'origine
seulement (pas dans les fissures). tension + pullV < 0 + confinement = triaxial.

**Cavité pressurisée (`confineFaces = bore`, 2D, 2026-08-14)** — pressurise les
SEULES faces extérieures d'origine dont le milieu est à moins de `boreSelectR`
de (`boreCX`, `boreCY`) (défauts : centre du bloc) : tunnel/forage sous pression
→ fissures radiales. Maillage troué : `make_unstructured_mesh.py tunnel W H R h`.
Montage type : `scenario = tension` + `pullV = 0` (mors immobiles = plaque tenue).
Cas de référence `configs/tunnel_bore.cfg` (granite, R = 10 mm, 40 MPa) :
rampe 100 µs → **53 joints rompus** (fissure diamétrale, amorçage ~30 MPa) ;
rampe 10 µs → **144 joints** (étoile radiale) — N croît avec ṗ, signature
d'obscuration du banc 6 Abaqus DP-DFH reproduite. ⚠️ Limite connue : la jauge
« achieved σ_xx » lit le cœur du bloc, sans signification en mode bore.

**Bilan B4 et confinement (corrigé 2026-08-14 soir)** — le travail de la
pression suiveuse est désormais comptabilisé (`confWork_`, ligne
« confinement » du résumé, imprimée seulement si `confiningPressure > 0` —
sorties des runs non confinés inchangées au bit). Les résidus [CHECK] à
~100 % des runs confinés/bore antérieurs à cette date étaient CE poste
manquant, pas une injection d'énergie.

**Arrêt post-rupture (`stopPeakDrop`, fdem3d tension, 2026-08-14)** —
0 = off (défaut) ; sinon, l'essai s'arrête PROPREMENT (hook `finished()` :
frame + history + summary) dès que la contrainte des mors retombe sous
(1 − stopPeakDrop) × pic, avec garde anti-bruit (pic > 1 MPa). Motivation :
le post-pic profond met toute la bande de cisaillement en contact de lèvres
— le prix D0 par pas — sans rien apporter à la mesure du pic. Décision
Fernando : « je n'ai pas besoin d'aller loin après la rupture ».

**Moniteur d'énergie runtime (E2 fiabilité, 2026-08-14)** —
`budgetAbortPct` (0 = off, défaut) : tous les 1024 pas, si le résidu B4
courant dépasse ce pourcentage de l'échelle (même définition que le résumé),
**arrêt PROPRE** via le hook `finished()` : dernière frame, dernière ligne
d'history et summary sont écrits, avec l'empreinte du hotspot (nœud le plus
rapide). Un run qui diverge laisse son autopsie au lieu de mégajoules de
débris. Typique : `budgetAbortPct = 5` en production, off pour les études
de diagnostic (E0) où l'on VEUT voir la divergence se développer.

> ⚠️ **À NE PAS ARMER SANS `energyBodyForces = on`** (mesuré le 2026-08-30).
> Le résidu que ce moniteur juge **ne contenait pas** le travail des forces
> volumiques : ni la pesanteur (aucun compteur n'existait) ni le tri des
> fragments (`brushWork_`, tenu hors bilan à dessein). Résultat mesuré sur
> `configs/fdem3d_percussion.cfg` + `gravity = 9.81` + `budgetAbortPct = 2` :
> **le moniteur ABORTE à t = 8,344·10⁻⁶ s** sur un résidu de **175 % de
> l'échelle** qui vaut, à 0,3 % près, le seul travail de la pesanteur — sur un
> run à **zéro joint rompu**. Avec `energyBodyForces = on`, même run, **aucun
> déclenchement**, résidu 0,25 %, `[OK]`. Repères `ebody_abort_defaut_3d` et
> `ebody_abort_on_3d` (tier `all`).

**Forces volumiques dans le bilan (`energyBodyForces`, 2026-08-30, 2D et 3D)** —
`off` (défaut) \| `on`. **La MESURE est inconditionnelle** : le travail de la
pesanteur (`gravWork_`, compteur créé à cette date — c'est le **septième poste**
du bilan d'ARMA 24-0952 éq. 3-7, le seul que rockim n'avait pas) et celui du tri
des fragments sont toujours calculés et **imprimés** au résumé, ligne
`forces vol.` ; ils ne touchent aucune force, donc la physique est
bit-identique dans les deux réglages. **Seule leur entrée dans `sumW` est
opt-in.** `off` : ils **tombent dans le résidu**, et le résumé le dit en toutes
lettres. `on` : ils entrent dans `sumW` **et dans l'échelle**, donc dans le
verdict `[OK|CHECK]` **et dans `budgetAbortPct`** — bilan à sept postes.
Mesures : résidu 1,05718e-06 → −2,015e-10 J (facteur **5 250**) sur la percussion
3D à T = 2e-5 ; 1,71064e-10 → −1,645e-13 J (facteur **1 040**) à T = 2e-6, où le
défaut affiche `[CHECK]` **à 116 % de l'échelle sur un run parfaitement sain**.
Repères `ebody_defaut_3d` / `ebody_on_3d` (tier `full`). Fiche :
[`chantier_imperial_2026-08-29/B10_bilan_energie_forces_volumiques.md`](etat_de_l_art/chantier/B10_bilan_energie_forces_volumiques.md).

> ⚠️ **Clé muette enregistrée en chemin** : en **3D**, `gravity` n'est lu que
> dans `placeTool()`, dont la première ligne est
> `if (scen_ == Scenario::TENSION) return;`. **Un deck de traction 3D qui pose
> `gravity` ne reçoit aucune pesanteur.** Le comportement n'a pas été changé —
> le corriger changerait la physique d'un deck existant — mais il n'est plus
> silencieux : un `AVERTISSEMENT` explicite est imprimé. Le solveur **2D**, lui,
> lit `gravity` dans `init()`, sans garde de scénario : **rupture de parité
> 2D/3D**, la sixième du registre.

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

### 5.10 Couplage hydro-mécanique (spec 004, 2D `fdem`, ajouté le 2026-08-19)

Modèle d'AbuAisha, Eaton, Priest & Wong (JPSE 154, 2017), le module HF du code
Y-Geo. Hypothèse centrale, qu'ils énoncent en toutes lettres : **la pression du
fluide est UNIFORME dans la cavité et les fissures** — pas de loi cubique, pas
de gradient, pas de leak-off, pas de pas de temps hydraulique propre. C'est un
fluide non visqueux, valable tant qu'on regarde le voisinage du puits.

Ce que ça apporte, et que `confineFaces = bore` ne savait pas faire — son
commentaire l'avouait, *« faces born from cracking receive nothing »* :
**la pression SUIT la fissure**.

| clé (défaut) | rôle |
|---|---|
| `hydro` (false) | active le couplage. Refusé hors des scénarios qui ont une cavité |
| `hydroSource` (bore) | `bore` = les faces extérieures dans `boreSelectR` autour de (`boreCX`, `boreCY`) ; `all` = toute la frontière extérieure |
| `hydroInjection` (rate) | `rate` = pompe à débit, la pression est une SORTIE ; `pressure` = pression imposée, pour les contrôles |
| `hydroRate` (0) | débit [m³/s par mètre d'épaisseur]. 20 l/s de l'article = `0.02` |
| `hydroPressure` (0) | pression imposée [Pa], mode `pressure` seulement |
| `hydroP0` (0) | pression de référence [Pa]. À 0 on travaille en pressions EFFECTIVES |
| `hydroRamp` (0) | rampe cosinus de la pompe [s]. Même forme que `confiningRamp` |
| `fluidBulk` (2.2e9) | K_f [Pa]. ⚠️ **l'article ne le donne pas** — 2,2 GPa est l'eau, c'est une hypothèse, et elle fixe toute la chronologie (§8) |
| `fluidDensity` (1000) | ρ_f0 [kg/m³] |

**Trois grandeurs, une seule variable d'état.** La *frontière mouillée* : les
faces reliées à la source par un chemin de joints rompus (recherche de
composante connexe sur les sommets, leur module 2). Le *volume* de cavité, par
décomposition locale — lacet de Green sur les seules faces source (contour fermé
par construction) plus, pour chaque fissure mouillée, son aire propre
L·(ouverture moyenne), forme employée par Lisjak et al. 2017. La *pression*, par
compressibilité linéaire `p = p0 + K_f·log(m / (V·ρ_f0))`, où la masse injectée
`m` — remplissage initial COMPRIS — est la variable d'état.

**Le chargement.** La pression pousse le solide à l'opposé de la normale
sortante, `−p·L·n/2` par nœud, en force suiveuse : exactement la forme de
`confiningForces()`.

> ⚠️ **Correctif du 2026-08-20 — le signe était inversé.** `hydroForces()`
> appliquait `+p n` : le fluide serrait la cavité au lieu de l'ouvrir, le forage
> produisait un breakout aligné sur σ′_h et rompait à 6,6 MPa au lieu de 12.
> La cause est une mauvaise lecture de leur éq. 7, `F = −(p/2)[y₂−y₁ ; x₂−x₁]`,
> dont le vecteur n'est pas orthogonal au segment : **la coquille est dans la
> seconde composante seule, pas dans le signe de tête**, et c'est ce signe de
> tête, le bon, qui avait été supprimé. Tranché par leur §3.1 (« negative sign
> to compressive stresses ») et par le maillage CCW de Y-Geo.
>
> **Le contrôle censé l'attraper existait et n'a rien vu** : H3
> (`parker_compare.py`) mesurait `max(y) − min(y)`, une valeur ABSOLUE, et a
> validé une interpénétration comme une ouverture. *Un contrôle de signe qui
> passe par une norme ne contrôle pas le signe.* Contrôle de remplacement :
> `bench_abuaisha/tools/hydro_sign_check.py` — même pression par les deux
> chemins de chargement, même déplacement de paroi POSITIF (mesuré : écart
> 0,000 % entre chemins, −2,0 % de Lamé).

### 5.11 Discontinuités préexistantes et catalogue microsismique (chantier f2, 2026-09-01)

Deux capacités ajoutées pour la §2.4 et les éq. 11-13 d'AbuAisha et al. 2017.
Toutes deux **opt-in strictes** : sans la clé, aucun chemin ne change, ce qui
est vérifié par diff **octet à octet** des sorties (28 fichiers, 0 différence)
entre le binaire d'avant et celui d'après sur `verify_fdem_voronoi_tension`.

| clé (défaut) | rôle |
|---|---|
| **`preBrokenJoints`** (—) | **sélecteur géométrique de fissures préexistantes** : `"x1 y1 x2 y2; x1 y1 x2 y2; ..."` en coordonnées de configuration (m). Tout joint dont le **milieu d'arête** tombe à moins de `preBrokenTol` d'un segment **et** dont l'orientation s'en écarte de moins de `preBrokenAngleDeg` naît rompu. Se compose avec `jointPrebrokenFrac` (union des deux populations) |
| `preBrokenTol` (−1) | tolérance de distance [m]. **Négatif = demi-longueur de l'arête courante**, ce qui suit la gradation du maillage — indispensable sur un maillage 3 mm en paroi / 300 mm au loin, où une tolérance uniforme ne sélectionnerait rien près du trou ou une bande large au loin |
| `preBrokenAngleDeg` (30) | écart d'orientation toléré entre l'arête et le segment. Sans ce filtre on ramasse les arêtes **transverses**, et la « discontinuité » devient un escalier de joints perpendiculaires qui scie la roche au lieu de la fendre |
| **`dampingViscous`** (0) | **amortissement nodal $-\mu v$** de leur éq. 9, $\mathbf{C} = \mu\mathbf{I}$ — la même force sur chaque nœud, indépendamment de sa masse. En kg/(m·s), multiplié par l'épaisseur pour rendre des newtons. Distinct de `dampingLocal` (Cundall sur $\lvert f\rvert$, **sans dimension**) et de `bulkViscosity` (Kelvin-Voigt de volume). Appliqué **une fois par groupe de liaison** ; son travail est compté avec celui de Cundall |
| **`dampingViscousScheme`** (explicit) | `explicit` = la forme de leur éq. 8, la force entre dans le résidu ; **borne $dt \le 2m/\mu e$**. `implicit` = division $v \mathbin{/}{=} 1 + dt\,c/m$, la forme des amortisseurs de Lysmer du même intégrateur : **inconditionnellement stable, aucune borne**. Refusé si `dampingViscous` n'est pas posée |
| **`hydroCavityClosure`** (false) | **referme le lacet de Green** des faces source en pontant les bouches de fissure débouchant en paroi. Sans lui l'éventail omet un coin d'aire $\approx aR/2$ par bouche : $V$ sous-estimé, donc $p$ **sur-estimée** |
| **`microseismic`** (false) | **catalogue microsismique**, éq. 11-13. Instrumente chaque joint : `tYield` (entrée en endommagement), `keYield` (énergie cinétique de ses **quatre** copies à cet instant) et `dKeMax` (maximum de l'accroissement jusqu'à la rupture). Écrit `fdem_seismic.csv` en fin de run |

**Ce que `preBrokenJoints` pose, et rien d'autre** : `pre = true`,
`bonded = false`, `D = 1`, `tBreak = 0`, `tInsert = 0`, `bmode = 4`. Il ne
touche **pas** `dead` : la fissure garde sa loi de joint, donc son **frottement
résiduel** `jointResidualMu` — c'est littéralement leur §2.4 (plan sans
cohésion, cisaillement par frottement pur). Elle est aussi **mouillée
naturellement** dès qu'une fissure hydraulique l'intersecte, et compte dans le
volume de cavité : le front mouillé ne teste que `!bonded && D >= wetDmin_`,
rien d'autre. À `t = 0` les copies sont co-localisées, donc sa contribution au
volume est **exactement nulle** et `hydroVol0_` est inchangé.

Trois garde-fous **refusent le run** plutôt que de produire un résultat
faux-mais-plausible : `jointResidualMu` non posé, `jointContactPenalty =
adaptive` (la raideur $(1-D)p_j$ s'annule à $D = 1$), `jointDeath = damage`
(toutes les pré-fissures mourraient au premier pas). Un quatrième, propre au
sélecteur : **un segment qui ne sélectionne aucun joint est une erreur**, pas
un avertissement. Et un avertissement bruyant si aucune pré-fissure n'a
d'extrémité scindée — la discontinuité serait alors cinématiquement **inerte**.

**Colonnes de `fdem_seismic.csv`** :
`jointId, x, y, tYield, keYield, dKeMax, tBreak, breakMode, magnitude`, une
ligne par joint **entré en endommagement**. `magnitude` applique leur éq. 13,
$M = \frac{2}{3}(\log_{10} E - 4{,}8)$, et vaut `nan` si `dKeMax` est nul.
L'article ne compte comme événement que la **rupture** : filtrer sur
`tBreak >= 0`. Les joints endommagés non rompus sont fournis parce qu'ils ne
coûtent rien et documentent l'incubation.

> ⚠️ **Pourquoi ceci est dans le solveur et pas dans un script.** Le maximum de
> l'éq. 13 se prend **pas à pas**. Mesure sur le banc AbuAisha : entre deux
> trames consécutives, **48 des 49 joints qui cassent** naissent, s'endommagent
> et rompent dans le **même intervalle**. Aucun post-traitement de trames ne
> peut donc reconstruire ni $E_{k,y}$ ni le maximum — il faudrait une trame par
> microseconde, soit ~350 Go par run, et le maximum resterait échantillonné.
>
> ⚠️ **Unités.** Le solveur 2D porte une tranche d'épaisseur `thickness` et
> `m_` vaut $\rho A\,e$ : l'énergie est donc en **joules pour cette tranche**,
> soit des joules par mètre de forage à la convention `thickness = 1`. La
> relation de Gutenberg attend des joules ; la magnitude hérite de cette
> convention plane et ne se compare à une magnitude de terrain qu'à ce titre.
> C'est une limite de la représentation 2D, pas du calcul.
>
> ⚠️ Les **pré-fissures sont exclues** du catalogue : nées à $D = 1$, elles
> n'ont jamais cédé et émettraient autant d'événements fantômes à $t_y = 0$
> portant l'énergie cinétique initiale du bloc.

### L'amortissement de leur équation 9 : deux lectures, et ce que la mesure en dit

Leur éq. 9 pose $\mathbf{C} = \mu\mathbf{I}$ et leur Table 1 donne
$\mu = 5{,}6\cdot10^5$ **kg/(m·s)**. Les deux ne s'accordent pas : ces unités
sont celles d'une **viscosité dynamique** — $\mu D$ rend bien des Pa — et non
celles d'un amortisseur nodal, qui serait en kg/s. D'où deux lectures :

- **(A) lecture des unités** — leur $\mu$ est la viscosité de Munjiza, et
  `bulkViscosity = 5.6e5` reproduit leur Table 1 **sans une ligne de code**.
  Obstacle : la borne diffusive $\rho h^2/4\mu$ écrase le pas de temps.
- **(B) lecture littérale de l'éq. 9** — c'est `dampingViscous`.

> **La mesure tranche, et contre (B).** La constante de temps d'un amortisseur
> nodal vaut $\tau = m/(\mu e)$. Sur le maillage du banc AbuAisha :
> $\tau = 1{,}7\cdot10^{-8}$ s pour un élément de 3 mm et
> $1{,}8\cdot10^{-9}$ s pour le plus fin — **contre un pas de temps de
> $2{,}1\cdot10^{-8}$ s**. À leur valeur, un nœud ne serait pas amorti, il
> serait **figé en un pas**. Mesuré aussi sur l'essai `t7` : le pas tombe d'un
> facteur 6,5 (la garde explicite $2m/\mu e$ prend la main) et le pic d'un
> essai de traction piloté en vitesse passe de 10,8 à **26,0 MPa** — parce
> qu'un $-\mu v$ taxe le mouvement imposé lui-même, ce qu'une viscosité de
> volume, aveugle à une translation uniforme, ne fait pas.
>
> `dampingViscous` reste utile à des valeurs **plus petites**, et c'est le seul
> moyen d'éprouver la forme littérale. Mais **la valeur de leur Table 1 n'est
> pas utilisable sous cette lecture** — ce qui est un argument de plus pour
> lire leur $\mu$ comme une viscosité de volume.

**Le schéma implicite : même réponse, 6,5 fois moins cher.** Mesuré sur `t8`,
à $\mu$ identique :

| | pas de temps | pas | pic |
|---|---|---|---|
| `explicit` | $2{,}12\cdot10^{-9}$ s | 283 415 | 25,982 MPa |
| `implicit` | $1{,}37\cdot10^{-8}$ s | **43 652** | 25,982 MPa |

Le pas redevient celui d'un run **sans amortissement du tout** : la borne
$2m/\mu e$ disparaît entièrement. Et les deux schémas donnent **le même pic à
$10^{-4}$ près**.

> Ce second point est le plus instructif : la montée du pic de 10,8 à 26 MPa
> n'est donc **pas** un artefact du traitement explicite. C'est le modèle. Un
> $-\mu v$ nodal taxe le mouvement imposé lui-même, et les deux schémas
> l'attestent également. L'implicite achète du pas de temps, **pas de la
> physique** — la constante de temps $\tau = m/(\mu e)$ reste une propriété du
> modèle, et un nœud dont $\tau$ est plus court que le pas reste figé, mais
> stablement.

### Le coin manquant du volume de cavité

Le volume de la cavité source est un **éventail** de triangles
(centroïde, $P_i$, $Q_i$) sur les faces du forage. Tant que l'anneau est fermé
l'éventail *est* l'aire du polygone. Dès qu'une fissure débouche en paroi, deux
faces voisines s'écartent d'une bouche $a$ et l'éventail perd le coin
(centroïde, $Q_i$, $P_{i+1}$), d'aire $\approx aR/2$.

Chiffré sur le deck du banc ($a = 89$ µm, $R = 50$ mm) : environ **0,6 MPa de
pression en trop, par bouche ouverte**. Au pic l'ouverture n'est que de
quelques micromètres et le biais est négligeable — **c'est un biais de
post-pic**, et il va dans le sens de la surpression résiduelle observée.

`hydroCavityClosure = true` referme le lacet en pontant les bouches (tri des
faces par angle polaire autour du centroïde, sens de parcours choisi par la
longueur des intervalles, et rejet signalé de tout intervalle de plus de trois
longueurs de face). Vérifié sur `t6` : écart **exactement nul** à contour
fermé, puis $9{,}4\cdot10^{-7}$ m²/m dès que les fissures débouchent, à
mécanique **rigoureusement identique** (même $t$, même `nBroken` à chaque
ligne — l'essai impose une pression nulle pour isoler le volume).

> ⚠️ Le défaut `false` est délibéré : il garde la comparabilité avec les sept
> calculs archivés du banc. La valeur `true` est la **plus juste des deux** —
> à poser pour toute étude qui lit le post-pic. Valide pour une cavité
> **étoilée** par rapport à son centroïde : vrai d'un forage circulaire.

**Contrôle de bon fonctionnement** (`tests_f2/`, ~20 s chacun) : `t1_seismic`
laisse le pic et le compte de ruptures **inchangés** (11,156 MPa, 15 joints)
tout en écrivant 133 événements — l'instrumentation ne perturbe pas la
mécanique ; `t2_prebroken` pose une entaille de 30 mm à mi-hauteur d'une
éprouvette de 60 mm et fait tomber le pic à **4,86 MPa**, sous la moitié du
témoin, comme l'exige la réduction de ligament plus la concentration en
pointe ; `t3` et `t4` vérifient que les garde-fous refusent bien le run.

### 5.12 Choc thermique de paroi (chantier f2, 2026-09-01)

Conduction transitoire **explicite sur le graphe des joints** — chaque joint
relie deux éléments, le maillage FDEM *est* le graphe de conduction — avec
condition de **Robin** (coefficient d'échange fini : la revue établit que
l'ébullition en film plafonne le flux, facteur 2,5-3 sous le Dirichlet idéal).
La contrainte thermique $-3K\alpha_T(T-T_{ref})\,\mathbf{I}$ entre au même site
que l'in situ : elle nourrit le VTU, la jauge de confinement et **le critère
d'insertion adaptative** — c'est par là que le froid casse. 2D `fdem`
seulement (gardé dans `main.cpp` : un deck 3D est **refusé**, pas ignoré).

| clé (défaut) | rôle |
|---|---|
| **`thermal`** (false) | arme le couplage. Refuse `law`, `bulkDamage`, brazilian/shpb |
| `thermalTemp` (**requis**) | température du fluide (le forçage — pas de défaut) |
| `thermalH` (**requis** > 0) | coefficient d'échange [W/m²K]. Dirichlet = h très grand, en le disant |
| `thermalFaces` (bore) | bore \| all \| top — sélection des faces refroidies |
| `thermalConduct` (3) / `thermalHeatCap` (800) / `thermalAlpha` (8e-6) / `thermalTref` (0) | propriétés |
| `thermalStart` (0) | la pompe à froid démarre ici |
| `thermalSpeedup` (1) / `thermalEvery` (1) | **accélération d'horloge thermique** : la mécanique sert de relaxation quasi-statique entre les incréments. Écrêté à la borne de stabilité $0{,}5\,\rho c V/\sum G$, en le disant. **Validation prescrite : diviser par deux ne doit rien changer** |
| `thermalCrackResist` (1) | conductance d'un joint **rompu** (1 = inchangé, 0 = isolant) — la rétroaction fissure→conduction de Yan & Jiao 2020, forme minimale |
| **`toolStart`** (0) | percussion : l'outil est **suspendu** jusque-là — le phasage choc thermique → percussion dans le même run |
| `dampingSwitchT` (−1) / `dampingLocalAfter` | bascule d'amortissement en cours de run (relaxation forte pendant le froid, valeur de production pour l'impact) |

Sorties : champ `temp` par élément dans les VTU, colonnes `thermTw` (T moyenne
de paroi) et `thermQ` (**énergie déposée** — la moitié manquante du bilan que
la littérature omet) dans l'historique, et un bilan de conservation au résumé.

**Validation (tests_f2/, faite le 2026-09-01) :**

| contrôle | résultat |
|---|---|
| conservation discrète $Q = \sum\rho cV\,\Delta T$ | **exacte, 1,4·10⁻¹³** (identité, pas une tolérance) |
| forme fermée $\sigma_{yy} = E\alpha_T\Delta T/(1-\nu)$, barre bloquée refroidie (`t9b`, adaptatif) | 5,33275 contre 5,3333 MPa : **0,011 %** |
| même essai en intrinsèque (`t9`) | −4,8 % — c'est la **complaisance de pénalité** documentée (~4-5 % sur E au facteur 20), mesurée ici par une voie indépendante |
| seuil : $\Delta T = 10$ K → σ < $f_t$ | **0 fissure** |
| seuil : $\Delta T = 40$ K → σ = 21 MPa ≫ $f_t$ | **21 joints rompus**, pic bloqué à ~$f_t$ |
| bit-identité, clé absente | 28 fichiers, 0 octet |
| gardes (clé orpheline, `law`, mode 3D) | trois refus propres, code 1 |

> ⚠️ Ce que `thermalSpeedup` n'achète pas : la validité quasi-statique. Le
> protocole est de vérifier que l'énergie cinétique reste négligeable devant
> l'énergie déposée, et que halver la vitesse ne change rien. Et la limite
> physique du modèle reste écrite : **un seul $\alpha_T$** — le désaccord
> quartz/feldspath, moteur réel de la microfissuration granulaire, demanderait
> un `phase.<nom>.thermalAlpha` (extension GBM naturelle, non écrite).

### 5.13 Schistosite : famille de plans paralleles (chantier f2, 2026-09-01)

`weakPlanes` pose une **famille de plans paralleles**, de pendage et
d'espacement donnes, dont les joints portent des proprietes **reduites**.
C'est le sens strict du mot schistosite, a distinguer du champ correle
anisotrope (`strengthCorrAngleDeg`) qui ne produit que des **taches**
allongees : la, un joint est faible parce qu'il est au mauvais *endroit* ;
ici parce qu'il est dans la bonne *direction* **et** sur un plan.

Selection a deux conditions, comme `preBrokenJoints` : distance du milieu
d'arete au plan le plus proche **et** ecart d'orientation. Sans la seconde on
ramasse les aretes transverses et le plan devient un escalier qui scie la
roche. Le facteur **multiplie** ce que le tirage de Weibull a deja pose :
texture et direction se composent.

| cle (defaut) | role |
|---|---|
| **`weakPlanes`** | `<pendage_deg> <espacement_m>` — arme la capacite |
| **`weakPlaneFactor`** (**requis**) | facteur sur `ft` et `cohesion` des plans |
| `weakPlaneAngleDeg` (30) | tolerance d'orientation des aretes retenues |
| `weakPlaneTol` (< 0 = demi-arete) | demi-largeur de la bande de selection ; le defaut suit la gradation du maillage |
| `weakPlaneOffset` (0) | decalage de la famille le long de sa normale |
| `weakPlaneFrictionDeg` (< 0 = inchange) | frottement impose sur les plans |
| `weakPlaneGf` (**follow**) | `follow` : `Gf` suit le facteur, donc **l_cz inchangee**. `keep` : `ft` baisse seule et l_cz explose en 1/f^2 — le plan cesse d'etre une entite distincte. Defaut inverse de `weibullScope`, deliberement |

Sorties : champ `weakPlane` par joint dans `fdem_joints_*.vtu` (arme seulement
quand la capacite l'est) — sans lui la famille n'est visible **nulle part**,
`ftScale` melangeant Weibull et schistosite.

**`CONTINUITE` : le chiffre a lire avant toute interpretation.** Rapport de la
longueur reellement affaiblie a la longueur ideale des plans (aire/espacement).
Sous 0,5, les plans sont des **troncons epars** et le run n'est pas comparable
aux autres pendages. Le solveur avertit.

**Validation (2026-09-01) :**

| controle | resultat |
|---|---|
| anisotropie de resistance, eprouvette de traction | pic **11,16 -> 4,37 MPa** (facteur **2,6**) ; minimum quand les plans barrent la traction, **retour a l'intact a 90 deg** (plans paralleles a la traction) — la courbe en U de Jaeger |
| le diagnostic `CONTINUITE` fait son office | 2 pendages sur 7 rejetes a 0,25 et 0,30 sur la petite eprouvette : la comparaison y melangeait orientation et **densite** |
| decks tunnel (100x100 m, 159 269 joints) | continuite **0,839 / 0,860 / 0,857** et **8740 / 8782 / 8750** joints affaiblis a 0/45/90 deg — ecart de densite **0,5 %**, les trois pendages sont strictement comparables |
| bit-identite, cle absente | 28 fichiers, 0 octet |
| garde : cle satellite orpheline / aucun joint selectionne | refus propres |

> ⚠️ **Reserves.** Le plan discret suit les aretes du maillage, donc en zigzag
> — meme approximation que pour toute fissure en FDEM, mais elle rend le plan
> un peu plus resistant qu'un plan lisse : l'effet mesure est une **borne
> inferieure**. Et la roche reste **elastiquement isotrope** : seule la
> resistance est anisotrope, pas E.

### 5.14 Schistosité pervasive de Lisjak : les trois briques (chantier f2, 2026-09-01/02)

Méthode de référence en FDEM pour une roche litée à l'échelle de l'ouvrage —
Lisjak, thèse Toronto 2013, ch. 5 ; Lisjak, Grasselli & Vietor 2014, *IJRMMS*
65:96-115 ; cas réel Lisjak et al. 2015, *TUST* 45:227-248. Elle remplace, pour
une schistosité **pervasive** (espacement sous-maille), l'approche à plans
discrets `weakPlanes` (§5.13), que son auteur a lui-même déclarée *« unsuitable
for field-scale models »*. Trois briques ; **aucune ne suffit seule**.

**Brique 1 — élasticité transversalement isotrope (sur le triangle).**
Cinq constantes : `E`, `nu` du deck = dans le plan du litage ; `beddingEperp`,
`beddingNuPerp`, `beddingGperp` = E', ν', G' ; direction `beddingDip`.
Complaisance plane condensée (déformation plane) dans le repère du litage :
S11 = (1−ν²)/E, S12 = −ν'(1+ν)/E', S22 = (1−Eν'²/E')/E', S66 = 1/G' ;
σ_zz = ν σ₁₁ + (Eν'/E') σ₂₂. Rotation faite *numériquement* (loi appliquée
aux trois déformations unité) — aucune convention de signe à se tromper. Le
triangle est co-rotationnel : le litage tourne avec l'élément. Les joints ne
portent **aucune** anisotropie élastique (verbatim de la source). Le budget de
pas de temps lit la vitesse d'onde maximale par balayage du tenseur de
Christoffel. Refuse `phases`, `law`, `neohookean` et **`thermal`** (la
contrainte thermique est écrite isotrope, −3Kα_TΔT I ; sous TI elle devrait
devenir −D:(α_TΔT 1) — non dérivé).

**Brique 2 — loi cohésive directionnelle (sur le joint).**
X(γ) = X_min + (X_max − X_min)·γ/90°, X ∈ {ft, c, G_Ic, G_IIc}, γ = angle
joint/litage. **Minimum à γ = 0** (joint parallèle au litage : la délamination),
maximum à γ = 90. φ constant. Le deck porte les **maxima** ; les clés
`beddingFtRatio`, `beddingCohRatio`, `beddingGfIRatio`, `beddingGfIIRatio`
= X_min/X_max ∈ ]0 ; 1]. Se compose avec Weibull, taille, `weakPlanes`.

> ⚠️ **Coquille de la thèse.** Le texte (p. 102) écrit *« maximum and minimum
> values, for γ = 0° and γ = 90° »* — c'est l'**inverse** de sa propre
> fig. 5.8a (droite de Min à γ = 0 vers Max à γ = 90, vignettes à l'appui) et
> de sa table 5.1. Vérifié sur le PDF le 2026-09-01. Implémenter la phrase
> donnerait une roche impossible à déliter.

Valeurs calibrées Opalinus (table 5.1) : E = 3,8 / E' = 1,3 GPa, ν = 0,35 /
ν' = 0,25, G' = 0,9 GPa ; ft 0,16→0,65 MPa (ratio 0,246), c 1→9 MPa (0,111),
G_Ic 0,4→7,0 J/m² (0,057), G_IIc 10→35 J/m² (0,286), φ = 22°.

**Brique 3 — maillage à arêtes alignées sur le litage** (hors solveur) :
`tunnel_schisto/make_tunnel_bedded_mesh.py W H hFine rFine hFar dip t rBed out`.
Cordes parallèles d'espacement t dans le disque r < rBed, **découpées en
Python** (bissection sur `isInside`, bande d'exclusion de 0,8 h le long de
la paroi, retrait 0,5 h) puis **plongées** dans la face roche par
`mesh.embed`. Le fragment OCC ne fait pas ce travail (il découpe sans
attacher), et un retrait fixe laisse des slivers là où une corde longe la
paroi — six versions ont été nécessaires, l'historique est dans le script.
Règle de Lisjak h ≈ t/3. Le lecteur de rockim jette les tags : aucun grain
parasite. Contrôle imprimé : **continuité** = longueur d'arêtes exactement
alignées (< 1°) / longueur des cordes dans la roche.

> ⚠️ **Le pas de temps est fixé par le plus petit triangle — et par lui seul.**
> En `mesh = file`, `pj = 4E/hmin` avec `hmin` **global** : un sliver de
> 20 mm gonfle la pénalité de *tous* les joints (×3,6) en plus de réduire la
> masse du pire nœud. Mesuré : dt divisé par 2,9 pour 17 triangles sur
> 304 000. Cible de qualité : diamètre inscrit minimal ≥ celui du maillage
> isotrope de référence (73 mm à h = 0,2).

**Deux réglages livrés (2026-09-02), trois pendages chacun :**

| réglage | t / h | triangles | d_inscrit min | continuité | dt | coût/pendage |
|---|---|---|---|---|---|---|
| **production** (`tunnel_lisjak*_4h.cfg`) | 0,60 / 0,20 m | 117 132–117 200 | **75,5–76,6 mm** (iso 72,9) | 1,005–1,010 | **3,22 µs** (iso 3,05) | **≈ 4 h 30** |
| convergence (`tunnel_lisjak*.cfg`) | 0,35 / 0,12 m | 303 654 | 44,3 mm | 1,015 | ~1,9 µs | ~17 h |

Le réglage production a la **même densité de maille que le run isotrope de
référence** (4 h 18) : la brique 3 ne coûte que les 10 % de triangles des
lignes de litage. t/D = 1/18 contre 1/30 chez Lisjak (tunnel de 3 m) ; le
réglage convergence est là pour vérifier que h/v n'en dépend pas.

| clé (défaut) | brique | rôle |
|---|---|---|
| **`beddingDip`** | 1+2 | pendage du litage [deg depuis x] — clé maîtresse, refusée seule |
| `beddingEperp`, `beddingNuPerp`, `beddingGperp` | 1 | E', ν', G' — **les trois ensemble** |
| `beddingFtRatio`, `beddingCohRatio`, `beddingGfIRatio`, `beddingGfIIRatio` (1) | 2 | X_min/X_max |

**Validation (tests_f2/, 2026-09-01/02) :**

| contrôle | résultat |
|---|---|
| bit-identité, clés absentes | 28 fichiers, 0 octet |
| loi γ neutre (ratios = 1) | **bit-identique** au témoin |
| TI à la limite isotrope, élastique **sans Cundall** (`t13c`) | écart **exactement nul** sur 2080 lignes (avec Cundall : 7,5·10⁻⁴, signe de la vitesse sur des arrondis à 10⁻¹⁶) |
| modules apparents, mesure tout intérieur, adaptatif (`t14c`) | en travers **1,5910** (cible 1,5910) : **+0,00 %** ; le long **4,3364** (cible 4,3307) : **+0,13 %** |
| idem en intrinsèque (`t14b`) | −2,0 % / −4,8 % = complaisance de pénalité, **en série et indépendante de la direction** (0,0127 / 0,0117 GPa⁻¹) |
| loi γ, ratios Lisjak, litage 0 vs 90° (`t16`) | pic **4,35 → 8,43 MPa** (rapport 1,94 ; Lisjak T_P/T_S ≈ 1,9) |
| maillage lité, fumée (h 0,6, t 1,8) | 1619 arêtes exactement alignées (9,4 % contre 0,7 % isotrope), **continuité 1,007** |
| gardes (TI+thermal, TI+law, `beddingDip` seule, 3D) | refus propres |

> Le message `WARNING: moins de 12 % des joints sont quasi parallèles au
> litage` de la brique 2 est le rappel de la brique 3 : sans maillage
> préconditionné, la délamination n'a pas de chemin continu.

### 5.15 État de contact des joints et lits bimodaux (chantier f2, 2026-09-02)

Deux ajouts issus de la revue `tunnel_schisto/REVUE_traversee_litage.md`
(pourquoi les fissures n'ont pas traversé les plans de litage). Tous deux
opt-in, bit-identiques clés absentes (28 fichiers, 0 octet).

**`writeJointState` (false) — solution S8.** Ajoute à `fdem_joints_*.vtu`
quatre champs par joint, moyennés sur les deux points d'intégration au pas
courant : `sigN` (traction normale, > 0 en traction), `tauS`, `dn`
(ouverture) et `contactState` : 0 ouvert (D ≥ 1, dn > 0), 1 fermé-glissant
(|τ| à la limite de frottement), 2 fermé-bloqué, 3 cohésif (D < 1), 4 mort
(relayé au contact général), 5 lié (non inséré). C'est le **diagnostic de
Renshaw & Pollard (1995)** aux plans délaminés : un plan ouvert ou glissant
ne transmet aucune traction — l'arrêt d'une fissure y est *mécanique* ; un
plan fermé-bloqué qui arrête quand même, c'est le rapport d'énergies
(He & Hutchinson 1989) qui bloque. Dépouillement :
`tunnel_schisto/tools/joint_state_stats.py out --dip β` (état des
plans-frontières, σ_n/τ par état, **événements de traversée**). Coût : deux
`double` par joint et un `int`, écritures privées au joint (sûres sous OpenMP),
aucune force ajoutée.

**`weakPlaneFactor2` + `weakPlaneFrac2` (+ `weakPlaneSeed`) — solution S5.**
Lits faibles et lits forts : une fraction `weakPlaneFrac2` des plans de
`weakPlanes` reçoit le facteur `weakPlaneFactor2`, tirée **par index de
plan** (hachage déterministe : les trois pendages partagent la séquence).
Chandler et al. 2016 (Mancos) : 5 lits sur 7 faibles, 2 sur 7 forts. Le champ
`weakPlane` du VTU vaut 2 sur les lits forts. Les deux clés vont ensemble ;
absentes, un seul facteur, comportement inchangé.

**Validation (tests_f2/, 2026-09-02) :** `t17_jointstate` — champs présents,
15 joints rompus tous à l'état 0 (ouvert, essai de traction), 1179 cohésifs,
σ_n dans ±3,2 MPa ; `t18_wp_bimodal` — 7 plans, 13 joints sur lits forts
(28,6 % demandés), pic 5,84 contre 5,64 MPa monomodal.

**Outils ajoutés (`tunnel_schisto/tools/`) :** `edz_sectors.py` (profondeur
d'enveloppe par secteur corrigée de la paroi — `edz_metrics.py` mesure des
demi-axes horizontal/vertical et est **aveugle à une ellipse à 45°** : h/v =
1,01 sur un losange dont le rapport le long / en travers vaut 1,77),
`tip_velocity.py` (S7, vitesse de pointe vs c_R depuis `tBreak`). Suite de
decks et arbre de décision : `tunnel_schisto/SUITE_solutions.md`.

> ⚠️ **S6, à écrire sur toute figure.** Sur le litage en mode I,
> ℓ_cz = EΓ/σ² = 179 mm pour une arête de 200 mm : la zone de process n'est
> **pas résolue**. Le rapport G_Ic 0,057 agit alors comme un rapport
> résistance × ouverture critique, non comme un rapport d'énergies au sens
> des critères de déflection. Le maillage fin (t = 0,35, h = 0,12) la résout
> à 1,5 élément près.

### 5.16 Polydispersité des grains et taille par phase (chantier f2, 2026-09-02)

Demande : « faire varier les tailles de grain, pas le maillage, en gardant
les proportions globales de chaque phase ». Deux clés opt-in, bit-identiques
clés absentes (deck GBM `calib_quick/q3_gbm_P050.cfg`, 8 fichiers, 0 octet).

| clé | défaut | rôle |
|---|---|---|
| `grainSizeSpread` | 0 | écart-type de ln(taille) des grains, dans [0 ; 1,5] (0,3 modéré, 0,6 fort). Exige `grainSeeding = random`. |
| `phase.<nom>.grainSize` | absent | taille cible [m] de la phase (affinité) ; avec plusieurs phases seulement. |

**Méthode (Tessellation::build).**

1. *Graines* : les N espacements s_i = s·L_i (L_i log-normale de moyenne 1,
   sd ln = `grainSizeSpread`) sont tirés **d'abord** puis placés du plus
   **grand** au plus petit, chacun accepté à ≥ 0,35 (s_i + s_j) de tous les
   précédents (addition séquentielle triée). Tirer la taille à chaque essai
   rejette surtout les grandes graines. N est divisé par E[L²] = exp(σ²)
   pour conserver l'aire moyenne.
2. *Cellules* : **diagramme de Laguerre à aires prescrites** A_i ∝ s_i²,
   normalisées à W·H. Pour des graines fixées il existe des poids, uniques à
   une constante près, réalisant exactement toute famille d'aires positives
   de somme W·H (Aurenhammer, Hoffmann & Aronov 1998) ; on les obtient par
   **Newton amorti sur le dual semi-discret** (Kitagawa, Mérigot & Thibert
   2019) : départ w = 0 (Voronoï, aucune cellule vide), Laplacien
   dA_i/dw_j = −ℓ_ij/(2 d_ij) assemblé depuis les arêtes (voisin identifié par
   le test de puissance au milieu de l'arête), système singulier régularisé
   par 11ᵀ/N et résolu en gradient conjugué, amortissement gardant
   min A ≥ ½ min(A(w₀), A_cible). C'est la construction de Bourne, Kok, Roper
   & Spanjer 2020 pour des grains de volumes donnés (Neper : Quey &
   Renversade 2018 ; Falco et al. 2017). Convergence : 3–8 itérations,
   < 0,1 %. Lloyd déplace les graines vers les centroïdes des cellules de
   puissance et les poids sont **résolus à nouveau** à chaque cycle.
3. *Phases* : sans `phase.<nom>.grainSize`, chemin historique (ordre mélangé,
   glouton sur le déficit d'**aire** — les fractions sont respectées par
   construction quelle que soit la distribution). Avec, les grains sont pris
   du plus grand au plus petit et vont à la phase de score maximal
   exp(−½ (ln(d_g/d_phase)/0,35)²) × déficit relatif ; une phase saturée ne
   recrute plus tant qu'une autre a du déficit. C'est une **affinité**, pas
   une contrainte : la taille réalisée par phase est dans le journal
   (moyenne en nombre ± sd, et **pondérée par l'aire**, la seule qui compte
   pour « où est la biotite »).

**Journal :** `[tess] laguerre: aires prescrites … atteintes à x % en n
itérations`, `[FDEM] POLYDISPERSITE : N grains, écart-type de ln(d_eq)
RÉALISÉ = …`, une ligne par phase (fraction d'aire vs cible, nombre, d_eq).
WARNING si le réalisé < 60 % de la demande (domaine trop petit pour la queue
de la distribution). `ROCKIM_TESS_DEBUG=1` imprime chaque itération de Newton.

**Validation (`calib_quick/_poly_*`, 20 × 40 mm, grains 3 mm, 3 phases) :**

| cas | demande | réalisé sd ln d_eq | grains | fractions (62/31/7) | Newton |
|---|---|---|---|---|---|
| a | 0,5, Lloyd 2 | **0,496** | 88 | 61,4 / 32,5 / 6,1 | 3 it, 0,000 % |
| c | 0,5, Lloyd 0, tailles par phase | 0,491 | 88 | 62,0 / 30,8 / 7,1 | 5 it |
| f | 0,8, Lloyd 2 | **0,802** | 60 | 59,6 / 30,5 / 9,9 (1 grain de biotite) | 5 it, 0,02 % |

Le réalisé est celui de l'**échantillon** tiré (corr(ln s, ln d) = 1,000 dans
le réplica), pas la valeur asymptotique : à 60 grains l'écart type
d'échantillon fluctue. Deux impasses documentées, à ne pas rejouer : Voronoï
ordinaire avec espacement par graine (0,5 → 0,16 : la médiatrice moyenne les
tailles des voisines), poids **fixes** (κ s_i)² (0,5 → 0,20–0,24 : les
interstices du Poisson-disc se partagent au périmètre). Réplicas Python :
`calib_quick/_lag_experiment*.py`, `_lag_newton.py` ; greffes : `_patch_*.py`.

> ⚠️ **Coût.** hmin suit le plus petit grain : à σ = 0,8, hmin passe de 0,20
> à 0,055 mm (delaunay intra-grain à 0,18 d) — dt ÷ 3,6. Une phase à petite
> fraction reçoit peu de grains quand σ est grand (biotite 7 % : 1 grain de
> 10 mm en f) : donner alors `phase.biotite.grainSize`. Le réseau hexagonal
> (`grainSeeding = hex`) est refusé avec `grainSizeSpread`.

Références : Aurenhammer F., Hoffmann F., Aronov B. (1998) *Algorithmica* 20,
61–76 ; Kitagawa J., Mérigot Q., Thibert B. (2019) *J. Eur. Math. Soc.* 21,
2603–2651 ; Bourne D.P., Kok P.J.J., Roper S.M., Spanjer W.D.T. (2020) *Phil.
Mag.* 100, 2677–2707 ; Quey R., Renversade L. (2018) *CMAME* 330, 308–333 ;
Falco S., Jiang J., De Cola F., Petrinic N. (2017) *Comput. Mater. Sci.* 136,
20–28.

### 5.16 bis Maillage intra-grain non structuré : `grainMeshRandom` (2026-09-02, `rockim_f2n.exe`)

Remarque de Fernando (14:00) : « le maillage GBM est structuré dans les
grains ». Vérifié : le Delaunay intra-grain (`grainMesh = delaunay`) place ses
points intérieurs sur un **réseau triangulaire** de pas h — les triangles sont
quasi équilatéraux et alignés : orientations d'arêtes intra-grain **R6 =
0,548, pic/creux 18,8** (cas 3, 113 grains), pire que le frontal de Gmsh
banni (0,34). Trois directions de fissure imposées à l'intérieur des grains.

`grainMeshRandom = true` (défaut false, bit-identique) : points intérieurs
par **Poisson-disc** dans le polygone (distance ≥ 0,75 h, marge 0,55 h aux
arêtes, densité de saturation ≈ celle du réseau), puis le même Delaunay.
Mesure sur le même deck : **R6 = 0,007, pic/creux 1,37**, angle minimal
médian 44° (réseau : 49°), +15 % d'éléments (6688 vs 5814), dt −15 %
(2,65 vs 3,11 ns) → coût ≈ +30 %. Exige `grainMesh = delaunay`. Planche :
`calib_quick/fig_tess_grainmesh.png` ; métriques : `calib_quick/_intra_metrics.py` ;
suite `fast` 44/44 (`suite_f2n.txt`). Même binaire : `mesh = file` accepté avec
`geometry = disc` (le fichier est pris tel quel comme disque, méplats compris ;
`discR_` = W/2 des clés W/H) — brésilien `calib_quick/bts_v070b.cfg` sur le
disque Gmsh `make_disc_mesh.py` (Ø40, méplats 2 × 20°, R6 0,09), fissure
amorcée au centre, BTS = k_band × σ_t nominal (k_band = jauge élastique du
solveur, 0,89).
**Tout deck GBM de calibration doit le poser** ; les résultats GBM antérieurs
(cas 3 « −20 % », campagnes de juillet/août) ont été obtenus sur ce maillage
structuré et sont à relire avec cette réserve. La subdivision des arêtes de
grain reste uniforme (points partagés entre grains voisins).

### 5.17 Clés de calibration triaxiale (chantier f2, 2026-09-02, `rockim_f2m.exe`)

Quatre ajouts opt-in issus de la critique adverse du plan de calibration
(`calib_quick/enquete/A5_critique_adverse_du_plan.md`, A3 « à ajouter au
code ») ; bit-identiques clés absentes sur deux decks (GBM tension et mors
`mesh = file`, 8 fichiers chacun, 0 octet).

| clé | défaut | rôle |
|---|---|---|
| `historyStrains` | false | scénario tension (mors **et** plateaux) : trois colonnes en fin de `history.csv`, `epsAx` = (ū_y haut − ū_y bas)/H, `epsLat` = (ū_x droite − ū_x gauche)/W, `epsVol` = epsAx + epsLat (déformation plane, ε_zz = 0), **traction positive**. Moyennes sur les copies de nœuds des quatre faces de la boîte (pondération par les éléments incidents). Sortie seule. C'est ce qui permet les seuils σ_ci/σ_cd par la méthode SBM (inversion de ε_v) avec le même opérateur que sur l'essai. Exige une boîte rectangulaire. |
| `gripsStopAfterPeak` (+ `gripsStopDelay` [s], 0) | false | miroir de `ucsStopAfterPeak` pour le montage à **mors** : le run s'arrête `gripsStopDelay` après le verrouillage du pic (`peakLocked`). Divise par 2–3 le coût des runs qui cassent beaucoup (146 → 530 s de post-pic profond sans valeur de mesure). |
| `stopPeakDrop` | absent (= 0,7 de chute) | fraction de chute sous le pic qui verrouille : `sigma < (1 − stopPeakDrop)·sigmaPeak` au lieu de l'historique `sigma < 0,3·sigmaPeak` (chute de 70 %, **inatteignable sous confinement** : sigma est la contrainte totale, bornée par σ₃). Vaut pour mors et plateaux ; dans ]0 ; 1[. Exemple triaxial : `stopPeakDrop = 0.3` + `gripsStopDelay = 2e-4` (Δε ≈ +0,1 % après la chute à 0,25 m/s sur 40 mm, la fenêtre de l'observable « chute »). |
| `weibullScope = lcz` | `strength` | troisième portée du Weibull de joint : G_f et G_II suivent **stat²**, donc ℓ_cz = E G_f/ft² est constante joint à joint — disperser la résistance sans changer la ductilité (règle du balayage `calib_quick`). La valeur est désormais validée aussi avec un bulk élastique (avant, `MatLaw::make` ne la vérifiait que sous `law` : une faute de frappe passait en silence). |
| — (avertissement) | — | `jointPenaltyFactor` posé avec `insertion = adaptive` → `[FDEM] WARNING: … INERTE …` (la pénalité effective des joints insérés est `insertionPenaltyFactor`, 4 E/h par défaut ; le 3D avertissait déjà, pas le 2D). |

**Validation (`calib_quick/_fk_*`, 2026-09-02) :** en-tête
`…,nInserted,nDamaging,epsAx,epsLat,epsVol` ; après consolidation à 50 MPa
(mors bloqués) epsAx = −1,6e-5, epsLat = −4,8e-4 (élastique :
[σ_xx(1−ν²) − νσ_yy(1+ν)]/E = −5,1e-4 ✓) ; `stopPeakDrop = 2` et
`weibullScope = xyz` refusés ; suite `fast` **44/44** (`suite_f2m.txt`).

## 6. Sorties

Tous les fichiers vont dans le dossier de sortie. Fréquences : VTU toutes les
`T/frames`, une ligne d'historique tous les ~1/2000 du run.

### 6.1 `history.csv` — colonnes par scénario (mode fdem)

| scénario | colonnes |
|---|---|
| percussion / shear | `t, toolFx, toolFy, toolX, toolY, toolVx, toolVy, work, toolKE, nBroken, nFrag, detachedVol, specificEnergy` + `eEl, eJnt, eGc, eFric, eCund, eLys` (V2/B4 : travaux cumulés par famille, signés — négatif = prélevé au solide) |
| tension (grips) | `t, gripFy, sigma, sigmaPeak, nBroken` |
| tension (platens) | + `epsPlaten, epsSpec, epsGauge, nBrokTen, nBrokShear, nFrag, confAchieved, peakLocked` |
| brazilian | `t, P, Pbot, drive, sigmaT, sigmaTpeak, nBroken, nFrag, sxxC, syyC, peakLocked` |
| shpb | `t, vDrive, epsM1, epsM2, sxxC, syyC, nBroken, nFrag, nInserted` |

Avec `hydro = on`, cinq colonnes s'AJOUTENT à celles du scénario :
`hydroP` (pression de puits [Pa] — c'est la courbe de leur fig. 11), `hydroVol`
(volume de cavité [m³/m]), `hydroMass` (la variable d'état [kg/m]), `hydroNWet`
(faces mouillées — elle décolle de sa valeur initiale dès qu'un joint rompt, et
le compte est exact : chaque joint livre ses deux lèvres) et `eHydro` (travail
du fluide sur le solide, poste séparé du bilan B4). ⚠️ **`nBroken` ne compte que
les joints ENTIÈREMENT rompus (D ≥ 1)** : l'amorçage réel est plus précoce, et
seul le champ `damage` des VTU joints le montre. Ni le nombre de joints insérés
ni l'endommagement maximal ne sont écrits dans `history.csv` — lacune connue.

fdem3d : `t, gripFz, sigma, sigmaPeak, nBroken` (tension) ; percussion/shear comme en
2D avec les trois composantes (+ `grpZ, grpVz` si `trackGroup`, + les six colonnes
énergie V2/B4). fem/fem3d/dem/dem3d : variantes proches (force outil,
travail, casse). ⚠️ `sigmaPeak` est un max glissant qui attrape la sonnerie
post-rupture : recalculer les pics depuis les courbes échantillonnées, ou lire le pic
verrouillé (`peakLocked`).

**Bilan d'énergie par sous-système (V2/B4, fdem + fdem3d, 2026-08-14).** Le
résumé de fin de run imprime un bloc `energy budget` : théorème
travail-énergie sur les nœuds, `KE(t) − KE(0) = Σ travaux par famille +
résidu`. Postes : éléments (−elWork, avec l'élastique stocké à la volée lu
sur le Cauchy stocké — invariants isotropes, exact en élastique, approché
sous law/caps), joints cohésifs (fissuration + stocké = −(jointWork −
dampWork)), dashpot, contact général (dont part frottement `gcFricWork`),
Cundall, frontières (amortisseurs + ressorts stockés), outil→solide,
platines, et **`integration`** : la correction leapfrog EXACTE
`f²dt²/2m` par nœud et par pas — les compteurs par famille lisent v⁻, le
théorème discret veut (v⁻+v⁺)/2 ; sur la percussion 2D grille ce poste vaut
+2443 J/m (forces de contact violentes du cas divergent connu) et sa prise
en compte fait passer le résidu de 91 % à **0,017 %**. Mesures de recette :
percussion 2D 0,017 %, percussion 3D 0,005 %, zeroload deux-corps −5e-24 J.
Le verdict (OK ≤ 1 % du flux BRUT échangé ; « zero machine » à charge
nulle) est verrouillé par l'extracteur `budget` de la suite
(zeroload_bench1_3d ≤ 1e-12, bench1_insert_impact ≤ 1 % de KE₀).
Instrumentation PURE : aucune trajectoire ne change (suite fast 12/12
bit-identique) ; les compteurs multi-threads se réduisent en ordre de
thread (même statut que dampWork_). Périmètre : percussion/impact complet ;
en quasi-statique les platines sont comptées à v imposée (approx O(dt)).

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
`tools/crater_metrics.py <run>` (V2/B3 : métriques de cratère depuis les VTU joints —
R_crater p95 dans la peau de surface, R_max, profondeur, aire cassée, volume
endommagé/détaché, fissures radiales cassées ET bras endommagés par secteur
angulaire ; multi-corps V1 géré : surface et fragments rapportés au CORPS IMPACTÉ,
l'insert n'est pas un débris ; V2/B5 : `--brush beta` (0.8) — volume brossable =
β × fragments détachés dont le centroïde est dans le bol du cratère, la masse
collectée du banc étant ρ × volume ; `--plot` vue de dessus, `--csv` export),
`tools/make_unstructured_mesh.py` (maillages simplexes non structurés uniformes via
Gmsh — `box3d W D H h out.msh [seed]` / `box2d W H h out.msh [seed]` — pour
`mesh = file` ; `pip install gmsh`).

**Post-traitement du couplage hydro** (`bench_abuaisha/tools/`, 2026-08-20) :
`hydro_sign_check.py <run_conf> <run_hydro>` (LE contrôle de signe, cf. §5.10),
`fig_controle_run.py <run>` (planche de diagnostic utilisable sur un run EN
COURS : pompe, volume, champ, marge à la rupture), `fig_b2.py` et `gif_b2.py`
(planche livrable et animation trois panneaux), `fig_postpic.py` (phase
post-pic, avec la limite de zone raffinée tracée — une aile qui touche ce cercle
est bornée par le MAILLAGE, plus par la physique), `fig_vitesse.py`
(trajectoires dans le champ de vitesse, format de leurs fig. 12-13),
`fig_ouverture.py` (**l'ouverture des fissures**).

⚠️ `fig_ouverture.py` reconstruit une donnée que **rockim n'écrit pas** : le
writer ne pousse qu'UNE lèvre par joint (`lines.push_back({J.a1, J.a2})`). La
seconde se retrouve sans toucher au solveur — les nœuds sont dédoublés par
élément selon `n = 3e + k`, les copies d'un même sommet sont confondues à
l'instant initial, et chaque arête géométrique est alors portée par exactement
deux éléments (284 124 arêtes internes pour 284 124 joints sur le maillage B2).
L'ouverture vaut ensuite (b − a)·n, la formule même dont le solveur se sert pour
son volume de cavité.

**Règle maison maillage (2026-08-11)** : le maillage de BASE de toute étude est
DÉSORDONNÉ (`mesh = file` non structuré, ou `mesh = voronoi` si le sujet est le
GBM). Les grilles régulières sont réservées aux vérifications qui en ont besoin
par construction ; le solveur imprime un WARNING si un scénario de fissuration
part sur `mesh = grid` (mesuré : la grille de Kuhn 3D intrinsèque part en cascade
énergétique en phase débris là où le même cas sur maillage non structuré est sain).
