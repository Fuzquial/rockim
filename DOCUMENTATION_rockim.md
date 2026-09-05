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
| `kpFactor` (1.0) | pénalité outil = facteur·E·t ; **accepté comme alias de `contactPenaltyFactor` en fem3d** depuis la revue w12 (auparavant ignoré en silence ; les deux clés avec des valeurs différentes sont refusées) | fem, fem3d (alias) |
| `contactPenaltyFactor` (1.0) | pénalité outil nœud-sphère/plan = facteur·E·h_min (`rockim_f2w12.exe`, §5.18 ; à 1 aucune opération, bit-identique) ; dt reste borné par 2√(m_min/kp), rapport ω_p/(2/dt) imprimé ainsi que le rapport **combiné** √((dt/CFL)² + (ω_p dt/2)²) (revue : le nœud de contact porte aussi la raideur d'élément), avertissement si le combiné > 0,5 ; rien d'imprimé en `toolContact = signorini` (kp inutilisé) | fem3d |

Contact général (débris, fdem/fdem3d) : `gcPenaltyFactor` (0.01 ×E·t — mou exprès),
`gcXi` (0.8), `gcRestitution` (0.2, quasi-plastique), `gcBirthTau` (1e-6 s, relaxation
de la pénétration de naissance) ; SHPB : `gcCell` (2·hDisc), `gcBoxMesh` (true),
`gcXwindow` (0.10 m autour du disque — le contact est le SEUL chemin de charge).

### 5.6 bis Loi de contact de l'OUTIL et bancs de la pompe (T0/T1, 2026-09-02)

*Section ajoutée le 2026-09-02. Les clés `toolContact`, `toolSignoriniRelax`,
`toolImpulseCap` et `cutterThick` existaient — écrites, argumentées en
en-tête de `FdemSolver.hpp` — et n'étaient documentées **nulle part**, ni dans
ce fichier, ni dans `CHANTIER_f2.md`, ni dans le README, ni dans la suite de
non-régression. Dette du principe VII, soldée ici avec les deux bancs qui
leur donnent un verdict.*

**Le problème que ces clés traitent.** Le rapport
[`coupe_pdc/RESULTATS_2026-08-18.md`](../rockim/coupe_pdc/RESULTATS_2026-08-18.md)
mesure, sur la coupe PDC 2D, une **pompe d'énergie dans le contact outil** :
77 286 J/m injectés dans le solide pour 189 J/m de travail de corps rigide
(**facteur 408**), des nœuds à **2 544 m/s** contre une borne physique de
2·v_outil = 20 m/s (**facteur 127**) — le tout pendant que le résidu du bilan
B4 affichait `[OK]` à 1,9e-10 %. **Le résidu ne peut pas voir une pompe logée
dans un canal COMPTÉ.**

La cause structurelle : en percussion l'outil est **libre** (il décélère,
réservoir fini) ; en coupe il est à **vitesse imposée** (réservoir infini).
Toute fuite de pénalité y devient non bornée. Et l'écrêtage de la voie
pénalité est **géométrique** (0,6 h) : il borne la *pénétration*, pas
l'*impulsion*, qui est la grandeur qui lance les nœuds.

| clé (défaut) | rôle | portée |
|---|---|---|
| **`toolContact`** (penalty) | penalty \| **signorini** = contact en **condition de VITESSE** (CD-Lagrange). Relation impulsion / saut de vitesse (lemme de viabilité de Moreau) : `si g > 0 alors r = 0 ; sinon 0 ≤ v ⊥ r ≥ 0`. La masse de rockim étant diagonale et l'outil un obstacle rigide, l'opérateur de Delassus `H = L M⁻¹ Lᵀ` est **diagonal et sphérique par nœud** : l'impulsion se calcule en **forme fermée**, sans système à résoudre, le solveur reste matrix-free. **`kp_` sort du budget de `computeStableDt()`** — le pas cesse de dépendre d'une pénalité arbitraire (gain mesuré sur ce deck : faible, la pénalité de joint domine ; le gain est sur la **physique**, pas sur le mur). Sources : Fekak, Brun & Gravouil (2017) ; Dureisseix, Greco & Brun, JTCAM (2024), Alg. 1 ; Ghesquière-Diérickx, Anciaux, Acary & Molinari, arXiv:2606.01355. Algèbre extraite dans `include/rockim/ToolSignorini.hpp` (la **géométrie** reste dupliquée dans le solveur, à dessein) | fdem, fdem3d |
| `toolSignoriniRelax` (0) | rattrapage façon Baumgarte de la pénétration **résiduelle**, dans [0, 1]. 0 = condition de vitesse **pure** : on annule l'approche, on ne résorbe pas la pénétration déjà acquise. Ce qu'il faut accepter en Signorini : une interpénétration résiduelle subsiste (les auteurs la quantifient à η = 0,43 % sur leur cas de référence) | idem |
| `toolImpulseCap` (0 = off) | écrêtage `\|Fc\| ≤ κ·2·\|v_outil\|·m/dt`. **Correct en direction, faux en formulation** (mesuré le 2026-08-18) : plafonne l'incrément **par pas**, alors que la borne physique porte sur l'impulsion de **toute** la collision — un nœud en contact soutenu prend 2v à chaque pas, seize pas donnent 311 m/s. Conservé comme instrument ; `toolContact = signorini` est le remède | fdem |
| `cutterThick` (0 = coin infini) | épaisseur du cutter derrière la face de coupe (**face de dégagement**). Sans elle la roche située derrière l'arête tombe dans un demi-espace infini et se fait labourer. Correctif physiquement juste, mais **ce n'était pas la cause** de la pompe (v2 et v3 rigoureusement identiques jusqu'à la trame 10) | fdem shear |

> ⚠️ **`chamferLen` / `chamferDeg` sont lus, validés, stockés — et
> n'interviennent NULLE PART dans le calcul du contact.** Un avertissement est
> imprimé à l'exécution. Anomalie consignée dans la spec `003-cutter-pdc-3d`.

**Les deux indicateurs, imprimés à CHAQUE run avec outil depuis cette date.**
Le rapport du 2026-08-18 concluait : « aucun des deux ne coûte quoi que ce
soit à calculer ; ils devraient être imprimés à chaque run avec outil ».
C'était resté lettre morte alors que **les deux nombres du facteur 408
étaient déjà imprimés tous les deux** — personne ne faisait la division.

```
[FDEM] injection outil   : <toolWork_> J/m vers le solide / <work_> J/m corps rigide = ratio <R>  [POMPE]
[FDEM] v nodale max      : <v> m/s = <X> x 2 v_outil (<2v> m/s)  [HORS BORNE]
```

Drapeau `[POMPE]` au-delà de ratio 2. Drapeau `[HORS BORNE]` au-delà de
**2 × la borne**, et non de la borne : 2·v_outil vaut pour la contribution du
**contact seul** (choc contre une masse infinie), or dans un continuum l'onde
émise se réfléchit sur une surface **libre** en doublant la vitesse
particulaire — un nœud de peau peut donc légitimement approcher le double.
La marge sépare sans ambiguïté (mesure T1 : pénalité 8,25, Signorini 1,08).
`vNodeMax_` est échantillonné tous les 1024 pas ; **limite assumée** : un pic
isolé plus bref passe entre deux mesures — la pompe, elle, est persistante.
Instrumentation **pure** : aucune force, aucune trajectoire ne change.

**Les deux bancs (~65 s au lieu des ~1 h 20 d'un run de coupe).**

`rockim selftest-toolcontact` (**T0, 0,1 s, sans maillage**) vérifie
l'**algèbre** du noyau réel en forme fermée — sept familles : condition de
Signorini (séparation ⇒ impulsion nulle, testée au seuil exact), théorème du
nœud au repos (repart **exactement** à v_outil : contact inélastique de
Moreau à relax = 0, jamais 2 v_outil), **invariance d'échelle** sur six
décades, absence d'adhésion + charge nulle ⇒ impulsion nulle, cap de Coulomb
sur l'impulsion (glissement *et* collage), **dissipativité** dans le repère de
l'outil sur un balayage (approche × glissement × µ), et effet réel de
`toolSignoriniRelax`. **Écart machine : 2,2e-16.** Le banc imprime en outre le
contraste pénalité sur les mêmes nombres : `F` écrêtée 7,24 MN/m → **373 m/s
en UN pas**, soit 18,7 × la borne, contre 10 m/s pour Signorini.

`configs/verify_fdem_toolcontact.cfg` (**T1**) est le banc de **raclage** :
bloc 4 × 2 mm, 334 éléments, cutter PDC à **10 m/s**. On teste à 10 m/s et non
à 1 parce que la borne est **sans échelle** (T0 cas C3) : c'est dix fois moins
cher *et* plus sévère. Le banc a des **dents** — le cas pénalité est là pour
**échouer** :

| | pénalité | signorini | |
|---|---|---|---|
| injection outil / corps rigide | **4,43** `[POMPE]` | **0,905** | ÷4,9 |
| v nodale max / 2 v_outil | **8,25** | **1,08** | ÷7,6 |
| joints rompus | 2 | 0 | |
| pic de force outil | 2,03 MN/m | 0,206 MN/m | ÷9,9 |
| résidu B4 | 1,7e-12 % `[OK]` | 2,6e-13 % `[OK]` | **aveugle dans les deux cas** |
| durée | 38 s | 22 s | |

> ⚠️ **Le cas pénalité est CHAOTIQUE hors `OMP_NUM_THREADS = 1`** : le même
> deck en multi-thread donne 7,31 et 8 joints rompus au lieu de 4,43 et 2. La
> suite fixe OMP = 1, où trois runs de recette sont **bit-identiques**. Le cas
> Signorini, lui, est reproductible.

> ⚠️ **CE QUE T1 EST, ET CE QU'IL N'EST PAS** (requalifié le 2026-09-02 soir
> après relecture du code et de `out_t1_*/history.csv` — la première rédaction
> lisait le pic de force comme un verdict, à tort).
>
> **T1 est un banc de POMPE. Ce n'est pas un banc de FORCE**, et il ne peut pas
> l'être : le bloc y est **LIBRE**. L'encastrement du fond n'est posé qu'en
> scénario `percussion` (`FdemSolver.cpp:3685-3686`, `flag_ = FIXED` sous
> `scen_ == Scenario::PERCUSSION`) ; en `shear` aucune condition ne retient le
> bloc, et le deck ne pose pas `absorbing`. Le bloc de 0,021 kg/m est donc
> **chassé par l'outil** : contact de 37,4 à 42,2 µs, soit 33 relevés sur 2 015
> (1,8 % du temps en Signorini, 6,4-8,6 % en pénalité), vitesse moyenne du bloc
> 10,2 m/s en fin de run, boîte translatée de 2,5 mm et tombée de 1,4 mm.
>
> **Le « pic divisé par 9,9 » n'est donc pas un signe de mollesse** : il compare
> deux pics dont l'un est produit par un run qui pompe. Sur la course engagée,
> les forces **moyennes** sont dans un rapport **×1,40** (1 698 contre
> 1 217 N/m), pas ×9,9.
>
> **Et zéro joint rompu est un fait de MATÉRIAU, pas de contact** : avec
> ft = 10,62 MPa et Gf = 152,9 J/m², ℓ_cz = E·Gf/ft² = **65 mm** pour un bloc de
> 4 mm, et l'ouverture de rupture dnF ≈ 2Gf/ft ≈ 30 µm est inatteignable en
> 4,8 µs de contact. Le banc ne *pouvait pas* casser. (227 joints insérés,
> D max = 0,062.)
>
> **Le verdict sur la force demande donc un autre banc** : bloc **tenu**
> (`shearSupport = fixedBottom`, §5.6 ter), `dampingLocal = 0` — le Cundall
> s'applique à l'impulsion de contact (`:7396-7404`) et absorbe 51 % du travail
> outil sur T1 — et un ℓ_cz à l'échelle du bloc. C'est T1h (élastique, 1 m/s,
> où la pénalité est une référence valide) puis T1b (qui casse). Voir le plan
> par étapes dans `panel_contact_2026-09-02/`.
>
> Enfin, ces bancs ne couvrent que le **canal outil** ; le canal **joint**
> (branche de compression, `jointContactPenalty = adaptive`) a sa propre
> famille, sans outil : `jointdeath_tension_2d` et `zeroload_jointdeath_2d`.

Repères de suite : `selftest_toolcontact`, `t1_toolcontact_penalty`,
`t1_toolcontact_signorini` (tier `fast`).

### 5.6 ter Instruments par canal, banc T0b et clés de montage de la coupe (étapes 1-3, 2026-09-02)

*Ajouté le 2026-09-02 au soir. Ces quatre briques répondent au constat du panel
de conception : le ratio d'injection outil est **≤ 1 par théorème** sous
`toolContact = signorini` (`toolWork_` lit v⁻ ; le ½mv² du premier toucher tombe
dans le poste « intégration »), donc il ne peut RIEN dire des canaux joint et
fragments — et le résidu B4 est aveugle à une pompe logée dans un canal compté.
Sans instruments par canal ni bloc tenu, aucune des étapes suivantes n'a de
critère de mort.*

#### Le banc T0b — « Signorini e = 0 est-il MOU ? » : non, et c'est démontré

Cas C8/C9 de `rockim selftest-toolcontact` (**+0,0 s**, aucun maillage). Une
chaîne de N masses lumpées + ressorts (granite Heilman) est frappée par le mur
de Signorini à 10 m/s, en appelant le **vrai** noyau `toolsig::impulse`.

Le soupçon à lever : le contact annule la vitesse d'approche (e = 0 **au
nœud**), donc un nœud isolé repart à v_outil et non à 2 v_outil (cas C2). Le
théorème de Saint-Venant dit que dans un **solide** la restitution n'est pas
portée par le nœud mais par l'**onde** : l'onde de compression se réfléchit en
traction sur l'extrémité libre, revient, et la barre se sépare à t = 2L/c en
repartant à **2v**.

| mesure | N = 100 | N = 400 | théorie |
|---|---|---|---|
| v_sortie | **19,816 m/s** | **19,936 m/s** | 2v = 20 |
| durée de contact | 46,558 µs | 46,569 µs | 2L/c = **46,511 µs** |
| perte d'énergie | **0,9997 %** | **0,2475 %** | 1/N = 1 % et 0,25 % |

La perte est **exactement** celle du premier toucher (½m₁v² = KE/N) et elle
**divise par 4,04** de N = 100 à N = 400 — donc en 1/N, pas un seuil. Le noyau
ne dissipe rien de ce que l'onde doit rendre. **Signorini n'est pas mou.**

> Le bilan d'énergie du banc **doit** inclure l'énergie élastique stockée : la
> barre repart en vibrant. Sans ce terme la « perte » vaut 2,34 % au lieu de
> 1 % à N = 100, avec un rapport 3,16 au lieu de 4 entre les deux maillages —
> la signature d'un terme oublié, pas d'une dissipation. Erreur commise puis
> corrigée à la première exécution.

**C9** rejoue la même barre à `dampingLocal = 0,05` : v_sortie tombe à
**16,89 m/s** (−15,5 %). Le Cundall d'`integrate()` s'applique à la force
totale, **impulsion de contact comprise** (`:7396-7404`) — c'est un amortisseur
de contact déguisé. D'où la règle : **`dampingLocal = 0` sur tout banc de
force**.

#### Les compteurs PAR CANAL (print-only, bit-neutres)

| sortie au résumé | ce qu'elle dit |
|---|---|
| `injection (trapeze)` | le travail outil compté en r·(v⁻+v⁺)/2, qui ne sous-compte pas. Le ratio de la ligne au-dessus est ≤ 1 **par théorème** en Signorini : il ne peut pas tirer |
| `canaux (travail +)` | somme des incréments **POSITIFS** par pas des joints et du contact général, en J/m et en % du travail outil. Un canal sain oscille autour de zéro, un canal qui pompe accumule |
| `noeuds profonds` | rejets `d < −capk` : un nœud plus profond que l'écrêtage QUITTE l'ensemble actif en silence — c'est une traversée que l'audit de pénétration ne voit pas (« 0,141 mm = exactement l'écrêtage ») |
| `contact outil` | fraction d'évaluations **COLLÉES** (cap de Coulomb non atteint). À µ = 0,8 et rake 20° le noyau colle : **44,5 % sur T1**, donc le calibreur F_v/F_h = tan(θ + atan µ), qui suppose le glissement, ne s'applique pas tel quel |

`vNodeMaxEvery` (défaut **1024**, inchangé) règle la cadence d'échantillonnage de
la vitesse nodale maximale ; poser **1** sur un banc, où un pic isolé ne doit pas
passer entre deux mesures.

> ⚠️ **Limite du « travail positif ».** Un ressort qui vibre accumule lui aussi
> du travail positif à chaque demi-cycle : cet indicateur ne distingue pas seul
> l'oscillation de la pompe. Il se lit **avec** `v nodale max` (borne physique
> dure) et avec le net (`eJnt`, `eGc` de l'history). C'est le couple qui tranche.
>
> ⚠️ **Limite du trapèze.** v⁺ y est *prédit* avant le Cundall et les Lysmer
> d'`integrate()` : sous amortissement le compteur SURESTIME l'injection (1,136
> sur T1 à `dampingLocal = 0,05`). Le lire à amortissement nul.

#### Les deux clés de montage

| clé (défaut) | rôle |
|---|---|
| **`shearSupport`** (none) | none \| **fixedBottom** = fond encastré en scénario de **coupe**. L'encastrement n'existait qu'en `percussion` (`:3685-3686`) : **aucune combinaison de clés existantes ne tenait un bloc de coupe** — `absorbing` ≠ none pose ressorts + Lysmer sur les DEUX flancs (donc sur la face d'ENTRÉE que l'outil engage, ~1,2e9 N/m par nœud contre 0,107 MN/m de pic outil), `absorbSpringR` est une clé unique lue pour les flancs ET le fond, `lateralRollers` impose u_x = 0 sur cette même face. Avertissement à l'exécution si `absorbing` ≠ none |
| **`cutterFloor`** (false) | true = un nœud sous la ligne d'arête (`rel.y() < 0`) n'est jamais en contact. **Le code faisait le contraire de son commentaire** : le test `dt2 < 0` exclut selon la direction de la face INCLINÉE, pas selon l'horizontale — pour `rel = (−0,20 ; −0,05)` mm à 20° on obtient dn = −0,205 et dt2 = **+0,021**, donc un nœud 0,05 mm SOUS l'arête est déclaré en contact. Bande labourée = `cutterThick·sin(rake)` : 0,171 mm sur T1, **0,855 mm sur v3, soit 84 % de la passe** |
| **`toolSignoriniGroup`** (false) | true = **une seule** impulsion par GROUPE de copies liées (M = Σm_i, F = Σf_i), redistribuée au prorata des masses, au lieu d'une décision par copie. Sous `insertion = adaptive`, `integrate()` intègre les groupes comme un seul nœud : l'opérateur de Delassus du solveur est diagonal PAR GROUPE (1/M), pas par copie. Séquentiel à dessein (instrument de mesure). Exige `toolContact = signorini` |

#### Ce que ces clés donnent sur T1 (mesuré, `OMP_NUM_THREADS = 1`)

| variante | injection (trapèze) | v_max / 2v | joints rompus | canal joints |
|---|---|---|---|---|
| Signorini seul (référence) | 1,1356 | 1,079 | 0 | 3,2 % |
| `+ toolSignoriniGroup` | 1,1354 | 1,077 | 0 | 3,2 % |
| `+ cutterFloor` | 1,1356 | 1,079 | 0 | 3,2 % |
| **`+ shearSupport = fixedBottom`** | **0,799** | **22,21 [HORS BORNE]** | **33** | **499 %** |

**Trois lectures, dans l'ordre d'importance :**

1. **Le bloc tenu change tout.** Contact soutenu (228 224 évaluations contre
   6 603, ×35), 33 joints rompus contre 0 — et **v_max à 22 × la borne** pendant
   que le canal **outil reste propre** (injection 0,80 < 1). C'est exactement ce
   que les compteurs par canal ont été construits pour voir, et que ni le ratio
   d'injection ni le résidu B4 ne montraient. **Signal fort que le canal JOINT
   pompe une fois l'outil propre** — à confirmer par T1b (étape 5), qui abaisse
   Gf pour que ℓ_cz tienne dans le bloc et pose `dampingLocal = 0`. Ce run-ci
   garde ℓ_cz = 65 mm ≫ 4 mm et un Cundall à 0,05 : ce n'est pas encore le banc.
2. **L'approximation par copie est négligeable** : 0,01 % sur l'injection,
   0,2 % sur v_max. Très en dessous du seuil de 5 % qui aurait rendu
   `toolSignoriniGroup` obligatoire. Le panel avait raison de dire que la
   linéarité sauve le cas d'un plan rigide touchant toutes les copies —
   c'est maintenant **mesuré**, plus supposé. La clé reste disponible pour les
   géométries où la compensation ne vaut plus.
3. **`cutterFloor` n'a AUCUN effet sur T1** — attendu : le bloc est libre et le
   contact dure 4,8 µs, le labourage sous l'arête ne se produit qu'en coupe
   engagée. Son verdict appartient à la porte v3 (étape 6), où la bande vaut
   0,855 mm.

Repères de suite : `selftest_toolcontact` (T0 + T0b), `t1_toolcontact_penalty`,
`t1_toolcontact_signorini` — références **inchangées** au bit près, l'ensemble de
ces ajouts étant de l'instrumentation et des clés opt-in.

### 5.6 quater Cutter PDC 3D — `toolShape = pdc` (2026-09-03)

**Ce que c'est.** La forme de cutter PDC pour le solveur 3D, afin de reproduire
Heilman et al., ARMA 24-0238 (Los Alamos, HOSS) : disque de 13 mm de diamètre,
2,5 mm d'épaisseur, arête chanfreinée, rigide, à vitesse imposée, sur un granite
Utah FORGE 40 × 30 × 20 mm. **La cinématique de coupe existait déjà** en 3D
(`scenario = shear` : `placeTool` pose un outil prescrit avec `cutDepth`,
`cutSpeed`, `toolX` ; le fond z = 0 est encastré ; `work_` et `toolFx..Fz` sont
dans `history.csv`). Il ne manquait que la forme — `toolShape` n'acceptait que
`sphere | flat | none`.

**Ce que ce n'est pas : l'extrusion du coin 2D.** En 2D le cutter est un coin
infini en profondeur testé par trois booléens sans rapport (`dn`, `dt2`,
`floorFlat`) et sa normale est **gelée** à `rakeNormal()` quelle que soit la
facette touchée (`FdemSolver.cpp:6259`). Ici le cutter est un **cylindre fini
chanfreiné, convexe**, et l'appartenance sort d'une **distance signée exacte**
dans le demi-plan méridien (ρ, s) : quatre demi-plans à gradient unitaire (face
de coupe, dos, tranche, chanfrein), `pen = −max(dᵢ)` à l'intérieur, normale de
la facette active. La surface de contact **s'élargit avec l'enfoncement** — c'est
l'argument de l'article pour le rôle de la garde arrière. Noyau et conventions :
`include/rockim/ToolPdc3d.hpp` (sans état de solveur, testable seul).

**Repère** : convention de `Tool.hpp:67-74` relevée y → z, pour qu'un même signe
de `backRakeDeg` soit la même physique en 2D et 3D. `tool_.x` est l'**arête de
coupe** (comme en 2D), n = (cos b, 0, sin b), u = (−sin b, 0, cos b), w = e₂.
**Signe de la garde** : à b < 0 le corps monte en arrière, rien sous l'arête
(orientation de `cut2d_v3ter_rake.cfg` et d'un vrai PDC) ; à b > 0 le dos
descend de t·sin b sous la ligne d'arête (0,855 mm à 20°) et laboure le
plancher — le solveur avertit, `cutterFloor = true` le rend à la roche.

| clé | défaut | rôle |
|---|---|---|
| `toolShape = pdc` | — | active la forme ; exige `scenario = shear` |
| `cutterDia` [m] | **requis** | diamètre du disque (13 mm) — `cutterLen` (clé 2D) **lève une erreur** |
| `cutterThick` [m] | **requis** | épaisseur (2,5 mm) ; sans elle le cutter serait un demi-espace qui piège le copeau |
| `backRakeDeg` | −20 | garde arrière, (−60, 60) |
| `chamferLen`, `chamferDeg` | 0, 45 | chanfrein sur la face de coupe ; **opérant** en 3D (contrairement au 2D) ; `chamferLen·tan(chamferDeg) < cutterThick` |
| `cutterFloor` | false | jamais de contact sous le niveau de l'arête ; `= flat` est **refusé** (getb le lirait FALSE) |
| `toolX`, `toolY`, `cutDepth`, `cutSpeed` | −2 mm, D/2, 4 mm, 10 m/s | position de l'arête et cinématique |

**Contact.** Voie pénalité : écrêtage en profondeur `0,6·h_el` (miroir 2D).
Voie Signorini : même noyau géométrique, impulsion inchangée. **Signorini est
la voie à utiliser** : à vitesse imposée la pénalité pompe (T1 2D : ×4,43), et
son écrêtage de 0,6 h atteint l'axe médian du disque dans l'anneau de l'arête,
où la normale saute de 90° (le chanfrein divise le saut par deux, il ne l'abolit
pas). Le solveur avertit fort en pénalité.

**Banc** : `rockim selftest-pdc3d`, < 1 s, sans maillage, entrée `selftest_pdc3d`
du tier fast. G1 pénétration = distance exacte au polygone méridien (20 000
tirages, 7·10⁻¹⁷ R) et sa variante qui **doit échouer** (noyau vif sur un solide
chanfreiné : 3·10⁻² R) ; G2 p + pen·N tombe sur la peau (3·10⁻¹⁶ R) ; G3 largeur
de contact à la surface libre 2√(R² − (d/cos b − R)²) = 7,180 mm pour l'article,
6,979 mm à b = 0 = le 2√(d(D−d)) du deck 2D, et la formule 2D qui **doit
échouer** à 20° (0,2 mm) ; G4 théorème du plancher ; G5 signe de la garde.
Sans la clé, la suite fast est bit-identique (principe VIII).

**Maillage** : `tools/make_cut3d_mesh.py` — entaille avec **jeu** (leçon 2D :
entaille exactement à la passe → 51 joints rompus avant la face verticale, 11
avec 0,284 mm de jeu) et zone fine en **couloir** le long de la course, pas en
boule. Deck de référence : `configs/cut3d_heilman.cfg` (phase 1, sans
confinement). À 0,5 mm : 18 596 tets, dt = 1,54·10⁻⁹ s, 324 000 pas pour 5 mm
de course à 10 m/s.

**Dette** : le confinement 3D est **scalaire** (`confiningPressure`), l'article
demande 30,29 MPa vertical + 17,85 horizontal (phase 2, ~145 lignes) ; les
frontières absorbantes combattent le confinement (jauge 91 → 67 %) et le 3D
n'avertit pas là où le 2D le fait ; les indicateurs de pompe (`injection
outil`, `v nodale max`) n'existent qu'en 2D.

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

### 5.18 Étude « briques constitutives » en fem3d (2026-09-03, `rockim_f2w2.exe` et suivants)

*Chantier `etude_lois_fem/` (proposition, phase A, bancs courts). Toutes les clés sont opt-in ; à clés
absentes, `history.csv` et `frames.csv` des 10 configs fem3d du dépôt sont bit-identiques au build de
référence `rockim_f2w0.exe` (`etude_lois_fem/bitid_w0` vs `bitid_w2`).*

**Briques du noyau `dpr` / `saksala` (`MatLaw.cpp`, classe `PlasticDamageLaw`)**

| clé (défaut) | rôle |
|---|---|
| `meridian = linear \| power` (linear) | méridien courbe pour σ₃ ≥ 0 : en compression triaxiale q = fc0 + `merB`·σ₃^`merN` (merB en MPa^(1−n), σ₃ en MPa ; défauts 56,59 / 0,538 = Red Bohus), fc0 = UCS du cône (2c cos φ/(1−sin φ), surcharge `merFc0`). Lu dans le plan (p, q) par la pseudo-σ₃ = −p − q/3 (exacte en triaxial de compression), résolu en σ₃ par Newton borné + bissection (la pente de σ₃ⁿ est infinie en 0), retour radial déviatorique inchangé. Côté tractif et forme de Lode identiques au cône linéaire ; continu en σ₃ = 0 |
| `compDamage = none \| crackband` (none), `compAc` (0,98), `compGIIc` (1e4 J/m² = 10 N/mm) | endommagement compressif ω_c = A_c (1 − exp(−b_c ε̄ᵖ)), b_c = fc0 h_e/G_IIc (h_e = ∛V₀), sur la partie spectrale négative de la contrainte effective (la brique ω_c de MH 2018, le d_c de CDP). État `MatState::Dc`, champ VTU `omegaC` |
| `dprCap` (false), puis `capP0`, `capH` (K) | avec `dprCap = true`, `dpr` lit `capP0`/`capH` (la compaction volumique de `saksala` sans sa viscosité) : l'ablation du cap change une clé. Derrière une clé d'activation parce que `fem3d_percussion_dpr.cfg` porte `capP0` en héritage du deck saksala et doit rester bit-identique |
| `erodeDc` (0 = off) | troisième canal d'érosion, sur ω_c / A_c NORMALISÉ (0,99 atteignable) ; exige `compDamage = crackband` |
| `dpApex` (false) | ⚠ **défaut du noyau d'origine** : au-delà de l'apex du cône (k − 3αp ≤ 0, traction hydrostatique > 18,8 MPa sur la carte Bohus, atteinte dès 0,03 % de déformation effective en traction) le retour radial à p fixe donne dλ > √J₂/G et le déviateur change de signe ; la part compressive fabriquée échappe au split unilatéral → injection d'énergie (B3 : 10¹⁵⁵ J ; A2b : wPlas −2·10⁴ J). `erodeD = 0,98` masquait le défaut. Avec `dpApex = true` : aucun retour DP au-delà de l'apex (la traction revient au cut-off de Rankine, OPTION-3 de la VUMAT DP-DFH). Promotion en défaut = décision Fernando |
| `dpTension = on \| off` (on) | `off` : le cône DP n'agit qu'en compression (p ≤ 0) — la vraie OPTION-3 de la VUMAT DP-DFH. Entre p = 0 et l'apex, la jambe tractive du cône (23,1 MPa effectifs en uniaxial) coulait dans l'espace effectif pendant que Rankine endommageait : sur la barre E1 cette plasticité tractive dissipait 3,3 Gf·A (revue adverse du 04/09). La référence R de l'étude pose `dpApex = true` et `dpTension = off` |
| `rankineDrive = strain \| stress` (strain) | ⚠ **défaut du noyau d'origine** : le cut-off de Rankine est piloté par la déformation principale maximale (l'en-tête disait « contrainte effective »). La dilatation de Poisson d'un élément comprimé (2νP/E = 7,5·10⁻⁴ sous 100 MPa latéraux, νσ/E = 4,7·10⁻⁴ sous l'UCS) dépasse k0 = ft/E = 1,2·10⁻⁴ et met D ≈ 0,85 partout SANS traction (matrice `C_T1_R_P100` : 100 % du bloc à D ≥ 0,5 avant l'impact). `stress` : κ = σ₁,eff/E — identique en traction uniaxiale (E1). La référence R pose `rankineDrive = stress` |
| `erodeWfrac` (0 = off) | spall sur l'ÉNERGIE de la bande : élément retiré (en traction nette) quand wDamT ≥ erodeWfrac · Gf/lc. ⚠ le seuil hérité `erodeD = 0,98` retire l'élément à κ ≈ 50 k0 où il porte encore 78 % de ft (D = 1 − k0/κ mesure la perte de raideur) : 80 % de Gf jamais dissipés (barre E1) — poser `erodeD = 2` et `erodeWfrac = 0,98` dans l'étude |
| compteurs (toujours calculés, sans effet) | `MatState::wPlas` (∫σ_nom : dε_p, contrainte nominale, cap compris), `wDamT` (∫Y_t dD), `wDamC` (∫Y_c dω_c), Y = ½ σ^± : C⁻¹ σ^± ; `eroCode` 1 spall / 2 broyage / 3 ω_c / 4 dfh |
| `rockim selftest-triax [csv]` | empreinte triaxiale (cibles analytiques à 10⁻⁷ % ; puissance +1,5 % au seul point σ₃ = 0, tangente verticale), deux contrôles qui doivent rater les données (corde −23 % à 50 MPa, pente 3,82 −50 % à 20), crack band en compression (aire adoucie × h = A_c G_IIc à 1,6 %) |

**Solveur `fem3d` (`Fem3dSolver.cpp`)**

| clé (défaut) | rôle |
|---|---|
| `confiningPressure` (0), `confiningRamp` (30e-6 s), `confineFaces = lateral \| all` (lateral), `topPressure` (0), `bottomPressure` (0), `confineGaugeTime` (2 rampes) | pression suiveuse sur les faces extérieures d'origine (copie de `Fdem3dSolver`), rampe cosinus ; `all` ajoute le fond (inerte si encastré) ; `topPressure` = boue sur z = H (faces sous l'outil comprises), `bottomPressure` = appui sous le bloc (équilibre d'une pression de dessus quand le fond est un amortisseur) ; les faces des éléments érodés et les faces y d'une tranche `symmetryY` ne reçoivent rien. Jauge cœur (moitié centrale) ⟨(σ_xx+σ_yy)/2⟩ (⟨σ_xx⟩ en tranche `symmetryY`) latchée à `confineGaugeTime`. ⚠ poser `absorbSpringFactor = 0` (avertissement imprimé sinon) |
| `mesh = file`, `meshFile` | maillage tétraédrique NON structuré (Gmsh MSH 2.2, même lecteur que fdem3d), translaté à l'origine, W/D/H redéfinis depuis la boîte ; hmin = plus petit diamètre inscrit 6V/A (⚠ pour un tet de Kuhn c'est 0,39 h : `dtFactor = 0,7` sur cette base équivaut au 0,3 sur l'arête de la grille) ; lc = ∛V₀. La grille de Kuhn (défaut) est STRUCTURÉE : la règle de la thèse impose `mesh = file` pour tout essai lu en faciès ou en énergie de bande ; générateur `etude_lois_fem/meshes/make_meshes.py` (Delaunay 3D + Netgen, champ de taille autour de l'impact) |
| `toolDelay` (0) | l'outil est gelé (pas de contact, pas d'intégration) tant que t < toolDelay |
| `symmetryY` (false) | tranche plane : v_y = 0 partout, pas de Lysmer ni de confinement sur les faces y |
| `activeNodes` (false) | masque des nœuds portés par ≥ 1 élément vivant (copie du fem 2D) : un nœud orphelin ne voit plus l'outil |
| `erodeDetMin` (0 = off), `erodeStrainMax` (0 = off) | soupape géométrique : det F < seuil (écrasement) ou ‖ε_Biot‖_F > seuil (étirement, le `strainCap` du fem 2D) ⇒ élément retiré et compté à part (`nEroGeo`, `V_eroGeo`, énergie élastique effacée `eRemoved`). Les deux sont nécessaires : la première seule laisse les orphelins étirer leurs voisins en cascade (banc B3) |
| `toolShape = blade` (shear seulement), `backRakeDeg` (20), `clearanceDeg` (10), `bladeHeight` (0,02) | lame rigide plane sur toute l'épaisseur ; `toolX` = la pointe, à `cutDepth` sous la surface ; corps = intersection des deux demi-espaces (coordonnées sr, sc ≥ 0 bornées par `bladeHeight`) ; face qui repousse : coupe au-dessus du niveau de la pointe, dépouille en dessous (la règle « face la plus proche » fait enjamber la matière : B4 à dépouille 60°) ; contact pénalité identique à la sphère. Validé : Fz/Fx = tan(rake) exact à µ = 0 sans rebond sur la dépouille ; ⚠ en tranche élastique la matière non enlevée rebondit sur la dépouille et annule la portance |
| `fieldStats` (false) | colonnes ajoutées EN FIN DE LIGNE de `history.csv` : `keBlock, wPlas, wDamT, wDamC, eRemoved, V_D09, V_wc05, V_eroLaw, V_eroGeo, detFmin, pMin, slatMean, nEroLaw, nEroGeo, nEroSpall, nEroCrush, nEroWc` ; champs VTU `omegaC, detF, erodedBy (1-4 loi, 5 soupape), sigLat`. Avec confinement : `confP, confAch, confWork` (avant les colonnes de champ) |
| `contactPenaltyFactor` (1 = bit-identique) | **raideur du contact outil** (2026-09-04 soir, `rockim_f2w12.exe`) : kp = facteur × E h_min par nœud (`toolContact`, sphère / poinçon / lame). Constat de Fernando sur le quart de bloc Delaunay : h_min est un **sliver** (0,11 mm pour des éléments de 0,5 mm), kp = E h_min = 8,6e6 N/m est 4 à 5 fois plus souple que la raideur nodale de la roche (E × 0,5 mm), interpénétration 0,07 mm en moyenne et 0,43 mm au nœud du pôle — autant que l'indentation. Stabilité **bornée par construction** : `computeStableDt` prend dt = dtFactor × min(h_min/c_P, 2√(m_min/kp)), donc ω_p/(2/dt) = ω_p dt/2 ≤ dtFactor (ω_p = √(kp/m_min), limite du schéma centré ω_p dt < 2, ou 2(√(1+ξ²) − ξ) avec l'amortissement `contactXi`). À l'init, ligne `[FEM3D] tool penalty kp = … omega_p/(2/dt) = …` ; **avertissement** (pas d'exception) si > 0,5 — atteint seulement avec dtFactor > 0,5 (0,7 sur les maillages Gmsh) quand la borne de pénalité domine la CFL. Banc `etude_lois_fem/bitid_w12/contact_kp{1,5,5_dt07}.cfg` (élastique, 20 µs) : facteur 5 → dt 1,41e-7 → 6,30e-8 s, rapport 0,3 = dtFactor, pic de force 4,48 → 7,89 kN ; avec dtFactor 0,7 : rapport 0,7 et `[FEM3D] WARNING`. **Revue (nuit du 4 au 5, même binaire w12 rebâti)** : (i) le ressort seul sous-estime la marge — un nœud de contact porte aussi la raideur d'élément (ω_el = 2 c_P/h_min, ω_el dt/2 = dt/CFL) et la borne de Rayleigh donne ω_tot ≤ √(ω_el² + ω_p²) ; la ligne imprime désormais `omega_el/(2/dt) = dt/CFL` et le **rapport combiné** √((dt/CFL)² + (ω_p dt/2)²), qui déclenche l'avertissement au-delà de 0,5 (à dtFactor ≤ 0,35 il reste < 0,5 quel que soit kp) : kp1 → 0,402, kp5 → 0,323, kp5_dt07 → **0,754** (ressort seul 0,7, élément 0,28) et `[FEM3D] WARNING: combined contact ratio … > 0.5` ; (ii) en `toolContact = signorini` (port de la session parallèle) rien n'est imprimé ni averti (kp inutilisé, vérifié `contact_sig.cfg`) ; (iii) `kpFactor` (clé du mode fem) est accepté comme alias — `contact_kpalias.cfg` (kpFactor = 5) donne le même history.csv que kp5 (6de5f5883d6081ca), `contact_kpconflit.cfg` (5 et 3) est refusé avec un message nommant les deux clés |
| `vtkCap` (false = VTU bit-identiques) | **champs VTU de la zone compactée** (revue w12) : ajoute aux frames les champs cellulaires `capPc` (= `MatState::pc`, Pa ; 0 sans cap) et `epsVpl` (= −tr ε_pl, compaction volumique plastique > 0, dilatance < 0), pour `cdp` (cdpCap) ET `dpr`/`saksala` (capP0) — la zone broyée par le cap n'était dans aucun champ (seul le résumé `max cap pc` la voyait). `history.csv` inchangé (vérifié : contact_kp1 avec/sans la clé, même hachage 70a48f23a409656d) |

Règle apprise sur ce chantier : **après toute modification d'un en-tête partagé (`MatState`), recompiler
TOUTES les unités** — un lien partiel garde le constructeur inline d'un vieux `.obj` (COMDAT) et corrompt
l'état en silence (ω_c lu à 1 avant le premier pas).

### 5.19 Loi `cdp` (Concrete Damaged Plasticity d'Abaqus) et essai triaxial continu (2026-09-04, `rockim_f2w11.exe`)

*Chantier `CONTINUUM/calib_bohus_triax/cdp_rockim/` (journal, bancs, decks matpoint, dossier `revue/`).
Tout est opt-in : à clés absentes, `history.csv`, `frames.csv` **et les `.vtu`** des 10 configs fem3d du
dépôt sont bit-identiques à `rockim_f2w10.exe` au même nombre de threads (`etude_lois_fem/bitid_w11/`,
`w10_omp14/` vs `w11_omp14/` : 9 configs prouvées identiques à OMP 14 au moment du rapport de revue,
`fem3d_shear` (150 s), `fem3d_shear_short` puis 7 configs à OMP 4 (`w10/` vs `w11/`) encore en cours dans
la chaîne de fond — comparer `hashes.txt` à la main quand elle a fini ; `comparaison_omp14.txt`,
`comparaison.txt`). Le nombre de threads change les hachages
(réduction OpenMP par thread) : seule la comparaison au même OMP fait foi. Binaire final : build du 04/09
14:54 (revue). Mise à jour du soir : la chaîne OMP 14 a fini, les **11** configs (`fem3d_shear` et
`fem3d_shear_short` comprises) sont identiques (`comparaison_omp14.txt`) ; le build crack-band de 15:16
(`cdpCompLength`, ci-dessous) est prouvé sur 3 configs à OMP 4 (`w11/hashes.txt` vs `w10/`, `comparaison.txt`).*

**Pourquoi.** La calibration CDP des triaxiaux Red Bohus s'est faite jusqu'ici dans Abaqus (`cdp_pente`,
`cdp_inverse`) ; le port rend la même loi dans rockim, au point matériel (quelques ms par confinement) et en
continu fem3d, pour calibrer sans jeton et pour comparer à armes égales avec les briques de §5.18.

**Conventions.** Tenseurs traction positive ; p̄ = −tr(σ̄)/3 (**compression positive**, comme `pbar` du port
dpdfh) ; q̄ = √(3/2 s̄:s̄) ; μ = G, K = module de compressibilité ; SI. Eigen rend les valeurs propres
croissantes : indice 2 = σ̂_max, indice 0 = σ̂_min.

**Surface, potentiel, retour (classe `CdpLaw`, `MatLaw.cpp`).**

| élément | formule |
|---|---|
| surface (effectif) | F = [q̄ − 3α p̄ + β⟨σ̂_max⟩ − γ⟨−σ̂_max⟩]/(1−α) − σ̄_c ; α = (fb0/fc0 − 1)/(2 fb0/fc0 − 1), γ = 3(1−Kc)/(2Kc−1), β = σ̄_c/σ̄_t (1−α) − (1+α) |
| pic triaxial de compression | q_pic = fc0_pic + m σ₃, m = (3α+γ)/(1−α) ; carte historique : m = 3,817 (126,6 / 202,9 / 317,5 / 412,9 / 508,3 MPa à 0/20/50/75/100) ; carte inverse (Kc 0,6061) : m = 6,751 (316,8 / 451,8 / 654,3 / 823,1 / 991,9). En nominal (r = 0, d = d_c, confinement piloté nominal) **toute** la courbe est la table σ_c(ε_c^pl) translatée de m σ₃ |
| potentiel (non associé) | G = √((ecc ft tan ψ)² + q̄²) − p̄ tan ψ ; ∂G/∂σ̄ = 1,5 s̄/√(a²+q̄²) + tan ψ/3 I. Hyperbolique, pas Mohr-Coulomb : ψ = 35° ⇒ −dε_v/dε₁ = tan ψ/(ρ − tan ψ/3) = 0,913 (ψ_MC équivalent 18,3°) |
| retour spectral, implicite | p̄ = p̄_tr + K tan ψ λ (**signe +** : la dilatance fait monter la pression effective, comme `pnew = pbar + xK tanp dlam` du port dpdfh) ; q̄(1 + 3μλ/√(a²+q̄²)) = q̄_tr (Newton monotone depuis q₀ = max(q̄_tr − 3μλ, 0), crochet [q₀, q̄_tr]) ; dε̂_i = λ[1,5 s̄_i/√(a²+q̄²) + tan ψ/3] ; r, β, σ̄_c à n+1 (Lee-Fenves 2001) ; F(λ) n'est que C0 → crochet par doublement puis regula falsi Illinois (bissection de garde), jamais de Newton pur ; état hydrostatique tractif (q̄_tr = 0) : retour purement volumique σ̂ = −p̄, P = (1−α)σ̄_c/(3α+β) |
| écrouissages | dε_t^pl = max(0, r dε̂_max), dε_c^pl = max(0, −(1−r) dε̂_min), r = Σ⟨σ̂⟩/Σ\|σ̂\| (clamp : dε̂_min > 0 près de l'axe hydrostatique ; ψ < 56,3° requis pour l'équibiaxial) ; en compression uniaxiale dε_c^pl = λ(ρ − tan ψ/3) = \|dε₁^pl\| (la variable de la table Abaqus) |
| endommagement (à chaque appel) | d = 1 − (1 − s_t d_c)(1 − s_c d_t), s_t = 1 − w_t r, s_c = 1 − w_c(1−r) ; σ_nom = (1−d) σ̄ ; d dépend de r, donc de la contrainte finale : ce n'est **pas** une variable d'état (un élément fissuré recomprimé affiche d = d_c) ; plafond 0,99 |

**Tables (parité Abaqus, « lecture B »).** Compression : `cdpHardening` (σ_c : ε_in) et `cdpCompDamage`
(d_c : ε_in), abscisses **fusionnées** (Abaqus accepte des abscisses différentes), chaque nœud converti à
l'init ε_c^pl = ε_in − d_c/(1−d_c) σ_c/E0, suite strictement croissante exigée (exception nommant le
point) ; σ_c, d_c linéaires en ε_c^pl, constantes après le dernier nœud ; σ̄_c = σ_c/(1−d_c). Carte
historique convertie : ε_c^pl = 0 / 8e-4 / 3,227e-3 / 1,0519e-2 → σ_c 50,6 / 126,6 / 60 / 10 MPa, σ̄_c 50,6 /
126,6 / 120 / 125 (l'adoucissement effectif est quasi nul : toute la chute nominale vient de d_c).
Traction : `cdpTension = gfi` (défaut rockim, = les decks du dépôt : σ_t = ft(1 − u/u_t0), u_t0 = 2Gf/ft =
2,353e-5 m) ou `table` (= type DISPLACEMENT, `cdpTensionTable` « σ:u_ck », σ_t(0) = ft à 1e-6 près) ; le
type STRAIN n'est **pas** couvert ; `cdpTensionDamage` (d_t : u_ck). Nœuds fusionnés {0} ∪ {u de σ_t} ∪ {u de
d_t} ∪ {u_t0} (≤ 16), c_k = d_t,k/(1−d_t,k) σ_t,k/E0 précalculés, conversion **à la volée par élément**
ε_t^pl,k = u_k/lc − c_k, interpolation linéaire **en ε_t^pl** (la lecture A — interpoler en u natif — rendrait
Gf exact et serait une clé opt-in ultérieure `cdpTensionInterp = u`, hors périmètre ; écart < 0,4 % à 1 mm).
Convention « lc partout » = Abaqus avec l0 = lc (identique aux decks mm à < 1 % ; en SI l0 = 1 m vaudrait
1000 lc). Plancher σ_t ≥ 1e-3 ft appliqué au **nominal** avant division par (1−d_t) (σ̄_t,res = 170 kPa,
β ≈ 650). lc = ∛V₀ du solveur (longueur caractéristique Abaqus des tétraèdres, V^(1/3) de mémoire du manuel).

**Énergie (banc (c)).** En traction uniaxiale, lecture B : W_tot × lc − lc × wDamT = Gf (identité exacte à
O(c_k), mesurée −0,007 % à 1 mm, −0,013 % à 2 mm) ; l'aire totale × lc vaut Gf + lc wDamT (+0,36 % à 1 mm,
+0,71 % à 2 mm) ; wDamT = ∫Y dd_t = 361 J/m³ indépendant de lc. La dissipation tractive est **plastique**
(ε_t^pl ≈ u_ck/lc) : le canal `erodeWfrac` est inerte en CDP (gardé, documenté).

**Gardes à l'init (lcMax du solveur, exception).** (G1) lcMax ≤ 2E0Gf/ft² (gfi ; 0,215 m ici) ou max
|Δσ_t/Δu| ≤ E0/lcMax (table) ; (G2) nœuds convertis croissants : u_{k+1} − u_k > lcMax (c_{k+1} − c_k) ;
(G3) snap-back **effectif** : σ̄_t = σ_t/(1−d_t) n'étant pas affine sur un segment, sa pente locale
[B(1−d) + Dσ]/(1−d)² est évaluée aux deux bouts de chaque segment et doit rester ≤ 0,9 (2μ + K tan ψ)/(1 +
tan ψ/3) = 75,4 GPa (Bohus) ; en compression ≤ 0,9 (3μ + 3αK tan ψ)/((1−α)(1 − tan ψ/3)) = 141,6 GPa.
Pente effective finale de la carte GFI + d_t 0,95 : 7,05 GPa × lc/mm → **lc < 10,7 mm** (banc (j) : 20 mm
lève l'exception, 1 mm passe). Bornes démontrées pour les chemins uniaxiaux ; le facteur 0,9 est prudent,
non prouvé suffisant hors uniaxial (une exception « no bracket / return mapping failed » sur un chemin mixte
signalerait ce point).

**Point aveugle (parité, pas une erreur).** L'adoucissement compressif est en déformation, non régularisé en
maillage ; remettre la carte à l'échelle en contrainte sans toucher aux abscisses raidit la chute nominale
post-pic entre les nœuds 2 et 3 : dσ/dε_total = −(s₂−s₃)/[(ε_in3−ε_in2) + (s₃−s₂)/E0] = −28,4 GPa (−0,37 E,
historique) contre −158 GPa (−2,0 E, carte inverse ×2,5027). À l'échelle de l'éprouvette (bande d'un élément
contre corps élastique) c'est un snap-back structurel : la « chute verticale » Abaqus (391,8 → 19,2 MPa en un
cadre) est attendue **à l'identique** dans rockim. Une régularisation (`cdpCompBand`) serait une clé opt-in
ultérieure. *(Soir du 2026-09-04 : cette clé existe, elle s'appelle `cdpCompLength` — voir « Crack band
en compression » ci-dessous ; à clé absente le point aveugle est conservé tel quel, parité oblige.)*

**Revue du 2026-09-04 (constats traités, tous opt-in ou sans effet à clés absentes).**

1. *Hétérogénéité* : `CdpLaw` ignorait `MatState::ftScale` (le facteur de Weibull que `Fem3dSolver`
   tire pour toute loi dès que `matWeibullM > 0`) — un deck `law = cdp` Weibull tournait homogène en
   silence. Désormais honoré comme par `dpr` (E1) : par élément s = ftScale, ft_loc = s·ft, a = ecc·ft_loc·tan ψ,
   table de traction σ_t × s et u_ck × κ (κ = 1/s en `weibullScope = strength` — Gf conservé — ou 1 en
   `strengthGf` — Gf × s), d_t(u) suit la même dilatation des abscisses ; la table de compression n'est pas
   touchée (comme la cohésion de dpr). Les gardes G1-G3 dépendent de (lc, s) : vérifiées à l'init à
   (lcMax, 1) et **re-vérifiées au premier appel de tout élément à s ≠ 1** (exception nommant lc et
   ftScale) ; G3 se durcit en s² (carte Bohus : lc < 10,7 mm / s²). À s = 1 l'arithmétique est identique
   (× 1,0 exact). Banc (l) : pic de traction 1,3 ft à −0,01 %, W lc − lc wDamT = Gf_loc (100 / 130), pic
   de compression inchangé, G3 levée au premier appel à s = 4 sur 1 mm. Deck `revue/fem3d_cdp_weibull_guard.cfg`
   (2×2×4, m = 3, OMP 14) : « `[fem3d] constitutive law threw in element 6 (lc = 0.005503 m, ftScale =
   1.657457) : cdp (G3) …` », code retour 1 ; sans `matWeibullM` le même deck tourne.
2. *Exceptions dans la région OpenMP* (`Fem3dSolver::elementForces`) : une exception levée par la loi
   (`cdp: no bracket / return mapping failed`, `crack-band` de dpr à ftScale > 1…) appelait
   `std::terminate` — mort sans message. Chaque thread capture la première exception de ses éléments
   (message, indice, lc, ftScale) dans son accumulateur et elle est relancée après la barrière avec
   l'en-tête `[fem3d] constitutive law threw in element N (lc = …, ftScale = …)`. Sans exception : aucun
   chemin nouveau (bit-identité).
3. *Viscosité* : le constat « les decks portent μ = 5e-5 s, rockim ne reproduit pas la surcontrainte
   visqueuse (+0,7-0,9 %) » est **rejeté sur le fond** : les decks sont `*Dynamic, Explicit` et
   **Abaqus/Explicit ignore le paramètre de viscosité** de `*CONCRETE DAMAGED PLASTICITY` (Keywords
   Reference, 5ᵉ donnée : « used … in Abaqus/Standard analyses. This parameter is ignored in
   Abaqus/Explicit ») — la parité avec les runs de référence est le défaut rate-indépendant. La clé opt-in
   `cdpViscosity` (s) est néanmoins fournie (régularisation de Duvaut-Lions de la CDP d'Abaqus/Standard :
   ε_pl,v += dt/(μ+dt)(ε_pl − ε_pl,v), d_v idem, σ = (1−d_v) D0 (ε − ε_pl,v) ; dt > 0 requis). Banc (k) :
   à μ = 5e-5 s et 0,75/s la surcontrainte stationnaire en compression uniaxiale n'est **pas** E0 μ ε̇ =
   2,9 MPa mais 18,1 MPa (formule exacte (m+1)c + μλ̇|(D0 n)_ax|, c = μλ̇(D0 n)_lat : la latérale visqueuse
   doit être compensée par un confinement inviscide que la pente m amplifie ; vérifiée à 2e-9 %, linéaire
   en vitesse, μ = 1e3 s → réponse élastique) ; sur le triaxial à 20 MPa de la carte historique à
   0,748/s : +9,9 MPa (+4,9 % du pic, 2,4 % de 404,8). Le chiffre du constat (0,7-0,9 %) était donc faux
   d'un facteur 6 même dans l'hypothèse Standard.
4. *Clés inconnues* : toute clé commençant par `cdp` absente de la liste connue est refusée à l'init
   (exception nommant la clé et la liste) — `Config` ignore les inconnues en silence. Banc (l) : `cdpKC`.
5. *Pilote latéral `Hold1D`* (bancs et matpoint, pas les solveurs) : le constat était réel et une
   première correction (crochet par le signe + bissection « après deux essais sans progrès ») **ne
   changeait rien** aux trois cas du constat (3 / 653 / 2947 pas non convergés, mêmes résidus) — une
   sécante qui progresse d'un iota remettait le compteur à zéro, et la recherche du signe ne doublait que
   dans la direction de la compliance, fausse quand le résidu n'est pas monotone. Sonde (balayage de
   2000 points à état restauré) : dans les trois cas le résidu a **une racine franche** (saut < 1e-8 Pa) mais
   n'est **pas monotone** en traction post-pic (1150 montées / 850 descentes ; monotone en équibiaxial).
   Algorithme actuel : (A) crochet par recherche **symétrique** en doublement depuis le pas de compliance
   endommagée (x₀ ± |dx| 2ᵏ), (B) regula falsi d'Illinois dans le crochet avec bissection de garde une
   itération sur trois, arrêt à |résidu| < 1e-3 Pa ou crochet réduit à 1e-16 relatif (discontinuité :
   échec honnête) ; l'état rendu est celui du meilleur point. Résultat (`revue/*.cfg`) : `bigsteps_ten`
   (10 pas de 3e-3 en traction) 0 non convergé, 311 appels, σ_ax = 8 500 Pa exact à ε = 0,03 (avant :
   76 kPa) ; `tiny_lc` (lc = 10 µm) 0/3000 ; `psi55_biax` (ψ 55,9°, équibiaxial) 0/3000 ; decks de
   calibration 0/20 000 (≈ 18 appels par pas, ≈ 90 ms par confinement).
6. *Verdict falsifiant du pilotage* : `selftest-cdp` rend 1 si un pas n'a pas convergé (verdict
   « pilotage latéral : 0 pas non convergé ») ; `matpoint` écrit une colonne `residu` (|résidu final| par
   ligne, Pa) et rend **2** si un pas n'a pas convergé (`[matpoint] WARN`), pour qu'une boucle Python le voie.
7. *Banc (j)* : le message des exceptions est inspecté — 20 mm doit lever « (G3) », 300 mm « (G1) »
   (> 2E Gf/ft² = 215 mm), une table d_t = 0,9 à u = 1 µm à lcMax = 2 mm « (G2) » (nœud converti
   décroissant), 1 mm passe.
8. *Plancher σ_t ≥ 1e-3 ft, artefact documenté* : une fois le dernier nœud de traction atteint, ρ =
   q̄/√(a²+q̄²) ≈ 0,3 et le potentiel hyperbolique coule dans les **trois** directions (dε̂_lat =
   λ(−ρ/2 + tan ψ/3) > 0 dès que ρ < 2 tan ψ/3 = 0,467) : l'élément fissuré gonfle latéralement
   (ε_lat remonte de −4,5e-3 à +1,9e-3 entre ε_ax 0,021 et 0,06 à 1 mm), sa déformation volumique
   plastique croît sans borne et wPlas dépasse Gf/lc. Parité Abaqus (même potentiel, même résidu de
   table) ; en fem3d un élément fissuré pousse ses voisins : garde-fou = `erodeD ≤ d_t,max` (0,95), à
   poser explicitement.
9. *Preuve de bit-identité* : les `.vtu` sont hachés (hachage combiné trié) avant leur suppression ;
   `fem3d_shear` rejouée à OMP 14 (150 s) et une copie courte `bitid_w11/fem3d_shear_short.cfg`
   (T 5e-4, 8 cadres, 3762 pas, mêmes chemins de code) à OMP 4 et 14 ; `BITID_OMP` / `BITID_TIMEOUT`
   dans `bitid_w11.sh` (défauts 4 / 185 s = comportement d'origine, sorties dans `<tag>_omp<N>/`).
10. *Logs bruts* : les `.log` matpoint déposés sont la sortie brute (l'avertissement « canal spall INERTE »
    en tête) ; `*_avant_revue.*` conservent les sorties du binaire de 13:39.

**Crack band en compression (`cdpCompLength`, 2026-09-04 soir, build 15:16 — opt-in, 0 = absent =
bit-identique).** Décision de Fernando (« le plus pertinent physiquement ») : la chute post-pic d'un
triaxial est **structurelle** (bande localisée), pas une propriété du point matériel ; la table CDP de
compression étant en déformation et non régularisée, l'énergie dissipée dans la bande dépendait de la
taille d'élément (le point aveugle ci-dessus). Avec `cdpCompLength = L_ref > 0` (m) :

- les abscisses ε_in des tables de compression (`cdpHardening` **et** `cdpCompDamage`, nœuds fusionnés)
  sont lues comme définies à la longueur de référence L_ref, et **seule la branche post-pic** est remise
  à l'échelle de l'élément (pic = premier nœud fusionné où σ_c nominal est maximal, ε_in,pic = 8e-4 sur
  la carte historique) : ε_in,loc = ε_in,pic + (ε_in − ε_in,pic) × L_ref/lc ; la branche pré-pic
  (écrouissage homogène) n'est pas touchée ;
- la conversion d'Abaqus s'applique **après** : ε_c^pl,k(lc) = ε_in,loc,k − c_k, c_k = d_c,k/(1−d_c,k) σ_c,k/E0
  (c_k ne dépend pas de lc) — faite à la volée par élément comme pour la traction (`CdpLaw::compNode`,
  `compState(ε_pl, lc)`) ; à clé absente `compState` lit les nœuds convertis à l'init (arithmétique
  identique : les 3 configs de bit-identité et `cdpCompLength = 0` explicite le prouvent) ;
- conséquence : le déplacement inélastique post-pic u_in = (ε_in,loc − ε_in,pic) lc = (ε_in − ε_in,pic) L_ref
  est **invariant**, donc G_c = ∫(σ_c − σ_res) du_in est la même quelle que soit la maille — crack band de
  Bažant-Oh en compression, même logique que `compGIIc` de la brique ω_c (§5.18). Sur la carte historique à
  L_ref = 2 mm : G_c = 933,1 J/m² (trapèze sur la table), 932,6 (quadrature de la lecture B, dε_in = dε_pl + dc) ;
- **gardes à l'init à lcMax** (la plus grande maille a les segments les plus courts ; compression
  indépendante de ftScale), levées seulement si la clé est > 0 : (G2c) nœuds convertis strictement
  croissants (carte historique, L_ref 1 mm, lc 8 mm : nœud 2 à 4,27e-4 < 8e-4 → exception) ; (G3) pente
  effective ≤ limC (141,6 GPa, existence du retour) ; **(G4) snap-back au point matériel** : la déformation
  totale ε = ε_c^pl + σ̄_c/E0 = ε_in,loc + σ_c/E0 doit croître avec ε_c^pl, soit −dσ̄_c/dε_c^pl ≤ E0
  localement (aux deux bouts de chaque segment, σ̄ = σ/(1−d) n'étant pas affine ; en corde
  −Δσ_c/Δε_in,loc ≤ E0, la condition d'Abaqus sur sa table) ; message nommant la pente, la corde, lc et
  L_ref, et la maille approximative à atteindre. G4 (E0) est plus stricte que G3 (limC ≈ 1,8 E0) : entre les
  deux le retour converge mais la réponse en déformation totale saute. **Écart au libellé de la spec** : le
  critère « −dσ_c/dε_pl ≤ E0 » (nominal sur la plastique) n'est pas le snap-back — il vaut 80,5 GPa > E0 à
  lc = 4 mm / L_ref 2 mm sur la carte historique alors que la corde en ε_in vaut 41,6 GPa et qu'il n'y a
  aucun retournement ; le banc à 1/2/4 mm demandé ne passerait pas — G4 implémente le critère exact.

Banc (m) de `selftest-cdp` (compression uniaxiale, carte historique, L_ref = 2 mm, lc = 1 / 2 / 4 mm ;
le point matériel restitue exactement la table, σ_ax = σ_c(ε_c^pl), ε_ax = ε_in,loc + σ/E0) : aire
adoucie post-pic [∫(σ − σ_res) dε + ((σ_pic − σ_res)² − (σ_end − σ_res)²)/(2E0)] × lc = 932,9 / 932,7 /
932,2 J/m² (= G_c lecture B à +0,03 / +0,002 / −0,05 %, trapèze brut à −0,02 / −0,05 / −0,10 %, tol 1 %),
invariance 0,99975 / 0,99925 (tol 1 %), branche pré-pic identique entre les trois lc à 4,3e-12 (2430 pas,
critère 1e-9) ; **doit échouer** : sans la clé l'aire × lc vaut 466,3 (1 mm) / 1865,3 (4 mm), rapport
4,000 (tol 2 %), et à 1 mm exactement G_c(L_ref)/2 (table lue à lc) ; `cdpCompLength = 0` explicite
bit-identique à la clé absente ; G4 levée sur une carte 126,6 → 10 MPa en 4e-4 (corde 291 GPa) à
L_ref = lc = 2 mm, la même carte passe à lc = 0,5 mm (segment étiré ×4 : 72,9 GPa < E0) ; G2c levée ;
`cdpCompLength < 0` refusé. `matpoint` : colonne `eps_in_c` (ε_c^pl + d_c/(1−d_c) σ_c/E0, table lue à
`mpLc`) en fin de ligne, u_in = (eps_in_c − ε_in,pic) × mpLc ; decks `cdp_rockim/matpoint_cdp_hist_compband_lc1.cfg`
/ `_lc4.cfg` (σ₃ = 0 et 20 MPa, mpStrainMax 0,05) et `check_compband.py` : u_in à l'**épuisement de la
table** 2,240e-5 m aux deux mailles (= (1,2e-2 − 8e-4) L_ref ; au-delà le plateau σ_res n'est pas de
l'adoucissement), G_c 932,9 / 932,2 J/m² aux deux confinements (identité de translation q − mσ₃ = σ_c),
branche pré-pic identique à 0, q_res = 10 + 3,817 × 20 = 86,35 MPa. Piège vu au passage : sous confinement
la déformation élastique **effective** σ/(1−d) atteint ≈ 0,017 à d_c = 0,92, la table à 1 mm (post-pic × 2)
n'est pas épuisée à ε = 0,03 — d'où 0,05. Ce que la calibration doit retenir : avec `cdpCompLength` la carte de compression se calibre
**à une longueur de bande** (L_ref = largeur physique de la bande de cisaillement de l'éprouvette, à
choisir — quelques mm sur un granite), et la pente post-pic vue par un maillage fem3d n'est plus un
artefact de maille ; la garde G4 dit jusqu'où l'on peut grossir les éléments sur une carte donnée.

**Calibration Red Bohus avec `cdpCompLength` (tâche 2, `CONTINUUM/calib_bohus_triax/cdp_rockim/calib_cdp_rockim.py`,
2026-09-04 — aucun code modifié, notes d'usage).** Partition homogène / structurel : Kc, ψ et la table pré-pic
(3 nœuds, plateau pendant l'ajustement) au point matériel (`matpoint`, 8 `least_squares` en parallèle, 84 s) ;
branche post-pic **construite** à L_ref = 2 mm (descente linéaire fc0 → σ_res sur u₁ = 2G/(fc0 − σ_res), d_c
linéaire jusqu'à 0,9). Carte « confinés » : Kc 0,6383, ψ 51,4°, σ_c = 112,2 / 257,6 @ 4,2e-4 / 328,2 @ 1,32e-3 /
132,2 @ 0,194 MPa, `cdpCompLength = 0.002` (`carte_confines_Wexc.cfg`) ; pics à ±5 % sur 20-100 MPa, UCS +159 %
(la droite CDP). Deux points d'usage mesurés sur ces cartes : (1) **l'énergie de bande n'est invariante en lc
qu'à O(c_res/Δε_in(lc))** — u_in = (eps_in_c − ε_in,pic) lc l'est à 1e-6 mm, mais G = ∫(q − q_res) du_in lu sur
la courbe dévie de −2 / −4 / −9 % à lc = 1 / 2 / 4 mm quand d_c,res = 0,9 (c_res = d/(1−d) σ_res/E = 1,5 % contre
Δε_in = u₁/lc), parce que la lecture B (linéaire en ε_pl) rend σ_c(ε_in) non affine ; à d_c,res = 0,5 la dérive
tombe à +0,05 / +0,1 / +0,2 % (banc (m) : 0,05 % avec c_res = 0,11 %) ; (2) sous confinement nominal tenu, la
contrainte effective vaut σ/(1 − d) : à d_c = 0,9 l'élément de bande porte 4-11 % de déformation élastique
dégradée (recouvrable, ∝ lc), qui s'ajoute à u_in dans la réponse d'éprouvette — d_c n'est donc pas neutre
même sans cycles. Un `fc0_pic` effectif lu sur une table à d_c,res = 0,9 vaut 10 σ_res, pas le pic nominal
(le constructeur et `[cdp]` l'affichent ainsi ; la formule q_pic = fc0_pic + mσ₃ du banc (a) s'entend à d_c = 0
au pic).

**Clés (`law = cdp`, défauts = carte historique Red Bohus des decks)**

| clé (défaut) | rôle / validation |
|---|---|
| `cdpDilationDeg` (35) | ψ ; 0 < ψ < 56 exigé (ψ > 0 pour le retour hydrostatique ; ψ ≥ 56,3° dégénère l'écrouissage équibiaxial — borne rockim plus stricte qu'Abaqus, sans effet sur la carte) |
| `cdpEcc` (0,1) | excentricité ; > 0 (sommet lisse, aucun cas apex) |
| `cdpFbFc` (1,16) | fb0/fc0 > 1 (0 < α < 0,5) |
| `cdpKc` (0,667) | 0,5 < Kc ≤ 1 ; Kc = rapport q_TE/q_TC à p̄ fixé, exact tant que σ̂_max < 0 sur les deux méridiens (banc (h) : 0,667002 ; Kc = 1 → 1) |
| `cdpWt` (0), `cdpWc` (1) | facteurs de récupération de raideur, dans [0, 1] |
| `cdpHardening` (« 50.6e6:0 126.6e6:0.0008 60e6:0.004 10e6:0.012 ») | σ_c [Pa] : ε_in, σ > 0, abscisses strictement croissantes, première = 0 |
| `cdpCompDamage` (« 0:0 0:0.0008 0.5:0.004 0.92:0.012 ») | d_c : ε_in, 0 ≤ d ≤ 0,99, première abscisse 0 avec d = 0 |
| `cdpTension` (gfi), `cdpTensionTable` (« σ:u_ck … », requis si table), `cdpTensionDamage` (« 0:0 0.95:2.35e-5 ») | voir Tables ; ft, Gf, E du bloc Material |
| `erodeD` (0,98), `erodeEpv` (1,5), `erodeDc` (0), `erodeWfrac` (0) | lus comme pour dpr : spall = d_t ≥ erodeD en traction nette (**avertissement** si erodeD > d_t,max : « canal spall INERTE ; poser erodeD ≤ 0,95 ») ; erodeEpv sur ε_c^pl ; erodeDc sur d_c/d_c,max normalisé (garde `compDamage` levée pour cdp) ; erodeWfrac inerte ; `eroCode` 1/2/3 |
| champs partagés | `MatState::D` = d_t (monotone : VTU `damage`, `V_D09`, canal spall), `Dc` = d_c, `epvEq` = ε_c^pl, `kappa` = ε_t^pl (VTU `kapDP`), sous-état `MatState::Cdp {epsTpl, epsCpl, d}` (d total, diagnostic) ; `wPlas`, `wDamT`, `wDamC` (incréments de d_t, d_c, pas de d) |
| syntaxe des tables | « v:x v:x … » séparés par des espaces, lecture stricte par jeton (virgule décimale ou jeton sans ':' refusés, message nommant clé et jeton) ; `Config` conserve les espaces internes |
| `cdpViscosity` (0 = rate-indépendant, bit-identique) | μ [s] de la régularisation de Duvaut-Lions (CDP d'Abaqus/**Standard** ; Abaqus/Explicit ignore ce paramètre, les decks de référence sont donc inviscides) ; dt > 0 requis (`mpDt` = Δε/ε̇ en matpoint) ; sous-état `MatState::Cdp {epsPv, dv}` |
| `matWeibullM`, `strengthCorrLength`, `weibullScope` (clés du solveur) | ftScale **honoré** par élément (ft, a, table de traction ; compression intacte) ; gardes G1-G3 re-vérifiées au premier appel à ftScale ≠ 1 (`MatState::Cdp::guarded`) |
| clé `cdp*` inconnue | **refusée** à l'init (exception nommant la clé et la liste connue) |
| `cdpCompLength` (0 = absent, bit-identique) | L_ref [m] du crack band en compression : abscisses ε_in des tables de compression définies à L_ref, branche **post-pic** remise à l'échelle L_ref/lc par élément (pré-pic intacte), conversion ε_in → ε_pl après ; u_in et G_c invariants en lc ; gardes G2c / G3 / **G4 (snap-back au point matériel, −dσ̄_c/dε_c^pl ≤ E0)** à lcMax ; ≥ 0 exigé |
| `cdpCap` (false = absent, bit-identique), `cdpCapP0` (obligatoire si cdpCap, Pa), `cdpCapH` (K = E/(3(1−2ν))) | **cap volumique de compaction** (`rockim_f2w12.exe`, voir ci-dessous) : sur la pression **effective** du prédicteur p̄ = −tr(σ̄_tr)/3, si p̄ > pc alors dev = (p̄ − pc)/(K + H), ε_pl −= dev/3 I, pc += H dev, σ̄_tr += K dev I, **avant** le retour de Lubliner ; pc initialisée à cdpCapP0 au premier appel (état partagé `MatState::pc`, résumé `max cap pc`) ; ne touche ni ε_t^pl/ε_c^pl ni d_t/d_c ; compté dans wPlas ; `cdpCapP0`/`cdpCapH` sans `cdpCap = true` refusés. **Revue** : le cap est imposé au **prédicteur seulement** (le retour de Lubliner qui suit est dilatant : p̄ de fin de pas = pc + K tan ψ Δλ, recappé au pas suivant) ; le seuil **nominal** vaut (1−d) pc (0,08 pc0 à d_c = 0,92) — `cdpCapP0` se calibre en pression effective ; la compaction seule n'érode jamais (`erodeEpv` lit ε_c^pl) ; champs VTU par `vtkCap = true` (§5.18) |

dt est ignoré à `cdpViscosity = 0` (défaut, = les runs Abaqus/Explicit de référence).

**Bancs (`rockim selftest-cdp [csv]`, code retour 0 seulement si tout passe ; log
`cdp_rockim/selftests/selftest_cdp_w11.log`).** (a) pics triaxiaux 0/20/50/75/100 = fc0_pic + m σ₃ à
−0,006 % (tol 0,2 %) et identité de translation sur toute la branche plastique à 4e-9 fc0 ; (b) équibiaxial
146,84 vs fb0 146,86 ; (c) traction : pic ft à −0,002 %, aire × lc = Gf à +0,35 % (1 mm) / +0,71 % (2 mm),
identité exacte à −0,007 / −0,013 %, ε_t^pl = ε^pl_zz à 5e-11 (critère 1e-8 ; r = 1 − O(résidu latéral/ft)) ; (d) décharges
(1−d_c)E0, (1−d_t)E0, E0 après refermeture à 1e-7 % ; (e) dilatance 0,9135 (ψ 35), 1,9777 (ψ 50), 0,9266
(ecc 0,5 à q̄ = 20 MPa) à 1e-12 ; (f) **doit rater** : −49,87 % à 20 MPa, −36,40 % à 100 ; (g) carte inverse
q(0) = 316,8, q(20) = 451,8 (+11,6 % sur 404,8 : 404,8 n'est atteint que par le rabattement structurel
r = 0,896 de l'éprouvette EF) ; (h) q_TE/q_TC = Kc ; (i) compression hydrostatique jamais plastique, traction
hydrostatique isotrope avec σ = (1−d)P à 1e-15, p̄(n+1) − p̄_tr = K tr(dε^p) > 0 ; (j) messages « (G3) » à
20 mm, « (G1) » à 300 mm, « (G2) » sur un nœud converti décroissant, 1 mm passe ; (k) `cdpViscosity` = 0
explicite bit-identique à la clé absente, surcontrainte de Duvaut-Lions exacte (18,10 MPa à μ 5e-5 s /
0,75/s), linéaire en vitesse, μ = 1e3 s → élastique, +4,9 % sur le triaxial 20 MPa si μ était actif ;
(l) ftScale 1,3 (deux portées), G3 par élément à ftScale 4, clé `cdpKC` refusée ; (m) crack band en
compression (`cdpCompLength` 2 mm à lc 1/2/4 mm : G_c invariant, pré-pic identique, sans la clé rapport 4,
gardes G4 / G2c) — 19 contrôles ajoutés le soir. 79 contrôles [OK] (60 à la revue) ; (n) cap volumique
`cdpCap` (ci-dessous) et (o) cap + cône de Lubliner — 15 contrôles le soir, **94 [OK]** ; revue de la nuit du 4 au 5 :
(n) cdpCapH absent = K et (o) refondu en trois variantes o1/o2/o3 (ci-dessous), **109 [OK] / 0 [FAIL]** avec le
`rockim_f2w12.exe` rebâti (`selftests/selftest_cdp_w12.log` du 05/09 00:06, 8,56 M appels, 0 pas non convergé).

**Cap volumique de compaction (`cdpCap`, 2026-09-04 soir, `rockim_f2w12.exe` — opt-in, false = absent =
bit-identique).** Constat de Fernando sur la percussion quart de bloc (`perc3d/`) : la CDP n'a **pas** de cap
de compaction et son méridien est une droite — sous l'insert (pression de contact ~ 8 GPa) la roche reste
élastique, aucune zone broyée, rebond de Hertz. Même brique que le cap de `dpr`/`saksala`
(`PlasticDamageLaw::stress`, « pressure cap », clés `capP0`/`capH`, activation `dprCap`), posée dans
`CdpLaw::stress` **entre le prédicteur élastique et le retour de Lubliner**, sur la pression **effective**
p̄ = −tr(σ̄_tr)/3 (compression positive) : si p̄ > pc (pc initialisée à `cdpCapP0` au premier appel, état
partagé `MatState::pc` — celui que le résumé `max cap pc` lit déjà), retour volumique dev = (p̄ − pc)/(K + H) > 0,
ε_pl −= dev/3 I (la trace de ε_pl **diminue** : compaction), pc += H dev, σ̄_tr += K dev I (p̄ redescend à pc + H dev,
avec H = `cdpCapH`, défaut K) ; le retour de Lubliner s'applique ensuite sur le prédicteur corrigé (décomposition
spectrale faite après le cap). **Choix documenté** : le cap ne touche ni ε_t^pl/ε_c^pl ni d_t/d_c — comme dans
`dpr`, la compaction a son propre écrouissage (pc) et n'alimente pas l'endommagement de la table de Lee-Fenves
(l'écrasement sous l'insert est une densification de pores, pas une fissuration ; un élément cappé garde sa
raideur nominale). **Conséquence (revue)** : la compaction seule n'érode **jamais** — `erodeEpv` lit ε_c^pl et le
canal crush aussi, que le cap ne nourrit pas ; en percussion c'est `erodeDetMin` (0,3 dans `perc3d/`) qui retire
un élément écrasé, et un seuil opt-in sur ε_v^pl reste à ajouter si Fernando veut un broyage « par compaction ». L'incrément de compaction est compté
dans `wPlas` (σ_nom : dε_pl, terme séparé pour laisser la ligne d'origine intacte). Aucun champ VTU ne portait pc ni ε_v^pl
(la zone compactée était invisible dans ParaView, seul le résumé `max cap pc` la voyait) : depuis la revue,
`vtkCap = true` (Fem3dSolver, opt-in, §5.18) ajoute les champs cellulaires `capPc` et `epsVpl` pour cdp ET dpr. Clés orphelines (`cdpCapP0`
ou `cdpCapH` sans `cdpCap = true`) **refusées** (rien d'appliqué en silence) ; `cdpCap = true` sans `cdpCapP0`
refusé ; clés ajoutées à la liste des `cdp*` connues. Réponse en compression hydrostatique : p = K|ε_v| jusqu'à
pc0, puis p = (K pc0 + K H |ε_v|)/(K + H), pente K H/(K + H) (= K/2 à H = K).

**Deux précisions de la revue (nuit du 4 au 5 septembre).** (i) Le cap est imposé au **prédicteur seulement**.
Contrairement à `dpr`, dont le retour du cône est purement déviatorique (p intact, cap exact en fin de pas), le
retour de Lubliner est **dilatant** (p̄ = p̄_tr + K tan ψ λ, signe +) : quand cap et cône agissent dans le même pas,
p̄ de fin de pas = pc + K tan ψ Δλ = pc + K (Δtr ε_pl + Δpc/H) > pc. Le dépassement est borné et recappé au
prédicteur suivant (lag d'un pas) ; il est mesuré pas à pas par le banc (o3) : **0,024 MPa au plus, 0,010 % de pc**
(tan ψ Δλ est petit à 1e-6 de déformation par pas), identité K (Δtr ε_pl + Δpc/H) vérifiée à 2,8e-7 Pa. Pas de
second retour volumique ajouté (il rouvrirait F : un retour de coin cap–cône exact serait une autre brique).
(ii) Le cap agit sur la pression **effective** p̄ = p_nom/(1−d) : dans une zone endommagée il se déclenche à une
pression nominale (1−d) pc — 35 MPa nominaux pour pc0 = 440 MPa à d_c = 0,92, et c'est ce qui fait frôler 256 MPa
effectifs au triaxial 20 MPa (σ₃/(1−d) = 95 MPa effectifs en fin de branche, pic à d_c 0,79). **`cdpCapP0` se calibre
donc en pression effective** ; sous l'insert (d_c élevé) le broyage démarre bien avant pc0 nominal.

Banc (n) (`selftest-cdp`, compression hydrostatique en déformation isotrope ε = −e I, e → 3 %, carte
historique, cap 440 MPa, H = K = 61,63 GPa) : 237 pas élastiques (p = K|ε_v| à 4e-14 %), 2763 pas cappés depuis
|ε_v| = 7,14e-3 (formule fermée à 7e-13 %, pc = pc0 + H ε_v^pl à 1e-12 %), pente dp/d|ε_v| = K/2 à 7e-13 %
(tol 0,5 %), fin p = pc = 2 993,6 MPa, ε_v^pl = 0,0414, wPlas = ∫p dε_v^pl à 0,027 %, d = ε_c^pl = ε_t^pl = 0,
ε_pl isotrope ; sans la clé p = K|ε_v| partout (4e-14 %), ε_pl = 0, pc = 0 ; `cdpCap = false` explicite
bit-identique ; **doit échouer** : `cdpCap = true` sans `cdpCapP0` → « cdp: cdpCap = true requires cdpCapP0 > 0 »,
`cdpCapP0` sans `cdpCap` refusé, `cdpCapP0 ≤ 0` refusé ; **revue** : `cdpCapH` absent → H = K exactement
(différence 0 Pa) et trace hydro bit-identique à `cdpCapH = K` explicite. Banc (o), **refondu à la revue** (la version
du soir, cap 440 seul, était vide : le cap ne s'activait jamais et l'inférence « n'agit que sur l'hydrostatique »
n'était pas fondée) — triaxial pilote du banc (a), 15 000 pas de 1e-6, sans/avec cap (H = K) :
**o1** σ₃ = 20 MPa, cap 440 : p̄ max **mesuré** 256,4 MPa (σ₃ effectif fin 95 MPa, d_c 0,79), pc reste à pc0 sur les
15 000 pas (jamais actif), q identique (écart 0) — invariance stricte d'une clé posée mais inactive ;
**o2** σ₃ = 20 MPa, cap **200** < 256 : activation au pas 12 361 (ε_ax −0,01236, **p_nom 56,9 MPa** mais d_c 0,72 →
p̄ = 200,05 MPa, post-pic : pic au pas 3 521), identité à 2,9e-11 avant, puis **doit diverger** : |Δq| max 0,775 MPa
(critère > 0,1), wPlas 991 664 → 1 026 990 J/m³, −tr ε_pl −0,00754 → −0,00663 (compaction supplémentaire
+9,1e-4), pc fin 249,8 MPa, d_c fin 0,790 → 0,783, q_pic 202,935 MPa des deux côtés (écart 0) ; **o3** σ₃ = 100 MPa,
cap **100** : le pilote part de σ = 0 (σ_lat tenue à −σ₃ dès le pas 1, σ_ax = λ tr ε ≈ −58 MPa), p̄ atteint pc0 au pas
540 exactement où q repasse par 0 (E dε = 77,7 kPa/pas) — en `matpoint` (confinement d'abord) ce serait le pas 0 ;
14 459 pas cappés dont **8 112 avec le cône actif** (Δε_c^pl > 0 : le couplage est exercé), q_pic 508,33 → 508,33 MPa
(+0,0007 %), ε_ax au pic −0,007887 → −0,008803, décalage −9,16e-4 = −ε_v^pl(cap)/3 à 0,03 % (pc au pic 269,4 MPa),
dépassement p̄ − pc max 0,0245 MPa (0,010 % de pc), identité K (Δtr ε_pl + Δpc/H) à 2,8e-7 Pa, fin pc 415,8 MPa,
−tr ε_pl −0,00283 → +0,00265 (la compaction l'emporte sur la dilatance), wPlas 1,19 → 2,17 MJ/m³, d_c fin 0,47 → 0,39. Deck `cdp_rockim/matpoint_cdp_hydro_cap.cfg` (`mpPath = hydro`) : première ligne cappée à
ε_v = −7,14e-3 (p = 440,04 MPa), fin p = pc = 2 993,6 MPa — mêmes chiffres que le banc.
Pilotage latéral (`Hold1D`, revue) : premier pas par la compliance **endommagée** (2(λ+μ) max(1−d, 0,02)),
crochet par recherche symétrique en doublement, regula falsi d'Illinois + bissection de garde, état
restauré à chaque essai ; 3,88 M appels, pire résidu 1e-3 Pa, 0 pas non convergé — et le verdict global
rend 1 si un pas ne converge pas (`selftest-triax` est intouché).

**Pilote générique `rockim matpoint <cfg> [out.csv]` (`matpointDrive`, toutes les lois).**

| clé (défaut) | rôle |
|---|---|
| `mpPath` (triax) | triax \| tension \| biaxial \| uniaxial (axe piloté z) ; triax/uniaxial : ε_zz décroissante, σ_xx = σ_yy = −σ₃ tenues ; tension : ε_zz croissante ; biaxial : ε_xx = ε_yy décroissantes, σ_zz = −σ₃ tenue ; **hydro** (`rockim_f2w12.exe`) : déformation isotrope imposée ε = −e I, e de 0 à `mpStrainMax` en `mpSteps` (aucun pilotage latéral, `mpSigma3`/`mpConfineFirst` ignorés), CSV propre `eps_iso, eps_v, p_nom, pc, eps_v_pl, d, wPlas, eroded` — les autres chemins sont intouchés |
| `mpSigma3` (« 0 ») | pressions [Pa] séparées par des espaces |
| `mpStrainMax` (0,03), `mpSteps` (3000), `mpLc` (1e-3 m), `mpDt` (1 s, lois visqueuses) | pas de déformation ε_max/n, lc du crack band / des tables CDP |
| `mpConfineFirst` (true) | consolidation isotrope à σ₃ (100 pas, p tenue par le pilote) puis phase axiale ; ε comptées depuis la fin de la consolidation |
| sortie | `sigma3, eps_ax, eps_lat, eps_vol, q, sig_ax, sig_lat, d_t, d_c, eps_c_pl, eps_t_pl, wPlas, wDamT, wDamC, residu, eps_in_c` (`eps_in_c` = ε_c^pl + d_c/(1−d_c) σ_c/E0, table lue à `mpLc`, ajoutée le soir du 2026-09-04 pour reconstruire u_in = (eps_in_c − ε_in,pic) × mpLc avec `cdpCompLength` ; autres lois : epvEq + Dc/(1−Dc)\|sig_ax\|/E0, proxy uniaxial) ; q = −sig_ax − σ₃ (compression positive ; tension : q = sig_ax − sig_lat) ; d_t = D, d_c = Dc, eps_c_pl = epvEq, eps_t_pl = kappa (cdp : exactement ses variables ; dpr : D, ω_c, ε_vp, κ de Rankine) ; `residu` = \|résidu final du pilotage latéral\| [Pa] de la ligne (> 1e-3 = pas non convergé) |
| code retour | 0, ou **2** si au moins un pas n'a pas convergé (`[matpoint] WARN`) — falsifiant pour une boucle de calibration |

Coût : ≈ 90 ms par confinement (4000 pas, ≈ 18 appels à la loi par pas). Decks : `cdp_rockim/matpoint_cdp_hist.cfg`,
`matpoint_cdp_inverse.cfg` ; cas de revue du pilote : `cdp_rockim/revue/{bigsteps_ten,tiny_lc,psi55_biax}.cfg`
(tous 0 pas non convergé). Un pas ≤ 1e-5 reste conseillé en traction post-pic pour la **précision** du chemin
(la convergence du pilote n'en dépend plus).

**Solveur `fem3d` : essai triaxial continu (`Fem3dSolver.cpp`, scénario `tension`, pullV < 0)**

| clé (défaut) | rôle |
|---|---|
| `pullDelay` (0 = bit-identique) | avant t = pullDelay les nœuds PRESCRIBED du mors supérieur sont traités comme **libres** (chargés par `topPressure`, amortis, le fond reste FIXED) ; à t = pullDelay ils redeviennent prescrits depuis leur position courante, rampe `pullRamp` comptée depuis pullDelay. Consolidation isotrope exacte (σ_xx = σ_yy = σ_zz = −P dans le tiers central) puis phase déviatoire. La pression de dessus agit sur ces nœuds tant qu'ils sont libres et est absorbée par le mors ensuite ; la colonne `sigma` (= \|F_grip\|/section) mesure alors le **déviateur** q, pas la contrainte totale (lire `sigZZmid`). Scénario tension seulement |
| `gripSection` (W·D) | section réelle pour σ = \|F_grip\|/section (éprouvette importée cylindrique) |
| `triaxStats` (false) | colonnes EN FIN DE LIGNE de `history.csv` (après `fieldStats`) : `sigZZmid, sigXXmid, sigYYmid` (moyennes pondérées par V₀ sur le tiers central `midEl_`), `epsAxMid, epsVolMid` (déformations de Biot moyennes), `epsAxGrip` = u_z moyen du mors / H. SI. Les ε comptent depuis t = 0 : soustraire la valeur à pullDelay pour la phase déviatoire |

Banc court (`cdp_rockim/bench_fem3d/`, élastique, grille 8×8×16, 20×20×40 mm, P = 50 MPa latéral + dessus,
pullDelay = 3 rampes, OMP 4, 7 s par run) : à la fin de la consolidation σ_xx / σ_yy / σ_zz = −49,97 /
−49,97 / −49,86 MPa (tol ±0,5), KE 3,5e-7 J, ε_vol = −P/K à 0,13 % ; phase axiale q/ε_ax = 77,64 GPa (E à
−0,03 %, 501 points) ; variante **qui doit échouer** (pullDelay = 0) : σ_zz = −31,3 MPa à la jauge de
confinement (écart 18,7 MPa > 10 ; −2νP = −29 MPa attendu à ε_zz = 0). `python check_bench.py out_delay out_nodelay`.

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
10. **Gardes des entrées (w20, 2026-09-05, chantier C du plan de robustesse : C3, C4, C6).** Trois familles
    d'erreurs *nommées*, levées avant le premier pas ou à cadence fixe, dans les six solveurs. Aucune ne change
    un flottant d'un deck valide : bit-identité w19 = w20 vérifiée sur quatre decks courts (fem3d cdp, fem3d dpr,
    fdem 2D voronoï, fdem3d grille — `etude_lois_fem/bitid_w20/comparaison_w20.txt`).
    - **C3 maillage** (`include/rockim/Guards.hpp`, appelé par `Fem3dSolver::buildMeshFile/finishMesh/buildMesh`,
      `Fdem3dSolver::buildMeshFile/buildFromTets`, `FdemSolver::buildMeshFile/buildFromTriangles`,
      `FemSolver::buildMesh`, `DemSolver/Dem3dSolver::init`) : (i) *nœud orphelin* — jamais référencé par un
      élément — `[rockim] error: mesh: N noeuds orphelins, premier : id <Gmsh>, (x, y, z) ; nettoyez le maillage
      (meshes/drop_orphans.py)`, contrôlé sur les coordonnées du fichier (avant translation) ; c'est la « broche
      fantôme » du 05/09 (point du champ de taille Gmsh, nœud 9 de tous les `T1_*.msh` non nettoyés et de
      `perc3d/Q1_c05.msh`, masse nulle épinglée FIXED sous le pôle d'impact, 5,9 kN sur 45,6). Plus **aucun**
      épinglage silencieux : seule exception, les nœuds de grille hors du cylindre en `geometry = cylinder`
      (fem3d), sans élément par construction, comptés et imprimés. (ii) *masse nodale nulle ou négative après
      lumping* (nœud, coordonnées, m). (iii) *élément dégénéré* : volume (aire) ≤ 0 après réparation
      d'orientation, ou < 1e-6 × médiane du maillage (sliver) — élément, ses nœuds et coordonnées, volume,
      médiane, seuil ; remplace « degenerate tet » sans identifiant. Le message « unknown node id » nomme
      désormais l'élément et le nœud.
    - **C4 NaN/Inf réel** : `checkFinite()` balaie **toutes** les composantes de u, v et f de **tous** les nœuds
      (x, v, f des particules en DEM) plus le travail et la force de l'outil, tous les `nanCheckEvery` pas
      (défaut **256**, 0 = off) et une dernière fois dans `finalize()`. Il remplace le détecteur *aveugle* de
      `Fem3dSolver` (`u_[0]`, nœud possiblement FIXED donc toujours fini) et l'échantillon E5 (~256 nœuds
      tous les 1024 pas) des fdem 2D/3D. Arrêt propre : `[rockim] error: <MODE> : NaN/Inf detecte au pas N
      (t = …) : noeud i (X0 = …), champ f composante x = nan, element voisin e`, écriture de
      `<outputDir>/ERROR.txt`, **code de retour 3** (les autres erreurs restent à 1). Coût mesuré : deck dpr court `C_T1_R_P000_court_grid` (55 296 tets, 4 209 nœuds, OMP 2, matrice v2 en parallèle), w20 : `nanCheckEvery = 256` 103,7 / 104,7 s contre `nanCheckEvery = 0` 104,7 / 102,3 s (+0,7 %, dans le bruit) ; `nanCheckEvery = 1` (balayage à CHAQUE pas) 107,6 s (+4 %) ; history.csv identique dans les trois cas (`bitid_w20/cost_*.log`).
      Limite connue : les énergies internes (joints, intégration) ne sont pas dans le balayage — un run de
      traction 2D à `dtFactor = 50` finit à code 0 avec `joints : nan J/m` alors que u, v, f sont finis
      (voir `bitid_w20/selftest_gardes/SELFTEST_gardes.md`).
    - **C6 clés par mode** : `hydro` et toute clé `hydro*` refusées hors `mode = fdem` (à côté de `thermal` et
      `bedding*` dans `main.cpp`) ; et **registre des clés par mode** `tools/keys_by_mode.json` → table
      compilée `include/rockim/KeysByMode.hpp`, tous deux GÉNÉRÉS par `tools/gen_keys_by_mode.py` (à relancer
      à chaque nouvelle clé, le build ne le fait pas). Méthode : lecture des `getd/getb/gets/geti/reqs/reqd/has`
      à clé littérale dans `src/` et `include/` ; propriétaire = fichier (`Fem3dSolver` → fem3d…, tout autre
      fichier = code partagé) ; clé lue par le code partagé ou par ≥ 2 solveurs = **commune**, jamais
      refusée ; lue par **un seul** solveur = clé de ce mode, refusée ailleurs : `cle 'X' sans effet en mode Y :
      cle du mode Z seulement`. Clés construites (`groupBond.<A>.<B>`, `groupPhase.<n>`, préfixe `cdp*`) hors
      registre, jamais refusées. Doute = commune (`ALWAYS_COMMON`). État au 05/09 : 388 clés, 238 communes,
      150 propres (fdem 112, fem3d 24 dont `quarterModel`, `kinematics`, `toolPulse*`, `bottomFree`, `probes`,
      `toolLockXY` ; fdem3d 7 dont `trackGroups`, `bulkViscosityGraded` ; fem 6 ; dem 1) ; essai à blanc
      `--scan-decks` sur 746 decks du dépôt : 0 refus.
    - **Banc falsifiant** `etude_lois_fem/bitid_w20/selftest_gardes/` (`run_gardes.py`, **30 / 30 PASS**) : le
      `T1_c05.msh` original est refusé (nœud 9), `T1_c05_clean.msh` passe ; tet plat et sliver à la main refusés
      (fem3d, fdem3d), triangle plat refusé (fdem) ; NaN → code 3 + ERROR.txt dans les six solveurs ; `hydro`
      en fem3d et fdem3d refusé ; clé fdem3d en fem3d, clé fem3d en fdem3d et en fdem refusées ; neuf témoins
      valides passent. Suite `verify_suite.py --tier fast` (OMP 1) : **TOUT PASSE (48/48)**, OMP 1, ~40 min sous charge (matrice v2 + bitid en parallele ; 22 min a vide), `bitid_w20/suite_fast_w20.log` + `.json` — aucun faux positif des gardes C3/C4/C6 sur les 48 tests.
    - Conséquence pour les références : le deck de bit-identité `bitid_w12/PQ_cdpI_P020_court.cfg` pointait
      `perc3d/Q1_c05.msh` (orphelin) — il est désormais **refusé** ; la chaîne cdp est rebasée sur
      `bitid_w20/PQ_cdpI_P020_court_clean.cfg` (`Q1_c05_clean.msh`) : w19 = w20, 6/6 IDENTIQUE, history e45a7cbc0609c9f8 / frames 5ddf25b84b9c15f1 / vtu a1e108561139734a, pic 7 791,04 N — history et frames sont ceux de la chaîne w14 → w19 sur le maillage à orphelin (la garde `m ≤ 0` du contact, w13, excluait déjà ce nœud), seul le hachage VTU change (un nœud de moins).
11. **La clé inconnue est une ERREUR (w21, 2026-09-05 — C1 et C7 du plan de robustesse ; décision de Fernando, 20:00 ;
    w22 le même soir : corrections D1 et D2 de la relecture adverse).**
    Jusqu'à w20, `Config` était une table sans suivi : une clé mal orthographiée, obsolète ou d'un autre mode était
    ignorée en silence (`fragBrushV` inerte dans 10 decks pendant des semaines, `kpFactor` en fem3d, `hydro` en 3D…).
    - **Mécanisme** (`include/rockim/Config.hpp`, `src/Config.cpp`) : chaque getter — `gets/getd/geti/getb/reqs/reqd/has`
      et `keysWithPrefix` (les clés rendues) — **marque la clé consommée** et note le défaut employé quand la clé est
      absente. Les solveurs copient `Config` par valeur (`Config cfg_`) : le stockage est **partagé** entre copies
      (`shared_ptr`), les lectures faites par le solveur, `Material::from`, `MatLaw`… sont donc vues par le deck que
      `main.cpp` audite. `keys()` (registre) ne marque rien. Lecture pure : aucun flottant touché, bit-identité w18 =
      w21 = w22 sur les 8 decks de `tools/bitid.py` (voir `etude_lois_fem/bitid_w21/`, `bitid_w22/`). **w22** : après
      l'audit et `config_effective.cfg`, `main.cpp` appelle `cfg.seal()` — les getters ne suivent plus rien (ni verrou ni
      chaîne du défaut) et lisent la table directement, comme avant w21 (`confineGaugeTime` est lu à **chaque pas** dans
      `FdemSolver::step` / `Fdem3dSolver::step` : le coût du suivi ne s'applique qu'à l'initialisation).
    - **Où** : `main.cpp`, **après `solver->init()`** (toutes les lectures conditionnelles de l'initialisation ont eu
      lieu) — `keyguard::enforce()` puis `keyguard::writeEffective()` (`include/rockim/KeyGuard.hpp`).
    - **Règle**, pour chaque clé du deck jamais consommée, dans cet ordre : (1) nom dans la table des **obsolètes**
      (`tools/obsolete_keys.json` → `kObsolete`) → `cle obsolete 'X' (ligne N du deck) : remplacee par 'Z'` ;
      (2) clé **du registre** : on consulte ses **lecteurs** (`kReaders`, w22 : les solveurs dont le fichier lit la clé,
      ou `shared` = code partagé — `main`, `Material`, `MatLaw`, `Tool*`… — et `ALWAYS_COMMON`). **Si le mode courant
      ou `shared` est lecteur → légitime, rien** — c'est la **lecture conditionnelle** : `capP0` n'est lu que si
      `dprCap = true`, `hydroStart` que si `hydro = on`, `jointResidualMu` que sous certains adoucissements, `gravity`
      qu'après init (`placeTool`) ; une clé que le solveur du mode lit dans un chemin du code existe et agit, la
      refuser serait un faux refus (le « réglage inerte » d'une telle clé reste possible sans message : prix assumé de
      zéro faux refus). **Sinon → `cle 'X' sans effet en mode Y : cle du mode Z seulement`** (un lecteur ; même
      sous-chaîne que le contrôle w20, désormais **fondu** ici — `keysbymode::check()` reste dans le header mais n'est
      plus appelé) **ou `… : cles des modes Z1, Z2`** (plusieurs lecteurs, aucun du mode). *w21 traitait toute clé à
      ≥ 2 lecteurs comme « commune = légitime » : `jointSoftening`, `insertion`, `bulkDamage`, `fragBrushV0`, `gc*`,
      `yan*`, `strainRate*`, `gravity` (fdem + fdem3d) passaient en fem3d sans un mot, `bond*`, `packing`,
      `particleRadius` hors dem, `kpFactor` hors fem/fem3d — trou D1 de la relecture adverse, fermé en w22 ; 0 deck
      existant sur 957 n'est touché (`tools/scan_decks.py`).* (3) clé d'une **famille dynamique** (`kDynamicPrefix` :
      `phase.<nom>.E`, `gb.<a>.<b>.<prop>`, `contactMu.<phase>`, `groupBond.<A>.<B>`, `groupPhase.<g>`,
      `groupVel.<g>`, `gauge.<g>`) non lue → `cle 'X' jamais lue : famille dynamique 'phase.<nom>...' — nom de phase /
      groupe declare nulle part ?` ; (4) sinon, **suggestion** : distance de Levenshtein ≤ 2 (ou même nom à la casse
      près) sur l'union registre + clés consommées, au plus 3, en précisant les modes où la suggestion agit si ce n'est
      pas le mode courant → `cle inconnue 'X' (ligne N du deck) : vouliez-vous dire 'Y' ?` ; (5) sinon `cle 'X' (ligne
      N du deck) inconnue de rockim`.
    - **Toutes** les clés fautives sont listées d'un coup (ordre du deck), puis : **`unknownKeys = error` (défaut)** →
      `[rockim] error: N cles du deck '…' refusees (unknownKeys = error) ; poser unknownKeys = warn pour continuer avec
      un avertissement :` + une ligne `- …` par clé, **code de retour 1** ; **`unknownKeys = warn`** → mêmes lignes en
      `[rockim] WARNING:` sur stderr, le run continue (vieux decks, campagnes en cours). Toute autre valeur est refusée.
      La garde `hydro*` hors fdem (C6) reste une erreur immédiate avant init ; en fdem, `keysWithPrefix("hydro")` n'est
      plus appelé pour qu'une clé `hydro*` mal orthographiée reste visible par l'audit.
    - **Ce qu'il faut regénérer** : `python tools/gen_keys_by_mode.py` à chaque nouvelle clé ou renommage (le build ne
      le fait pas) — il produit `tools/keys_by_mode.json` (registre : `modes`, `common`, `prefixes_common`,
      `prefixes_dynamic`, `obsolete`, `readers`) et `include/rockim/KeysByMode.hpp` (`kTable`, `kKnown`,
      `kDynamicPrefix`, `kCommonPrefix`, `kObsolete`, `kReaders` w22). **Un renommage de clé = une entrée dans
      `tools/obsolete_keys.json` dans le même commit** (ne jamais en retirer : les vieux decks doivent rester
      diagnosticables). État au 05/09 : 389 clés connues, 239 communes, 150 propres à un mode, 7 familles dynamiques,
      1 obsolète (`fragBrushV` → `fragBrushV0`) ; lecteurs : 85 `shared`, 112 fdem seul, 85 fdem + fdem3d, 24 fem3d
      seul, 18 les six solveurs, 11 dem + dem3d…
    - **C7 — `<outputDir>/config_effective.cfg`** (format w22) : écrit à la fin de l'initialisation, **tri stable par
      clé** : lignes **actives** = toutes les clés du deck sauf les fautives — `cle = valeur   # deck ligne n`
      (consommée) ou `cle = valeur   # deck ligne n ; non lue a l'initialisation (lecture conditionnelle ou apres init)`
      (légitime mais non consommée, ex. `gravity`, `hydroStart` sans `hydro`) — et défauts du code consommés en lignes
      **commentées** `# cle = valeur   (defaut)` (`; autres defauts lus : …` si deux lectures à défauts différents, ex.
      `absorbSpringR`) ; en-tête : deck, mode, exe, politique, comptes, et la liste des clés refusées/averties (retirées
      des lignes actives). **Le fichier est un deck REJOUABLE** : ses clés actives sont celles du deck d'origine moins
      les refusées → mêmes résultats (banc `selftest_cles` partie C : `history.csv` bit-identique au rejeu sur les six
      modes). *w21 écrivait les défauts en clés actives : au rejeu, les gardes « satellite orpheline » (`has()`) les
      voyaient posées — fem3d rc 1 `tensionShearRetention requires tensionDamage = fixed`, fdem rc 1 `dampingLocalAfter
      est posee sans dampingSwitchT` — défaut D2 de la relecture, corrigé en w22.* La console imprime
      `[rockim] cles : N consommees (n du deck, m au defaut), k du deck non lues [dont f fautives] -> …/config_effective.cfg`
      (actives = n + k − f, commentées = m) ; une clé lue plus tard dans `step()` (ex. `gravity` dans `placeTool`)
      compte parmi les « non lues » et reste active sans être marquée consommée : limite connue, E16 du plan.
    - **Balayage statique des decks** : `python tools/scan_decks.py [--fix-obsolete]` applique la règle sans lancer
      rockim (registre + lecteurs + obsolètes ; les clés dynamiques et les lectures conditionnelles ne sont pas
      jugeables) sur `rockim_f2` + `CONTINUUM/calib_bohus_triax/cdp_rockim` → `tools/scan_decks_<date>.md`. w22 : il
      modélise aussi les **gardes pré-init** de `main.cpp` (`thermal`/`beddingDip` hors fdem, `law` hors
      fem3d/fdem/fdem3d, `phases`/`mesh = voronoi` hors fdem/fdem3d, `mesh = file` hors fdem/fdem3d/fem3d, `hydro*`
      hors fdem, mode inconnu) et met à part les fichiers **sans clé de solveur** (decks `matpoint` à clés `mp*`, cartes
      matériau `law = cdp` + constantes, decks temporaires de selftest : ils ne passent jamais par l'audit) ; il ne voit
      **pas** les gardes de maillage C3 (nœud orphelin, tet plat : il faudrait lire le `.msh`, ex. `tunnel_schisto/
      S4_ratios1.cfg`). Le 05/09 : 957 fichiers, 881 decks de solveur, 7 `fragBrushV` corrigés en place (C0 ; aussi les 7
      de `rockim_f2_wt`), 8 clés inconnues (`dfhPsiVar`, `dfhPsi0`, `dfhKPsi`, `dfhPsiMax` dans
      `bench_impact/configs/impact3d_dpdfh*.cfg` : la dilatance variable ψ(p) de `vumat_hole.f` n'existe **pas** dans le
      `dpdfh` de rockim, qui lit `dfhPsiDeg` — à trancher par Fernando), 1 garde pré-init (`tests_f2/tg_3d.cfg` :
      `thermal` en fdem3d), **0 clé « sans effet en mode » par la règle des lecteurs**, 0 faute de frappe, 400 clés
      dynamiques, 76 fichiers sans clé de solveur.
    - **Banc falsifiant** `etude_lois_fem/bitid_w21/selftest_cles/` (`run_cles.py`, `SELFTEST_cles.md`) : faute de
      frappe (`contactMuu` → suggestion `contactMu`), obsolète (`fragBrushV` → `fragBrushV0`), inventée
      (`zorglubFactor` → inconnue de rockim), les trois d'un coup (3 lignes, rc 1), `unknownKeys = warn` (WARNING ×3,
      rc 0, `config_effective.cfg` écrit), clé d'un autre mode (`trackGroups` en fem3d, régression w20), clé
      conditionnelle légitime (`hydroStart` en fdem sans hydro, `capP0` sans dprCap : rc 0), famille dynamique
      (`phase.granit.E` sans phase granit : rc 1), politique invalide (`unknownKeys = maybe` : rc 1), casse
      (`ContactMu`), six témoins valides ; **w22** : clés à plusieurs lecteurs hors mode (`jointSoftening`, `insertion`,
      `bulkDamage`, `fragBrushV0` en fem3d → 4 refus « cles des modes fdem, fdem3d » ; `packing`, `particleRadius` en
      fem3d → « dem, dem3d » ; `kpFactor` en fdem → « fem, fem3d » ; les mêmes 4 avec `unknownKeys = warn` → rc 0 et
      4 WARNING) ; **partie C** : rejeu des six témoins depuis leur `config_effective.cfg` → rc 0 et `history.csv`
      bit-identique ; partie B : les 20 tests les plus courts du tier fast (rc 0 + `config_effective.cfg` cohérent :
      actives = n + k − f, commentées = m).

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

### 5.20 Percussion fem3d : contact de Signorini, impulsion de force, quart de bloc, garde des nœuds sans masse (2026-09-04/05)

| clé (défaut) | effet |
|---|---|
| `toolContact = penalty \| signorini` (penalty) | fem3d : contact outil en IMPULSION (port de la branche A1 de fdem3d, noyau `ToolSignorini.hpp`) pour la sphère et le poinçon plat ; par nœud, v* = v + (dt/m) f, condition de Signorini sur le gap prédit, r_n = m (relax·pen/dt − v_n), cap de Coulomb sur l'impulsion tangentielle, report en force r/dt ; `toolSignoriniRelax` (0). La lame garde la pénalité. Banc de Hertz (quart de bloc élastique, 16 J) : Signorini = pénalité ×10 = ×100 à 1 % ; la pénalité ×1 (défaut, kp = E h_min avec h_min = sliver du Delaunay) donnait −18 % de force et +12 % d'enfoncement |
| `contactPenaltyFactor` (1) | kp = facteur × E h_min ; pulsation de pénalité et rapport à 2/dt imprimés, avertissement au-delà de 0,5 |
| `quarterModel` (false) | plans x = 0 et y = 0 symétriques (v_x = 0 sur x = 0, v_y = 0 sur y = 0 seulement), ni pression suiveuse ni Lysmer sur ces faces ; insert sur l'arête (0, 0), masse d'outil à diviser par 4 dans le deck |
| `toolPulseForce` (0), `toolPulseTable`, `toolPulseImpedance` (0) | impulsion de force imposée sur la masse de l'outil, F_z = −toolPulseForce × a(t), a(t) linéaire par morceaux dans la table "t:a t:a …" (0 hors table) — le protocole tige + `*Cload amplitude` des decks Abaqus (P_cdpQ_v11 : 21 318,4 N de crête sur le quart, 185 µs). Appliquée à l'intégration seulement : `toolFz`, `peakF_`, `work_` restent la force de contact ; résumé : ∫F dt et ∫F v dt. Avec `impactSpeed = 0` et `toolGap` petit  Avec `toolPulseImpedance = Z` (N s/m, ρcA de la tige) : source à impédance, F = 2 F_inc − Z v_bouton, `toolPulseForce` = amplitude de l'onde incidente (deck à *Cload au sommet d'une tige absorbée par des dashpots Z_d : F_inc = Cload × Z/(Z + Z_d)). Une force pure (Z = 0) pousse l'insert à travers la roche (vérifié : 8 mm, 13 m/s)|
| garde nœuds sans masse | `toolContact()` ignore les nœuds de masse nulle (nœud 0-D d'un maillage Gmsh jamais référencé par un tétraèdre, épinglé FIXED à la lecture). Avant : la pénalité lui appliquait kp × pénétration = BROCHE FANTÔME sous le pôle quand le point du champ de taille Gmsh coïncide avec le point d'impact (tous les maillages T1 de l'étude des briques : 5,9 kN à P = 0 sur 45,6, 1,9 J de ressort) ; Signorini divisait par sa masse. Outils `etude_lois_fem/meshes/check_orphans.py`, `drop_orphans.py` (→ `*_clean.msh`). Change les résultats des seuls maillages à nœud orphelin |

#### Viscosité de volume `bulkViscosity = b1 b2` (2026-09-05, `rockim_f2w15.exe`, opt-in)

| clé (défaut) | effet |
|---|---|
| `bulkViscosity = b1 b2` (absent) | fem3d : viscosité de volume à la Abaqus/Explicit (*Analysis User's Guide*, « Bulk viscosity » ; défauts Abaqus 0,06 et 1,2). Clé absente **ou** `0 0` : rien n'est ajouté, aucune colonne, dt inchangé (bit-identique, prouvé `etude_lois_fem/bitid_w15/`). Active (b1 > 0 ou b2 > 0) : contrainte visqueuse isotrope ajoutée aux forces internes, colonne `wBulk` en **dernière** position de `history.csv`, ligne de résumé, facteur sur le dt |

**Formules (par élément, à chaque pas).** Taux de déformation volumique ε̇_vol = (J_{n+1} − J_n)/(dt J_{n+1}), J = det F déjà calculé dans `elementForces` (mémoire `Elem::Jbv` du J précédent, lue et mise à jour seulement quand la clé est active). Longueur d'élément L_e = `lc` = V0^{1/3} (celle de la crack band ; Abaqus prend sa « characteristic element length », pour un C3D4 une longueur du même ordre — on garde V0^{1/3} pour n'avoir qu'une seule longueur par élément dans le code). Célérité dilatationnelle NON endommagée c_d = √((λ + 2µ)/ρ) = `Material::cP()`.

- pression linéaire p₁ = b1 ρ c_d L_e ε̇_vol (agit dans les deux sens : amortit la sonnerie derrière un front, en compression comme en traction) ;
- pression quadratique p₂ = ρ (b2 L_e ε̇_vol)² **seulement en compression volumique** (ε̇_vol < 0 ; nulle sinon : le choc, pas la détente) ;
- contrainte visqueuse σ_bv = q I avec q = p₁ − [ε̇_vol < 0] p₂, donc q du signe de ε̇_vol : en compression rapide elle **ajoute de la compression**, en expansion rapide de la traction — la pression visqueuse s'oppose au taux volumique. Elle est ajoutée à la contrainte de la loi **pour les forces internes seulement** : `svm`, `pm`, `szz`, `sigLat`, `wEl`, les jauges triaxiales et les champs VTU restent la contrainte matériau (comme Abaqus, où la bulk viscosity n'apparaît pas dans S).
- dissipation `wBulk` = Σ_el V0 q ε̇_vol dt ≥ 0 par construction (q ε̇_vol = b1 ρ c_d L_e ε̇² + [ε̇<0] ρ b2² L_e² |ε̇|³) ; cumulée, imprimée dans le résumé (`dissipation wBulk = … J`) et écrite en dernière colonne de `history.csv` ; réduction OpenMP dans l'ordre des threads (déterministe à nombre de threads fixé, comme les autres compteurs).
- pas de temps : Abaqus réduit le dt stable de l'élément du facteur √(1 + ξ²) − ξ, ξ = b1 − b2² L_e ε̇_vol/c_d ; `computeStableDt` applique la **part linéaire** ξ = b1 à la CFL (0,9418 pour b1 = 0,06 ; 0,7440 pour 0,3) — la part quadratique dépend du taux courant (pas de dt adaptatif dans rockim) et n'est **pas** portée ; la borne du ressort de contact 2√(m_min/kp) n'est pas touchée (elle peut rester la borne active : en scénario tension kp = E h_min borne aussi dt, comportement hérité — poser `contactPenaltyFactor` < 1 si l'on veut lire le facteur sur dt).

**Ce qui n'est pas porté.** Le terme quadratique de ξ dans le dt ; la longueur caractéristique exacte d'Abaqus (V/A_max pour les tétraèdres) ; la viscosité de volume des éléments 2D (`fem`), des cohésifs et de `fdem3d` — clé fem3d seulement ; l'énergie `ALLVD` d'Abaqus n'a d'équivalent que `wBulk` (la viscosité de contact et l'amortissement de Cundall restent comptés ailleurs).

**Banc falsifiant** (`etude_lois_fem/bitid_w15/selftest_bv/`, `SELFTEST_bv.md`) : bloc élastique 6 × 6 × 24 mm, échelon de vitesse −2 m/s (tension, `pullV < 0`, `pullRamp = 0`, damping 0) ; sonnerie = écart-type détendancé de la jauge de tête / du tiers central derrière le front. `0 0` bit-identique à sans clé ; `0.06 1.2` : sonnerie de tête −23 %, wBulk 8,1·10⁻⁵ J > 0, dt × 0,94180 ; `0.3 1.2` : −62 % (tête) / −46 % (tiers central), wBulk 2,2·10⁻⁴ J, dt × 0,74403 ; le témoin sans clé au même dt réduit ne baisse que de 6 % (tête) : c'est la viscosité qui amortit, pas le dt. Terme quadratique seul en rampe lente : wBulk = 0 en traction (1,8·10⁻³⁸ J) et > 0 en compression (porte compression vérifiée) ; sans chargement wBulk = 9·10⁻²⁴ J. Aucun scénario fem3d n'impose un mouvement à volume constant (cisaillement pur, rotation rigide) : le cas (v) n'est couvert que par le cas sans chargement et par la construction (q ne dépend que de det F, invariant par rotation, = 1 en cisaillement pur).

#### Cinématique `kinematics = hencky` (2026-09-05, `rockim_f2w16.exe`, opt-in)

| clé (défaut) | effet |
|---|---|
| `kinematics = biot \| hencky` (biot) | fem3d : mesure de déformation et couple conjugué. Clé absente **ou** `biot` : chemin historique, bit-identique (prouvé `etude_lois_fem/bitid_w16/comparaison_w16.txt`). `hencky` : déformation logarithmique ε = ln U, contrainte de la loi lue comme Cauchy co-rotationnelle, forces internes sur la configuration courante — le couple conjugué d'Abaqus/Explicit `nlgeom`. Toute autre valeur est refusée. Aucune colonne ni sortie nouvelle |

**Pourquoi.** Constat 1 de `cdp_rockim/audit_crush/SYNTHESE_audit_broyage.md` : sur le même maillage et la même carte CDP, rockim porte 20 % de moins qu'Abaqus au pic d'un impact d'insert ; la couche de contact est à J = 0,84-0,94 (20-30 % de déformation), là où la mesure de Biot et la mesure logarithmique diffèrent de 5-19 % sur la force transmise.

**Chemin historique (`biot`).** Par tétraèdre F = I + ∂u/∂X, rotation R par trois itérations de Newton (R ← ½(R + R⁻ᵀ) depuis F√3/‖F‖ ; exacte à ~10⁻⁷ près à 20 % de déformation), déformation de **Biot** ε = sym(Rᵀ F) − I, la loi rend σ_c dans le repère co-rotationnel, forces nodales f_a = −V0 (R σ_c) ∂N_a/∂X (premier Piola-Kirchhoff P = R σ_c, volume et gradients de **référence**) : le couple conjugué est (contrainte de Biot, déformation de Biot).

**Chemin `hencky`** (`include/rockim/Fem3dKinematics.hpp`, appelé dans `elementForces`) — décomposition **spectrale** de C = Fᵀ F = Σ λ_i² n_i n_iᵀ (`Eigen::SelfAdjointEigenSolver` 3×3, exacte au conditionnement près) :

- U = Σ λ_i n_i n_iᵀ, U⁻¹ = Σ (1/λ_i) n_i n_iᵀ, R = F U⁻¹ (polaire **exacte**, plus d'itération) ;
- ε = ln U = Σ ln(λ_i) n_i n_iᵀ, passée à la loi à la place de la déformation de Biot (interface `MatLaw::stress` inchangée : la loi ne sait pas quelle mesure elle reçoit ; `epsP` des lois plastiques devient une déformation plastique logarithmique) ;
- σ_c rendue par la loi = contrainte de **Cauchy** dans le repère co-rotationnel, σ = R σ_c Rᵀ ;
- forces internes sur la configuration **courante** : P = J σ F⁻ᵀ = J R σ_c U⁻¹ (F⁻ᵀ = R U⁻¹), f_a = −V0 P ∂N_a/∂X — équivalent à −V (R σ_c Rᵀ) ∂N_a/∂x avec V = J V0 et ∂N/∂x = F⁻ᵀ ∂N/∂X, mais avec les gradients de référence déjà stockés : une décomposition et trois produits 3×3 par élément ;
- garde-fou : det F ≤ 10⁻⁹ ou λ_min ≤ 10⁻⁹ (élément dégénéré, ou échec du solveur propre) → l'élément reprend le chemin historique (R = I, ε de Biot), exactement comme aujourd'hui ;
- viscosité de volume (`bulkViscosity`) : q I est ajoutée à σ_c **avant** P, dans les deux cinématiques (en hencky c'est une pression de Cauchy, comme dans Abaqus) ; OpenMP : tout est local à l'élément, aucune course.

**Ce qui est approché.** Abaqus/Explicit intègre le **taux** de déformation (D dt) dans le repère tourné (Green-Naghdi) — une déformation logarithmique *incrémentale*, égale à ln U seulement pour un trajet **coaxial** (directions principales de U fixes dans le repère matériel : traction/compression uniaxiale, sphère sous l'insert au pôle). rockim calcule ln U **total** à chaque pas : exact pour un trajet coaxial, différent d'Abaqus dès que les directions principales tournent dans la matière (cisaillement fini) ; le contrôle hyperélastique (ln U = mesure totale, pas de dérive d'intégration) est la contrepartie. Les deux cinématiques coïncident au premier ordre : à λ = 1,002 elles diffèrent de 0,2 %.

**Ce qui n'est pas porté / inchangé.** Masse lumpée sur V0 ; Lysmer et ressorts d'absorption (géométrie initiale) ; contact outil (géométrie courante des nœuds, déjà) ; `lc` = V0^{1/3} de la crack band (la longueur de référence, pas la longueur courante) ; **pas de temps** : la CFL est calculée une fois sur la géométrie initiale (`computeStableDt`) et n'est **pas** recalculée — sous forte compression la célérité apparente sur la configuration courante augmente, la marge est celle de `dtFactor` ; soupapes d'érosion : `erodeDetMin` lit det F (inchangé), `erodeStrainMax` lit ‖ε‖ de Frobenius — en hencky c'est ‖ln U‖ (0,3 ⇔ λ ≈ 0,74 en uniaxial au lieu de 0,70 en Biot) ; sorties `svm`, `pm`, `szz`, `sigLat`, champs VTU et jauges triaxiales (`sigZZmid`, `sigXXmid`, `sigYYmid`) = contrainte de la loi, donc **Cauchy** en hencky (Biot sinon), `epsAxMid`/`epsVolMid` = ln U (ε_vol = ln J exact) ; `wEl` = ½ σ_c : (ε − ε_p) par unité de volume de référence (nominal, pas l'énergie exacte) ; la force de mors et `toolFz` restent des forces vraies. Modes `fem` (2D), `fdem3d` : non concernés.

**Banc à formes fermées** (`etude_lois_fem/bitid_w16/selftest_hencky/`, `SELFTEST_hencky.md`) : bloc élastique 4 × 4 × 8 mm (768 tets), E 77,66 GPa, ν 0,29, `scenario = tension`, `gripLateralFree = true` (mors libres latéralement, faces libres : état uniaxial homogène à contrainte latérale nulle), rampe cosinus 200 µs puis ±1,8 m/s, Cundall 0,7, quasi statique (ρ c v = 29 MPa contre 15,5 GPa). Force de mors / A0 attendue : biot = E (λ − 1) ; hencky = λ^(−2ν) E ln λ (Cauchy E ln λ sur l'aire courante A0 λ^(−2ν)) ; à λ = 0,8 les deux diffèrent de 27 %, à 1,2 de 18 %, à 1,002 de 0,2 %. Mesuré (`selftest_hencky/resultats.md`) : chaque cinématique à **−0,01 %** de SA forme fermée (RMS 0,02 % sur le trajet), la mauvaise forme rejetée (hencky contre la forme biot +27,0 % à 0,8, −18,0 % à 1,2) ; `sigZZmid` (contrainte de la loi) = E ln λ à −0,01 % et σ_mors/|sigZZmid| = λ^(−2ν) (1,1382 à 0,8) : Cauchy sur l'aire courante, vérifié par deux jauges indépendantes ; petites déformations : hencky/biot −0,215 % ; clé absente ≡ `biot` octet pour octet, w15 ≡ w16 sans clé. Rotation rigide et objectivité : aucun scénario fem3d ne les impose — test unitaire `selftest_hencky/test_kinematics.cpp` (23/23) sur `Fem3dKinematics.hpp` : F = Q → ln U = 0, R = Q, P = 0 à 1e-16 ; F' = Q F → ε inchangée, P' = Q P. Coût : ×1,7 sur une loi élastique (la décomposition spectrale domine), moindre avec cdp. Deck percussion `cdp_rockim/perc3d/PQ2_cdpH_pulseZ2_abqmesh_hencky.cfg` (écrit à 06:28 ; lancé ensuite dans la file de nuit avec `bulkViscosity = 0.06 1.2` ajoutée à 07:05 — verdict ci-dessous).

#### Verdict du 2026-09-05 : les deux clés reproduisent la courbe d'Abaqus/Explicit

Impact d'insert R 7,94 mm, quart de bloc, carte CDP historique, tige à impédance, P = 0, sur le maillage Abaqus converti à l'identique (`cdp_rockim/perc3d/P_cdpQ_v11_block.msh`, 160 156 C3D4 ; référence `percussion_bohus/P_cdpQ_v11`, 42,5 kN à 161 µs ; conventions `perc3d/nuit_lib.py` : origine au premier contact F × 4 > 0,2 kN, pic sur la force filtrée 10 µs, forces × 4 ; bilan `perc3d/bilan_nuit.md`). Pic de force de l'insert complet : chemin historique **34,2 kN** (−20 %) → `bulkViscosity = 0.06 1.2` seule **38,3 kN** (bulbe V(d_c > 0,5) gelé au pic comme Abaqus, 149 → 168 mm³ contre 173 → 189 ; wBulk 0,009 J : c'est la sonnerie des éléments broyés qui est amortie, pas une dissipation) → `bulkViscosity = 0.06 1.2` + `kinematics = hencky` **42,5 kN** à 148 µs (Abaqus 42,5 à 161), enfoncement max 0,829 mm (0,838), résiduel 0,493 (0,516), impulsion 8,84 N s (8,97), bulbe 154 → 177 mm³. Force à enfoncement 0,5/0,6/0,7/0,75/0,8 mm : Abaqus 22,0/29,3/35,1/36,6/39,8 kN, historique 19,7/26,5/30,4/32,6/34,2, bv 21,0/27,7/33,1/35,3/37,3, bv + hencky 23,2/30,0/36,5/39,2/41,6. Les variantes suppression (`noero` 34,2), contact (`pen10` 32,8), vitesse d'arrivée (`gap1um` 34,5) n'expliquaient rien ; la loi est à parité au point matériel (`cdp_rockim/audit_crush/SYNTHESE_audit_broyage.md`). L'écart de 20 % était pour moitié la viscosité de volume et pour moitié la cinématique (Biot + forces sur V0 contre logarithmique + Cauchy courante, couche de contact à J 0,84-0,94) : deux ingrédients d'Abaqus/Explicit, aucune propriété du matériau.

**Règle.** Toute comparaison fem3d avec Abaqus/Explicit (même maillage, même carte) se fait avec **les deux clés** `bulkViscosity = 0.06 1.2` (défauts Abaqus) **et** `kinematics = hencky` (le couple conjugué de `nlgeom`), binaire `rockim_f2w16.exe` ou plus récent ; sans elles l'écart de pic attendu en régime broyé est de l'ordre de −20 %. Les défauts restent intacts (clés absentes = chemin historique, bit-identique). Sensibilités sur base bv (`bv_mu05` μ 0,5, `bv_Gf200` G_f × 2) : en cours au 05/09 matin, à lire dans `bilan_nuit.md`.

#### Outil bloqué latéralement `toolLockXY = true` (2026-09-05, `rockim_f2w17.exe`, opt-in)

| clé (défaut) | effet |
|---|---|
| `toolLockXY` (false) | fem3d, **percussion** : l'outil rigide garde v_x = v_y = 0 et x = `toolX`, y = `toolY` **exactement** à tous les pas (`Tool3::integrate` : v.x = v.y = 0 après la mise à jour par F, x += dt·0). La force de contact latérale est toujours calculée et comptée (`toolFx` de `history.csv`, `peakF_`, ligne de résumé) mais ne déplace pas l'outil — l'équivalent de `*Boundary RPB1 1,2` sur le nœud de référence du bouton des decks Abaqus. Elle ne travaille pas (`work_` = −F·v dt avec v_x = v_y = 0, comme un RP bloqué). Scénario **shear** (lame, sphère traînée : outil piloté en déplacement, v_x = cutSpeed) : clé **ignorée** avec `[FEM3D] WARNING`. Clé absente : chemin bit-identique (prouvé `etude_lois_fem/bitid_w17/comparaison_w17.txt`, w14 → w17) |
| résumé « tool lateral » | ligne de console ajoutée (aucun fichier) : x, y de fin, départ, max \|dx\|, \|dy\|, \|Fx\|, \|Fy\| sur tous les pas — l'observable de la dérive, qui n'était visible que dans `toolX` de l'historique et `toolY` des images |

**Pourquoi.** Contre-expertise `cdp_rockim/perc3d/DECHARGE_abaqus_vs_rockim.md` § 6 : dans le quart de bloc l'insert posé sur l'arête (0, 0) **dérive** vers les plans de symétrie (toolX −0,055 mm au pic, −0,18 mm à 320 µs, −0,81 mm à 800 µs) — le quart de sphère ne voit qu'un quart de la pression de contact, rien n'équilibre F_x, F_y (≈ 0,17 F_z au pic sur le deck court), et la vitesse x, y de l'outil est libre ; Abaqus bloque le RP. Banc `etude_lois_fem/bitid_w17/selftest_lock/` (`SELFTEST_lock.md`, PASS 5/5) : sans clé x ≠ 0 dès le premier contact (−1,15 µm à 57 µs de contact, croissance ~t³) ; avec la clé x = y = `0` sur 2102/2102 pas et 15/15 images, \|Fx\| jusqu'à 1,48 kN toujours compté ; contrôle shear : avertissement et l'outil avance de cutSpeed·t exactement. **Constat à relire** : la clé n'est pas neutre (pic filtré 10 µs −14 %, travail −6 % sur le deck court) alors que la dérive n'y est que de 1 µm — c'est la *vitesse* latérale de l'outil libre (10⁻² m/s ≫ `contactVreg` = 10⁻³ m/s) qui fixe la direction du frottement des nœuds du pôle ; avec la clé elle est nulle, comme pour le RP d'Abaqus. Deck de comparaison à l'identique, NON lancé : `perc3d/PQ2_cdpH_pulseZ2_abqmesh_bvh_nrb.cfg` (= bv + hencky + `absorbing = all` + `toolLockXY = true`, `bottomFree` retiré).

#### Sondes de point matériel `probes = x,y,z ; x,y,z ; …` (2026-09-05, `rockim_f2w18.exe`, opt-in)

| clé (défaut) | effet |
|---|---|
| `probes = x,y,z ; x,y,z ; …` (absente) | fem3d : suit des POINTS MATÉRIELS du bloc. Mètres, repère du bloc recalé [0,W]×[0,D]×[0,H] (dessus z = H ; `mesh = file` est translaté à l'origine à la lecture). Clé absente : **aucun test par élément, bit-identique** (prouvé `etude_lois_fem/bitid_w18/comparaison_w18.txt`) ; clé présente : `history.csv`, `frames.csv`, `*.vtu` inchangés, seul `<outputDir>/probes.csv` en plus. Points séparés par `;`, composantes par `,` (virgule décimale = erreur nommée). |

**Localisation (à l'init, après le maillage).** Pour chaque point, le tétraèdre qui le contient : coordonnées barycentriques λ_a = ∂N_a/∂X · (x − X0[n0]) (a = 1..3, les gradients déjà stockés), λ_0 = 1 − Σ, dedans si toutes ≥ −10⁻⁹ ; premier élément trouvé (point sur une face partagée : l'un des deux, ordre du maillage, déterministe). Aucun tétraèdre : **centroïde le plus proche** avec `[FEM3D] WARNING: probe k lies in no tetrahedron (OUTSIDE the block …)` (ou « inside the box but in no element »). Plusieurs sondes peuvent partager un élément (message `probes j and k share element`, colonnes identiques). Console : une ligne par sonde (élément, centroïde, lc, distance point-centroïde).

**Fichier `probes.csv`** (une ligne par ligne d'historique, même cadence et mêmes instants que `history.csv`, ~2000 lignes, flush à chaque ligne, précision 9 chiffres) : `t`, puis par sonde k les colonnes `p<k>_<champ>` :

| colonnes | contenu |
|---|---|
| `elem` | indice du tétraèdre (0-based, ordre du maillage) |
| `sxx syy szz sxy syz sxz` | contrainte **de la loi** σ_c (Pa, traction > 0) : celle des forces internes, **avant** la viscosité de volume (`bulkViscosity` n'y entre pas, comme pour `svm`/`pm`/`szz`), dans le repère **co-rotationnel** (Biot : contrainte de Biot ; `kinematics = hencky` : Cauchy co-rotationnelle) ; cisaillements tensoriels σ_xy |
| `exx eyy ezz exy eyz exz` | déformation **passée à la loi** : Biot sym(Rᵀ F) − I, ou ln U en hencky ; cisaillements tensoriels ε_xy (= γ_xy/2) |
| `p q` | p = −tr σ/3 (**compression > 0**), q = √(3/2 s:s) — le plan des méridiens |
| `s1 s2 s3` | contraintes principales **décroissantes** (s1 ≥ s2 ≥ s3, traction > 0) |
| `trEpl eplEq epvEq` | tr(ε_pl) (dilatance > 0, compaction < 0), ε_pl équivalente √(2/3 dev ε_pl : dev ε_pl), et `epvEq` propre à la loi (dpr/saksala : déformation viscoplastique équivalente ; cdp : ε_c^pl) |
| `dt dc pc` | d_t = `MatState::D` (dpr/saksala : D de Rankine ; cdp : d_t), d_c = `MatState::Dc` (compDamage crackband ω_c ; cdp : d_c), pc = pression de cap courante (**0 si pas de cap**) |
| `detF eroded` | det F (J) et drapeau de suppression (0/1). Élément supprimé : σ = 0, ε = dernière valeur calculée (soupape det F : celle du pas précédent), les variables d'état restent à leur dernière valeur |
| `epsTpl epsCpl dcdp` | **cdp seulement** : ε_t^pl, ε_c^pl de Lee-Fenves et d combiné |
| `Dv1 Dv2 Dv3` | **dpdfh seulement** : les trois endommagements directionnels D_i |

Champs sans objet pour une loi (elastic : tout ce qui est plastique ; dpr sans cap : `pc` ; dpr : `dc` = 0 sans `compDamage`) : écrits **0**. Composantes de la loi ≠ contrainte globale : pour un élément qui tourne beaucoup (lèvre du cratère) les composantes co-rotationnelles ne sont pas celles du repère du bloc — p, q, s1-s3 sont invariants et ne souffrent pas de ce choix.

**Coût.** Sans la clé : un booléen par appel de `elementForces`. Avec : un test d'indice par élément et une copie de deux matrices 3×3 pour les éléments sondés ; un thread par élément dans la région OpenMP (aucune course). Localisation O(n_el × n_sondes) une fois.

**Banc falsifiant** (`etude_lois_fem/bitid_w18/selftest_probes/`, `SELFTEST_probes.md`, `analyse_probes.py` → `resultats.txt`, PASS 14/14) : (i) bloc élastique 8³ mm sous confinement isotrope 30 MPa (tension, `pullV = 0`, `pullDelay > T`, `topPressure` = `confiningPressure`, `gripLateralFree`) : p = 29,990 MPa (−0,03 %), q = 5·10⁻¹¹ MPa, ε_ii = −P/3K à 0,03 % ; **contrôle qui doit échouer** sans `topPressure` : p = 20,0 = 2P/3, q = 30,0 = P, rejeté ; (ii) compression élastique 4×4×8 mm `triaxStats` : `p1_szz` = `sigZZmid` à 0,21 % max (0,01 % moyen) sur 2415 lignes ; (iii) sonde à (4, 4, 12) mm hors du bloc de 8 mm : WARNING et élément le plus proche dans la couche du dessus (centroïde z = 7,75 mm) ; (iv) trois sondes à 10⁻⁷ m : même élément, 26 colonnes identiques sur 3335 lignes. Figures : `etude_lois_fem/matrice_v2/fig_probes.py` (une figure par point : plan p-q avec le méridien du deck — dpr linéaire ou puissance, apex, cap ; cdp effectif —, σ_zz-ε_zz, d_t/d_c/ε_pl depuis le contact ; 3 runs au plus, PDF + PNG) ; ligne à ajouter aux decks D : `matrice_v2/probes_D.txt` (cinq points : pôle −1 et −5 mm, lèvre, bord du bulbe, champ lointain).

#### Endommagement de traction à direction figée `tensionDamage = scalar | fixed` (2026-09-05, `rockim_f2w19.exe`, opt-in)

**Clés.** `tensionDamage = scalar` (défaut = chemin actuel, bit-identique clé absente) `| fixed` ; `tensionShearRetention = β` ∈ [0, 1] (défaut 1, lu seulement avec `fixed`, orpheline refusée). Lois `dpr` et `saksala` (noyau `PlasticDamageLaw`) ; **refusée** pour `saksala2011` (port fidèle : endommagement piloté par la déformation viscoplastique, sans bande (f_t, G_f, l_c)), `cdp` (tables propres) et `dpdfh` (déjà directionnel). Faute de frappe = exception. Un binaire antérieur à w19 **ignore la clé en silence** et rend le scalaire.

**Pourquoi.** Demande de Fernando (05/09) : le d_t de Rankine est un scalaire qui dégrade **toute** la partie tendue du tenseur — un élément fissuré en x perd aussi sa raideur en y (banc (a) : 0,1 E après d = 0,9). Le *fixed crack model* (Rashid 1968 ; Rots & Blaauwendraad 1989, « smeared crack » à direction fixe) remplace ce scalaire par un endommagement **par direction**, repère figé au premier dépassement de f_t, raideur perdue **normalement** au plan de fissure et **conservée parallèlement** — la géométrie de la loi `dpdfh` (D_i dans le repère figé `Dfh::eul`), avec la **cinétique de la bande de fissuration de Rankine existante** (f_t, G_f, l_c) et non l'obscuration.

**Formules** (`src/MatLaw.cpp`, `PlasticDamageLaw::fixedCrackUpdate` et l'assemblage nominal ; état `MatState::Fcm {nAct, d[3], kap[3], eul[3]}`). Tenseur pilote T = σ_eff/E (`rankineDrive = stress`) ou ε − ε_p (`strain`), la **même grandeur** que le cut-off scalaire. *Amorçage* : la plus grande valeur propre de T restreinte au **complément orthogonal** des directions déjà figées dépasse k_0 = f_t/E → la direction est figée : 1re direction = repère principal complet (n1 = vecteur propre max, n2, n3 = les deux autres, main droite), stocké en angles d'Euler ZYX (`fcmEuler`/`fcmFrame`, mêmes formules que `dfhk::eulr`/`reul`) ; 2e = rotation de (n2, n3) autour de n1 vers le vecteur propre max du 2 × 2 restreint (θ = ½ atan2(2 T_23, T_22 − T_33)) ; 3e = n3 seul. *Croissance* : κ_i = max(κ_i, n_iᵀ T n_i), d_i = 1 − (k_0/κ_i) exp(−(κ_i − k_0)/k_f), k_f = G_f/(l_c f_t) − k_0/2 — **identique** à la formule scalaire (ouverture u_i = κ_i l_c), monotone. *Contrainte nominale* : S = Rᵀ σ_eff R (colonnes de R = n_i) ; S_ii ← (1 − d_i) S_ii si S_ii > 0 (fissure ouverte ; **composante normale compressive intacte** : unilatéral, l'équivalent du split spectral du scalaire) ; S_ij ← (1 − β max(d_i^ouv, d_j^ouv)) S_ij avec d^ouv = d si la fissure est ouverte (S_ii > 0), 0 sinon — β = 1 : convention min(f_i, f_j) de la VUMAT DP-DFH (fissure ouverte = pas de cisaillement transmis, fermée = transfert total) ; β < 1 : *shear retention factor* de Rots ; σ_nom = R S_n Rᵀ. `compDamage = crackband` : ω_c reste sur la partie spectrale négative de σ_eff, soustraite après l'opérateur figé (sans d_i : (1 − ω_c) σ⁻ + σ⁺ comme avant). *Compteurs* : wDamT += ½ ⟨S_ii⟩²/E dd_i sur les trois directions (le cisaillement retenu n'est pas compté) ; wDamC, wPlas inchangés. *Sorties* : `s.D = max_i d_i` (history, VTU `damage`, V_D09, canal spall `erodeD`), `s.kappa = max_i κ_i` (`kapDP`) ; **champs VTU `dFix1..3`** ajoutés au mode courant (stats / vtkCap honorés) et **colonnes `p<k>_dFix1..3`** dans `probes.csv`. Plasticité DP, cap, méridien, apex, érosion (`erodeD`, `erodeEpv`, `erodeDc`, `erodeWfrac`), soupapes : inchangés.

**Ce qui change dans le code** (originaux dans `etude_lois_fem/bitid_w19/orig/`) : `include/rockim/MatLaw.hpp` (sous-état `Fcm`, fiche), `src/MatLaw.cpp` (`BrickOpts::fixedCrack/shearRet`, branche `if (!br_.fixedCrack)` autour du Rankine scalaire — texte du scalaire inchangé —, branche `if (br_.fixedCrack)` avant l'assemblage nominal, `fixedCrackUpdate`, validation des clés dans `MatLaw::make`, `fixedCrackSelftest`), `src/main.cpp` (`rockim selftest-fixed [out.csv]`), `include/rockim/Fem3dSolver.hpp` + `src/Fem3dSolver.cpp` (`fixedDam_`, branche VTU, colonnes probes). `rockim matpoint` hérite de la clé par `MatLaw::make` (chemin `tension` avec `tensionDamage = fixed`).

**Bit-identité** (`etude_lois_fem/bitid_w19/comparaison_w19.txt`, `bitid_w19.sh`, OMP 2) : (a) `PQ_cdpI_P020_court` (cdp) sans clé, w19 vs sorties conservées de w18 : history e45a7cbc0609c9f8, frames 5ddf25b84b9c15f1, vtu[4] 180a5a99311bbe44, cmp 6/6, pic 7791,04 N — chaîne w14 → w19 ; (b) **deck court dpr** `bitid_w19/C_T1_R_P000_court_grid.cfg` (= C_T1_R_P000 à T 100 µs, 2 images, maillage interne 24 × 24 × 16 mm, 55 296 tets : le noyau touché) sans clé, w18 vs w19 : e23cbec3a88c2e6f / 08d05ce828a8e2af / bd067aa893e231ac, cmp 6/6, pic 24 623,9 N ; (c) `tensionDamage = scalar` explicite ≡ sans clé (cmp 6/6) ; (d) contrôle `fixed` : 6/6 DIFFÉRENT, pic 23 348,9 N (−5,2 %), wDamT 0,166 vs 0,122 J, 2686 érodés (canal `erodeWfrac`) contre 0, damage ≥ 0,9 sur 17 187 vs 23 234 tets, dFix1 > 0 sur 50 429 tets, deux directions sur 33 094, trois sur 8288, `damage` = max(dFix) à 0,0 ; 76,6 s contre 102,3 s (la rotation R S Rᵀ est moins chère que la décomposition spectrale du scalaire).

**Banc falsifiant** (`etude_lois_fem/bitid_w19/selftest_fixed/`, `SELFTEST_fixed.md`, `rockim selftest-fixed`, **33 PASS / 0 FAIL**) au point matériel, carte Red Bohus dpr (`dpApex`, `dpTension = off`, `rankineDrive = stress`), pilotage en déformation (état uniaxial ε = e n nᵀ − ν e (I − n nᵀ)) : (a) traction x jusqu'à d = 0,9, décharge, traction y : pente σ_yy/ε_yy = **1,000 E** (fixed) contre **0,0997 E = (1 − D) E** (scalar) — la réponse scalaire est **rejetée** par le critère E à 1 %, c'est le test qui sépare les deux ; même cinétique en uniaxial (5·10⁻¹⁶), σ_xx nominal ≤ 3·10⁻¹⁷ f_t pendant la traction y, d1 inchangé bit pour bit ; puis y s'amorce à E ε_yy = f_t (1,02 k_0), n2 = y, d1 inchangé ; (b) traction x puis compression x : E dans les deux (−27,000 MPa = E ε) ; (c) traction à 45° : eul = (45,000, 0, 0)°, n1 = (0,7071, 0,7071, 0), σ_45 = R_z σ_x R_zᵀ à 2,6·10⁻¹⁵ ; (d1) tel qu'écrit (x maintenu, ε_yy croissant, σ_zz = 0) : σ_xy = 0 par symétrie (axes = repère), y s'amorce à S_22 = f_t (1,006), d1 0,5 → 0,763 par couplage de Poisson, σ_yy nominal fin identique dans les deux modèles (8,435 MPa) ; (d2) **donnée** : x maintenu à d = 0,5 puis traction croissante selon 45° — cisaillement transmis sur le plan de fissure σ_xy^nom/σ_xy^eff = **0,242** (β = 1, = 1 − max(d1 0,758, d2 0,505)) et **0,924** (β = 0,1), direction principale nominale désalignée de **13,9°** (β = 1) et 13,4° (β = 0,1) de l'effective (31,7°) : le *stress locking* de Rots — β règle le cisaillement transmis, pas le désalignement ; scalar 0,183 ; (e) énergie dissipée en traction : ∫σ dε × l_c/G_f = 0,9999 et wDamT × l_c/G_f = 1,0008 dans les deux modèles.

**Points à relire.** (1) `wDamT` **somme** les trois directions (spec) : un élément fissuré dans trois directions rend jusqu'à 3 G_f/l_c et le canal `erodeWfrac` se déclenche plus tôt qu'en scalaire (deck court : 2686 contre 0) — garder, ou passer le critère sur max_i w_i ; (2) le repère est celui du pas de dépassement (comme la DP-DFH) : à un pas explicite près ; (3) directions 2 et 3 : cherchées dans le complément orthogonal (pas de re-rotation de n1) ; (4) les angles eul2/eul3 d'une fissure unique ne sont pas contraints (n2, n3 = base quelconque du plan propre) — lire `dFix1` et n1 ; (5) deck de démonstration `matrice_v2/D_T1_fixed_P000.cfg` (= C_T1_R_P000 + `tensionDamage = fixed` + sondes de `probes_D.txt`), **non lancé**.
