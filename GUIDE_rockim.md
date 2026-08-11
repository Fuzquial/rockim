# Guide pratique de rockim

*Mini-code C++ de simulation d'impact et de coupe sur roche — FEM / DEM / FDEM,
2D et 3D, multi-lois. Guide utilisateur en français ; le `README.md` (anglais)
reste la référence détaillée des modèles et des chiffres de vérification.*

---

## 1. Où se trouve rockim

| Emplacement | Contenu |
|---|---|
| `phd_geothermie\FDEM\rockim\` (base OneDrive, **l'archive de référence**) | `rockim_gbm.tar.gz` (source complète), `FICHE_rockim.md` (historique + verdicts), figures, scripts |
| `simulations\FDEM\rockim\` (dossier de travail, miroir) | copie identique |
| `Downloads\rockim_gbm.tar.gz` | copie de commodité ; `Downloads\rockim.tar.gz` = l'archive d'ORIGINE avant extensions |

L'arborescence de la source, une fois extraite :

```
rockim/
  src/            le code (un .cpp par solveur + MatLaw.cpp + Tessellation.cpp)
  include/rockim/ les en-têtes (modèles documentés en tête de fichier)
  configs/        toutes les configs : démos + vérifications (verify_*.cfg)
  tools/          rockim_gui.py (interface), bayes_bench.py (banc bayésien),
                  export_abaqus.py (maillage+champ -> .inp mm-t-s-MPa),
                  plot_results.py, make_gif.py
  README.md       référence des modèles (EN) ; GUIDE_rockim.md = ce guide
  CMakeLists.txt
```

## 2. Compiler

**Windows / MSVC (testé)** — depuis le dossier `rockim/` extrait :

```bash
cmd /c '"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" && cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX /I include /I ..\eigen-3.4.0 src\*.cpp /Fe:rockim.exe'
```

Seule dépendance : **Eigen** (headers seuls). Si absent :
`curl -L -o eigen.tar.gz https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz && tar -xzf eigen.tar.gz` puis pointer `/I` dessus.

**CMake (Linux/Mac/Windows)** : `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j` — Eigen est téléchargé automatiquement si non installé, OpenMP détecté s'il existe.

## 3. Lancer

```bash
./rockim.exe configs/fem3d_percussion_saksala2011.cfg mon_dossier_sortie
```

Premier argument : la config ; second (optionnel) : le dossier de sortie.

**Interface graphique** (tkinter + matplotlib, rien d'autre) :

```bash
python tools/rockim_gui.py
```

Trois zones : configs (édition + sauvegarde), lancement (exe, dossier, threads
OpenMP, arrêt, **bouton « Suite de vérification »** qui enchaîne les 11 verify_*
et affiche PASS/FAIL), et tracés (courbe F–δ, historiques, coupe médiane) pour
tout dossier `out_*`.

**Autotests des lois portées** (rejouent les chemins des harnais VUMAT) :

```bash
./rockim.exe selftest-saksala2011
./rockim.exe selftest-dpdfh
```

**Variables d'environnement utiles** : `OMP_NUM_THREADS` (nombre de threads ;
`1` = résultats bit-identiques au build série ; N fixé = déterministe) ;
`ROCKIM_PROF=1` (profil par pas en fin de run, mode fdem).

## 4. Écrire une config

Format `clé = valeur`, une par ligne, `#` = commentaire. **Unités SI partout**
(m, s, Pa, kg — PAS le mm-MPa d'Abaqus) et **point décimal obligatoire** (`0,5`
est rejeté avec un message nommant la clé — protection locale FR). Les clés
inconnues sont ignorées, mais les valeurs invalides et les combinaisons
incohérentes (ex. `phases` sans `mesh = voronoi`) arrêtent le run avec un
message explicite.

### 4.1 Bloc commun

| clé | rôle (défaut) |
|---|---|
| `mode` | `fem` \| `fem3d` \| `dem` \| `dem3d` \| `fdem` \| `fdem3d` |
| `scenario` | `percussion` \| `shear` \| `tension` (traction ; **pullV < 0 = compression uniaxiale**) — les trois disponibles en 2D ET en 3D |
| `T`, `frames` | durée physique [s], nombre de frames VTK |
| `W, H` (+ `D` en 3D) | dimensions du bloc [m] ; `thickness` en 2D |
| `nx, ny` (+ `nz`) | découpage du maillage |
| `seed` | graine du maillage (jitter, Voronoï) |
| `dtFactor` | fraction du pas critique (0.15–0.3 typique) |
| `dampingLocal` | amortissement de Cundall (0.05 dynamique ; ~0.1 en quasi-statique — 0.7 biaise les mesures de contrainte, cf. README) |
| `absorbing` | `none` \| `sides` \| `all` (frontières de Lysmer) |

### 4.2 Matériau (partagé par tous les modes)

`rho, E, nu, ft, cohesion, frictionDeg, Gf, gfShearFactor` — validation
stricte (E, rho, ft, cohesion, Gf > 0 ; nu ∈ [0, 0.5) ; frictionDeg < 89).

### 4.3 Outil et chargement

| clé | rôle |
|---|---|
| `toolRadius`, `toolMass`, `impactSpeed` | insert sphérique/disque libre (percussion) |
| `toolShape` | 2D : `disc` \| `flat` ; 3D : `sphere` \| `flat` (poinçon cylindrique à bout plat, rayon `toolRadius`) ; `cutDepth`, `cutSpeed` pour la coupe |
| `contactMu`, `contactXi` | frottement et amortissement du contact outil |
| `pullV`, `pullRamp` | vitesse du mors (traction/compression) ; **rampe cosinus** [s] — sans rampe, le transitoire casse au mors quel que soit le matériau |
| `gripLateralFree` | mors sans frottement latéral (essais uniaxiaux propres) — rampe et mors libres disponibles dans TOUS les modes à scénario tension |

## 5. Les lois de comportement (mode `fem3d`)

Sélection par `law = ...` ; la percussion 3D multi-lois se fait en changeant
UNE ligne. Fissuration = endommagement lissé + érosion (pas de joints).

| loi | modèle | clés propres |
|---|---|---|
| `elastic` | élasticité linéaire | — |
| `dpr` | Drucker-Prager (calé Mohr-Coulomb triaxial) + endommagement de Rankine **régularisé crack-band** (Gf/taille d'élément), rate-indépendant | `erodeD` (0.98), `erodeEpv` (1.5) |
| `saksala` | version SIMPLIFIÉE type Saksala : Perzyna (sur-contrainte linéaire) + cap en pression durcissant + même endommagement | `saksalaEta` [Pa·s], `capP0`, `capH` |
| `saksala2011` | **portage fidèle ligne à ligne de `VUMATS\saksala\vumat_saksala_2011.f90`** (Saksala, IJNAMG 2011) : cône DP non associé, Rankine modifié, cap parabolique raccordé, consistance viscoplastique bi-viscosité, dégradation de cohésion au confinement, écrouissage log du cap, dommage exponentiel, fermeture unilatérale, coin de Koiter. SANS érosion ni régularisation (comme la loi publiée). Vérifié à 8×10⁻¹⁴ contre la référence Fortran | `skBetaDP, skCres, skHdp, skSdp, skSmr, skAt, skBetaT, skPp0, skPtr0, skDcap, skWcap, skNd` — **défauts = Table I du papier** (en Pa) ; il suffit donc de poser le bloc matériau (E=60e9, nu=0.2, ft=13e6, cohesion=37.5e6, frictionDeg=30, rho=2600) |
| `dpdfh` | **portage fidèle de `VUMATS\dfh\vumat_kstdfh.f` — LA loi DP-DFH de la thèse** (Shariati/Saadati/Hild) : DP en compression seule (OPTION-3, retour radial + apex, ψ non associé), **endommagement DFH anisotrope** (3 D_i, repère principal FIGÉ au premier amorçage, seuils = 3 tirages Weibull triés par hash spatial 64 bits DÉTERMINISTE, σ_k = σw(Zeff/V_el)^(1/m), V_el = lc³), croissance par intégrateur racine-cubique de l'obscuration (S, k·c, plancher 1/V_el), nominal (1−D_i) en traction + gating cisaillement (fissure fermée = plein transfert). Sans régularisation ; suppression OFF par défaut (lit de débris). Vérifié à **4,7×10⁻¹²** contre la référence Fortran ifx (`rockim selftest-dpdfh`, 4 chemins dont repère gelé hors axes, tirages bit-identiques) | `dfhBetaDeg, dfhDCoh, dfhPsiDeg, dfhWeibullM, dfhSigW, dfhZeff, dfhK, dfhS, dfhDeld` — **défauts = carte Red Bohus** (§2 de phd/CONTINUUM.md, en SI : β=51,7°, d=153,3 MPa, ψ=15°, m=24, σw=120 MPa @ 1 mm³, k=0,38, S=4π/3) ; poser seulement E=52e9, nu=0.25, rho=2620 (ft/cohesion/frictionDeg du bloc matériau sont IGNORÉS par cette loi). `matWeibullM`/`ftScale` multiplie σw (extension champ, neutre à 1) |

Endommagement unilatéral partout (`dpr`/`saksala` : split spectral — un élément
endommagé porte encore en compression ; c'est le lit de débris — l'équivalent
du choix DELD=1e9 des VUMAT percussion).

**Hétérogénéité** : `matWeibullM = m` tire un facteur de résistance Weibull de
moyenne 1 par élément (mécanisme FIELD du VUMAT ; le papier utilise m = 3) —
tirage indépendant, ou **champ corrélé 3D** avec `strengthCorrLength`
(+ `strengthCorrLengthB`/`strengthCorrAngleDeg` pour une foliation inclinée,
`fieldSeed` indépendant du maillage). Champ visualisable (`ftScale` dans les VTU).

**Géométrie** : `geometry = cylinder` découpe un cylindre Ø min(W,D) dans le
grid (surface courbe absorbante) ; `meshMirror` (défaut on) = Kuhn miroité en
damier (diagonales alternées — sinon les fissures s'alignent sur LA diagonale
globale) ; `meshJitter` désordonne les nœuds.

## 6. Le mode GBM (`mesh = voronoi` — FDEM 2D **et** 3D)

Grains de Voronoï + phases minérales + joints cohésifs classés
(intra-grain / homophase / hétérophase), mêmes clés en 2D (`mode = fdem`)
et en 3D (`mode = fdem3d`, tessellation `Tessellation3` : semis HCP jitté ou
Poisson, clipping par demi-espaces, contraction d'arêtes, maillage tet par
éventails de faces partagés). L'essentiel :

```
mesh = voronoi
grainSize = 0.008          # diametre moyen [m]
grainSeeding = random      # hex | random (Poisson, isotrope — recommande)
lloydIters = 2
refineLevels = 0           # 2D : 0..4 (x4/niveau) ; 3D : 0..2 (x8/niveau)
vertexMergeFrac = 0.25     # EN 3D : 0,25 recommande (0,12 en 2D) — les
                           # eventails des faces minces controlent le dt

phases = quartz feldspar biotite
phase.quartz.fraction = 0.33
phase.quartz.E = 94e9      # toute propriete materiau surchargeable par phase
...
gbAlphaTen = 0.5           # joints de grains = moyenne des phases x alpha
gbAlphaCoh = 0.5           # (gbAlphaGf, gbAlphaE, gbAlphaFric idem)
gbHeteroFactor = 0.8       # malus des joints entre mineraux differents

jointWeibullM = 6          # resistances de joints statistiques (moyenne 1)
strengthCorrLength = 0.008 # 0 = tirage independant ; >0 = champ correle
strengthCorrLengthB = ...  # + strengthCorrAngleDeg : bandes orientees (foliation ;
                           # en 3D : plan de texture incline autour de l'axe y)
fieldSeed = 555            # graine du CHAMP, independante du maillage
```

Le résumé du run donne les fractions atteintes et la **fraction
intergranulaire** de la casse ; ParaView colore par `grain`, `phase`,
`ftScale`, `type` et `damage` des joints.

## 7. Sorties et post-traitement

| fichier | contenu |
|---|---|
| `history.csv` | force outil, position/vitesse, travail, énergies, casse — à chaque ~1/2000 du run |
| `fem3d_XXXX.vtu` / `fdem_XXXX.vtu` (+ `_joints`) | frames ParaView : vonMises, pressure, damage (ω_t), kapDP, epvEq, ftScale, eroded / joints |
| `*_final_*.csv`, `summary` (stdout) | états finaux scriptables + bilans (dont bilan d'énergie outil) |

Scripts fournis (racine archive `FDEM\rockim\`) : `fig_fp.py` (courbes F–δ),
`fig_fem3d.py` (comparaison de lois en coupe), `fig_cyl.py` (vue oblique +
coupes cylindre), `fig_voronoi.py`, `fig_weibull.py` (GBM 2D). La GUI trace
les mêmes choses en un clic.

**Export vers Abaqus** (`tools/export_abaqus.py`) — pour la validation
croisée rockim ↔ Abaqus/Explicit+VUMAT sur maillage IDENTIQUE :

```bash
python tools/export_abaqus.py out_mon_run maillage.inp
```

Lit la frame 0 (fem3d direct ; fdem3d : noeuds dupliqués SOUDÉS — un
maillage Voronoï/GBM devient un C3D4 classique avec un ELSET par phase
minérale), écrit `*NODE/*ELEMENT` en **mm-t-s-MPa** (×1000), NSET
bas/haut, et si un champ `ftScale` existe (matWeibullM) : CSV des
centroïdes + `*INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=1` nodal — le
même champ corrélé σw(x) injectable dans le mécanisme FIELD des VUMAT.
Autocontrôle : volume total des tets réordonnés (jacobien C3D4 positif).

## 8. Vérifications — à relancer après toute modification

`configs/verify_*.cfg` (ou le bouton de la GUI) : tension FDEM 2D (bit-repère
−1,29039 %), onde de barre FEM, tension DEM 2D et 3D, tension FDEM3D grille
(−2,6 %), tension Voronoï 2D (+9…+12 % = tortuosité, normal) et **Voronoï 3D**
(même bande −5…+25 %), compression DP 3D (analytique, ±5 %), tension ft 3D,
sur-contrainte de Perzyna aux deux vitesses (±25 %, linéarité ~2), et
`selftest-saksala2011` (8×10⁻¹⁴ vs Fortran).

## 9. Pièges connus

- **Virgule décimale** → erreur explicite (voulu) ; corriger la config.
- Un `.exe` fraîchement écrasé peut être verrouillé quelques secondes
  (antivirus) — relancer.
- `meshMirror = false` restitue l'ancien maillage fem3d à l'identique.
- Les paramètres des démos sont des ordres de grandeur NON calibrés ; la
  calibration est le rôle du banc bayésien (`tools/bayes_bench.py`).
- rockim est le **bac à sable** (exploration, calibration, figures) — la
  production 3D calibrée reste Abaqus + VUMAT.
