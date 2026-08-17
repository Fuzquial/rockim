# Guide pratique de rockim

*Mini-code C++ de simulation d'impact et de coupe sur roche — FEM / DEM / FDEM,
2D et 3D, multi-lois. Guide utilisateur en français ; le `README.md` (anglais)
reste la référence détaillée des modèles et des chiffres de vérification.*

---

## 1. Où se trouve rockim

Le dépôt vit en **bundles git** (pas de remote pour l'instant) : le bundle le
plus récent fait foi, on le clone pour travailler.

| Emplacement | Contenu |
|---|---|
| `phd_geothermie\FDEM\rockim\` (base OneDrive, **l'archive de référence**) | bundle le plus récent + `HANDOFF_*.md` (état des chantiers) + docs et rapports ; versions superseded dans `old\` ; `rockim_gbm.tar.gz` = la lignée GBM séparée |
| clone de travail : `simulations\FDEM\rockim\rockim_p1\` (branche `article-exact`) | déplacé là le 2026-08-17 (il était dans `Downloads`) ; `simulations` n'est **pas** synchronisé OneDrive, donc pas de fichier verrouillé pendant la compilation |

Cloner un bundle (se placer DANS le dossier du bundle) :

```bash
git clone rockim_2026-08-17.bundle -b article-exact rockim_p1
```

**Repartir de zéro sur une machine neuve — procédure vérifiée de bout en bout
le 2026-08-17** (clone vierge → compilation → suite à 15/15) :

1. récupérer le bundle dans `phd_geothermie\FDEM\rockim\` et le cloner comme
   ci-dessus. Le clone ne contient **ni exécutable, ni maillage, ni sortie** :
   uniquement 12 sources, 17 en-têtes, les configs, les outils et les docs ;
2. placer **Eigen 3.4.0 À CÔTÉ du clone**, pas dedans : le script de compilation
   le cherche en `..\eigen-3.4.0`. C'est le seul piège de l'installation, et il
   s'est manifesté le 2026-08-17 quand le clone a été déplacé sans Eigen ;
3. `pip install gmsh numpy matplotlib scipy` pour les outils Python ;
4. compiler — **107 s** mesurées. Sur une machine dont Visual Studio n'est pas
   au chemin par défaut, `build_chk.cmd` échoue (chemin de `vcvars64.bat` en
   dur) : passer par **CMake**, qui trouve Eigen seul ;
5. **régénérer les maillages** : les `.msh` ne sont PAS dans le dépôt (156 Mo,
   exclus par `.gitignore` car reproductibles). La ligne de commande de
   génération figure en tête de chaque config concernée ;
6. valider avant tout calcul :
   `python tools\verify_suite.py --exe C:\chemin\ABSOLU\vers\rockim.exe` — 15
   tests. Le chemin **doit être absolu**, le runner changeant de répertoire de
   travail.

L'arborescence du dépôt :

```
rockim/
  src/            le code (un .cpp par solveur + MatLaw.cpp + Tessellation*.cpp)
  include/rockim/ les en-têtes (modèles documentés en tête de fichier)
  configs/        toutes les configs : démos + bancs + verifications
  meshes/         maillages Gmsh generes (mesh = file)
  tools/          verify_suite.py (LA suite de non-regression), rockim_gui.py,
                  make_unstructured_mesh.py (maillages multi-corps),
                  crater_metrics.py, export_abaqus.py, bayes_bench.py, ...
  specs/          chantiers spec-kit (spec + plan + tasks)
  .specify/memory/constitution.md   LA discipline du projet (a lire d'abord)
  README.md       référence des modèles (EN) ; GUIDE_rockim.md = ce guide
  DOCUMENTATION_rockim.md           référence complète clés/sorties
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
OpenMP, arrêt, **bouton « Suite de vérification »** qui lance la suite
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

## 4. Pas à pas : une simulation d'impact multi-corps de A à Z

*Exemple fil rouge : l'impact d'un insert sphérique en carbure de tungstène
(R = 11 mm, 8 m/s) sur un bloc de granite 120×120×120 mm — le « banc » du
pilier performance, monté le 2026-08-14.*

### 4.1 Le modèle = le maillage (deux corps, zéro joint entre eux)

Depuis V1, un « modèle » multi-corps est un maillage Gmsh où chaque **volume
physique** devient un corps : joints cohésifs À L'INTÉRIEUR de chaque corps,
**aucun joint entre corps** — l'interaction passe par le contact général.
Le générateur fourni fait tout :

```bash
python tools/make_unstructured_mesh.py bench1 0.12 0.12 0.12 0.011 0.00005 0.0045 0.0028 meshes/mon_banc.msh 1
```

soit `bench1 W D H R gap h hIns out.msh [seed]` :

| argument | rôle | effet sur le coût |
|---|---|---|
| `W D H` | bloc de roche [m] | volume → nombre de tets |
| `R` | rayon de l'insert sphérique [m] | — |
| `gap` | jeu initial insert/roche [m] | fixe l'instant du contact : **t = gap / v** |
| `h`, `hIns` | taille d'élément roche / insert [m] | tets ∝ 1/h³ ; **dt ∝ h_min** — le coût total explose en ~1/h⁴ |
| `seed` | graine du maillage non structuré | reproductibilité |

Ordres de grandeur mesurés (ce banc) : h = 4,5/2,8 mm → 82 k tets ;
h = 2/2 mm → 842 k tets. Le résumé du générateur donne le h inscrit
min/med/max — c'est le **min** qui fixera le pas de temps.

### 4.2 Dimensionner la fenêtre T AVANT tout (leçon payée 4 h)

Trois temps à poser sur une enveloppe avant d'écrire la config :

1. **arrivée au contact** : t₀ = gap / v (ex. 0,05 mm / 8 m/s = 6,25 µs) ;
2. **durée du contact** (Hertz, ordre de grandeur) :
   t_c ≈ 2,87 (m² / (R E*² v))^(1/5) — pour l'insert WC de 79 g sur granite :
   **~90 µs**. E* combiné : 1/E* = (1−ν₁²)/E₁ + (1−ν₂²)/E₂ ;
3. **ce qu'on veut voir** : chargement seul → T ≈ t₀ + 15 µs suffit ;
   cycle complet avec REBOND → **T ≳ t₀ + t_c + marge** (ici 120 µs).

Le piège historique (banc P1 v1, 2026-08-14) : gap 0,5 mm et T = 20 µs
→ contact à 62,5 µs, **jamais atteint** — 4 h de calcul d'approche à vide.
D'où la règle suivante.

### 4.2 bis La règle de FAISABILITÉ : le cas peut-il casser ? (leçon payée 2 h 30)

**Trois nombres à poser avant tout run de fissuration.** Dimensionner T ne suffit
pas : un cas peut tourner jusqu'au bout, fermer son bilan d'énergie à 1e-5 %, et
ne rien casser — non par bug, mais parce que le matériau ne *peut pas* casser à
l'échelle où on le sollicite.

La grandeur qui décide est la **longueur de la zone cohésive** :

```
ℓ_cz ≈ E · Gf / ft²
```

C'est la taille du processus de fissuration. Elle doit être encadrée :

```
2 · dx  <  ℓ_cz  <  a
```

- **borne haute** — `ℓ_cz < a`, avec `a` la taille de la zone CHARGÉE (pour un
  contact sphère/plan, le rayon de Hertz `a = √(R·δ)`, **pas** l'étendue du champ
  de contrainte visible, cinq fois plus grande) : si la zone cohésive ne tient
  pas dans la zone chargée, aucune fissure ne peut localiser. L'endommagement
  s'étale et plafonne à quelques pourcents ;
- **borne basse** — `ℓ_cz > 2·dx` : il faut au moins deux éléments pour résoudre
  le processus, sinon la rupture se concentre sur un élément et dépend du
  maillage.

**Le contre-exemple, mesuré le 2026-08-17.** Banc percussion, `ft = 10 MPa`,
`Gf = 70 J/m²`, `E = 50 GPa` → **ℓ_cz = 35 mm** pour un rayon de contact de
**1,45 mm** : rapport 24. Résultat sur six runs : **4 joints rompus sur 158 423**,
restitution 0,90, et la fissuration consommant 0,55 % de l'énergie contre 12 %
pour l'amortissement numérique. Ni le maillage (facteur 6 balayé), ni
l'hétérogénéité de Weibull (m = 6 puis 24), ni l'insertion adaptative, ni la
vitesse d'impact n'y changeaient quoi que ce soit — et frapper plus fort était
exclu : la roche était **déjà sollicitée à 49 fois sa résistance en traction**.

Avec `ft = 87 MPa`, `Gf = 150 J/m²` → ℓ_cz = 1,00 mm, et un maillage à 0,46 mm
dans la zone de contact (ℓ_cz/dx = 2,2) : **731 joints rompus**, cratère localisé
à l'échelle du contact, fissuration à parité avec l'amortissement.

**Corollaire sur (ft, Gf).** Les deux ne se choisissent pas séparément : seul leur
rapport `Gf/ft²` fixe ℓ_cz. Un matériau plus TENACE a une zone cohésive plus
GRANDE, donc plus facile à résoudre. À `ft` élevé (échelle de l'indentation), il
faut donc un `Gf` élevé pour rester calculable — c'est contre-intuitif mais c'est
l'algèbre. Et l'ouverture critique de rupture d'un joint vaut
`dnF ≈ 2,59 · Gf/ft` (facteur `1/∫f(D)dD` en adoucissement de Yan) : c'est une
propriété du **matériau**, indépendante du chargement, ce qui explique pourquoi
taper plus fort ne rachète jamais un mauvais couple.

### 4.3 La règle du smoke test — TOUJOURS

Avant tout run > 30 min : le MÊME cas (mêmes clés, même physique) sur un
maillage réduit, en multithread, quelques minutes. Verdicts à cocher :

- le contact démarre à t₀ prévu (`grpFz` décolle) ;
- ça casse (`nBroken` > 0) si le cas doit casser ;
- bilan B4 : résidu < 1 % (en pratique ~1e-5 %), dashpot et Cundall ≤ 0 ;
- le wall time du smoke, extrapolé (× N_tets × N_pas), confirme le budget.

### 4.4 La config, bloc par bloc

```
mode = fdem3d              # solveur FDEM 3D (joints cohesifs + contact)
scenario = percussion
mesh = file                # maillage Gmsh multi-corps
meshFile = meshes/mon_banc.msh
T = 1.2e-4                 # cf. 4.2 !
frames = 6                 # frames VTU (chaque frame 842k = ~360 Mo en texte)

rho = 2650                 # bloc materiau = defauts globaux (la roche)
E = 50e9
nu = 0.25
ft = 10e6
cohesion = 25e6
frictionDeg = 40
Gf = 70
gfShearFactor = 10

phases = rock insert       # une phase par volume physique HOMONYME du .msh
phase.rock.fraction = 0.5  # fractions obligatoires mais SANS effet ici
phase.insert.fraction = 0.5
phase.insert.rho = 14500   # l'insert : carbure de tungstene, resistances
phase.insert.E = 600e9     # hors d'atteinte -> il reste elastique
phase.insert.ft = 400e6
phase.insert.cohesion = 800e6
phase.insert.Gf = 400

jointPenaltyFactor = 20

toolShape = none           # PAS d'outil rigide : l'insert est MAILLE
groupVel.insert = 0 0 -8   # vitesse initiale du corps 'insert' [m/s]
trackGroup = insert        # suivi -> colonnes grpZ/grpVz/grpFx..Fz/grpSzz

dampingLocal = 0.05
absorbing = all            # frontieres de Lysmer
dtFactor = 0.15

contact = potential        # potentiel de Munjiza (conservatif) | penalty
gcActivation = adaptive    # activation adaptative (Fukuda) : x2 et plus
```

Rappels : **unités SI** (m, s, Pa, kg), **point décimal** obligatoire.

### 4.5 Lancer

```bash
set OMP_NUM_THREADS=1      # 1 = bit-identique, comparable aux refs
rockim.exe configs\mon_banc.cfg out_mon_banc
```

`OMP_NUM_THREADS` non posé = tous les cœurs (perf, mais plus de
bit-identité : l'ordre de sommation flottante change). Deux précautions :

- **budget** : coût ≈ N_tets × N_pas × (µs/tet/pas de la machine) ;
  N_pas = T / dt et dt ≈ dtFactor × h_min / c_max (c = √(E/ρ) du matériau
  le plus raide — souvent l'insert). Repères mesurés (Core Ultra, MSVC,
  potentiel+adaptatif) : ~1,2 µs/tet/pas à 18 threads ;
- **l'exe d'un run en cours est verrouillé par Windows** : pour continuer à
  compiler pendant qu'un run tourne, lancer le run sur une COPIE de l'exe
  (ou compiler sous un autre nom, `/Fe:rockim_dev.exe`).

### 4.6 Suivre un run en cours

`history.csv` s'écrit en continu dans le dossier de sortie — les colonnes
utiles d'un impact : `t`, `grpZ/grpVz` (cinématique du corps suivi),
`grpFz` (la F-δ en direct), `grpSzz` (jauge sous l'insert), `nBroken`,
`eEl/eJnt/eGc/eFric/eCund/eLys` (bilan B4 par sous-système). Avancement =
dernier `t` / T. Attention : valeurs imprimées en 6 chiffres significatifs —
les positions évoluent « en escalier » de 1 µm, c'est l'arrondi d'écriture,
pas la physique.

### 4.7 Dépouiller

1. **Le résumé stdout d'abord** : résidu B4 (doit être ≪ 1 %), dashpot
   [OK, dissipative], joints cassés (traction/cisaillement), `potential
   stats` (paires, tGrid/tLoop), fragments, wall time.
2. **Restitution** : e = |vz sortie| / |vz entrée| du corps suivi
   (repère : e = 0,71 sur bench1, 2D comme 3D).
3. **Cratère** : `python tools/crater_metrics.py out_mon_banc` (rayon,
   profondeur, fissures radiales, bras endommagés — multi-corps géré).
4. **ParaView** : frames `fdem3d_XXXX.vtu` (+ `_joints` : `damage`, `type`).
5. Courbes : la GUI (`tools/rockim_gui.py`) ou matplotlib sur `history.csv`.

### 4.8 Archiver

Jalon validé → `git bundle create rockim_<date>.bundle --all` → copie dans
`phd_geothermie\FDEM\rockim\` (l'ancien bundle part dans `old\`), docs à
jour DANS LE MÊME COMMIT (constitution §7).

## 5. Écrire une config — référence des clés

Format `clé = valeur`, une par ligne, `#` = commentaire. **Unités SI partout**
(m, s, Pa, kg — PAS le mm-MPa d'Abaqus) et **point décimal obligatoire** (`0,5`
est rejeté avec un message nommant la clé — protection locale FR). Les clés
inconnues sont ignorées, mais les valeurs invalides et les combinaisons
incohérentes (ex. `phases` sans `mesh = voronoi`) arrêtent le run avec un
message explicite.

### 5.1 Bloc commun

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

### 5.2 Matériau (partagé par tous les modes)

`rho, E, nu, ft, cohesion, frictionDeg, Gf, gfShearFactor` — validation
stricte (E, rho, ft, cohesion, Gf > 0 ; nu ∈ [0, 0.5) ; frictionDeg < 89).

### 5.3 Outil et chargement

| clé | rôle |
|---|---|
| `toolRadius`, `toolMass`, `impactSpeed` | insert sphérique/disque libre (percussion) |
| `toolShape` | 2D : `disc` \| `flat` ; 3D : `sphere` \| `flat` (poinçon cylindrique à bout plat, rayon `toolRadius`) ; `cutDepth`, `cutSpeed` pour la coupe |
| `contactMu`, `contactXi` | frottement et amortissement du contact outil |
| `pullV`, `pullRamp` | vitesse du mors (traction/compression) ; **rampe cosinus** [s] — sans rampe, le transitoire casse au mors quel que soit le matériau |
| `gripLateralFree` | mors sans frottement latéral (essais uniaxiaux propres) — rampe et mors libres disponibles dans TOUS les modes à scénario tension |

## 6. Les lois de comportement (mode `fem3d`)

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

## 7. Le mode GBM (`mesh = voronoi` — FDEM 2D **et** 3D)

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

## 8. Sorties et post-traitement

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

## 9. Vérifications — à relancer après toute modification

**La** suite de non-régression est `tools/verify_suite.py` (~48 repères) :

```bash
python tools/verify_suite.py --exe rockim.exe --tier fast
```

`--tier fast` = 12 tests (~50 s) : selftests des lois (Saksala 2011, DP-DFH),
conservation du potentiel 2D/3D, intégrale de Yan, barre FEM, tensions
DEM/Voronoï, et les **zeroload** (charge nulle : 0 joint cassé, travail de
contact = 0 exact — ils attrapent ce que rien d'autre ne voit).
`--update-refs` re-baseline les repères. **Les réfs sont PAR PLATEFORME** :
MSVC et libstdc++ divergent à graine égale (Voronoï/Weibull) — ne jamais
comparer un pic MSVC à une réf Linux, ni des pics à threads différents.
Constat 2026-08-14 : la suite fast passe 12/12 sous MSVC sans re-baseline.

## 10. Pièges connus

- **Virgule décimale** → erreur explicite (voulu) ; corriger la config.
- **Fenêtre T vs gap** : vérifier t_contact = gap/v ≪ T AVANT de lancer
  (cf. §4.2 — 4 h payées le 2026-08-14 pour un banc qui ne touchait jamais).
- Un `.exe` fraîchement écrasé peut être verrouillé quelques secondes
  (antivirus) — relancer. Et l'exe d'un **run en cours** est verrouillé tout
  du long : compiler sous un autre nom pendant ce temps (`/Fe:rockim_dev.exe`).
- `git clone <bundle>` échoue si on n'est pas dans le dossier du bundle.
- OneDrive verrouille des fichiers pendant la compilation : ne jamais compiler
  depuis `phd_geothermie\`. Le clone de travail vit dans `simulations\`, qui est
  un dossier local ordinaire.
- `meshMirror = false` restitue l'ancien maillage fem3d à l'identique.
- Les paramètres des démos sont des ordres de grandeur NON calibrés ; la
  calibration est le rôle du banc bayésien (`tools/bayes_bench.py`).
- rockim est le **bac à sable** (exploration, calibration, figures) — la
  production 3D calibrée reste Abaqus + VUMAT.
