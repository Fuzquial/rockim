# WP6 — μ de contact dégradé (`contactResidualMu`) et campagne pulvérisation × coulomb

Date : 2026-08-28. Conçu en session conteneur (cartographie 4 lecteurs + 2 plans
concurrents + revue adverse), sur la base de Yang et al., IJRMMS 206 (2026)
106660 (granite Kuru Grey, données Mines Paris-PSL / Aising).

## 1. État des lieux — ce qui EXISTE déjà (vérifié, commits à l'appui)

Le modèle de pulvérisation du papier a DEUX ingrédients. Le premier est déjà
dans rockim, le second n'y est qu'à moitié :

| ingrédient Yang 2026 | état rockim | où |
|---|---|---|
| éq. 3-4 : D(δm) linéaire, δm = h_e·ε_vm, plafond Dmax, σ ← Cd(1−D)σ | **FAIT** (WP1, commit 792a9613, 2026-08-22 ; validation 3 branches) | `bulkDamage = yang`, 2D `FdemSolver.cpp:3333` + 3D `Fdem3dSolver.cpp:2197` ; défauts = Table 1 |
| base néo-hookéenne | **FAIT** | `bulkModel = neohookean` (composable) |
| compteur « pulverized elements » (leur fig. 18) | **FAIT** | `nPulv_` dans history.csv |
| dissipation d'endommagement ventilée | **FAIT** | `bdWork_` |
| friction résiduelle **de joint** (pic tanφ → résiduel par f(D)) | **FAIT** (2026-08-25) | `jointResidualMu`, 3D `Fdem3dSolver.cpp:2546-2560` |
| DIF de taux de déformation (leurs éq. 1-2) | **FAIT**, opt-in | `strainRateDIF = yang \| yang-fig2` |
| friction glissante 0,18 **au contact** post-rupture (perte d'appui sous l'insert, mobilité des fragments) | **MANQUANT** — c'est WP6 | — |

Decks existants : `impact_pulv_a/b.cfg` (réplique 4 corps + `bulkDamage = yang`),
antérieurs au correctif `jointShearRange = coulomb` (2026-08-27). **Aucun deck
ne combine encore pulvérisation et coulomb** ; les gardes sont composables
(`bulkDamage` n'exclut que les lois MatLaw, `Fdem3dSolver.cpp:272-276`).

## 2. Pourquoi le contact OUTIL est le point d'insertion décisif

Chaîne causale (revue adverse, vérifiée dans le code) :
1. il n'existe **jamais** de joint outil↔roche → `jointResidualMu` n'agit pas là ;
2. la voie potential saute toute paire dont le joint vit
   (`Fdem3dSolver.cpp:3018-3023` : « le joint vivant porte la paire ») ;
3. sous l'insert en compression avec `jointDeath = separation` (défaut), le
   joint ne meurt pas → le relais contact **n'existe pas** (piège documenté
   `:2652-2656`).

Donc l'unique chemin par lequel un μ dégradé pilote la **perte d'appui sous
l'insert** (le mécanisme du rebroussement non monotone du rebond) est le
contact outil : 3D `:3504/3519` (tanh pénalité) et `:3578` (cap Signorini) ;
2D `:4995/5004/5017` et `:5106`, lookup par `elemOf_[i]`. Les sites de contact
général (5 en 3D : 3147/3444/3504/3519/3578 ; 6 en 2D) portent l'autre moitié
(éjection des débris) — utile, mais secondaire pour le rebond.

## 3. Design retenu

- **Clé `contactResidualMu`** (défaut −1 = absent = bit-identique, même
  convention que `jointResidualMu`). Quand un nœud appartient à un élément
  **pulvérisé** (`bdD >= bdDmax`), le contact qui l'implique utilise ce μ au
  lieu de `contactMu`.
- **Bascule BINAIRE au franchissement de Dmax**, pas de rampe en D : le papier
  applique la sliding friction au matériau *rompu* (post-fracture), et la
  discontinuité est le geste déjà accepté au relais joint mort → contact.
  Une rampe dès D>0 contaminerait la force de pénétration en pleine charge
  (δ0 = 0,014 mm est franchi quasi immédiatement sous l'insert).
- **Gardes** : `contactResidualMu` sans source de D (`bulkDamage = off`) →
  throw ; `contactResidualMu > contactMu` → **warning bruyant** (pas throw,
  cohérent E3/E6) ; bannière `[FDEM(3D)] contactResidualMu = ...` au démarrage.
- **Sites CSV** : les compteurs/canaux vont dans les **4 sites par solveur**
  (en-têtes + lignes des scénarios tension ET percussion : 2D
  6215/6242/6293/6322 ; 3D 4051-4058 et 4070/4159) — la moitié oubliée casse
  silencieusement les repères zeroload.
- **Preuve bit-identique** : protocole du correctif coulomb — runs de
  référence 2D + 3D sans la clé, comparés au bit près avant/après.

À signaler sans le corriger en douce (constitution I) : asymétrie latente
`meanTensionCapFactor` — gardé `!law_` en 2D (`FdemSolver.cpp:3329`) mais pas
en 3D (`Fdem3dSolver.cpp:2194`). À trancher explicitement un jour.

## 4. Validation

- **Test 0 (analytique, cellule unique)** : traction uniaxiale pilotée ;
  cible fermée δm = 1,065e-4 m → D = 0,9 exactement (recalculée et confirmée
  en revue). Critère : exact.
- **Repère suite rapide** : jumeau zeroload (aucun contact → bit-identique
  clé posée/absente) + micro-impact 2 s vérifiant que μ passe à la valeur
  résiduelle sur éléments pulvérisés (ratio Ft/Fn à l'interface outil,
  nouvelle colonne history).
- **Banc `bench_pulverisation`** (3 leviers identifiés en revue, aucun des
  deux plans initiaux ne les avait tous) :
  1. variante `jointDeath = damage` pour isoler la contribution roche/roche ;
  2. `gcRestitution` FIXÉ dans les decks puis balayé (défaut 0,2 : il écrase
     le cap tangentiel à la détente et contamine e — personne ne le fixait) ;
  3. mesures : Ft/Fn interface outil, pénétration résiduelle max (raideur de
     pénalité non dégradée = interpénétrations possibles sur matière à 10 %
     de raideur), énergie cinétique des fragments éjectés.
  Critères PASS : perte d'appui (F sous l'insert à p égal, avec vs sans),
  non-monotonie de e sur 3 vitesses (leur fig. 9d), tendance nPulv (fig. 18).

## 5. Campagne de couplage (dès maintenant, sans code neuf)

`impact_pulv_coulomb.cfg` = `impact_pulv_a.cfg` + `jointShearRange = coulomb`
(+ `contactResidualMu = 0.18` quand WP6 est codé — pour St Anne, la valeur
calcaire de Solidity est 0,6 ; 0,18 est la valeur granite Table 1).
Prédiction inscrite : noyau broyé (coulomb) + perte d'appui (μ dégradé)
→ e non monotone en énergie ; sans WP6, coulomb seul ne suffit pas à faire
chuter e à haute énergie.

## 6. Limites assumées

- `bulkDamageDelta0/DeltaF` scalaires globaux : un deck multi-minéral
  (`phases`) garde des seuils uniques (contrairement à `crushCapP_[phase]`
  et au Weibull `ftScale`). Hétérogénéité des seuils = chantier séparé.
- Voie `law = dpdfh` non couverte (dégradation unilatérale en traction,
  DCAP codé en dur `MatLaw.cpp:732/903/1008` — un `dfhDmax` serait un autre
  chantier ; le fdem3d percussion n'en a pas besoin).
- Raideur de pénalité non dégradée : mesurée au banc (pénétration résiduelle),
  pas corrigée dans WP6.

## 7. Effort

Code ~0,5-1 j (11 sites + clé + gardes + colonnes) ; preuve bit-identique
0,5 j ; banc 0,5 j + runs. Le plus gros du chantier « pulvérisation » était
déjà fait (WP1-WP5) — WP6 est la dernière pièce mécanique.
