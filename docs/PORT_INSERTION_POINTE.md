# Port de la branche `insertion-pointe` dans `g0`

**Date** : 2026-09-05 · **Branche source** : `insertion-pointe` (2ead636, 4 commits) ·
**Cible** : `g0` · **Regle appliquee** : toute capacite arrive en cle **OPT-IN a
defaut bit-identique** (croissance par addition, principe VIII).

## 1. Ce que la branche apportait, et pourquoi c'etait vraiment orphelin

Les quatre commits de `insertion-pointe` n'ont jamais ete fusionnes et **aucune
de leurs six cles n'existait dans `f2`** (verifie en comparant les litteraux des
getters de configuration de `FdemSolver.cpp`, `Fdem3dSolver.cpp`, `MatLaw.cpp`
entre chaque branche et `g0` : `joint-handoff`, `dif-intrinseque` et `main`
rendent l'ensemble vide, `insertion-pointe` rend les six ci-dessous). Deux decks
de `f2` posaient deja ces cles et etaient donc **refuses** par `rockim_f2w21`
(garde C1 « la cle inconnue est une ERREUR ») :
`bench_impact/configs/impact3d_dpdfh.cfg` et `..._gros.cfg`.

| commit | apport | cles |
|---|---|---|
| `28c6ec2` | insertion preferentielle en POINTE de fissure (2D + 3D) | `insertionTipFactor`, `insertionTipDamage` |
| `f8b379a` | 3 decks de calibration du facteur (tunnel EDZ) | — |
| `c2c2d42` | dilatance VARIABLE psi(pbar) du DP-DFH + sortie de l'endommagement au `.vtu` | `dfhPsiVar`, `dfhPsi0`, `dfhKPsi`, `dfhPsiMax` |
| `2ead636` | `insertion = none` : le continuum pur (2D + 3D) | valeur nouvelle de `insertion` |

## 2. Ce qui a ete porte

### 2.1 `insertionTipFactor` / `insertionTipDamage` — l'enveloppe relachee en pointe

**Mesure qui la motive** (relevee dans la branche) : l'insertion adaptative ne
PROPAGE que 43,7 % de ses ruptures — le reste NUCLEE en terrain vierge — contre
58,9 % pour le schema intrinseque a loi de joint identique (chiffre CORRIGE le 2026-09-06 : « 56,8 % » etait une coquille de recopie ; la mesure d'origine, BILAN_insertion_adaptative.md l. 102/127/154, donne 58,9 %). Cause : la
contrainte est moyennee sur deux CST, ce qui ECRASE la singularite de pointe
(2,7 elements par zone cohesive de mode I) ; la facette devant une pointe ne se
distingue plus d'une facette quelconque de l'anneau plastique.

**Ce qui est fait** : une facette dont un sommet porte deja un joint insere et
endommage (`D >= insertionTipDamage`) voit son enveloppe **DIVISEE** par
`insertionTipFactor` ; une facette en terrain vierge garde l'enveloppe
nominale. Le choix de relacher la POINTE plutot que de penaliser la NUCLEATION
n'est pas cosmetique : penaliser la nucleation x1,6 faisait passer
`verify_fdem_voronoi_tension` de 15 joints rompus a ZERO — sans fissure
preexistante tout est terrain vierge, donc le facteur relevait la resistance
MACROSCOPIQUE de 60 % et aurait invalide le calage GBM Red Bohus.

- Defaut `insertionTipFactor = 1` : le facteur vaut 1 partout, aucun test
  supplementaire n'est evalue (`tipBias` est faux), **chemin d'origine au bit
  pres**. `insertionTipDamage` (defaut 0,5, valeur de la branche) n'est LU que
  lorsque le facteur est > 1.
- Gardes : `insertionTipFactor < 1` refuse (erreur nommee), `insertionTipDamage`
  hors [0,1] refuse, `insertionTipFactor > 1` sans `insertion = adaptive`
  refuse (en intrinseque tous les joints existent deja : il n'y a rien a
  inserer).
- Diagnostic : en fin de run, `insertions en POINTE : N propagations /
  M nucleations (x % de propagation)`, avec les deux references mesurees.
- Sites : `FdemSolver::init`, `FdemSolver::insertionSweep` (les DEUX branches,
  OpenMP et sequentielle), `FdemSolver::finalize` ; miroir exact en 3D
  (`Fdem3dSolver`), ou la pointe est un FRONT mais le test reste porte par les
  trois SOMMETS de la facette.

### 2.2 `insertion = none` — le continuum pur

Aucun joint n'existe ni ne peut naitre : les copies de noeuds restent liees pour
toujours (liaison rigide = elements finis a noeuds partages, exactement). C'est
le mode a employer quand la fissuration est portee par la **loi de volume**
(`law = dpdfh`, `saksala`, `bulkDamage`...) et non par des cohesifs.

**Pourquoi la cle existe.** On obtenait ce comportement en posant des
resistances de joint inatteignables (`ft = 1e12`). C'est une bombe a
retardement : mesure du 2026-08-25 sur l'impact 3D DP-DFH — un element
totalement endommage (D vers DCAP, 1 % de raideur, aucune suppression) se
distord, sort une contrainte aberrante, franchit **meme une enveloppe a 1e12**,
et les 89 424 joints qui s'activent alors portent `dnE = ft/pj = 8 cm` :
l'energie passe de 53 J a **-89 GJ en 10 microsecondes**. Avec
`insertion = none` le balayage n'existe pas, donc le piege non plus.

Trois consequences cablees, toutes portees :
1. les joints sont **lies** comme en adaptatif (`adaptive_ || noJoints_` avant
   `buildBindingTables()`) — sinon ils restent intrinseques et cassent, l'exact
   contraire du continuum voulu ;
2. les groupes lies **integrent comme UN noeud** (`integrate()`) — sinon les
   copies bougent independamment et le maillage se comporte comme un NUAGE de
   simplexes libres (bug du 2026-08-25 : outil freine de 0,1 J seulement,
   460 J d'energie elastique nee de rien) ;
3. le ressort de penalite des joints **sort du budget de pas de temps** : il
   n'est jamais evalue, il ne contraint donc pas la stabilite.

### 2.3 `dfhPsiVar` — la dilatance variable psi(pbar) du DP-DFH

Forme de `vumat_hole.f` / `vumat_kstdfh_psivar_phicap.f` l. 336-339 :

```
psi = clamp(dfhPsi0 - dfhKPsi * pbar[MPa], 0, dfhPsiMax)
```

avec les valeurs de la these (`dfhPsi0 = 160,345` deg, `dfhKPsi = 0,213793`
deg/MPa, `dfhPsiMax = 51,7` deg = `dfhBetaDeg`, donc ecoulement ASSOCIE tant que
le plafond mord). La constante #5 de la carte (`dfhPsiDeg`, 15 deg) est alors
MORTE. `dfhPsiVar` absente ou 0 : `psiDeg` fixe, chemin d'origine bit-identique.

### 2.4 Sortie de l'endommagement DP-DFH au `.vtu` (2D seulement)

`law = dpdfh` calcule trois endommagements directionnels dans un repere fige
(SDV 4-6 de la VUMAT) mais rien ne sortait : les `.vtu` ne portaient que les
contraintes. On ecrit desormais `dfhD = max(D1,D2,D3)` — exactement ce que
lisent les extracteurs du banc 6 (leur SDV 2) — et `dfhTini`, l'instant du
premier amorcage. **Ajout de sortie pur**, aucune trajectoire ne change ; en 2D
il prend la forme de deux entrees de la carte `vtk::ScalarField` de `f2`.

**RESERVE de la relecture adverse (2026-09-06).** La branche d'origine ecrivait
ces deux champs dans les DEUX solveurs (`git show insertion-pointe:src/FdemSolver.cpp`
l. 5635 et `:src/Fdem3dSolver.cpp` l. 3421) ; le port ne les a mis qu'en 2D —
`grep dfhD src/Fdem3dSolver.cpp` ne rend rien. Le titre de cette section disait
« 2D et 3D » : il est corrige. Reste a porter en miroir dans
`Fdem3dSolver::writeFrame()`, avec une passe `bitid` de confirmation.

### 2.5 Decks

- `tunnel_edz/configs/tunnel_tip13.cfg`, `tunnel_tip16.cfg`, `tunnel_tip20.cfg` :
  les trois decks de calibration du facteur (1,3 / 1,6 / 2,0) sur le tunnel EDZ,
  repris tels quels.
- `bench_impact/configs/impact3d_dpdfh.cfg` et `..._gros.cfg` : reparés (ils
  chargent maintenant, cf. §5).

## 3. Ce qui n'a PAS ete porte

- `build_dfh.cmd`, `build_dfh3.cmd`, `build_tip.cmd` : trois scripts de
  compilation `cl` ad hoc de la branche. `g0` nait avec CMake (decision 9) et
  `tools/build.ps1` ; ajouter trois `build_*.cmd` de plus serait du bruit.
- Le message d'annonce du bloc de liaison des joints a ete rendu HONNETE : il
  disait `adaptive insertion: N bonded edges` meme sous `insertion = none`. Il
  dit maintenant `insertion = none: N bonded edges` dans ce cas (sortie seule).

## 4. Preuves

### 4.1 Bit-identite globale

`python tools/bitid.py --exe build/rockim.exe --threads 4` contre l'ancre de
naissance `tools/bitid_refs.json` : **8/8 IDENTIQUE**, avant comme apres le
port. Rapport : `results/bitid_apres_port_insertion_pointe.json` (le nom
`bitid_g0-0.2.0_insertion_pointe.json` cite ici jusqu'au 2026-09-06 n'a jamais
existe, et il n'y a pas de tag `g0-0.2.0`).

**Attention** : les 8 decks de l'ancre couvrent `cdp`, `dpr`, `saksala2011`, les
joints cohesifs et le contact de Signorini — **aucun ne pose `law = dpdfh`**.
L'ancre ne prouve donc rien sur le portage de psi(p) ; d'ou le banc §4.2.

### 4.2 Banc court psi(p) — `tests_f2/psivar/`

`python tests_f2/psivar/check_psivar.py --exe build/rockim.exe --ref ../rockim_f2/rockim_f2w21.exe`
(deux points materiels DP-DFH, ~1 s chacun, aucun maillage) :

- **[A] bit-identite** : la trace du TEMOIN (`dfhPsiVar` absente) produite par
  `build/rockim.exe` est **rigoureusement celle** du binaire d'avant le portage
  `rockim_f2w21.exe` — meme SHA-256
  (`341fac848845623c2db4a702897f4efab6a797842d843964b6cf83274c997580`). C'est la
  preuve que l'ancre ne pouvait pas donner.
- **[B] la cle AGIT** : a `sigma3 = 0` l'ecart vaut `+1,04e-1`.
- **[C] psi decroit avec p** : `psi(p)` coupe les 15 deg fixes du temoin a
  `pbar = (160,345 - 15)/0,213793 = 679,8 MPa`. Le banc mesure le `pbar` au
  seuil (premier pas ou les deux traces divergent) et exige un **changement de
  signe** de l'ecart de deformation volumique de part et d'autre :

| sigma3 [MPa] | pbar au seuil | psi(p) | eps_vol TEMOIN | eps_vol ESSAI | ecart | attendu |
|---:|---:|---:|---:|---:|---:|---|
| 0 | 88 MPa | 51,7 deg | 1,3602e-02 | 1,1773e-01 | +1,04e-01 | psi > 15 : ecart > 0 |
| 200 | 434 MPa | 51,7 deg | 6,9094e-03 | 9,5055e-02 | +8,81e-02 | psi > 15 : ecart > 0 |
| 300 | 608 MPa | 30,5 deg | 3,5630e-03 | 2,2129e-02 | +1,86e-02 | psi > 15 : ecart > 0 |
| 400 | 781 MPa | 0,0 deg | 2,1664e-04 | -1,0977e-02 | -1,12e-02 | psi < 15 : ecart < 0 |
| 500 | 954 MPa | 0,0 deg | -3,1298e-03 | -1,3084e-02 | -9,95e-03 | psi < 15 : ecart < 0 |

  Les `pbar` mesures collent a la prediction analytique du seuil DP
  (`pbar = 1,7297 sigma3 + 88,4 MPa` pour `beta = 51,7` deg,
  `dcoh = 153,3 MPa`) : 88 / 434 / 607 / 780 / 954 MPa.
  **VERDICT : tous les criteres passes.**

- **[D] la variante qui DOIT echouer**, `--falsify` : le meme deck d'essai avec
  `dfhPsiVar = 0` rend un ecart max de `+0,000e+00` et une trace **bit-identique**
  au temoin, donc le critere B echoue comme il doit. L'option etait annoncee dans
  l'en-tete du script et dans `mp_dpdfh_psivar_on.cfg` mais **n'existait pas** ;
  elle a ete ecrite a la relecture adverse du 2026-09-06.

### 4.3 Banc court insertion — `tests_f2/insertion_pointe/`

`python tests_f2/insertion_pointe/check_tip.py --exe build/rockim.exe --threads 2`
(quatre compressions uniaxiales 2D courtes sur maillage voronoi, 2 399
elements, arret automatique apres le pic) :

| deck | dt [s] | joints inseres | rompus | pic axial [MPa] | sha256(history.csv) |
|---|---:|---:|---:|---:|---|
| `ins_tip_ref` (cle absente) | 1,46441e-08 | 1290 | 334 | 51,0807 | `d2bbb36ea38acc57…` |
| `ins_tip_1` (`= 1` pose) | 1,46441e-08 | 1290 | 334 | 51,0807 | `d2bbb36ea38acc57…` |
| `ins_tip_16` (`= 1,6`) | 1,46441e-08 | 1370 | 239 | 51,0807 | `3e9ed03de36f1065…` |
| `ins_none` | 2,51169e-08 | — | — | 84,7594 (pas de rupture) | `03fff7971e34e783…` |

- **[A] defaut NEUTRE** : cle absente et `insertionTipFactor = 1` rendent un
  `history.csv` **bit-identique**.
- **[B] la cle AGIT** : le facteur 1,6 change la trace et imprime
  `305 propagations / 1065 nucleations`.
- **[C] l'amorcage est INTACT** : le pic passe de **51,0807 a 51,0807 MPa**,
  a l'affichage pres. C'est exactement la propriete que le choix « relacher la
  pointe » visait, et que « penaliser la nucleation » n'avait pas (+60 % sur la
  resistance macroscopique). La variante qui doit echouer est donc la
  penalisation de la nucleation : elle deplacerait ce pic.
- **[D] continuum pur** : `insertion = none` annonce le mode, n'insere aucun
  joint, et tourne a **dt x1,715** (le ressort des joints est sorti du budget).

## 5. Les deux decks `impact3d_dpdfh*` reparés

Ils posaient les six cles et etaient refuses par `w21`. Ils chargent maintenant
sous `build/rockim.exe` (garde C1 : `140 cles consommees, 1 du deck non lue`
— `verifyFt`, lecture conditionnelle legitime). Deux autres reparations ont ete
necessaires :

1. **le maillage n'existait nulle part.** `meshes/impact3d_dfh.msh` et
   `meshes/impact3d_gros.msh` ne sont ni dans `g0`, ni dans `rockim_f2`, et les
   `.msh` ne sont pas versionnes. La **commande exacte de regeneration** est
   desormais dans l'en-tete de chaque deck. Elle ne reproduit pas au tetraedre
   pres les comptes annonces par la version d'origine (98 858 et 143 451), dont
   les parametres n'ont jamais ete consignes : le deck le dit et donne les
   comptes reellement obtenus (105 498 et 150 535 tetraedres).
2. **le noeud orphelin.** gmsh ecrit le point du champ de taille (le point
   d'impact lui-meme) comme un noeud 0-D jamais reference par un tetraedre ; le
   lecteur de maillage le refuse (broche fantome de masse nulle sous le pole,
   cf. `ETAT_DES_LIEUX_2026-09-05.md`). Les outils qui le nettoient,
   `etude_lois_fem/meshes/drop_orphans.py` et `check_orphans.py`, n'avaient pas
   ete repris a la naissance de `g0` (seuls les `*.py` de la RACINE de
   `etude_lois_fem/` l'avaient ete) alors que le message d'erreur du solveur les
   NOMME : ils sont ajoutes, et les decks pointent sur `*_clean.msh`.

Les decks n'ont **pas** ete lances (T = 1,5e-4 s sur 105 k / 150 k tetraedres :
c'est un run de campagne, pas un banc). Seule leur INITIALISATION a ete jouee,
sur une copie a `T = 3e-7` s.
