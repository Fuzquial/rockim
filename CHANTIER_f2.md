# Chantier f2 — cinq capacités, et un crash préexistant corrigé

**2026-09-01. F. Uzquiano, Mines Paris – PSL.**
Cinq capacités ajoutées au solveur FDEM 2D pour AbuAisha et al. 2017
(JPSE 154, 100-113) : `preBrokenJoints` (leur §2.4), `microseismic`
(leurs éq. 11-13), `dampingViscous` (leur éq. 9) et `hydroCavityClosure`
(leur éq. 4), et le **couplage thermo-mécanique complet** — conduction transitoire sur le graphe des joints, Robin, contrainte thermique dans l'élément, phasage `toolStart` et bascule d'amortissement (voir DOCUMENTATION §5.12 : conservation exacte à 1,4e-13, forme fermée E αT ΔT/(1−ν) reproduite à 0,011 %). Toutes **opt-in, inertes par défaut**.

> ### ⚠️ Un crash préexistant, trouvé en chemin et corrigé
>
> **`hydro = on` avec `insertion = intrinsic` corrompait le tas.**
> `nVert_` n'était écrit que par `buildBindingTables()`, qui ne tourne qu'en
> mode adaptatif ; en intrinsèque il restait à 0, `updateWetBoundary()`
> allouait `wetV` à **un** élément puis y écrivait pour tous les sommets du
> maillage. Reproduit sur trois combinaisons (voronoï/grille × débit/pression),
> `0xC0000374` ou `0xC0000005`, **avant même le message `[HYDRO]`** — et
> reproduit à l'identique sur le binaire **de référence**, donc antérieur à ce
> chantier.
>
> **Ce n'est pas académique** : le deck `e5_intrinsic_iso12.cfg`, préparé pour
> trancher le dernier écart structurel avec Y-Geo, pose précisément
> `insertion = intrinsic` **et** `hydro = on`. Il aurait planté au démarrage.
>
> Correctif : poser `nVert_` dès la construction du maillage, où `vOf_` vient
> d'être rempli. `buildBindingTables()` le recalcule à l'identique, donc aucun
> chemin existant ne change — **vérifié bit-identique, 28 fichiers**.

---

## 0. Ce qu'est ce dossier

Copie de `rockim_f1` (mini-dépôt de partage du 31 août), **moins** les
2 Go de sorties de campagne de `calib_triax3d/`. C'est `f1` et non `rockim_p1`
qui sert de base : son `FdemSolver.cpp` fait 7357 lignes contre 6319, il est
donc nettement en avance.

> ⚠️ **Conséquence à retenir.** L'audit du 2026-09-01 portait sur `rockim_p1`
> et concluait qu'aucune pré-rupture n'existait. **C'est faux sur `f1`** : la
> machinerie `applyPrebrokenPopulation()` / clé `jointPrebrokenFrac` /
> champ `Joint::pre` y est datée du 31 août. Ce chantier n'a donc **rien créé
> de neuf** côté pré-rupture : il a ajouté un **sélecteur géométrique** à une
> mécanique déjà écrite, déjà gardée, et qui avait déjà résolu le piège
> principal (la reconstruction des tables de liaison).

Build : `.\build_f2.cmd [nom.exe]` — MSVC + les en-têtes Eigen de
`..\rockim\eigen-3.4.0`. Pas de CMake sur cette machine.

---

## 1. Les deux capacités

Référence complète des clés : `DOCUMENTATION_rockim.md` §5.11.

### `preBrokenJoints` — les discontinuités préexistantes de leur §2.4

```
preBrokenJoints   = 0.015 0.06 0.045 0.06 ; 0.02 0.09 0.05 0.09
preBrokenTol      = -1        # < 0 = demi-arête locale (défaut)
preBrokenAngleDeg = 30        # écart d'orientation toléré
jointResidualMu   = 0.5       # OBLIGATOIRE
gcSurfaceRefresh  = eager     # fortement recommandé
```

Leur joint préexistant n'est **pas** un élément cohésif affaibli : c'est un
élément de fracture **déjà rompu** — géométrie plane, ouverture
infinitésimale, aucune cohésion, cisaillement par le frottement résiduel de
leur éq. 2. C'est exactement l'état que pose la fonction.

Le sélecteur ajoute **deux** conditions à la sélection : distance du **milieu
d'arête** au segment, et **écart d'orientation**. La seconde n'est pas
cosmétique — sans elle on ramasse les arêtes transverses qui croisent le
segment, et la « discontinuité » devient un escalier de joints
perpendiculaires qui scie la roche au lieu de la fendre.

### `microseismic` — leurs éq. 11-13

```
microseismic = true
```

Instrumente chaque joint (`tYield`, `keYield`, `dKeMax`) et écrit
`fdem_seismic.csv` :
`jointId, x, y, tYield, keYield, dKeMax, tBreak, breakMode, magnitude`.

**Filtrer sur `tBreak >= 0`** pour retrouver leur définition de l'événement
(la rupture) ; les joints endommagés non rompus sont fournis parce qu'ils ne
coûtent rien et documentent l'incubation.

---

## 2. Ce qui est prouvé, et comment

### Bit-identité sans les clés — **preuve forte, faite**

Même jeu de données (`configs/verify_fdem_voronoi_tension.cfg`,
`OMP_NUM_THREADS = 1`) joué par le binaire **d'avant** et celui **d'après**,
puis comparaison **SHA-256 fichier par fichier** de tout l'arbre de sortie :

> **28 fichiers, 0 différence d'octet.** Et pas de `fdem_seismic.csv` quand la
> clé n'est pas posée.

C'est la preuve de niveau « diff d'arbre », plus forte que « la suite passe » :
elle porte sur `history.csv`, toutes les trames VTU et les CSV de fin.

### Le chemin `jointPrebrokenFrac` existant est intact — **preuve faite**

C'est l'édition la plus risquée du chantier : pour que les messages de garde
nomment la clé réellement armée, ils sont passés d'un littéral à une
concaténation avec une variable. Contrôle dédié (`tests_f2/t5_fracpath.cfg`,
`jointPrebrokenFrac = 0.02`) :

> **28 fichiers bit-identiques**, et **journaux identiques ligne pour ligne**
> (26 lignes, hors temps machine) entre le binaire d'avant et celui d'après.

La concaténation rend donc les textes identiques **au caractère près** quand
seule `jointPrebrokenFrac` est posée, et le bloc de journal propre au
sélecteur géométrique ne s'imprime pas.

### Suite de non-régression `fast` — **binaire de référence : 44/44**

Passée sur `rockim_base.exe` (MSVC, `OMP_NUM_THREADS = 1`) : **TOUT PASSE
(44/44)** contre les références **Linux embarquées**, sans avoir besoin d'une
baseline de plateforme. C'est un résultat en soi : le dossier prévenait que
macOS ARM en perdait 10 sur 44 et recommandait de créer une baseline avant
toute modification. **MSVC n'en perd aucun** — la baseline n'était pas
nécessaire ici.

**Et sur le binaire modifié : 44/44 également.** Zéro échec, contre les mêmes
références Linux embarquées. Passée sur `rockim_f2c.exe` ; le binaire final
`rockim_f2d.exe` n'en diffère que par l'option de schéma implicite, dont
l'inertie a été prouvée bit-identique sur 28 fichiers, chemin explicite
compris — le verdict s'y transporte donc.

| binaire | contenu | suite `fast` |
|---|---|---|
| `rockim_base.exe` | copie de `f1`, aucune modification | **44/44** |
| `rockim_f2c.exe` | + les 4 capacités + le correctif de crash | **44/44** |
| `rockim_f2d.exe` | + l'option de schéma implicite | bit-identique à `f2c` |

### Les capacités font ce qu'elles annoncent — **fait**

| essai (`tests_f2/`) | attendu | mesuré |
|---|---|---|
| `t1_seismic` | l'instrumentation ne perturbe **pas** la mécanique | pic **11,156 MPa**, **15** ruptures — identiques au témoin ; 133 événements écrits |
| `t2_prebroken` | une entaille de 30 mm dans une éprouvette de 60 mm doit faire **chuter** le pic | 8 joints sélectionnés, pic **4,86 MPa** contre 11,156 |
| `t3_garde_vide` | segment hors du maillage → **refus** | refus, code 1, message nommant les trois clés à vérifier |
| `t4_garde_mu` | `preBrokenJoints` sans `jointResidualMu` → **refus** | refus, code 1, message nommant `preBrokenJoints` |
| `t5_fracpath` | le chemin `jointPrebrokenFrac` **existant** reste intact | 28 fichiers bit-identiques, journaux identiques ligne pour ligne |
| `t6_cavclose` | volume : correction **nulle** à contour fermé, non nulle dès qu'une bouche s'ouvre | écart initial **0,000**, puis $9{,}4\cdot10^{-7}$ m²/m ; mécanique **rigoureusement identique** (même $t$, même `nBroken` sur 2080 lignes) |
| `t7_damping` | l'amortissement nodal freine, et la garde de pas de temps mord | garde $2m/\mu e$ prend la main (pas ÷ 6,5) ; pic 10,8 → **26,0 MPa** |
| `t8_damping_impl` | le schéma implicite : même physique, sans la borne de pas | pas de temps **rendu à sa valeur sans amortissement**, 283 415 → **43 652 pas** ; pic **identique à $10^{-4}$ près** |
| `t8_garde_orpheline` | `dampingViscousScheme` sans `dampingViscous` → **refus** | refus, code 1 |
| `t9`/`t9b_thermo` | forme fermée $\sigma_{yy} = E\alpha_T\Delta T/(1-\nu)$, conservation | **0,011 %** en adaptatif (−4,8 % en intrinsèque = complaisance de pénalité, mesurée par une voie indépendante) ; conservation **exacte** 1,4·10⁻¹³ |
| `t10_thermo_cracks` | $\Delta T = 40$ K → $\sigma$ = 21 MPa ≫ $f_t$ | **21 joints rompus** ; à 10 K : **0** |
| `tg_orphan`/`tg_law`/`tg_3d` | gardes thermiques | trois refus propres, code 1 |

L'essai `t6` mérite un mot sur sa construction : il arme l'hydro mais impose
une **pression nulle**, donc aucune force hydraulique. La mécanique est alors
celle du témoin à la ligne près, et le volume rapporté est la **seule** chose
qui peut différer — c'est ce qui rend la comparaison concluante.

L'essai `t8` conforte le précédent d'une façon qui n'était pas acquise : les
schémas explicite et implicite donnent **le même pic à $10^{-4}$ près**. La
montée de 10,8 à 26 MPa n'est donc pas un artefact du traitement explicite,
c'est le modèle. Et le schéma implicite rend au pas de temps sa valeur d'un run
sans amortissement — 43 652 pas au lieu de 283 415, soit **6,5 fois moins
cher pour la même réponse**. C'est le schéma à préférer dès qu'on explore des
$\mu$ élevés.

Et `t7` livre un résultat qui dépasse la vérification de code : voir
« l'amortissement : deux lectures » dans `DOCUMENTATION_rockim.md` §5.11.
À la valeur de leur Table 1, la constante de temps nodale $\tau = m/(\mu e)$
vaut $1{,}7\cdot10^{-8}$ s sur une maille de 3 mm, **contre un pas de temps de
$2{,}1\cdot10^{-8}$ s** : un nœud ne serait pas amorti, il serait figé en un
pas. La lecture littérale de leur éq. 9 n'est donc pas utilisable **à leur
valeur** — argument de plus pour lire leur $\mu$ comme la viscosité de volume
de Munjiza, que `bulkViscosity` porte déjà.

Le pic de `t2` mérite un mot : le ligament passe de 60 à 30 mm, donc la
section nette doublerait la contrainte et prédirait ~5,6 MPa ; on mesure 4,86,
un peu plus bas — l'écart est la concentration en pointe d'entaille. La
direction **et** l'ordre de grandeur sont ceux qu'on attend.

### Cohérence du catalogue — **fait**

Sur `t1` : 133 joints entrés en endommagement, dont **15 rompus** — le compte
exact que le solveur affiche par ailleurs. $t_y \le t_f$ sur **tous** (0
violation). Magnitudes des ruptures : **−5,11 à −4,66**.

> Leur figure 17 porte des magnitudes **de l'ordre de −4 à −5**. Le nôtre
> tombe dans la même bande sur un problème pourtant très différent (éprouvette
> de traction, pas forage). Ce n'est **pas** une validation de la physique —
> c'est une vérification que l'éq. 13 et les unités sont bien implantées.

---

## 3. Ce qui n'est PAS prouvé, et qu'il faut vérifier au premier run du banc

**Le mouillage d'une pré-fissure par une fissure hydraulique qui l'intersecte
est établi par LECTURE DU CODE, pas par un run.** Le raisonnement est solide —
le front mouillé ne teste que `!bonded && D >= wetDmin_`, rien d'autre, et
l'état posé est identique à celui de `jointPrebrokenFrac` qui, lui, tourne
depuis le 31 août — mais aucun calcul hydraulique ne l'a exercé ici : les
maillages du banc (`hf_bore.msh`, `parker_crack.msh`) vivent dans
`rockim_p1/meshes` et n'ont pas été copiés.

**À contrôler explicitement au premier run de B5**, dans cet ordre :

1. `hydroNWet` doit **rester à 105** (les faces du forage) tant que la fissure
   n'a pas atteint la discontinuité — si la pré-fissure était happée par
   `boreSelectR`, elle serait pressurisée dès $t = 0$ et l'essai serait faux ;
2. au moment où la fissure atteint la discontinuité, `hydroNWet` doit sauter
   **d'un coup** de deux fois le nombre de joints de la chaîne — c'est la
   « perméabilité infinie » de leur §2.4, la relaxation étant un point fixe
   global et non incrémentale ;
3. `hydroVol` doit croître en conséquence.

**Second point non prouvé** : `preBrokenJoints` n'a jamais été exercé sur un
maillage **gradué**. Le défaut `preBrokenTol < 0` (demi-arête locale) a été
conçu pour cela, mais il n'a tourné que sur un Voronoï à grain quasi uniforme.
Lire le compte `nGeo` du journal et le comparer à la longueur attendue de la
chaîne.

---

## 4. Les pièges rencontrés, pour qui reprendra

**Le piège principal, résolu par `f1` avant nous.** En insertion adaptative,
les copies de nœuds d'un joint `bonded` sont **liées cinématiquement** et
`grpsOfVert_` est un **cache**. Poser `bonded = false` après la construction
des tables ne délie rien : le joint tourne, mais ses quatre copies étant
soudées, $dn \equiv 0$ et il ne transmet ni n'oppose rien. **La défaillance
est silencieuse** : la pré-fissure devient indiscernable d'une arête intacte
incassable. D'où l'appel à `buildBindingTables()` après la pose, déjà présent.

**Corollaire à ne pas oublier** : la pré-rupture doit être posée **après** la
ligne qui force `bonded = true` sur tous les joints en mode adaptatif — sinon
elle est effacée sans un mot, et **seulement en mode adaptatif**, si bien
qu'un essai en intrinsèque passerait et masquerait le bug.

**Une chaîne isolée est cinématiquement inerte.** Une pré-fissure dont aucune
extrémité n'est scindée a ses deux paires de copies dans le même groupe : elle
ne peut ni s'ouvrir ni glisser. C'est `nFree`, pas `nPre`, qui mesure la
population active à $t = 0$ — d'où l'avertissement ajouté.

**Le trou du cache de contact.** Les faces mortes du contact général sont
estampillées sur `nBroken_`, que les pré-fissures **n'incrémentent jamais** (et
c'est voulu : une fissure préexistante est une condition initiale, pas une
réponse au chargement). La première pré-fissure qui meurt par séparation, s'il
n'y a eu aucune vraie casse avant, resterait donc invisible du contact
général. `gcSurfaceRefresh = eager` supprime ce trou — le code avertit, et
pour un banc à discontinuités préexistantes il faut le poser.

**Ne PAS « corriger » cela** en faisant compter les pré-fissures dans
`nBroken_` : cela fausserait toutes les courbes de fissuration et casserait la
bit-identité des runs sans la clé.

---

## 5. Ce qui reste à faire

- Exercer les deux capacités sur un **cas hydraulique** (§3 ci-dessus).
- Outil de maillage `make_bore_joint_mesh.py` (forage **et** joint dans le
  même maillage, position et orientation libres) pour leur figure 8, et
  `make_joint_swarm.py` pour le semis de leur figure 14 — `preBrokenJoints`
  accepte déjà plusieurs segments, donc un semis se déclare par une chaîne.
- Post-traitement `micro_seismic.py` : lire `fdem_seismic.csv`, filtrer sur
  `tBreak >= 0`, tracer le semis et l'histogramme des magnitudes au format de
  leurs figures 16-17.
- Décider si le catalogue doit aussi porter les **événements de glissement**
  sans rupture (leur seconde famille) : rien ne les enregistre aujourd'hui, et
  c'est un chantier distinct des éq. 11-13.

---

## 7. Les trois briques de Lisjak — schistosite pervasive (2026-09-01/02)

La revue de litterature (tunnel_schisto/BIBLIO_anisotropie_fdem.md, 4 references
verifiees) a etabli que l approche a plans discrets `weakPlanes` (section 5.13
de la doc) est celle que son auteur a abandonnee a l echelle de l ouvrage. La
methode de reference — Lisjak, these Toronto 2013, ch. 5 — tient en trois
briques ; DOCUMENTATION_rockim.md section 5.14 en donne la reference complete.

| brique | ou | etat |
|---|---|---|
| 1. elasticite transversalement isotrope | triangle (`setupBeddingElastic`) | faite, validee en forme fermee |
| 2. loi cohesive directionnelle | joint (`applyBeddingCohesive`) | faite, validee (bit-identite neutre, anisotropie 1,94) |
| 3. maillage a aretes alignees sur le litage | `tunnel_schisto/make_tunnel_bedded_mesh.py` | fait, continuite 1,015-1,018 sur les trois pendages |

### Ce qui est prouve

| essai (`tests_f2/`) | attendu | mesure |
|---|---|---|
| bit-identite f2g -> f2h -> f2i, cles absentes | 0 octet | 28 fichiers, 0 different (x2) |
| `t15_gamma_neutre` (ratios = 1) | bit-identique | 28 fichiers, 0 different |
| `t13c_nocund_*` TI a la limite isotrope, elastique, SANS Cundall | ecart nul | ecart EXACTEMENT nul sur 2080 lignes |
| `t13b` idem AVEC Cundall | — | 7,5e-4 : le signe de la vitesse amplifie des arrondis a 1e-16 (pas un bug) |
| `t14c_ti_*` modules apparents, mesure tout interieur, adaptatif | 1,5910 / 4,3307 GPa | 1,5910 (+0,00 %) / 4,3364 (+0,13 %) |
| `t14b_ti_*` idem intrinseque | — | -2,0 % / -4,8 % = complaisance de penalite, en serie, independante de la direction |
| `t16_gamma_0/90` ratios Lisjak | anisotropie | pic 4,35 -> 8,43 MPa (rapport 1,94 ; Lisjak T_P/T_S ~ 1,9) |
| maillage lite, fumee (h 0,6 / t 1,8) | cordes plongees | 1619 aretes exactement alignees, continuite 1,007 |
| maillages lites reels (h 0,12 / t 0,35, 0/45/90 deg) | continuite ~1 | 303 580 / 303 907 / 304 052 triangles, continuite 1,015 / 1,016 / 1,018 |
| dry-run `tunnel_lisjak45` sur le maillage lite | briques 1+2 armees, pas de WARNING | TI armee (rapport apparent 2,56), 455 595 joints, 19,3 % quasi paralleles — le WARNING de la brique 2 s eteint : la brique 3 est validee PAR la 2 |
| gardes : TI+thermal, TI+law, TI+phases, beddingDip seule, 3D | refus | refus propres |
| suite `fast` sur `rockim_f2i.exe` | 44/44 | **44/44** |

### Deux pieges trouves en chemin, a ne pas reproduire

1. **Coquille de la these de Lisjak (p. 102)** : « maximum and minimum values,
   for gamma = 0 and 90 » est l INVERSE de sa fig. 5.8a et de sa table 5.1.
   Verifie sur le PDF. Le minimum est a gamma = 0 (joint parallele au litage).
2. **Gmsh/OCC** : `occ.fragment` DECOUPE les cordes a la cavite mais ne les
   ATTACHE PAS a la face (0 frontiere, tous les troncons sans face) ; il faut
   `mesh.embed` des troncons dont le milieu est dans la roche. Et la metrique
   « aretes a moins de 10 deg » ne bouge pas (la triangulation contrainte
   appauvrit le voisinage de la ligne dans sa direction) : compter les aretes
   a moins de 1 deg et rapporter leur longueur aux cordes (continuite).

### Ce qui n est PAS fait

- thermal + TI : refuse (la contrainte thermique est isotrope, -3 K alphaT dT I ;
  sous TI elle devrait devenir -D:(alphaT dT 1), non derive).
- Lysmer/absorbing sous TI : lit encore Material::cP() isotrope (approximation,
  decks tunnel en absorbing = none).
- alphaT par phase, et la loi gamma par phase : hors perimetre.
- Aucun run de production : tunnel_schisto/tunnel_lisjak00/45/90.cfg attendent
  validation. Cout mesure au dry-run : dt = 1,06e-6 s, 662 000 pas,
  ~35 h par run, ~106 h pour les trois, sequentiels.

| binaire | ajout | suite fast |
|---|---|---|
| `rockim_f2g.exe` | weakPlanes + continuite | 44/44 |
| `rockim_f2h.exe` | briques 1 et 2 | bit-identique a f2g |
| **`rockim_f2i.exe`** | + verrou thermal/TI | **44/44** |

### 7.1 Le pas de temps, les slivers, et le reglage a 4 h 30 (2026-09-02)

Question posee : « comment on est passe de 4 h a 35 h ? juste en ajoutant de
la schistosite ? » Reponse mesuree : NON. Les briques 1 et 2 ne coutent rien
(la borne CFL — 1,57 us sous TI — ne mord pas ; c est le budget de raideur
nodale qui fixe dt). Le x8,3 venait du MAILLAGE : x2,86 elements (h = t/3 a
t = 0,35) et x2,89 pas, dont ~x1,6 dus a 17 slivers sur 304 000 triangles.

Mecanisme, verifie a 10 % pres par la formule dt = 0,4 sqrt(m_i/(K_i+2 k_c)) :
en `mesh = file`, `pj = 4E/hmin` avec hmin GLOBAL (src/FdemSolver.cpp:2295),
donc un seul sliver gonfle la penalite de TOUS les joints (x3,6) ET reduit la
masse du pire noeud. Les ressorts de joints font 73-85 % de K_i. En adaptatif
la boucle de budget est inconditionnelle (les joints lies comptent, a pj =
4E/h au lieu de 20E/h en intrinseque — c est le « dt x2 » de l adaptatif).

Six versions de l outil de maillage pour tuer les slivers :
  v1 fragment disque+cordes -> extremites sur le cercle sans face (isotrope)
  v2 fragment cordes -> decoupe sans attacher ; mesh.embed regle ca
  v3 idem : sommets imposes sur la paroi -> 17 slivers (d_min 20 mm)
  v4 decoupe python + retrait fixe 0,5 h -> 109 slivers DANS le vide du retrait
  v5 retrait h/sin(theta) -> 6 slivers restants la ou une corde LONGE la paroi
  v6 bande d exclusion 0,8 h le long de la paroi -> d_min 76 mm > isotrope 73
Diagnostic decisif a chaque etape : position des triangles sous 40 mm par
rapport a la paroi (tous a 3-10 cm), et le fait que la variante « retrait plus
large » donnait un maillage strictement identique (donc pas les extremites).

Resultat : reglage PRODUCTION t = 0,60 / h = 0,20 = meme densite que le run
isotrope de reference. Trois maillages 0/45/90 : 117 200 / 117 132 / 117 174
triangles, d_min 75,5 / 76,2 / 76,6 mm, continuite 1,005 / 1,007 / 1,010.
Dry-run : dt = 3,22e-6 s (> 3,05 isotrope), 217 730 pas -> ~4 h 30 par
pendage, 13 h 30 les trois. Decks : tunnel_schisto/tunnel_lisjak00|45|90_4h.cfg
(diff hors commentaires : meshFile + beddingDip, rien d autre).
Reglage CONVERGENCE t = 0,35 / h = 0,12 conserve (tunnel_lisjak45.cfg,
~17 h) pour verifier que h/v ne depend pas de t.

Piste solveur NON ecrite (demande accord) : `jointPenaltyLocal` pour
mesh = file, h local par joint comme le fait deja le mode Voronoi —
decouplerait dt du pire element global. Non necessaire au reglage a 4 h 30.

### 7.2 Premier run de production : tunnel lite 45 deg, methode de Lisjak (2026-09-02)

Deck tunnel_schisto/tunnel_lisjak45_4h.cfg, rockim_f2i.exe, 14 threads,
lance 01:23, termine 05:34 : 15 053 s = 4 h 11 de mur (estimation 4 h 30 ;
l extrapolation a 35 min faite a 3 % du run etait fausse — la phase
d equilibrage ne coute rien, la fissuration coute tout).

Resultats (sortie out_lisjak45_4h, temoin isotrope out_tun_corr_th4 de p1) :
| mesure                                   | lite 45      | isotrope     |
|------------------------------------------|--------------|--------------|
| joints rompus                            | 18 339       | 26 524 (-31 %) |
| part cisaillement (solveur)              | 64 %         | 58 %         |
| EDZ rayon p95 (edz_metrics)              | 10,75 m      | 20,33 m      |
| demi-axes h / v (edz_metrics)            | 9,54 / 9,41  | 19,39 / 14,36 |
| profondeur p95 - paroi LE LONG du litage | 6,0 m        | 14,5 m (memes secteurs) |
| profondeur EN TRAVERS du litage          | 3,4 m        | 12,2 m       |
| rapport le long / en travers             | **1,77**     | 1,19 (bruit h/v 1,35) |
| deplacement max (reins)                  | 3,42 m       | 0,99 m (Wang 0,347) |
| blocs mono-element (block_sizes)         | 93 %         | 77 %         |
| 5 plus gros blocs [m2]                   | 2,8 1,4 0,8 0,4 | 55,8 19,1 16,4 11,4 |
| selection MC conjuguee +-[18,42] deg     | 1,14 (hors 0,81) | 1,65 (hors 0,94) |
| selection par orientation absolue        | ~1 partout   | ~1 partout   |

Lectures :
1. H1 (le grand axe tourne avec le pendage) : OUI sur la profondeur
   d enveloppe par secteur (1,77 contre 1,19), et visible a l oeil (losange a
   45 deg, halos sigma_1/sigma_3 de meme forme). ATTENTION : edz_metrics
   mesure des demi-axes HORIZONTAL/VERTICAL, aveugles a une ellipse a 45 deg
   (h/v = 1,01 ici). Pour les pendages 0 et 90 l outil convient ; pour 45 il
   faut la profondeur par secteur (script inline du 2026-09-02, a ranger).
2. Orientation des fissures : PAS de selection (rompu/offre ~ 1 dans tous les
   secteurs, comme l isotrope). L EDZ est SATUREE — tout casse — donc le test
   « les fissures piquent a 45 deg » ne discrimine pas ici. La selection
   conjuguee de Mohr-Coulomb est affaiblie (1,14 contre 1,65).
3. La roche litee casse MOINS LOIN (EDZ deux fois moins profonde, -31 % de
   joints) mais PLUS FIN (93 % de blocs mono-element, aucun bloc > 3 m2
   contre 56 m2 en isotrope) et CONVERGE ENORMEMENT (3,4 m aux reins : la
   peau comminuee s ecoule dans le vide — effondrement des piedroits).
4. Reserves : (a) un seul run, une graine ; (b) loi ou chemin de maillage
   non separes (temoin ratios = 1 non lance) ; (c) les RAPPORTS d Opalinus
   (c/9, GIc/17) appliques a une matrice deja tres faible donnent
   c_litage = 0,09 MPa — la severite de ces rapports est un choix a
   discuter, une sensibilite est necessaire ; (d) Wang ne modelise pas de
   litage : le 3,4 m contre 0,35 m n est pas une comparaison a l article.

### 7.3 Solutions a la remarque « rien n a traverse le litage » (2026-09-02)

Revue : tunnel_schisto/REVUE_traversee_litage.md (13 agents, 8 references
confirmees par relecture, ~45 lues). Verdict : avec Gamma_i/Gamma_b = 0,057 et
sigma_b/sigma_i = 4,07, l arret sur plan est PREDIT par He-Hutchinson (seuil
1/4 a 90 deg, 0,55-0,9 en oblique) et Parmigiani-Thouless ; les rapports de
tenacite MESURES dans les shales valent 0,10-0,30 en G, jamais 0,057 (constante
de calibration de Lisjak, dependante de son maillage). In situ, la fissure
isolee s arrete au litage (Nussbaum 2011, Chandler 2016), mais la ZONE traverse
par trois voies absentes du run : anisotropie de contrainte (Bure, lambda 1,3
suffit), front 3D, flambage des plaques avec le temps.

Applique EN AJOUT (rockim_f2j.exe, bit-identique cles absentes 28/0) :
- S8  writeJointState : sigN, tauS, dn, contactState par joint (diagnostic de
      Renshaw-Pollard) + tools/joint_state_stats.py (etats des plans-frontieres,
      evenements de traversee). Valide t17.
- S5  weakPlaneFactor2 / weakPlaneFrac2 / weakPlaneSeed : lits faibles et lits
      forts tires par plan. Valide t18 (13 joints sur lits forts / 80).
- S6  calcul : l_cz mode I sur le litage = 179 mm pour une arete de 200 mm ->
      NON RESOLUE ; en travers 556 mm (resolue) ; mode II resolu. Le 17,5 agit
      en resistance x ouverture critique, pas en energie. Remede : maillage fin.
- S1/S2/S3/S4/S7/S9a : 13 decks tunnel_schisto/S*.cfg, chacun a variable
  unique par rapport a tunnel_lisjak45_4h.cfg (+ writeJointState), dry-runs
  passes. S9a exige jointFrictionScaled = 0 (exclusif de jointResidualMu).
- tools/edz_sectors.py (edz_metrics est aveugle a 45 deg), tools/tip_velocity.py.
- SUITE_solutions.md : ordre S4a, S4b, S8, S1(0,30), S3(1,3) ; run_solutions.ps1.
Rien de lance. Non couvert : 3D (S11), temps/pressions interstitielles (S10),
tenacites de l argilite cible (aucune mesure par orientation trouvee).

| `rockim_f2j.exe` | + writeJointState (S8) + lits bimodaux (S5) | **44/44** (2026-09-02 09:40) |

### 7.4 Calibration rapide Red Bohus et polydispersite des grains (2026-09-02)

Dossier `calib_quick/` (README a jour). Trois cas de moins de 5 min a
sigma3 = 50 MPa, bulk elastique, joints de la sonde 4 : homogene 709 MPa
(+18 %), Weibull m = 8 identique (709), GBM alpha = 1 515 MPa (-14 %). Le
maillage Gmsh FRONTAL (algo 6) est BANNI par Fernando (pavage quasi
equilateral, R6 = 0,34 : les cas 1-2 n y cassaient pas) -> Delaunay algo 5 +
taille bruitee (R6 = 0,035). Balayage de resistance a Gf FIXE = artefact
(l_cz x 4-16, pics PLUS hauts) ; a l_cz constant (Gf ~ ft^2) : x0,5 -> 547,
x0,35 -> 496, x0,25 -> 439 MPa, monotone. Regle : parametrer (ft, c, l_cz).

Polydispersite (rockim_f2l.exe) : `grainSizeSpread`, `phase.<nom>.grainSize`.
Trois tentatives, la troisieme est la bonne — voir DOCUMENTATION 5.16 :
Voronoi a espacement par graine (0,5 -> 0,16), Laguerre a poids fixes
(0,5 -> 0,20), **Laguerre a aires prescrites par Newton amorti** (0,5 ->
0,496 ; 0,8 -> 0,80 ; fractions conservees ; 3-8 iterations). Bit-identique
cles absentes (deck GBM, 8 fichiers). Greffes en ajout : `voronoiCells` recoit
des poids et rend le Laplacien (optionnels), la branche polydisperse de la
graine est SEPAREE du Poisson-disc historique (restaure verbatim).

| `rockim_f2l.exe` | + polydispersite Laguerre-Newton + taille par phase | **44/44** (2026-09-02 13:05, suite_f2l.txt) |
| `rockim_f2m.exe` | + historyStrains (epsAx/epsLat/epsVol), gripsStopAfterPeak + stopPeakDrop, weibullScope = lcz, avertissement jointPenaltyFactor inerte (DOC 5.17) | **44/44** (2026-09-02 14:20, suite_f2m.txt) |
| `rockim_f2n.exe` | + grainMeshRandom (Delaunay intra-grain NON structure : R6 0,55 -> 0,007, DOC 5.16 bis) + mesh = file accepte pour geometry = disc (bresilien sur disque Gmsh) | **44/44** (2026-09-02 15:55, suite_f2n.txt) |

Remarque de Fernando 14:00 : « le maillage GBM est structure dans les grains ».
Confirme et mesure (R6 0,548 sur les aretes intra-grain, pire que le frontal
banni) ; corrige en AJOUT par grainMeshRandom. Les resultats GBM anterieurs
(cas 3, juillet, aout) sont a relire avec cette reserve.

## 8. Handoff

Etat complet, plan de runs A-E, recuperation de la revue en cours et regles :
**HANDOFF_2026-09-02.md** (a lire en premier dans une nouvelle session).
Memoire de l agent : project-rockim-schistosite-lisjak.md.
