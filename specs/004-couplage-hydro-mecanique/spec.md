# Feature Specification: Couplage hydro-mécanique 2D

**Feature Branch**: `004-couplage-hydro-mecanique`

**Created**: 2026-08-19 — **révisée le même jour** après lecture de l'article
d'AbuAisha transmis par F. Uzquiano.

**Status**: Draft — spécification à valider, rien n'est implémenté

**Input**: Donner à rockim un solveur hydraulique couplé permettant de
reproduire l'essai de fracturation hydraulique près du puits d'AbuAisha,
Eaton, Priest & Wong (JPSE 154, 2017, 100–113) — c'est-à-dire faire ce que le
dossier `bench_abuaisha/` déclare aujourd'hui **hors de portée**.

---

## 0. Correction de la première rédaction

> La version initiale de cette spec, écrite avant d'avoir l'article, posait
> comme cible la formulation de **Lisjak et al. 2017** : réseau de cavités et
> de canaux, écoulement de Darcy, loi cubique, sous-cyclage hydraulique.
>
> **Ce n'est pas ce que fait AbuAisha.** Leur module HF d'Y-Geo repose sur une
> hypothèse bien plus simple, qu'ils énoncent explicitement (leur §2.2) :
>
> > *For each time step the fluid pressure inside the cavities/fractures is
> > assumed **constant**, i.e. **no flow due to hydraulic gradients** is
> > considered, i.e. similar to an **inviscid flow**.*
>
> Pas de loi cubique, pas de gradient de pression, pas de leak-off — ils
> l'écrivent aussi en introduction : *« the FDEM–HF in its current form assumes
> inviscid flow restricted to fractures and it does not account for leakoff
> into the rock formation »*.
>
> Conséquence : **le chantier est environ deux fois plus court que je ne
> l'avais estimé**, et la formulation de Lisjak devient une extension
> ultérieure, pas le point de départ.

---

## 1. Contexte : ce qui manque, écrit dans le code lui-même

`bench_abuaisha/README.md` pose le diagnostic. Notre architecture mécanique
**est** la leur : triangles Delaunay élastiques, éléments cohésifs, mode I de
Hillerborg, mode II en slip-weakening, couplage mixte elliptique. Ce qui manque
est le **fluide**. Et la conséquence tient dans un commentaire du source de
rockim, au sujet de `confineFaces = bore` :

> *faces born from cracking receive nothing*

**La pression ne suit pas la fissure.** Le confinement actuel s'applique aux
faces extérieures **d'origine seulement** — juste pour une cellule triaxiale,
qui a une membrane ; faux pour une fracturation hydraulique, où le fluide
pénètre la fissure qu'il vient d'ouvrir et c'est cette pression qui la propage.

---

## 2. Le modèle d'AbuAisha, dans le détail

Trois modules en boucle (leur Fig. 1) : le solveur mécanique FDEM, un
**calculateur de volume de cavité**, et un **modèle de pompe**.

### 2.1 L'hypothèse centrale : fluide non visqueux

La pression est **uniforme dans toute la cavité connectée** et instantanée. Il
n'y a donc ni réseau d'écoulement, ni résistance, ni pas de temps hydraulique
propre. Leur justification : *« appropriate as long as we are investigating
near wellbore behaviour where convection effects are dominant »*.

C'est une limite assumée, et elle borne le domaine de validité — mais c'est
aussi ce qui rend le chantier abordable.

### 2.2 Volume de cavité, par le théorème de Green

$$V = \frac{1}{2}\oint x\,\mathrm{d}y - y\,\mathrm{d}x \tag{4}$$

évalué sur **toutes** les frontières mouillées, anciennes et nouvelles.

> **Détail d'implémentation qui simplifie tout** : discrétisée, cette intégrale
> est une somme sur les segments de frontière,
> $V = \tfrac12\sum (x_1 y_2 - x_2 y_1)$, et elle est **indépendante de l'ordre
> de parcours** tant que les orientations sont cohérentes. Il n'y a donc
> **aucun contour à assembler ni à ordonner** — c'est la formule du lacet
> sommée face par face. rockim oriente déjà ses faces extérieures.

### 2.3 Suivi des frontières mouillées

Le module *« tracks the newly created wet boundaries by checking their
**connection with the initial source of fluid** »*. Une face devient mouillée
si elle est reliée au forage par un chemin de joints ouverts.

> **Ce que rockim possède déjà** : c'est une recherche de composante connexe,
> exactement la machinerie de `computeFragments()` et de l'union-find de
> l'insertion adaptative — mais parcourue sur les joints **rompus/morts** au
> lieu des joints vivants. Et `confineFaces = bore` sait déjà désigner les
> faces du forage, c'est-à-dire la source.

### 2.4 Compressibilité et pompe

$$K_f = -V\frac{\mathrm{d}p}{\mathrm{d}V} = \rho_f\frac{\mathrm{d}p}{\mathrm{d}\rho_f} \tag{5}$$

$$p = p_0 + K_f \log\!\left(\frac{m}{V\rho_{f0}}\right) \tag{6}$$

où $m$ est la masse injectée, intégrée du débit de pompe. Une seule variable
d'état scalaire par cavité : $m$. La pression s'en déduit à chaque pas.

### 2.5 Chargement du solide

Force nodale sur une frontière mouillée définie par les nœuds 1 et 2 :

$$\mathbf{F}_{p12} = -\frac{p}{2}\begin{bmatrix} y_2 - y_1 \\ x_2 - x_1\end{bmatrix} \tag{7}$$

C'est une **force suiveuse** sur la configuration courante — la machinerie de
`confiningForces()` au signe et à l'assiette près.

### 2.6 Joints préexistants

*« if at any point the rock joint is intersected by a fluid-driven fracture,
the fluid pressure percolation is distributed **evenly** over the entire
discontinuity »* — le joint entier devient mouillé d'un coup, à la même
pression. Conséquence assumée par les auteurs : concentration de contrainte aux
**pointes** du joint, d'où l'amorçage préférentiel qu'ils observent.

---

## 3. Exigences fonctionnelles

- **FR-001** — Une **cavité fluide** est définie par l'ensemble des faces
  mouillées connectées à une source.
- **FR-002** — L'appartenance à la cavité se propage par **connexité** à
  travers les joints rompus, réévaluée quand la topologie change.
- **FR-003** — Le volume suit le théorème de Green, sommé face par face, sur la
  configuration **courante**.
- **FR-004** — La pression suit la compressibilité linéaire (éq. 6) à partir de
  la masse injectée.
- **FR-005** — La pression charge les faces mouillées par l'éq. 7, en force
  suiveuse.
- **FR-006** — La pompe impose soit un **débit**, soit une **pression**.
- **FR-007** — Un joint préexistant intersecté est mouillé **en entier**.
- **FR-008** — Le travail du fluide entre au bilan B4 dans un poste **séparé**
  (`hydroWork_`). *Après quatre pompes d'énergie découvertes en août, un canal
  de force non instrumenté est exclu.*
- **FR-009** — Tout est **opt-in** : sans `hydro = on`, aucune trajectoire ne
  change d'un bit.
- **FR-010** — Le solveur imprime à chaque frame : volume de cavité, pression,
  masse injectée, nombre de faces mouillées. *Une cavité qui se remplit sans
  qu'on puisse la regarder est ingouvernable.*

### Clés proposées

```
hydro           = off | on       # defaut off, bit-identique
hydroSource     = bore           # faces sources (reutilise confineFaces)
fluidBulk       = 2.2e9          # K_f [Pa]
fluidDensity    = 1000           # rho_f0 [kg/m3]
hydroP0         = 0              # p_0 [Pa]
hydroInjection  = rate | pressure
hydroRate       = 1e-5           # [m3/s/m]
hydroPressure   = 12e6           # [Pa] si pressure
hydroWetJoints  = whole | partial  # FR-007 : joint entier ou face par face
```

---

## 4. Échelle de validation, et reproduction des figures

### 4.0 Correction : Parker est DÉJÀ validé, et sans hydro

> Ma première rédaction plaçait l'ouverture de Parker comme « le barreau
> décisif » du couplage. **C'est faux**, et `bench_abuaisha/README.md` le disait
> déjà : *« B1 est sans hydro par construction — pression uniforme imposée, et
> résistances mises à des valeurs irreal high pour interdire toute
> propagation »*. Ce cas teste une **charge suiveuse en élasticité pure**, pas
> le fluide.
>
> Il a été couru le 2026-08-18 et **il passe** :
> `parker_vs_theorie_2maillages.png` donne **−5,65 %** au centre sur le
> maillage grossier au plateau, **+3,16 %** sur le raffiné à 6 ms — celui-ci
> n'ayant pas atteint son plateau, il est en dépassement dynamique. La forme
> colle à l'ellipse sur toute la demi-longueur.

L'échelle repart donc plus bas, et le premier vrai test du fluide est un
**pont de non-régression**.

### 4.1 Les barreaux du module

| # | essai | critère |
|---|---|---|
| **H1** | charge nulle hydro | pompe à zéro : 0 joint cassé, `hydroWork` nul, volume constant |
| **H2** | volume de Green | **FAIT** — forage circulaire : 5,6e-07 d'écart à l'aire du polygone maillé |
| **H3** | **pont Parker** | rejouer le cas Parker en `hydro = on, hydroInjection = pressure` et retrouver **les mêmes −5,65 % / +3,16 %** que par `confiningPressure`. Physiquement le même problème : tout écart est un bug de cavité, de connexité ou de signe |
| **H4** | compressibilité | en `hydroInjection = rate`, la pression doit suivre $p_0 + K_f\log(m/Vho_0)$ pendant que la fissure s'ouvre |
| **H5** | connexité | une fissure qui se propage doit **mouiller de nouvelles faces** et faire croître le volume. C'est précisément ce que `confineFaces = bore` ne sait pas faire |

H3 à H5 testent le fluide ; H1 et H2 le protègent.

### 4.2 Reproduction des figures de l'article

**Le problème aux limites (leur Fig. 6)** : forage nu de 0,1 m dans un bloc de
8 × 8 m, déformation plane, épaisseur unité. Maillage Delaunay de **3 mm** dans
une zone de 0,8 × 0,8 m² autour du puits, grossissant jusqu'à 0,3 m au loin.
Formation granitique à 1800 m.

**Contraintes (leur §3.2)**, en convention « négatif = compression », et
**effectives** (leur pression de formation de 28 MPa est déjà retranchée) :

| état | $\sigma'_h$ | $\sigma'_H$ |
|---|---|---|
| isotrope | $-4{,}6$ MPa | $-4{,}6$ MPa |
| anisotrope | $-4{,}6$ MPa | $-6{,}8$ MPa |

**Procédure** : étape géostatique sans le puits jusqu'à équilibre, puis
excavation, puis injection à **Q = 20 l/s** (choisi sous le seuil de ~100 l/s
au-delà duquel la pression de rupture augmente).

> **Ce que rockim a déjà pour ça** : `insituSh`/`insituSv` pour la
> pré-contrainte et `excavRelease` pour l'excavation, tous deux **validés
> contre Kirsch à 1,7–2,1 %** dans l'étude tunnel. Leur excavation à eux passe
> par une réduction de module puis un retrait d'éléments ; la nôtre par
> convergence-confinement — même état initial, même état final.

| # | figure | ce qu'on compare | cible chiffrée |
|---|---|---|---|
| **F1** | **Fig. 11(a)**, cas de référence sans joint | **pression de rupture**, anisotrope | **12 MPa** analytique (leur éq. 10), **~12,5** chez eux |
| **F2** | idem, isotrope | pression de rupture | **14,2 MPa** |
| **F3** | **Fig. 9** | champ de $\sigma_{11}$ effectif juste avant amorçage | qualitatif + Kirsch |
| **F4** | **Fig. 7** | **faciès** à $t = 1{,}2$ et $1{,}26$ ms | isotrope → étoile radiale complexe ; anisotrope → **bi-aile**, puis branchement et incurvation vers $\sigma_H$ |
| **F5** | **Fig. 11**, 9 cas à joint | décalage du seuil selon orientation et distance | **+5,2 %** en longitudinal et oblique ; **~0 %** en transverse |
| **F6** | **Fig. 11**, post-pic | plateau de propagation | **5,5 MPa** (la contrainte lointaine effective), **indépendant de la distance du joint** |
| **F7** | **Fig. 10, 12, 13** | trajectoires selon la distance du joint (0,3L / 0,6L / L) | l'effet du joint s'estompe avec la distance ; à $L$, quasi identique au cas sans joint |
| **F8** | **Fig. 15** | formation très fracturée : 320 joints/m², longueur moyenne 3,1 cm | faciès à $t = 1{,}36$ ms |

**F1 et F2 sont les barreaux quantitatifs** — deux nombres, une solution
analytique. **F6 en est un troisième**, et il est plus discriminant qu'il n'en
a l'air : le plateau post-pic doit valoir la contrainte lointaine **et ne pas
dépendre** de la distance du joint. Le reste est qualitatif.

**Hors de portée, et il faut le dire** : leur Fig. 19 (essai de terrain
Montney) demande des données de champ ; leurs Fig. 16-17 (microsismicité) sont
un post-traitement de la masse et de la vitesse nodales — que rockim sort déjà,
mais dont le semis suivra le trajet de fracture, donc ne coïncidera pas au-delà
de l'amorçage.

**Réserve de méthode, posée avant les runs.** Leur modèle est non visqueux :
la pression est uniforme et instantanée dans toute la cavité. Un faciès qui
diffère du leur ne prouvera donc **ni** que notre mécanique est fausse, **ni**
que la leur l'est — la FDEM est dépendante du maillage, ils le disent
eux-mêmes et emploient une discrétisation aléatoire pour limiter le biais
d'orientation. Les nombres (F1, F2, F6) départageront ; les images
illustreront.

## 5. Décisions en attente

1. **2D seulement, ou 2D et 3D ?** La constitution (III) impose les deux de
   front ; la feuille de route dit « hydro 2D d'abord » ; AbuAisha est 2D.
   *Recommandation : invoquer l'exception documentée et écrire le 2D seul.*
2. **Une cavité ou plusieurs ?** Leur modèle suppose une source unique. Un
   fragment détaché qui emporte du fluide, ou deux fissures qui se rejoignent,
   demandent une gestion multi-cavités (fusion/scission). *Recommandation :
   une cavité pour commencer, structure prévue pour plusieurs.*
3. **Extension Lisjak ensuite ?** Le modèle non visqueux ne peut pas rendre le
   plateau post-pic ni le leak-off. La formulation de Lisjak (réseau, loi
   cubique, sous-cyclage) reste spécifiée dans l'historique git de ce fichier
   et pourra devenir `hydroModel = inviscid | network`.

---

## 6. Estimation révisée

| étape | contenu | effort |
|---|---|---|
| A | cavité : connexité aux sources, volume de Green, sortie | 1,5 j |
| B | compressibilité + pompe + chargement éq. 7 + bilan B4 | 1,5 j |
| C | **V1, V2, V3 (Parker)** | 1 j |
| D | joints préexistants mouillés en entier (V7) | 1 j |
| E | seuils de rupture, faciès (V5, V6, V8) | 1,5 j |

**Environ 6 jours**, contre la dizaine estimée avec le modèle de Lisjak. Et
**V3 tombe dès l'étape C** : on saura si la formulation est juste après trois
jours, pas après dix.

---

## 7. Sources

- **AbuAisha, Eaton, Priest & Wong**, *Hydro-mechanically coupled FDEM
  framework to investigate near-wellbore hydraulic fracturing in homogeneous
  and fractured rock formations*, **JPSE 154 (2017) 100–113**. Classé dans la
  bibliographie.
- **Parker (1981)**, *The mechanics of fracture and fatigue*, p. 33 — la
  solution analytique de V3, reprise en annexe A de l'article ci-dessus.
- **Fjaer et al. (2008)**, *Petroleum related rock mechanics* — le seuil de
  l'éq. 10.
- **Lisjak, Kaifosh, He, Tatone, Mahabadi & Grasselli**, Computers and
  Geotechnics 81 (2017) 1–18 — la formulation **visqueuse** d'Y-Geo/Irazu, pour
  l'extension ultérieure. Ses équations sont reprises en accès libre par
  Shandilaya & Roshankhah, Stanford Geothermal Workshop 2026 (classé dans la
  bibliographie).
