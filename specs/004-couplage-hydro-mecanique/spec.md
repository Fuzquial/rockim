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

## 4. Échelle de validation

Chaque barreau réfutable, et **chacun doit tomber avant le suivant**.

| # | essai | référence | critère |
|---|---|---|---|
| **V1** | charge nulle hydro | — | `hydro = on`, pompe à zéro : 0 joint cassé, `hydroWork` nul, volume constant |
| **V2** | volume d'une cavité connue | géométrie | un forage circulaire non déformé : le théorème de Green doit rendre $\pi R^2$ |
| **V3** | **ouverture sous pression uniforme** | **Parker (1981), éq. A.1** | **le barreau décisif** — voir ci-dessous |
| **V4** | Kirsch avec pression interne | solution fermée | leur B3, déjà outillé (`tunnel_edz/tools/kirsch_check.py`) |
| **V5** | pression de rupture, anisotrope | **12 MPa** (leur éq. 10), ~12,5 chez eux | leur B2 |
| **V6** | pression de rupture, isotrope | **14,2 MPa** | idem |
| **V7** | décalage dû à un joint préexistant | **+5,2 %** longitudinal et oblique, ~0 transverse | leur B5 |
| **V8** | faciès à l'amorçage | étoile radiale (isotrope) / bi-aile (anisotrope) | leur B4 |

### V3 — le contrôle décisif, entièrement spécifié

C'est **le seul point de tout leur article confronté à une solution fermée**,
et leur annexe A le donne complètement.

**Solution de Parker (1981), p. 33**, ouverture d'une discontinuité sous
pression uniforme en déformation plane :

$$w(x) = \frac{2\sigma'(1-\nu^2)}{E}\sqrt{\ell^2 - x^2} \tag{A.1}$$

avec $\sigma' = p - \sigma_n$ la contrainte effective d'ouverture, $\ell$ la
demi-longueur, $p$ la pression interne.

**Le cas** : domaine 8 × 8 m, discontinuité de 1,5 m ($\ell = 0{,}75$ m),
$p = 12$ MPa, $\sigma_H = 15$ MPa, $\sigma_v = 10$ MPa, donc $\sigma' = 2$ MPa.
Maillage Delaunay de 0,003 à 0,3 m. Résistances portées à des valeurs
*« irreal high »* pour interdire toute propagation : **c'est de l'élasticité
pure sous charge suiveuse**, rien d'autre n'est testé.

$$w(0) = \frac{2 \times 2\cdot10^6 \times (1-0{,}2^2)}{45\cdot10^{9}} \times 0{,}75
       = 6{,}4\times10^{-5}\ \text{m} = \mathbf{0{,}064\ mm}$$

> ⚠️ **Coquille de l'article**, déjà repérée dans `bench_abuaisha/README.md` :
> l'annexe écrit « E = 45 MPa ». C'est **45 GPa** — sans quoi l'ouverture
> vaudrait 64 m. Leur figure A.21 lit bien 0,065 mm.

**Le maillage et la config de ce cas sont déjà écrits** dans `bench_abuaisha/`
(`make_crack_mesh.py`, `parker.cfg`, `parker_check.py`), avec le dédoublement
des lèvres par le greffon `Crack` de gmsh et l'astuce de la contrainte nette
(appliquer 2 MPa dans un milieu non précontraint plutôt que 12 contre 10, pour
éviter des lèvres plaquées à $t = 0$).

### V5/V6 — le seuil analytique

$$p^{HF} = \sigma'_H - 3\sigma'_h + f_t \tag{10}$$

(Fjaer et al. 2008). Pour leur cas : $-6{,}8 - 3\times(-4{,}6) + 5 = 12$ MPa,
et leur code donne ~12,5.

### Paramètres de leur Table 1

| | valeur |
|---|---|
| $E$, $\nu$, $\rho$ | 35 GPa, 0,27, 2500 kg/m³ |
| $f_t$, $c$, $\varphi_i = \varphi_f$ | 5 MPa, 24 MPa, 38° |
| $G_{Ic}$, $G_{IIc}$ | 10 N/m, 80 N/m |
| pénalités $p_n$, $p_t$, $p_f$ | 350, 35, 175 GPa·m |
| amortissement $\mu$ | $5{,}6\times10^{5}$ kg/m/s |

> **Deux remarques sur ces valeurs.** (i) $G_{Ic} = 10$ J/m² et $f_t = 5$ MPa
> avec $E = 35$ GPa donnent $\ell_{cz} = EG_f/f_t^2 = 14$ mm, pour une maille de
> 3 mm et un forage de 100 mm : **la règle maison est respectée des deux
> côtés**, ce que `bench_abuaisha/README.md` avait déjà relevé. (ii) Leur
> amortissement est **massif** — c'est une relaxation dynamique assumée, le
> problème étant quasi-statique. À rapprocher du balayage B1 en cours : le
> même paramètre, deux régimes opposés, deux réponses opposées.

---

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
