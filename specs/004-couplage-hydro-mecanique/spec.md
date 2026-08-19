# Feature Specification: Couplage hydro-mécanique 2D

**Feature Branch**: `004-couplage-hydro-mecanique`

**Created**: 2026-08-19

**Status**: Draft — spécification à valider, rien n'est implémenté

**Input**: Donner à rockim un solveur hydraulique couplé, sur le modèle du
solveur d'Y-Geo/Irazu (Lisjak et al. 2017), pour reproduire l'essai de
fracturation hydraulique près du puits d'AbuAisha et al. (2017) — c'est-à-dire
faire ce que le dossier `bench_abuaisha/` déclare aujourd'hui **hors de
portée**.

---

## 1. Contexte : ce qui manque, écrit dans le code lui-même

Le dossier `bench_abuaisha/README.md` pose le diagnostic sans détour. Notre
architecture mécanique **est** celle de leur code Y-Geo : triangles Delaunay
élastiques, éléments cohésifs, mode I de Hillerborg, mode II en
slip-weakening. Ce qui manque est le **fluide**. Et la conséquence tient en une
phrase, présente en toutes lettres dans le source de rockim, au sujet de
`confineFaces = bore` :

> *faces born from cracking receive nothing*

**La pression ne suit pas la fissure.** Le confinement actuel est une pression
suiveuse appliquée aux faces extérieures **d'origine seulement** — choix
délibéré et juste pour une cellule triaxiale, qui a une membrane. Mais une
fracturation hydraulique est exactement le phénomène inverse : le fluide
pénètre la fissure qu'il vient d'ouvrir, et c'est cette pression-là qui la
propage.

D'où le partage acté dans `bench_abuaisha/` : cinq essais reproductibles sans
hydro (tout ce qui se joue **avant** l'amorçage), et hors de portée le plateau
post-pic, le branchement, l'incurvation, l'interaction fissure-joint. Cette
spécification vise précisément ce qui est hors de portée.

---

## 2. Le modèle retenu, et pourquoi celui-là

**Référence de formulation** : Lisjak, Kaifosh, He, Tatone, Mahabadi &
Grasselli, *A 2D, fully-coupled, hydro-mechanical, FDEM formulation for
modelling fracturing processes in discontinuous, porous rock masses*,
**Computers and Geotechnics 81 (2017) 1–18**. C'est le solveur d'Y-Geo, devenu
celui d'Irazu — donc **le même code que celui dont AbuAisha se sert**.

Ses équations sont reprises et explicitées dans un article d'accès libre que
j'ai pu lire intégralement : Shandilaya & Roshankhah, *Fluid
Injection-Induced Fracture Evolution and Breakdown*, **Stanford Geothermal
Workshop 2026**. Toutes les formules ci-dessous en viennent, et sont
attribuées à Lisjak et al. par leurs auteurs.

Ce choix n'est pas seulement pratique : c'est le **même modèle que la cible**.
Une divergence de résultat sera donc imputable à notre implémentation ou à nos
paramètres, jamais à un désaccord de formulation. C'est ce qui rend la
validation falsifiable.

### 2.1 Topologie : cavités et canaux

Le réseau d'écoulement **se superpose au maillage mécanique**, il n'en crée pas
un second.

- une **cavité virtuelle** par sommet du maillage ;
- un **canal d'écoulement** par élément cohésif, reliant les cavités de ses
  deux extrémités ;
- une **cavité physique** optionnelle et explicite, pour le forage : c'est le
  volume dans lequel la pompe injecte.

> **Ce que rockim possède déjà, et qui rend le portage direct.** En FDEM les
> nœuds sont dupliqués : « le sommet » n'est pas un nœud mais un **groupe** de
> copies. Or ce groupe existe déjà et il est tenu à jour — c'est l'union-find
> de l'insertion adaptative (`copiesOfVert_`, `grpsOfVert_`, `vOf_`,
> `rebindVertex()`). Une cavité = un groupe. Et quand un joint s'insère, le
> sommet se scinde : **le réseau hydraulique se scinde avec lui, gratuitement**.
> C'est la même remarque que pour le tétraèdre à pression nodale de la fiche
> D2, et elle vaut ici encore davantage.

### 2.2 Ouverture hydraulique

Pour un canal, l'ouverture aux deux extrémités vaut

$$a_0 = a_i + o_0, \qquad a_1 = a_i + o_1$$

où $a_i$ est l'**ouverture initiale** (celle d'un joint intact — un milieu
poreux n'a pas une conductivité nulle) et $o$ l'ouverture **mécanique** lue au
point d'intégration correspondant. rockim les a déjà : la loi de joint calcule
$d_n$ à ses deux points d'intégration en 2D.

Deux garde-fous de la littérature :
- **ouverture résiduelle** $a_r$ : plancher sous lequel l'ouverture hydraulique
  ne descend pas, même joint fermé ou en compression ;
- **ouverture seuil** $a_t$ : plafond au-delà duquel on cesse d'appliquer la
  loi cubique (au-delà, l'écoulement n'est plus de Poiseuille).

Valeurs de l'article Stanford, pour mémoire : $a_i = a_r = 10^{-5}$ m,
$a_t = 4{,}3\times10^{-4}$ m.

### 2.3 Volume et masse de fluide d'une cavité

Chaque cavité reçoit **la moitié** du volume de chacun de ses canaux :

$$V_c = \sum_{j \ni c} L_j \frac{a_0 + a_1}{2} \cdot \frac{1}{2}$$

et la masse initiale, avec un modèle de compressibilité linéaire :

$$m^{t_0} = S^{t_0}\, V_c\, \rho_f \left(1 + \frac{P^{t_0}}{K_f}\right)$$

$S$ = saturation, $\rho_f$ = masse volumique, $K_f$ = module de compressibilité.

### 2.4 Écoulement : loi de Darcy, résistance par loi cubique

Débit massique dans un canal entre deux cavités de pressions $u_0$ et $u_1$ :

$$q = \frac{dm}{dt} = -\frac{u_1 - u_0 + \rho_f g (y_1 - y_0)}{R}$$

La résistance $R$ suit la **loi cubique**, intégrée le long d'un canal dont
l'ouverture varie linéairement — et l'intégrale a une forme fermée :

$$R = 12\nu \int_{x_1}^{x_0}\frac{\mathrm{d}x}{a(x)^3}
    = \frac{6\nu\,(a_0 + a_1)}{(a_0 a_1)^2}\, L$$

$\nu$ = viscosité cinématique. Noter que $R \to \infty$ quand une ouverture
tend vers zéro : le plancher $a_r$ n'est pas cosmétique, il évite une division
par zéro autant qu'il porte une physique.

### 2.5 Intégration en temps et mise à jour de pression

Euler explicite avant, sur la masse :

$$\Delta m_0 = -\Delta m_1 = f(S)\,\frac{u_1-u_0+\rho_f g(y_1-y_0)}{R}\,\Delta t_h,
\qquad f(S) = S^2(3-2S)$$

Le facteur $f(S)$ éteint progressivement l'écoulement dans une cavité non
saturée — c'est ce qui permet de simuler un front d'invasion plutôt qu'un
milieu partout plein.

Puis la pression, par compressibilité :

$$P^{t} = \begin{cases}
P^{t-1} + K_f\,\dfrac{\Delta m}{\rho_f V_c^{t}} & \text{si } S^t = 1\\[2mm]
0 & \text{si } 0 \le S^t < 1
\end{cases}$$

Une cavité non saturée est donc **à pression nulle** : elle se remplit d'abord,
elle ne pousse qu'ensuite.

### 2.6 Couplage mécanique

Deux sens, et c'est ce qui fait le « fully coupled » :

- **fluide → solide** : la pression de cavité s'applique aux lèvres du joint,
  en force suiveuse sur la configuration courante — exactement la machinerie de
  `confiningForces()`, mais sur les faces **internes** et avec une pression
  **par cavité** au lieu d'une pression globale ;
- **solide → fluide** : l'ouverture mécanique change $a$, donc $R$ (en $a^{-3}$)
  et $V_c$. Une fissure qui s'ouvre devient conductrice ; une fissure qui se
  ferme s'étrangle.

**Schéma décalé, avec sous-cyclage.** Le pas mécanique est fixé par l'onde
(0,9 ns sur `indent3d_grad`) ; le pas hydraulique est fixé par la diffusion et
vaut typiquement des ordres de grandeur de plus. On avance donc $N$ pas
mécaniques par pas hydraulique, avec $N = \Delta t_h/\Delta t$ posé par
configuration. L'article Stanford trace d'ailleurs ses résultats en « temps
hydraulique », qui va jusqu'à 10 s.

---

## 3. Exigences fonctionnelles

- **FR-001** — Un réseau hydraulique est construit sur le graphe des joints :
  une cavité par sommet (groupe de l'union-find), un canal par joint.
- **FR-002** — Le réseau **suit la topologie mécanique** : quand un joint
  s'insère et qu'un sommet se scinde, les cavités se scindent, et le fluide de
  la cavité mère se répartit au prorata des volumes.
- **FR-003** — L'ouverture hydraulique vaut $a = a_i + o$, bornée par $a_r$ et
  $a_t$.
- **FR-004** — La résistance suit la loi cubique sous sa forme fermée à
  ouverture linéaire.
- **FR-005** — La pression suit la compressibilité linéaire, et une cavité non
  saturée est à pression nulle.
- **FR-006** — La pression charge les lèvres en **force suiveuse**, sur la
  configuration courante.
- **FR-007** — Le couplage est **décalé avec sous-cyclage** : `hydroSubSteps`
  pas mécaniques par pas hydraulique.
- **FR-008** — Une **cavité physique** (le forage) peut être déclarée, alimentée
  soit à **débit imposé**, soit à **pression imposée**.
- **FR-009** — Le travail du fluide sur le solide entre au bilan d'énergie B4,
  dans un poste **séparé** (`hydroWork_`). *Après quatre pompes d'énergie
  découvertes en août, un canal de force non instrumenté est exclu.*
- **FR-010** — Tout est **opt-in** : sans `hydro = on`, aucune trajectoire ne
  change d'un bit.

### Clés de configuration proposées

```
hydro          = off | on        # defaut off, bit-identique
fluidDensity   = 1000            # rho_f [kg/m3]
fluidBulk      = 2.2e9           # K_f [Pa]
fluidViscosity = 1.0e-6          # nu [m2/s] (cinematique)
apertureInit   = 1e-5            # a_i [m]
apertureRes    = 1e-5            # a_r [m], plancher
apertureMax    = 4.3e-4          # a_t [m], plafond de validite cubique
hydroSubSteps  = 1000            # pas mecaniques par pas hydraulique
hydroGravity   = 0               # terme rho g (y1 - y0)
boreInjection  = rate | pressure
boreRate       = 1e-5            # [m3/s/m] si rate
borePressure   = 12e6            # [Pa] si pressure
boreVolume     = 0               # volume de la cavite physique
```

---

## 4. Échelle de validation — chaque barreau réfutable

L'ordre est celui de la difficulté croissante, et **chaque barreau doit tomber
avant de monter au suivant**.

| # | essai | référence | critère |
|---|---|---|---|
| **V1** | charge nulle hydro | — | `hydro = on`, pression uniforme partout, aucun gradient : **0 joint cassé**, débit nul machine, `hydroWork` nul |
| **V2** | Kirsch avec pression interne | solution fermée | déjà outillé (`tunnel_edz/tools/kirsch_check.py`), c'est le B3 d'AbuAisha |
| **V3** | **ouverture d'une discontinuité sous pression uniforme** | **Parker (1981), w(0) = 0,0640 mm** | c'est le **B1** du dossier abuaisha, *le seul point de tout leur article confronté à une solution fermée*. Maillage et config **déjà écrits** |
| **V4** | diffusion 1D dans un canal unique | solution de diffusion | vérifie $R$, $K_f$ et le sous-cyclage indépendamment de la mécanique |
| **V5** | pression de rupture, isotrope | **14,2 MPa** (AbuAisha) | leur B2 |
| **V6** | pression de rupture, anisotrope | **12 MPa** (leur éq. 10) | idem, avec $\lambda \ne 1$ |
| **V7** | décalage dû à un joint préexistant | **+5,2 %** | leur B5 |
| **V8** | faciès à l'amorçage | étoile radiale / bi-aile | leur B4 |

**V1 et V3 sont les deux barreaux décisifs** : le premier parce qu'un contrôle à
charge nulle a attrapé trois bugs sur le seul chantier du contact par potentiel ;
le second parce qu'il confronte à une **solution analytique**, ce qui est rare.

Cerise : l'article Stanford donne un contrôle supplémentaire indépendant —
pression de rupture de **92 MPa** contre l'analytique de Hubbert & Willis, avec
$G_{Ic} = 50$ et $G_{IIc} = 500$ J/m² sur un granite à $E = 55$ GPa.

---

## 5. Ce qu'il faut décider avant de coder

1. **Milieu poreux ou fractures seules ?** Lisjak couple aussi la **porosité de
   la matrice** (leak-off dans les triangles). AbuAisha en a besoin pour le
   plateau post-pic. Le faire d'emblée double le chantier ; ne pas le faire
   interdit un de leurs résultats. *Ma recommandation : fractures seules
   d'abord, matrice ensuite, avec la clé prévue dès maintenant.*
2. **2D seulement, ou 2D et 3D ?** La constitution (III) impose les deux de
   front. Mais la feuille de route dit « hydro 2D d'abord », et AbuAisha est un
   cas 2D. *Ma recommandation : invoquer l'exception documentée — l'objet
   n'existe pas encore en 3D — et écrire le 2D seul, en gardant la structure
   portable.*
3. **La cavité physique du forage.** Faut-il un volume de puits explicite avec
   sa compressibilité (donc un stockage, donc une montée en pression réaliste),
   ou une pression imposée en paroi ? *Le premier est nécessaire pour reproduire
   une courbe de pression de puits ; le second suffit pour un seuil d'amorçage.*

---

## 6. Ce qui manque encore

**L'article d'AbuAisha lui-même.** Le dossier `bench_abuaisha/` en cite les
figures, les équations et les valeurs, donc quelqu'un l'a lu — mais le PDF n'est
pas dans la bibliographie. Il sera nécessaire pour V5 à V8 : leur éq. 3
(couplage mixte elliptique), leur éq. 10 (pression de rupture analytique) et
leurs Tables de paramètres.

> AbuAisha, Eaton, Priest & Wong, *Hydro-mechanically coupled FDEM framework to
> investigate near-wellbore hydraulic fracturing in homogeneous and fractured
> rock formations*, **J. Petrol. Sci. Eng. 154 (2017) 100–113**.

**Lisjak et al. 2017** serait un confort, pas une nécessité : l'article Stanford
en donne les équations utiles. Il trancherait en revanche les deux points ouverts
du §5 — porosité de matrice, et traitement de la cavité physique.

---

## 7. Estimation

| étape | contenu | effort |
|---|---|---|
| A | structures (cavités, canaux) + construction sur l'union-find | 2 j |
| B | solveur hydraulique seul (V1, V4) | 2 j |
| C | couplage fluide → solide + bilan d'énergie (V2, V3) | 2 j |
| D | cavité physique, pompe, sous-cyclage (V5, V6) | 2 j |
| E | scission dynamique du réseau à l'insertion (V7, V8) | 2-3 j |

Une dizaine de jours pour l'ensemble, mais **V3 (Parker) tombe dès l'étape C** —
c'est-à-dire qu'on saura si la formulation est juste avant d'avoir tout écrit.
