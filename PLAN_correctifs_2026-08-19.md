# PLAN DE CORRECTIFS — manques physiques et numériques de rockim

*Rédigé le 2026-08-19. **Proposition à valider : rien n'est implémenté.** Chaque
fiche donne le manque mesuré, ce que dit la littérature, la correction proposée
avec son API, le critère qui la réfuterait, et le risque de non-régression.*

---

## 0. Synthèse

| # | chantier | classe | effort | ce que ça débloque |
|---|---|---|---|---|
| **A1** | Contact outil en **condition de vitesse** (CD-Lagrange par nœud) | numérique | **court** | supprime la pompe ×408 ; toute comparaison quantitative en coupe/indentation |
| **A2** | Contact général **non lisse** (Nonsmooth Newmark, QP sur ensemble actif) | numérique | gros | découple le pas de temps de la raideur de contact (×1000 rapporté) ; supprime la dérive d'énergie |
| **A3** | **Cap de raideur cohésive** à faible endommagement | numérique | court | évite que dt s'effondre quand D → 0 (prérequis de A2) |
| **B1** | `dampingLocal = 0` en dynamique + balayage falsifiant | physique | *campagne* | rend défendable **tout** bilan énergétique |
| **B2** | Deux observables anti-pompe imprimées à chaque run | numérique | **très court** | détecte la classe de défaut que le résidu B4 ne voit pas |
| **C1** | Mass scaling **sélectif d'Olovsson** (blocs 4×4 inversés une fois — la duplication des nœuds le rend possible) | numérique | moyen | dt ×3 à ×10 **sans ajouter de masse** ni toucher la translation rigide |
| **C2** | Checkpoint / restart | numérique | moyen | campagnes de plusieurs jours |
| **C3** | VTU binaire | numérique | court | 360 Mo/frame → ~40 ; l'E/S cesse de dominer |
| **C4** | `MatState` par loi | numérique | court | ×4-5 mémoire ; débloque les gros maillages |
| **D1** | Loi de volume **avec cap** pour l'indentation | physique | *config* | supprime les 1,75 GPa de pression de contact non physiques |
| **D2** | Tétraèdre à **pression nodale moyennée** (Bonet & Burton ; COMPTet en repli) | physique | **moyen** | lève le verrouillage volumique de la zone broyée, sans changer l'élément |
| **E** | Lot de dette silencieuse (8 points) | numérique | **demi-journée** | supprime une famille entière de faux résultats |
| **F** | Validation contre les articles + recalibration | *campagne* | gros | passe de « vérifié » à « validé » |

---

## A. LE CONTACT — le verrou principal

### A1. Contact outil : passer d'une condition en pénétration à une condition en VITESSE

**Le manque, mesuré.** L'écrêtage actuel est *géométrique* : il borne la
pénétration à 0,6 h sans référence à la masse nodale ni au pas de temps. Il ne
borne donc pas l'**impulsion**, qui est la grandeur qui lance les nœuds. Mesure
du 2026-08-18 (coupe PDC) : l'outil injecte **77 286 J/m** dans le solide pour
**189 J/m** de travail de corps rigide — rapport **408** — et lance des nœuds à
**2 544 m/s**, soit **127 fois** la borne physique 2·v_outil. Le correctif
`toolImpulseCap` posé le même jour ramène le rapport à 57, mais sa formulation
est fausse : il plafonne l'incrément **par pas**, alors que la borne physique
porte sur l'impulsion de **toute la collision**. Un nœud en contact soutenu
prend 20 m/s à chaque pas et seize pas donnent les 311 m/s mesurés.

**Ce que dit la littérature.** C'est exactement le schéma **CD-Lagrange**
(Fekak, Brun, Gravouil et al. 2017 ; exposé pédagogique dans Dureisseix et al.,
*JTCAM* 2024, accès libre). Le contact n'y est pas pénalisé : il est écrit comme
une relation **impulsion / saut de vitesse** (lemme de viabilité de Moreau),

```
si g > 0        alors r = 0
sinon           0 ≤ v ⊥ r ≥ 0
```

et la dynamique réduite s'écrit `v = v_libre + H·r`, où `H = L M⁻¹ Lᵀ` est
l'opérateur de Delassus.

> **Le point décisif pour rockim** : l'article note que *« due to the properties
> of the lumped mass matrix M, this leads to a Delassus operator which is
> diagonal, definite positive and spherical per node »*. La masse de rockim est
> **déjà diagonale**, et le contact outil est un contact **nœud contre obstacle
> rigide** : `L` est donc local à chaque nœud, `H` est un scalaire, et
> l'impulsion se calcule **en forme fermée, nœud par nœud**. Aucun système à
> résoudre, aucune itération : le solveur reste matrix-free.

**La correction proposée.** Dans `toolContact()`, remplacer la force de pénalité
par une impulsion de Signorini résolue par nœud :

```
1. v_libre,i  = v_i + (dt/m_i)·f_i          (vitesse libre, forces déjà sommées)
2. g_{n+1}    = gap prédit du noeud i au pas suivant
3. si g_{n+1} > 0 :  r_i = 0
   sinon           :  v_rel = (v_libre,i - v_outil)·n
                      r_i   = max(0, -m_i·v_rel)      <- forme fermée, H = 1/m_i
4. f_i += r_i·n/dt        (report en force pour que le compteur toolWork reste juste)
5. frottement : cap de Coulomb sur l'impulsion tangentielle (bipotentiel,
   Fekak 2017), pas sur la force
```

Clé : `toolContact = penalty | signorini` (défaut `penalty`, bit-identique).

**Ce que ça change dans le budget du pas de temps.** `kp_ = E·t` sort du calcul
de `computeStableDt()` en mode `signorini` : la borne ressort-masse perd le
terme `nExtra·k_contact`. Sur `indent3d_grad`, ce terme n'est pas dominant, mais
sur un cas SHPB ou une coupe il l'est.

**Critère de falsification, binaire.** Après correctif, sur la config `cut_v3`
inchangée par ailleurs :
- **aucun nœud au-dessus de 2·v_outil** (aujourd'hui : 2 544 m/s pour 20 attendus) ;
- **rapport injection / travail de corps rigide ≤ 2** (aujourd'hui 408) ;
- contrôle à charge nulle : 0 joint cassé, travail de contact exactement nul ;
- suite `fast` bit-identique clé éteinte.

Si le rapport ne descend pas sous 2, l'hypothèse « la pompe est dans l'outil »
est réfutée et il faut rouvrir le diagnostic.

**Risque.** L'interpénétration résiduelle ne disparaît pas (l'article la
quantifie par η = −min g / max|g| ≈ 0,4 % sur son cas 0D) : on échange une
force fausse contre une petite violation de contrainte, ce qui est le bon
échange. Le pic de force **va baisser** : il faut vérifier qu'il ne tombe pas
sous la cible Heilman de 3,08 MN/m comme l'a fait la combinaison EPFL + écrêtage
(0,798 MN/m, quatre fois trop bas).

---

### A2. Contact général : Nonsmooth Newmark sur l'ensemble actif

**Le manque.** Deux problèmes distincts, tous deux mesurés : la déclaration
tardive des faces neuves porte un facteur **6 à 8** sur toutes les observables,
et le correctif évident (déclarer plus tôt) **aggrave** — l'injection suit le
*nombre* de basculements cohésif/contact, pas leur retard. La pénalité
adaptative `k⁻ = k⁺(D)` divise l'injection par 14 mais ses auteurs la
qualifient eux-mêmes de *« diagnostic rather than a definitive remedy »*.

**Ce que dit la littérature.** Le remède est dans l'article **compagnon** de
celui déjà cité : Ghesquière-Diérickx, Anciaux, Acary & Molinari,
*The semi-explicit nonsmooth Newmark time integrator for robust unilateral
contact in dynamic fragmentation simulations*, **arXiv:2606.01355** (mai 2026,
accès libre) — le même problème que le nôtre : éléments finis + CZM extrinsèque
+ contacts denses entre fragments.

Leur schéma **NSN** :
- le volume et la fissuration sont intégrés **explicitement** (Newmark β = 0,
  γ = 1/2), le contact **implicitement** (Euler) ;
- le contact est un **programme quadratique convexe** sur le seul ensemble
  actif : `min_{p≥0} ½ pᵀW'p + pᵀb`, avec `W'` l'opérateur de Delassus modifié ;
- l'ensemble actif ne contient que les nœuds à gap prédit négatif, donc **le
  coût suit le nombre de contacts, pas la taille du système** ;
- **le pas de temps cesse de dépendre de la raideur de contact** : gain de trois
  ordres de grandeur rapporté sur leur barre impactée endommagée
  (10⁻¹⁰ s en pénalité contre 10⁻⁷ s en NSN).

**La correction proposée, en deux temps.**

1. **Étape courte** : appliquer A1 (forme fermée par nœud) au contact
   nœud-face **sur les faces exposées à un corps rigide ou à un plan** — même
   algèbre, `H` reste diagonal.
2. **Étape longue** : pour le contact déformable-déformable (nœud-face et
   élément-élément du potentiel), `H` n'est plus diagonal (un nœud en contact
   avec une facette couple 4 nœuds). Il faut alors un **Gauss-Seidel projeté**
   sur l'ensemble actif, quelques dizaines d'itérations par pas, parallélisable
   par blocs de contacts disjoints. Clé `contact = penalty | potential |
   nonsmooth`.

**Critère de falsification.** Leur propre banc : barre impactée avec interfaces
cohésives, où la solution non lisse et la pénalité doivent converger vers la
même réponse quand la pénalité tend vers l'infini — sauf que NSN doit y arriver
avec un pas 100 à 1000 fois plus grand. Plus, chez nous : le contrôle à charge
nulle, la conservation du potentiel (ΔKE/KE₀ doit rester à 10⁻¹²), et la
percussion 2D de référence.

**Risque et arbitrage à trancher avant de commencer.** C'est un chantier de
plusieurs semaines qui touche l'intégrateur. Il ne doit pas démarrer avant que
A1 et B1 aient été mesurés : si A1 suffit à ramener les faciès dans le domaine
physique, A2 devient un chantier de fond, pas d'urgence.

---

### A3. Cap de raideur cohésive à faible endommagement

**Le manque.** À D → 0, la sécante de pénalité vaut `p_J = 20E/h` — c'est elle
qui domine la borne ressort-masse et donc le pas de temps. En insertion
adaptative on est déjà à `4E/h`, mais le problème demeure à chaque insertion.

**Littérature.** C'est le §2.4 de arXiv:2606.01355 : ils modifient la loi de
traction-séparation pour **plafonner la raideur cohésive** et éviter la
restriction prohibitive de pas quand l'endommagement approche zéro.

**Proposition.** Clé `jointStiffnessCap` (défaut désarmé) : borner `p_J` par
`κ·m_min/(A·dt²)` de sorte que le joint n'impose jamais un pas plus petit qu'une
fraction fixée du pas élastique. À chiffrer avant de coder : sur `indent3d_grad`
la borne ressort vaut 5,97e-9 s contre 1,16e-8 s pour la CFL, donc les joints
coûtent déjà un **facteur 1,9** sur le pas — c'est exactement ce que ce cap
récupérerait.

---

## B. AMORTISSEMENT ET COMPTABILITÉ ÉNERGÉTIQUE

### B1. `dampingLocal` : la valeur par défaut est probablement fausse en dynamique

**Le manque, mesuré trois fois.** Sur `impact3d_ultra`, l'amortissement de
Cundall dissipe **0,402 J contre 0,165 J** pour les joints — 2,4 fois la
physique. Sur `indent2d_yan`, 56,7 J/m contre 5,1 — **onze fois**. Sur
`indent3d_grad` (19/08), 0,0914 J contre 0,0035 — **vingt-six fois**. La valeur
0,05 est héritée d'un banc antérieur et n'a jamais été justifiée.

**Ce que dit la littérature, et qui tranche.** L'amortissement local de Cundall
(1987) est un outil de **mise à l'équilibre et de quasi-statique** : *« for
compact particle models, non-zero local damping can be used to establish
equilibrium and quasi-static deformation simulation. However, in dynamic
analysis, the energy in the model system should be dissipated based on the
viscous damping of the contact model »*. En FDEM, la dissipation dynamique doit
venir du **dashpot de joint**, des **frontières de Lysmer** et, si on la veut,
de la **viscosité de volume de Yan**.

**Proposition — et c'est le chantier n° 1.** Ce n'est pas du code, c'est une
campagne de trois runs sur `indent3d_grad`, qui est désormais notre cas de
référence :

| run | `dampingLocal` | ce qu'on regarde |
|---|---|---|
| D0 | 0,05 (référence, déjà fait) | 26× la fissuration |
| D1 | 0,01 | le faciès change-t-il ? le pic ? |
| D2 | **0** | la dissipation restante est-elle portée par joints + Lysmer + dashpot ? |

**Critère de décision.** Si le pic de force et le compte de fissures varient de
moins de 10 % entre D0 et D2, alors `dampingLocal` ne fait que **manger de
l'énergie sans changer la physique** : on le met à 0 en dynamique et tous les
bilans deviennent lisibles. S'ils varient beaucoup, c'est que le calcul
s'appuyait sur lui pour rester stable — et c'est un tout autre problème, à
traiter par A1/A2.

**Note de coût.** Chaque run ≈ 5 h. Les trois peuvent s'enchaîner en une nuit.

### B2. Imprimer les deux observables qui exposent une pompe

**Le manque.** Le résidu B4 est **structurellement aveugle** à une injection
logée dans un canal comptabilisé : il boucle à 10⁻¹⁰ % pendant que la physique
est détruite. C'est arrivé **quatre fois** dans la campagne d'août.

**Proposition (très court, aucune physique touchée).** Ajouter au résumé de fin
de run, dès qu'un outil ou un contact existe :

```
[FDEM] injection outil / travail de corps rigide : 408.0   <-- ALERTE si > 2
[FDEM] vitesse nodale max / 2 v_outil            : 127.2   <-- ALERTE si > 1
```

Et les verrouiller par un repère de la suite sur la percussion 2D.

---

## C. PAS DE TEMPS, MÉMOIRE, COÛT

### C1. Mass scaling — la duplication des nœuds nous ouvre la bonne solution

*Fiche révisée le 2026-08-19 après lecture d'Olovsson, Simonsson & Unosson
(IJNME 63, 2005, 1436–1445) et de Tkachuk & Bischoff (Comput. Mech. 52, 2013,
563–570), transmis par F. Uzquiano.*

**Le manque.** Aucun mass scaling n'existe dans le code. Or les éléments qui
étranglent le pas font **moins de 0,6 % du volume**. C'est la condition pour
raffiner sous le contact sans payer en 1/h⁴.

**Deux familles, et ce que dit vraiment la littérature.**

*Le mass scaling conventionnel* alourdit les éléments lents. Il est trivial mais
il **ajoute de la masse**, donc modifie l'inertie là où il agit.

*Le mass scaling SÉLECTIF* (SMS) ajoute à la matrice de masse un terme
$\boldsymbol\lambda$ tel que $\tilde{\mathbf M} = \mathbf M + \boldsymbol\lambda$
avec la contrainte **$\boldsymbol\lambda\,\ddot{\mathbf u}_r = 0$ pour toute
accélération de corps rigide** : la masse translationnelle de l'élément est
*exactement* préservée, seules les fréquences élevées sont abaissées. Olovsson
en donne deux formes :

- **Méthode I** : $\boldsymbol\lambda = \alpha\,\mathbf k$ (proportionnelle à la
  raideur). Les modes propres sont **préservés** et les fréquences deviennent
  $\tilde\omega_i^2 = \omega_i^2/(1+\alpha\omega_i^2)$ : à
  $\alpha = 100/\omega_{\max}^2$, la fréquence maximale chute d'un ordre de
  grandeur pendant que le bas du spectre bouge peu. **Mais** $\mathbf k$ change
  à chaque pas en grandes transformations : les auteurs écrivent eux-mêmes que
  c'est *« not practical in geometrically nonlinear applications »* et qu'ils
  *« cannot afford to invert $\tilde{\mathbf M}$ every time step »*. Écartée.
- **Méthode II** : $\boldsymbol\lambda$ **algébrique**, ne dépendant que du type
  d'élément et de sa masse, donc **constante** sous déformation et rotation. Pour
  l'hexaèdre à 8 nœuds, $\boldsymbol\lambda_{8\times8} = \frac{\beta m^e}{56}
  [\,7\ \text{sur la diagonale},\ -1\ \text{hors diagonale}\,]$, sans couplage
  entre $x$, $y$ et $z$. C'est la forme disponible dans **LS-DYNA et RADIOSS**.
  Gain mesuré sur leur poutre : 12 900 pas → ~360.

**Ce que la duplication des nœuds change, et c'est décisif.** L'objection
standard au SMS est que $\tilde{\mathbf M}$ devient non diagonale, donc qu'il
faut résoudre un système à chaque pas. **En FDEM à nœuds dupliqués, aucun nœud
n'est partagé entre deux éléments** : la matrice assemblée est donc
**bloc-diagonale par élément**, un bloc $4\times4$ par direction, et ces blocs
sont **constants**. On les inverse **une fois à l'initialisation**, en forme
fermée. Le pas de temps coûte alors un produit $4\times4$ par élément — quelques
pourcents — et **le solveur reste matrix-free**. La transposition au tétraèdre
est immédiate :
$\boldsymbol\lambda_{4\times4} = \frac{\beta m^e}{12}[\,3\ \text{diag},\ -1\,]$,
qui double la masse nodale à $\beta = 1$ et annule exactement la translation
rigide (somme des lignes nulle).

> **Réserve à porter au dossier.** Tkachuk & Bischoff notent que la méthode
> algébrique *« changes eigenmodes and increases rotational inertia of the
> element »*. rockim est un code de **fragmentation** : ses débris tournent, et
> la formulation co-rotationnelle est là pour ça. Alourdir l'inertie de rotation
> des éléments est donc un effet à **mesurer**, pas à supposer négligeable — sur
> la vitesse de rotation des fragments détachés et sur la restitution.
> Ils relèvent aussi que *« the factor defining the amount of SMS has no clear
> meaning »* : $\beta$ devra être calibré, pas choisi.

**Proposition, en deux niveaux.**

```
massScale = none | conventional | selective     # defaut none, bit-identique
massScaleBeta   = 0.0        # selective : le beta d Olovsson (methode II)
massScaleMinDt  = 0.0        # conventional : cible de pas
massScaleMaxAdd = 0.01       # refus au-dela de 1 % de masse ajoutee
```

Le mode `conventional` est le repli simple (une demi-journée) ; le mode
`selective` est la vraie solution (deux à trois jours, l'inversion des blocs
$4\times4$ à l'init et un produit par élément dans `integrate()`).

⚠️ **En insertion adaptative**, les groupes liés s'intègrent comme un nœud
unique : les blocs ne sont plus par élément mais **par groupe**. Ils restent
locaux et constants tant que le groupe existe, et se réinversent exactement là
où `rebindVertex()` recalcule déjà les groupes. C'est le seul point délicat.

**Garde-fous non négociables** : audit imprimé (masse ajoutée en %, éléments
touchés, dix pires), refus au-delà de `massScaleMaxAdd` (erreur, pas
avertissement), et l'énergie cinétique supplémentaire portée au bilan.

**Critères de falsification.** Sur `indent3d_grad` : (i) pic de force et compte
de fissures à moins de 2 % de la référence ; (ii) temps mur divisé par ≥ 2 ;
(iii) **vitesse de rotation moyenne des fragments détachés** à moins de 5 % —
c'est le test qui vise la réserve de Tkachuk ; (iv) translation rigide
strictement inaltérée (test dédié : bloc lancé sans force, la vitesse ne doit
pas bouger d'un bit).

### C2. Checkpoint / restart

Sérialiser `X0_, u_, v_, f_, m_, flag_`, l'état des joints (D, slip, omax,
smax, dead, bonded, dn0, propriétés tirées), l'état des lois (`MatState`),
l'état de l'outil, les compteurs d'énergie et le RNG. Format binaire simple,
un fichier par run, écrit toutes les N frames. Clés `checkpointEvery`,
`restartFrom`.

**Critère** : un run coupé à mi-parcours et redémarré doit donner une trajectoire
**bit-identique** à un run continu. C'est le seul test qui vaille.

### C3. VTU binaire

Une frame de 842 k tétraèdres pèse **~360 Mo en ASCII**. Encodage base64 en
`appended` avec `header_type="UInt64"` : facteur ~8 en taille, facteur ~20 en
temps d'écriture. Clé `vtkBinary` (défaut true après validation, ParaView lit
les deux).

### C4. `MatState` par loi

700 octets par élément quelle que soit la loi, alors qu'un bulk élastique n'en
a besoin d'aucun. Union discriminée ou allocation par la loi. Gain ×4-5 : c'est
ce qui a tué le cas 842 k tets sur la machine à 7 Go.

---

## D. PHYSIQUE DU VOLUME

### D1. L'indentation exige une loi AVEC CAP — c'est un choix de config, pas du code

**Le manque, mesuré le 18/08 puis confirmé le 19/08.** Sous un bouton
sphérique, le confinement est triaxial et le cône de Mohr-Coulomb **s'ouvre
sans borne** : la pression de contact atteint **1,75 GPa** et τ_max 1 700 MPa
sur une roche à 235 MPa de résistance en compression. La loi fait ce qu'on lui
demande, mais ce n'est pas physique. Le run `indent3d_grad` du 19/08 le
confirme : σ₁ = −150 MPa dans un noyau, von Mises 250 MPa, **aucune fissure
médiane ni radiale**, un simple disque broyé de 1,5 mm de profondeur.

**Proposition.** Rejouer `indent3d_grad` à l'identique avec `law = saksala`
(Perzyna + cap en pression durcissant) puis `law = dpdfh` (la loi calibrée de la
thèse, avec son cap et son endommagement anisotrope). Aucun développement : les
deux lois sont portées et vérifiées contre leur référence Fortran.

**Ce qu'on attend, et qui falsifierait.** Avec un cap, la pression de contact
doit plafonner à l'ordre de la résistance en compression (quelques centaines de
MPa, pas 1,75 GPa), et le faciès doit faire apparaître le cône hertzien et les
fissures médianes que la littérature d'indentation décrit. Si le faciès reste
un disque plat, le problème n'est pas dans la loi de volume mais dans la loi de
joint ou dans le contact.

### D2. Tétraèdres composites non bloquants

**Le manque.** Le tétraèdre linéaire à déformation constante **verrouille** en
quasi-incompressible — et la zone broyée sous l'outil, plastifiée et confinée,
est exactement ce régime. La raideur y est donc surestimée par construction,
ce qui fausse la force de contact et l'énergie stockée.

**Littérature — et la solution la moins chère n'est pas celle que je croyais.**
*Fiche révisée le 2026-08-19 après lecture de Bonet & Burton (CNME 14, 1998,
437–449), transmis par F. Uzquiano, et récupération du manuscrit COMPTet.*

**Option retenue : le tétraèdre à PRESSION NODALE MOYENNÉE (Bonet & Burton).**
L'idée tient en quatre équations et les auteurs annoncent qu'elle
*« can be implemented into any existing dynamic code with very minimal
additional expense »* — ce qui est vrai :

$$V_a = \sum_{e\ni a}\tfrac14 V^{(e)},\qquad v_a = \sum_{e\ni a}\tfrac14 v^{(e)},
\qquad J_a = \frac{v_a}{V_a}$$
$$p_a = \kappa\,\frac{v_a - V_a}{V_a},\qquad
\bar p^{(e)} = \tfrac14\sum_{a=1}^{4} p_a,\qquad
\mathbf T^{(e)}_{\text{vol},a} = \bar p^{(e)}\,v^{(e)}\,\nabla N_a^{(e)}$$

Le diagnostic qu'ils posent est exactement le nôtre : le tétraèdre standard
impose la conservation du volume **par élément**, soit $m$ contraintes pour
$3n$ degrés de liberté, et comme $m/n > 5$ en maillage tétraédrique, le
verrouillage est inévitable. En moyennant la pression **par nœud**, on ne
pose plus que $n$ contraintes : le verrouillage disparaît.

> **Ce qui rend cette option particulièrement adaptée à rockim** : la
> « pression nodale » est une somme sur les éléments **partageant le nœud**. En
> FDEM les nœuds sont dupliqués, donc la question devient : *sur quoi
> assemble-t-on ?* La réponse est déjà écrite dans le code — sur les
> **groupes liés** de l'union-find (`grpsOfVert_`, `copiesOfVert_`). Le volume
> nodal se scinde alors **exactement quand le nœud se scinde**, c'est-à-dire
> quand un joint s'insère. La machinerie existe, il n'y a rien à inventer.
>
> Coût : une passe supplémentaire sur les éléments (accumulation des $v^{(e)}$
> par groupe) avant la passe de forces. Estimation : +15 à 25 % sur
> `elementForces()`, à comparer au facteur ~3 de dt qu'on gagnerait ailleurs.

**Option de référence, plus lourde : COMPTet** — *A non-locking composite
tetrahedron element for the combined finite discrete element method*, Lei,
Rougier et al., **Engineering Computations 33(7) 2016, 1929–1956**,
DOI 10.1108/EC-09-2015-0268. Élément composite à 10 nœuds assemblé à partir de
**huit tétraèdres linéaires**, intégration **sélective** (réduite sur la partie
volumique, complète sur la déviatorique), sur la décomposition multiplicative de
Munjiza. C'est la solution du code **HOSS** (LANL), donc écrite *pour la FDEM*.
**Le manuscrit accepté est en accès libre sur le dépôt LANL** et est désormais
classé dans la bibliographie — il n'est plus à chercher.

**Recommandation révisée.** Commencer par **Bonet & Burton**, qui ne change ni
l'élément, ni le maillage, ni les joints, ni le pas de temps, et qui se greffe
sur l'union-find existant. COMPTet reste la référence FDEM-native, mais il
change l'élément — donc le maillage, les joints, la masse et le pas — et ne se
justifie que si l'ANP échoue.

**Critères de falisification (ANP).** (i) une barre élasto-plastique en
compression à $\nu \to 0{,}5$ plastique doit cesser de sur-raidir ; (ii) sur
`indent3d_grad`, la pression de contact doit **baisser** si le verrouillage
contribuait aux 1,75 GPa ; (iii) contrôle à charge nulle inchangé ; (iv) veiller
au damier de pression, mode parasite connu des tétraèdres à pression nodale —
si des oscillations spatiales apparaissent, passer aux formulations stabilisées.

**Réserve honnête.** L'ANP vise le **quasi-incompressible**. Le granite élastique
($\nu = 0{,}2$) ne l'est pas ; c'est la **zone plastifiée** sous l'outil qui
l'est (écoulement Mohr-Coulomb à faible dilatance, $\psi = 5^\circ$). L'effet
attendu est donc **local à la zone broyée** — ce qui est précisément là où la
force de contact se joue, mais il ne faut pas en attendre un changement global.

---

## E. LE LOT DE DETTE SILENCIEUSE — une demi-journée, huit points

Tous partagent le même mode de panne : **le solveur accepte un réglage qui ne
fait rien**, ou observe un état qu'il ne peut pas voir.

| # | défaut (vérifié dans les sources le 19/08) | correctif |
|---|---|---|
| E1 | `ftScale` n'agit que dans `dpdfh` et `saksala2011` : sous `dpr` et `saksala` le champ de Weibull est **écrit dans les VTU sans aucun effet** | appliquer le facteur, ou **refuser la config** |
| E2 | `toolShape` est validé puis **écrasé par `DISC`** en percussion 2D : un poinçon plat et un disque donnent des runs bit-identiques | honorer la clé (elle l'est déjà en 3D) |
| E3 | `chamferLen` / `chamferDeg` lus, stockés, **jamais employés** | implémenter le chanfrein ou refuser la clé |
| E4 | `pairKey(i,j) = (i<<32)|j` **non trié** en `dem3d` : perte silencieuse de l'historique tangentiel | trier les indices |
| E5 | garde NaN sur `u_[0]`, qui peut être un nœud FIXED donc toujours fini : **détecteur aveugle** | tester une norme globale échantillonnée |
| E6 | la trace de l'écrêtage affiche « 0 m/s » (`tool_.v` pas encore assigné) : **le journal ment** | déplacer l'impression |
| E7 | `detachedVol` compte l'insert comme fragment détaché (5 446 mm³) | exclure les groupes physiques outil |
| E8 | `difT` / `edotIns` **non exportés** : impossible d'auditer la population insérée, là où un attracteur s'était logé | ajouter aux VTU de joints |

**Règle commune** : un réglage qui ne fait rien doit **lever une erreur**, pas
être ignoré en silence.

---

## F. CE QUI N'EST PAS DU CODE

1. **Reproduire un résultat de Yan** (UCS ou brésilien de l'article 2023) et
   **un de Yang** (un des sept critères d'impact) : les deux capacités sont
   implémentées et vérifiées en interne, **aucune n'est validée** contre son
   article.
2. **Recalibration Red Bohus de zéro**, avec `jointWeibullM` cette fois —
   l'ingrédient manquant identifié après l'échec de l'émulateur GP.
3. **Données expérimentales d'Aising** (Mines Paris) : bloquant pour toute
   validation d'impact, et personne ne peut le coder.

---

## G. ORDRE D'EXÉCUTION PROPOSÉ

*Ordre révisé le 2026-08-19 : C1 s'alourdit (mais rapporte bien plus) et D2
s'allège nettement, une fois les trois articles lus.*

```
   nuit 1        B1  balayage dampingLocal (3 runs, AUCUN code)
   jour 1        E   lot de dette silencieuse + B2 observables anti-pompe
   jour 2-3      A1  contact outil en condition de vitesse   <-- LE correctif
   jour 3        D1  rejouer indent3d_grad en saksala puis dpdfh (aucun code)
   jour 4-5      D2  tetraedre a pression nodale moyennee (Bonet & Burton)
   jour 6-8      C1  mass scaling selectif d Olovsson, methode II
   semaine 2     C2  checkpoint/restart, C3 VTU binaire, C4 MatState
   ensuite       A3  cap de raideur cohesive
   fond          A2  contact non lisse (NSN)   |  COMPTet si l ANP echoue
```

D2 remonte avant C1 pour une raison de méthode : **les deux touchent la zone
broyée**, et si le mass scaling est posé d'abord, on ne saura plus attribuer un
changement de pression de contact au verrouillage ou à l'inertie ajoutée. Une
correction à la fois, mesurée contre le même cas de référence.

La logique : **B1 et E ne coûtent presque rien et rendent tout le reste
lisible**. A1 est le seul correctif qui débloque une comparaison à l'expérience.
D1 ne demande aucun développement et tranchera si le faciès plat vient de la loi
de volume. A2 et D2 sont des chantiers de fond qu'il serait imprudent de lancer
avant d'avoir mesuré l'effet des précédents.

---

## H. ARTICLES

**Accès libre, déjà consultés :**
- Ghesquière-Diérickx, Molinari & Anciaux, *Stability of Extrinsic Cohesive-Zone
  Model with Penalty-Based Contact...*, arXiv:2511.14323.
- Ghesquière-Diérickx, Anciaux, Acary & Molinari, *The semi-explicit nonsmooth
  Newmark time integrator...*, arXiv:2606.01355.
- Dureisseix et al., *Explicit dynamics and non-smooth interface behaviors*,
  JTCAM 2024 (algorithme CD-Lagrange complet).
- Fekak, Brun, Gravouil et al. (CD-Lagrange original), HAL hal-01852099.

**Reçus de F. Uzquiano le 2026-08-19, lus et intégrés aux fiches C1 et D2 :**
- **Olovsson, Simonsson & Unosson**, *Selective mass scaling for explicit finite
  element analyses*, IJNME **63** (2005) 1436–1445 → fiche C1 réécrite
  (méthodes I et II, contrainte $\lambda\ddot u_r = 0$).
- **Tkachuk & Bischoff**, *Variational methods for selective mass scaling*,
  Comput. Mech. **52** (2013) 563–570 → réserve sur l'inertie de rotation et
  sur l'absence de sens physique du facteur.
- **Bonet & Burton**, *A simple average nodal pressure tetrahedral element…*,
  CNME **14** (1998) 437–449 → fiche D2 réécrite, l'ANP devient l'option n° 1.

**Récupéré par mes soins le même jour :**
- **Lei, Rougier et al.**, COMPTet, Engineering Computations 33(7) 2016 —
  manuscrit accepté **en accès libre sur le dépôt LANL**, classé dans
  `Documents\bibliographie\` sous
  `Lei_Rougier_2016_COMPTet_non-locking_composite_tetrahedron_FDEM_EngComput33.pdf`.

**Reste introuvable, et ce n'est pas bloquant :**
- **Cundall (1987)** sur l'amortissement local. La règle dont j'ai besoin —
  *l'amortissement local sert à la mise à l'équilibre et au quasi-statique ; en
  dynamique la dissipation doit venir de l'amortissement visqueux du modèle de
  contact* — est reprise à l'identique dans la documentation Itasca et dans la
  littérature FDEM récente. On citera une reprise, ou on s'en passe : le
  balayage B1 tranchera par la mesure, pas par l'autorité.

**Aucun article ne conditionne plus le démarrage.** A1, B1, B2, C1, D1, D2 et E
sont spécifiés avec ce qui est en main.
