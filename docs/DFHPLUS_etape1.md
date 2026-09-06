# DFH+ étape 1 — la loi `dfhplus` : le périmètre de DP-DFH sur un cadre neuf

**2026-09-06, branche `g0`.** Compte rendu chiffré de la première étape du chantier DFH+ :
reconstruire la loi de la thèse sur une **énergie libre postulée unique**, sans encore ajouter
aucun mécanisme nouveau, et **prouver que le cadre tient** avant d'y greffer quoi que ce soit.

Fondement : `CONTINUUM/loi_dfh_plus/loi_DFH_plus.pdf` et `.tex` (66 p., 146 équations, deux
relectures adverses). Instrument : `rockim thermobench`, écrit la veille, dont les deux contrôles
sont `elastic` (doit passer) et `dpdfh` (doit échouer).

---

## 1. Le résultat central, en une ligne

| | `dpdfh` | `dfhplus` |
|---|---|---|
| `rockim thermobench <loi>` | **ÉCHEC** (7 172 violations) | **PASS** (0 violation) |
| symétrie majeure de la tangente élastique | 0 / 3 574 | **0 / 2 078** |
| dissipation ≥ 0 | **7 172 / 59 483 (12,1 %)** | **0 / 120 000** |
| dont violations significatives (> 1 % de Gf/ℓc) | **1 124** | **0** |
| pire déficit absolu | **−15 240 J/m³ = 17,6 % de Gf/ℓc** | **−4,67 × 10⁻³ J/m³ = 3,8 × 10⁻⁶ % de Gf/ℓc** |
| réduction élastique | 0 / 96 000 (7,7e−16) | 0 / 96 000 (3,3e−15) |
| objectivité (invariants) | 0 / 480 000 (5,1e−14) | 0 / 480 000 (2,3e−14) |
| continuité σ(ε), à temps gelé | 0 / 120 000 | 0 / 120 000 |
| incréments EXCLUS du verdict de dissipation | **60 517 / 120 000** | **0** |

Mêmes 12 000 tirages, même graine 20260906, mêmes six familles de chemins, même carte Red Bohus.
Le facteur sur le pire déficit absolu est **3,3 × 10⁶**.

La colonne « incréments exclus » est aussi importante que les autres : `dpdfh` n'expose pas son
énergie libre, le banc devait l'estimer par sonde de décharge et écarter la moitié des incréments
(contaminés ou non relâchés). `dfhplus` **expose** son énergie libre (§4) : le test 2 devient
**exact et complet**, sans estimateur ni exclusion.

### Le même verdict avec le MÊME instrument

Une loi qui expose son énergie libre est jugée par un instrument différent : la comparaison
ci-dessus n'est donc pas tout à fait à armes égales. `rockim thermobench dfhplus --probe` force
l'estimateur par sonde de décharge, celui de `dpdfh` :

| avec la SONDE DE DÉCHARGE (instrument identique) | `dpdfh` | `dfhplus --probe` |
|---|---|---|
| dissipation < 0 | **7 172 / 59 483** | **1 / 79 826** |
| dont significatives | **1 124** | **0** |
| pire déficit absolu | −15 240 J/m³ (17,6 % Gf/ℓc) | −403 J/m³ (0,63 % Gf/ℓc) |
| incréments contaminés | 49 650 | 40 089 |
| incréments non relâchés | 10 867 | **85** |

7 172 → 1, et la seule violation résiduelle reste *sous* le seuil de significativité. La chute des
« non relâchés » de 10 867 à 85 dit la même chose autrement : la contrainte de `dfhplus` reste
alignée sur la compliance vierge pendant une décharge, parce qu'elle dérive d'un potentiel.

---

## 2. L'énergie libre postulée

Notation de `loi_DFH_plus.tex`. `\epse` = déformation élastique, `D_i` = les trois endommagements
de traction portés par le repère **figé** {n₁, n₂, n₃}, λ et G les coefficients de Lamé.

**Décomposition spectrale** (Miehe et al. 2010) de la déformation élastique :

```
\epse = Σ_a ε_a  m_a ⊗ m_a ,      \epse^± = Σ_a <ε_a>_±  m_a ⊗ m_a
```

{m_a} est le repère **principal courant** de `\epse`. Il tourne avec le chargement et n'a aucune
raison de coïncider avec le repère figé : c'est exactement là que DP-DFH perdait la positivité.

**Tenseur d'intégrité** porté par le repère figé, et **intégrité volumique** (moyenne harmonique) :

```
A   = Σ_i (1 − D_i)  n_i ⊗ n_i ,      g(m) = m · A m = Σ_i (1 − D_i) (n_i·m)²
g_v = 3 / [ 1/(1−D_1) + 1/(1−D_2) + 1/(1−D_3) ]
```

```
┌──────────────────────────────────────────────────────────────┐
│  ρψ(\epse, D₁, D₂, D₃, {n_i}) =                              │
│      (λ/2) [ g_v <tr \epse>₊² + <tr \epse>₋² ]                │
│    + G  A : (\epse⁺)²                                        │
│    + G  ‖\epse⁻‖²                                            │
└──────────────────────────────────────────────────────────────┘
```

La deuxième ligne se lit aussi `G Σ_a g(m_a) <ε_a>₊²` (identité `A : (\epse⁺)² = Σ_a g(m_a)
<ε_a>₊²`). Écrite sous la forme `A : (\epse⁺)²`, elle est manifestement objective et se dérive
**sans jamais dériver les vecteurs propres**.

**Pourquoi la moyenne harmonique pour g_v.** C'est le couplage en *série* : une dilatation isotrope
doit être transmise à travers les trois familles de facettes, les compliances s'ajoutent, et le
ressort casse dès qu'un maillon casse. C'est la seule moyenne qui rende à la fois `g_v = 1−D` pour
`D₁=D₂=D₃=D` (réduction au cas isotrope) **et** `g_v → 0` dès qu'un seul `D_i → 1` — sans quoi la
traction uniaxiale saturée garderait une raideur résiduelle que la VUMAT n'a pas.

**Trois propriétés par construction :**

1. **Réduction élastique exacte.** À D = 0 : `g_v = 1`, `A = I`, les deux moitiés spectrales se
   recollent, `ρψ = (λ/2)(tr \epse)² + G \epse:\epse`. Mesuré : **3,3e−15** sur 96 000 états.
2. **Unilatéralité naturelle.** Seule la partie positive est dégradée. Aucun interrupteur, aucun
   test de signe de contrainte.
3. **Symétrie majeure de la tangente élastique.** σ est un gradient, la tangente est une hessienne.
   Mesuré : **7,3e−11** — le bruit machine des différences finies.

---

## 3. Ce qu'on en dérive

### 3.1 La contrainte

```
σ = λ [ g_v <tr \epse>₊ + <tr \epse>₋ ] I
  + P⁺ : [ G (A \epse⁺ + \epse⁺ A) ]
  + 2 G \epse⁻
```

`P⁺ = ∂\epse⁺/∂\epse` est la projection spectrale positive (Daleckii–Krein) : dans le repère propre
de `\epse`, `(P⁺:T)_aa = H(ε_a) T_aa` et `(P⁺:T)_ab = θ_ab T_ab` avec
`θ_ab = (<ε_a>₊ − <ε_b>₊)/(ε_a − ε_b) ∈ [0,1]`, prolongé par H à la dégénérescence. `P⁺` a la
symétrie majeure (c'est la hessienne de `\epse ↦ ½‖\epse⁺‖²`), ce qui autorise cette écriture.

**Vérifiée contre les différences finies** sur ρψ, 6 états × 6 composantes, banc (1) du selftest :
**écart relatif 2,0 × 10⁻⁹** (pas h = 1e−10, précision de la DF centrée).

### 3.2 Les forces motrices

```
Y_i = (λ/2) [ g_v² / (3 (1−D_i)²) ] <tr \epse>₊²  +  G ‖\epse⁺ n_i‖²   ≥ 0
```

**Positivité inconditionnelle** : deux carrés multipliés par λ > 0 et G > 0. Comme la cinétique
d'obscuration est monotone (`dD_i ≥ 0`), la dissipation d'endommagement `Σ_i Y_i dD_i` est positive
**sans aucune condition** — c'est ce que DP-DFH ne pouvait pas tenir (§T:B de
`loi_DFH_plus.tex` : contre-exemple en compression uniaxiale, `Y < 0` dès l'amorçage).

**Borne** : `g_v ≤ 3(1−D_i)` toujours, donc `g_v²/(1−D_i)² ≤ 9`. Y_i reste **borné** quand D_i → 1,
contrairement au `Y_i^o ~ (1−D)^{−1/2}` de la variante R du document (facteur 7 à D = 0,98).

Vérifiées analytique-contre-différences-finies, banc (2) : écart **0** (les états testés ont soit
D = 0, soit D plafonné ; les cas intermédiaires passent par le potentiel exposé) ; **min Y_i = 0**,
jamais négatif.

### 3.3 La plasticité, et sur quelle contrainte

La contrainte effective du cadre est **celle qui rend le squelette sain** :

```
\bar\sigma := ∂(ρψ)/∂\epse |_{D=0} = \Cel : \epse
```

Trois raisons, et une conséquence :

* c'est une **dérivée de la même énergie libre**, prise à endommagement nul : elle n'est pas
  postulée à côté du potentiel, elle en sort ;
* c'est la **seule mesure de contrainte indépendante de D**, donc la seule qui donne une surface de
  charge dont la position ne bouge pas quand la roche se fissure. Physiquement : le frottement se
  joue sur les contacts intacts, pas sur l'aire perdue ;
* elle est **linéaire en `\epse`**, donc le critère `f = q − p tan β − d`, le retour radial et
  l'apex de `dpdfh` s'y transposent **sans aucune modification**.

**Conséquence mesurée** : tant qu'aucune direction n'est amorcée, `dfhplus` et `dpdfh` donnent la
**même** réponse plastique. Compression uniaxiale et triaxiaux, pic de q :

| σ₃ | `dfhplus` | `dpdfh` | écart |
|---|---|---|---|
| 0 MPa | 265,259 MPa | 265,259 MPa | 2,4e−10 |
| 20 MPa | 309,078 MPa | 309,078 MPa | 1,7e−10 |
| 50 MPa | 374,807 MPa | 374,807 MPa | 2,5e−10 |
| 100 MPa | 484,355 MPa | 484,355 MPa | 2,6e−10 |

### 3.4 L'écrêtage de dilatance (clé `dfhpPsiClamp`, défaut actif)

Avec une énergie libre unique, la dissipation vaut `D = σ^nom : d\epsp + Σ_i Y_i dD_i` : c'est la
contrainte **nominale**, pas l'effective, qui est appariée à `d\epsp`. Avec l'écoulement DP non
associé `d\epsp = dγ [ (3/(2\bar q)) \bar s + (tan Ψ / 3) I ]` :

```
D_A / dγ = q_proj − tan(Ψ) p^nom ,
q_proj := (3/(2\bar q)) \bar s : s^nom ,   p^nom := −tr(σ^nom)/3
```

qui devient négatif dès que l'endommagement est fortement anisotrope. On impose donc
`tan Ψ ≤ max(q_proj, 0)/p^nom` quand `p^nom > 0` : c'est l'écrêtage d'admissibilité
`eq:T:clamp` / `eq:D1:psibound` de `loi_DFH_plus.tex`, transposé au cône, avec un
**prédicteur-correcteur** (la contrainte nominale MOYENNE du pas, pas celle du début, faute de quoi
D_A est annulée au premier ordre et rendue négative par le second).

---

## 4. L'ajout au banc : l'énergie libre exposée

Deux méthodes virtuelles sont **ajoutées** à `MatLaw` (croissance par addition, principe VIII —
défaut « non exposée », zéro effet sur les lois existantes, `dpdfh` comprise) :

```cpp
virtual bool   hasFreeEnergy() const                          { return false; }
virtual double freeEnergy(const Matrix3d&, const MatState&)   { return 0.0; }
virtual int    damageForces(const Matrix3d&, const MatState&, double*) { return 0; }
```

Quand la loi les expose, le banc **bascule automatiquement** de l'estimateur par sonde à la valeur
exacte. C'est le remède annoncé dans les *issues* du banc, et il supprime d'un coup les trois
limites L1 (part récupérable seule), L2 (incréments contaminés) et L3 (incréments non relâchés).

**Deux discriminants ont dû être ajoutés au banc** pour que ses verdicts veuillent dire quelque
chose sur une loi C¹ par morceaux et pilotée par le temps. Les deux sont documentés dans
`src/ThermoBench.cpp` et publiés dans le compte rendu.

**(a) Quadrature du test 2.** Le trapèze commet, à chaque *coin* de la réponse (une valeur propre de
ε qui change de signe alors que `g_v < 1`, l'endommagement qui plafonne), une erreur en O(h²) que
rien ne distingue d'une dissipation négative. Un incrément flagué à 8 sous-pas est donc rejoué à 32 :
l'erreur de quadrature est divisée par 16, une vraie violation ne bouge pas. Sur `dfhplus` :
**58 incréments flagués à 8 sous-pas se sont tous effacés à 32** ; pire cas relatif au trapèze
grossier −8,1e−2, à 8 sous-pas −1,2e−3, à 32 sous-pas sous la tolérance. Le rapport 65 ≈ 8² entre
grossier et fin est la signature exacte de la quadrature.

**(b) Temps gelé du test 5.** `r = ‖dσ‖/(M‖dε‖)` rapporte une variation de *contrainte* à un
incrément de *déformation*. Un mécanisme piloté par le **temps** — l'obscuration de Denoual–Hild,
`dx = (S λ)^{1/3} k c dt` — fait tomber σ **sans que ε bouge** : r y est arbitrairement grand *et*
invariant au raffinement (subdiviser divise dt d'autant), et pourtant σ(ε) n'a pas la moindre
discontinuité. Le banc rejoue donc chaque incrément à temps gelé (dt × 1e−6) et fait porter le
verdict sur cette valeur-là.

> **Ceci corrige une lecture du compte rendu du banc.** Il y était écrit que les 84 flags de
> `dpdfh` étaient « une discontinuité de la loi à saturation de l'endommagement en traction, pas un
> transitoire de pas de temps », au motif que le raffinement ne les effaçait pas. **C'est faux** :
> le raffinement ne pouvait pas les effacer, puisqu'il divise dt en même temps que dε. À temps
> gelé, **les 84 flags de `dpdfh` et les 87 de `dfhplus` disparaissent tous**, et le pire r tombe à
> **0,99999 ≤ 1** — exactement la borne élastique — pour les deux lois. Les deux réponses σ(ε) sont
> continues ; ce qui saute, c'est la relaxation temporelle, et c'est la physique voulue.

Le contrôle `elastic` reste PASS (0 violation sur les cinq tests) après ces deux ajouts, et
`dpdfh` reste ÉCHEC pour la seule raison qui compte : 7 172 dissipations négatives, dont 1 124
significatives.

---

## 5. Ce que `dfhplus` change par rapport à `dpdfh`, chiffré

Banc (6) de `rockim selftest-dfhplus`, carte Red Bohus, ℓc = 1 mm.

### Traction uniaxiale (pilotage libre latéralement)

| grandeur | `dfhplus` | `dpdfh` | écart |
|---|---|---|---|
| pic | 112,529 MPa | 112,529 MPa | **−5,3e−13 %** |
| ε au pic | 1,449e−3 | 1,449e−3 | identique |
| aire sous la courbe | 0,082148 MJ/m³ | 0,081961 MJ/m³ | **+0,23 %** |
| σ à D = 0,98 | **0,0269 MPa** | 1,094 MPa | ÷41 |

Le **pic est identique** parce qu'il est atteint à l'amorçage, où D = 0 et où les deux lois sont la
même élasticité : la contrainte équivalente d'énergie `σ_i^eq = √(2 E Y_i / c_ν)` est calibrée
(`c_ν = (2/(1+ν))[ν(1−2ν)/6 + ½]`, 0,8067 à ν = 0,29) pour valoir **exactement** `\bar\sigma_ii` en
traction uniaxiale libre à D = 0. L'écart apparaît sur la **branche adoucissante** :

* En cours d'endommagement, `dpdfh` rend `σ₁₁ = (1−D) \bar\sigma₁₁` ; `dfhplus` rend
  `σ₁₁/(E ε) = [ν g_v + (1−D)] / (1+ν)`, plus **raide** de `ν (g_v − (1−D))/(1+ν)`. Cet excès est
  nul à D = 0 et à D = 1, maximal vers D ≈ 0,65 où il vaut **+6,0 % de E en module sécant**, soit
  **+17 % relatifs** sur la contrainte. C'est le prix du couplage volumique : le terme λ ne peut pas
  être dégradé direction par direction sans casser la symétrie majeure.
* En fin de course, l'inverse : à D = 0,98 la moyenne harmonique a déjà tué `g_v` alors que la
  VUMAT garde `(1−D) \bar\sigma` sur une contrainte effective qui, elle, continue de croître — d'où
  1,09 MPa contre 0,027. `dfhplus` **relâche plus complètement**.
* Le bilan énergétique est à **+0,23 %** : les deux écarts se compensent presque exactement.

### Contrainte latérale résiduelle

Le découpage spectral laisse, en traction uniaxiale endommagée, `σ₂₂ = (Eν/(1+ν))(g_v − 1) ε < 0`
(artefact connu du *split* de Miehe, rigoureusement nul à D = 0). À D = 0,5 il vaut −5,6 % de E·ε.

### Cisaillement

La VUMAT **postule** un facteur `min(f_i, f_j)` sur les composantes de cisaillement. Ici le facteur
**tombe du potentiel** et vaut la moyenne **arithmétique** des intégrités des deux directions
(`A : (\epse⁺)²` évalué sur un cisaillement pur dans le plan (1,2) donne
`½[(1−D₁)+(1−D₂)]`). Ce n'est pas un choix, c'est une conséquence.

### Chemin non coaxial (traction jusqu'à saturation, puis rotation du repère + cisaillement)

| grandeur | `dfhplus` | `dpdfh` |
|---|---|---|
| ‖σ‖ max | 113,37 MPa | 116,47 MPa (−2,7 %) |
| travail total | 0,3114 MJ/m³ | 0,2767 MJ/m³ (+12,5 %) |
| D final | 0,9999 | 0,9999 |
| écrêtage de dilatance actif | 0 fois | — |

C'est la famille où DP-DFH concentrait les cinq pires déficits du banc (−15,2 / −15,0 / −9,9 /
−9,8 / −9,8 kJ/m³) : `dfhplus` y a **0 violation** et y dépense **12,5 % de travail en plus** — la
différence n'est pas cosmétique, c'est l'énergie que le cadre ancien se rendait à lui-même.

---

## 6. L'exposant de vitesse 3/(m+3) : la physique survit au changement de cadre

C'est le contrôle qui dit si l'obscuration de Denoual–Hild est encore la même physique une fois
branchée sur `Y_i` au lieu de la contrainte normale du repère figé. Traction uniaxiale pilotée en
déformation, 7 vitesses de 10³ à 10⁹ /s, régression de log(σ_pic) sur log(ε̇).

**La loi d'échelle ne vaut que dans le régime de fragmentation MULTIPLE**, `λ_frag V_el ≫ 1` : en
dessous, l'élément casse sur *un* défaut (plancher `xlam·Vel < 1` de la VUMAT) et le pic est
rate-indépendant. Le banc ne régresse donc que sur les points où `λ_frag V_el > 10`, soit
`σ_pic > σ_w 10^{1/m}`, et où la branche adoucissante a bien été atteinte.

| m | cible 3/(m+3) | `dfhplus` | écart | `dpdfh` | écart | points |
|---|---|---|---|---|---|---|
| 6 | 0,33333 | **0,33327** | −0,02 % | 0,33326 | −0,02 % | 4 |
| 12 | 0,20000 | **0,19890** | −0,55 % | 0,19882 | −0,59 % | 6 |
| 24 | 0,11111 | **0,11066** | −0,41 % | 0,11063 | −0,44 % | 6 |

**L'exposant est conservé au demi-pour-cent près, et `dfhplus` est aussi près de la cible que
`dpdfh`.** C'était attendu : `λ_frag ∝ (σ^eq)^m` avec `σ^eq ≡ \bar\sigma` sur ce chemin, et
l'intégrateur en racine cubique est repris mot pour mot.

---

## 7. Le contrôle qui doit échouer

`dfhpPsiClamp = false` désarme l'écrêtage d'admissibilité. Chemin : traction jusqu'à saturation
selon x, puis **compression uniaxiale dans un repère tourné de 40°** — la roche est fissurée
normalement à x et on l'écrase obliquement, le cône travaille avec un endommagement fortement non
coaxial. 371 incréments plastiques.

| variante | incréments à dissipation plastique < 0 | pire |
|---|---|---|
| carte de la thèse, Ψ = 15°, **écrêtage actif** | **0 / 371** | 0 |
| carte de la thèse, Ψ = 15°, écrêtage désarmé | **0 / 371** | 0 |
| écoulement **associé**, Ψ = β = 51,7°, écrêtage actif | 252 / 371 | **−386 J/m³** (0,39 % Gf/ℓc) |
| écoulement **associé**, Ψ = β = 51,7°, écrêtage désarmé | 252 / 371 | **−8 127 J/m³** (8,1 % Gf/ℓc) |

Lecture, en trois points et sans arrondir les angles :

1. **Sur la carte de la thèse l'écrêtage ne mord jamais** : sur la surface de charge
   `\bar q/\bar p > tan β = 1,266` alors que `tan Ψ = 0,268`. La dissipation plastique y est
   positive sur les 371 incréments, avec ou sans la clé. La clé est une **assurance**, pas un
   correctif ; c'est une bonne nouvelle et il faut la dire comme telle.
2. **En écoulement associé la clé sert, et le facteur est 21** sur le pire déficit (8 127 → 386
   J/m³). C'est la falsification demandée : on a exhibé un régime où désarmer la clé dégrade
   mesurablement l'admissibilité.
3. **Elle ne suffit pas en associé.** Les 252 incréments restants ont `q_proj < 0` : c'est le terme
   **déviatorique** de `σ^nom : d\epsp` qui est négatif, pas le terme de dilatance, et aucun
   écrêtage de Ψ ne peut le corriger. Le déficit résiduel (0,39 % de Gf/ℓc) reste **sous** le seuil
   de significativité du banc, mais le point est ouvert : voir §9.

---

## 8. Les clés

Les **neuf clés matériau de `dpdfh` sont reprises à l'identique** — `dfhBetaDeg`, `dfhDCoh`,
`dfhPsiDeg`, `dfhWeibullM`, `dfhSigW`, `dfhZeff`, `dfhK`, `dfhS`, `dfhDeld` — plus `dfhPsiVar`,
`dfhPsi0`, `dfhKPsi`, `dfhPsiMax`. **Les cartes sont interchangeables : `law = dpdfh` →
`law = dfhplus` suffit.** Le tirage de Weibull utilise le **même hachage spatial** que la VUMAT
(`kst_seed` translittéré, pas appelé — `dpdfh` ne doit être ni modifié ni rendu dépendant du
nouveau fichier) : pour un même x₀ les deux lois tirent exactement les mêmes seuils.

Deux clés **propres** au cadre neuf :

| clé | défaut | effet |
|---|---|---|
| `dfhpPsiClamp` | `true` | écrêtage d'admissibilité de la dilatance (§3.4). `false` = la variante falsifiante. |
| `dfhpVolInteg` | `harmonic` | forme de `g_v` : `harmonic` (couplage en série), `min` (borne inférieure, non lisse), `none` (`g_v = 1` : ablation qui isole le couplage volumique, celui que le contre-exemple de `loi_DFH_plus.tex` §T:B accuse). |

Les deux sont enregistrées dans `tools/keys_by_mode.json` / `KeysByMode.hpp` (catégorie `shared`).

---

## 9. Ce qui n'est PAS fait, et ce qui reste ouvert

1. **Les mécanismes B et C ne sont pas codés.** L'étape 1 pose le cadre sur le périmètre connu.
   L'ouverture axiale directionnelle et la pulvérisation avec plafond viendront s'**ajouter** au
   potentiel — c'est le sens même du travail : il n'y a plus qu'un objet à modifier.
2. **La dissipation plastique déviatorique en écoulement associé** (§7, point 3). Le cadre garantit
   `Σ Y_i dD_i ≥ 0` inconditionnellement, mais `σ^nom : d\epsp ≥ 0` n'est garanti que par
   l'écrêtage, qui ne couvre que la dilatance. Un traitement complet demanderait de projeter la
   direction d'écoulement — donc de renoncer à la parité exacte du retour DP avec `dpdfh`. À
   trancher par Fernando : la parité vaut-elle plus que l'admissibilité inconditionnelle en associé,
   sachant que la carte de la thèse est à Ψ = 15° et n'atteint jamais le régime ?
3. **Le test 3 en mode `--ref` n'a toujours pas d'objet exact.** `dfhplus` ne *doit pas* rendre
   `dpdfh` à 1e−12 : ce n'est pas la même loi, et §5 chiffre pourquoi. La réduction qui est vérifiée,
   et qui est la bonne, est la réduction **élastique** (0 / 96 000 à 3,3e−15). Si l'on veut un test
   `--ref` qui morde, il faudra le poser sur `dfhpVolInteg = none` contre une variante de référence,
   ou accepter que le critère soit un écart chiffré et non une identité.
4. **La convexité de ρψ n'est pas établie.** `A : (\epse⁺)²` avec A anisotrope n'est pas
   manifestement convexe en `\epse`. Le banc ne teste ni la symétrie mineure, ni la convexité, ni la
   normalité, ni la stabilité de Drucker. Ces tests 6, 7… s'ajouteront au même banc, et ils
   dépendent de la forme retenue pour ρψ — c'est-à-dire maintenant qu'ils peuvent être écrits.
5. **La population « élastique » du test 1 est plus petite pour `dfhplus`** (2 078 contre 3 574) :
   la signature d'état irréversible du banc voit `wDamT`, que `dfhplus` alimente et `dpdfh` non.
   0 violation sur 2 078 à 7,3e−11 reste concluant, mais le chiffre n'est pas comparable tel quel.
6. **Aucun run de structure.** `dfhplus` n'a jamais tourné sur un maillage. La bit-identité 8/8 le
   confirme d'ailleurs : aucun deck existant ne l'utilise. Un banc à l'échelle de l'élément fini
   (objectivité au maillage, énergie de bande) reste à faire avant tout usage en percussion.

---

## 10. Reproduire

```
powershell -ExecutionPolicy Bypass -File tools\build.ps1     # build/rockim.exe
set OMP_NUM_THREADS=4
build\rockim.exe selftest-dfhplus  out_dfhplus.csv           # bancs (1)-(7), ~1 s
build\rockim.exe thermobench dfhplus                          # PASS, 0 violation
build\rockim.exe thermobench dfhplus x.csv --probe            # meme instrument que dpdfh
build\rockim.exe thermobench dpdfh                            # ECHEC, 7172
build\rockim.exe thermobench elastic                          # PASS (controle)
python tools\bitid.py --exe build\rockim.exe --threads 4      # 8/8 IDENTIQUE
```

Fichiers : `src/MatLawDfhPlus.cpp` (la loi, avec l'énergie libre en en-tête),
`include/rockim/MatLawDfhPlus.hpp`, `src/DfhPlusBench.cpp` (les sept bancs),
`include/rockim/MatLaw.hpp` (sous-état `MatState::Dfhp`, virtuelles `hasFreeEnergy` /
`freeEnergy` / `damageForces`), `src/ThermoBench.cpp` (chemin exact + les deux discriminants),
`src/MatLaw.cpp` (une branche `else if (kind == "dfhplus")` dans la fabrique — **`DpDfhLaw` n'est
pas touchée**).
