# L'insertion des joints dans rockim : critère, loi, objectivité

Note de fond, 7 octobre 2026. Lecture seule de `src/` et `include/` : aucun fichier de code n'a été modifié.
Les numéros de ligne renvoient à `src/FdemSolver.cpp` (2D) et `src/Fdem3dSolver.cpp` (3D) dans l'état du
dépôt à cette date, c'est-à-dire après la correction du contact du 7 octobre (`contactCandidates`,
`gcBirth = offset`, `potForceExact` présents dans le source). La référence principale pour la comparaison
est Yan, Zheng et Wang, IJRMMS 169 (2023) 105439 (PDF fourni par le doctorant, lu pages 3 à 8, et texte
extrait). Figure 1 (lois comparées, en ouverture) : `figures/insertion/fig_lois_insertion.pdf` (script `fig_lois_insertion.py`).

Conventions de la note : traction positive, σ_n > 0 en ouverture ; p_j raideur de joint par unité d'aire
(Pa/m) ; « adaptatif » désigne `insertion = adaptive`, « intrinsèque » `insertion = intrinsic`.

---

## 0. Résumé

Le critère de rockim est celui de Yan et al. (éq. 7-8) à la lettre : moyenne arithmétique des contraintes
des deux triangles voisins, projetée sur l'arête en configuration courante, testée à chaque pas, après le
calcul des contraintes et avant celui des forces de joint, sans plafond du nombre d'insertions. La loi
après insertion reproduit celle de Yan (courbe en z de Munjiza, endommagement elliptique, énergie G_I
exacte en traction pure) à trois différences près, qui ne sont pas des détails : rockim garde une branche
élastique de raideur `insertionPenaltyFactor` · E/h (4 E/h par défaut) sur laquelle le joint naît décalé
de δ_n0 = min(σ, f_t)/p_j ; la composante qui n'a pas déclenché l'insertion reste sur cette branche au lieu
d'être figée à sa valeur d'insertion multipliée par f(D) comme chez Yan ; la décharge part vers une origine
décalée de −δ_n0, ce qui borne la raideur sécante à p_j alors que celle de Yan est non bornée.

Avant insertion, le continu est exact (nœuds partagés), et la résistance apparente de l'adaptatif dépasse
celle de l'intrinsèque de 13 à 16 % sur les essais de Yan. La cause n'est pas établie ; deux mécanismes
sont plausibles et testables : la traction intrinsèque est une traction nodale locale, la contrainte
adaptative une moyenne sur deux triangles qui écrête les concentrations ; et les interfaces intrinsèques,
partout présentes, ajoutent une souplesse (E_eff = E/(1 + α/p_f), α = 1,239) qui localise. Le critère est
en contrainte seule, sans condition d'énergie ni de résolution de la zone cohésive : la propagation d'une
pointe dépend donc de h quand la zone cohésive n'est résolue que par 2 à 3 éléments, ce qui est le cas du
tunnel en mode I.

Le facteur de pointe ne relâche que le déclenchement ; le joint né à f_t/k doit ensuite remonter jusqu'à
f_t sur un ressort de 4 E/h. Il installe donc devant la fissure une bande de joints intrinsèques souples,
ce qui explique qu'il ramène la part de propagation au niveau de l'intrinsèque (57,4 % contre 58,9 %) au
prix de 1,85 fois plus de joints rompus. Le défaut de contact corrigé le 7 octobre gonflait la
fragmentation de cette variante plus que celle de la référence ; la part de propagation de `tip16` n'a pas
été remesurée avec le contact corrigé.

Trois incohérences de code sont relevées en 2D (section 5.1) : le seuil de cisaillement du balayage ignore
`jointFrictionMobilised = damage` alors que l'activation l'applique ; le décalage δ_n0 n'est pas inversé
sur la parabole de Guo (le 3D l'est depuis le 11 septembre) ; en compression, sous `jointElastic =
parabolic`, la traction de naissance vaut 2σ en 2D et n'est pas exacte non plus en 3D.

---

## 1. Le critère d'insertion, ligne par ligne

### 1.1 Où il se place dans le pas de temps

`FdemSolver::step()` (l. 5428-5470) enchaîne :

1. `elementForces()` : contraintes des triangles à la configuration x_n et forces nodales des éléments ;
2. `bodyForces()` ;
3. `insertionSweep()` si `adaptive_ && !noJoints_` (l. 5449) ;
4. `jointForces()` (l. 5464), qui calcule déjà la traction des joints nés au point 3 ;
5. `generalContact()`, `toolContact()`, puis l'intégration (groupes de nœuds liés).

Il n'y a donc pas de pas de retard entre le franchissement et l'existence du joint : le joint né au pas n
transmet sa traction au pas n, sur la configuration x_n. L'ordre est celui de Yan (§2.4 : « In each
calculation step, stresses in all triangular elements are first calculated »). Le 3D suit le même ordre
(l. 4267 et 4278).

Le dépassement est d'au plus un incrément de contrainte par pas. En charge monotone,

$$
0 \le \sigma_n^{(n)} - f_t \le \Delta\sigma_\text{pas} \approx E\,\dot\varepsilon\,\Delta t ,
$$

négligeable sur le tunnel quasi statique, de l'ordre de quelques pour cent de f_t sous un impact
(ε̇ = 10⁴ s⁻¹, E = 60 GPa, Δt = 10⁻⁹ s donnent 0,6 MPa). Sous la loi par pénalité, ce dépassement est
perdu à l'activation (δ_n0 lit min(σ, f_t), §2.1) : c'est la seule discontinuité temporelle de traction
qui subsiste en mode I. Sous `jointTSL = camacho`, il est conservé (t_ins inclut le dépassement).

Toutes les arêtes qui franchissent l'enveloppe au même pas sont insérées ensemble (liste `hits`, triée par
indice pour le déterminisme, l. 4467-4479), sans réévaluation après chaque activation et sans plafond. C'est
le comportement de Pandolfi et Ortiz (2002) selon `tunnel_edz/REVUE_insertion_adaptative.md` §1 ; je n'ai
pas relu l'article.

### 1.2 La contrainte de facette

Repère de l'arête (l. 4279-4285) : avec a_1, b_1 les deux copies du premier sommet (côté A, côté B), a_2, b_2
celles du second,

$$
P = \tfrac12\big(\mathbf x_{a_1} + \mathbf x_{b_1}\big), \quad
Q = \tfrac12\big(\mathbf x_{a_2} + \mathbf x_{b_2}\big), \quad
\mathbf e = \frac{Q-P}{\lVert Q-P\rVert}, \quad \mathbf n = (e_y,\,-e_x),
$$

n sortante de l'élément A. Tant que l'arête est liée, les copies coïncident et P, Q sont les sommets
courants.

Contrainte moyenne (l. 4294-4303, clé `facetAverage`, lue l. 3786) :

$$
\bar{\boldsymbol\sigma} = w_A\,\boldsymbol\sigma_A + w_B\,\boldsymbol\sigma_B, \qquad
(w_A, w_B) = \begin{cases} (\tfrac12, \tfrac12) & \texttt{arith} \text{ (défaut)}\\[2pt]
\big(\tfrac{V_A}{V_A+V_B}, \tfrac{V_B}{V_A+V_B}\big) & \texttt{volume}\end{cases}
$$

σ_A et σ_B sont les tenseurs constants des deux CST (contrainte nominale de la loi de volume, dégradée
s'il y a endommagement de volume). Projection (l. 4304-4307) :

$$
\sigma_n = \mathbf n\cdot\bar{\boldsymbol\sigma}\,\mathbf n, \qquad \tau = \mathbf e\cdot\bar{\boldsymbol\sigma}\,\mathbf n .
$$

Sous `facetAverage = max` (l. 4330-4352), le couple (σ_n, τ, f_s) du triangle le plus chargé, au sens de
max(σ_n/(d_T f_t), |τ|/f_s), remplace la moyenne. Le 3D dispose en plus de `facetAverage = nodal`
(l. 3427-3438 de `Fdem3dSolver.cpp`) : la traction réellement transmise, obtenue par partition des forces
nodales, à la manière de Camacho et Ortiz ; ce mode n'existe pas en 2D.

### 1.3 L'enveloppe, le DIF et le facteur de pointe

Facteurs dynamiques (l. 4319-4327) : avec ε̇ = ½(ε̇_A + ε̇_B) par défaut, ou deux taux de facette
(ε̇_eq, γ̇) sous `facetAverage = volume` ou `facetRate = tensor`,

$$
d_T = \mathrm{DIF}_t(\dot\varepsilon_\text{eq}), \qquad d_C = \mathrm{DIF}_c(\dot\gamma),
$$

`include/rockim/YangDif.hpp` : DIF_t = 0,95 + 0,41 ε̇^{a_t} borné à [1 ; 1,85] (a_t = 0,07 ou 0,1707 sous
`yang-fig2`), DIF_c = 0,77 + 0,56 ε̇^{a_s} borné à [1 ; 1,84]. Le terme frottant n'est pas amplifié (choix
de Yang et al. 2025, rappelé l. 4309-4313).

Résistance au cisaillement (l. 4328-4330) :

$$
f_s = \max\!\Big(0,\; d_C\,c_J + \tan\varphi_J\,M(\sigma_n)\Big), \qquad
M(\sigma_n) = \begin{cases} \langle -\sigma_n\rangle & \texttt{jointShearEnvelope = yan} \text{ (défaut, éq. 8 de Yan)}\\
-\min(\sigma_n, f_t) & \texttt{yang} \text{ (éq. 1 de Yang 2025, Guo éq. 2.24)}\end{cases}
$$

où f_t, c_J, φ_J sont les propriétés du joint (`assignJointProps`, l. 2700-2800 : propriétés de phase,
facteurs inter-granulaires, surcharges par paire, tirage de Weibull ou champ corrélé). La coupure de
l'enveloppe `yang` lit le f_t statique.

Facteur de pointe (l. 4256-4270 et 4356-4358) : un sommet est une pointe s'il porte un joint non lié
d'endommagement D ≥ D_tip (`insertionTipDamage`, 0,5). Une arête dont l'un des deux sommets est une pointe
voit son enveloppe divisée par k = `insertionTipFactor` :

$$
\kappa = \begin{cases} 1/k_\text{tip} & \text{si un sommet de l'arête est une pointe}\\ 1 & \text{sinon.}\end{cases}
$$

Critère (l. 4365-4370, clé `insertionCriterion`, lue l. 3834) :

$$
\texttt{or} \text{ (défaut)} : \quad \sigma_n \ge \kappa\,d_T f_t \quad\text{ou}\quad |\tau| \ge \kappa\,f_s ,
$$

$$
\texttt{elliptic} : \quad \Phi_F = \left(\frac{\langle\sigma_n\rangle}{\kappa\,d_T f_t}\right)^2
+ \left(\frac{\tau}{\kappa\,f_s}\right)^2 \ge 1
\qquad (\texttt{jtsl::phiInsert}).
$$

Hystérésis (l. 4371-4384, `insertionHoldSteps` = n_h, défaut 1) : le joint n'est inséré qu'après n_h pas
consécutifs au-dessus du seuil ; le compteur `nPhi` est remis à zéro dès un pas sous le seuil. Le
dépassement accumulé pendant ces n_h pas est lui aussi perdu à l'activation.

### 1.4 Ce que le critère ne fait pas

- Il ne lit que des contraintes. Aucune condition d'énergie (taux de restitution, longueur de zone
  cohésive résolue) n'intervient ; l'énergie n'entre que par la loi après insertion.
- Il ne distingue pas une arête en pointe d'une arête en terrain vierge, sauf sous le facteur de pointe.
- Il est local à l'arête : aucun lissage sur un voisinage de rayon ℓ_cz.
- Il ne sélectionne pas de direction : la normale est celle de l'arête, et la cohérence directionnelle
  mesurée est au niveau du hasard (§3.4).

---

## 2. La loi du joint après insertion

### 2.1 L'activation (`activateJoint`, l. 4513-4594)

Dans l'ordre :

1. `J.bonded = false`, horodatage `tInsert`.
2. DIF figé (`stampDif`, l. 4217-4231) : f_t ← d_T f_t, G_I ← d_T G_I, c ← d_C c, G_II ← d_C G_II, puis
   recalcul des longueurs (`setJointLengths`). L'énergie de rupture suit le DIF.
3. Décalage d'ouverture (l. 4538) :
   $$\delta_{n0} = \frac{\min(\sigma_n, f_t)}{p_j}, \qquad p_j = \texttt{insertionPenaltyFactor}\cdot\frac{E}{h}\ (4\,E/h).$$
   Un joint né en traction a δ_n0 = f_t/p_j = δ_nE : il naît au pic. Un joint né en cisaillement garde sa
   traction normale sous-critique ou de compression.
4. Glissement initial (l. 4546-4552) :
   $$f_s^\text{now} = \begin{cases} c & \texttt{jointFrictionMobilised = damage}\ (D = 0)\\ c + \tan\varphi\,M(\sigma_n) & \text{sinon}\end{cases}, \qquad
   \tau_0 = \operatorname{clamp}(\tau, -f_s^\text{now}, f_s^\text{now}), \qquad s_0 = -\tau_0/p_j .$$
   La traction d'essai p_j(δ_t − s_0) vaut τ_0 à glissement géométrique nul.
5. Sous `jointTSL = camacho` (l. 4584-4591) : δ_n0 = 0, s_0 = 0, tampon `jtsl::stampInsertion` (§2.6).
6. Union-find aux deux sommets (`rebindVertex`) : les copies se regroupent par composantes connexes de
   l'éventail sur les arêtes encore liées, ce qui reproduit la fig. 7 de Yan (dédoublement des nœuds
   1, 1', 2, 2', 2'').

Le facteur de pointe n'entre pas dans δ_n0 ni dans f_s^now : un joint né en pointe avec σ_n ≈ f_t/k
naît sur sa branche élastique, au-dessous du pic, avec sa résistance pleine (§4.2).

### 2.2 Ouverture, glissement, points d'intégration

À chaque point k (deux sommets par défaut, deux points de Gauss sous `jointQuadrature = midedge`),

$$
\delta_n = (\mathbf x_{b_k} - \mathbf x_{a_k})\cdot\mathbf n + \delta_{n0}, \qquad
\delta_t = (\mathbf x_{b_k} - \mathbf x_{a_k})\cdot\mathbf e .
$$

δ_n est donc une ouverture effective, décalée ; l'ouverture géométrique o de Yan vaut δ_n − δ_n0.

### 2.3 Branche normale, adoucissement `yan` (alias `munjiza`), l. 6474-6540

Longueurs (`setJointLengths`, l. 3695-3720) :

$$
\delta_{nE} = \frac{f_t}{p_j}, \qquad
\delta_{nF} - \delta_{nE} = \begin{cases} G_I/(f_t I) & \texttt{jointDeltaC = exact}\ (\text{défaut})\\ \max(2\delta_{nE},\,3G_I/f_t) & \texttt{solidity}\\ 3G_I/f_t - \delta_{nE} & \texttt{guo}\end{cases}, \qquad
s_F = \frac{G_{II}}{c\,I},
$$

avec I = ∫₀¹ f(D) dD = 0,386307 pour (a, b, c) = (0,63 ; 1,8 ; 6). f(D) est l'éq. 11 de Yan, identique à
la courbe en z de Munjiza 1999 (`include/rockim/YanSoftening.hpp`).

Endommagement (l. 6495-6503, branche `plastic` par défaut) :

$$
r_n = \frac{\langle \delta_n - \delta_{nE}\rangle}{\delta_{nF} - \delta_{nE}}, \qquad
r_s = \frac{|s_p|}{s_F}, \qquad
D = \min\Big(1,\ \max_{t'\le t}\sqrt{r_n^2 + r_s^2}\Big),
$$

où s_p est le glissement plastique du retour radial (initialisé à s_0, §2.1). Sous `jointShearUnload =
origin`, r_s = (s_max − s_E)/s_F avec s_E = f_s(σ_geo)/p_j le glissement au pic, forme littérale de
l'éq. 14 de Yan lue depuis le pic.

Traction normale en ouverture (δ_n ≥ 0), avec o_max la plus grande ouverture effective atteinte :

$$
\sigma_\text{env} = \min\big(p_j\,o_\text{max},\ f(D)\,f_t\big) \quad(\text{linéaire}), \qquad
\sigma_\text{env} = \min\big(f_t(2r - r^2),\ f(D)\,f_t\big),\ r = o_\text{max}/\delta_{nE} \quad(\texttt{parabolic}),
$$

$$
\sigma = \sigma_\text{env}\,\frac{\delta_n}{o_\text{max}} ,
$$

et, sous `jointSecantRatchet = on`, la sécante σ_env/o_max est rendue non croissante. En compression :

$$
\sigma = k^-\,\delta_n, \qquad k^- = p_j\ (\text{linéaire}),\ 2p_j\ (\texttt{parabolic}),\ \text{multiplié par } (1-D) \text{ sous } \texttt{jointContactPenalty = adaptive}.
$$

Un terme visqueux (`jointXi`, ou η_n de l'option B) s'ajoute à σ sans entrer dans le critère de rupture.

### 2.4 Branche tangentielle et frottement, l. 6600-6760

Cap de cisaillement :

$$
\tau_\text{lim} = f(D)\,c + \mu_\text{eff}\,M(\sigma_\text{crit}), \qquad
\mu_\text{eff} = \begin{cases} \tan\varphi & \text{défaut (frottement non dégradé)}\\ f(D)\tan\varphi & \texttt{jointFrictionScaled = 1}\ (\text{éq. 10 littérale})\\ \mu_r + (\tan\varphi - \mu_r) f(D) & \texttt{jointResidualMu} = \mu_r\end{cases}
$$

σ_crit vaut σ (visqueux compris, défaut) ou sa part élastique (`jointViscInCriterion = off`). Sous
`jointFrictionMobilised = damage`, le terme frottant est multiplié par D (`jtsl::shearCap`).

- `plastic` (défaut) : τ_tr = p_j(δ_t − s_p), τ = clamp(τ_tr, ±τ_lim), retour radial s_p += (τ_tr − τ)/p_j,
  et D mis à jour par l'ellipse. La décharge suit la pente p_j et garde le glissement acquis.
- `origin` : τ = min(p_j s_max, τ_lim)·(δ_t − s_0)/s_max, sécante vers l'origine figée s_0 ; non
  conservative sans cliquet (création ½(k_2 − k_1)s² par cycle, rapport §2.2).
- `solidity` : transcription de Y3Dfd.c, sans mémoire d'endommagement ; écartée des calculs longs.

### 2.5 Branche linéaire historique (`jointSoftening = linear`, défaut du code)

σ = min[(1 − D) p_j δ_n, enveloppe triangulaire de pic f_t et d'ouverture finale δ_nF] ; D est mis à jour
par la sécante, D = 1 − σ_env/(p_j δ_n) ; en cisaillement, D = |s_p|/s_F. Pas d'ellipse : D est le plus
grand des deux. Avec k_I = 2, l'aire de la branche descendante vaut G_I.

### 2.6 Variante `jointTSL = camacho` (§2.4 de la note de septembre 2026)

`include/rockim/JointTsl.hpp`, cablée l. 4584-4591 et 6286-6296, 6511-6533, 6680-6727. Loi initialement
rigide, raccordée à la traction réellement transmise :

$$
\beta = \frac{f_s^\text{now}}{d_T f_t}, \qquad
t_m^\text{ins} = \sqrt{\langle\sigma_n\rangle^2 + \tau^2/\beta^2}, \qquad
G_C = \begin{cases} G_I & \texttt{jointMixLaw = none}\ (\text{défaut})\\ G_I + (G_{II} - G_I)\,m^\eta & \texttt{bk},\ m = \frac{r^2}{\beta^2 + r^2},\ r = |\tau|/\langle\sigma_n\rangle\end{cases}
$$

$$
\delta_m^f = \frac{2G_C}{t_m^\text{ins}}, \quad
\delta_m = \sqrt{\langle\delta_n\rangle^2 + \beta^2\delta_s^2}, \quad
t_m = t_m^\text{ins}\Big(1 - \frac{\delta_m}{\delta_m^f}\Big), \quad
D = \frac{\delta_m^\text{max}}{\delta_m^f}, \quad
t_n = \frac{t_m}{\delta_m}\langle\delta_n\rangle, \quad t_s = \frac{t_m}{\delta_m}\beta^2\delta_s .
$$

Branche montante (`jointTSLRise` = ρ, défaut 10⁻³) : le joint naît au sommet d'une branche fictive de
longueur δ_m0 = ρ δ_m^f, décalée dans la direction de la traction tamponnée ; la raideur de décharge est
bornée par k_0 = t_ins/δ_m0 = t_ins²/(2ρG_C), qui entre au budget du pas de temps. À ρ = 0 la loi est
instable en explicite (mesure du 11 septembre : sécantes jusqu'à 549 p_j). La compression est portée par
p_j δ_n ; le frottement est un curseur de Coulomb écrêté au cap (`jfric::capReturn`).

L'adoucissement camacho est linéaire, et non la courbe en z.

### 2.7 Comparaison terme à terme avec Yan et al. 2023

| élément | Yan 2023 | rockim adaptatif (défauts du code) | verdict |
|---|---|---|---|
| contrainte de facette | moyenne des deux triangles, projetée (§2.4) | idem, `arith` | identique |
| critère, éq. 7 | σ ≥ f_t ou τ ≥ f_s | σ_n ≥ f_t ou \|τ\| ≥ f_s | identique (rockim prend \|τ\|, ce que Yan sous-entend) |
| f_s, éq. 8 | −σ tan φ + c si σ < 0, c sinon | `jointShearEnvelope = yan` | identique |
| ordre dans le pas | contraintes puis insertion | idem, puis forces de joint au même pas | identique |
| continuité à l'insertion | « la contrainte du joint est celle de la facette avant insertion » | δ_n0 = min(σ, f_t)/p_j, s_0 = −τ_0/p_j | identique en traction et en cisaillement non écrêtés ; le dépassement σ − f_t est perdu dans les deux |
| cas 1, déclenché en traction | σ = f(D) f_t ; τ = f(D) τ_ins | σ sur la courbe en z ; τ sur la branche élastique p_j depuis τ_ins, jusqu'au cap τ_lim | diffère : chez Yan τ ne peut que décroître, chez rockim il peut croître jusqu'à c f(D) + frottement |
| cas 2, déclenché en cisaillement | τ = f(D) f_s ; σ = f(D) σ_ins | τ au cap par retour radial ; σ sur la branche élastique p_j depuis σ_ins, jusqu'à f_t | diffère de la même façon pour σ |
| compression, éq. 9 | σ = p_f o | σ = p_j δ_n, p_j = 4 E/h | même forme ; Yan ne chiffre pas p_f pour l'adaptatif |
| D, éq. 12, 14, 16 | o/o_t, \|s\|/s_t, √((o/o_t)² + (\|s\|/s_t)²), o et s depuis l'insertion | r_n identique pour un joint né en traction ; r_s lu sur le glissement plastique (`plastic`) ou depuis le pic (`origin`) | identique en traction pure ; écart d'ordre s_E/s_F en cisaillement |
| o_t, s_t, éq. 13 et 15 | G_I = f_t o_t I, G_II = c s_t I | identiques sous `jointDeltaC = exact` | identique ; `solidity` (configuration d'impact) dissipe 1,159 G |
| f(D), éq. 11 | a = 0,63, b = 1,8, c = 6,01 dans le texte, 6,0 dans la fig. 5 | 0,63 ; 1,8 ; 6,0 | identique à la coquille près |
| irréversibilité | D = D_max (Fukuda) | D ne décroît jamais | identique |
| décharge normale, éq. 17 | σ = f(D_max) T_s o/o_max, sécante vers o = 0 | sécante vers δ_n = 0, soit o = −δ_n0 | diffère : raideur de Yan non bornée quand o_max → 0, celle de rockim bornée par p_j (fig. 1c) |
| décharge tangentielle, éq. 18 | sécante vers s = 0 | `plastic` : pente p_j, glissement conservé ; `origin` : sécante vers s_0 | `origin` ≈ éq. 18 ; `plastic` (défaut) diffère |
| frottement dans le cap | (c − σ tan φ) f(D) : le frottement s'efface avec f(D) | frottement plein (défaut) ; f(D) tan φ sous `jointFrictionScaled = 1` | le tunnel pose `jointFrictionScaled = 1` : identique |
| contact après insertion | contact activé entre les deux triangles dès l'insertion (§2.2) | paire exclue du contact tant que le joint vit ; la compression passe par p_j | diffère ; le texte de Yan ne dit pas comment éq. 9 et contact se partagent la compression |
| nœuds, fig. 7 | dédoublement local | union-find sur l'éventail | identique |
| amortissement | 2μD dans la contrainte, μ = 7,6·10³ Pa·s | Cundall local par défaut ; `bulkViscosity` pour l'éq. 6 | diffère par défaut |

T_s n'est pas défini dans le texte de Yan ; je le lis comme f_t (ou σ_ins dans le cas 2).

### 2.8 Comparaison avec la loi intrinsèque de Munjiza, Wang et Y-Geo

Munjiza (2004) et Wang et al. (2024, fig. 2, éq. 2) : joints présents dès l'origine, branche montante
parabolique,

$$
\sigma = \Big[2\frac{o}{o_p} - \Big(\frac{o}{o_p}\Big)^2\Big] f_t,\quad 0 \le o \le o_p, \qquad
o_p = \frac{2h f_t}{P_f}, \qquad
\sigma = f(D)\,f_t,\quad D = \frac{o - o_p}{o_r - o_p}\ \text{au-delà},
$$

et l'ellipse de l'éq. 2 de Wang écrite depuis le pic. La pente à l'origine vaut 2f_t/o_p = P_f/h ; la
sécante au pic f_t/o_p = P_f/(2h). Avec P_f = 1 000 GPa et E = 10 GPa (Wang, tableau 1), la raideur
initiale vaut 100 E/h et la sécante au pic 50 E/h. Le calcul intrinsèque de contrôle de rockim sur le
tunnel prend p_f = 20 : 20 E/h en branche linéaire, soit 2,5 fois plus souple en sécante que Wang. La
formule o_p = 2h f_t/P_f est celle de Munjiza reprise par Y-Geo ; je l'ai vérifiée par une recherche web
(citation dans un article arXiv), pas dans le livre de Munjiza, que je n'ai pas lu.

Y-Geo (thèse de Lisjak 2013, lue, §3.3) : branche montante linéaire par pénalité p_f, D depuis le pic
(éq. 3.5, 3.8, 3.9), adoucissement en mode II vers un frottement résiduel f_r = σ_n tan φ_f, avec
G_IIc = ∫(τ − f_r) ds (éq. 3.11), et trois points d'intégration par élément de fissure.

rockim intrinsèque : σ = p_j δ_n, p_j = p_f E/h avec p_f = 20 et h le plus petit diamètre inscrit par
défaut ; `jointElastic = parabolic` redonne la parabole de Guo (éq. 2.31), pente 2p_j à l'origine.

Les deux familles ne diffèrent que par la branche avant pic : identique après le pic, absente (Yan) ou
remplacée par un décalage (rockim adaptatif) avant.

### 2.9 Énergie dissipée

Sur un trajet monotone jusqu'à rupture, avec décharge sécante, toute l'aire sous la courbe est dissipée.

| loi | traction pure | cisaillement pur | trajet proportionnel d'angle θ dans (o/o_t, s/s_t) |
|---|---|---|---|
| Yan 2023 | G_I (éq. 13) | G_II en traction ; (c + μ\|σ\|) s_t I = G_II(1 + μ\|σ\|/c) en compression, frottement compris | G_I cos θ + G_II sin θ |
| rockim adaptatif, `yan`, `exact` | G_I exactement (joint né au pic, aire de la seule branche descendante) | G_II au cap c f(D) ; plus le frottement, plein par défaut et non borné en durée | G_I cos θ + G_II sin θ, plus l'écart s_E/s_F |
| rockim, `jointDeltaC = solidity` | 3I G_I = 1,159 G_I | 1,159 G_II | idem × 1,159 |
| rockim, `jointShearRange = coulomb` | inchangé | plage divisée par f_s/c : G_II quelle que soit la pression | idem |
| Munjiza, Wang (intrinsèque, parabole) | G_I + ⅔ f_t o_p = G_I + (4/3) h f_t²/P_f | G_II + ⅔ f_s s_p | idem + branches avant pic |
| rockim intrinsèque linéaire | G_I + f_t²/(2p_j) = G_I + h f_t²/(2 p_f E) | G_II + f_s²/(2p_j) | idem |
| camacho, sans `bk` | ½ t_ins δ_m^f = G_I exactement, dépassement compris | G_I (et non G_II) | G_I |
| camacho, `bk` | G_I | G_II | G_C de Benzeggagh-Kenane, figé à l'insertion |

Trois lectures découlent du tableau. Le mode mixte elliptique (Yan comme rockim) ne dissipe ni G_I ni G_II
ni une loi de mélange calibrée : sur un trajet à 45°, il dissipe 0,707(G_I + G_II), et le maximum sur θ
vaut √(G_I² + G_II²), supérieur à G_II. Le surplus avant pic de l'intrinsèque est négligeable : sur le
tunnel (h = 0,22 m, f_t = 0,6 MPa, P_f = 1 000 GPa, G_I = 20 J/m²), (4/3) h f_t²/P_f = 1,1·10⁻¹ J/m²,
soit 0,5 % de G_I. Enfin `jointTSL = camacho` sans `jointMixLaw = bk` dissipe G_I en cisaillement pur,
20 fois moins que G_II pour le granite de Kuru : la combinaison n'est pas refusée par le code.

---

## 3. Raideur et résistance avant insertion

### 3.1 Raideur

En adaptatif, une arête liée n'a pas de joint : les copies de nœuds sont intégrées par groupe (liaison
rigide), ce qui est exactement l'élément fini à nœuds partagés. La bande pesante le montre : tassement
−6,00866·10⁻⁴ m contre −6,0087·10⁻⁴ m analytique (rapport, r05a). Le pas de temps ne paie pas les joints
liés, sauf sous camacho avec branche montante (k_0 budgété pour toutes les facettes).

En intrinsèque, chaque joint ajoute en série une souplesse 1/p_j = h/(p_f E). Un modèle de joints en série
donne

$$
M_\text{eff} = \frac{M}{1 + \alpha/p_f}, \qquad \frac{c_\text{eff}}{c} = \Big(1 + \frac{\alpha}{p_f}\Big)^{-1/2}, \qquad \alpha = 1{,}239
$$

(ajusté sur huit points du balayage en pénalité de l'onde dans une barre, r04 ; α dépend de la qualité du
maillage par la définition de h). À p_f = 20 : E_eff = 0,942 E, mesuré 0,95 E sur le tunnel.

Une conséquence peu remarquée : les joints insérés et non rompus de l'adaptatif portent p_j = 4 E/h, cinq
fois plus souple que l'intrinsèque à 20 E/h. Dans un anneau entièrement inséré, le même modèle donne
E_eff = E/(1 + 1,239/4) = 0,76 E. Sur le tunnel, 10,8 % des joints insérés restent à D < 0,05 (REVUE
§2) : l'EDZ adaptative est une bande plus souple que la même zone intrinsèque.

### 3.2 Résistance : l'écart de 13 à 16 %

Sur les essais de Yan, rockim retrouve la hiérarchie de l'article mais pas sa convergence : le pic
brésilien intrinsèque plafonne à 4,44 MPa (500E) contre 5,30 MPa en adaptatif, soit 16 % de moins, et la
résistance uniaxiale à 44,3 MPa contre 51,1 MPa, soit 13 % (rapport r05a). Yan écrit que la famille
intrinsèque converge vers l'adaptatif. La raideur, elle, converge dès 100E. Le déficit est donc de
résistance, pas de raideur. Le rapport l'attribue à une accumulation d'endommagement sous-critique dans les
joints intrinsèques ; ce n'est pas démontré. Quatre mécanismes, non exclusifs, sont compatibles avec le
code.

1. Contrainte moyennée contre traction locale. Le joint intrinsèque lit une traction nodale, p_j fois le
   saut de déplacement entre copies, aux deux sommets : c'est la partition des forces nodales, qui garde
   la concentration de contrainte locale. Le critère adaptatif lit la moyenne de deux tenseurs constants.
   Pour un champ régulier, cette moyenne approche la contrainte au milieu de l'arête à O(h²) près ; pour
   un champ concentré (appuis du brésilien, pointe de fissure), elle l'écrête. L'adaptatif paraît alors
   plus résistant ; l'intrinsèque n'a pas de résistance manquante, il a moins de lissage.
2. Souplesse localisante. Les interfaces intrinsèques sont des bandes molles préexistantes ; Falk,
   Needleman et Rice (2001) font de l'espacement des surfaces cohésives une longueur caractéristique. La
   pénalité finie concentre la déformation sur les joints les plus chargés avant le pic.
3. Le chemin. L'intrinsèque offre toutes les arêtes à la fois, et un joint intrinsèque qui atteint le pic
   s'adoucit sans attendre ; l'adaptatif ne crée de surface que là où la moyenne franchit l'enveloppe.
4. L'endommagement sous-critique. En adaptatif aussi, 4,7 % des joints dépassent D = 0,01 avant la
   première rupture du brésilien, et σ_xx au centre tombe à 0,700 fois la solution de Hondros au dernier
   pas avant rupture (r05a). L'argument ne distingue donc pas les deux schémas tant qu'il n'est pas chiffré
   dans les deux.

Les tests de la section 5.2 séparent ces mécanismes.

### 3.3 Dépendance au maillage et absence de critère d'énergie

Un critère en contrainte sur une moyenne de CST est objectif tant que la zone cohésive est résolue :
devant la pointe, la traction est bornée par f_t et varie sur ℓ_cz, que la moyenne voit si ℓ_cz couvre
plusieurs éléments. Quand ℓ_cz ≈ h, la moyenne voit un champ quasi singulier écrêté à l'échelle h, et la
ténacité apparente devient de l'ordre de f_t √h : la propagation dépend du maillage. Sur le tunnel,
ℓ_cz = E G_I/f_t² = 0,556 m et h = 0,22 m : 2,5 éléments en mode I, 7 en mode II (code : « 2,7 éléments
par zone cohésive de mode I », `FdemSolver.hpp` l. 708). La règle d'objectivité h < ℓ_cz/2 est tenue de
justesse ; la littérature des CZM demande en général au moins trois à cinq éléments dans la zone cohésive
(valeur de mémoire, non relue dans ce travail). Le rejeu du 6 octobre sur maillage réduit (h = 0,40 m)
viole la règle et ne vaut que comme tendance.

Le brésilien adaptatif ne converge pas : le pic croît de 25 % quand la maille passe de 1,5 à 0,5 mm (4,73 ;
5,12 ; 5,30 ; 5,93 MPa), alors que ℓ_cz y est résolue par 20 à 70 éléments (r05a). Ce n'est donc pas la
résolution de la zone cohésive. Deux suspects restent : p_j = 4 E/h, qui raidit les joints insérés quand h
diminue, et le critère près des appuis, où le champ est concentré à l'échelle du maillage.

Nucléation contre propagation : un critère à seuil par facette, sans désordre, répond à un champ
uniformément au seuil par une insertion en tapis. C'est le comportement attendu du schéma fondateur
(Camacho-Ortiz, Pandolfi-Ortiz, d'après les volets `biblio_insertion/`, que je n'ai pas revérifiés) : ce
n'est pas une erreur d'implémentation.

### 3.4 Ce que nos mesures montrent

| mesure | adaptatif | intrinsèque | source |
|---|---|---|---|
| joints insérés allant à rupture (D ≥ 0,999) | 76,4 % (homogène), 75,9 % (Weibull m = 6) | sans objet | REVUE §2 |
| joints insérés restés à D < 0,05 | 10,8 % et 11,4 % | sans objet | REVUE §2 |
| angle médian entre arêtes rompues voisines | 53,1° | 54,2° | BILAN §4 ; tirage aléatoire 55,5° |
| arêtes voisines alignées à moins de 30° | 23,2 % | 22,2 % | BILAN §4 ; aléatoire 20,7 % |
| part de propagation, production | 43,7 % | 58,9 % | r05c, tab-tun-intr |
| part de propagation, rejeu du 6 octobre (h = 0,40 m) | 29,2 % | 52,5 % | r05c, tab-tun-rejeu |
| joints rompus, production | 21 259 | 36 527 (+72 %) | r05c |
| blocs d'un seul élément | 75,7 % | 84,9 % | r05c |
| plus gros blocs libres | 37 m² | 318, 120, 114, 104, 68 m² | r05c |
| pic brésilien (essais de Yan) | 5,30 MPa | 4,44 MPa à 500E (−16 %) | r05a |
| résistance uniaxiale | 51,1 MPa | 44,3 MPa à 500E (−13 %) | r05a |

Le critère tire au bon moment : trois joints insérés sur quatre vont à rupture complète. La différence
entre schémas n'est pas l'orientation des fissures, au niveau du hasard dans les deux. Elle tient à la
capacité de l'intrinsèque à fermer des contours autour de grands blocs, par ses 159 269 interfaces
préexistantes. Cette lecture reste conditionnelle à l'injection d'énergie par le contact, maximale en
intrinsèque (18,9 MJ/m de travail net du contact pour 2,66 MJ/m d'énergie cohésive sur le rejeu, avant
correction).

### 3.5 Ce que dit la littérature des CZM extrinsèques

- Papoulia, Sam et Vavasis (2003) ; Sam, Papoulia et Vavasis (2005, EFM 72) : une loi initialement rigide
  activée sur un test au pas de temps présente une discontinuité temporelle ; conséquences nommées :
  oscillations et non-convergence en temps, que raffiner Δt ne corrige pas. Le remède est une condition de
  continuité en temps. Je n'ai lu ni l'un ni l'autre : la recherche web n'a retrouvé que le titre d'une
  communication ASCE 2002 de même contenu (« division by zero and nonconvergence in time ») et l'article de
  2005 ; le contenu est repris de `biblio_insertion/biblio_patho.md`, qui le marque partiellement
  [MEMOIRE]. rockim satisfait cette condition à l'insertion par δ_n0 et s_0, sauf pour le dépassement
  σ − f_t ; camacho la satisfait dépassement compris et borne la raideur par la branche montante, qui est
  précisément la régularisation que la loi rigide n'a pas.
- Molinari et al. (2007, IJNME 69) : en fragmentation, un léger aléa de maillage améliore la convergence en
  énergie jusqu'à deux ordres de grandeur ([VERIFIE] dans `biblio_patho.md`).
- Zhou et Molinari (2004, IJSS 41) : résistance de Weibull par facette normalisée par le volume effectif,
  qui réduit la dépendance au maillage ([VERIFIE] dans le même volet).
- Seagraves et Radovitzky (2010, chapitre de « Dynamic Failure of Materials and Structures ») ; Radovitzky,
  Seagraves, Tupek et Noels (2011, CMAME) : hybride Galerkine discontinu et cohésif extrinsèque, qui évite
  à la fois la souplesse intrinsèque et le changement de topologie en cours de calcul. Non lus ici ; le
  volet bibliographique les marque [MEMOIRE].
- Papoulia, Vavasis et Ganguly (2006) : maillages pinwheel, qui assurent la convergence des trajets
  restreints aux arêtes. Non lu.
- Falk, Needleman et Rice (2001) : la loi intrinsèque modifie les propriétés élastiques et fait de
  l'espacement des surfaces une longueur du calcul ([VERIFIE] dans `biblio_patho.md`, non relu ici).

Aucune référence trouvée ne compare les motifs (blocs contre nuage) entre intrinsèque et extrinsèque en
FDEM (`biblio_fdem.md`, point 5).

---

## 4. Le facteur de pointe

### 4.1 Origine

Ajouté le 24 août 2026 (rockim_p2, branche `insertion-pointe`, porté en 3D), opt-in, bit-identique à k = 1.
La motivation écrite dans le code (`FdemSolver.hpp` l. 704-721) est un déficit de propagation, pas un excès
de fracturation : l'adaptatif ne propageait que 43,7 % de ses ruptures contre 56,8 % pour l'intrinsèque
(le rapport et le BILAN donnent 58,9 % pour la même comparaison ; l'écart de chiffre n'est pas expliqué).
Le diagnostic retenu est que la moyenne de deux CST écrase la singularité de pointe, si bien que la facette
devant une pointe ne se distingue pas d'une facette de l'anneau plastique. Le choix de relâcher la pointe
plutôt que de pénaliser la nucléation vient d'une mesure : pénaliser la nucléation de 1,6 fait passer
`verify_fdem_voronoi_tension` de 15 joints rompus à zéro, c'est-à-dire relève la résistance
macroscopique de 60 % et invalide le calage GBM Red Bohus.

### 4.2 Ce qu'il fait réellement

Le facteur ne divise que le seuil de déclenchement. À l'activation, δ_n0 = min(σ_n, f_t)/p_j ≈ f_t/(k p_j)
et τ_0 ≈ f_s/k : le joint naît sur sa branche élastique, à 1/k du pic, et doit remonter jusqu'à f_t sur un
ressort de raideur 4 E/h (fig. 1b, courbe grise). En pointe, l'adaptatif devient donc localement un schéma
intrinsèque, avec une pénalité cinq fois plus souple que celle du calcul intrinsèque de contrôle. Le facteur
agit par deux canaux : il insère plus tôt, et il installe devant la fissure une bande souple qui concentre la
déformation, le mécanisme de localisation attribué à l'intrinsèque en §3.2. Cela explique que la part de
propagation rejoigne celle de l'intrinsèque (57,4 % à k = 1,6 contre 58,9 %), et pas au-delà.

La définition de la pointe est isotrope et large. Est pointe tout sommet portant un joint non lié de
D ≥ 0,5, y compris les sommets intérieurs d'une fissure déjà rompue. Toutes les arêtes de l'éventail de ces
sommets sont relâchées, y compris celles qui flanquent la fissure. Le facteur favorise donc la ramification
et la granulation le long des lèvres autant que l'avancée en pointe ; c'est cohérent avec la mesure : la part
de blocs d'un seul élément passe de 75,7 à 88,3 % et les joints rompus sont multipliés par 1,85 à temps
égal, R⁹⁵ de l'EDZ +31 % sur le rejeu.

Le compteur interne (`nProp_`, `nNuc_`, l. 4472-4477) compte une insertion à un sommet de pointe ; la part
de propagation du rapport (outil `nucleation_vs_propagation.py`) compte les ruptures qui prolongent un joint
déjà rompu. Ce ne sont pas les mêmes observables.

### 4.3 Lien avec le défaut de contact corrigé le 7 octobre

`docs/rapport_guide/ENQUETE_CONTACT.md` établit que presque toute l'injection d'énergie du contact sur le
tunnel vient de paires de triangles à sommet commun : tant qu'aucun joint de l'éventail n'est mort, aucun
des deux n'est candidat au contact (`act_` ne contient que les faces extérieures et les lèvres des joints
morts), leurs coins se recouvrent sans force de rappel quand les joints de l'éventail glissent ou
s'écrasent (65 % des joints meurent en compression), puis la paire naît d'un coup avec ce recouvrement et
le potentiel de Munjiza crée p∫(φ_A + φ_B) dA sans travail (gcBirthTau = 1 µs < Δt = 5,8 µs).

Le facteur de pointe multiplie les joints insérés dans les éventails, donc ces recouvrements : sur le rejeu,
le travail net du contact valait 1,94 MJ/m pour 0,96 MJ/m d'énergie cohésive avec k = 1,6, contre 0,027 et
0,43 MJ/m sans. Le diagnostic 4 concluait déjà que « le facteur de pointe exposait le biais sans en être la
cause ». La correction (`correction_contact/IMPLEMENTATION.md`, TESTS.md) sur `tip16_mu0` (T = 0,12 s,
p = 10¹⁰ N/m) : `contactCandidates = vertex` fait tomber l'énergie de naissance de 1,45·10⁶ à 16 J/m, le
travail du contact de +1,31·10⁶ à −8,8·10⁴ J/m, la vitesse nodale maximale de 24,7 à 2,5 m/s ; les joints
rompus passent de 4 889 à 4 317 (−12 %) et l'énergie cohésive de 393 à 317 kJ/m (−19 %). Une partie de la
fragmentation de `tip16` était donc entretenue par l'énergie créée. La part de propagation, les blocs et
l'EDZ de `tip16` et de la référence adaptative n'ont pas été remesurés avec le contact corrigé : les chiffres
de la section 4.2 restent ceux d'avant correction. Le faciès dépend encore de la raideur de contact
(−24 % de joints rompus entre p = 10⁹ et 10¹¹ N/m, corrigé).

---

## 5. Ce qui n'est pas maîtrisé, et propositions

### 5.1 Incohérences de code relevées (2D surtout)

1. Seuil de cisaillement et `jointFrictionMobilised = damage` (2D). Le balayage teste |τ| ≥ c + tan φ
   M(σ_n) sans regarder `fricMob_` (l. 4328, 4421), alors que `activateJoint` écrête τ_0 à c seul sous cette
   clé (l. 4546-4552). Un joint né en compression-cisaillement perd μ|σ_n| de traction au pas de sa
   naissance, le saut que la continuité devait supprimer. Le 3D est cohérent (l. 3381-3390 : seuil c seul
   sous `fricMob_`), ce qui insère en revanche en compression dès |τ| ≥ c, bien avant Mohr-Coulomb.
   Correction décrite : en 2D, reprendre l'expression du 3D dans les deux branches du balayage, ou faire le
   choix inverse dans les deux solveurs (seuil Mohr-Coulomb, frottement mobilisé à l'insertion par un D
   initial non nul). À trancher, car les deux choix changent l'instant d'insertion en compression.
   Test : joint isolé en compression-cisaillement, traction tangentielle avant et après l'activation égale
   à 10⁻¹² près.
2. δ_n0 sous `jointElastic = parabolic` (2D). `activateJoint` 2D pose δ_n0 = s/p_j ; la branche
   parabolique rend alors f_t(2r − r²) > s pour un joint né sous le pic (+25 % de f_t à s = f_t/2). Le 3D
   l'inverse depuis le 11 septembre (l. 3606-3616). Correction décrite : porter les dix lignes du 3D.
3. Compression sous `parabolic` (2D et 3D). En compression la pente est 2p_j ; avec δ_n0 = s/p_j, le 2D
   rend 2s ; le 3D, qui inverse la parabole, rend 2f_t(1 − √(1 + |s|/f_t)), soit −0,83 f_t pour s = −f_t et
   −4,6 f_t pour s = −10 f_t. Correction décrite : pour s < 0, δ_n0 = s/(k⁻), k⁻ étant la pente de
   compression réellement appliquée ((1 − D) compris sous `jointContactPenalty = adaptive`). Test :
   traction de naissance pour s ∈ {−10 f_t, −f_t, 0, f_t/2, f_t} sous `linear` et `parabolic`.
   Ces deux défauts sont dormants : aucune configuration de référence ne combine `parabolic` et
   `adaptive`.
4. DIF continu et adaptatif. Sous `strainRateDIFArm = continuous`, f_t et donc δ_nE sont recalculés à chaque
   pas, alors que δ_n0 est figé à l'insertion. Un joint né au pic se retrouve dans l'adoucissement si le
   taux baisse (r_n = (δ_n0 − δ_nE)/plage > 0 sans ouverture), ou sous le pic s'il monte. Combinaison non
   refusée, non utilisée. Correction décrite : refuser `continuous` sous `adaptive`, ou recalculer δ_n0 en
   proportion de δ_nE.
5. Camacho sans `jointMixLaw = bk` dissipe G_I en cisaillement (§2.9). Correction décrite : avertissement
   à la lecture des clés, ou G_C = G_II pour un joint déclenché en cisaillement pur.
6. Écarts de parité 2D/3D sans justification écrite : `facetAverage = nodal` n'existe qu'en 3D ; le
   compteur d'hystérésis est tenu en permanence en 3D, seulement si n_h > 1 en 2D (sans effet sur le
   résultat) ; la définition de la pointe porte sur deux sommets en 2D, trois en 3D (cohérent).

### 5.2 Ce qui n'est pas compris

1. L'écart de 13 à 16 % entre le plateau intrinsèque et l'adaptatif sur les essais de Yan.
2. La non-convergence en maillage du brésilien adaptatif (+25 % de 1,5 à 0,5 mm).
3. La part de propagation et le motif en blocs après correction du contact, avec et sans facteur de pointe.
4. La convergence en temps de l'insertion : jamais testée sur un faciès (seul le bilan d'énergie l'a été).
5. Le rôle de la souplesse des joints insérés intacts (4 E/h, E_eff ≈ 0,76 E dans un anneau inséré) dans la
   granulation de l'EDZ.
6. Le partage traction-cisaillement à l'insertion : 94 % de ruptures en cisaillement en intrinsèque sur le
   rejeu, 63 % en adaptatif ; on ne sait pas lequel est physique pour Hutou Beishan.
7. L'énergie dissipée en mode mixte (G_I cos θ + G_II sin θ sous l'ellipse), jamais confrontée à un essai.

### 5.3 Propositions classées

Priorité 1, cohérence et mesures sans changer la physique par défaut.

1. Corriger 5.1-1, 5.1-2 et 5.1-3 derrière la règle de croissance par addition (défauts bit-identiques
   quand la combinaison concernée n'est pas armée), avec le banc d'un joint isolé décrit.
2. Instrumenter l'insertion : par joint inséré, le mode de déclenchement (traction, cisaillement, les deux),
   le dépassement (σ_n − d_T f_t)/f_t, le saut de traction à la naissance et l'énergie perdue
   ½(σ_n − f_t)²/p_j ; sortie dans `fdem_joints_*.vtu` et en histogramme de fin de run. Sans cela, la
   continuité temporelle n'est vérifiée que sur des cas jouets.
3. Remesurer le tunnel avec le contact corrigé (`contactCandidates = vertex`, `gcBirth = offset`), sur le
   maillage de production : référence adaptative, k = 1,3 et 1,6, intrinsèque ; part de propagation, blocs,
   EDZ, à temps égal. Tant que ce n'est pas fait, aucun chiffre du facteur de pointe n'est utilisable.

Priorité 2, séparer les mécanismes de l'écart adaptatif-intrinsèque (brésilien et compression de Yan, à la
maille de 0,75 mm de l'article).

4. Critère adaptatif sur la traction nodale : porter `facetAverage = nodal` du 3D au 2D, puis comparer le
   pic à `arith`, `volume` et `max`. Si le pic adaptatif tombe vers 4,4 MPa, le mécanisme 1 de §3.2 est
   confirmé.
5. Intrinsèque à souplesse nulle : intrinsèque à 500E avec `jointQuadrature = midedge` contre `vertex`, et
   comptage de D > 0,01 avant le pic dans les deux schémas. Mesure le mécanisme 4.
6. Balayage `insertionPenaltyFactor` ∈ {4, 10, 40, 100} sur le brésilien adaptatif aux trois mailles. Si la
   dérive de +25 % disparaît à p_j élevé, la cause est la raideur des joints insérés.
7. Convergence en temps : Δt, Δt/2, Δt/4 sur le brésilien et sur le tunnel court, avec n_h = 1 et 3.

Priorité 3, changer la manière d'insérer (opt-in, à valider sur les bancs précédents).

8. Pointe orientée : ne compter comme pointe qu'un sommet portant exactement un joint rompu (extrémité
   de fissure), et ne relâcher que les arêtes de l'éventail dont la normale est proche de la direction
   principale de traction au sommet. Attendu : même part de propagation qu'à k = 1,6, moins de granulation
   le long des lèvres.
9. Naissance en pointe au pic : sous facteur de pointe, faire naître le joint au pic de sa branche
   (résistance du joint ramenée à la traction de déclenchement, f_t ← σ_n, G_I conservé) plutôt que sur un
   ressort de 4 E/h. Sépare l'effet du déclenchement anticipé de celui de la bande souple.
10. Camacho avec `jointMixLaw = bk` et branche montante sur le tunnel : continuité dépassement compris,
    énergie G_C exacte et indépendante du maillage. C'est la formulation la plus proche de la littérature
    des CZM extrinsèques (§3.5).
11. Désordre géométrique (jitter des nœuds, Molinari 2007) à résistances homogènes, puis champ corrélé avec
    normalisation par volume effectif (`jointSizeEffect`), qui n'était pas activée dans le balayage de
    Weibull du tunnel.
12. Critère non local ou en énergie pour la nucléation : contrainte moyennée sur un disque de rayon ℓ_cz/2,
    ou condition de taux de restitution pour une arête hors pointe. Plus lourd ; à n'envisager qu'après 4
    à 9.

---

## 6. Sources

Lues pour cette note : `src/FdemSolver.cpp` (insertionSweep, facetFrame, activateJoint, stampDif,
refreshDif, setJointLengths, assignJointProps, step, jointForces, lecture des clés), `src/Fdem3dSolver.cpp`
(insertionSweep, activateJoint), `include/rockim/JointTsl.hpp`, `YangDif.hpp`, `YanSoftening.hpp`,
`FdemSolver.hpp` (commentaires du facteur de pointe), `docs/notes/BILAN_insertion_adaptative.md`,
`docs/notes/AUDIT_loi_adaptative_2026-09-11.md`, `tunnel_edz/REVUE_insertion_adaptative.md`, `docs/notes/CHANGES_YAN.md`,
`docs/SPEC_loi_note_2026.md` (début), `docs/rapport_guide/ENQUETE_CONTACT.md`,
`docs/rapport_guide/correction_contact/{IMPLEMENTATION,TESTS}.md`, sections r02, r05a, r05c du rapport,
`biblio_insertion/` (passages cités), Yan et al. 2023 (PDF pages 3 à 8 et texte), Wang et al. 2024 (texte
extrait et figure 2), thèse de Lisjak 2013 (§3.3).

Non lus, ou lus seulement de seconde main :

- Munjiza 2004 (livre) et Munjiza, Owen et Bicanic 1999 : non disponibles dans `bibliographie/` ; la forme
  o_p = 2h f_t/P_f est vérifiée par une recherche web seulement ([résultat ASCE Y-Geo](https://ascelibrary.com/doi/10.1061/%28ASCE%29GM.1943-5622.0000216),
  [Munjiza 1995](https://www.doi.org/10.1108/02644409510799532)) ; arxiv.org est bloqué par le proxy.
- Camacho et Ortiz 1996, Pandolfi et Ortiz 2002, Falk, Needleman et Rice 2001, Papoulia, Sam et Vavasis
  2003, Sam et al. 2005, Seagraves et Radovitzky 2010, Radovitzky et al. 2011 : absents de
  `bibliographie/` ; repris des volets `biblio_insertion/` avec leurs marques [VERIFIE] ou [MEMOIRE]. La
  recherche web n'a pas retrouvé l'article 2003 sous ce titre
  ([publications de Papoulia](https://math.uwaterloo.ca/~papoulia/publications.html),
  [communication 2002](https://invenio.itam.cas.cz/record/4823?ln=en)).
- Fukuda et al. 2021 : non lu ; cité par Yan pour l'irréversibilité et l'activation du contact.
- Guo 2014 : non relu ; équations 2.24, 2.28, 2.30-2.33 reprises du rapport (r02).
- Mahabadi et al. 2012 (Y-Geo) : non lu ; loi Y-Geo prise dans la thèse de Lisjak.
- La note de travail de septembre 2026 (« Lois constitutives proposées pour un FDEM hybride à insertion
  adaptative ») : introuvable dans les deux dépôts ; connue par la SPEC, l'AUDIT et les commentaires du code.
