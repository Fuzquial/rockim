# PLAN — la loi de joint, d'abord mathématique, puis numérique (12/09/2026, 1 h)

Demande de Fernando : « corrige tous les problèmes mathématiques, écris les équations
correctement, regarde les vraies équations, puis passe à l'implémentation numérique ».
Ce document est le plan soumis au conseil avant tout code. Il ne lance rien.

Sources lues cette nuit, dans l'ordre d'autorité :
1. **Solidity, `src/Y3Dfd.c`, fonction `Sigma_tau`** (GitHub `ImperialCollegeLondon/solidity-solver-open`,
   LGPL) — le code que Yang et al. 2026 ont fait tourner. Transcrit verbatim (§1.2).
2. **Yang et al. 2026**, IJRMMS 206:106660, §2 (p. 2-4) : renvoie la loi de joint à la thèse de
   Guo 2014, ajoute le DIF (éq. 1-2) et l'endommagement de volume (éq. 3-4). Table 1 (p. 6).
3. **Guo 2014**, thèse, §2.3.3 (éq. 2.24, 2.30, 2.31) — PDF > 100 Mo, texte non extractible ici ;
   fiche `BIB/guo2014.md`. La transcription 1 fait foi (c'est son implémentation).
4. **Yan et al. 2023**, IJRMMS 169:105439, §2.5, éq. 9-18 (z-curve f(D), moteur mixte,
   décharge sécante) — via `include/rockim/YanSoftening.hpp` et les commentaires du solveur ;
   le PDF n'est pas dans la bibliographie (à ingérer).
5. **rockim_g1**, `src/Fdem3dSolver.cpp` l. 4351-4660 (3D) et `src/FdemSolver.cpp` l. 5940-6400
   (2D, miroir déclaré exact) — la loi telle qu'elle est codée, toutes branches.

---

## 1. Les équations

Notations, par point d'intégration k d'une facette : δ_n ouverture (> 0 en traction),
δ_s glissement (vecteur dans le plan en 3D, scalaire en 2D), pj = pf·E/h raideur de pénalité
(rockim) ou pe/el (Solidity), ft, c, tanφ résistances, G_I, G_II énergies de rupture,
f(D) la z-curve de Munjiza (a = 0,63, b = 1,8, c = 6 ; f(0) = 1, f(1) = 0, ∫₀¹f = 0,386).

### 1.1 rockim, telle que codée (deck v2 : `munjiza` + `parabolic` + `yang` + `coulomb` + `origin`)

**Normal.** Élastique parabolique (Guo 2.31), o_max = max historique de δ_n :

    δ_n ≥ 0 :  σ = min( env_E(o_max), f(D)·ft ) · δ_n / o_max,     env_E(o) = ft·(2r − r²), r = o/δ_E, δ_E = ft/pj
    δ_n < 0 :  σ = 2·pj·δ_n                                                     (pente de la parabole à l'origine)

Décharge en mode I = sécante à l'origine passant par (o_max, enveloppe) : **éq. 17 de Yan**.

**Moteur d'endommagement (éq. 16 de Yan), ratchet** D ← max(D, √(r_n² + r_s²)) avec
r_n = (δ_n − δ_E)/(δ_F − δ_E) et, selon la branche de cisaillement,

    plastic :  r_s = |s_p| / s_F                              s_F = max(2 s_E, 3 G_II / c) (« cohesion ») 
    origin  :  r_s = (s_max − s_E) / den,   den = s_F   ou   den = max(2 s_E, s_F·c/f_s)   (« coulomb »)

avec le seuil élastique de glissement lu **à la contrainte normale courante** :

    s_E = ( c + tanφ·F(σ_n) ) / pj,    F(σ) = −min(σ, ft)  (enveloppe `yang`, cut-off en traction)
    f_s = c + tanφ·max(0, −σ_n)         (résistance de Mohr-Coulomb courante)

**Cap de cisaillement** τ_lim = f(D)·c + tanφ·F(σ_n)  (cohésion endommagée, frottement entier).

**Cisaillement, branche `plastic`** (défaut) : élasto-plasticité à retour radial

    τ_tr = pj·(δ_s − s_p) ;   si |τ_tr| > τ_lim :  τ = τ_lim·τ_tr/|τ_tr|,   s_p += (τ_tr − τ)/pj ;  sinon τ = τ_tr

**Cisaillement, branche `origin`** (éq. 18 de Yan transcrite) : sécante à l'origine

    s_eff = δ_s − s₀,  s_max = max historique de |s_eff|
    τ = min( pj·s_max, τ_lim(σ_n) ) · s_eff / s_max

**Mort** : D ≥ 1 sur deux points sur trois (`majority`), quel que soit le signe de δ_n (`damage`) ;
la facette passe au contact par potentiel, naissance par rampe (`gcBirth = ramp`, force nulle sur
l'offset de naissance, relaxé sur gcBirthTau = 1 µs) ou par continuité de force (`penalty`).

### 1.2 Solidity, `Sigma_tau`, telle que codée (verbatim, par point d'intégration)

    op = 2·el·ft/pe ;  ot = max(2 op, 3 G_I/ft)
    σ_tmp = 2·o·ft/op  si o < 0  (= pe·o/el, pénalité linéaire), sinon 0⁺
    f_s   = c − tanφ·σ_tmp  si σ_tmp ≤ 0  ;  c sinon                      ← MC à la pression COURANTE
    sp = 2·el·f_s/pe ;  st = max(2 sp, 3 G_II/f_s)                        ← plage de mode II divisée par f_s
    D  : (o−op)/ot, (|s|−sp)/st, ou √ des carrés si les deux dépassent ; z = f(D) ; **D lu sur o et s COURANTS**
    σ  = pe·o/el (o<0) ;  ft·z·(2ô − ô²) (0≤o≤op, ô = o/op) ;  ft·z (o>op)
    τ  = z·f_s·(2ŝ − ŝ²) (|s|≤sp, ŝ = |s|/sp) ;  z·f_s (|s|>sp) ;  (− dpefm·σ en compression, dpefm = 0)
    nfail > 1 (deux points sur trois) → élément joint retiré, deux faces de contact créées.

Trois faits, lus dans le code et absents des articles :
- **D n'a pas de mémoire** : z est recalculé sur l'ouverture et le glissement courants. La loi de
  Solidity est **élastique non linéaire réversible** jusqu'à la mort (le joint « guérit » en
  décharge). Aucune dissipation avant la rupture ; G_I, G_II ne sont libérées qu'au retrait.
- **Aucun frottement dans le joint** (dpefm = 0) : le frottement n'existe qu'au contact, après la mort.
- La pression n'entre dans τ qu'au second ordre : τ = pe·s/el − (pe·s)²/(4 el² f_s) ; la
  raideur initiale pe/el ne dépend pas de σ_n.

### 1.3 Ce qui est mathématiquement faux, et où

Condition à satisfaire par toute loi cohésive (thermodynamique des interfaces, Ortiz-Pandolfi
1999, Lisjak-Grasselli 2014 §3) : il existe une énergie libre W(δ_n, δ_s ; α) avec α les variables
internes (D, s_p), traction = ∂W/∂δ, et la dissipation Φ = −∂W/∂α · α̇ ≥ 0 le long de tout
chemin. Corollaire : la partie « élastique » doit dériver d'un potentiel, ∂σ/∂δ_s = ∂τ/∂δ_n.

**(M1) rockim `origin` : violée au premier ordre.** τ = k·s_eff avec k = τ_lim(σ_n)/s_max dès que
le cap est actif ; ∂τ/∂δ_n = s_eff·tanφ·∂F/∂δ_n·(1/s_max) ≠ 0 = ∂σ/∂δ_s. Cycle fermé à glissement
s fixé : comprimer (k monte de k₁ à k₂, aucun travail de cisaillement), glisser en retour
(travail rendu ½k₂s²), décomprimer, recharger (travail payé ½k₁s²) : **gain ½(k₂ − k₁)s² > 0 par
cycle**, indépendant de dt et du maillage. Sous l'insert, σ_n oscille de ±100 MPa autour de
−200 MPa : τ_lim = 30 + 1,85·200 = 400 MPa varie de ±46 %. C'est la pompe mesurée en 2D (V0 vs
V19) et en 3D (référence vs B1) ; les conventions Solidity n'en changent que l'amplitude
(B4 : −0,1 J ; B4a −8,6 J ; B4c −12 J ; référence −16 J).

**(M2) Solidity : violée au second ordre** dans le même régime (f_s entre dans la courbure de la
parabole de cisaillement et dans sp). Gain par cycle O(ŝ³·Δf_s·sp). Leur code n'a aucun bilan
d'énergie (`Y3Dsd.c` : différences centrées, sans amortissement) ; l'écart n'est mesuré nulle
part chez eux. À quantifier chez nous par le banc du §3 : si l'amplitude est < 1 % du travail de
l'outil, on l'accepte comme « la loi publiée » ; sinon on ne la transcrit pas.

**(M3) rockim `plastic` : conforme.** W = ½·2pj⟨−δ_n⟩² + W_n⁺(δ_n ; D) + ½pj(δ_s − s_p)², le cap
n'entre que dans le critère, Φ = τ·ṡ_p ≥ 0 (retour radial le long de τ, τ_lim ≥ 0), et la
composante d'endommagement dissipe car D est un ratchet et f(D) décroît. C'est la
plasticité de Mohr-Coulomb non associée standard sur une interface. B1 le confirme (joints
+0,24 J, leapfrog 0,38 J sur 90 µs).

**(M4) `jointShearRange = coulomb` est attaché à `origin` par une garde**, alors que la
normalisation 3 G_II/f_s(σ_n) de Solidity ne concerne que le **moteur** r_s (une longueur de
référence), pas la raideur. Sur la branche `plastic`, r_s = |s_p| / max(2 s_E, s_F·c/f_s) est
dissipatif (D ratchet). La garde est une erreur de transcription, pas une nécessité.

**(M5) Le relais joint → contact en compression.** À la mort, le joint porte F = 2pj·A·⟨−δ_n⟩
(jusqu'à 3 à 12 kN par facette, 0,9 à 3,7 MN cumulés sur le banc). La rampe fait naître le
contact à force nulle : le chemin d'effort est coupé pendant ~1 µs, les deux éléments se
détendent, puis le contact les rattrape par pénalité sur des nœuds de 6·10⁻⁷ kg. Ce n'est pas
une loi, c'est une discontinuité de force : elle nourrit le terme leapfrog Σf²dt²/2m (+7 J sur la
référence, +1,7 J encore sur B4 sain, 20 % du travail). Solidity impose la **continuité de
force** (Y3Did.c l. 915-964 : pénalité de la paire recalée sur F_joint/F_contact, bornée
[0,01 ; 3]). rockim la porte sous `gcBirth = penalty` mais l'applique aussi aux paires **sans**
joint mort (facteur 1 sur un recouvrement déjà formé) : c'est le 1 J de la naissance piston/bit.

**(M6) La sécante non croissante est la forme dissipative de l'éq. 18.** Si l'on veut garder la
décharge à l'origine de Yan (pas de glissement résiduel), la condition Φ ≥ 0 impose k(t) =
min(k(t⁻), τ_lim(σ_n)/s_max) — un endommagement de cisaillement d_s = 1 − k/pj monotone. La
pression ne peut alors que **réduire** la raideur, jamais l'augmenter. Mais cette forme n'a pas de
frottement cyclique (une interface comprimée qui glisse et revient ne dissipe rien) ; Solidity
n'en a pas non plus dans le joint. `plastic` en a.

---

## 2. Décisions de loi proposées

D1. **Branche de référence = `plastic`** (M3) avec la plage coulomb rattachée (M4) : la loi
« Mohr-Coulomb élasto-plastique + endommagement z-curve à ratchet », qui est ce que Guo/Yang
décrivent en mots et ce que la thermodynamique impose. Différence assumée avec Solidity :
frottement de joint avant la mort (le leur est nul) et endommagement irréversible (le leur guérit).
Les deux différences vont dans le sens dissipatif.

D2. **`origin` reste disponible mais non conservatif déclaré** : impression d'un AVERTISSEMENT
au démarrage (« ressort paramétrique, crée de l'énergie sous σ_n variable, disqualifié pour
l'impact »), et refus sous `budgetAbortPct` sauf clé explicite `allowNonConservative = true`.
Principe VIII (croissance par addition) : rien n'est retiré.

D3. **Nouvelle valeur `jointShearUnload = damage`** = sécante non croissante (M6), opt-in, pour
qui veut la décharge à l'origine de Yan sans la pompe. Elle sert de contrôle : même géométrie que
`origin`, dissipative par construction.

D4. **Relais** : `gcBirth = penalty` réservé aux paires nées d'un joint mort (continuité de
force, Solidity) ; les paires sans joint mort naissent à force nulle sur un offset **figé**
(pas de relaxation gcBirthTau, qui injecte à son tour). Nouvelle valeur `gcBirth = relay`.

D5. **L'endommagement de volume de Yang (éq. 3-4)** : σ = (1 − D)σ̄ avec D linéaire en δ_m
ratcheté ; dérive d'un potentiel (1 − D)W_NH avec D monotone → dissipatif. Rien à corriger.
Le DIF (éq. 1-2) multiplie des résistances (critères), pas des raideurs → rien à corriger.

---

## 3. Le banc falsifiant, avant tout run

`tools/joint_cycle.cpp` : un pilote point-matériel qui appelle **les mêmes fonctions** que le
solveur (à extraire de `processJoint` dans un en-tête `JointLaw.hpp`, comme `YanSoftening.hpp`
l'a été pour `yan_point.cpp`). Chemin imposé (δ_n, δ_s) en quatre segments :
compression à s = 0 → glissement s₁ à σ_n = σ₁ → compression jusqu'à σ₂ à s fixé → retour à
s = 0 à σ₂ → retour à σ₁. Sortie : W = ∮(σ dδ_n + τ·dδ_s).
- **Doit passer** : `plastic` (W ≥ 0), `damage` (W ≥ 0), pour trois amplitudes (ŝ = 0,3 ; 0,9 ; 1,5).
- **Doit échouer** : `origin` (W < 0, et W = −½(k₂−k₁)s₁² à 1 % près, la formule de M1).
- Mesure M2 : la loi de Solidity transcrite telle quelle (parabole en s/sp) — W < 0 attendu, à
  chiffrer.
Coût : secondes. Entre dans `verify` (tier fast) avec le deck 2D `origin` qui DOIT aborter sous
la borne KE (le banc à échec obligatoire de la règle du 25/08).

---

## 4. Implémentation numérique, dans cet ordre, une seule fois

N1. `JointLaw.hpp` : extraction sans changement de la loi (bit-identité 8/8 exigée avant toute
    modification — c'est le refactor qui porte le risque).
N2. `jointShearUnload = damage` (D3) et levée de la garde coulomb sur `plastic` (M4, D1).
N3. Avertissement + refus `origin` sous la borne (D2).
N4. `gcBirth = relay` (D4). Mesure : banc s = 2,5 `plastic`, `ramp` vs `relay` : terme leapfrog
    et KE bound, puis test dt/2 (le terme doit être divisé par 2, sinon il ne vient pas de dt).
N5. Banc court s = 2,5 sous la loi D1, 300 µs : les sept critères de Yang, avec le banc de
    §3 archivé à côté. Puis seulement la décision du run s = 1 (coût à chiffrer, ETAT §9).
Hors périmètre ici : estimateur de facette adaptatif (`max`), coût, maillage s = 1.

## 5. Ce que ce plan ne sait pas encore

- ~~Le texte exact des éq. 17-18 de Yan 2023 (PDF absent)~~ **Lu le 13/09 (PDF déposé)** : l'éq. 18
  publiée est bien τ = (c − σ·tanφ)·f(D_max)·|s|/s_max avec σ la contrainte normale COURANTE — la sécante
  paramétrique est dans l'article ; `jointSecantRatchet` est « Yan amendé » (AUDIT §8).
- La grandeur de M2 chez Solidity : mesurée au §3, pas supposée.
- Si `plastic` fracture assez : B1 n'a rompu que 4 joints en 90 µs ; le banc 300 µs en cours
  (`yang2026_bench_s25_plastic`, 60 % à 1 h 05) répondra avant le conseil ou juste après.
