# Correction du contact par potentiel de rockim : ce que dit la littérature

Date : 2026-10-07. Lecture seule de `/home/user/rockim` (rien n'y a été modifié). Point de départ :
`rockim/docs/rapport_guide/ENQUETE_CONTACT.md` (défauts a à d). Textes locaux convertis par pdftotext dans
`correction_contact/txt/`. Les fichiers commencent par `._` : utiliser `ls -a` ou le motif `.??*`.

Étiquettes de vérification :
- `[LU]` : passage lu dans le PDF local, page citée ;
- `[VAULT]` : repris d'une note du vault (`phd_geothermie/FDEM/rockim/biblio_insertion/*.md`) qui déclare avoir lu le texte intégral. Je ne l'ai pas relu aujourd'hui : arxiv.org est bloqué par le proxy ;
- `[RÉSUMÉ]` : connu seulement par le résumé ou un extrait de moteur de recherche, sans le texte intégral ;
- `[DÉRIVÉ]` : dérivation faite ici, contrôlée numériquement dans l'enquête (`loop2.cpp`), sans source publiée trouvée ;
- `[MÉMOIRE]` : souvenir de lecture, non revérifié. À contrôler avant toute citation.

Accès réseau : le proxy refuse arxiv.org, sciencedirect, springer, osti, sustech et e-bookshelf. Je n'ai eu
que les extraits de recherche web. **Le livre de Munjiza (2004), Munjiza & Andrews (2000), Munjiza & Andrews
(1998), Munjiza-Knight-Rougier (2011, 2015) et Mahabadi et al. (2012, Y-Geo) ne sont PAS dans la
bibliographie locale** et n'ont pas pu être lus (voir §7).

---

## 1. (a) Paires candidates : qui doit voir qui ?

### 1.1 Ce que fait la lignée Munjiza

**Y / Y-Geo / Irazu, contact « intrinsèque » où tout le maillage est discret dès t = 0.**
- « In actual implementation, the separation of adjacent element edges is assumed in advance through the
  topology of adjacent elements being described by different nodes. Thus **no two elements share any
  nodes** – the continuity between elements is enforced through the penalty function method »
  (AbuAisha et al. 2015, *FDEM by Y-Geo: an overview*, MIC vol. 5, ch. 7, p. 67) `[LU]`.
- « Contact interaction forces are calculated **between all pairs of elements that overlap in space** »
  (Tatone & Grasselli 2015, IJRMMS 75, §2.2, p. 57 `[LU]` ; même phrase chez Lisjak 2013, thèse, §4.3.1,
  p. 64 `[LU]`, et chez AbuAisha 2015, p. 70 `[LU]`). La détection NBS (Munjiza & Andrews 1998) porte sur
  tous les éléments.
- Le rôle du contact pour les éléments liés : « a finite stiffness is required for the crack elements …
  represented by the normal, tangential and fracture penalty values, pn, pt and pf, **for compressive,
  shear and tensile loading conditions**, respectively » (Lisjak 2013, §3.3.2, p. 31 `[LU]` ; AbuAisha 2015,
  p. 67 `[LU]`). Dans Y-Geo, la COMPRESSION à travers un joint intact est donc portée par la pénalité de
  contact pn, et non par le joint. Les paires adjacentes (arête commune liée) sont des paires de contact
  ordinaires. Les paires qui ne partagent qu'un sommet aussi, puisqu'aucun nœud n'est partagé.

  Conséquence : dans Y, une paire « à sommet commun » ne peut pas naître avec un recouvrement profond. Elle
  est vue dès que les deux triangles s'interpénètrent.

**Lignée 3D Imperial (Xiang, Guo, Latham : Y3D puis Solidity), activation au moment de la rupture.**
- Guo 2014 (thèse, §2.3.3.2, pp. 72-74, éqs. 2.39-2.46) `[LU]`. Quand un joint N1N2N3-N4N5N6 casse, le
  couple face-à-face est formé (éq. 2.39). Pour CHACUN des 6 nœuds du joint, le groupe de **tous les
  tétraèdres qui contiennent ce nœud** est ensuite construit :
  $\mathrm{Group}_k=\{\mathrm{tet}_i \mid N_k\in\mathrm{tet}_i\}$, $k=1..6$ (éqs. 2.40-2.45). La détection
  est lancée « in and between these six groups » (éq. 2.46). Exclusion : « contacting couples that still
  have joint elements connecting tetrahedral elements are excluded ». Le texte précise : « node-node
  contact, edge-edge contact, edge-face contact ». **C'est exactement l'anneau de sommets** proposé dans
  l'enquête (§4-1a), activé à la mort du joint.
- Même procédure, en version publiée : Guo, Xiang, Latham & Izzuddin (manuscrit UCL, local
  `Manuscript_UCL_deposit.pdf`, §« contact », l. 488-530 du texte) `[LU]` : « the contact detection in
  the continuum region (no fractures) is only activated after new factures are formed ».
- Dans ce modèle, le joint intact porte lui-même la compression (Guo 2014, p. 64 : « This compressive
  stress in the joint element increases as the adjacent tetrahedral elements penetrate into each other
  … which functions similarly to the contact force ») `[LU]`.

**FDEM adaptative (CZM extrinsèque) : MultiFracS (Yan), codes GPGPU de Fukuda et Mohammadnejad.**
- Yan et al. 2023 (IJRMMS 169, 105439, §2.2, p. 5) `[LU]` : « Both intrinsic and extrinsic cohesive zone
  model-based approaches need to **carefully consider the timing for the activation of contact force
  calculation** as indicated by Fukuda et al. [70]. In this paper, once a cohesive element is dynamically
  inserted, the contact between the triangular elements connected by this cohesive element is activated. »
  ⚠️ Pris à la lettre, Yan n'active que la paire d'arête, donc PAS les voisins par le sommet. C'est le défaut
  de rockim.
- **Fukuda et al. 2021**, « Modelling of dynamic rock fracture process using the FDEM with a novel and
  efficient contact activation scheme », IJRMMS 138, 104645 `[RÉSUMÉ]`. Points du résumé :
  1. l'ACAA de Mohammadnejad et al. (2020, CPM) n'active comme candidats que les éléments de bord et ceux
     voisins des cohésifs **rompus**. Elle produit un « **spurious fracturing mode** » (fragmentation non
     physique) aux barres de Hopkinson ;
  2. « it is proven that spurious fracturing is caused by **unphysical movements of mesh due to missing
     contacts during shear softening of cohesive elements**, rather than from sudden activation of contact
     force calculations » ;
  3. remède (semi-ACAA) : activer le contact des éléments AUTOUR des cohésifs dès que leur fonction
     d'adoucissement en cisaillement franchit un seuil ; « when the threshold is set to around unity, the
     spurious fracturing mode can be overcome » ;
  4. l'activation brute de tous les éléments (BCAA) est jugée coûteuse, et même « physically
     unreasonable » pour Mohammadnejad 2020 `[RÉSUMÉ]`.

  **C'est le diagnostic de l'enquête, publié** : contacts manquants pendant le glissement ou l'écrasement
  des joints vivants, d'où mouvements non physiques et fragmentation parasite. Non lu : la définition exacte
  de « around » (arête, anneau de sommets, rayon ?) et de la grandeur seuil.
- HOSS (Knight et al. 2020, CPM 7:765, §2.2-2.3, p. 768) `[LU]` : SCA (smooth contact) couplé à l'UCZM,
  « concentrates the efforts of resolving contact only around the areas of interest, i.e., where fracture
  processes within the material are ongoing ». Le critère géométrique n'est pas donné dans ce papier.

### 1.2 Solution recommandée pour rockim (défaut a)

La littérature converge sur ceci : **la paire doit être candidate AVANT de pouvoir se recouvrir**. L'activation
doit avoir lieu au début de l'endommagement, pas à la mort du joint (Fukuda 2021), et elle doit porter sur
l'**anneau des sommets**, pas seulement sur la paire d'arête (Guo 2014). Seules les paires reliées par un
joint qui porte lui-même la compression restent exclues.

Définition du jeu actif (2D ; en 3D, remplacer triangle par tétraèdre) :

$$
\mathcal{C}=\Big\{(A,B)\ :\ A\in\mathcal{E}_{\text{act}}\ \lor\ B\in\mathcal{E}_{\text{act}},\ \ \mathrm{bbox}(A)\cap\mathrm{bbox}(B)\neq\emptyset,\ \ \nexists\,J\ \text{vivant et lié entre }A\text{ et }B\Big\}
$$

$$
\mathcal{E}_{\text{act}}=\mathcal{E}_{\partial}\ \cup\ \bigcup_{J\,:\,D_J>D^\star\ \text{ou}\ J\ \text{inséré}}\ \ \bigcup_{n\in\mathrm{nœuds}(J)}\mathrm{Ring}(n),\qquad \mathrm{Ring}(n)=\{E\mid n\in E\}
$$

Ici $D^\star$ est le seuil d'activation. Avec un CZM extrinsèque comme celui de rockim, prendre l'insertion
elle-même comme déclencheur ($D^\star=0$), ce qui revient au « seuil ≈ 1 » de Fukuda sur la fonction
d'adoucissement. `Ring` est calculé sur la topologie INITIALE (Guo : « same initial coordinates »), parce que
la scission des sommets duplique les nœuds.

Pseudo-code (FdemSolver.cpp, à la place de `:7134-7139` ; même chose en 3D) :

```cpp
// à l'insertion d'un joint (activateJoint) — PAS à sa mort :
for (int n : jointNodesInitial(jI))            // 4 nœuds en 2D (6 en 3D), topologie initiale
    for (int e : ringInitial[n])               // tous les éléments contenant ce sommet
        if (!inAct[e]) { inAct[e] = 1; actElems.push_back(e); }   // persistant, événementiel
// recherche de paires (chaque gcActEvery pas) : elems = actElems (+ bord) ;
// exclusion inchangée : paire liée par un joint VIVANT qui porte la compression.
// supprimer le rafraîchissement « stepCount_ % 8 » de deadList_ (:7421) : la mise à jour est faite
// au fil de l'eau, à l'insertion.
```

Coût : l'anneau compte environ 6 éléments par sommet, soit quelques dizaines d'éléments par joint inséré.
Ce coût est borné par la zone endommagée (ACAA/semi-ACAA), bien inférieur à la BCAA (tous les éléments).

Garde-fou (diagnostic, recommandé) : à la naissance d'une paire, si $|S_0| > \alpha\, v_{\max}\,\Delta t\,L$,
compter $e(S_0)$ dans un poste « énergie de naissance » et le neutraliser (§3). Après correction de la
détection, ce poste doit rester voisin de zéro.

---

## 2. (c) Forces nodales du potentiel : gradient exact pour des triangles déformables

### 2.1 Formule de Munjiza
- Potentiel (Munjiza 2004 ; écrit chez Tatone & Grasselli 2015, éq. 8, p. 58 `[LU]`) :
  $\varphi_i(P)=p_n\min\{3A_1/A,\,3A_2/A,\,3A_3/A\}$ si $P\in E_i$, 0 sinon. Autrement dit
  $\varphi_A(x)=\min_a 3\lambda_a(x)$ (normalisé, $p$ mis en facteur), et vaut 1 au centroïde.
- Force (Tatone éqs. 3-4 et 7, p. 57) `[LU]` :
  $\mathbf f_c=\int_{S}(\nabla\varphi_c-\nabla\varphi_t)\,dA$. Forme de bord (Xiang, Munjiza, Latham &
  Guises 2009, Eng. Comput. 26:673, éqs. 6-7, p. 676) `[LU]` :
  $\mathbf f_c=\sum_i\sum_j\oint_{\Gamma_{\beta_{c_i}\cap\beta_{t_j}}}\mathbf n\,(\varphi_{c_i}-\varphi_{t_j})\,d\Gamma$,
  et par arête $\mathbf f_{c,\text{edge}}=\frac{1}{u^2}\,\mathbf u\!\int_0^{L}p\,\varphi(v)\,dv$ (Xiang
  2009, éq. 7 ; la normale est obtenue en tournant le vecteur d'arête $\mathbf u$). Répartition sur les
  nœuds : la force sur l'arête va aux deux nœuds de l'arête, la réaction va aux nœuds du triangle cible
  (Tatone p. 58, fig. 1e) `[LU]`.
- Munjiza & Andrews 2000 annoncent un algorithme « incorporating contact kinematics **preserving energy
  balance** » `[RÉSUMÉ]`. Cette preuve repose sur des MOUVEMENTS RIGIDES : $\mathbf f$ est l'opposé du
  gradient de $E$ par rapport à la translation relative (l'enquête le confirme, travail sur boucle rigide
  = 2e-16). Non vérifié sur le texte : si Munjiza 2004 (§2) traite le cas déformable.

### 2.2 Gradient exact `[DÉRIVÉ]`
Énergie : $E=p\int_{S}(\varphi_A+\varphi_B)\,dA$ avec $S=A\cap B$. Le potentiel est matériel
($\varphi_A=\psi(\lambda^A(x))$, transporté de façon affine par $A$). On fait varier le nœud $\mathbf x_a$ de
$A$, à $B$ fixe. Le transport de Reynolds donne, avec $\mathbf v_A(x)=\sum_a\lambda_a\dot{\mathbf x}_a$,
$\partial_t\varphi_A=-\nabla\varphi_A\cdot\mathbf v_A$, et $\varphi_A=0$ sur $\partial A$ :

$$
\frac{\partial E}{\partial\mathbf x_a}=p\Big[\oint_{\partial A\cap B}\varphi_B\,\lambda_a\,\mathbf n_A\,ds\;-\;\int_S\lambda_a\nabla\varphi_A\,dA\Big]
$$

On intègre par parties ($\nabla\lambda_a$ est constant sur $A$ et $\varphi_A=0$ sur $\partial A$) :

$$
\boxed{\;\mathbf f_a=-\frac{\partial E}{\partial\mathbf x_a}
=\underbrace{p\Big[-\oint_{\partial A\cap B}\varphi_B\,\lambda_a\,\mathbf n_A\,ds+\oint_{\partial B\cap A}\varphi_A\,\lambda_a\,\mathbf n_B\,ds\Big]}_{\mathbf f_a^{\text{Munjiza}}}\;-\;p\,I_A\,\nabla\lambda_a,\qquad I_A=\int_S\varphi_A\,dA\;}
$$

La formule est symétrique pour les nœuds $b$ de $B$ ($-p\,I_B\nabla\lambda_b$, $I_B=\int_S\varphi_B\,dA$).
Avec $\nabla\lambda_a=-\boldsymbol\nu_a/h_a=-\ell_a\boldsymbol\nu_a/(2|A|)$, où $\boldsymbol\nu_a$ est la
normale sortante et $\ell_a$ la longueur de l'arête opposée à $a$. En 3D, même formule :
$\varphi=\min 4\lambda_a$, $I_A=\int_V\varphi_A\,dV$, $\nabla\lambda_a=-\mathcal A_a\boldsymbol\nu_a/(3V)$.

Propriétés : $\sum_a\nabla\lambda_a=0$ et $\sum_a\mathbf x_a\times\nabla\lambda_a=0$. Le terme correctif est
donc auto-équilibré (force totale et moment inchangés : la 3e loi tient toujours) et ne travaille que dans les
modes de déformation (dilatation de $A$ : travail $p\,I_A\,d\ln|A|$). Ordre relatif : $\delta/h$. Ce résultat
coïncide avec la vérification aux différences finies de l'enquête (écart 5e-10, boucle 1e-16).

**Référence publiée** : aucune trouvée qui écrive ce terme pour des triangles déformables. Pistes non lues :
- **Cai, Gao, Ai & Zhi 2024**, « A 2D energy-conserving contact model for the combined finite-discrete
  element method (FDEM) », Comput. Geotech. 166, 105972 `[RÉSUMÉ]` : titre et extrait seulement (énergie de
  contact fonction de la région recouverte). C'est la première chose à lire ;
- Feng 2021, CMAME 373:113493, théorie ECC (« energy-conserving contact ») `[RÉSUMÉ]` : corps rigides
  uniquement ;
- Lei, Rougier, Euser & Munjiza 2020, « A smooth contact algorithm for FDEM » (CPM) `[RÉSUMÉ]` : la force
  est le gradient d'un potentiel nodal lissé, mais je n'ai pas vu de preuve de conservation en déformable.

### 2.3 Implémentation (PotentialContact.hpp::pairForce)
```cpp
// après le calcul Munjiza fA[3], fB[3] (inchangé) :
double IA = integralPhi(S, triA), IB = integralPhi(S, triB);  // découpe de S par les médianes
                                                            // (φ affine par morceaux : 3 sous-triangles par élément)
for (a=0..2) { Vec2 g = gradLambda(triA, a);  fA[a] -= p*sc*IA*g; }
for (b=0..2) { Vec2 g = gradLambda(triB, b);  fB[b] -= p*sc*IB*g; }
// test : boucles fermées sur les ddl nodaux (loop2) -> travail < 1e-12 ; selftest-potential2d
// étendu à des triangles déformés (aujourd'hui rigides seulement).
```
Priorité faible : W_extra pèse 0,3 % de l'injection mesurée. La correction est néanmoins nécessaire pour que
le bilan B4 ferme à la précision machine.

---

## 3. (b) Naissance du contact et recouvrement initial

### 3.1 Solutions publiées
| Méthode | Équation | Énergie créée | Source |
|---|---|---|---|
| Rampe de pénalité à géométrie fixe (rockim `gcBirth=ramp` ; Abaqus « gradual interference fit ») | $f=s(t)\,f(S)$, $s:0\to1$ | $\int e\,ds=e(S_0)$, quelle que soit τ : **injecte toujours** | Abaqus : le recouvrement est « resolved … gradually », il faut garder l'énergie cinétique faible « by applying the fit slowly » `[RÉSUMÉ]` ; enquête §2.1-4 |
| Déplacement sans déformation | on déplace les nœuds esclaves pour annuler le recouvrement | 0, mais la géométrie est modifiée | Abaqus « strain-free adjustment » ; LS-DYNA `IGNORE=0` « move the penetrating node outwards » `[RÉSUMÉ]` |
| **Offset figé avec cliquet** | $F=k\,(d-d_0)$ si $d>d_0$ ; sinon $d_0\leftarrow d$, $F=0$ | 0 (force continue, énergie continue) | LS-DYNA `IGNORE=1`, *CONTROL_CONTACT carte 4 : « the initial penetration (d0) is stored … F = k·(d − d0) … If the current penetration is less than the reference penetration, the reference penetration is updated » `[RÉSUMÉ]` (lsdyna.ansys.com, « Initial Penetrations ») |
| Continuité de raideur cohésif → contact | $k^-=k^+(D)$ | supprime $\Delta H=\tfrac12(\alpha-1)k_0\delta_0\delta_1$ par traversée de δ = 0 | Ghesquière-Diérickx, Molinari & Anciaux 2025, arXiv:2511.14323 `[VAULT]`. Ils la jugent « outil de diagnostic, pas remède définitif », car elle autorise plus d'interpénétration |
| Contact non régularisé | Newmark non lisse semi-explicite : β = 0, γ = ½ pour la partie lisse, Euler implicite sur les impulsions de contact | pas de pénétration stockée ; restitution contrôlée par un coefficient e | Ghesquière-Diérickx et al. 2026, arXiv:2606.01355 `[RÉSUMÉ]` ; la conclusion de 2511.14323 est que le contact par pénalité « n'est pas viable » pour la fragmentation à long terme `[VAULT]` |
| Activation précoce (évite la naissance profonde) | voir §1.2 | ≈ 0 | Fukuda 2021 `[RÉSUMÉ]`, Guo 2014 `[LU]` |
| « Re-scaling » de pénalité de Y3D, « contact birth » de Munjiza | — | — | **NON TROUVÉ** : aucune trace dans les sources accessibles. Je ne peux ni confirmer ni décrire ces mécanismes. |

### 3.2 Recommandation pour rockim
1. **Cause** : la détection selon §1.2. Les naissances se font alors à $|S_0|\sim v\,\Delta t\,L$ et la rampe
   devient inutile (`gcBirth` = pénalité pleine).
2. **Filet de sécurité, version conservative de l'offset LS-DYNA adaptée au potentiel**. Le $e(S)$ de Munjiza
   n'est pas quadratique en une pénétration scalaire. On décale donc le NIVEAU du potentiel (équivalent
   exact de $d-d_0$) :

$$
E_{\text{eff}}=p\int_{S}\big(\varphi_A+\varphi_B-c_0\big)_+\,dA,\qquad c_0^{\text{naissance}}=\max_{S_0}(\varphi_A+\varphi_B),\qquad c_0\leftarrow\min\big(c_0,\max_{S}(\varphi_A+\varphi_B)\big)
$$

   - $E_{\text{eff}}(S_0)=0$ : aucune énergie n'est créée à la naissance. Force continue : la région active
     croît depuis l'aire nulle. Le cliquet n'agit que quand $E_{\text{eff}}=0$, donc sans saut. $c_0$ décroît
     jusqu'à 0 quand la paire se sépare : on retrouve alors le contact standard ;
   - le domaine $\{\varphi_A+\varphi_B>c_0\}$ est un polygone : $\varphi_A+\varphi_B$ est affine sur chaque
     cellule de la découpe par les médianes, et on obtient le domaine par un découpage de plus ;
   - version simple (non lisse) : $E_{\text{eff}}=\max(0,\,e(S)-e_0)$ avec cliquet $e_0\leftarrow\min(e_0,e(S))$.
     L'énergie est exacte, mais la force saute de $0$ à $|\nabla e(S_0)|$ au franchissement. Erreur O(Δt) par
     franchissement, du même type que le saut de raideur de Ghesquière : acceptable seulement si ces
     naissances sont rares (après 1.) ;
   - ⚠️ la forme actuelle `sc = 1 − aRef/A` reste non conservative même avec `aRef` figé : la force appliquée
     est $sc\,\nabla e$ alors que le gradient de $sc\cdot e$ contient en plus $e\,\nabla sc$. Si on garde ce
     facteur, il faut figer `aRef` (pas de relaxation à géométrie fixe) ET ajouter $-e\,(a_{\text{ref}}/A^2)\nabla A(S)$.
3. **Bilan** : poste « énergie neutralisée à la naissance » $=\sum e(S_0)$, imprimé. Avertissement si
   `gcBirthTau < dt` (la rampe se réduit à une marche).
4. **Transition joint → contact** (`gcBirth=relay`) : imposer la continuité de raideur $k^-=k^+(D)$ au moment
   où le joint meurt (Ghesquière 2025 `[VAULT]`). Sinon, chaque oscillation autour de δ = 0 pompe
   $\tfrac12(\alpha-1)k_0\delta_0\delta_1$.

```cpp
// naissance d'une paire (isNewH) :
H.c0 = (R.area > areaTol) ? maxPhiSum(S) : 0.0;   // niveau d'offset
birthNeutralized_ += energy(S, /*c0=*/0) ;         // e(S0), pour le bilan
// à chaque pas :
H.c0 = std::min(H.c0, maxPhiSum(S));               // cliquet (ne remonte jamais)
forces = pairForceLevel(S, H.c0);                  // ∫ sur {φA+φB > c0}, + terme -p I ∇λ (§2)
```

---

## 4. Stabilité : pas critique avec pénalité, choix de p, ordre de l'erreur d'énergie

**Pas critique.** Pour la différence centrale : $\Delta t\le 2/\omega_{\max}$. La borne de Gershgorin est
conservative : $\omega_{\max}^2\le\max_i \sum_j|K_{ij}|/m_i$ (O'Sullivan & Bray 2004, Eng. Comput. 21:278
`[MÉMOIRE]`). Raideur de contact Munjiza par unité de longueur d'arête qui pénètre de δ dans un triangle de
hauteur $h$ : $\varphi\approx3\delta/h$, d'où $k_n\approx 3p/h$ par potentiel (enquête : environ
$3pL/h$ par nœud). Règle pratique, nœud de masse $m$ engagé dans $n_c$ contacts :

$$
\Delta t\le\frac{2}{\sqrt{\big(k_{\text{el}}+n_c\,3pL/h_{\min}\big)/m}}\times\text{facteur de sécurité}
$$

Guo 2014 (§2.3.5, pp. 78-79, éqs. 2.56-2.60) `[LU]` :
$\Delta t=\min\{\Delta t_{\text{FEM}}\sim0.1h\sqrt{\rho/E},\ \Delta t_{\text{DEM}}\sim\tfrac{\pi}{5}\sqrt{m/k}\}$
(Tsuji 1993), avec $k=E\,h$, d'où $\Delta t_{\text{DEM}}\approx0.2h\sqrt{\rho/E}$.

**Choix de p** (sources lues) :
| Source | Valeur |
|---|---|
| Mahabadi 2012, repris par Lisjak 2013 (thèse §4.4, p. 71) `[LU]` | $p_n=10E$, $p_t=1E$, $p_f=5E$ |
| Tatone & Grasselli 2015, §3.4, p. 64 `[LU]` | procédure : balayer chaque pénalité sur 3-4 ordres de grandeur en UCS et retenir le plateau ; mortier : $100E$, « the largest possible values that did not necessitate a reduction in the time step » |
| Guo 2014, éq. 2.28, p. 67 `[LU]` | joint : $E\le p_0\le10E$ |
| Deng et al. 2022 (EFM 276) `[LU]` | joint $P_f=30E$ ; $P_n$ calibré selon Deng 2021 (EFM 242:107459) |
| Ghesquière 2025 `[VAULT]` | $k^-\approx\alpha E/h$, $\alpha\sim10^1$-$10^4$, référence 10 E/h |
| Munjiza 2004, « p = 10 à 100 E » | **NON VÉRIFIÉ** (livre inaccessible) |

**Ordre de l'erreur d'énergie** `[MÉMOIRE]`. Le saute-mouton (Verlet) est symplectique : pour une force lisse,
l'énergie modifiée est conservée et l'erreur reste bornée en O(Δt²) sans dérive (Hairer, Lubich & Wanner,
*Geometric Numerical Integration*, ch. IX). Le contact par pénalité est seulement C⁰ à l'ouverture et à la
fermeture (raideur discontinue). Chaque franchissement coûte alors une erreur O(Δt) de signe non contrôlé,
et ces erreurs s'accumulent avec le nombre de basculements. Ghesquière 2025 `[VAULT]` donne la forme exacte
de ce saut pour un ressort à raideur discontinue : $\Delta H=\tfrac12(\alpha-1)k_0\delta_0\delta_1$. En
pratique, ils recommandent $\Delta t\lesssim0.2\,\Delta t_c$ quand les commutations sont fréquentes.
L'enquête mesure un compteur de biais $D\approx-\tfrac12\sum\Delta x\,K\,\Delta x<0$, proportionnel à Δt :
c'est cohérent. Le défaut (a) est d'une autre nature, O(1) et indépendant de Δt : la naissance profonde.

---

## 5. (d) Bilan d'énergie avec vitesses imposées (rouleaux) en saute-mouton

Identité exacte du saute-mouton : $\mathbf v^{n+1/2}=\mathbf v^{n-1/2}+\Delta t\,M^{-1}\mathbf F^n$. Pour
chaque ddl :

$$
\tfrac12M\big(v^{n+1/2}\big)^2-\tfrac12M\big(v^{n-1/2}\big)^2=\Delta t\,F^n\,\bar v^{\,n},\qquad \bar v^{\,n}=\tfrac12\big(v^{n-1/2}+v^{n+1/2}\big)
$$
$$
=\Delta t\,F^n v^{n-1/2}+\frac{(\Delta t\,F^n)^2}{2M}\quad(\text{terme « biasW_ » de rockim})
$$

Cette identité n'est vraie que si $\mathbf F^n$ est la force qui produit réellement le saut de vitesse. Sur un
ddl bloqué, la vitesse n'est pas modifiée ($\Delta v=0$) : la force effective y est
$F^n_{\text{eff}}=F^n+R^n=M\,\Delta v/\Delta t=0$, où $R^n$ est la réaction du rouleau. Ordre correct
(Belytschko, Liu & Moran, *Nonlinear FE for Continua and Structures*, §6.2 « energy balance » : les
réactions des ddl imposés font partie de $W_{\text{ext}}$ `[MÉMOIRE]`) :

1. assembler $\mathbf F^n$ (interne + contact + joints + charges) ;
2. **appliquer les contraintes à la force** : sur un ddl bloqué, $R^n=-F^n+M\,a^n_{\text{imp}}$ (avec
   $a_{\text{imp}}=0$ pour un rouleau fixe), puis $F^n\leftarrow M a^n_{\text{imp}}$ ;
3. calculer les postes de bilan (biais $(\Delta tF)^2/2M$, travaux) avec la force CONTRAINTE ;
   $W_R\mathrel{+}=\Delta t\,R^n\bar v^{\,n}$, qui vaut 0 pour un rouleau immobile ;
4. mettre à jour $v$, puis $x$.

Mieux, et sans biais : compter chaque poste au point milieu, $W_k\mathrel{+}=\Delta t\,\mathbf F_k^n\cdot\bar{\mathbf v}^{\,n}$,
avec $\bar{\mathbf v}^{\,n}$ calculé APRÈS l'application des contraintes. Alors
$\sum_kW_k=\Delta KE_{1/2}$ à la précision machine : le terme `biasW_` disparaît et chaque poste devient
lisible séparément (déjà proposé dans l'enquête, §4-4).

```cpp
// integrate(), chemin groupé (:9148-9163) et non groupé (:9258-9273) :
Vec2 F = Fint[i] + Fc[i] + Fj[i] + Fext[i];
if (flag_[i0] == ROLLERX) { R.x() = -F.x(); F.x() = 0.0; }   // AVANT tout comptage
if (flag_[i0] == ROLLERY) { R.y() = -F.y(); F.y() = 0.0; }
Vec2 vOld = v[i], vNew = vOld + dt*F/M;
Vec2 vBar = 0.5*(vOld + vNew);
Wc += dt*Fc[i].dot(vBar); Wj += dt*Fj[i].dot(vBar); ...;   // Fc, Fj : composantes bloquées
                                                            // à exclure (vBar.x = 0 sur ROLLERX : automatique)
Wroll += dt*R.dot(vBar);                                    // = 0 pour un rouleau fixe
v[i] = vNew;
// test : smoke T = 0,03 s -> résidu B4 = 0 à 1e-9 près (aujourd'hui -136 005 J/m).
```

---

## 6. Synthèse par défaut

| Défaut | Solution recommandée | Référence principale | Statut |
|---|---|---|---|
| (a) paires à sommet commun invisibles | candidats = anneau des sommets de tout joint inséré ou endommagé, activé à l'insertion ou au début de l'adoucissement ; dans Y, tous les éléments sont candidats | Guo 2014 éqs. 2.39-2.46 ; Fukuda et al. 2021 (semi-ACAA, « missing contacts during shear softening ») ; Y-Geo « all pairs » | Guo, Y-Geo `[LU]` ; Fukuda `[RÉSUMÉ]` |
| (b) rampe gcBirth (τ < dt), énergie de naissance | supprimer la rampe à géométrie fixe ; offset de niveau $c_0$ avec cliquet (analogue de LS-DYNA IGNORE=1) ; $k^-=k^+(D)$ au relais joint → contact | LS-DYNA IGNORE=1 ; Ghesquière 2025 | `[RÉSUMÉ]` / `[VAULT]` ; adaptation au potentiel `[DÉRIVÉ]` |
| (c) forces non gradient | $\mathbf f_a=\mathbf f_a^{\text{M}}-p\,I_A\nabla\lambda_a$ | aucune publiée trouvée ; à vérifier dans Cai et al. 2024 | `[DÉRIVÉ]` + contrôle numérique de l'enquête |
| (d) biais leapfrog avant le rouleau | contraindre F avant le comptage ; travail au point milieu $\bar v$ | BLM §6.2 | `[MÉMOIRE]` + identité exacte ci-dessus |

## 7. Non vérifié (à faire si l'accès revient)
- Munjiza 2004 (livre) : le §2 traite-t-il le cas déformable ? valeur recommandée de p ? existe-t-il une
  « contact birth » ? **Non lu.**
- Munjiza & Andrews 2000 (IJNME 49:1377) : connu par le résumé seulement (« preserving energy balance »).
  Munjiza & Andrews 1998 (NBS) : non lu.
- Munjiza, Knight & Rougier 2011 (*Computational Mechanics of Discontinua*) et 2015 (*Large Strain FDEM*) :
  non lus.
- Mahabadi et al. 2012 (Y-Geo, IJG) : connu par la citation de Lisjak seulement.
- Irazu (Mahabadi 2016), Solidity : pas de description du traitement des paires adjacentes trouvée (le
  papier Solidity local, DEM7 2016, ne donne que la formule de force).
- « Re-scaling Y3D » et « initial interference removal » propres à la FDEM : rien trouvé.
- Fukuda 2021, Mohammadnejad 2020 : définition exacte du voisinage activé et du seuil (résumés seulement).
- Ghesquière 2025 (2511.14323) : repris des notes du vault, non relu aujourd'hui. 2606.01355 : résumé seulement.
- Cai et al. 2024 (2D energy-conserving FDEM contact) : titre et extrait seulement. **C'est la source à lire
  en priorité pour (c).**
- BLM §6.2 et Hairer-Lubich-Wanner : de mémoire.
