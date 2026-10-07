# Enquête : injection d'énergie par le contact potentiel (FDEM 2D, tunnel `tip16`)

Date : 2026-10-07. Lecture seule du dépôt `/home/user/rockim` : aucun fichier n'y a été modifié.
Les preuves viennent d'une COPIE instrumentée du code, compilée dans le scratchpad
(`enquete_contact/src_copy`, binaire `enquete_contact/bld/rockim`). L'instrumentation ne fait que lire :
le `gcWork_` imprimé est identique au bit près à celui du binaire d'origine. Elle s'appuie sur deux bancs
autonomes (`loop2.cpp`, `blocks.cpp`) qui incluent le `PotentialContact.hpp` du dépôt tel quel.
Tous les calculs : 2 fils, `nice -n 10`. Maillage : `meshes/tunnel_hs_red.msh` du dépôt (22 730
triangles, l'original, dt = 5,83e-6 s). La campagne diag_contact avait utilisé le maillage régénéré de
23 004 triangles : les chiffres absolus diffèrent, la signature non.

## 0. Résultat en une phrase

Presque toute l'injection vient d'une seule source. Des paires de triangles qui ne partagent qu'un
SOMMET se recouvrent en profondeur tant que le contact ne les voit pas, parce qu'aucun des deux
n'appartient encore au jeu actif. Le jour où l'un d'eux y entre, la paire naît d'un coup avec ce
recouvrement : le potentiel de Munjiza crée alors l'énergie p∫(φA+φB)dA sans qu'aucun travail ne l'ait
payée, puis la restitue en écartant les blocs. La rampe `gcBirth` ne peut rien contre ce saut, parce que
gcBirthTau (1 µs) est plus court que dt (5,8 µs). Cette énergie est proportionnelle à p à recouvrement
donné. Le résidu B4 est un problème indépendant : c'est la correction leapfrog comptée sur les ddl bloqués
des rouleaux.

## 1. Mesure : décomposition exacte du travail normal du contact (copie instrumentée)

Pour chaque paire, à chaque pas, on calcule l'énergie EXACTE de la fonctionnelle de Munjiza,
e = p ∫_S (φA + φB) dA (intégrale exacte par découpe de S selon les médianes des deux triangles,
`potx.hpp`), multipliée par le facteur `sc` réellement appliqué. Le travail compté
`Wn = Σ f·v_old dt` (la ligne 7299) se décompose alors exactement en :

  Wn = −ΣdE_géom + W_extra + D,   avec   −ΣdE_géom = naissances + ΣdE_sc − disparitions − U_fin

- dE_géom : variation de e à sc fixé, c'est-à-dire le travail conservatif ;
- dE_sc : variation de sc à géométrie fixée, c'est-à-dire de l'énergie créée ou détruite sans travail ;
- W_extra : travail de la part non gradient des forces nodales de Munjiza (§3.2) ;
- D : erreur d'intégration en temps du compteur.

Contrôle de fermeture : 1e-6 à 1e-9 J/m sur tous les calculs.

| run (T) | dt (s) | Wn = gcWork_ | ΣdE_sc (naissance en ≤ 3 pas) | dont paires à sommet commun, élément entrant dans le jeu actif | naissances de ce type / aire moy. | max par paire | W_extra | D | U_fin |
|---|---|---|---|---|---|---|---|---|---|
| `tip16_mu0` (0,25 s) | 5,83e-6 | **8,20e6** | **9,09e6** | 9,09e6 | 4 977 / 4,6e-6 m² | — | −2,6e4 | −7,6e5 | 1,1e5 |
| `tip16_mu0_gcpen` (0,25 s) | 5,83e-6 | 9,60e6 | 0, mais 1,09e7 comptés dans « naissances autres paires » | — | 5 031 / 5,1e-6 m² | — | −3,4e4 | −1,10e6 | 1,7e5 |
| `tip16_mu0` (0,12 s) | 5,83e-6 | 1,31e6 | 1,452e6 | **1,452e6 (100 % à sommet commun)** | 1 281 / 2,9e-6 m² | 3,2e5 J/m | −4,4e3 | −1,12e5 | 2,7e4 |
| `tip16_mu0_dt05` (0,12 s) | 2,92e-6 | 1,13e6 | 1,200e6 | 1,200e6 (100 %) | 1 175 / 2,6e-6 m² | 2,5e5 J/m | −4,2e3 | −4,5e4 | 1,8e4 |
| `tip16_mu0_pot01` (0,12 s) | 5,83e-6 | 1,32e5 | 1,79e5 | 1,79e5 (100 %) | 1 082 / 3,6e-6 m² | 3,9e4 J/m | −2,0e2 | −8,3e3 | 3,8e4 |

En comparaison, les naissances « normales » (paires déjà visibles qui commencent à se recouvrir,
environ 1e5 par run) ont une aire moyenne de 1e-7 à 3e-7 m² et créent moins de 1e3 J/m au total. Les
paires de joint mort pèsent 1e3 à 2e4 J/m. Le relais joint vers contact n'est donc pas la source.
Journaux : `enquete_contact/runs/{mu0b,gcpen,tip16_mu0_t12,tip16_mu0_dt05_t12,tip16_mu0_pot01_t12}.log`.

## 2. Cause principale : naissance à recouvrement profond des paires « à sommet commun »

### 2.1 Mécanisme (fichier:ligne, dépôt d'origine)

1. **Jeu de candidats** : `FdemSolver.cpp:7134-7139`. Les éléments candidats sont ceux qui possèdent au
   moins une arête dans `act_` :
   ```cpp
   for (const auto& be : act_)
       if (emark[be.elem] != epoch) { emark[be.elem] = epoch; elems.push_back(be.elem); }
   ```
   `act_` contient le pool extérieur activé plus les faces des joints MORTS (`deadList_`,
   `:7421-7426`, rafraîchi seulement quand `nBroken_` change et que `stepCount_ % 8 == 0`). Un triangle
   intérieur dont tous les joints sont vivants (liés, ou insérés et cassés mais pas morts) n'est jamais
   candidat.
2. **Exclusion** : `:7237`. Seules les paires reliées par un joint VIVANT sont exclues, puisque le joint
   porte leur interaction. Deux triangles qui ne partagent qu'un sommet n'ont pas de joint entre eux :
   AUCUNE loi ne gère leur recouvrement tant que l'un d'eux reste hors du jeu actif. Quand les joints
   insérés de l'éventail autour du sommet glissent ou s'écrasent (65 % meurent en compression, décharge
   `origin`, facteur de pointe), les coins de ces voisins par le sommet se chevauchent sans aucune force
   de rappel, sur des centaines de pas.
3. **Naissance** : `:7250-7283`. Quand un joint de l'éventail meurt, l'élément entre dans `elems`, et la
   paire est évaluée pour la première fois avec un recouvrement S0 accumulé (aire moyenne 4,6e-6 m², soit
   20 à 30 fois une naissance normale, avec des cas extrêmes à 3,2e5 J/m pour une seule paire) :
   ```cpp
   if (isNewH) H.aRef = R.area;            // sc = 0 au pas de naissance
   else H.aRef *= relax_;                  // relax_ = exp(-dt/gcBirthTau), :1428
   sc = std::max(0.0, 1.0 - H.aRef / R.area);
   ```
   Avec gcBirthTau = 1e-6 s et dt = 5,83e-6 s, relax_ = exp(−5,8) = 0,003 : sc passe de 0 à 0,997 en UN
   pas, à géométrie figée. L'énergie e(S0) = p∫_{S0}(φA+φB) apparaît sans travail (ligne ΣdE_sc), puis la
   paire la rend en se séparant. C'est ce que `gcWork_` compte comme « injecté ». Sous `gcBirth = penalty`
   (`:7262-7271`), une paire sans joint reçoit `fac = 1.0` (fnJ = 0), donc la pénalité pleine dès le pas
   de naissance : même énergie, comptée cette fois dans « naissances ». C'est pourquoi gcpen ne change
   presque rien. Sur ce maillage, il augmente même l'injection (9,6e6 contre 8,2e6), et le −22 % du
   diagnostic 2 est du même ordre que la variabilité de faciès.
4. La rampe elle-même, même lente (τ ≫ dt), ne serait pas neutre : faire croître sc à géométrie fixée
   crée ∫ e dsc = e(S0). Elle ne fait qu'étaler la restitution.

### 2.2 Cohérence avec les quatre diagnostics

- **∝ p** : S0 est fixé par la mécanique des joints et de l'éventail, pas par p (l'élément est hors
  contact pendant que le recouvrement grandit). Donc e(S0) ∝ p. Mesuré : 1,79e5 → 1,45e6 J/m de p = 1e9
  à 1e10 (×8,1, contre ×8,6 dans la campagne). L'écart à ×10 vient de recouvrements un peu plus grands à
  p faible (aire moyenne 3,6e-6 contre 2,9e-6 m²). Le ralentissement à ×4,5 entre 1e10 et 1e11 s'explique
  par le dt réduit et des recouvrements plus petits, qui naissent plus tôt quand le faciès devient plus
  violent.
- **Effet de dt (−41 %)** : pas de dépendance intrinsèque. À T = 0,12 s, dt/2 retire 17 % du terme de
  naissance (1,45e6 → 1,20e6), en raison de 8 % de naissances en moins (1 281 → 1 175) et d'une aire
  moyenne plus faible. Le compteur ajoute sa propre erreur D ∝ dt (§4). Le reste du −41 % sur 0,25 s
  vient d'un faciès différent (cohésif −37 % dans la campagne), et le nombre de naissances tardives suit
  la fragmentation. La disparition de la rampe joue aussi : relax_ = exp(−2,9) = 0,055 à dt/2, ce qui
  reste un saut.
- **Cliquet −16 %, `plastic` +114 %** : ces clés changent la cinématique des joints vivants de l'éventail
  (glissement ou écrasement avant la mort), donc S0. Elles agissent sur la cause sans en être une.
- **Le frottement masque** : il dissipe la restitution de e(S0) (glissement des coins sortants).

## 3. Pistes examinées

### 3.1 Piste 1, schéma d'intégration et cohérence du compteur : **contribution de signe opposé**
La force est bien évaluée à x_n (`step()` : forces puis `integrate()`), comme les autres postes. Tous
les postes comptent f·v_old·dt (point droit) et le biais global est rendu par `biasW_ = Σ|F|²dt²/2M`
(`:9148`). Le bilan B4 est donc cohérent poste par poste au second ordre près. Pour un canal raide,
l'erreur par canal vaut D ≈ −½ Σ Δx·K·Δx < 0. Mesurée : −7,6e5 J/m (base), divisée par 2,5 avec dt/2.
Le compteur fait donc paraître le contact plus dissipatif qu'il n'est. Banc `blocks` (rebond de deux
blocs déformables) : Wc = −840 J pour un KE de 5 300 J à p = 1e10, divisé par 2 avec dt, alors que
l'énergie mécanique exacte ne varie que de +0,17 %.
**Erreur de commentaire** : `FdemSolver.cpp:9914-9916` (« un petit résidu positif est le biais O(dt) du
compteur ») et le commentaire du selftest `:10620-10630` (« biais systématique POSITIF ») ont le signe
inverse. Le biais est négatif, si bien qu'un `gcWork_` positif est une injection réelle et même
sous-estimée.

### 3.2 Piste 2, le potentiel est-il un gradient : **défaut réel, contribution négligeable ici**
- La 3e loi et le moment sont exacts à la précision machine (`loop.cpp` : somme des forces 1e-16,
  moment 6e-17).
- Pour les mouvements RIGIDES relatifs, le champ est conservatif : travail sur boucle 2e-16. C'est
  l'énergie E = p∫_S(φA+φB)dA (dérivation : F_A = p∮(φA−φB)n se réécrit en −∂E/∂d).
- Pour des triangles DÉFORMABLES, il ne l'est pas. Les forces nodales de `pairForce`
  (`PotentialContact.hpp:180-214`, répartition λ_a(X) des charges de bord) valent
  −∂E/∂x_a **+ p (∫_S φA dA) ∇λ_a** pour les nœuds de A, et de même pour B. Ce terme en plus est
  auto-équilibré (Σ∇λ_a = 0, moment nul) : c'est une « pression » interne fictive, dont le travail vaut
  p·I_A·d(ln aire_A), sans potentiel. Vérification (`loop2.cpp`) : écart de 16 % de la force nodale
  maximale avec −dE/dq aux différences finies, et travail sur boucle de 2,9e-4 (en unités de p). La force
  corrigée « Munjiza − p I ∇λ » coïncide avec −dE/dq à 5e-10 près, et son travail sur boucle tombe à
  1e-16.
- Mesuré dans le tunnel : W_extra = −2,6e4 J/m, soit 0,3 % de l'injection. Son ordre de grandeur relatif
  est δ/h ≈ σ/(6p) ≈ 1e-4. Ce n'est pas la cause.
- Le selftest `selftest-potential2d` (`:10485-10650`) n'emploie que des corps RIGIDES : il est aveugle
  à ce défaut.
- La non-dérivabilité du min (cassures sur les médianes) ne crée pas de travail sur cycle (boucles
  rigides ≈ 0) : écartée.

### 3.3 Piste 3, détection et activation : **CAUSE (§2)**
Point précis : ni la fréquence de recherche ni la marge (`gcActEvery`, `gcActMargin`) ne sont en
cause. Le défaut est la définition du jeu de candidats PAR ARÊTES ACTIVES (`:7134`), qui laisse hors
contact les voisins par le sommet d'éléments intérieurs. Le cache `deadList_` rafraîchi « % 8 » sur
`nBroken_` (`:7421`) ajoute un retard, mais 100 % de l'énergie mesurée vient des paires à sommet commun,
aucune des paires de joint mort.

### 3.4 Piste 4, stabilité du pas : **écartée**
`computeStableDt()` budgète bien la raideur du potentiel : `:5222`
`if (contactPot_) kContact = std::max(kContact, std::max(potP_, potKt_));`, avec `nExtra = 2` (`:5202`).
La remarque de coupe_pdc (« kpGC_ seulement en SHPB ») est périmée pour le mode potential. La raideur
nodale réelle de Munjiza (≈ 3pL/h par nœud et par contact) peut dépasser 2p d'un facteur 2 à 4 sur des
éléments allongés, mais dtFactor = 0,2 donne ω·dt ≈ 0,4 à 0,8, loin de 2. Le pas est d'ailleurs commandé
par les ressorts de joint (jointPenaltyFactor = 20), sauf à p = 1e11. Banc `blocks` : énergie bornée
pour p de 1e9 à 1e11, aucune divergence.

### 3.5 Piste 5, amortissement de contact, rampe gcBirth, écrêtage : **la rampe est en cause, le reste est écarté**
- `potXi = 0` par défaut et `contactMu = 0` dans les runs `mu0` : aucun terme visqueux ni tangentiel dans
  le canal normal mesuré.
- La rampe (`:1428`, `:7281-7283`) : τ = 1 µs < dt, donc une seule marche (§2.1-3). Mesuré :
  ΣdE_sc = 9,09e6 J/m, pour 0,2 J/m après le 3e pas.
- `pairForce` rend `false` pour un triangle inversé : la paire disparaît, ce qui ne peut que perdre de
  l'énergie (disparitions mesurées : 1e3 J/m), donc écarté.

### 3.6 Piste 6, bilan B4 : **le résidu est la correction leapfrog des rouleaux, sans lien avec le contact**
`integrate()`, chemin groupé : `:9148`, `bias += F.squaredNorm()*dt*dt/(2M)` est calculé AVANT
`:9163`, `if (flag_[i0] == ROLLERX) vn.x() = 0.0;`. La part F_x²dt²/2M de l'énergie cinétique est
inscrite dans `biasW_` puis annulée par le rouleau : elle n'arrive jamais dans l'énergie cinétique.
Preuve : compteur `rollerBias` = Σ ½M(dt F_x/M)² sur les nœuds ROLLERX.

| run | résidu B4 | KE retirée par les rouleaux (mesurée) |
|---|---|---|
| smoke T = 0,03 s (pas de contact) | −136 005 | 136 005 |
| `tip16_mu0` T = 0,25 s | −947 419 | 947 419 |
| `tip16_mu0` T = 0,12 s | −527 664 | 527 664 |
| `tip16_mu0_dt05` T = 0,12 s | −263 954 | 263 954 |

L'égalité est exacte et ∝ dt (F_x²dt²/2M par pas, avec 1/dt pas). Ce terme vient de la réaction σ0 = 5
MPa sur les rouleaux latéraux : « 13,6 J/m par pas » d'après le commentaire de `red_tip16.cfg`. Il
explique la lecture 3 du diagnostic 3 (résidu divisé par 2 avec dt, −8,5e5 à dt × 0,62).

## 4. Correctifs PROPOSÉS (non appliqués) et tests de validation

1. **Contact des voisins par le sommet** (cause principale). Plusieurs options, à combiner :
   a. étendre `elems` (`:7134`) à l'anneau de sommets de chaque élément actif, ou y inclure tout
      élément portant un joint inséré non lié. Les paires à sommet commun naissent alors à recouvrement
      quasi nul ;
   b. naissance sans création d'énergie : si une paire naît avec une aire supérieure à un seuil
      (≈ v_max·dt·L), mémoriser e0 = e(S0) et ne faire agir que le potentiel EN EXCÈS, par exemple
      `aRef` figé (pas de relaxation) pour ces paires, ou retrait de e0 sous forme d'offset. Imprimer
      l'énergie ainsi neutralisée ;
   c. au minimum, une ligne dans le bilan : « énergie créée à la naissance des paires » (le ΣdE_sc
      instrumenté), et un avertissement quand gcBirthTau < dt, car la rampe dégénère en marche.
   **Test** : rejouer `tip16_mu0` (et `red_tip16`). Critères : ΣdE_sc + naissances < 1 % de l'énergie
   cohésive, et `gcWork_` ≤ 0 (contact dissipatif ou neutre) pour p = 1e9, 1e10 et 1e11. Un banc
   minimal : éventail de 6 triangles autour d'un sommet, joints insérés en glissement sous compression,
   puis mort d'un joint. Il faut que l'énergie mécanique KE + U_el + U_c ne saute pas à la mort.
2. **Résidu B4 des rouleaux** (`:9148/9163`, et le chemin non groupé `:9258` puis `:9273`) : appliquer la
   contrainte ROLLERX à F (F.x = 0) AVANT le calcul de `bias`, ou retrancher ½M(vn.x)² de `bias`.
   **Test** : le smoke T = 0,03 s doit donner un résidu nul à la précision machine (aujourd'hui
   −136 005 J/m).
3. **Forces nodales non gradient** (`PotentialContact.hpp:pairForce`) : retrancher p·I_A·∇λ_a (nœuds
   de A) et p·I_B·∇λ_b (nœuds de B), avec I = ∫_S φ dA calculée par la découpe selon les médianes
   (`potx.hpp:integrals/correct`, environ 40 lignes). **Test** : `loop2` (boucles sur ddl nodaux, travail
   inférieur à 1e-12), et étendre `selftest-potential2d` à des triangles DÉFORMABLES.
4. **Commentaires** `:9914-9916` et `:10620-10630` : corriger le signe du biais du compteur (négatif).
   Le cas échéant, compter le travail du contact en trapèze, ½(v_n + v_{n+1})·f_n·dt, pour un poste
   lisible sans biais O(dt).

## 5. Fichiers de l'enquête (scratchpad `enquete_contact/`)
- `potx.hpp` : énergie exacte p∫(φA+φB), I_A, I_B, force corrigée.
- `loop.cpp`, `loop2.cpp` : 3e loi, boucles fermées, comparaison avec −dE/dq.
- `blocks.cpp` : deux blocs CST déformables, rebond et pressage, Munjiza ou force corrigée, balayage de p
  et dt.
- `patch.py`, `patch2.py` et le bloc ajouté en ligne : instrumentation de la copie (`src_copy/`,
  `bld/rockim`).
- `runs/*.cfg|*.log` : decks (seuls `meshFile`, `frames` et `T` modifiés) et journaux.
