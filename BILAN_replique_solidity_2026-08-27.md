# Réplique Imperial College — bilan au 2026-08-27

Branche `joint-handoff` (worktree `rockim_p4`). Objectif : reproduire l'impact
à insert unique de Yang, Xiang, Naderi, Wang, Aising, Ugarte et Latham
(IJRMMS 191, 2025) **avec leur physique**, et non avec celle de rockim.

---

## 1. Le fait nouveau : leur solveur est public

`github.com/ImperialCollegeLondon/solidity-solver-open` — LGPL-3.0, C,
17 000 lignes, format `.Y3D`, poussé le 2026-03-31. C'est bien la lignée
Munjiza décrite par la thèse de Guo (2014) et les articles de Yang *et al.*

**Conséquence de méthode.** Jusqu'ici les conventions manquantes étaient
« déduites des articles ». Elles sont désormais **relevées dans le code**,
fichier et ligne cités. Plusieurs points réputés impossibles à répliquer faute
de source ne le sont plus.

Carte de lecture utile :

| fichier | contenu |
|---|---|
| `Y3Dfd.c` | `CauchyTet4` (loi de volume + viscosité), `Sigma_tau` l. 1078 (loi de joint), `S_N_direction` l. 878 (repère local, points d'intégration), `Yfd3TET4JOINT` l. 1382 (DIF) |
| `Y3Did.c` | `force_friction_cal` l. 899 (frottement, naissance d'un contact sur joint mort) |
| `Y3Dsd.c` | intégration en différences centrées — **sans aucun amortissement** |
| `Y3Drd.c` | lecture des propriétés matériau (`d1peks` = viscosité = propriété 0) |

⚠️ Piège de lecture : `d1nbp` / `dbkp` ressemblent à un mécanisme de naissance
de contact dans la signature de `Yfd3TET4JOINT`, mais c'est une **contre-pression
hydraulique** (20 MPa en dur dans `Y3Drd.c` l. 1422, utilisée par `Y3Dhydpre.c`).

---

## 2. Ce que la source CONFIRME — rockim était déjà juste

Rien à changer sur ces points :

- loi de volume `T = (μ/J)B + [(λ lnJ − μ)/J]I + η·D` (l. 716) → `bulkModel = neohookean` ;
- viscosité `dpeks·D`, et rockim écrit `2μD` → `bulkViscosity = η/2 = 2000` pour leur η = 4000 ;
- courbe z a = 0,63, b = 1,8, c = 6 (l. 1088-1090) ;
- couplage elliptique `SQRT(tmp1²+tmp2²)` (l. 1136) ;
- branche élastique parabolique `(2r−r²)` (l. 1274) et raideur **double** en compression (l. 1265) ;
- 3 points d'intégration aux **milieux d'arêtes**, poids A/6 sur chacun des deux nœuds (l. 900-917 et 940) ;
- DIF appliqué à f_t **et** G_I, c **et** G_II par le même facteur (l. 1448-1456) ;
- **aucun amortissement dans l'intégrateur** : la viscosité de volume est leur seule dissipation hors joints et frottement.

**Désactivé chez eux** : le DIF lui-même (`dpeftdif = R1`) et l'endommagement
diffus du volume (`df = R0`, le `(1−df)` ne s'appliquant **pas** au terme
visqueux) — ce dernier recoupant ce que leur article de 2026 dit du calcaire
et du grès.

---

## 3. Ce qui a été implémenté — six clés, toutes opt-in

Principe VIII respecté : **aucun défaut n'a changé**, chaque capacité est une
valeur nouvelle.

| clé | valeurs | défaut | source |
|---|---|---|---|
| `jointDeltaC` | `exact` \| `guo` \| **`solidity`** | `exact` | Y3Dfd.c l. 1099 |
| `jointFailRule` | `any` \| **`majority`** | `any` | Y3Dfd.c l. 1175 |
| `strainRateDIFArm` | `insertion` \| `envelope` \| **`continuous`** | `insertion` | Y3Dfd.c l. 1448-1456 |
| `gcBirth` | `ramp` \| **`penalty`** | `ramp` | Y3Did.c l. 915-964 |
| `gcBirthPenMin` / `gcBirthPenMax` | bornes | 0.01 / 3.0 | idem |
| `strainRateFilter` | `exponential` \| **`none`** | `exponential` | Y3Dfd.c l. 1448 |

Plus une valeur de config, sans code : **`potTangentFactor = 1.4286`**, pour que
`kt/kn = potTangentFactor/potPenaltyFactor` vaille les **2/7** de leur
`ktss = 2.0/(7.0)*penalty` (Y3Did.c l. 1017).

### Points de conception à retenir

- **`jointFailRule = majority` implique un endommagement PAR POINT.** Chez eux
  `z` est une variable locale de la boucle d'intégration ; rockim ne portait
  qu'un scalaire `J.D` par joint, déjà le *max* des points. Compter les points
  sur un max n'aurait rien voulu dire. D'où `Joint::Dk[]`.
- **`strainRateDIFArm = continuous` interdit d'appliquer le DIF en place** :
  le facteur se composerait à chaque pas. `snapBase()` garde les valeurs de
  base, `refreshDif()` reconstruit. Contrôle falsifiant : `difmed = 1,53036`,
  sous le plafond 1,85 de `yang-fig2`.
- **`gcBirth = penalty` est l'inverse de `ramp`.** `ramp` part de force NULLE
  et monte sur `gcBirthTau` ; `penalty` cale la pénalité de la paire pour que
  la force soit **continue**. Sous l'insert, où les joints meurent comprimés,
  la rampe était une perte de portance.
- **`strainRateFilter = none` et `strainRateDIFArm = continuous` se répondent.**
  La justification du filtre était « trop bruité pour **figer** un DIF dessus » —
  argument qui vise le gel. En continu, un pic ne dure qu'un pas.
- **Gardes d'exclusivité** : `gcBirth = penalty` exige `contact = potential` et
  refuse `gcBirthTau` ; `strainRateFilter = none` refuse `strainRateTau` ;
  `jointFailRule = majority` exige `jointQuadrature = midedge`. Aucune clé
  inerte silencieuse.

### Refactorisation

`setJointLengths()` centralise `dnE`/`dnF`/`slipF`, auparavant dupliqués en
**trois sites** qui avaient déjà divergé une fois.

---

## 4. Vérification

**Suite fast : 40 tests, 36 PASS.** Les 4 échecs (`fdem_voronoi_tension`,
`jointdeath_tension_2d`, `bulkmodel_neohooke_2d`, `jointquad_midedge_2d`) sont
la baseline MSVC déjà documentée dans `SUITE_full_MSVC_2026-08-25.md`.

**Preuve de neutralité.** Un binaire a été compilé depuis `HEAD` dans un
worktree jetable et passé sur la même suite : sur les **29 tests d'origine,
zéro différence de valeur physique**, et les mêmes 4 échecs aux mêmes
décimales.

**Nouveaux repères — 10 en tout.**

| tier fast | valeur |
|---|---|
| `jointdeltac_solidity_2d` | err −2,33345 %, 24 cassés |
| `jointfailrule_majority_2d` | err −2,78505 % |
| `dif_continuous_2d` | edot 7,65775, dif 1,53036 |
| `gcbirth_ramp_2d` / `gcbirth_penalty_2d` | err −1,67766 → −1,85267 % |
| `srfilter_none_2d` | edot 4,02148, dif 1,46994 |
| + 4 contrôles à charge nulle | 0 cassé, dampWork ≤ 0 |

| tier full | valeur |
|---|---|
| `jointdeltac_solidity_3d` | err −1,10893 %, **0 cassé** |
| `jointfailrule_majority_3d` | err −2,25444 % |
| `gcbirth_penalty_3d` | err −4,75889 %, 118 paires calées |
| `srfilter_none_3d` | edot 0,0746516, dif 1,21328 |
| `gcbirth_penalty_percussion_2d` | facteur 1,03077, 260 paires, résidu −0,994861 J/m |

### Trois pièges de repère, consignés

1. **`jointdeltac_solidity_3d` casse 0 joint, pas 200.** Vérifié : `guo` seul
   donne déjà 0 (−1,11074 %), `exact` en donne 200 (−1,27035 %). Les
   conventions à k_I ≈ 3 rendent le joint assez ductile pour survivre au
   déplacement imposé de cet essai. Ce n'est pas une régression.
2. **`gcbirth_penalty_2d` ne teste PAS le rééchelonnement.** En traction pure
   les joints meurent sans charge à relayer (`fDeath = 0`), le facteur retombe
   à 1 et le repère ne mesure que la *suppression de la rampe*. Il faut un
   indenteur : d'où `gcbirth_penalty_percussion_2d`, où 4 joints meurent
   **comprimés** en relayant 651 kN/m.
3. **`gcbirth_penalty_3d` ne verrouille pas un écart** mais le *fait* que le
   mécanisme s'arme (118 paires) : en traction 3D les deux modes donnent le
   même `err_pct`.

### Un garde-fou qu'il ne faut pas perdre

Le relevé de naissance par volume/aire ne servait pas qu'à adoucir : il
empêche une **injection d'énergie** sur les paires nées en recouvrement
(mesure historique **+936 J/m sans relevé**). `gcBirth = penalty` le supprime,
donc le **résidu du bilan devient le contrôle obligatoire** — le solveur
l'écrit lui-même : *« en mode penalty tout positif est une injection »*.
Mesure sur la percussion 2D : **−0,9174 J/m** en `ramp`, **−0,994861** en
`penalty`. Négatif donc dissipatif, et même légèrement plus. **À revérifier
sur tout nouveau cas.**

---

## 5. Le run de réplication `impact_imperial.cfg`

Lancé le 2026-08-26 à 17:01, 14 threads, sortie `out_imperial`.
**État au 2026-08-27 01:26 : t = 104,6 µs sur 550, soit 19,0 %.**

### Ce qui est validé — trois recoupements indépendants

| contrôle | prédiction hors calcul | observé |
|---|---|---|
| instant d'impact du piston | 0,2 mm / 10,66 m/s = **18,8 µs** | énergies démarrent à 19,2 µs |
| arrivée de l'onde à l'insert | 0,242 m / 5048 m/s → **66,7 µs** | l'insert décolle 65,9–67,8 µs |
| masse de l'outil | quantité de mouvement 25–60 µs | **1,2182 ± 0,0004 kg** = le bit seul, à 0,03 % |

Et **σ_zz max à la jauge = 205,7 MPa**, dans leur fourchette 200–260. Premier
critère sur sept : **validé**.

### Les sept critères à cet instant

| critère | mesure | cible | statut |
|---|---|---|---|
| contrainte max au bit | **205,7 MPa** | 200–260 | ✔ **dans** |
| vitesse d'indentation | pic 11,0 m/s (insert) | 9,40–9,85 | trop haut (transitoire) |
| vitesse de rebond | — | 6,87–7,10 | pas de rebond |
| rebond / indentation | — | 0,72–0,73 | idem |
| profondeur d'indentation | 0,294 mm | 1,45–1,60 | 19 % |
| fissure radiale max | 6,61 mm | 20,2–24,5 | 30 % |
| rayon de cratère | 4,02 mm | 10,0–12,1 | 37 % |

Fait notable : **la fissuration court en avance sur la pénétration** (cratère à
37 %, radiale à 30 %, enfoncement à 19 %).

### Le résultat qui compte

**239 joints rompus, et sur la dernière frame dépouillée : 85 en traction,
0 en cisaillement.**

C'est la même signature que le run de référence adaptatif (841 traction /
3 cisaillement, 0,36 %). **Les six conventions portées ne changent rien au mode
de rupture.** La chaîne causale du déficit est complète :

> pas de cisaillement → pas de surfaces broyées → pas de surfaces qui glissent
> → pas de frottement.

Chez eux ce poste porte **32 J sur 49,3** (65 %). Ici : **0,58 J** contre
11,8 J de fissuration.

⚠️ **Réserve.** On est à 0,29 mm sur 1,53. La pression sous l'insert va encore
monter, et c'est elle qui active le cisaillement dans un critère de
Mohr-Coulomb. **Conclure maintenant serait conclure sans avoir atteint le
régime où le broyage se produit.**

### Coût

Temps par frame : 12,5 → 26 → 43,5 → 56 → 84 → 110 → **103 min**. Plafonne
vers 105 min. Soit **~58 h** pour aller à 100 %.

Mais l'enfoncement atteint 1,53 mm vers **t ≈ 275 µs, soit 50 %** — donc
**~20 h** suffisent pour acquérir quatre des cinq critères manquants. Les 38 h
suivantes n'affinent que le cratère et les radiales.

**Cadence anormale, non expliquée** : 128 ms/pas pour 43 k tétraèdres sur
14 threads. Suspicion : en schéma intrinsèque les 81 797 joints existent tous,
et `potentialContact()` interroge `jointOfPair_` pour chaque paire candidate à
chaque pas uniquement pour découvrir qu'un joint vivant la porte. À investiguer
— l'optimisation serait sans effet sur la physique.

---

## 6. Corrections faites en dépouillant — toutes des chiffres plausibles et faux

Consignées parce que chacune aurait pu passer inaperçue.

1. **`imp_lib.broken()` utilisait `damage >= 0.999`.** Depuis
   `jointFailRule = majority`, `J.D` est le *max* sur les points : un seul
   point à 1 met `damage` à 1 sans tuer le joint. **206 facettes comptées
   contre 88 réelles** — facteur 2,3 qui gonflait radiale et cratère.
   **Corrigé** : critère `tBreak >= 0`.
2. **Déplacement du bit ≠ pénétration de l'insert.** Le bit se fait comprimer
   par le piston à son sommet pendant que l'insert n'a pas bougé : 0,102 mm lu
   là où il y en avait 0,004. Toutes les figures lisent maintenant `z_insert`.
3. **Enveloppe étoilée entièrement artefactuelle.** Deux causes cumulées : le
   centre pris sur 7 facettes rompues était faux de **0,67 mm** (Rayleigh
   p = 0,000, « anisotrope »), et le rayon **max** par secteur est corrélé
   **+0,79** au nombre de facettes du secteur. Recentré sur la zone de
   processus, au 90ᵉ percentile, découpage adaptatif : **p = 0,83, isotrope**.
4. **Force-pénétration, trois erreurs successives** : signe (vitesse en
   « descente positive » avec le `+g` d'un axe vers le haut → −65 kN) ; masse
   bit+insert appliquée à l'accélération du bit seul (l'insert à 9·10⁵ m/s²
   pèse 58 kN à lui seul) ; et l'idée fausse que jauge et Newton devaient
   coïncider — le transit de l'onde dure 48 µs pour un événement de 100, l'outil
   n'est **jamais** en équilibre quasi statique.
5. **Bannière `gcBirth` jamais imprimée** : placée avant la lecture de la clé.
   Troisième occurrence de ce piège d'ordre sur ce projet.

---

## 7. Écarts restants avec leur code

**Physique** : aucun écart connu et non porté.

**Discrétisation** :
- maillage 42 882 éléments contre leurs **230 788** ;
- leurs tétraèdres de volume sont des **T10 quadratiques avec F-bar**
  (Xiang *et al.* 2009) là où rockim est en T4. Leur élément **joint**, lui,
  est bien le TET4JOINT à 6 nœuds que rockim reproduit.

**Confirmé absent chez eux** : `jointResidualMu` — leur `dpefm` vaut `0.0` en
dur (l. 1091). Un joint rompu ne porte aucun cisaillement ; tout le frottement
vient du contact. La clé reste désactivée.

---

## 8. Prochain pas recommandé

1. **Laisser le run jusqu'à ~50 %** (~20 h), puis arrêter et dépouiller : cinq
   critères et un verdict solide sur le mode de rupture.
2. **Tester `jointFrictionScaled = 1`.** Le solveur avertit au démarrage que
   `jointShearUnload = origin` avec `jointFrictionScaled = 0` rend le
   glissement **réversible, sans hystérésis** — le joint ne dissipe rien en
   cisaillement, par construction. C'est la forme littérale de leur article,
   c'est une clé de config, et ça vise exactement le déficit.
3. **Instrumenter les 128 ms/pas** avant toute autre réplique fidèle.
