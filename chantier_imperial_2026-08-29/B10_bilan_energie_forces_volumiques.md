# B10 — Fermer le bilan d'énergie : les deux forces VOLUMIQUES y entrent
# `energyBodyForces = off | on`

*Chantier du 2026-08-30. Suite du contre-audit §9, verdicts **M-7** (poste
gravitaire absent) et sa conséquence sur `brushWork_`. Opt-in, défaut
bit-identique, mesuré.*

---

## 1. Le constat, et il est plus grave que ce que j'en avais écrit

Le contre-audit relevait que **le septième poste du bilan d'énergie n'existe
pas** : `bodyForces()` applique la pesanteur **sans aucun compteur de travail**,
alors qu'ARMA 24-0952 pose explicitement l'**énergie potentielle gravitaire**
dans ses éq. 3-7, et que `gravity = 9.81` est posé dans **20 des 22 decks** de
`bench_impact/configs`, `impact_imperial.cfg` compris.

**J'avais écrit que le défaut était « structurel, pas numérique », l'ordre de
grandeur étant ~10⁻⁴ J contre ~49 J injectés. C'est vrai pour l'énergie
cinétique. C'est FAUX pour le résidu — et le résidu est justement ce que le
verdict et `budgetAbortPct` jugent.** Ma propre magnitude « invisible » était
la bonne réponse à la mauvaise question : la pesanteur est négligeable devant
l'énergie **injectée**, et dominante devant l'énergie **non expliquée**.

Mesure sur `configs/fdem3d_percussion.cfg` (36 000 tets, 70 000 joints,
144 000 nœuds, T = 2·10⁻⁵ s, `gravity = 9.81`) :

| | valeur |
|---|---|
| travail de la pesanteur, mesuré | **1,05738·10⁻⁶ J** |
| **résidu B4 du même run** | **1,05718·10⁻⁶ J** |

**Le résidu ÉTAIT le travail non compté de la pesanteur — à 0,02 % près.**
Autrement dit : sur ce run, tout ce que le bilan attribuait à l'erreur
d'intégration était en réalité un poste manquant. Une fois compté, le résidu
tombe à **−2,015·10⁻¹⁰ J**, soit un facteur **5 250**.

### 1.0 bis Et le garde-fou COUPE un run sain à cause de ça — mesuré

Ce n'est pas une inquiétude théorique. Même deck, `T = 9·10⁻⁶ s`
(1 105 pas, donc **un** passage du moniteur, qui tourne tous les 1 024 pas),
`budgetAbortPct = 2` :

| | `energyBodyForces = off` (défaut) | `= on` |
|---|---|---|
| **ENERGY ABORT** | **OUI, à t = 8,34411·10⁻⁶ s** | **non** |
| résidu au déclenchement | 2,36997·10⁻⁹ J = **175,1 % de l'échelle** | — |
| travail de la pesanteur au même instant | **2,37808·10⁻⁹ J** | — |
| le run va-t-il au bout ? | **NON, coupé à 93 %** | oui, t = 9·10⁻⁶ s |
| résidu final | (run coupé) | **−9,849·10⁻¹² J = 0,25 %** `[OK]` |

> **Le résidu qui déclenche l'arrêt est, à 0,3 % près, exactement le travail
> non compté de la pesanteur.** Le run est parfaitement sain — **zéro joint
> rompu**, aucune instabilité, le hotspot est à `|v| = 7,95·10⁻⁵ m/s. Le
> garde-fou censé attraper une divergence attrape **un poste manquant dans sa
> propre comptabilité.**

C'est la réponse chiffrée à la question laissée ouverte par le contre-audit
(« un garde-fou peut-il couper un run sain ? ») : **oui, et voici le cas.**
Cela lève aussi le seul obstacle à armer `budgetAbortPct` sur les decks de
réplique — geste que le contre-audit réclamait et que je ne pouvais pas
recommander tant que le bilan était incomplet.

### 1.1 Le second poste : le tri des fragments

`brushWork_` était tenu **hors de `sumW`** à dessein, et la raison écrite dans
`FdemSolver.hpp` est bonne :

> « après la journée du 2026-08-18 passée à traquer des pompes logées dans des
> canaux comptabilisés, en fabriquer une serait l'erreur à ne pas commettre. »

**Elle ne protège pourtant pas de ce qu'elle craint.** Le théorème
travail-énergie ne distingue pas une force « physique » d'une force
« numérique » : toute force appliquée aux nœuds fait un travail, et ce travail
est **soit dans `sumW`, soit dans le résidu**. Hors bilan, le travail du tri
tombe **entièrement** dans le résidu — où `budgetAbortPct` peut couper un run
**sain** sur un artefact purement numérique. **Le mesurer et le montrer protège ;
le cacher, non.**

---

## 2. Ce que fait la clé

`energyBodyForces = off | on`, **défaut `off`, dans les DEUX solveurs.**

**La MESURE est inconditionnelle et imprimée** — elle lit `v_` et écrit un
scalaire, elle ne touche **aucune force**, donc la physique est bit-identique
dans les deux cas. **Seule l'entrée dans `sumW` est opt-in.**

| | `off` (défaut) | `on` |
|---|---|---|
| `gravWork_` mesuré | oui | oui |
| `brushWork_` mesuré | oui (inchangé) | oui |
| ligne `forces vol.` au résumé | oui, avec **« ces J tombent dans le résidu ci-dessous »** | oui, avec **« DANS le bilan »** |
| dans `sumW` et dans l'échelle | **non** | **oui** |
| pèse sur le verdict `[OK\|CHECK]` et sur `budgetAbortPct` | non | **oui** |

Le défaut affiche en outre une bannière dès que `gravity > 0` est posé :
*« energyBodyForces = off (defaut) : gravity = … est pose, et le travail de la
pesanteur N ENTRE PAS dans sumW — il tombera dans le residu B4. Le resume
chiffre ce qu il y verse. »* — application directe de la règle maison, *une
capacité active et muette est indiscernable d'une capacité inerte*.

### 2.1 Convention du compteur

`gravWork_ += Σ F·v dt` avec `F = −ρ V₀ g ẑ / 4` par nœud (miroir 2D :
`/3`, selon `−ŷ`). **Même convention que tous les autres compteurs** — `v` lu au
moment de l'application (leapfrog), donc **positif quand la pesanteur injecte de
l'énergie cinétique**. La réduction OpenMP (`reduction(+ : gw)`) n'introduit
aucune écriture concurrente sur `f_`.

---

## 3. Preuve de bit-identité

Deux runs de percussion 3D identiques au seul `energyBodyForces` près,
sorties complètes comparées ligne à ligne :

```
diff <(grep -v "energyBodyForces|forces vol|residu|wall time" off.txt) \
     <(grep -v "energyBodyForces|forces vol|residu|wall time" on.txt)
  -> AUCUNE DIFFERENCE
```

Énergie cinétique finale **0,0145049 J** dans les deux cas, joints rompus 0/70 000
dans les deux cas, tous les postes du bilan identiques. **Les trois seules lignes
qui changent sont la bannière, la ligne de ventilation et le résidu** — c'est-à-dire
exactement la comptabilité, et rien de la physique.

`energyBodyForces` **absent** du deck et `= off` explicite donnent des sorties
**identiques au caractère près**.

---

## 4. Une clé MUETTE trouvée en chemin, et enregistrée

En cherchant où lire la clé, j'ai trouvé que **`gravity` n'est lu, en 3D, que
dans `placeTool()`** — dont la première ligne est
`if (scen_ == Scenario::TENSION) return;`.

**Conséquence : un deck de traction 3D qui pose `gravity` ne reçoit AUCUNE
pesanteur, et rien ne le lui disait.** Le solveur 2D, lui, lit `gravity` dans
`init()`, sans garde de scénario : **c'est une rupture de parité 2D/3D**, la
sixième du registre.

**Je n'ai pas changé le comportement** — le corriger changerait la physique d'un
deck existant, et ce n'est pas à moi d'en décider. **Mais il cesse d'être
silencieux** : un avertissement explicite est imprimé, vérifié :

```
[FDEM3D] AVERTISSEMENT : gravity = 9.81 est pose, mais le scenario est TENSION
et le 3D ne lit `gravity` que dans placeTool(), qui sort avant sur ce scenario.
LA PESANTEUR EST INERTE ICI. (Le solveur 2D, lui, l applique : rupture de parite
2D/3D, enregistree le 2026-08-30.)
```

Et la clé `energyBodyForces` elle-même est lue dans `init()`, **pas** dans
`placeTool()` : une clé validée dans une branche que le scénario n'atteint pas
est une clé silencieusement ignorée, et une faute de frappe y passerait inaperçue.

---

## 5. Parité 2D/3D — et le 2D dit autre chose, ce qui est instructif

**Le changement est porté des deux côtés**, à l'identique : même clé, même
défaut, même message d'erreur, même ventilation, même entrée dans `sumW` et dans
`gross`, aux deux endroits qui les calculent (`checkEnergyAbort()` et le résumé).

**Mesure 2D** (`configs/fdem_percussion.cfg`, `T = 5·10⁻⁵ s`, `gravity = 9.81`) :

| | valeur |
|---|---|
| travail de la pesanteur | **4,10104·10⁻³ J/m** |
| résidu au défaut | **−5,65299·10⁻² J/m** (0,0096 % de l'échelle, `[OK]`) |
| résidu, poste compté | **−6,06310·10⁻² J/m** (0,0103 %, `[OK]`) |

> **Et c'est une nuance qu'il faut garder.** En 2D la pesanteur ne fait que
> **7 %** du résidu, là où en 3D elle en faisait **99,98 %**. La raison est
> simple et elle disqualifie toute généralisation : **le banc 2D casse des
> joints**, donc il a un vrai résidu physique d'intégration, tandis que le banc
> 3D mesuré n'en casse **aucun** (0/70 000) — son résidu y était presque pur
> bruit, et la pesanteur pouvait donc le dominer.
>
> **Le facteur 5 250 n'est pas une propriété du correctif ; c'est une propriété
> du run.** Ce que le correctif garantit, dans les deux cas, c'est que le poste
> cesse d'être attribué à l'erreur d'intégration.

Je venais d'écrire, au lot 4 §3 entrée 8, que *« c'est le registre des ruptures
de parité qui est l'avantage »*. En créer une septième en corrigeant la
comptabilité du seul 3D aurait été le contre-exemple parfait.

**Preuve de non-régression globale** : suite complète au tier par défaut après
le changement — **`[suite] TOUT PASSE (44/44)`**. Aucune référence de la suite
existante n'a bougé, dans aucun des deux solveurs.

---

## 6. Réserves honnêtes

1. **Le facteur 5 250 vaut pour CE run — et le 2D le démontre** (§5) : sur un
   banc qui casse des joints, la part gravitaire tombe à **7 %** du résidu.
   Ce run 3D-ci ne casse **aucun** joint (0/70 000) : son résidu est presque pur
   bruit d'intégration, donc la pesanteur y domine d'autant plus facilement.
   **La mesure à refaire est celle du deck de réplique St Anne**, sur un run qui
   fissure — je ne peux pas la faire ici (~40 h).
2. **`brushWork_` n'est pas exercé** par cette mesure : le tri n'était pas armé
   (`0 J`). Son entrée dans `sumW` est donc **implémentée et lue, mais non
   mesurée sous charge**. À exercer sur un run qui arme `fragBrushStart`.
3. **Le poste gravitaire n'est pas l'énergie POTENTIELLE**, c'est le **travail**
   de la pesanteur. Les deux sont opposés au signe près sur une trajectoire, mais
   ARMA écrit un bilan en énergie potentielle : si l'on veut la comparer terme à
   terme à leur éq. 3-7, il faudra convertir, et dire lequel des deux on publie.
4. **Je n'ai pas armé `budgetAbortPct` sur les decks de réplique.** C'est le
   geste qui rendrait l'avantage revendiqué au lot 4 §3 entrée 2 réellement actif
   — mais il change ce qu'un run produit (il peut l'interrompre), et cette
   décision est celle de l'utilisateur. **Elle devient sûre maintenant** : avec
   `energyBodyForces = on`, le garde-fou ne jugera plus un run sur un poste
   manquant.

---

## 7. Contrôles de non-régression

Ajoutés à `tools/verify_suite.py`, tier `fast` :

| repère | ce qu'il verrouille |
|---|---|
| `ebody_defaut_3d` (tier `full`) | `= off` : `gravwork = 1,71228e-10 J` et résidu **1,71064e-10 J**, soit **116 % de l'échelle → verdict `[CHECK]`** sur un run où **rien n'est cassé**. **Le repère verrouille le DÉFAUT, y compris le fait qu'il soit faux.** |
| `ebody_on_3d` (tier `full`) | `= on`, même run : `gravwork` **identique**, résidu **−1,64459e-13 J**, facteur **1 040**. Le `gravwork` identique est la preuve que la clé ne change **que** la comptabilité, jamais une force. |
| `ebody_abort_defaut_3d` (tier `all`) | **le garde-fou DÉCLENCHE** au défaut (`abort = 2`), sur `gravwork = 2,37808e-09 J`. |
| `ebody_abort_on_3d` (tier `all`) | **le garde-fou NE DÉCLENCHE PAS** clé armée, et `gravwork = 2,70127e-09 J` — valeur **plus grande** parce que ce run-là **va jusqu'au bout**. L'écart entre les deux références *est* la preuve que l'autre a été interrompu. |

**Mécanisme ajouté à la suite pour ce dernier : `absent:<métrique>`.** Sans lui,
un contrôle non obligatoire dont la métrique manque est simplement « toléré » —
il **passe que le phénomène se produise ou non**, donc il ne prouve rien.
Vérifier qu'un garde-fou *ne* se déclenche *pas* demande d'exiger le vide, et
la suite ne savait pas le faire.

**Ce que ces repères NE prouvent pas** : que `brushWork_` entre correctement dans
`sumW` sous charge. Aucun repère de la suite n'arme le tri des fragments ; ce
trou est nommé ici et reste ouvert.
