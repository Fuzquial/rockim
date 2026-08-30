# A11 — La raideur tangentielle du contact entre dans le budget de pas de temps (3D)

*Chantier du 2026-08-29. Clé `dtBudgetTangential`. **Défaut `off` : bit-identique.***

---

## 1. La source

> **Xiang, J., Munjiza, A., Latham, J.-P., Guises, R. (2009)**, « On the
> validation of DEM and FEM/DEM models in 2D and 3D », *Engineering
> Computations*, **26**(6), 673-687, DOI 10.1108/02644400910975469, **p. 677**.

C'est la **publication d'origine de la loi de frottement d'Imperial** (ses
équations 8-9), identifiée au cours de la mission — voir
[`biblio_insertion/2026-08-29_lot2c_frottement_tangentiel.md`](../fiches/2026-08-29_lot2c_frottement_tangentiel.md) §3ter.
Elle porte, verbatim, un avertissement que ses propres auteurs qualifient
d'alarmant :

> « with the larger of the two time steps, **the errors become significant**.
> However, using the larger time step, the calculation of FEM/DEM **with zero
> friction is fairly stable**. This **somewhat alarming conclusion** suggests
> that **in order to reduce the numerical error for calculation of tangential
> forces, the smaller time step is required**. »

**Le calcul des forces tangentielles exige un pas de temps plus petit que le cas
sans frottement.** Aucune autre source du corpus Imperial ne le dit.

## 2. Ce que le code faisait avant

`src/Fdem3dSolver.cpp`, `computeStableDt()` :

    dtMin = min(dtMin, 2 sqrt( m_[i] / (K[i] + (toolSig_ ? 0 : nExtra * kp_)) ))

Seule `kp_` — la raideur de contact de l'outil — entrait dans le budget.
**`potKt_`, la raideur tangentielle du contact par potentiel, n'y entrait pas.**

Or **le solveur 2D la prend depuis longtemps** (`src/FdemSolver.cpp`,
`computeStableDt()` : `if (contactPot_) kContact = max(kContact, max(potP_, potKt_));`).
**C'était une rupture de parité 2D/3D, du côté qui porte l'impact.**

## 3. LE PIÈGE D'UNITÉS — et pourquoi le miroir du 2D est faux

**Ne jamais recopier le `max(potP_, potKt_)` du 2D dans le 3D.**

| | 2D (`FdemSolver.cpp:834-835`) | 3D (`Fdem3dSolver.cpp:696, 701`) |
|---|---|---|
| `potP_` | `potPenaltyFactor · maxE · thk_` → **N/m** | `potPenaltyFactor · maxE` → **Pa** |
| `potKt_` | `potTangentFactor · maxE · thk_` → **N/m** | `potTangentFactor · maxE · hmin_` → **N/m** |

En 2D les deux sont homogènes et le `max` a un sens. **En 3D, `potP_` est en
pascals** : le comparer à une raideur en N/m n'a aucun sens, et sur un maillage à
`hmin_ = 1 mm` le `max` aurait choisi `potP_`, **mille fois trop grand**, divisant
le pas de temps par **≈ 32 sans la moindre erreur visible**.

**Seule `potKt_` entre**, et elle a bien les mêmes unités que `kp_` et que `K[i]`.

Le journal du solveur imprimait déjà l'indice — « p = 5e+10 **Pa**, kt = 2.8572e+08
**N/m** » — mais personne ne l'avait lu comme un avertissement.

## 4. Ce que le code fait maintenant

`src/Fdem3dSolver.cpp`, `computeStableDt()` :

    double kContact = toolSig_ ? 0.0 : kp_;
    if (contactPot_) {
        if (dtTangential_) kContact = max(kContact, potKt_);
        else if (potKt_ > kContact) → AVERTISSEMENT au journal
    }

Clé `dtBudgetTangential = on | off`, **défaut `off`**, lue dans le bloc
`if (contactPot_)` de la configuration, avec exception sur toute autre valeur.

**Le risque n'est pas muet.** Sous `off`, si la configuration est justement une de
celles où la clé changerait quelque chose, le solveur **avertit**, chiffre le
facteur, cite la source et nomme la clé à poser.

## 5. Preuve de bit-identité au défaut — mesurée

Configuration `verify_fdem3d_tension.cfg`, `contact = potential`,
`potTangentFactor = 1.0` (le **défaut**) :

| | `dt` |
|---|---|
| `dtBudgetTangential = off` | **1,30191e-08 s** |
| `dtBudgetTangential = on` | **1,30191e-08 s** |

**Identiques.** C'est structurel, pas une coïncidence : `kp_ = maxE · hmin_`
(l. 565) et `potKt_ = potTangentFactor · maxE · hmin_` (l. 701) sont **le même
produit** quand le facteur vaut 1. Le `max` ne peut rien changer.

## 6. Mesure de l'effet — et elle corrige une affirmation du lot 4

Même configuration, `potTangentFactor = 1.4286` — **la valeur que posent les decks
d'impact** (`bench_impact/configs/impact_imperial.cfg:379`) :

| schéma d'insertion | `off` | `on` | écart |
|---|---|---|---|
| **intrinsèque** | 1,30191e-08 s | 1,27315e-08 s | **−2,21 %** |
| **adaptative** | 1,92829e-08 s | 1,83836e-08 s | **−4,66 %** |

> **JE CORRIGE MON PROPRE BILAN.** Le
> [lot 4](../fiches/2026-08-29_lot4_bilan_rockim.md) §2.6 écrivait :
> « dès que le frottement travaille sous l'insert, **la marge de stabilité n'est
> pas celle qu'on croit** ». **C'est exagéré.** L'effet mesuré est de **2 à 5 %**
> sur le pas de temps, pas un précipice.
>
> La raison est structurelle et vaut d'être retenue : le budget nodal est dominé
> par `K[i]`, la raideur des **ressorts de joint**, devant laquelle le terme de
> contact `nExtra · kContact` pèse peu. L'écart **double** en insertion
> adaptative — précisément parce que les ressorts de joint de la phase liée
> sortent du budget, comme le dépôt l'avait déjà établi le 2026-08-07 — mais il
> reste sous 5 %.
>
> **Le correctif reste justifié** : par la parité 2D/3D, par la source publiée,
> et parce qu'un budget de stabilité doit compter ce qui raidit le système.
> Il ne l'est **pas** par une urgence : aucun run n'explosait à cause de ça.

## 7. Contrôles de non-régression ajoutés

Dans `tools/verify_suite.py`, plus une métrique `dt` (regex sur le journal) :

| contrôle | ce qu'il verrouille |
|---|---|
| `dtbudget_tangential_defaut_3d` | le **défaut** vaut 1,30191e-08 s à `potTangentFactor = 1.4286` — donc la clé absente ne fait rien |
| `dtbudget_tangential_on_3d` | la **clé armée** vaut 1,27315e-08 s — donc elle fait quelque chose |

**Les deux forment un couple. L'un sans l'autre ne prouverait rien** : le premier
seul ne dirait pas que la capacité existe, le second seul ne dirait pas que le
défaut est sauf.

Les deux passent. Un contrôle 2D existant (`dif_yang_litteral_2d`) a été rejoué :
inchangé.

## 8. Ce que ce changement NE fait PAS

* Il **ne touche pas le 2D**, qui prenait déjà la raideur tangentielle.
* Il **ne change aucun défaut** : tous les decks existants tournent à l'identique
  tant qu'ils ne posent pas la clé.
* Il **ne dit pas quelle valeur de `potTangentFactor` est juste.** Le rapport
  `k_t/k_n = 2/7` des decks d'impact vient d'un code qui n'est pas celui
  d'Imperial, et **aucune valeur de k_t n'est publiée**, sur huit sources
  ([lot 2c](../fiches/2026-08-29_lot2c_frottement_tangentiel.md) §3ter).
  C'est l'objet de l'action A3, pas de celle-ci.
* Il **ne traite pas `potP_`**, dont la place éventuelle dans un budget de
  stabilité 3D reste une question ouverte — elle n'est pas une raideur nodale.

## 9. Réserve

Les deux valeurs de référence des contrôles sont mesurées **sur cette machine**
(g++, Linux, `OMP_NUM_THREADS = 1`). Le dépôt a déjà rencontré des écarts
Linux/MSVC de l'ordre de 0,1 % sur d'autres métriques
(`tools/verify_suite.py`, commentaire l. 100). La tolérance retenue est **1e-5**
en relatif ; si un jour ces deux contrôles échouent sous MSVC d'un écart de cet
ordre, c'est la tolérance qu'il faut relever, pas la mesure qu'il faut refaire.
