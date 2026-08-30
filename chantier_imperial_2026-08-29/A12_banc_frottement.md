# A12 — Le banc analytique du frottement de rockim
# **SPÉCIFIÉ, NON IMPLÉMENTÉ** — et mon estimation d'effort était fausse

*Chantier du 2026-08-29. Cette fiche existe pour que le travail soit prêt, et
pour dire honnêtement pourquoi il n'a pas été fait ce soir.*

---

## 1. Ce qu'il faut construire, et pourquoi

**Le chemin tangentiel du contact de rockim n'a AUCUN contrôle.** Vérifié :
`tools/verify_suite.py` contient bien des contrôles de frottement —
`jointdeath_friction_2d`, `jointResidualMu`, `jointFrictionScaled` — mais **tous
portent sur le frottement de JOINT**. La branche
`Ft -= potKt_·dt·vt` puis écrêtage de Coulomb (`Fdem3dSolver.cpp:3310-3315`,
miroir 2D `FdemSolver.cpp:4733`) n'est exercée par aucun test.

Or c'est cette branche qui gouverne le poste que l'ARMA 24-0952 mesure à
**64,9 % de l'énergie d'impact** sur St Anne.

## 2. Le banc publié — configuration exacte

> **Xiang, J., Munjiza, A., Latham, J.-P., Guises, R. (2009)**, « On the
> validation of DEM and FEM/DEM models in 2D and 3D », *Engineering
> Computations* **26**(6), 673-687, **§3 p. 677**, « Verification for FEM/DEM ».

Un rectangle posé sur un plan horizontal, lancé à une vitesse initiale
horizontale, ralenti par le frottement avec la base, s'arrête à une distance

    L = v_i^2 / (2 mu g)                                            (leur eq. 10)

Configuration publiée, complète :

| grandeur | valeur |
|---|---|
| côté du carré | **l = 0,05 m** |
| masse volumique | **2650 kg/m³** |
| coefficient de frottement glissant | **µ = 0,5** |
| module d'Young | **E = 1,0×10⁹ Pa** |
| pas de temps testés | **1,0×10⁻⁷ s** et **1,0×10⁻⁸ s** |

Et **le résultat qualitatif est aussi précieux que le résultat métrique** :

> « with a small time step, the numerical results are in excellent agreement with
> the theoretical values. It is worth noting that with the **larger** of the two
> time steps, **the errors become significant**. »

**Un bon banc doit reproduire les deux** : l'accord au pas fin **et la dégradation
au pas grossier**. Un test qui ne vérifierait que le premier ne distinguerait pas
une implémentation correcte d'une implémentation qui ne dépend pas du pas.

## 3. POURQUOI CE N'EST PAS FAIT — et je corrige mon estimation

J'avais porté A12 à **effort FAIBLE** dans le [lot 4](../biblio_insertion/2026-08-29_lot4_bilan_rockim.md)
§6 et dans la [synthèse](../SYNTHESE_etat_de_l_art_2026-08-29.md) §6 étape 10,
en supposant qu'il suffisait d'écrire un deck. **C'est faux.**

**rockim n'a pas de scénario capable d'exprimer ce banc.** Les scénarios du
solveur 2D sont `percussion | shear | tension | brazilian | shpb`
(`src/FdemSolver.cpp:78-85`). Aucun ne décrit un corps **libre**, posé sur un
plan, **lancé à une vitesse initiale**. Il n'existe par ailleurs aucune clé de
vitesse initiale de corps : `v_` est initialisé à zéro
(`src/FdemSolver.cpp:2039`) et n'est ensuite écrit que par les pilotes de bord.

Le coût réel d'un scénario `slide` :

    grep -c "scen_ == Scenario::" src/FdemSolver.cpp   ->  43 points de branchement
    repartis sur 17 fonctions : init, setupBoundaries, setupConfinement,
    placeTool, toolContact, computeStableDt, integrate, step, writeFrame,
    historyHeader, historyRow, finalize, finished, buildFromTriangles, ...

**Estimation honnête : 1 à 2 jours**, et surtout **un risque de régression réel**
sur une suite de 97 contrôles — pour héberger un seul banc. C'est disproportionné,
et c'est contraire à la discipline du dépôt, qui veut des changements minimaux,
opt-in et bit-identiques.

## 4. Les trois voies, avec leur coût et leur risque

| voie | ce que c'est | coût | risque | ce qu'on obtient |
|---|---|---|---|---|
| **V1 — scénario `slide`** | un scénario neuf dans le solveur de production | 1-2 j | **élevé** : 43 points de branchement, suite à 97 contrôles | le banc publié **à l'identique**, citable dans le manuscrit |
| **V2 — pilote autonome** | un `tools/slide_point.cpp` sur le modèle de `tools/yan_point.cpp` | 0,5-1 j **+ refactor** | moyen | un test du chemin tangentiel, **mais pas le banc publié** |
| **V3 — ne rien faire** | consigner le manque | 0 | — | le trou reste, et il est nommé |

**La voie V2 a un obstacle que la fiche doit dire** : `tools/yan_point.cpp` marche
parce que la loi de joint vit dans un en-tête réutilisable
(`include/rockim/YanSoftening.hpp`), si bien que le pilote « exercise the very
same functions the solver calls » — c'est écrit dans l'en-tête même. **La loi
tangentielle, elle, est enfouie dans la boucle de maillage de
`potentialContact()`** : un pilote ne pourrait pas l'appeler sans qu'on l'extraie
d'abord dans un en-tête. Le refactor est faisable et propre — c'est exactement ce
que le dépôt a fait pour `YangDif.hpp` et `YanSoftening.hpp` — mais il s'ajoute au
coût.

## 5. Recommandation

**V2, en deux temps**, et dans cet ordre :

1. **extraire** la loi tangentielle dans `include/rockim/TangentialContact.hpp`,
   à raideur et coefficient explicites, sans état de solveur — refactor **à
   comportement bit-identique**, verrouillé par la suite existante ;
2. **écrire** `tools/slide_point.cpp` qui intègre un bloc rigide de masse
   ρ·l²·(épaisseur) sous gravité avec **cette fonction-là**, et compare la
   distance d'arrêt à `v_i²/(2µg)` **aux deux pas de temps publiés**.

On perd la géométrie du banc — ce ne sera pas un rectangle maillé — mais on gagne
un test **falsifiable, exact, et qui exerce le code expédié**. Et le manuscrit
peut citer la source pour la loi et pour la solution analytique, en disant
honnêtement que la vérification porte sur le point matériel.

**V1 reste la seule voie qui reproduise le banc publié tel quel.** Elle se
justifierait si le scénario `slide` servait à autre chose ensuite.

## 6. Ce qui est acquis dès maintenant

Même non implémenté, ce travail a produit trois choses :

1. **La constatation, vérifiée, qu'aucun contrôle n'exerce le chemin tangentiel
   de contact** — c'est un trou nommé, avec sa preuve.
2. **La configuration publiée du banc**, complète, prête à être posée.
3. **La correction d'une estimation fausse de mon propre bilan** : A12 n'est pas
   un effort faible. Les deux documents qui la portaient sont corrigés.
