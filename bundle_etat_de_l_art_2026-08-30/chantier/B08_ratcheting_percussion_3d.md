# B08 — Le « diffuse ratcheting » porté sur la PERCUSSION 3D
# et l'avertissement de lecture que la mesure a imposé

*Chantier du 2026-08-30. Verdict **M-8** du contre-audit §9. Diagnostic pur :
aucune force n'est touchée, aucun défaut n'est modifié.*

---

## 1. Ce qui manquait, et où

rockim possédait depuis le 2026-08-25 un diagnostic de la pathologie propre à la
**pénalité intrinsèque** : la part des joints portant déjà un `D > 0,01` **alors
qu'aucun n'a rompu**. Chaque joint cède un peu, la structure perd de la raideur,
et la résistance apparente baisse **sans qu'une seule fissure soit apparue**.

**Il était enfermé dans `if (scen_ == Scenario::BRAZILIAN)` (`src/FdemSolver.cpp:6770`)
et n'existait qu'en 2D** — `grep "gDfrac_" src/Fdem3dSolver.cpp` → **0**.
Autrement dit : la pathologie de l'insertion intrinsèque était diagnostiquée sur
un essai de **calibration 2D**, et **pas du tout sur la percussion**, c'est-à-dire
pas sur le cas comparé à Imperial.

> **C'est la jauge qui manquait à B1 et B5.** Sans elle, on mesure l'injection
> d'énergie de l'insertion intrinsèque sans pouvoir dire **quelle fraction des
> joints ratchette**.

## 2. Ce qui a été ajouté

Un relevé `scanSubCriticalDamage()` en 3D, figé :

* **en percussion** — au **pic d'effort d'outil**, l'instant où la structure est
  la plus chargée (analogue exact du pic brésilien) ;
* **en traction / cisaillement** (pas d'outil) — au dernier pas où `nBroken_ == 0`,
  soit juste avant la première fissure, **miroir exact du 2D**.

Les valeurs de **fin de run** sont imprimées en plus, ainsi que le nombre de
balayages, pour que le coût soit lisible. Le seuil `0,01` est tenu **en dur des
deux côtés** : un seuil configurable rendrait la comparaison dimensionnelle
muette — le mode de panne que la règle de parité existe pour éviter.

**Le relevé est inconditionnel et imprimé.** Il ne touche aucune force ; il lit
`J.D` et écrit trois scalaires. Aucune clé n'est introduite : une jauge opt-in
serait une jauge qu'on oublie d'armer.

## 3. La mesure — et elle tranche

`configs/fdem3d_percussion.cfg`, `T = 2·10⁻⁵ s`, 36 000 tets, 70 000 joints,
**`OMP_NUM_THREADS = 1`** (convention de la suite — voir §5 pour pourquoi ce
détail compte). **Seul le schéma d'insertion change**, la pénalité est la même :

| | joints > D = 0,01 **au pic** | D moyen | **pic d'effort** | joints rompus |
|---|---|---|---|---|
| insertion **intrinsèque** | **0,1571 %** | 8,099·10⁻⁴ | **8 729 N** | 0 / 70 000 |
| insertion **adaptative** | **0,0486 %** | 2,498·10⁻⁴ | **9 495 N** | 0 / 70 000 |
| rapport | **× 3,24** | × 3,24 | **− 8,1 %** | — |

*(À 4 fils : × 3,20 et − 8,8 %. La conclusion ne dépend pas du nombre de fils ;
seul le troisième chiffre bouge.)*

> **Trois fois plus de joints ratchetant, et un pic d'effort inférieur de 8,8 %,
> sans qu'AUCUN joint n'ait rompu de part et d'autre.** C'est exactement la
> signature annoncée : la structure est plus molle sous insertion intrinsèque,
> non parce que quelque chose s'est fissuré, mais parce que **chaque joint cède
> un peu**.

## 4. J'ai testé la causalité, et la mesure m'a REFUTÉ

Corrélation n'est pas mécanisme. J'ai donc posé le test qui devait trancher :
**si le ratcheting vient de la souplesse de la pénalité, raidir la pénalité doit
le réduire.** À schéma constant (intrinsèque), même deck :

| `jointPenaltyFactor` | dt | joints > D = 0,01 | D moyen | pic d'effort |
|---|---|---|---|---|
| **20** | 8,15·10⁻⁹ s | 0,1571 % | 8,099·10⁻⁴ | 8 729 N |
| **60** | 5,45·10⁻⁹ s | 0,1929 % | 1,187·10⁻³ | 8 557 N |
| **180** | 3,34·10⁻⁹ s | **0,2100 %** | **1,334·10⁻³** | 8 629 N |

**La jauge MONTE de 34 % (et le D moyen de 65 %) pour une pénalité multipliée par
NEUF.** L'hypothèse est réfutée, et proprement : la tendance est monotone sur
trois points, et elle l'est aussi à 4 fils (+36 %).

**L'explication tient en une ligne de code.** `D` est un déplacement
**normalisé** : `J.dnE = J.ft / J.pj` avec `J.pj = pf·E/h`
(`src/Fdem3dSolver.cpp:1813` et `:1588`), donc `dnE = ft·h/(pf·E)` **rétrécit
quand `pf` monte**. La **même ouverture physique se lit comme plus
d'endommagement.**

### 4.1 Ce que ce contre-test établit, et qui est plus fort que l'hypothèse

Le **pic d'effort reste PLAT à 1,2 % près** (8 729 / 8 557 / 8 629 N) sur un
facteur 9 de pénalité — alors qu'il bouge de **8,1 %** entre les deux schémas.

> **Les deux grandeurs bougent ENSEMBLE d'un schéma à l'autre, et SÉPARÉMENT
> d'une pénalité à l'autre.** C'est précisément ce qui distingue l'effet de
> schéma, réel, de l'artefact de normalisation. Sans ce contre-test, le ×3,2 du
> §3 aurait pu n'être qu'un effet de `dnE`.

## 5. ⚠️ L'avertissement de lecture, et il vise B1 directement

> **`gDfrac_` NE SE COMPARE PAS D'UNE PÉNALITÉ À L'AUTRE.**
>
> **Conséquence pour B1** — qui consiste précisément à poser
> `jointPenaltyFactor ≈ 9,6` au deck de réplique, contre 26,32 : **comparer la
> part de joints ratchetant avant et après ce changement n'a aucun sens.** Le
> **pic d'effort** et le **partage d'énergie**, eux, se comparent.

L'avertissement est écrit **dans `include/rockim/Fdem3dSolver.hpp`**, avec ses
trois points de mesure, pour qu'on ne puisse pas lire la jauge sans lui.

## 5 bis. Une sensibilité au NOMBRE DE FILS, trouvée en calibrant les repères

Mes premières mesures ont été faites à 4 fils ; la suite impose
`OMP_NUM_THREADS = 1`. Les repères ont donc **échoué**, et l'écart n'était pas de
l'arrondi :

| observable, insertion intrinsèque `pf = 20` | 1 fil | 4 fils | écart |
|---|---|---|---|
| **pic d'effort d'outil** | 8 728,85 N | 8 566,23 N | **1,9 %** |
| joints > D = 0,01 au pic | 0,157143 % | 0,155714 % | **0,9 %** |
| D moyen au pic | 8,0993·10⁻⁴ | 8,0809·10⁻⁴ | **0,23 %** |

C'est **reproductible**, donc ce n'est pas du hasard : c'est la non-associativité
des sommations OpenMP, amplifiée. **`peakF_` est un `max` sur un signal
oscillant** — la moindre divergence de sommation change *quel maximum local est
attrapé*, ce qui transforme un écart au dernier chiffre en écart de 2 %.

> **Et le résultat inattendu est là : la jauge de ratcheting est PLUS ROBUSTE que
> le pic d'effort.** `dfrac` est même **identique au chiffre près** entre 1 et
> 4 fils sur l'adaptatif (0,0485714 %) et sur `pf = 60` (0,192857 %). Des deux
> observables, c'est celle-ci qu'il faut citer.

**Ce que cela ne remet pas en cause** : les deux conclusions tiennent aux deux
nombres de fils — l'effet de schéma vaut × 3,24 (1 fil) et × 3,20 (4 fils), et
l'artefact de pénalité + 34 % (1 fil) et + 36 % (4 fils). **Le rapport 3,2 est
donc 4 à 5 fois le bruit de fil sur le pic, et bien davantage sur la jauge.**

**Conséquence pratique** : les tolérances des trois repères sont posées
**au-dessus** de ce bruit (± 200 N sur le pic, ± 3·10⁻³ point sur `dfrac`), sinon
ils échoueraient sur toute machine à autre nombre de cœurs — exactement le mode
de panne que les références « périmées » de la suite avaient déjà connu sous MSVC.

## 6. Coût

**118 balayages** en intrinsèque, **77** en adaptatif, sur 2 455 pas — le pic
d'effort est monotone, donc le nombre de relevés est borné par le nombre de pas
où l'effort progresse. Temps de paroi 41,0 s contre 25,3 s : l'écart vient du
schéma d'insertion, pas du relevé (l'adaptatif a moins de joints actifs).

## 7. Contrôles de non-régression

Références à **`OMP_NUM_THREADS = 1`**, tolérances posées au-dessus du bruit de
fil du §5 bis (± 200 N sur le pic, ± 3·10⁻³ point sur la jauge) :

| repère | ce qu'il verrouille |
|---|---|
| `ratchet_intrinseque_3d` (`full`) | 0,157143 %, D moyen 8,0993·10⁻⁴, pic 8 728,85 N, **0 joint rompu** |
| `ratchet_adaptatif_3d` (`full`) | 0,0485714 %, D moyen 2,4981·10⁻⁴, pic 9 495,44 N, **0 joint rompu** |
| `ratchet_penalite_artefact_3d` (`all`) | **verrouille l'AVERTISSEMENT, pas la capacité** : à `pf = 180` la jauge vaut 0,2100 % et le pic 8 629,0 N. Mon hypothèse a été réfutée par cette mesure ; ce repère la garde réfutée. |

Le couple des deux premiers est le test : un seul ne prouverait rien.

## 8. Réserves honnêtes

1. **Aucun joint ne rompt sur ce banc** (0/70 000 dans les cinq runs). Le
   diagnostic est donc mesuré dans le régime où il est le plus propre —
   ratcheting pur, sans fissuration pour le masquer. **Sur un run qui fissure,
   la lecture sera plus difficile**, et le rapport ×3,2 n'est pas transposable.
2. **Ce n'est pas le banc de réplique.** `configs/fdem3d_percussion.cfg` est un
   banc générique ; `bench_impact/configs/impact_imperial.cfg` exige
   `meshes/impact_fidele_r10.msh` et `_s15.msh`, **absents du dépôt**. Les
   chiffres ci-dessus ne valent donc pas pour St Anne.
3. **Le pic d'effort n'est pas isolé du schéma par autre chose que ce test.**
   L'insertion adaptative modifie aussi le pas de temps et le jeu de contacts
   actifs ; je montre que la pénalité seule ne déplace pas le pic, ce qui rend
   l'explication par le ratcheting la plus simple — **pas la seule possible**.
   Et il faut y ajouter le §5 bis : **le pic est bruité à 1,9 % par le seul
   nombre de fils**, donc les 8,1 % ne valent que ~4× ce bruit. La jauge de D,
   elle, est stable à 0,9 %, et c'est sur elle que le rapport × 3,24 repose.
4. **Le message est en français en 3D, en anglais en 2D** (`sub-critical damage
   at peak`). Asymétrie cosmétique, non corrigée pour ne pas casser un
   extracteur existant ; notée ici.
