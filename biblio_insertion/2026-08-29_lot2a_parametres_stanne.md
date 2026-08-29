# LOT 2a — Les paramètres de St Anne chez Imperial, sur sources primaires
# et DEUX CORRECTIONS à des conclusions en vigueur au dépôt

*Fiche du 2026-08-29, écrite dès réception des PDF. Sources lues de première main
dans cette session (extraction texte PyMuPDF ; l'outil PDF natif du conteneur est
hors service, pdftoppm absent).*

**Urgence de cette fiche** : deux conclusions actuellement inscrites au dépôt sont
fausses pour le calcaire St Anne, et l'une d'elles inviterait à dégrader un deck
qui est juste. Elle est écrite avant le mémo de formulation complet pour cette
seule raison.

---

## 0. Les sources, et leur statut

| code | fichier fourni | référence | statut |
|---|---|---|---|
| **A** | `ARMA_paper_energy_XW_final.pdf` | **ARMA 24-0952**, Yang X., Xiang J., Naderi S., Wang Y., Latham J.-P. (Imperial) ; Aising J., Gerbaud L. (Mines Paris – PSL) ; Ugarte I. (Drillco Tools), « Where does the energy go in percussion drilling? FDEM's answer », 58e US Rock Mech./Geomech. Symposium, Golden, Colorado, 23-26 juin 2024 | **LU** intégralement |
| **B** | `main.pdf` | **Yang X., Xiang J., Naderi S., Wang Y., Aising J., Ugarte I., Latham J.-P.**, « Multi-criteria validation of hi-fidelity numerical model of impact breakage », *IJRMMS* **191** (2025) 106125 | **LU** en partie (tableaux, DIF, retrait des fragments) |
| **C** | `Manuscript_UCL_deposit.pdf` | **Guo L., Xiang J., Latham J.-P., Izzuddin B.**, « A generic computational model for three-dimensional fracture and fragmentation problems of quasi-brittle materials », manuscrit déposé (UCL Discovery) | **LU** §2.1 à §2.4 |

Étiquettes : **[LU]** = recopié du texte, citation à l'appui · **[INFÉRÉ]** =
déduction, raisonnement donné · **[OUVERT]** = question non tranchée par ces
sources.

---

## 1. LE TABLEAU DE ST ANNE, tel qu'Imperial le publie

**[LU]** Deux tableaux indépendants, deux articles, mêmes auteurs. Ils concordent
sur tout ce qu'ils ont en commun.

### 1.1 Source B, Table 4, p. 6 (page imprimée de la revue)

| paramètre | St Anne limestone | Rhune sandstone |
|---|---|---|
| Density (kg/m³) | 2731 | 2670 |
| Young's Modulus (GPa) | 57 | 37 |
| Poisson's Ratio | 0,31 | 0,2 |
| Energy Release Rate Mode I (J·m⁻²) | 12 | 40 |
| Energy Release Rate Mode II (J·m⁻²) | 800 | 1400 |
| Tensile Strength (MPa) | 7,0 | 7,5 |
| Cohesive Strength (MPa) | 18,8 | 33,6 |
| Internal friction coefficient | 1,0 | 1,0 |
| **Sliding friction coefficient** | **0,6** | **0,6** |

### 1.2 Source A, Table 1, p. 4 du PDF — le même, PLUS deux lignes

| paramètre | St Anne limestone | Rhune sandstone |
|---|---|---|
| Density (kg/m³) | 2731 | 2670 |
| Young's Modulus (GPa) | 57 | 37 |
| Poisson's Ratio | 0,31 | 0,2 |
| **Mass Damping Coefficient** | **4000** | **5000** |
| **Penalty Number (GPa)** | **3000** | **1800** |
| G_I (J/m²) | 12 | 40 |
| G_II (J/m²) | 800 | 1400 |
| Tensile Strength (MPa) | 7,0 | 7,5 |
| Cohesive (MPa) | 18,8 | 33,6 |
| Internal friction coefficient | 1,0 | 1,0 |
| **Sliding friction coefficient** | **0,6** | **0,6** |

**La source A est la seule en notre possession qui publie la PÉNALITÉ et
l'AMORTISSEMENT pour St Anne.** L'unité du coefficient d'amortissement de masse
n'est **pas imprimée** ; je ne la devine pas. **[OUVERT]**

### 1.3 Une coquille à ne pas recopier — source B, Table 5, p. 6

**[LU]** Le tableau acier/carbure de la source B est inconsistant avec celui de
l'article de 2026 (Table 1, déjà en fiche) :

| | source B, Table 5 « Steel » | source B, Table 5 « Carbide » | Yang 2026, Table 1 |
|---|---|---|---|
| ρ (kg/m³) | 15 250 | 7850 | carbure 15 250, acier 7850 |
| E (GPa) | 600 | 2000 | carbure 600, acier 200 |
| ν | 0,2 | 0,29 | carbure 0,2, acier 0,29 |

**[INFÉRÉ]** Les colonnes de la Table 5 sont **interverties** (l'acier n'a ni
ρ = 15 250 ni E = 600 ; le carbure, si), et le « 2000 » est vraisemblablement un
« 200 » fautif. Raisonnement : les trois lignes concordent avec la Table 1 de
2026 une fois les colonnes échangées, et le carbure de tungstène est
notoirement le plus dense et le plus raide des deux. **Quiconque transcrit cette
table telle quelle dans un deck obtient un outil à l'envers.**

---

## 2. CORRECTION 1 — le frottement 0,18 NE S'APPLIQUE PAS à St Anne

### Ce que le dépôt affirme aujourd'hui

[`yang2026_pulverisation.md`](yang2026_pulverisation.md) §3, dernier paragraphe :

> « ECART AVEC LE DECK DE ROCKIM : `contactMu = 0.6` y est pose GLOBALEMENT, et
> 0,18 n'est atteint que via `contactResidualMu` [...]. Le frottement roche/roche
> du deck est donc **3,3 fois trop eleve** sur l'immense majorite des contacts. »

### Ce que disent les sources primaires

**[LU]** Le coefficient de frottement glissant de **St Anne est 0,6** — source B
Table 4 p. 6, et source A Table 1 p. 4, deux fois indépendamment. Celui de
l'acier et du carbure est **0,6 également** (source B, Table 5).

Le **0,18 appartient au granite Kuru Grey**, dans l'article de 2026, et à lui
seul. Il arrive avec un modèle de pulvérisation que les articles de 2025 ne
contiennent pas (§5 ci-dessous).

### Verdict

**La conclusion du §3 de la fiche 2026 ne se transporte pas à St Anne.** Sur la
roche de la thèse, `contactMu = 0.6` **est exactement la valeur d'Imperial**. Le
deck de rockim n'a aucun écart à corriger sur ce point.

Ce qui reste vrai de ce §3 : que chez Imperial le frottement glissant est une
propriété **par matériau**, et que **la règle de combinaison pour une paire de
matériaux différents n'est publiée nulle part**. Ici la question ne se pose pas —
tous les matériaux sont à 0,6 — ce qui explique peut-être qu'ils n'aient jamais
eu à l'écrire. **[INFÉRÉ]**

> **À ne pas faire** : passer `contactMu` de 0,6 à 0,18 sur un cas St Anne au
> motif de l'article de 2026. Ce serait importer la calibration d'un granite sur
> un calcaire, ce que la fiche 2026 §4 interdit déjà pour d'autres raisons.

---

## 3. CORRECTION 2 — la pénalité de 3000 GPa était bien réelle

### L'historique de l'erreur

1. `DOCUMENTATION_rockim.md` §5.4 avait dérivé une pénalité de joint d'Imperial
   valant **≈ 52,6 E**, à partir d'un « 3000 GPa ».
2. La **CORRECTION 1 du 29/08**
   ([`BILAN_interference_2026-08-29.md`](../BILAN_interference_2026-08-29.md),
   §C1) a déclaré cette dérivation fausse : « **Raideur de joint** :
   `Spring_Stiffness = 900e9 Pa` (`mat.txt` l. 10), soit 15 E sur leur granite à
   60 GPa — et **NON les 52,6 E** que `DOCUMENTATION_rockim.md` §5.4 avait
   derives. Le « 3 000 GPa » de cette derivation est `D1PEPE`, la penalite
   d ELEMENT/CONTACT ».
3. Cette réfutation reposait sur la lecture de `/home/user/solidity`, dont la
   **CORRECTION 2 de la même nuit** a établi que ce n'est pas le code d'Imperial.

### Ce que dit la source primaire

**[LU]** Source A, Table 1, p. 4 : **Penalty Number = 3000 GPa** pour St Anne,
dont le module d'Young est **57 GPa**.

    3000 / 57 = 52,6

**La dérivation d'origine était juste. C'est la correction qui était fausse.**
Pour Rhune : 1800 / 37 = **48,6 E**. Les deux tournent autour de **50 E**.

### La tension que ça révèle — et elle est de premier ordre

**[LU]** Source C, p. 12, équation 9, énonce la recommandation des mêmes auteurs :

> « to maintain a balance between accuracy and computational efficiency, the
> value of the penalty term p₀ is usually chosen as » **E ≤ p₀ ≤ 10E**

avec, juste avant (éq. 8 et son commentaire, p. 12) :

> « ideally the value of the penalty term p₀ should be large enough so that the
> **extra elasticity introduced into the domain by the joint elements can be
> negligible** (Klein et al., 2001; Turon et al., 2007). However, a larger
> penalty term may cause numerical stability problems (Schellekens and de Borst,
> 1993), which usually requires smaller time-steps »

**Donc : la recommandation publiée du groupe est E ≤ p₀ ≤ 10E, et leur propre
pratique sur l'impact est ≈ 50 E — cinq fois au-dessus du haut de leur plage.**

**[INFÉRÉ]** L'explication la plus économique est que la plage E-10E de la
source C vise des cas quasi-statiques (brésilien, compression polyaxiale) où le
pas de temps n'est pas le facteur limitant du coût, tandis que l'impact exige une
raideur de joint bien plus haute pour ne pas amollir la roche sous le taillant.
Ce n'est **pas** écrit dans les sources : c'est une hypothèse, à confirmer.
**[OUVERT]**

### Ce que ça change pour rockim — et ça n'est pas ce qu'on croyait

Le `jointPenaltyFactor` de rockim vaut **20** par défaut. Il n'est donc pas
« au-dessus de la plage publiée » comme je l'avais écrit en lisant la seule
source C : il est **entre les deux** — deux fois au-dessus de la recommandation
générale, **2,6 fois EN DESSOUS de la pratique d'Imperial sur l'impact**.

Et la CORRECTION 1 du 29/08 avait mesuré que le levier pénalité pèse
**4,4 fois** le levier schéma d'insertion (§C2.3). Si la pénalité d'Imperial est
2,6 fois la nôtre sur le paramètre le plus sensible du problème, **c'est la
première chose à tester**, avant toute discussion sur l'insertion.

**Réserve honnête** : la Table 1 de la source A donne **un seul** « Penalty
Number » par roche, sans dire s'il désigne la pénalité de JOINT (le p₀ de
l'équation 6 de la source C) ou la pénalité de CONTACT, ou les deux. La source C
n'emploie « penalty » qu'au sens du joint. **[OUVERT]** — mais l'ordre de
grandeur, ≈ 50 E, tient quel que soit le référent.

---

## 4. LE BILAN D'ÉNERGIE SUR ST ANNE — la réponse chiffrée

**[LU]** Source A, §4, pp. 5-7. Piston à **9,41 m/s**, énergie d'entrée
**49,3 J**. Quatre stades :

| stade | fenêtre | ce qui s'y passe |
|---|---|---|
| **I** | 0 → 104 µs | le piston frappe le taillant. Énergie cinétique du piston tombée à ≈ **5 J** (11 % de l'initiale) ; celle du taillant montée à **40,9 J** |
| **II** | 104 → 267 µs | indentation et fissuration. **C'est là que tout se dissipe.** En fin de stade : énergie de rupture **1,66 J**, énergie de frottement **32,0 J**. Déformation élastique stockée : **13,9 J** dans le taillant, **4,23 J** dans la roche |
| **III** | 267 → 524 µs | rebond : ≈ **20 J** d'élastique redeviennent de la cinétique de taillant. Gain cinétique de la roche : minime |
| **IV** | 524 → 631 µs | le taillant, plus rapide que le piston au rebond, **le repercute** |

**Bilan final, St Anne** : frottement **32,0 J**, rupture **1,3 J**, soit
**2,6 % à la fissuration** et **64,9 % au frottement**.
**Rhune** : **2,46 %** et **62,2 %**.

**Tendance avec l'énergie d'impact** (Fig. 7, p. 6) : la part de rupture reste
**constante à ≈ 2,6 %** ; celle du frottement monte de **30 % (3 J) à 46 %
(22 J)** puis plafonne vers **65 %** pour le calcaire, **70 %** pour le grès
au-delà de 67 J.

Conclusion de l'article, p. 7 : « approximately **2.4% to 2.6%** of the energy is
utilized for crack propagation, while approximately **30% to 70%** of the energy
is used for friction between fragments ».

### Leurs postes d'énergie, et l'aveu qui va avec

**[LU]** Source A, p. 3, équations 3 à 7 : cinétique, déformation élastique
stockée, **contact** (« the elastic energy stored due to the **overlap** between
separated tetrahedral elements [...] a type of **recoverable** elastic energy »),
**rupture** (« includes elastic energy **and** plastic energy for joint
elements »), **frottement**, et potentielle gravitaire.

Puis, mot pour mot :

> « The **damping energy and numerical error** are obtained **by subtracting the
> above-mentioned energies from the total energy**. »

**Leur bilan n'est donc pas fermé indépendamment** : amortissement et erreur
numérique forment un résidu, non une mesure. **C'est un point où rockim, dont le
bilan d'énergie est fermé et contrôlé, est en avance sur la source.** À porter au
crédit du dépôt dans le lot 4.

### L'aveu le plus utile de tout le papier

**[LU]** Source A, p. 5, sur le stade II :

> « for the model used in this study, **rock fragments smaller than 1 mm cannot
> be further fractured**. However, in experiments, large fragments greater than
> 1 mm in size are indeed generated, and smaller fragments are also produced.
> [...] Therefore, the **friction energy output by FDEM simulation should include
> some of the fracture energy**, which creates a fracture area that is smaller
> than 1 mm. »

et, sur la cause mécanique :

> « Because of the relatively **sharp tetrahedral elements** in this simulation.
> When the bit impacts downwards, a **large friction force** will be generated,
> leading to an increase in friction energy consumption. »

**Traduction** : leur 65 % de frottement est gonflé par deux artefacts qu'ils
nomment eux-mêmes — le plancher de maillage à 1 mm, qui interdit de fissurer plus
fin et renvoie au frottement l'énergie qui aurait dû créer de la surface, et
l'angularité des tétraèdres, qui majore la force de frottement. **Le 2,6 % de
rupture est donc un plancher, pas une mesure.** Toute comparaison rockim/Imperial
sur la répartition d'énergie doit porter cette réserve.

---

## 5. DATATION DU MODÈLE DE PULVÉRISATION

**[LU]** Recherche exhaustive sur les trois articles de 2025 en notre possession
(sources A et B, plus `main_5.pdf` = Yang *JRMGE* 2025 et `main_3.pdf` = Naderi
*JRMGE* 2025) des termes : `pulveri`, `damage factor`, `load-bearing`,
`stiffness degradation`, `delta_m`, `D1PEM`.

**Aucune occurrence**, sauf une mention descriptive du mot « pulverised » dans la
description d'un cratère (`main_5.pdf`, p. 2), sans modèle derrière.

**Le modèle de pulvérisation est absent de tout le corpus 2025.** Il apparaît
dans l'article de 2026. **Le mécanisme a moins d'un an**, et il est né sur le
granite. **[INFÉRÉ, à partir d'un ABSENT bien étayé]**

---

## 6. CE QUE LA SOURCE A APPORTE ENCORE

**[LU]** p. 3, §3.1 — la géométrie exacte du banc, à reproduire :

* piston : cylindre **acier**, longueur **260 mm**, diamètre **26,5 mm** ;
* taillant : cylindre **acier**, longueur **265 mm**, diamètre **30 mm** ;
* insert **carbure hémisphérique, rayon 8,51 mm**, au bas du taillant ;
* plaque de chargement (« weight on bit »), conservée pour la distribution de
  masse **même si le WOB n'est pas appliqué dans cette étude** ;
* éprouvette : cylindre **150 mm de haut, 250 mm de diamètre** ;
* maillage : **1 mm** au centre, sur une demi-sphère de **25 mm de diamètre** ;
  **2 mm** jusqu'à **50 mm** de diamètre ; surface de l'insert à **0,7 mm** ;
* UCS : St Anne **179 MPa**, Rhune **165 MPa**.

**[LU]** p. 2 — **le code est nommé** : dans l'ARMA jumeau (24-0788, fichier
`ARMA_Impact_Final.pdf`, p. 1), les mêmes auteurs écrivent « our in-house rock
fracture software, **SOLDITY** [sic, pour Solidity] that leverages a hybrid FDEM
method ». C'est la première source en notre possession qui relie explicitement
les travaux d'impact de Yang/Latham/Xiang à **Solidity**.

**[LU]** p. 3 — la loi de joint y est rappelée sous une forme **différente** de
celle de la source C. L'équation (2) de la source A donne le cut-off
**explicitement en deux branches** :

    f_s = c − σ_n tan φ   si σ_n < f_t
    f_s = c − f_t tan φ   si σ_n > f_t

là où la source C (éq. 5, p. 11) écrit `f_s = c − σ_n tan φ` seule, en disant que
le cut-off est « **automatically guaranteed** » parce que σ_n ne peut pas
dépasser f_t. Les deux se rejoignent, mais **la source A est la forme à
implémenter** : elle est explicite et ne suppose rien du plafonnement de σ_n.
La figure 1 de la source A est créditée « (Liwei Guo, 2014) » : la thèse reste
la référence de formulation du groupe.

---

## 7. CONSÉQUENCES IMMÉDIATES POUR ROCKIM

Par ordre d'urgence, et sans anticiper le lot 4 :

1. **Ne pas toucher à `contactMu = 0.6`** sur les cas St Anne. C'est la valeur
   d'Imperial. (§2)
2. **Tester la pénalité de joint à ≈ 50 E**, contre les 20 actuels. C'est le
   paramètre le plus sensible du problème selon notre propre mesure du 29/08, et
   nous sommes 2,6 fois sous la pratique publiée. (§3)
3. **Vérifier l'amortissement de masse** : Imperial en met (4000 pour St Anne),
   unité non publiée. Savoir ce que rockim met, et sous quelle convention. (§1.2)
4. **Reprendre la comparaison d'énergie** avec les vrais chiffres : 2,6 % rupture,
   65 % frottement, et la réserve du plancher de maillage. (§4)
5. **Ne pas transcrire la Table 5 de la source B** sans intervertir les colonnes.
   (§1.3)

---

## 8. CE QUE CETTE FICHE NE TRANCHE PAS

* Le **couplage endommagement → contact** : hors du corpus 2025. Il faut
  l'article *IJRMMS* 206 (2026) en PDF. **[OUVERT]**
* L'**algorithme tangentiel de frottement** : la source C p. 16 le réduit aux
  trois mêmes phrases que la thèse de Guo (Coulomb, glissement si f_tan > µN).
  Toujours non publié. **[OUVERT]**
* La **règle pour une paire de matériaux différents** : sans objet ici (tout à
  0,6), non publiée ailleurs. **[OUVERT]**
* Le référent exact du **« Penalty Number »** (joint, contact, ou les deux).
  **[OUVERT]**
* L'**unité** du coefficient d'amortissement de masse. **[OUVERT]**
