# LOT 2b — LE COUPLAGE ENDOMMAGEMENT → CONTACT : verdict sur source primaire
# et le DIF, avec une coquille corrigée par contrôle de continuité

*Fiche du 2026-08-29. Source lue de première main dans cette session, sur le PDF
fourni par F. Uzquiano.*

**Objet** : trancher la question qui a piégé la session du 29/08 et que la
[CORRECTION 2](../BILAN_interference_2026-08-29.md) avait laissée ouverte —
**où, exactement, l'endommagement se couple-t-il au contact chez Imperial ?**

---

## 0. La source, référence complète et vérifiée

> **Yang, X., Xiang, J., Naderi, S., Wang, Y., Aising, J., Ugarte, I., Latham, J.-P.**
> « **High-fidelity modelling of fragmentation and pulverisation in hard granite
> under percussion loading: a FDEM-based approach** », *International Journal of
> Rock Mechanics and Mining Sciences*, **206** (2026) 106660.

Affiliations imprimées p. 1 : Yang, Xiang, Naderi, Latham — Dept of Earth Science
and Engineering, Imperial College London ; Wang — Resource Geophysics Academy,
Imperial ; **Aising — Dept of Geosciences, Mines Paris - PSL University,
Fontainebleau** ; Ugarte — Drillco Tools S.A, Santiago.

Le titre annoncé le 2026-08-29 à partir d'une fiche ResearchGate est **confirmé
mot pour mot**. L'ordre des auteurs de l'en-tête de
[`yang2026_pulverisation.md`](yang2026_pulverisation.md) plaçait Wang en
troisième position ; l'article imprime **Yang, Xiang, Naderi, Wang, Aising,
Ugarte, Latham**. Corrigé ici.

---

## 1. LA RÉPONSE, en une phrase

**Le seul couplage (1−D) que l'article écrit est sur la CONTRAINTE DE L'ÉLÉMENT.
La pénalité de contact n'est jamais mentionnée : le mot « penalty » n'apparaît
pas une seule fois dans tout l'article.**

Le frottement, lui, est un **paramètre matériau calibré, constant** — pas une
fonction de D. Le lien entre endommagement et frottement existe dans le
*discours* des auteurs, jamais dans une équation.

---

## 2. CE QUI EST ÉCRIT — le modèle d'endommagement (§2.2, pp. 3-4)

**[LU]** Équations (3) et (4), p. 4 :

    sigma = { sigma_NH,                        0 <= eps <= eps_d
            { C_d eps = (1 - D) sigma_barre,   eps > eps_d

    D = { 0,                                                     delta_m <= delta_m^0
        { delta_m^f (delta_m^max - delta_m^0)
          / [delta_m^max (delta_m^f - delta_m^0)],     delta_m^0 < delta_m <= delta_m^f
        { D_max,                                                 delta_m > delta_m^f

Légende donnée p. 4, verbatim : « σ represents the stress, ε represents the
strain, ε_d represents the **damage initiation strain**, σ_NH represents the
**Neo-Hookean stress**, σ̄ represents the **effective stress**, and D represents
the **damage factor**. **C_d is a material-dependent constant.** »

**[LU]** Domaine d'application, p. 4, verbatim :

> « The damage model employed in this study is **assigned to finite elements
> (tetrahedral elements) and is not intended for joint elements**. »

et la division du travail, même page :

> « **joint-element failure represents meso-scale crack initiation and
> propagation**, while the **damage treatment assigned to tetrahedral elements
> represents local over-crushing and pulverisation-like stiffness degradation**. »

**Conclusion partielle [LU]** : D dégrade la **raideur/contrainte de l'élément
fini**. Rien d'autre n'est écrit sous forme d'équation.

---

## 3. CE QUI N'EST PAS ÉCRIT — et c'est le point décisif

### 3.1 La pénalité de contact : **ABSENTE**

**[ABSENT — recherche exhaustive]** Termes cherchés sur l'intégralité du texte
extrait (14 pages) : `penalty`, `contact stiffness`, `d_fact`, `1 - D`,
`contact penalty`.

| terme | occurrences |
|---|---|
| `penalty` | **0** |
| `contact stiffness` | **0** |

**L'article ne dit nulle part que la pénalité de contact, ou la raideur normale
de contact, est multipliée par (1−D) ou dégradée de quelque manière que ce soit.**

C'est la réfutation directe de l'hypothèse que la session du 29/08 avait tirée de
la *forme* du code local : `penalty = penalty*d_fact`. Cette forme n'a **aucun
répondant dans l'article**. Elle peut exister dans le vrai Solidity — je n'en
sais rien et l'article ne permet pas d'en décider — mais **elle n'est pas
documentée**, et rien ne permet de l'attribuer à Imperial.

### 3.2 Le frottement : un paramètre constant, pas une fonction de D

**[LU]** Table 1, p. 6 — le coefficient de frottement glissant est une **ligne de
propriété matériau**, au même titre que la densité :

| paramètre | granite Kuru Grey | carbure | acier |
|---|---|---|---|
| **Sliding friction coefficient** | **0,18** | **0,6** | **0,6** |

**Une seule valeur par matériau. Aucune ligne « frottement résiduel », aucune
ligne dépendant de D.** À noter aussi : **cette table ne contient ni pénalité ni
amortissement de masse**, contrairement à celle de l'ARMA 24-0952 pour St Anne.

**[LU]** §3.2, p. 6, verbatim :

> « The sliding friction coefficient is treated as an **effective post-fracture
> contact parameter** controlling the mobility and ejection of generated
> fragments, rather than as a direct representation of compacted powder friction.
> For example, **a lower sliding friction coefficient allows rock fragments to be
> ejected more easily after fracture**, reducing their support beneath the bit
> and consequently altering the bit penetration process »

C'est la description d'un **choix de calibration** — mettre le frottement du
granite bas — pas d'un mécanisme couplé.

### 3.3 La phrase ambiguë, et pourquoi je ne la surinterprète pas

**[LU]** Deux formulations, p. 4 puis p. 11, disent la même chose :

> p. 4 : « By applying larger deformations to elements at the impact centre and
> **reducing sliding friction**, the model reproduces experimentally observed
> bit-rock interactions and fragment ejection patterns. »
>
> p. 11 : « By assigning large deformations to the elements at the impact centre
> and simultaneously reducing **their** sliding friction coefficient, the model
> successfully captures key physical responses of hard granite »

Le « **their** » de la p. 11 renvoie grammaticalement aux « elements at the impact
centre » : lu au pied de la lettre, il suggère un frottement réduit **par
élément**, donc corrélé à l'endommagement.

**Mais trois éléments du même article tirent dans l'autre sens :**

1. Table 1 donne **une seule valeur globale** par matériau, sans ligne
   conditionnelle ;
2. §2.2 p. 4 présente le « **modified sliding friction treatment** » comme
   l'ingrédient *ajouté pour le granite*, par contraste avec St Anne et Rhune
   qui n'en avaient pas besoin — or St Anne est à 0,6 et le granite à 0,18 :
   **le « traitement » peut n'être que ce changement de valeur** ;
3. §6 p. 12 le décrit encore au singulier : « the sliding friction coefficient
   between joint elements after failure serves a role **analogous to the residual
   strength** of the failed material ».

**[VERDICT — AMBIGU, NON TRANCHÉ PAR LA SOURCE]** L'article **ne contient aucune
équation du traitement de frottement**. Il est impossible de décider, sur ce
texte, entre « μ réduit par élément en fonction de D » et « μ du granite
simplement calibré bas, globalement ». **Je ne tranche pas, et personne ne
devrait trancher sur cette base.**

C'est exactement la distinction que la CORRECTION 2 exigeait : ici, une
**intention** est décrite (« combined effect of damage evolution and frictional
contact behaviour », résumé p. 1) sans que le **mécanisme** soit publié.

---

## 4. CE QUE ÇA CHANGE POUR ROCKIM

Le WP7 de rockim avait été motivé par la lecture du code local, puis re-motivé
par l'article. Le bilan honnête, maintenant :

| ce que rockim fait | statut au regard de la source |
|---|---|
| dégradation de la contrainte d'élément par (1−D) | **CONFORME** — c'est l'équation (3), le seul couplage publié |
| couplage continu (1−D) sur le **frottement** | **NON RÉFUTÉ, NON CONFIRMÉ** — l'article décrit l'effet, jamais la formule. Le garder est légitime **à condition de le documenter comme un choix de rockim**, pas comme une réplication d'Imperial |
| couplage sur la **pénalité de contact** | **NON DOCUMENTÉ chez Imperial.** Ne pas l'implémenter en se réclamant d'eux. Si on le fait, c'est notre hypothèse, et elle porte notre nom |

**La « moitié manquante » diagnostiquée le 29/08 — la pénalité de contact non
couplée à D — n'est pas un retard de rockim sur Imperial. C'est une extrapolation
faite à partir d'un code qui n'est pas le leur.** Le canal d'effondrement de
portance existe bien chez eux, mais il passe par **(1−D) sur la contrainte
d'élément** (éq. 3) plus **un frottement bas** — pas, à notre connaissance
publiée, par la raideur de contact.

---

## 5. LE DIF — et une coquille de l'article de 2025, corrigée par le calcul

**[LU]** Équations (1) et (2), p. 3 de l'article 2026 :

    DIF_compression = { 1,                        eps_point <= 5,0e-6
                      { 0,77 + 0,56 eps_point^0,07,  5,0e-6 < eps_point <= 1,0e4
                      { 1,84,                     1,0e4 < eps_point

    DIF_tension     = { 1,                        eps_point <= 5,0e-6
                      { 0,95 + 0,41 eps_point^0,17,  5,0e-6 < eps_point <= 1,0e2
                      { 1,85,                     1,0e2 < eps_point

**[LU]** Où ils s'appliquent, verbatim p. 3 :

> « In the FDEM models, **DIF_Compression is applied to cohesion and G_II**, and
> **DIF_Tension is applied to tension strength and G_I**. According to the theory
> proposed by Zhao, the influence of the **internal friction coefficient** on
> strain rate **is considered to be insignificant**. »

Donc : le DIF touche **quatre** paramètres de joint (c, G_II, f_t, G_I) et
**laisse le frottement interne intact**. C'est net, et directement implémentable.

### La coquille

> **CRÉDIT — LE DÉPÔT L'AVAIT TROUVÉE AVANT MOI, ET MIEUX.**
> `include/rockim/YangDif.hpp` porte, en en-tête, l'analyse complète de cette
> coquille, faite le **2026-08-18**, donc **onze jours avant cette fiche et sans
> disposer de l'article de 2026**. Elle va plus loin que mon contrôle :
> * elle vérifie les **DEUX** raccords, pas seulement le haut : avec 0,07 la loi
>   vaut 1,1245 en 5e-6 /s (elle devrait valoir 1) et 1,5160 en 1e2 /s (elle
>   devrait valoir 1,85) — « elle ne se raccorde a aucune de ses bornes » ;
> * elle **relève l'exposant sur la figure 2(b) de l'article** et trouve
>   **0,1707**, qui raccorde EXACTEMENT les deux bornes (1,0010 et 1,8499) —
>   « deux raccords simultanes avec un seul parametre : ce n est pas une
>   coincidence » ;
> * elle confirme que l'éq. 2 (compression) est bien à 0,07 par sa figure 2(a) ;
> * elle **mesure une conséquence physique** que je n'avais pas vue : le saut de
>   22 % en 1e2 /s est un **attracteur** en insertion extrinsèque — un joint qui
>   franchit le seuil voit sa résistance bondir et cesse de s'insérer, si bien
>   que la population insérée s'empile juste sous 1e2. Mesure du 2026-08-18 :
>   médiane 99,36 /s (max 99,9988) avec l'exposant littéral, contre 40,22 /s
>   avec 0,1707 ;
> * elle **verrouille les deux variantes** par les contrôles `dif_yang_litteral_2d`,
>   `dif_yang_fig2_2d`, `dif_yang_fig2_plateau_2d` et `dif_yang_fig2_3d` de
>   `tools/verify_suite.py`.
>
> **Et l'article de 2026 valide cette inférence** : il imprime 0,17 là où le
> dépôt avait dérivé 0,1707 à partir de la seule figure. C'est une prédiction du
> dépôt confirmée par une source publiée un an plus tard. **Cela vaut d'être
> écrit dans le manuscrit.**
>
> Ce que cette fiche-ci ajoute, et rien de plus : la **confirmation sur source
> primaire** que l'exposant imprimé en 2026 est bien 0,17, ce qui clôt le doute.


**[LU]** L'article de 2025 (*IJRMMS* **191**, 106125, p. 3, éq. 3) imprime le
**même** DIF en traction avec l'exposant **0,07** au lieu de **0,17**.

**[INFÉRÉ — contrôle de continuité, calcul fait]** Les deux lois doivent raccorder
leur plateau à la borne haute de leur plage. Vérification :

| loi | valeur au raccord | plateau annoncé | écart |
|---|---|---|---|
| compression, exposant 0,07, à ε̇ = 1,0e4 | **1,8371** | 1,84 | 0,003 ✔ |
| traction, exposant **0,17**, à ε̇ = 1,0e2 | **1,8470** | 1,85 | 0,003 ✔ |
| traction, exposant 0,07, à ε̇ = 1,0e2 | **1,5160** | 1,85 | **0,334 ✘** |

**L'exposant de la loi en traction est 0,17. Le « 0,07 » de l'article de 2025 est
une coquille**, et elle est démontrable sans rien lire d'autre : avec 0,07 la loi
saute de 1,52 à 1,85 au raccord, soit une discontinuité de 22 %.

> **Conséquence pratique** : quiconque implémente le DIF d'après l'article de
> 2025 obtient une loi discontinue et sous-estime la résistance dynamique en
> traction de ~18 % sur toute la plage haute. **Utiliser l'article de 2026.**

---

## 6. LE RETRAIT DES FRAGMENTS, confirmé et daté

**[LU]** §2.1, p. 3 :

> « rock fragments are identified based on the **failure of joint elements**. As a
> result, some simulated fragments may be **completely surrounded by failed joint
> elements but still mechanically constrained** by adjacent intact rock blocks,
> meaning they cannot be considered free debris in a physical sense. In this
> study, after the impact, an **upward velocity and acceleration were assigned to
> all elements initially identified as rock fragments**. After running the
> simulation for a short period, the **displacement of each element was evaluated
> to determine whether it should be removed** »

L'algorithme détaillé est renvoyé à Yang et al. réf. 28 — c'est l'article
*IJRMMS* **191** (2025), qui en donne les trois étapes et les valeurs
(vitesse initiale 2,5 mm/s, accélération opposée à l'impact). Les deux sources
concordent.

**C'est un post-traitement, pas une érosion en cours de calcul.** Aucun élément
n'est supprimé pendant la simulation.

---

## 7. POURQUOI LE MODÈLE N'A PAS ÉTÉ APPLIQUÉ À ST ANNE — dit par eux

**[LU]** §2.2, p. 4, verbatim — c'est la justification la plus importante de
l'article pour nous :

> « For **St Anne limestone and Rhune sandstone** under hemispherical insert
> impact, previous validated FDEM studies were able to reproduce the main
> fragmentation characteristics **without introducing the additional damage model
> or modified sliding friction treatment**. This is consistent with experimental
> observations, where **severe pulverisation and intense fragment ejection were
> not prominent**. In contrast, both experimental evidence and previous modelling
> studies indicate that **Kuru Grey granite exhibits significantly more severe
> local pulverisation**, fragment ejection, and nonlinear rebound behaviour. »

**Imperial dit noir sur blanc que le modèle de pulvérisation n'est PAS nécessaire
pour St Anne, et que leurs simulations de St Anne ont réussi sans lui.**

C'est la confirmation, par la source, de ce que la fiche de 2026 §4 avait déduit,
et c'est plus fort : ce n'est pas seulement que la calibration granite ne se
transporte pas au calcaire — c'est que **le mécanisme entier est superflu sur le
calcaire**, de leur propre aveu.

> **Conséquence pour rockim** : sur un cas St Anne, activer le modèle de
> pulvérisation, c'est ajouter un mécanisme qu'Imperial juge inutile pour cette
> roche. S'il change les résultats, ce n'est pas une amélioration de fidélité,
> c'est un artefact. **Le banc de référence St Anne doit tourner SANS.**

---

## 8. LA MÉTHODE DE CALIBRATION, telle qu'ils la décrivent

**[LU]** §3.2, p. 5-6. Ce qui a été calibré, faute de données : **G_I, G_II, les
paramètres du modèle d'endommagement, et le coefficient de frottement glissant**.

Cas de calibration : trois vitesses de piston — **5,68, 9 et 13 m/s** — choisies
« because they cover low, intermediate, and high impact-energy conditions and
capture the **nonlinear variation of bit rebound velocity** ».

Évalué contre : vitesse d'indentation, vitesse de rebond, profondeur
d'indentation, masse de fragments, longueur de radiale, rayon de cratère,
morphologie de cratère, motifs de fissures observés au **scanner CT**, et
observations d'éjection à la **caméra rapide**.

Et l'avertissement, p. 6, verbatim :

> « G_I, G_II, the damage-model parameters, and the sliding friction coefficient
> were **not calibrated to fit a single output variable**. Instead, they were
> calibrated as an **integrated parameter set** to reproduce the overall bit-rock
> interaction process. »

Déjà relevé dans la fiche de 2026 §6 ; confirmé.

---

## 9. CE QUI RESTE OUVERT APRÈS CETTE FICHE

| question | statut |
|---|---|
| l'**équation** du traitement de frottement (μ constant ou μ(D) ?) | **NON PUBLIÉE.** Ni ici, ni dans les trois autres sources Imperial en main |
| l'**algorithme tangentiel** de frottement (ressort k_t ? régularisation ?) | **NON PUBLIÉ** — quatre sources concordent pour ne rien en dire : thèse Guo §2.3.4, manuscrit UCL p. 16, ARMA 24-0952 p. 3, cet article |
| la **règle pour une paire de matériaux différents** | **NON PUBLIÉE** |
| la **pénalité de joint** pour le granite Kuru Grey | **ABSENTE de la Table 1** — seule celle de St Anne et Rhune est publiée (ARMA 24-0952) |
| l'unité de C_d, et sa valeur | **NON PUBLIÉE** — « a material-dependent constant », sans plus |
| ε_d, la déformation d'amorçage de l'endommagement | **NON PUBLIÉE** en tant que telle — la Table 1 ne donne que δ_m^0, δ_m^f et D_max |

Le dernier candidat publié pour combler le frottement reste **Xiang, Latham &
Farsi (2017), « Algorithms and Capabilities of Solidity »**, Springer Proc. Phys.
188, ch. 16.
