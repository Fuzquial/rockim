# LOT 2c — L'ALGORITHME TANGENTIEL DE FROTTEMENT, enfin trouvé
# Xiang, Latham & Farsi (2017), « Algorithms and Capabilities of Solidity »

*Fiche du 2026-08-29. Source lue de première main sur le PDF fourni par
F. Uzquiano, extraction PyMuPDF, équations relues en préservant la mise en page
(l'extraction linéaire les brouille : elles sont composées en colonnes).*

---

## 0. Aveu préalable

J'avais annoncé à F. Uzquiano un **mauvais pronostic** sur cette pièce : « ce
chapitre porte sur l'empilement et l'interaction de formes complexes […] la
probabilité qu'il contienne l'algorithme tangentiel est modeste ». J'avais tort.
**Il le contient, explicitement, en deux équations.** C'était le dernier trou du
lot 2 et il est comblé.

Leçon de méthode : le titre d'un article prédit mal son §2. La section
« Governing equations » d'une communication de conférence rappelle souvent le
noyau du code entier, y compris ce que les articles d'application passent sous
silence.

---

## 1. La source

> **Xiang, J., Latham, J.-P., Farsi, A.** (2016/2017), « Algorithms and
> Capabilities of Solidity to simulate interactions and packing of complex
> shapes », *The 7th International Conference on Discrete Element Methods*
> (DEM7), Dalian, Chine, 1-4 août 2016, **paper number G010111** ; publié dans
> *Springer Proceedings in Physics* **188**, ch. 16.

Affiliation imprimée p. 1 : **Applied Modelling and Computation Group, Department
of Earth Science and Engineering, Imperial College London**. Le code y est nommé
**Solidity** dès le résumé : « a FEMDEM code, **Solidity**, is used to more
accurately capture the influence of complex shape ».

---

## 2. L'ALGORITHME, tel qu'il est publié (§2.2, p. 4)

### 2.1 Le contexte : la pénalité, c'est le NORMAL

**[LU]** p. 3, verbatim, première phrase du §2.2 :

> « In FEM/DEM, a **penalty function method** is employed to calculate the
> **normal** contact force when two particles are in contact. »

puis, éq. (3) p. 4, la force de contact distribuée de Munjiza (réf. [5]) :

    f_c = SUM_i SUM_j INTEGRALE_{beta_ci ∩ beta_tj} ( ... ) n d(Gamma)

> « Munjiza showed that **integration over finite elements was equivalent to
> integration over finite element boundaries** »

**[INFÉRÉ]** Ceci renforce la conclusion du [lot 2b](2026-08-29_lot2b_couplage_endommagement_contact.md) :
chez eux « penalty » désigne la **force normale de contact**, et le frottement
est un mécanisme **séparé**, traité par les équations ci-dessous. Il n'existe
donc, dans toute la littérature en notre possession, aucun endroit où la pénalité
serait couplée à l'endommagement.

### 2.2 LE FROTTEMENT — équations (4) et (5), p. 4

**[LU]** Attribution, verbatim :

> « **Xiang et al (2009) developed further the FEMDEM method by taking account of
> the sliding friction force.** The well-known classic Coulomb type friction was
> implemented and described as follows, »

**Équation (4)** — force tangentiale, régime adhérent :

    f_t = − k_t δ_t − η v_t

**Équation (5)** — bascule en glissement. Condition, verbatim : « **If f_t is
bigger than the friction force obeying the Coulomb-type friction law**,
f_t ≥ µ f_n, **the particles slide over each other and the tangential force is
calculated using the total normal contact force f_n** » :

    f_t = − µ f_n

**Légende, verbatim p. 4** : « where **η is the coefficient of viscous
dissipation**, **f_t is the tangential elastic contact force** and **v_t is the
tangential relative velocity** » ; « where **µ is the coefficient of sliding
friction** ».

> **Réserve d'extraction.** Les équations sont composées en colonnes ; une
> extraction linéaire les rend illisibles. J'ai relu les lignes (4) et (5) par
> extraction **positionnelle** (mots triés par abscisse à ordonnée constante),
> qui restitue littéralement `f = − k δ − ⟨η⟩ v` et `f = − µ f`. Les indices `t`
> et `n` sont composés en indice et se rattachent sans ambiguïté aux symboles
> nommés dans la légende. **δ_t n'est pas explicité dans la légende** : c'est le
> déplacement tangentiel relatif, seule lecture compatible avec `k_t δ_t` ayant
> la dimension d'une force. **[INFÉRÉ, dimensionnel]**

---

## 3. CE QUE ÇA ÉTABLIT

### 3.1 Il y a bien un ressort tangentiel

**Le frottement d'Imperial n'est PAS un Coulomb rigide.** C'est un contact
tangentiel **régularisé** :

* une **raideur tangentielle k_t** produisant une force élastique proportionnelle
  au déplacement tangentiel relatif ;
* un **amortisseur visqueux η** sur la vitesse tangentielle relative ;
* un **plafond de Coulomb** µ f_n qui prend le relais dès que la force élastique
  le dépasse.

C'est exactement la famille à laquelle appartient le `k_t = 2/7` de rockim.
**La question « rockim a-t-il le droit d'avoir un ressort tangentiel ? » est
réglée : oui, Imperial en a un.**

### 3.2 La chaîne d'attribution est bouclée

Trois sources se répondent enfin :

1. **Guo (2014), thèse, §2.3.4** — « A Coulomb friction law **was implemented
   into the three-dimensional FEMDEM code by Dr Jiansheng Xiang** », sans autre
   détail (déjà en fiche : [`guo2014_s234_235_contact_integration.md`](guo2014_s234_235_contact_integration.md) §1) ;
2. **ce chapitre** — « **Xiang et al (2009)** developed further the FEMDEM method
   by taking account of the sliding friction force », suivi des équations ;
3. donc **[INFÉRÉ, avec un fort degré de confiance]** la publication d'origine
   est **Xiang, J., Munjiza, A., Latham, J.-P. (2009)**, « Finite strain, finite
   rotation quadratic tetrahedral element for the combined finite-discrete
   element method », *IJNME* **79**(8), 946-978, DOI 10.1002/nme.2599.

**Le trou de publication que quatre sources laissaient béant n'en était pas un :
il était juste ailleurs que là où on le cherchait.** Il était dans l'article de
2009 sur l'élément, et rappelé dans une communication DEM de 2016 — jamais dans
les articles d'impact ni dans la thèse.

---

## 3ter. LA SOURCE PRIMAIRE EST TROUVÉE — le frottement d'Imperial est publié

*Ajout du 2026-08-29, dernier de la session. L'hypothèse formulée au §3bis
(lecture 2) est CONFIRMÉE.*

### La référence

> **Xiang, J., Munjiza, A., Latham, J.-P., Guises, R. (2009)**, « **On the
> validation of DEM and FEM/DEM models in 2D and 3D** », *Engineering
> Computations: International Journal for Computer-Aided Engineering and
> Software*, **26**(6), 673-687. Emerald. **DOI 10.1108/02644400910975469**.

Affiliations p. 673 : Xiang, Latham, Guises — Department of Earth Science and
Engineering, **Imperial College London** ; Munjiza — Department of Engineering,
**Queen Mary, University of London**.

### La phrase qui referme la chaîne d'attribution

**[LU]** p. 677, verbatim :

> « **In this paper, we develop further the FEM/DEM method by taking account of
> the sliding friction force.** The well-known classic Coulomb-type friction is
> implemented and described as follows, »

C'est **mot pour mot** la phrase que le chapitre DEM7 recopie en l'attribuant à
« Xiang et al (2009) ». **L'attribution visait bien cet article-ci**, et non
l'IJNME 79 éliminé au §3bis. La chaîne est close :

    Guo 2014, these §2.3.4  : « implemented [...] by Dr Jiansheng Xiang », sans detail
    Xiang, Latham & Farsi 2017 : recopie les equations, attribue a « Xiang et al (2009) »
    Xiang, Munjiza, Latham & Guises 2009, Eng. Comput. 26(6) : LA PUBLICATION D'ORIGINE

### Les équations, dans leur forme d'origine

**[LU]** p. 677, équations (8) et (9) :

    f_t = k_t delta_t - eta v_t                            (8)
    si |f_t| > mu |f_n| :   f_t = mu f_n                   (9)

Légende verbatim : « where **η is the coefficient of viscous dissipation**, f_t
is the **tangential elastic contact force** and v_t is the **tangential relative
velocity** » ; « where µ is the **coefficient of sliding friction** ».

> **⚠️ Divergence de signe entre les deux publications.** L'article de 2009 écrit
> `f_t = +k_t δ_t − η v_t` ; le chapitre de 2017, qui le recopie, écrit
> `f_t = −k_t δ_t − η v_t`. **Les mêmes auteurs, deux signes.** Seule la forme de
> 2009 est cohérente avec sa propre équation (9), qui pose `f_t = µ f_n` sans
> signe négatif. **[INFÉRÉ]** La convention de 2009 est la bonne, celle de 2017
> ajoute un signe pour orienter la force contre le glissement — ce qui est
> physiquement équivalent mais formellement différent. **Toute transcription doit
> choisir explicitement et le dire.**

### Et le côté DEM (particules sphériques), pour mémoire

**[LU]** p. 675, éq. (3b) — modèle ressort-amortisseur-patin de Cundall & Strack :

    f_t,ij = min{ k_t delta_t,ij - eta_t v_t,ij ,  mu |f_n,ij| t_ij }

Même loi, écrite en une seule ligne par un `min`. « η_n et η_t sont les
coefficients d'amortissement visqueux de contact **normal et tangentiel** ».

### CE QUI RESTE NON PUBLIÉ, définitivement

**k_t et η n'ont toujours AUCUNE valeur publiée.** Huit sources dépouillées, la
publication d'origine comprise : la **forme** est publiée, les **valeurs** ne le
sont pas. **La recherche s'arrête ici** — c'est un résultat, pas un échec :

> **Le `k_t = 2/7` de rockim ne peut être ni conforme ni divergent : il n'existe
> aucun nombre publié auquel le comparer. C'est un choix propre à rockim, à
> documenter et à justifier en son nom.** Position parfaitement tenable dans un
> manuscrit, et désormais adossée à une bibliographie exhaustive qui l'établit.

### DEUX CADEAUX que cet article fait à rockim

**1. Un banc de vérification analytique du frottement, prêt à l'emploi.**

**[LU]** §3, p. 677 — « Verification for FEM/DEM ». Un rectangle lancé sur un
plan horizontal avec une vitesse initiale, dont la distance d'arrêt vaut

    L = v_i^2 / (2 mu g)                                   (10)

Configuration publiée : côté **l = 0,05 m**, masse volumique **2650 kg/m³**,
coefficient de frottement **µ = 0,5**, module d'Young **E = 1,0×10⁹ Pa**, deux
pas de temps **Δt = 1,0×10⁻⁷ s** et **1,0×10⁻⁸ s**.

**C'est un contrôle de non-régression que rockim peut ajouter tel quel à
`tools/verify_suite.py`** : solution analytique fermée, configuration complète,
et il teste exactement le chemin tangentiel. **Recommandation : l'ajouter.**

**2. Un avertissement de stabilité qu'on ne devinerait pas.**

**[LU]** p. 677, verbatim :

> « with the larger of the two time steps, **the errors become significant**.
> However, using the larger time step, the calculation of FEM/DEM **with zero
> friction is fairly stable**. This **somewhat alarming conclusion** suggests
> that **in order to reduce the numerical error for calculation of tangential
> forces, the smaller time step is required**. »

**Le calcul des forces tangentielles exige un pas de temps plus petit que le cas
sans frottement.** Les auteurs le qualifient eux-mêmes d'« alarmant ». Aucune
autre source du corpus ne le dit, et c'est directement actionnable : le budget de
pas de temps de rockim (`FdemSolver.cpp:3105-3112`) compte la pénalité de joint
et la parabole — **compte-t-il le chemin tangentiel ?** À vérifier au lot 4.

---

## 3bis. L'ARTICLE DE 2009 EST ÉLIMINÉ — et l'attribution du chapitre est cassée

*Ajout du 2026-08-29, après lecture du PDF fourni par F. Uzquiano.*

**[LU]** Xiang, J., Munjiza, A., Latham, J.-P. (2009), « Finite strain, finite
rotation quadratic tetrahedral element for the combined finite-discrete element
method », *IJNME* **79**(8), 946-978, a été dépouillé intégralement (33 pages).

**Il ne contient PAS l'algorithme de frottement.** Recherche exhaustive sur le
texte extrait :

| terme | occurrences |
|---|---|
| `friction` | **0** |
| `tangential` | **0** |
| `sliding` | **0** |
| `Coulomb` | **0** |
| `stick` | **0** |

Les deux seules occurrences de `penalty` sont : une phrase d'introduction, p. 947
— « Once elements in contact are detected, a **penalty function method is
employed to calculate the normal contact force** when two particles are in
contact [16] » — et l'entrée [16] de la bibliographie, qui est Munjiza & Andrews
(2000). L'article porte **exclusivement sur l'élément tétraédrique quadratique**
et sur le *locking* du tétraèdre linéaire.

### L'attribution du chapitre DEM7 ne tient pas

Le chapitre écrit « **Xiang et al (2009)** developed further the FEMDEM method by
taking account of the sliding friction force », puis donne les équations (4)-(5).
Sa référence **[5]** est bien cet article. **Il ne décrit aucun frottement.**

**[INFÉRÉ]** Deux lectures possibles :

1. l'attribution est **approximative** — les auteurs renvoient à l'article de
   référence de leur élément 3D plutôt qu'à la publication du frottement, qui
   n'existerait pas ;
2. « Xiang et al (2009) » désigne en réalité la référence **[4]** du même
   chapitre : **Xiang, J., Munjiza, A., Latham, J.-P., Guises, R., « On the
   validation of DEM and FEM/DEM models in 2D and 3D », *Engineering
   Computations*, **26**(6), 673-687**. Le chapitre la date de 2008, mais le
   volume 26 n° 6 paraît en **2009** — d'où l'ambiguïté. Un article de
   **validation** est un endroit plausible pour décrire une implémentation de
   frottement.

**~~PROCHAINE PIÈCE À DEMANDER~~ : *Engineering Computations* 26(6), 673-687. — OBTENUE ET DÉPOUILLÉE, voir §3ter ci-dessus. L'hypothèse était juste.**
C'est la dernière hypothèse identifiée. Si elle ne donne rien non plus, la
conclusion est ferme : **k_t et η n'ont jamais été publiés**, et le chapitre DEM7
reste la seule source au monde pour la FORME de leur loi tangentielle.

### Un homonyme de plus, à ne pas confondre

**[LU]** L'article de 2009 publie en revanche une valeur de viscosité :
p. 969, cas de vérification du cube « jelly-like » — ρ = 1000 kg/m³,
λ = 10,8 kPa, µ = 7,22 kPa, et « a **viscosity constant responsible for damping**
of **1,0×10⁴ Pa·s** ».

**Ce η n'est PAS celui du frottement.** C'est le η **constitutif** du modèle
néo-hookéen (`T = … + ηD`, éq. 1 du manuscrit UCL), une viscosité de volume ; le
η de l'équation (4) du chapitre DEM7 est un **amortisseur tangentiel de
contact**. Deux objets, même lettre — comme les deux `D` du [lot 2b](2026-08-29_lot2b_couplage_endommagement_contact.md) §1.

Et la valeur 1,0×10⁴ Pa·s vaut pour un matériau de gelée à µ = 7,22 kPa :
**elle n'est transférable à aucune roche.** Elle ferme la question « existe-t-il
une valeur publiée de η constitutif ? » — oui, une, sans usage pour nous.

---

## 4. CE QUI RESTE NON PUBLIÉ

**[ABSENT]** Même dans ce chapitre :

> ### ⚠️ B9 — LA LISTE NOMINATIVE, arrêtée le 2026-08-30
> *Le contre-audit (§9, M-11) a relevé que **quatre nombres différents** ont
> circulé pour deux affirmations voisines : « quatre sources muettes » ici même
> (§6), « cinq » au lot 4 §1.1, « six » au lot 4 §2.6, « sept » ci-dessous et
> « huit » ailleurs. **Un décompte qui varie de 4 à 8 ne peut pas aller au
> manuscrit.** Voici la liste, nommée une fois pour toutes ; partout ailleurs on
> renvoie ici plutôt que de recompter.*
>
> **Les sources dépouillées pour la loi de frottement de contact — les huit :**
>
> | # | source | ce qu'elle donne sur le frottement |
> |---|---|---|
> | 1 | **Xiang, Munjiza, Latham & Guises**, *Eng. Comput.* **26**(6) (2009) 673-687 | **la loi elle-même**, éq. 8-9 p. 677 : `f_t = k_t δ_t − η v_t`, plafonnée à `µ f_n` |
> | 2 | **Xiang, Latham & Farsi** (2017), chapitre DEM7 « Algorithms and Capabilities of Solidity » | la loi reprise, éq. 3-5 p. 4 — **convention de signe opposée**, physiquement équivalente |
> | 3 | **Guo, Xiang, Latham & Izzuddin**, manuscrit UCL | contact et potentiel, éq. 14-18 — **rien sur k_t ni η** |
> | 4 | **Munjiza** (2004), *The Combined Finite-Discrete Element Method* | le potentiel, la z-curve — **rien sur k_t ni η** |
> | 5 | **Yang et al.**, *IJRMMS* **191** (2025) 106125 | DIF, fragments — **muette** |
> | 6 | **Yang et al.**, *IJRMMS* **206** (2026) 106660 | pulvérisation, éq. 3-4 — **muette** |
> | 7 | **ARMA 24-0952** | bilan d'énergie, éq. 3-7 — **muette** |
> | 8 | **ARMA 24-0788** | impact — **muette** |
>
> **Trois affirmations, trois nombres, et ils ne sont plus à recompter :**
>
> | affirmation | nombre | sur quelles sources |
> |---|---|---|
> | **la valeur de `k_t`** n'est publiée nulle part | **8 sources**, dont **2 donnent la loi sans sa valeur** | les huit ci-dessus |
> | **la valeur de `η`** n'est publiée nulle part | **8 sources** | idem |
> | **la règle de paire** (deux matériaux différents) n'est traitée nulle part | **8 sources** | idem |
>
> ⚠️ **Et la règle de paire n'est pas un silence coupable** : dans **tous** leurs
> essais publiés, les corps en contact sont du **même matériau** — la question ne
> se pose pas. À écrire ainsi, jamais comme un reproche.
>
> *L'IJNME **79**(8) (2009) 946-978 a été **ÉLIMINÉ** du corpus (§3bis) : il ne
> porte pas cette loi. Il n'entre donc dans aucun des trois décomptes.*

| grandeur | statut |
|---|---|
| **la valeur de k_t**, ou sa relation à la pénalité normale | **non donnée.** C'est le paramètre que rockim fixe à 2/7. ~~Sept sources muettes~~ → **8 sources, liste nominative ci-dessus (B9)**, l'IJNME 79 (2009) exclu — voir §3bis |
| **la valeur de η** (dissipation visqueuse tangentielle) | **non donnée.** Ne pas confondre avec le η CONSTITUTIF, dont une valeur est publiée (1,0×10⁴ Pa·s, sur une gelée) — §3bis |
| la **règle pour une paire de matériaux différents** | **toujours absente.** Ici tous les corps sont du même matériau : la question ne se pose pas, ce qui explique peut-être qu'elle ne soit jamais traitée |
| l'existence d'un **état de glissement mémorisé** (stick/slip avec historique de δ_t) | **non dit.** L'équation (4) suppose un δ_t, donc un état, mais sa remise à zéro au décollement n'est pas décrite |

**Recommandation, révisée le 2026-08-29** : l'IJNME 79 (2009) est **ÉLIMINÉ**
(§3bis). Le dernier candidat est **Xiang, Munjiza, Latham & Guises,
*Engineering Computations* **26**(6), 673-687**. S'il ne donne rien, on déclare
k_t et η **non publiés** et le `2/7` de rockim devient un choix documenté en
propre — position parfaitement tenable dans un manuscrit.

---

## 5. AUTRES ACQUIS DE CE CHAPITRE

### 5.1 Solidity a un mode CORPS RIGIDE

**[LU]** §2.1, p. 3, sous l'intertitre « **2. Rigid body solid** », les équations
de la dynamique d'un corps rigide :

    m_p,i v̇_p,i = f_d,i + m_p,i g + SUM_j ( f_cn,ij + f_ct,ij )
    I_p,i θ̈_p,i = SUM_j T_c,ij

> « where m_p,i and I_p,i are **mass and moment of inertia** of the particle i
> […] **T_c,ij is contact force torque**. For multiple interactions, the
> interparticle forces and torques are summed for k_i elements interacting with
> particle i. »

Deux enseignements :
1. **Solidity offre les deux formulations**, déformable (§2.1.1, éq. 1-4) et
   rigide (§2.1.2). Le `Tool3` rigide de rockim n'est donc pas une entorse :
   Imperial fait pareil quand la déformation du corps n'intéresse pas ;
2. la décomposition **f_cn + f_ct** (normal + tangentiel) est explicite dans
   l'équation du mouvement, ce qui confirme le §2.1 ci-dessus.

Un terme **f_d,i** apparaît aussi — une force d'amortissement — sans être défini.
**[OUVERT]**

### 5.2 Performances et parallélisation

**[LU]** §2.3, pp. 4-5 :

> « Recently Xiang **optimized the contact detection algorithm** in the FEMDEM
> and **parallelized the code using OpenMP** […] a **speedup of 6.5 on 8 threads
> and 9 on 12 threads** in a 3D deformable rock deposition with 288 particles.
> The **runtime was reduced by half after the code was optimized**. In future,
> the OpenMP implementation will be **redesigned with a hybrid MPI and OpenMP** »

**[INFÉRÉ]** À la date de ce chapitre (2016), Solidity est **parallélisé en
mémoire partagée seulement** (OpenMP), MPI restant un projet. Cela recoupe le
constat de la thèse de Guo (2014), qui décrivait un code **série** à 2,62 s/pas
pour 1,77 M d'éléments. À porter au lot 4 : c'est une donnée sur ce que le code
d'Imperial peut et ne peut pas faire, à comparer honnêtement à rockim.

### 5.3 Un jeu de paramètres à NE PAS transférer

**[LU]** §3.1, p. 5 : pour l'étude d'empilement, « the particles are modelled as
**rigid bodies** with density of **2.56 g/cm³** and **Coulomb coefficients of
friction from 0.2-1.0**, together with a **damping coefficient of 0.6** ».

**Ce sont des paramètres d'empilement de solides de Platon et d'Archimède, pas
d'impact.** Ne rien en tirer pour le forage percussif. Mentionné ici uniquement
pour que personne ne les recopie par erreur en croyant tenir une calibration.

---

## 6. ÉTAT DU LOT 2 APRÈS CETTE FICHE

Les six éléments que le brief demandait de reconstituer :

| élément | statut | source |
|---|---|---|
| loi de joint cohésive (modes I, II, mixte, adoucissement, plages) | **ÉTABLI** | Guo et al., manuscrit UCL, éq. 1-13, pp. 9-14 |
| algorithme de pulvérisation et couplage au contact | **ÉTABLI**, y compris ses limites | Yang et al. 2026, éq. 3-4 p. 4 ; [lot 2b](2026-08-29_lot2b_couplage_endommagement_contact.md) |
| DIF et son armement | **ÉTABLI**, coquille de 2025 corrigée | Yang et al. 2026, éq. 1-2 p. 3 |
| contact : potentiel, détection, pénalité, **frottement** | **ÉTABLI** | manuscrit UCL éq. 14-18 ; **ce chapitre, éq. 3-5 p. 4** |
| règle pour une paire de matériaux différents | **NON PUBLIÉE — et pour cause : tous leurs corps sont du même matériau** | **8 sources, liste nominative au §4 (B9)** |
| retrait des fragments et masse détachée | **ÉTABLI** | Yang et al. 2025 *IJRMMS* §2.3 ; Yang et al. 2026 §2.1 |
| calibrations publiées et transférabilité | **ÉTABLI** | [lot 2a](2026-08-29_lot2a_parametres_stanne.md) |

**Le lot 2 est complet**, à une exception près — la règle de combinaison du
frottement pour deux matériaux différents — et deux valeurs numériques
manquantes (k_t et η), toutes trois documentées comme non publiées.
