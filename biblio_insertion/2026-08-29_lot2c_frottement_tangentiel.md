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

## 4. CE QUI RESTE NON PUBLIÉ

**[ABSENT]** Même dans ce chapitre :

| grandeur | statut |
|---|---|
| **la valeur de k_t**, ou sa relation à la pénalité normale | **non donnée.** C'est le paramètre que rockim fixe à 2/7. Aucune source en main ne dit ce que fait Imperial |
| **la valeur de η** (dissipation visqueuse tangentielle) | **non donnée** |
| la **règle pour une paire de matériaux différents** | **toujours absente.** Ici tous les corps sont du même matériau : la question ne se pose pas, ce qui explique peut-être qu'elle ne soit jamais traitée |
| l'existence d'un **état de glissement mémorisé** (stick/slip avec historique de δ_t) | **non dit.** L'équation (4) suppose un δ_t, donc un état, mais sa remise à zéro au décollement n'est pas décrite |

**Recommandation** : le prochain et dernier candidat pour k_t et η est
**Xiang, Munjiza & Latham (2009), IJNME 79(8), 946-978**. C'est la source
d'origine désignée par ce chapitre. Si elle ne les donne pas, on déclare le
paramètre non publié et le `2/7` de rockim devient un choix documenté en propre.

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
| règle pour une paire de matériaux différents | **NON PUBLIÉE** | quatre sources muettes |
| retrait des fragments et masse détachée | **ÉTABLI** | Yang et al. 2025 *IJRMMS* §2.3 ; Yang et al. 2026 §2.1 |
| calibrations publiées et transférabilité | **ÉTABLI** | [lot 2a](2026-08-29_lot2a_parametres_stanne.md) |

**Le lot 2 est complet**, à une exception près — la règle de combinaison du
frottement pour deux matériaux différents — et deux valeurs numériques
manquantes (k_t et η), toutes trois documentées comme non publiées.
