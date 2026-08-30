# LOT 1 — Bibliographie annotée : le FDEM d'Imperial College
# (Y → VGW/VGeST → Solidity), impact et forage percussif

*Fiche du 2026-08-29, en réponse au [brief de mission](../MISSION_etat_de_l_art_2026-08-29.md)
§2. Rédacteur : session Claude « état de l'art ». Commanditaire : F. Uzquiano.*

---

## 0. AVERTISSEMENT DE MÉTHODE — lire avant d'utiliser une seule ligne de cette fiche

### 0.1 Ce qui a été matériellement possible dans cette session

Le conteneur de cette session **n'a pas d'accès sortant au web au-delà du moteur
de recherche**. Toute tentative de récupérer une page ou un PDF — y compris des
pages entièrement publiques — est refusée par la politique de sortie réseau :

    imperial.ac.uk        → EGRESS_BLOCKED
    spiral.imperial.ac.uk → EGRESS_BLOCKED
    solidityproject.com   → EGRESS_BLOCKED
    sciencedirect.com     → EGRESS_BLOCKED
    arxiv.org             → EGRESS_BLOCKED
    orchyd.eu             → EGRESS_BLOCKED
    doi.org, crossref, openalex, semanticscholar → EGRESS_BLOCKED

(vérifié par sondage direct ; le journal du mandataire enregistre
`gateway answered 403 to CONNECT` pour chacun.)

**Conséquence dure : AUCUN article n'a été lu en plein texte dans cette
session.** Ce qui suit est construit sur les métadonnées et les fragments de
résumé que le moteur de recherche restitue. C'est suffisant pour établir
*quelles pièces existent, qui les a écrites, où elles sont et ce qu'elles
coûtent à obtenir* — c'est exactement le livrable 1 demandé. Ce n'est **pas**
suffisant pour établir une équation, une valeur numérique ou un algorithme.

### 0.2 Étiquettes épistémiques employées ici

| étiquette | signification |
|---|---|
| **[MÉTA]** | référence bibliographique (auteurs, titre, revue, volume, année, pages, DOI) relayée par le moteur. Fiable sur l'identité de la pièce, à revérifier au téléchargement. |
| **[RÉSUMÉ]** | contenu tiré d'un résumé ou d'un fragment de résumé relayé par le moteur. Dit de quoi parle la pièce ; ne dit pas ce qu'elle démontre ni comment. |
| **[INFÉRÉ]** | déduction explicite du rédacteur à partir de [MÉTA]/[RÉSUMÉ] ou de faits déjà acquis au dépôt. Le raisonnement est donné. |
| **[SUPPOSÉ]** | hypothèse de travail, non étayée. À vérifier avant tout usage. |
| **[LU-DÉPÔT]** | lu de première main, mais dans une pièce **déjà présente au dépôt** (fiches antérieures, PDF fournis par F. Uzquiano lors de sessions précédentes). |

**Il n'y a pas une seule occurrence de « [LU] » dans cette fiche pour une source
Imperial.** C'est voulu. Le lot 2 ne pourra être écrit qu'après réception des
pièces listées au §5.

### 0.3 Rappel du piège de la session précédente

`/home/user/solidity` **n'existe même pas dans ce conteneur** (vérifié :
`ls: cannot access '/home/user/solidity'`). Il n'a donc pas pu servir, ni en
bien ni en mal. Le rappel reste utile pour la suite : ce dossier est une
transposition d'équipe, jamais une preuve de ce que fait Imperial
(cf. [`BILAN_interference_2026-08-29.md`](../BILAN_interference_2026-08-29.md)
CORRECTION 2).

---

## 1. LE PAYSAGE — trois lignées FDEM à ne jamais confondre

Le FDEM n'est pas un code, c'est une famille. Trois branches descendent du même
tronc (Munjiza, fin des années 1980-1990) et **publient des formulations
voisines mais non identiques**. Citer l'une pour l'autre est précisément
l'erreur que ce dépôt s'interdit.

| lignée | code(s) | acteurs | où |
|---|---|---|---|
| **Imperial College + QMUL** | Y2D/Y3D → **VGW/VGeST** → **Solidity** | Munjiza, Latham, Xiang, Guo, Obeysekara, Lei, Farsi, Naderi, Yang | Londres |
| **Los Alamos (LANL)** | **HOSS** | Munjiza (co-auteur), Rougier, Knight | États-Unis |
| **Toronto** | **Y-Geo**, puis **Irazu** (Geomechanica) | Grasselli, Mahabadi, Lisjak, Tatone | Canada |

**[INFÉRÉ]** Munjiza est co-auteur des trois branches à des époques
différentes : sa signature sur un article **n'établit pas** que l'article décrit
Solidity. Exemple concret rencontré dans cette revue :

* Rougier, Knight, Broome, Sussman & Munjiza, *IJRMMS* **70** (2014) 101-108,
  « Validation of a three-dimensional finite-discrete element method using
  experimental results of the split Hopkinson pressure bar test » **[MÉTA]**.
  Le titre attire (validation dynamique 3D, granite, SHPB — exactement notre
  sujet). Mais les auteurs sont l'équipe **LANL/HOSS**. **Cette pièce ne dit
  rien de Solidity.** À lire pour le fond de la méthode, jamais à citer comme
  « ce que fait Imperial ».

Même réserve pour toute la littérature Y-Geo/Irazu (Lisjak, Mahabadi, Tatone),
souvent en accès plus facile et plus explicite sur les équations : elle décrit
**la loi de Munjiza telle que Toronto l'a implémentée**. Utilisable pour
comprendre la forme générale ; jamais comme preuve du contenu de Solidity.

---

## 2. GÉNÉALOGIE DU CODE D'IMPERIAL

**[MÉTA + RÉSUMÉ]**, recoupé sur plusieurs pages :

1. **Y** — le code de référence de Munjiza, publié **avec le livre de 2004**
   (listing source inclus). Base 2D.
2. **Y2D / Y3D** — bibliothèques FDEM libres issues de Y ; Y3D est attribué à
   Munjiza (2004) et **Xiang et al. (2009)**.
3. **VGW — Virtual Geoscience Workbench** — projet EPSRC de 5 ans, ~0,8 M£,
   **2004-2009**, co-dirigé par **J.-P. Latham** (Imperial) et **A. Munjiza**
   (Queen Mary). Objectif déclaré : rendre les codes Y, 2D et 3D, disponibles en
   **source ouverte** pour les géosciences. Atelier de lancement le
   **30 mars 2009**, RSM Building, Imperial. Distribution : **SourceForge**
   (`sourceforge.net/projects/vgw/`).
4. **VGeST — Virtual Geoscience Simulation Tools** — la plateforme FDEM séparée
   issue de VGW.
5. **Solidity** — **relancé sous ce nom en 2016** avec la mécanique du contact
   et la modélisation de la rupture. Site : `solidityproject.com`. Décrit comme
   « Open Source general purpose, two and three dimension finite element –
   discrete element solid mechanics code », avec « fracture and fragmentation
   **without the need to seed** » **[RÉSUMÉ, verbatim du site relayé]**.

> **[INFÉRÉ — important pour le lot 3]** « fracture and fragmentation *without
> the need to seed* » est une formule marketing qui, prise au mot, décrit
> l'**insertion intrinsèque** : les joints sont partout dès l'origine, donc
> l'utilisateur n'a rien à « semer ». Elle est **compatible** avec ce que la
> session du 29/08 avait observé sur la transposition locale (zéro occurrence
> d'`adaptive`/`extrinsic`). Mais une phrase de page d'accueil n'est pas une
> preuve d'algorithme : à confirmer sur Guo 2014/2016 (§4.B ci-dessous).

> **[INFÉRÉ]** « Open Source » est revendiqué depuis VGW (2009, SourceForge).
> **Je n'ai pas pu vérifier si le Solidity actuel est réellement téléchargeable,
> ni sous quelle licence.** C'est une question à trancher — elle décide si une
> comparaison ligne à ligne rockim/Imperial est possible un jour, ou si tout
> devra passer par les publications. Voir §5, item **T1**.

---

## 3. BIBLIOGRAPHIE ANNOTÉE

Légende d'accès : **[LIBRE]** téléchargeable sans abonnement ·
**[OA?]** vraisemblablement en accès ouvert, à confirmer ·
**[SPIRAL?]** post-print probablement déposé au dépôt institutionnel d'Imperial ·
**[PAYWALL]** abonnement requis → accès Mines Paris-PSL ·
**[LIVRE]** achat ou bibliothèque.

### A. Fondations Munjiza — la formulation de base

| # | référence | ce que ça apporte | accès |
|---|---|---|---|
| A1 | **Munjiza, A., Andrews, K.R.F., White, J.K. (1999)**, « Combined single and smeared crack model in combined finite-discrete element analysis », *Int. J. Numer. Meth. Engng* **44**(1), 41-57 **[MÉTA]** | **La loi de joint cohésif d'origine** (mode I), calée sur les courbes σ-ε expérimentales du béton en traction **[RÉSUMÉ]**. C'est l'ancêtre direct de tout ce qui suit. | [PAYWALL] Wiley |
| A2 | **Munjiza, A., Andrews, K.R.F. (2000)**, « Penalty function method for combined finite-discrete element systems comprising large number of separate bodies », *IJNME* **49**(11), 1377-1396 **[MÉTA]** | **Le potentiel de contact.** Force de contact *distribuée*, cinématique préservant le bilan d'énergie **[RÉSUMÉ]**. C'est l'équation (2.49) que Guo cite (déjà en fiche au dépôt). | [PAYWALL] Wiley |
| A3 | **Munjiza, A., Andrews, K.R.F. (1998)**, « NBS contact detection algorithm for bodies of similar size », *IJNME* **43**(1), 131-149 **[MÉTA]** | **La détection de contact** NBS (*No Binary Search*), coût total ∝ N **[RÉSUMÉ]**. Pertinent pour le coût, pas pour la physique. | [PAYWALL] Wiley |
| A4 | **Munjiza, A., Latham, J.-P., Andrews, K.R.F. (2000)**, « Detonation gas model for combined finite-discrete element simulation of fracture and fragmentation », *IJNME* **49**(12), 1495-1520 **[MÉTA]** | Modèle de gaz de détonation. Hors sujet impact ; utile si le sujet dérive vers le tir. | [PAYWALL] Wiley |
| A5 | **Munjiza, A. (2004)**, *The Combined Finite-Discrete Element Method*, Wiley, 352 p., ISBN 978-0-470-02017-3 **[MÉTA]** | **Le livre de référence, avec le listing source de Y.** Seule source où l'on peut lire *à la fois* la formulation et le code qui la réalise, publiée et citable. | [LIVRE] |
| A6 | **Munjiza, A., Knight, E.E., Rougier, E. (2011)**, *Computational Mechanics of Discontinua*, Wiley, ISBN 978-1-119-97301-0 **[MÉTA]** | Suite du précédent, orientée 3D et parallélisme. | [LIVRE] |
| A7 | **Munjiza, A., Knight, E.E., Rougier, E. (2015)**, *Large Strain Finite Element Method: A Practical Course*, Wiley, 486 p., ISBN 978-1-118-40530-7 **[MÉTA]** | Cinématique grandes déformations du tétraèdre FDEM. Pertinent pour comprendre `CauchyTet4` et ses équivalents. | [LIVRE] |

### B. Le noyau 3D d'Imperial — la ligne la plus importante pour nous

| # | référence | ce que ça apporte | accès |
|---|---|---|---|
| **B1** | **Xiang, J., Munjiza, A., Latham, J.-P. (2009)**, « Finite strain, finite rotation quadratic tetrahedral element for the combined finite-discrete element method », *IJNME* **79**(8), 946-978, DOI 10.1002/nme.2599 **[MÉTA]** | **L'élément.** Motivé explicitement par le *locking* du tétraèdre linéaire, qui « dégrade sérieusement la précision » **[RÉSUMÉ]**. Fonde Y3D. | [PAYWALL] Wiley · [SPIRAL?] |
| **B2** | **Guo, L. (2014)**, *Development of a three-dimensional fracture model for the combined finite-discrete element method*, thèse de doctorat, Imperial College London **[MÉTA]** — handle Spiral annoncé `10044/1/28974` **[MÉTA, à vérifier]** | **LA pièce maîtresse.** C'est le modèle de rupture 3D de Solidity, écrit par son auteur. **F. Uzquiano l'a déjà** ; quatre sections sont en fiche au dépôt (§2.3.4-5, §2.4, ch. 5, ch. 6). **Les sections qui manquent sont listées au §6.** | [LIBRE] via Spiral — thèses Imperial en accès ouvert **[INFÉRÉ]** |
| **B3** | **Guo, L., Latham, J.-P., Xiang, J. (2015)**, « Numerical simulation of breakages of concrete armour units using a three-dimensional fracture model in the context of the combined finite-discrete element method », *Computers & Structures* **146**, 117-142 **[MÉTA]** | **La version revue par les pairs du modèle de rupture 3D.** Points forts annoncés : modèle dynamique transitoire explicite, interaction multi-corps 3D avec fissuration, collisions et fracturation d'objets fragiles de forme complexe **[RÉSUMÉ]**. **C'est la référence citable à la place de la thèse dans le manuscrit.** | [PAYWALL] Elsevier · [SPIRAL?] |
| **B4** | **Guo, L., Xiang, J., Latham, J.-P., Izzuddin, B. (2016)**, « A numerical investigation of mesh sensitivity for a new three-dimensional fracture model within the combined finite-discrete element method », *Engineering Fracture Mechanics* **151**, 70-91 **[MÉTA]** | **Sensibilité au maillage : taille ET orientation.** Fissure unique en traction, flexion trois points (maillages structurés) ; disque comprimé diamétralement sous divers angles (maillages non structurés) **[RÉSUMÉ]**. **[INFÉRÉ] C'est l'article correspondant au §2.4 de la thèse, dont [`guo2014_s24_maillage.md`](guo2014_s24_maillage.md) est déjà la fiche.** URL ScienceDirect sans `/abs/` → **[OA?]**. | [OA?] Elsevier · [SPIRAL?] |
| **B5** | **Guo, L., Xiang, J., Latham, J.-P., Izzuddin, B. (2020)**, « A generic computational model for three-dimensional fracture and fragmentation problems of quasi-brittle materials » **[MÉTA — revue, volume et pages non établis]** | Modèle générique : durcissement pré-pic, adoucissement post-pic, transition continu→discontinu, interaction explicite des fragments, code C/C++ maison ; validé par essais brésiliens et compression polyaxiale **[RÉSUMÉ]**. **PDF déposé librement : `discovery.ucl.ac.uk/id/eprint/10217490/1/Manuscript_UCL_deposit.pdf`.** | **[LIBRE]** dépôt UCL |

> **[INFÉRÉ]** B5 est signé des quatre mêmes auteurs que B4 et décrit « un code
> C/C++ maison » : c'est très probablement Solidity, quelques années après la
> thèse. Comme le manuscrit est **librement déposé**, c'est la pièce la moins
> chère à obtenir de tout le bloc B — et probablement la plus riche en équations
> accessibles. **À télécharger en premier.**

### C. Impact et forage percussif — la cible directe

| # | référence | ce que ça apporte | accès |
|---|---|---|---|
| **C1** | **Yang, X., Xiang, J., Latham, J.-P., Naderi, S., Wang, Y. (2023)**, « 3D numerical modelling of insert-rock interaction and piston-bit-rock interaction in percussive drilling using FDEM », 84e conf. annuelle EAGE, juin 2023, DOI 10.3997/2214-4609.202310940 **[MÉTA]** | **Le premier jalon de la série.** Résumé : la fissuration se sépare en deux temps — **amorçage et propagation des fissures de cisaillement en phase de CHARGE**, puis **propagation des radiales et des sub-horizontales qui tournent vers la surface libre en phase de DÉCHARGE** **[RÉSUMÉ, quasi-verbatim]**. | [PAYWALL] EarthDoc/EAGE |
| **C2** | **Yang, X., Xiang, J., Naderi, S., Wang, Y., Latham, J.-P., Aising, J., Gerbaud, L., Ugarte, I. (2024)**, « Where Does the Energy Go in Percussion Drilling? FDEM's Answer », 58e US Rock Mechanics/Geomechanics Symposium (ARMA 24), Golden, Colorado, juin 2024 **[MÉTA]** | **LE BILAN D'ÉNERGIE.** Pour le calcaire St Anne **et** le grès Rhune, sur un impact unique : **2,4-2,6 % de l'énergie seulement va à la propagation des fissures**, et **30 à 70 % est consommée par le FROTTEMENT ENTRE FRAGMENTS** **[RÉSUMÉ]**. | **[LIBRE]** — PDF public sur le site ORCHYD : `orchyd.eu/wp-content/uploads/2024/07/ARMA_Impact_Final.pdf` |
| **C3** | **Yang, X., Xiang, J., Latham, J.-P., Naderi, S., Wang, Y. (2025)**, « Cracking and fragmentation in percussive drilling: Insight from FDEM simulation », *J. Rock Mech. Geotech. Eng.* **17**(10), 6095-6110 **[MÉTA]** | **Granite Kuru Grey**, la même roche que l'article 2026. Fissures radiales, latérales, inclinées ; zones broyée et fissurée ; **coalescence de cratères adjacents produisant davantage d'éclats** **[RÉSUMÉ]**. | **[LIBRE]** — JRMGE est **intégralement en accès ouvert** |
| **C4** | **Yang, X., Xiang, J., Naderi, S., Wang, Y., Aising, J., Ugarte, I., Latham, J.-P. (2025)**, « Multi-criteria validation of hi-fidelity numerical model of impact breakage: towards next generation percussion drill simulation », *IJRMMS* **191**, 106125 **[MÉTA]** | **Les 7 critères de validation**, calcaire St Anne et grès Rhune. Deux aspects visés : perte d'énergie du taillant et morphologie de la rupture **[RÉSUMÉ]**. | [OA?] — URL ScienceDirect **sans** `/abs/`, et projet H2020 → mandat d'accès ouvert **[INFÉRÉ]** |
| **C5** | **Yang, X., Xiang, J., Wang, Y., Naderi, S., Aising, J., Ugarte, I., Latham, J.-P. (2026)**, *IJRMMS* **206**, 106660 — pulvérisation du granite Kuru Grey **[LU-DÉPÔT]** | **Déjà dépouillé** : [`yang2026_pulverisation.md`](yang2026_pulverisation.md). δ_m est une longueur ; 0,18 est le frottement de la roche ; le mécanisme dépend de la roche. **Ne pas refaire.** | déjà en possession (PDF) |
| **C6** | **Naderi, S., Latham, J.-P. (2025)**, « Optimised hammer drilling bit design using artificial neural networks trained by FDEM-generated data », *J. Rock Mech. Geotech. Eng.* **17**(11) **[MÉTA]** | **[INFÉRÉ]** Un réseau de neurones entraîné sur données FDEM suppose un **grand plan d'expériences FDEM** : c'est la pièce la plus susceptible de contenir un tableau complet de paramètres et de configurations de deck. | **[LIBRE]** — JRMGE, accès ouvert |
| **C7** | « A Study of Rock Breakage Under Extreme Submerged Confining Pressure: Can DTH Hammer Drilling Deliver? », *Rock Mech. Rock Eng.* (2025), DOI 10.1007/s00603-025-04626-1 **[MÉTA — liste d'auteurs non établie ; groupe Imperial/ORCHYD selon le contexte]** | Rupture sous pression de confinement immergée extrême — le régime *fond de puits* profond. Directement pertinent pour la géothermie profonde. | [PAYWALL] Springer · **[LIBRE?]** une version « Can DTH Hammer Drilling Deliver? » est annoncée sur orchyd.eu |
| C8 | **Xiang, X.(?), Naderi, S., Latham, J.-P. et al. (2024)**, « Destruction of Rock Microstructure: An Experimental and Numerical Modelling Study of High-Pressure Water Jet Rock Cutting Under Subsurface Confining Pressure Conditions », ARMA 24 **[MÉTA]** | Volet jet d'eau d'ORCHYD. Hors percussion, mais même équipe, même code. | [PAYWALL] OnePetro · [LIBRE?] orchyd.eu |

### D. Le projet ORCHYD — la mine à ciel ouvert

**[MÉTA + RÉSUMÉ]** ORCHYD (« Novel Drilling Technology Combining Hydro-Jet and
Percussion for ROP Improvement in deep geothermal drilling »), **H2020, projet
n° 101006752**, fiche CORDIS `cordis.europa.eu/project/id/101006752`.

* **Coordinateur : ARMINES / Mines ParisTech.** Partenaires : **Imperial College
  London**, SINTEF, Université du Pirée, China University of Petroleum,
  Drillstar. **L. Gerbaud** (Mines Paris) est co-auteur de C2.
* Objectif : porter la vitesse d'avancement en roche dure de 1-2 m/h à 4-10 m/h,
  jusqu'à 6 km. Résultat annoncé : **+170 % de vitesse dans le granite du
  Sidobre** avec assistance HPWJ **[RÉSUMÉ]**.
* **WP6 = destruction de la roche par percussion.** **WP5 = jet d'eau.**
* Livrables publics repérés :
  * **D7.1** « Report on optimization of the simultaneous interactions during
    HPWJ-percussive drilling »,
    `orchyd.eu/wp-content/uploads/2025/02/D7.1_...pdf` **[LIBRE]** ;
    renvoie à **D6.2** pour la configuration percussion.
  * page **`orchyd.eu/repository/`** — le catalogue des livrables **[LIBRE]**.

> **[INFÉRÉ] — et c'est peut-être le point le plus opérationnel de cette fiche.**
> ORCHYD est un projet **H2020 coordonné par ARMINES/Mines ParisTech**, c'est-à-dire
> **l'institution de F. Uzquiano**. Deux conséquences :
> 1. **Mandat d'accès ouvert H2020** : toutes les publications à comité de
>    lecture issues d'ORCHYD *doivent* être déposées en libre accès (voie dorée
>    ou verte). Les C1-C5 et C7-C8 sont donc probablement récupérables
>    **gratuitement**, sur `orchyd.eu`, sur Spiral, ou via CORDIS.
> 2. **Les livrables de WP6 (rapports internes détaillés, non publiés en revue)
>    peuvent être accessibles en interne à Mines Paris.** Un rapport de livrable
>    contient typiquement bien plus de détail numérique qu'un article : decks,
>    tableaux de paramètres, études de sensibilité. **C'est la piste la plus
>    prometteuse et personne ne l'avait ouverte.**

### E. Le logiciel lui-même

| # | pièce | statut |
|---|---|---|
| **T1** | `solidityproject.com` — site officiel : pages *Technology* (`/3d-femdem/`, `/fracture-modelling/`, `/particle-shape-library/`), *Applications*, **`/recent-publications/`** | **[LIBRE]**, mais **inaccessible depuis ce conteneur**. La page *Recent Publications* est le catalogue faisant autorité de ce que le groupe publie : **à relever en entier**. |
| **T2** | `sourceforge.net/projects/vgw/` — **Virtual Geoscience Workbench**, source ouverte | **[LIBRE]** — **la seule source de CODE publiée et citable de cette lignée** à laquelle je puisse renvoyer avec certitude. Millésime 2009 : antérieur au modèle de rupture 3D de Guo (2014). Utile pour le contact et le noyau, pas pour la loi de joint 3D actuelle. |
| T3 | **Xiang, J., Latham, J.-P., Farsi, A. (2017)**, « Algorithms and Capabilities of Solidity to Simulate Interactions and Packing of Complex Shapes », in *Proc. 7th Int. Conf. on Discrete Element Methods (DEM7, Dalian, août 2016)*, Springer Proc. in Physics **188**, ch. 16 **[MÉTA]** | [PAYWALL] Springer. **[INFÉRÉ]** Un chapitre intitulé « Algorithms and Capabilities » est le document le plus proche d'une **description d'architecture** publiée. Prioritaire pour le lot 3. |
| T4 | **Latham, J.-P., Munjiza, A., Garcia, X., Xiang, J., Guises, R. (2008/2009)**, « The Virtual Geoscience Workbench, VGW: Open Source tools for discontinuous systems », *Particuology* **[MÉTA — auteurs et pagination à confirmer]** | [PAYWALL] Elsevier. Article de présentation de la plateforme. |

### F. Le contexte hors-Imperial (à connaître, à ne pas confondre)

**[MÉTA]**, sans annotation détaillée — ces pièces servent à comprendre la
famille, pas à documenter Imperial :

* **LANL/HOSS** : Rougier et al. *IJRMMS* 70 (2014) 101-108 (SHPB) ; Knight,
  Rougier et al.
* **Toronto/Y-Geo** : Mahabadi et al., « Y-Geo: New Combined Finite-Discrete
  Element Numerical Code for Geomechanical Applications », *Int. J. Geomech.*
  **12**(6) ; Mahabadi et al., « Y-GUI ... incorporating material
  heterogeneity » ; Lisjak & Grasselli (revues FDEM roche).
* **Branche chinoise / GPGPU** : Fukuda et al., *RMRE* (2019/2020), simulateur
  3D FDEM parallélisé GPGPU ; Liu et al., *IJNAMG* (2020), multi-GPU CUDA ;
  « An Improved GPU-Parallelized 2D/3D Elastoplastic-Damage-Fracture Joint
  Framework », *RMRE* (2023).
* **Contexte forage percussif non-FDEM** : Saksala (modélisation continue,
  bit-rock, Kuru granite) ; Hustrulid ; les essais de chute sur granite de Kuru
  (`pmc.ncbi.nlm.nih.gov/articles/PMC5179971/`, **[LIBRE]** — c'est
  vraisemblablement la source expérimentale que Yang et al. utilisent pour
  valider ; **à récupérer, c'est gratuit**).

---

## 4. CE QUI N'EST PAS PUBLIÉ — les trous, nommés

**[INFÉRÉ]** à partir de ce que les fiches déjà au dépôt ont établi de première
main et de la structure de la bibliographie ci-dessus :

1. **L'algorithme de frottement tangentiel.** La thèse de Guo y consacre **trois
   phrases** et attribue l'implémentation à J. Xiang, **sans la publier**
   (`guo2014_s234_235_contact_integration.md` §1, **[LU-DÉPÔT]**). Aucune des
   pièces repérées ci-dessus ne promet de combler ce trou. **[SUPPOSÉ]** T3
   (« Algorithms and Capabilities ») est le meilleur candidat.
2. **Le couplage endommagement → contact** (raideur normale ? frottement ? les
   deux ?). L'article C5 décrit l'**objectif** — « severe local stiffness
   degradation, loss of load-bearing capacity » — mais la session précédente n'a
   trouvé **aucune équation publiée** qui dise *où* le couplage se branche.
   **C'est le trou le plus coûteux** : c'est exactement la question qui a piégé
   la session du 29/08.
3. **Le retrait des fragments et la mesure de la masse détachée.** Aucune pièce
   repérée n'en fait son sujet. **[SUPPOSÉ]** décrit incidemment dans C4 (la
   masse de fragments est un des 7 critères) ou dans un livrable ORCHYD.
4. **Le DIF.** Aucune pièce Imperial repérée ne le traite frontalement.
   **[SUPPOSÉ]** soit il n'y en a pas dans leur modèle d'impact, soit il est
   mentionné en passant dans C3/C4/C5. À vérifier avant d'affirmer quoi que ce
   soit — **ne pas conclure de son absence bibliographique qu'il est absent du
   code**.
5. **Le maillage adaptatif en cours de calcul.** Le brief (§4) le suppose dans la
   thèse de Guo. **[INFÉRÉ]** Le titre de B4 est « mesh **sensitivity** », pas
   « mesh adaptivity » : les deux termes ne disent pas la même chose, et la
   fiche `guo2014_s24_maillage.md` décrit une **étude de sensibilité à h
   fixe par maillage**, pas un raffinement adaptatif. **Il se peut qu'il n'y
   ait pas de maillage adaptatif du tout.** À trancher — voir §6, demande **G3**.

---

## 5. LISTE DE TÉLÉCHARGEMENT, PAR ORDRE DE VALEUR

*Pour F. Uzquiano. Colonne « coût » : ce qu'il en coûte réellement d'obtenir la
pièce.*

### Rang 1 — gratuit et à très fort rendement (à faire d'abord, ~30 min)

| # | pièce | pourquoi (une ligne) | coût |
|---|---|---|---|
| 1 | **B5** — Guo et al. 2020, PDF `discovery.ucl.ac.uk/id/eprint/10217490/1/Manuscript_UCL_deposit.pdf` | Le modèle de rupture 3D d'Imperial, en manuscrit libre : la source d'équations la moins chère du lot. | gratuit |
| 2 | **C2** — ARMA 2024, PDF `orchyd.eu/wp-content/uploads/2024/07/ARMA_Impact_Final.pdf` | Le bilan d'énergie sur *votre* roche (St Anne) : 2,4-2,6 % aux fissures, 30-70 % au frottement entre fragments. | gratuit |
| 3 | **C3** — JRMGE 17(10) 6095-6110, `sciencedirect.com/science/article/pii/S1674775525001477` | Revue en accès ouvert intégral, granite Kuru, fissures et éclats : la formulation y est forcément rappelée. | gratuit |
| 4 | **T1** — page `solidityproject.com/recent-publications/` (+ `/technology/`) | Le catalogue faisant autorité : il fermera les trous de cette bibliographie en une seule visite. | gratuit |
| 5 | **D** — `orchyd.eu/repository/`, tous les livrables **WP6** et **D6.2** | Un livrable de projet contient les decks et les tableaux qu'un article comprime ; et vous êtes chez le coordinateur. | gratuit |
| 6 | **C6** — Naderi & Latham 2025, JRMGE 17(11) | Un réseau de neurones entraîné sur FDEM implique un plan d'expériences documenté : tableaux de paramètres. | gratuit |
| 7 | **B2** — thèse Guo 2014 sur Spiral (`10044/1/28974`, à vérifier) | Vous l'avez déjà, mais la version Spiral permet de citer un handle stable dans le manuscrit. | gratuit |
| 8 | essais de chute sur granite de Kuru, `pmc.ncbi.nlm.nih.gov/articles/PMC5179971/` | La source expérimentale probable de leur validation granite : indispensable pour juger de leur accord. | gratuit |

### Rang 2 — accès Mines Paris-PSL, valeur élevée

| # | pièce | pourquoi (une ligne) | éditeur |
|---|---|---|---|
| 9 | **B3** — Guo, Latham, Xiang 2015, *Comput. Struct.* **146**, 117-142 | Le modèle de rupture 3D revu par les pairs : la référence à citer dans la thèse, et le texte des équations. | Elsevier |
| 10 | **T3** — Xiang, Latham, Farsi 2017, *Springer Proc. Phys.* **188**, ch. 16 | Le seul document publié qui promette une description d'**architecture** de Solidity — dont, peut-être, le frottement. | Springer |
| 11 | **C4** — Yang et al. 2025, *IJRMMS* **191**, 106125 | Les 7 critères de validation sur **St Anne**, votre roche : le protocole de comparaison à reproduire. | Elsevier (tenter l'OA d'abord) |
| 12 | **B4** — Guo et al. 2016, *Eng. Fract. Mech.* **151**, 70-91 | Sensibilité taille **et orientation** de maillage : la version citable du §2.4 déjà en fiche. | Elsevier (tenter l'OA d'abord) |
| 13 | **B1** — Xiang, Munjiza, Latham 2009, *IJNME* **79**(8), 946-978 | L'élément tétraédrique quadratique : dit si leur souplesse vient de l'élément ou des joints. | Wiley |
| 14 | **A1** — Munjiza, Andrews, White 1999, *IJNME* **44**(1), 41-57 | La loi de joint cohésif d'origine, dont tout le reste est une variante. | Wiley |
| 15 | **A2** — Munjiza & Andrews 2000, *IJNME* **49**(11), 1377-1396 | Le potentiel de contact : la seule façon de vérifier ce que « pénalité » veut dire chez eux. | Wiley |
| 16 | **C1** — Yang et al., EAGE 2023, DOI 10.3997/2214-4609.202310940 | La chronologie charge/décharge des fissures : cisaillement à la charge, radiales à la décharge. | EAGE |
| 17 | **C7** — *RMRE* 2025, DOI 10.1007/s00603-025-04626-1 | Rupture sous confinement immergé extrême : le régime réel du fond de puits géothermique. | Springer |

### Rang 3 — livres, si le budget le permet

| # | pièce | pourquoi |
|---|---|---|
| 18 | **A5** — Munjiza 2004, *The Combined FDEM*, Wiley | Contient le **listing source de Y** : formulation *et* code, publiés et citables. La seule pièce qui puisse jouer le rôle que `/home/user/solidity` ne peut pas jouer. |
| 19 | **A7** — Munjiza, Knight, Rougier 2015, *Large Strain FEM* | La cinématique grandes déformations du tétraèdre FDEM. |
| 20 | **A6** — Munjiza, Knight, Rougier 2011, *Comput. Mech. of Discontinua* | Le 3D et le parallélisme. |

**Ne pas acheter avant d'avoir fait le rang 1 :** il est très possible que les
pages `solidityproject.com/technology/` et les livrables ORCHYD rendent une
partie du rang 2 inutile.

---

## 6. DEMANDES PRÉCISES SUR LA THÈSE DE GUO (2014)

*Le brief autorise à demander la thèse **par extraits**, une section à la fois.
Quatre sections sont déjà en fiche ; voici les suivantes, par ordre de besoin.*

| # | section demandée | pourquoi j'en ai besoin | pour quel lot |
|---|---|---|---|
| **G1** | **§2.2 (ou la section « fracture model » / « joint elements »), en entier** — la définition de la loi de joint cohésive : contrainte normale et tangentielle en fonction des ouvertures, forme de la courbe d'adoucissement, définition des plages, critère de rupture mode I / mode II / mixte, **avec les numéros d'équation** | C'est **le** trou central du lot 2. Sans ça, je ne peux rien affirmer sur la loi cohésive d'Imperial. | **lot 2** |
| **G2** | **§2.3.1 à §2.3.3** — l'insertion des joints : à quel moment ils sont créés, dans quels corps, quelle **pénalité élastique** leur est donnée, et **si les auteurs corrigent le module d'Young pour compenser la souplesse artificielle** | Question directe du lot 3, et question ouverte du dépôt depuis le 28/08. | **lot 3** |
| **G3** | **la table des matières complète**, plus le **titre de chaque section des chapitres 3 et 4** | Pour trancher s'il existe un **maillage adaptatif** (le brief le suppose ; B4 s'intitule « mesh *sensitivity* », ce qui n'est pas la même chose), et pour repérer les validations dynamiques transformables en bancs. | **lots 3 et 5** |
| **G4** | toute section où apparaissent les mots **« damage »**, **« friction coefficient »**, **« erosion »**, **« fragment »** hors §2.3.4 — même une simple liste de numéros de page ferait l'affaire dans un premier temps | Pour savoir si le couplage endommagement→contact et le retrait des fragments existent dans la thèse, avant de les chercher ailleurs. | **lots 2 et 4** |

**Le plus urgent est G1.** Si vous ne devez en envoyer qu'un seul, envoyez G1.

Un raccourci possible : **B5 est en accès libre** et signé de Guo. Si sa
section formulation reprend celle de la thèse — ce qui est probable — G1 devient
inutile et vous économisez l'extraction. **Ouvrez B5 d'abord.**

---

## 7. UNE CORRECTION À PORTER SUR UNE FICHE EXISTANTE

[`yang2026_pulverisation.md`](yang2026_pulverisation.md) **§5**, dernier
paragraphe, affirme :

> « Le code Solidity le realise par `penalty *= min(1-D_i, 1-D_j)` et
> `mu = mud*d_fact` (Y3Did.c l. 995, 1044, 1263-1265) : raideur de contact ET
> frottement suivent (1-D) EN CONTINU. »

Cette fiche a été écrite **avant** la CORRECTION 2 du 29/08. Ces trois lignes
attribuent à Solidity ce qui a été lu dans `/home/user/solidity`, lequel n'est
pas le code d'Imperial et dont le facteur d'endommagement est câblé à zéro.
**Le reste de la fiche (§1 à §4, §6) porte sur l'article : il tient.** Seul ce
paragraphe est à requalifier.

Correctif appliqué : un encart d'avertissement a été inséré en tête du §5 de
cette fiche, renvoyant ici et à la CORRECTION 2. Le corps du texte est conservé
(on ne réécrit pas l'historique).

---

## 8. ÉTAT DU LOT 1 ET SUITE

**Livrable 1 : fait**, dans les limites posées au §0 — 20 pièces référencées,
classées, avec statut d'accès et liste de téléchargement ordonnée.

**Ce que le lot 1 ne peut pas donner et qu'il ne prétend pas donner** : une
seule équation d'Imperial. Le lot 2 est **bloqué** tant qu'au moins une des
pièces suivantes n'est pas en main : **B5** (gratuit), **B3**, ou l'extrait
**G1** de la thèse.

**Recommandation d'enchaînement** : faire le rang 1 de la liste (§5) — c'est
gratuit, ça prend une demi-heure, et il est plausible que ça débloque le lot 2
en entier sans dépenser un euro ni extraire une seule page de la thèse.
