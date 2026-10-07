# Anisotropie (schistosité / stratification) en FDEM — bibliographie

**Revue du 2026-09-01.** Six angles de recherche indépendants + vérification
adverse des références principales (consigne : *réfuter par défaut* si
l'existence ou le mécanisme ne peut pas être confirmé).

> **Statut de vérification.** Seules les entrées marquées ✅ ont été soumises à
> une vérification indépendante (existence de l'article + conformité du
> mécanisme décrit). Les autres ont été rapportées par les agents de recherche
> mais **non re-vérifiées** : à confirmer vous-même avant citation.

---

## A. Le noyau — la lignée Lisjak / Grasselli / Vietor ✅

C'est la référence pour l'argilite d'Opalinus (Mont Terri) et l'anisotropie en
FDEM. C'est aussi ce que Wang et al. 2024 citent pour « le rôle clé des plans
de stratification dans l'EDZ ». **Les quatre entrées ci-dessous sont vérifiées.**

**A1. Lisjak Bradley, A. (2013, dépôt 2014) — thèse de doctorat, University of
Toronto.**
*Investigating the Influence of Mechanical Anisotropy on the Fracturing
Behaviour of Brittle Clay Shales with Application to Deep Geological
Repositories.*
→ **Accès libre :** https://utoronto.scholaris.ca/bitstreams/49c1c1ee-83bc-4b7b-8fcf-74cc096e343a/download
→ **Déjà téléchargée en local** (56 Mo) — voir §Fichiers locaux plus bas.

> **C'est LE document à lire.** Il contient les trois briques du modèle, et
> surtout la justification écrite de l'abandon de l'approche à plans discrets :
> *« due to computational constraints associated with the mesh resolution, the
> technique was found unsuitable for field-scale models »*.
> Les chapitres 4, 5 et 6 correspondent respectivement à A2, A3 et à
> l'application. Le résumé méthodologique de l'auteur : anisotropie traitée par
> *« (i) smearing the transversely isotropic elastic deformation […] (ii)
> preconditioning the triangular mesh along the bedding plane direction, and
> (iii) developing a directional cohesive fracture formulation »*.

**A2. Lisjak, Tatone, Grasselli & Vietor (2014)** — *Rock Mech. Rock Eng.*
47(1), 187–206.
*Numerical Modelling of the Anisotropic Mechanical Behaviour of Opalinus Clay
at the Laboratory-Scale Using FEM/DEM.*
→ https://doi.org/10.1007/s00603-012-0354-7
*L'approche **discrète** (échelle labo) — celle qui a ensuite été abandonnée
pour l'échelle de l'ouvrage. À noter : elle incluait **déjà** l'élasticité
transversalement isotrope.*

**A3. Lisjak, Grasselli & Vietor (2014)** — *Int. J. Rock Mech. Min. Sci.*
65, 96–115.
*Continuum–discontinuum analysis of failure mechanisms around unsupported
circular excavations in anisotropic clay shales.*
→ https://doi.org/10.1016/j.ijrmms.2013.10.006
*L'approche **pervasive** (« smeared ») appliquée à une excavation circulaire.
C'est le pivot méthodologique. Contient la phrase sur l'absence d'anisotropie
élastique dans la pénalité des joints.*

**A4. Lisjak, Garitte, Grasselli, Müller & Vietor (2015)** — *Tunn. Undergr.
Space Technol.* 45, 227–248.
*The excavation of a circular tunnel in a bedded argillaceous rock (Opalinus
Clay): Short-term rock mass response and FDEM numerical analysis.*
→ https://doi.org/10.1016/j.tust.2014.09.014
*Le cas d'ouvrage réel (galerie FE, Mont Terri, tunnel 3 m) confronté aux
mesures in situ. C'est le seul cas d'ouvrage confirmé, et la source des ordres
de grandeur d'espacement (t ≈ 10 cm pour un tunnel de 3 m).*

### Compléments de la même équipe *(non re-vérifiés)*

- **Lisjak, Figi & Grasselli (2014)** — *J. Rock Mech. Geotech. Eng.* 6(6),
  493–505. *Fracture development around deep underground excavations: Insights
  from FDEM modelling.* → https://doi.org/10.1016/j.jrmge.2014.09.003
  **(accès libre — PDF déjà en local)**
- **Lisjak, Tatone, Mahabadi, Grasselli, Marschall, Lanyon, de La Vaissière,
  Shao, Leung & Nussbaum (2016)** — *Rock Mech. Rock Eng.* 49(5), 1849–1873.
  *Hybrid Finite-Discrete Element Simulation of the EDZ Formation and
  Mechanical Sealing Process Around a Microtunnel in Opalinus Clay.*
  → https://doi.org/10.1007/s00603-015-0847-2
- **Mahabadi, Kaifosh, Marschall & Vietor (2014)** — *J. Rock Mech. Geotech.
  Eng.* 6(6), 591–606. *Three-dimensional FDEM numerical simulation of failure
  processes observed in Opalinus Clay laboratory samples.*
  → https://doi.org/10.1016/j.jrmge.2014.10.005 **(accès libre)**

---

## B. FDEM anisotrope — autres équipes *(non re-vérifiées sauf mention)*

**B1. Deng, Liu, Huang, Pan & Wu (2022)** — *Comput. Geotech.* 142, 104535.
*FDEM numerical modeling of failure mechanisms of anisotropic rock masses
around deep tunnels.* → https://doi.org/10.1016/j.compgeo.2021.104535
✅ *existence confirmée ; mécanisme **en doute** — les agents ne s'accordent pas
sur pervasif vs plans discrets. À trancher vous-même : c'est le cas le plus
proche du vôtre (tunnel profond, masse rocheuse anisotrope).*

**B2. Liu, Liu, Deng, Huang & Xie (2023)** — *Eng. Fract. Mech.* 289, 109359.
*A new phenomenological anisotropic tensile failure criterion and its
application in FDEM simulations.*
→ https://doi.org/10.1016/j.engfracmech.2023.109359

**B3. Liu, Liu, Huang, Hu, Bo, Yuan & Xie (2022)** — *Rock Mech. Rock Eng.*
55(12), 7765–7789. *Direct Tensile Test and FDEM Numerical Study on Anisotropic
Tensile Strength of Kangding Slate.*
→ https://doi.org/10.1007/s00603-022-03036-x
*Ardoise — donc une vraie schistosité, pas seulement un litage sédimentaire.*

**B4. Li, Chapman, Faramarzi & Metje (2024)** — *Rock Mech. Rock Eng.* 57,
2385–2405. *The Analysis of the Fracturing Mechanism and Brittleness
Characteristics of Anisotropic Shale Based on Finite-Discrete Element Method.*
→ https://doi.org/10.1007/s00603-023-03672-x

**B5. Zhang, Qiu, Jiang, Zheng, Xie & Fang (2025)** — *Rock Mech. Rock Eng.*
58, 4139–4158. *Study of the Mechanical Characteristics and Crack Evolution of
Layered Rocks Using Voronoi Block-Based Finite-Discrete Element Method.*
→ https://doi.org/10.1007/s00603-024-04372-w
*Voronoï + FDEM — proche de votre machinerie GBM existante.*

**B6. Sun, Liu, Grasselli & Tang (2020)** — *Comput. Geotech.* 117, 103237.
*Simulation of thermal cracking in anisotropic shale formations using the
combined finite-discrete element method.*
→ https://doi.org/10.1016/j.compgeo.2019.103237
*⚠️ Croise directement votre chantier thermique de ce matin.*

**B7. He, Liu & Deng (2020)** — *Adv. Mater. Sci. Eng.* 2020, 8793214.
*Investigation of the Anisotropic Characteristics of Layered Rocks under
Uniaxial Compression Based on the 3D Printing Technology and the Combined
Finite-Discrete Element Method.*
→ https://doi.org/10.1155/2020/8793214 **(accès libre)**

**B8. Seyed Ghafouri, Aboayanah, Abdelaziz & Grasselli (2022)** — ARMA
2022-0310, 56ᵉ US Rock Mechanics Symposium. *Numerical Investigation of the
Influence of Bedding Plane Thickness and Friction on Cracking Pattern and
Mechanical Behavior of Shale Under Unconfined Loading Condition Using FDEM.*
*Traite explicitement l'**épaisseur** des plans de litage — votre question
d'espacement.*

### École chinoise (Yan / Jiao, MultiFracS — les auteurs de votre article)

**B9. Hu, Yan, Jiao, Wang, Jia & Wang (2024/2025)** — *Eng. Fail. Anal.*
*Analysis and countermeasures of asymmetric failure in layered surrounding rock
tunnels based on FDEM: A case study.*
→ https://doi.org/10.1016/j.engfailanal.2024.109049
*⚠️ **Tunnel en roche stratifiée, rupture asymétrique** — c'est exactement votre
question, par l'équipe même de Wang et al. 2024. Priorité haute.*

**B10. Guo, Yan, Zhang, Xu, Wang & Jiao (2024)** — *Comput. Geotech.* 165,
105883. *Mechanical analysis of toppling failure using FDEM: A case study for
soft-hard interbedded anti-dip rock slope.*
→ https://doi.org/10.1016/j.compgeo.2023.105883

**B11. Shi, Zheng, Kong, Luo & Chen (2022)** — *Eng. Fract. Mech.* 272, 108718.
*Study on the failure mechanism in shale-sand formation based on hybrid
finite-discrete element method.*
→ https://doi.org/10.1016/j.engfracmech.2022.108718

**B12. Wang, Gan, Wang, Ma, Yan, Benson, Wang & Elsworth (2024)** —
*Geomech. Geophys. Geo-energ. Geo-resour.* 10, 71. *Propagation and complex
morphology of hydraulic fractures in lamellar shales based on FDEM.*
→ https://doi.org/10.1007/s40948-024-00788-4

---

## C. Ubiquitous joint / continuum équivalent *(non re-vérifiées)*

L'origine théorique de l'approche pervasive, et ses usages en continuum.

**C1. Jaeger, J.C. (1960)** — *Geological Magazine* 97(1), 65–72.
*Shear failure of anisotropic rocks.*
*L'article fondateur du plan de faiblesse. C'est la courbe en U.*

**C2. Sainsbury, B.L. & Sainsbury, D.P. (2017)** — *Rock Mech. Rock Eng.* 50,
1507–1528. *Practical Use of the Ubiquitous-Joint Constitutive Model for the
Simulation of Anisotropic Rock Masses.*
→ https://doi.org/10.1007/s00603-017-1177-3
*Les pièges pratiques du modèle ubiquitaire — utile pour savoir ce qu'on
achète et ce qu'on perd.*

**C3. Leng, Wang, Sheng, Chen & Li (2021)** — *Front. Earth Sci.* 9, 744900.
*An Enhanced Ubiquitous-Joint Model for a Rock Mass With Conjugate Joints and
Its Application on Excavation Simulation of Large Underground Caverns.*
→ https://doi.org/10.3389/feart.2021.744900 **(accès libre)**

**C4. Vazaios, Swan & Stewart (2022)** — Tunnelling Association of Canada.
*The Implementation of the Ubiquitous Joint Model to Capture Small Scale
Defects in Sedimentary Rock Masses for Underground Excavations.*
*Traite précisément l'argument d'échelle : défauts sous-maille.*

---

## D. DEM / bonded-particle transversalement isotrope *(comparaison, non re-vérifiées)*

Comment la communauté **discrète non-FDEM** traite le même problème.

**D1. Park, B. & Min, K.-B. (2015)** — *Int. J. Rock Mech. Min. Sci.* 76,
243–255. *Bonded-particle discrete element modeling of mechanical behavior of
transversely isotropic rock.*
→ https://doi.org/10.1016/j.ijrmms.2015.03.014

**D2. Park, Min, Thompson & Horsrud (2018)** — *Int. J. Rock Mech. Min. Sci.*
110, 120–132. *Three-dimensional bonded-particle discrete element modeling of
mechanical behavior of transversely isotropic rock.*
→ https://doi.org/10.1016/j.ijrmms.2018.07.018

**D3. Rasmussen, L.L. & Min, K.-B. (2023)** — *Int. J. Rock Mech. Min. Sci.*
170, 105518. *Developments to the Bonded Block Modeling technique for Discrete
Element simulation of transversely isotropic rocks.*
→ https://doi.org/10.1016/j.ijrmms.2023.105518

**D4. Gao, F. (2013)** — thèse, Simon Fraser University. *Simulation of failure
mechanisms around underground coal mine openings using discrete element
modelling.* *(approche UDEC Trigon — l'analogue polygonal)*

---

## E. Revues récentes *(non re-vérifiées — utiles pour le panorama)*

**E1. Li, Feng, Su, Chen & Guo (2025)** — *Int. J. Coal Sci. Technol.* 12, 99.
*A review of computational approaches for simulating fracturing mechanisms in
layered rock formations.*
→ https://doi.org/10.1007/s40789-025-00836-8 **(accès libre)**
*⚠️ Lue intégralement par un agent mais **non soumise à la vérification
indépendante**. Panorama a priori exactement sur le sujet.*

**E2. Ghorbani, Shahinfar & Taheri (2025)** — *Deep Resources Engineering*,
100219. *A review of the geological characterization, classification, modeling,
and case studies of anisotropic rock masses.*
→ https://doi.org/10.1016/j.deepre.2025.100219 **(CC-BY)**

---

## F. Hors littérature revue par les pairs

**F1. Geomechanica Inc. — « Modelling approaches for layered rocks ».**
→ https://www.geomechanica.com/blog/modelling-approaches-layered-rocks/
*Billet technique de l'éditeur d'**Irazu**, le code FDEM commercial issu de la
lignée Y-Geo/Grasselli. **Non revu par les pairs**, mais c'est la doc de ceux
qui ont écrit A1–A4 : utile pour savoir ce qui a survécu en production.*

---

## Fichiers déjà téléchargés en local

Dans le répertoire de travail temporaire de la session :

| fichier | contenu |
|---|---|
| `lisjak_thesis.pdf` (56,7 Mo) | **A1 — la thèse intégrale** |
| `lisjak_thesis.txt` / `lisjak.txt` (840 ko) | texte extrait, cherchable |
| `lisjak2014_jrmge.pdf` (833 ko) | Lisjak, Figi & Grasselli 2014, accès libre |

Chemin :
`C:\Users\FUZQUI~1\AppData\Local\Temp\claude\C--Users-fuzquianoalricabi-simulations\7e7c5111-23f3-4e3a-87dc-8684a34311a8\scratchpad\`

⚠️ Ce répertoire est **temporaire**. Dites-le-moi si vous voulez que je copie la
thèse quelque part de durable (`OVERLEAF\`, la base OneDrive, ou ici).

---

## Ordre de lecture suggéré

1. **A1, chapitres 5 et 6** — le pivot méthodologique et l'abandon du discret.
   Si vous ne lisez qu'une chose, c'est ça.
2. **A3** — la formulation pervasive sur excavation circulaire.
3. **A4** — le cas d'ouvrage réel, les ordres de grandeur d'espacement.
4. **B9 (Hu/Yan 2024)** — tunnel stratifié asymétrique, par l'équipe de votre
   article de référence.
5. **A2** — pour voir ce que faisait l'approche discrète *avant* l'abandon,
   et pourquoi elle incluait quand même l'élasticité TI.

---

## Ce que je n'ai PAS trouvé

Aucune référence confirmée n'emploie l'option **(C)** — champ aléatoire à
corrélation anisotrope — pour représenter une schistosité. C'est une absence de
preuve, pas une preuve d'absence, mais aucun candidat n'a survécu à la
vérification. Le seul usage voisin trouvé est **Deng, Liu & Lu (2022)**,
*Comput. Geotech.* 154, 105138, qui utilise une distribution de Weibull pour
l'**hétérogénéité** — pas pour l'anisotropie.
→ https://doi.org/10.1016/j.compgeo.2022.105138

Et aucune référence confirmée ne combine **plans discrets + élasticité
isotrope**, la configuration du deck actuel.
