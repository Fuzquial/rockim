# Les sources — références complètes et statut d'accès

**Aucun PDF n'est inclus dans ce dossier.** Ces articles sont sous droits
(Elsevier, Springer, Wiley, EAGE, ARMA) et leur redistribution serait illicite.
Chaque entrée porte son DOI ou son identifiant pour que le lecteur en obtienne sa
propre copie.

## Les neuf sources primaires dépouillées

| code | référence | accès |
|---|---|---|
| **UCL** | Guo, L., Xiang, J., Latham, J.-P., Izzuddin, B., « A generic computational model for three-dimensional fracture and fragmentation problems of quasi-brittle materials » — **la formulation, éq. 1-20** | manuscrit déposé, **UCL Discovery** eprint 10217490, **libre** |
| **XLF** | Xiang, J., Latham, J.-P., Farsi, A., « Algorithms and Capabilities of Solidity to simulate interactions and packing of complex shapes », *Springer Proc. Phys.* **188**, ch. 16 (DEM7, Dalian 2016), paper G010111 — **la loi de frottement** | DOI 10.1007/978-981-10-1926-5_16 · aussi *Simulations in Bulk Solids Handling* ch. 8, DOI 10.1002/9783527835935.ch8 |
| **XMLG** | Xiang, J., Munjiza, A., Latham, J.-P., Guises, R., « On the validation of DEM and FEM/DEM models in 2D and 3D », *Engineering Computations* **26**(6) (2009) 673-687 — **publication d'origine du frottement, éq. 8-9 p. 677** | DOI 10.1108/02644400910975469 |
| **IJ191** | Yang, X., Xiang, J., Naderi, S., Wang, Y., Aising, J., Ugarte, I., Latham, J.-P., « Multi-criteria validation of hi-fidelity numerical model of impact breakage », *IJRMMS* **191** (2025) 106125 — **St Anne + Rhune, 7 critères** | Elsevier · projet H2020, accès ouvert probable |
| **IJ206** | Yang, X., Xiang, J., Naderi, S., Wang, Y., Aising, J., Ugarte, I., Latham, J.-P., « High-fidelity modelling of fragmentation and pulverisation in hard granite under percussion loading: a FDEM-based approach », *IJRMMS* **206** (2026) 106660 — **la pulvérisation** | Elsevier |
| **JR** | Yang, X., Xiang, J., Latham, J.-P., Naderi, S., Wang, Y., « Cracking and fragmentation in percussive drilling: Insight from FDEM simulation », *J. Rock Mech. Geotech. Eng.* **17**(10) (2025) 6095-6110 — **Kuru Grey** | **JRMGE, accès ouvert intégral** |
| **ANN** | Naderi, S., Wang, Yuyang, Yang, X., Xiang, J., Pain, C., Heaney, C., Gerbaud, L., Velmurugan, N., Latham, J.-P., « Optimised hammer drilling bit design using artificial neural networks trained by FDEM-generated data », *JRMGE* **17**(11) (2025) — **modèle 2D** | **JRMGE, accès ouvert intégral** |
| **A952** | Yang, X., Xiang, J., Naderi, S., Wang, Y., Latham, J.-P., Aising, J., Gerbaud, L., Ugarte, I., « Where does the energy go in percussion drilling? FDEM's answer », **ARMA 24-0952**, 58e US Rock Mech./Geomech. Symposium, Golden, juin 2024 — **le bilan d'énergie sur St Anne** | OnePetro |
| **A788** | Gerbaud, L., Velmurugan, N., Aising, J., Chambres, C. (Mines Paris – PSL) ; Naderi, S., Latham, J.-P., Xiang, J., Yang, X. (Imperial), « A Study of Rock Breakage under Extreme Submerged Confining Pressure: Can DTH Hammer Drilling Deliver? », **ARMA 24-0788** — **confinement jusqu'à 130 MPa ; nomme « SOLIDITY »** | OnePetro · **PDF public sur orchyd.eu** |

## Sources de contexte

* **ORCHYD**, H2020 n° 101006752, coordonné par **ARMINES / Mines ParisTech**.
  Livrables **D6.1** et **D6.4** — modèle continu *Concrete Damaged Plasticity*
  sur granite **Red Bohus**. Catalogue : `orchyd.eu/repository/`.
* **Dumoulin, S. et al.**, « Three-dimensional numerical study of DTH bit-rock
  interaction with HPWJ downhole slotting », *Rock Mechanics Bulletin*,
  DOI 10.1016/j.rockmb.2024.100169 — SINTEF, ARMINES et Imperial, Red Bohus.
* **Le code** : `github.com/ImperialCollegeLondon/solidity-solver-open`,
  **LGPL-3.0**, C, ~17 000 lignes, format `.Y3D`. Dépôt public d'Imperial College
  London. ⚠️ **Ce n'est pas la version qui a produit l'article de 2026** — son
  facteur d'endommagement y est câblé à zéro.
* **Guo, L. (2014)**, thèse de doctorat, Imperial College London, « Development of
  a three-dimensional fracture model for the combined finite-discrete element
  method » — dépôt **Spiral**.

## À récupérer

| pièce | pourquoi |
|---|---|
| Guo, Xiang, Latham & Izzuddin (2016), *Eng. Fract. Mech.* **151**, 70-91 | la seule étude de sensibilité au maillage du modèle |
| livrables **WP6** d'ORCHYD (dont D6.2) | un livrable contient les decks qu'un article comprime |
| Munjiza, Andrews & White (1999), *IJNME* **44**(1), 41-57 | l'origine des constantes a = 0,63, b = 1,8, c = 6,0 de la courbe d'adoucissement |
