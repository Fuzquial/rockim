# FICHE — Yang X., Xiang J., Naderi S., Wang Y., Aising J., Ugarte I., Latham J.-P.
# « High-fidelity modelling of fragmentation and pulverisation in hard granite
# under percussion loading: a FDEM-based approach », IJRMMS 206 (2026) 106660.
# (Titre et ordre des auteurs verifies sur l'article le 2026-08-29, lot 2b.)
# Lue le 2026-08-29 sur le PDF fourni par F. Uzquiano. Source PRIMAIRE.

## 1. La question qui bloquait : delta_m est-il une longueur ou une deformation ?

**C'est une LONGUEUR.** Texte p. 4, verbatim : « To calculate the damage
factor D, the displacement components are used. Specifically, delta_m
represents the current effective DISPLACEMENT component, delta_m^max
represents the effective displacement max value, and delta_m^0 represents the
effective displacement min value. »

Leurs equations (3) et (4) :

    sigma = { sigma_NH,              0 <= eps <= eps_d
            { C_d eps = (1-D) sigma_barre,   eps > eps_d

    D = { 0,                                                  delta_m <= delta_m^0
        { delta_m^f (delta_m^max - delta_m^0)
          / [delta_m^max (delta_m^f - delta_m^0)],   delta_m^0 < delta_m <= delta_m^f
        { D_max,                                              delta_m > delta_m^f

CONSEQUENCE : `dm = hEl * sqrt(2/3) * ||dev eps||` de rockim (Fdem3dSolver.cpp
l. 2274) est DIMENSIONNELLEMENT CORRECT. L'hypothese d'une erreur d'unite,
formulee le 29/08 a partir du deck d'exemple BST.Y3D de Solidity
(D1PEM0 = 5,0e-3), est REFUTEE : ce deck est un cas SANS RAPPORT avec l'essai
d'impact. La dependance au maillage est reelle mais ASSUMEE par les auteurs,
qui la discutent explicitement (p. 4 : l'element minimal de 1 mm est bien plus
gros que la poudre mesuree a 0,035 mm ; descendre a la taille de poudre
multiplierait le nombre d'elements par plus de 23 000).

## 2. Table 1 — proprietes du granite Kuru Grey, du carbure et de l'acier

| parametre | granite | carbure | acier |
|---|---|---|---|
| masse volumique [kg/m3] | 2626 | 15250 | 7850 |
| module d'Young [GPa] | 60 | 600 | 200 |
| Poisson | 0,24 | 0,2 | 0,29 |
| G_I [J/m2] | 50 | | |
| G_II [J/m2] | 1000 | | |
| resistance a la traction [MPa] | 10,98 | | |
| cohesion [MPa] | 29,84 | | |
| coefficient de frottement interne | 1,85 | | |
| **coefficient de frottement glissant** | **0,18** | **0,6** | **0,6** |
| delta_m^0 (effective displacement min) | 0,014 | | |
| delta_m^f (effective displacement max) | 0,4 | | |
| D_max | 0,9 | | |

Les unites de delta_m ne sont PAS imprimees dans la table. Elles sont
necessairement des MILLIMETRES : 0,014 mm = 1,4e-5 m et 0,4 mm = 4,0e-4 m,
soit exactement les valeurs du deck de rockim. La transcription du depot est
donc JUSTE.

## 3. LE POINT LE PLUS IMPORTANT — 0,18 est le frottement de la ROCHE

> **AVERTISSEMENT AJOUTE LE 2026-08-29 (LOT 2a) — VALABLE POUR LE GRANITE SEUL.**
> Le 0,18 est le frottement du granite Kuru Grey, et de lui seul. Sur le calcaire
> St Anne, Imperial publie **0,6** — deux fois, independamment : IJRMMS 191 (2025)
> Table 4 p. 6, et ARMA 24-0952 Table 1 p. 4. L'acier et le carbure y sont aussi a
> 0,6. Le dernier paragraphe de ce §3, qui conclut que le deck de rockim est
> « 3,3 fois trop eleve », NE VAUT PAS pour St Anne : sur la roche de la these,
> `contactMu = 0.6` est exactement la valeur d'Imperial. Voir
> [la fiche du lot 2a](2026-08-29_lot2a_parametres_stanne.md) §2.

Table 1 range le coefficient de frottement glissant comme une propriete
MATERIAU : 0,18 pour le granite, 0,6 pour le carbure et 0,6 pour l'acier.
Le texte precise (p. 6) : « The sliding friction coefficient is treated as an
effective POST-FRACTURE contact parameter controlling the mobility and
ejection of generated fragments » ; et p. 12 : « the sliding friction
coefficient between joint elements AFTER FAILURE serves a role analogous to
the residual strength of the failed material ».

Autrement dit : chez eux, tout contact ROCHE/ROCHE glisse a 0,18. Comme la roche intacte
est liee par ses joints, un contact roche/roche n'existe qu'APRES rupture :
0,18 est donc le frottement de TOUT contact roche/roche, pas seulement de la
matiere pulverisee. Le 0,6 est reserve au carbure et a l'acier.

ECART AVEC LE DECK DE ROCKIM : `contactMu = 0.6` y est pose GLOBALEMENT, et
0,18 n'est atteint que via `contactResidualMu` sur les elements a D = Dmax
(115 elements sur 18 123 dans P1). Le frottement roche/roche du deck est donc
3,3 fois trop eleve sur l'immense majorite des contacts.

## 4. LE MECANISME DEPEND DE LA ROCHE — et St Anne n'est pas le granite

p. 9, §4.3, verbatim : « When impacted with a hemispherical insert, rock
failure in **St Anne limestone and Rhune sandstone occurs during the REBOUND
stage of the bit, whereas in granite, it develops during the PENETRATION
stage**. This difference is crucial for verifying whether the bit-rock
interaction for the numerical model accurately reflects the actual bit-rock
interaction mechanisms. »

p. 11 : « below the crater, neither the SHEAR CRACK ZONE typically observed in
Rhune sandstone, nor the DOWNWARD-EXTENDING MEDIAN CRACKS commonly seen in
St Anne limestone, are observed [dans le granite]. »

p. 12 : « For St Anne limestone and Rhune sandstone, ROCK REBOUND provides the
dominant driving force for chipping. In contrast, for granite, the outward
ejection of fragmented rock between the insert and the crater exerts dynamic
forces on the surrounding rock, which drives the continued propagation of side
cracks. »

CONSEQUENCE CAPITALE : le modele de pulverisation ET le frottement 0,18 ont
ete developpes pour une roche ou il n'y a NI zone de cisaillement NI fissures
medianes, et ou la rupture se produit A LA PENETRATION. Le calcaire St Anne
possede les deux et rompt AU REBOND. Appliquer la calibration granite au
calcaire, c'est importer un modele construit pour supprimer precisement les
mecanismes que le calcaire exhibe. L'intuition initiale de F. Uzquiano — « les
fissures radiales interviennent en phase de rebond » — est CONFIRMEE par
l'article, pour SA roche.

## 5. Ce que le modele de pulverisation vise, en leurs termes

> **AVERTISSEMENT AJOUTE LE 2026-08-29 (LOT 1).** Le dernier paragraphe de ce
> §5 attribue a Solidity un couplage lu dans `/home/user/solidity` (Y3Did.c).
> Ce dossier N'EST PAS le code d'Imperial et son facteur d'endommagement est
> cable a zero : la phrase « le code Solidity le realise par ... » est donc
> une INFERENCE d'architecture, pas un fait etabli sur Imperial. Voir
> **[MISE A JOUR DU 2026-08-29, LOT 2b] L'ARTICLE EST DESORMAIS DEPOUILLE SUR CE
> POINT : le mot « penalty » n'apparait PAS UNE SEULE FOIS dans les 14 pages.
> Le seul couplage (1-D) publie porte sur la CONTRAINTE D'ELEMENT (eq. 3). La
> penalite de contact n'est jamais dite degradee. La derniere phrase de ce §5 —
> « rockim n'a que la moitie tangentielle » — est donc SANS OBJET : il n'y a pas
> de moitie normale documentee chez Imperial.**
> [Fiche du lot 2b](2026-08-29_lot2b_couplage_endommagement_contact.md) §3.
> [CORRECTION 2](../BILAN_interference_2026-08-29.md) et
> [la fiche du lot 1](2026-08-29_lot1_bibliographie_imperial.md) §7. Le reste
> de cette fiche porte sur l'ARTICLE (source primaire) et tient sans reserve.

p. 4 : « severe local stiffness degradation, LOSS OF LOAD-BEARING CAPACITY,
and fragment-support reduction beneath the insert » ; « By applying larger
deformations to elements at the impact centre and REDUCING SLIDING FRICTION,
the model reproduces experimentally observed bit-rock interactions and
fragment ejection patterns. »
p. 11 : « the combined effects of local stiffness degradation, loss of fragment
support, fragment ejection, and post-fracture contact interactions beneath the
insert ».

L'effondrement de PORTANCE est donc explicitement vise par le papier. Le code
Solidity le realise par `penalty *= min(1-D_i, 1-D_j)` et `mu = mud*d_fact`
(Y3Did.c l. 995, 1044, 1263-1265) : raideur de contact ET frottement suivent
(1-D) EN CONTINU. rockim n'a que la moitie tangentielle, en tout-ou-rien.

Le modele est assigne aux ELEMENTS tetraedriques et « is NOT intended for
joint elements » (p. 4) — comme dans rockim.

## 6. Reperes chiffres pour la validation (granite, piston)

* Chronologie : « rock pulverisation in granite occurs within approximately
  150 us after impact, followed by the development of side cracks and the
  formation of surface chippings » (p. 11).
* Bilans d'energie du bit (masse 1,509 kg), p. 11 :
  - piston 9 m/s : penetration 5,62 m/s, rebond 4,65 m/s, 23,84 -> 16,31 J,
    perte 7,53 J (31,6 %) ; e = 0,827
  - piston 11 m/s : penetration 6,83 m/s, rebond 2,66 m/s, 35,20 -> 5,34 J,
    perte 29,86 J (84,8 %) ; e = 0,389
  Le REBOND S'EFFONDRE entre 9 et 11 m/s : c'est le point de bascule non
  lineaire que le modele de pulverisation existe pour reproduire.
* Fig. 18 : nombre d'elements pulverises 100 (5 m/s) -> 1030 (13 m/s).
* Fig. 9-10 : profondeur d'indentation 0,6 mm (5 m/s) -> 1,6-2,1 mm (13) ;
  masse de fragments 0,1 -> 2,1 g ; longueur de radiale 5 -> 15 mm ; rayon de
  cratere 4 -> 15 mm.
* Maillage : 1 mm dans un DIAMETRE de 25 mm (donc R 12,5 mm), 2 mm jusqu'a
  50 mm de diametre, 10 mm au bord ; insert a 0,7 mm ; 230 788 elements,
  dt = 2,5e-9 s. Eprouvette cylindre 250 x 150 mm.
* Ils declarent n'avoir calibre AUCUN parametre isolement : « GI, GII, the
  damage-model parameters, and the sliding friction coefficient were NOT
  calibrated to fit a single output variable. Instead, they were calibrated as
  an INTEGRATED PARAMETER SET » (p. 6). Reprendre 0,18 sans reprendre G_I,
  G_II et les seuils de delta_m est donc un demembrement du jeu.
