# État de l'art — le FDEM d'Imperial College et le bilan de rockim
**Dossier constitué les 29-30 août 2026. F. Uzquiano, Mines Paris – PSL.**

---

## Ce que contient ce dossier

Une reconstitution, **sur sources primaires**, de la manière dont le solveur FDEM
d'Imperial College (**Solidity**) traite l'impact et le forage percussif — puis le
bilan de ce qui manque à **rockim** pour en faire autant.

**Neuf sources primaires** ont été dépouillées : la totalité de la littérature
d'impact du groupe (IJRMMS 191 et 206, JRMGE, deux ARMA 2024), sa source de
formulation (le manuscrit de Guo, Xiang, Latham & Izzuddin), le chapitre DEM7 qui
porte la loi de frottement, la publication d'origine de celle-ci
(*Engineering Computations* 26(6), 2009), et deux livrables ORCHYD.

**Les PDF ne sont pas dans ce dossier.** Ce sont des articles sous droits
(Elsevier, Springer, Wiley, EAGE) : les redistribuer serait illicite.
[`SOURCES.md`](SOURCES.md) donne la référence complète et le DOI de chacun, plus
son statut d'accès.

---

## Dans quel ordre lire

| # | fichier | pour qui |
|---|---|---|
| 1 | [`SYNTHESE_etat_de_l_art_2026-08-29.md`](SYNTHESE_etat_de_l_art_2026-08-29.md) | **commencer ici.** Présente les cinq lots, distingue *établi / inféré / ouvert*, et se termine par un plan de travail |
| 2 | [`chantier/CONTRE_AUDIT_corrections.md`](chantier/CONTRE_AUDIT_corrections.md) | **lire juste après.** Un contre-audit adversarial a noté le bilan à **52/129**. Ce document porte les corrections et **le plan applicable** |
| 3 | [`fiches/`](fiches/) | le détail, lot par lot, chaque affirmation avec sa page et son équation |
| 4 | [`chantier/`](chantier/) | les changements apportés au code, un document par action |
| 5 | [`MISSION_etat_de_l_art_2026-08-29.md`](MISSION_etat_de_l_art_2026-08-29.md) | le cahier des charges d'origine — **sa section 0 contient une prémisse fausse**, voir ci-dessous |
| 6 | [`BILAN_interference_2026-08-29.md`](BILAN_interference_2026-08-29.md) | le bilan antérieur auquel les fiches renvoient — il porte les deux CORRECTIONS dont la seconde a produit la prémisse fausse |

---

## Avertissement de lecture — deux réserves qui comptent

**1. Ce dossier se corrige lui-même, et il faut lire les corrections.**
Deux passes de vérification adversariale ont été conduites : un second lecteur,
ignorant du premier, retournait à la page ou à la ligne de code et cherchait la
surinterprétation.

| ce qui a été vérifié | verdicts | confirmés |
|---|---|---|
| le dépouillement des **sources** | 673 | **614 — 91 %** |
| le **bilan de rockim** (lot 4) | 129 | **52 — 40 %** |
| dont son **plan de travail** | 17 | **2 — 12 %** |

Les documents d'origine sont conservés avec un avertissement en tête ; **en cas de
contradiction, `chantier/CONTRE_AUDIT_corrections.md` fait foi.** Le motif est
net et vaut d'être retenu : ce qui est passé par une relecture adverse tient, ce
qui n'y est pas passé ne tient pas, et la pièce la plus fausse était celle qui ne
s'appuyait sur aucune source — le plan.

**2. La section 0 du cahier des charges est fausse.**
Elle affirme que le dossier `solidity` est « une transposition partielle faite par
l'équipe ». C'est **le dépôt public d'Imperial College London**
(`ImperialCollegeLondon/solidity-solver-open`, **LGPL-3.0**, C, 17 000 lignes).
Ce qui reste vrai, et qui est l'avertissement utile : **ce n'est pas la version
qui a produit l'article de 2026** — son facteur d'endommagement y est câblé à
zéro. Détail dans [`chantier/A03_resourcer_attributions.md`](chantier/A03_resourcer_attributions.md).

---

## Ce que le dossier établit — en dix lignes

* **La formulation d'Imperial est reconstituée**, équation par équation, avec page
  et numéro : loi de joint cohésive, DIF, contact et potentiel de Munjiza,
  frottement tangentiel, modèle de pulvérisation, retrait des fragments, bilan
  d'énergie, insertion et maillage.
* **Trois choses ne sont publiées nulle part**, et c'est un résultat, pas une
  lacune de recherche : la raideur tangentielle **k_t**, la dissipation visqueuse
  **η**, et la **règle de combinaison du frottement pour deux matériaux
  différents**. Huit sources muettes.
* **Leur partage d'énergie sur St Anne** : 2,6 % à la fissuration, 64,9 % au
  frottement. **Mais ils écrivent eux-mêmes que ce 65 % est gonflé** par leur
  plancher de maillage à 1 mm et par l'angularité de leurs tétraèdres. Le 2,6 %
  est un plancher, pas une mesure.
* **Leur modèle de pulvérisation est superflu sur le calcaire**, de leur propre
  aveu — il a été développé pour le granite.
* **Leur coefficient de frottement glissant n'est pas une propriété physique** :
  c'est pour partie un correctif d'angularité de maillage, pour partie un
  contraste voulu qui amorce les radiales. Le même granite est calibré deux fois
  différemment par la même équipe à un an d'écart.
* **La radiale est le seul des sept critères qu'Imperial rate**, et la cause est
  structurelle : leurs fissures suivent les faces d'éléments, donc en zigzag,
  donc plus courtes que la réalité.
* **rockim a une loi de joint complète et fidèle**, et un DIF conforme jusqu'aux
  points d'application.
* **Un résultat du dépôt mérite publication** : l'exposant 0,1707 de la loi de
  vitesse de déformation, dérivé le 2026-08-18 de la seule figure d'un article de
  2025 qui en imprimait un autre — et **imprimé 0,17 par l'article de 2026**, un
  an plus tard.

---

## Ce qui reste ouvert

Trois pièces à récupérer, aucune bloquante :
*Eng. Fract. Mech.* **151** (2016) 70-91 (la seule étude de sensibilité au
maillage) · les livrables **WP6 d'ORCHYD** · *IJNME* **44**(1) 1999 (l'origine des
constantes de la courbe d'adoucissement).

Et trois questions qui se règlent avec le clone de `solidity-solver-open` sous la
main : a-t-il un commit correspondant aux articles ? le facteur d'endommagement
est-il absent de toute son histoire ? une note de version dit-elle que le code
public diffère de celui des articles ?
