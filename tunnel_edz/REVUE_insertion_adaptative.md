# Pourquoi l'insertion adaptative diffuse : revue de la littérature et mesures

Revue commandée le 2026-08-23 après le constat du tunnel de Hutou Beishan :
à lois, maillage, domaine et conditions aux limites comparables, rockim
rompt **28 000** joints là où Wang et al. (2024, schéma intrinsèque) en
comptent **~10 000**, et produit un **nuage de fissures courtes** au lieu de
macro-fractures découpant des blocs. La même signature était apparue au banc
d'impact (`bench_impact/BILAN_fissures_radiales.md`). La question posée :
**le critère de déclenchement et la manière d'insérer sont-ils en cause ?**

Sources : six volets bibliographiques (généalogie Camacho-Ortiz, pathologies
publiées, école Paulino/maillages, codes FDEM, applications roches) et deux
mesures nouvelles faites sur nos propres runs.

---

## 1. Le critère de rockim est conforme à l'état de l'art — vérifié ligne à ligne

L'audit de `insertionSweep()` (FdemSolver.cpp:1931-2015) confronté à la
littérature donne un verdict net sur trois points.

**L'évaluation de la traction.** rockim moyenne les deux tenseurs élémentaires
voisins et les projette sur l'arête en configuration déformée. Camacho & Ortiz
(1996) reconstruisent la traction par **partition nodale des forces**,
`T = (F⁺ − F⁻)/L` (forme vérifiée dans Falk, Needleman & Rice 2001, éq. 2.9-2.11) ;
Pandolfi & Ortiz (2002) évaluent **par facette**. Notre moyenne à deux éléments
est mécaniquement équivalente à la première et plus simple que la seconde :
**conforme**.

**L'enveloppe.** `sigma_n >= ft` OU `|tau| >= c + tanφ·(terme de frottement)`,
avec `jointShearEnvelope = yan` par défaut — c'est l'**éq. 8 de Yan et al.
(2023)** elle-même, le frottement tombant à zéro dès la traction. (L'option
`yang` implémente l'éq. 1 du papier d'impact ; le deck tunnel utilise bien le
défaut, donc suit Yan.) **Conforme à l'article que nous reproduisons.**

**L'absence de bridage.** rockim insère au même pas *toutes* les arêtes qui
franchissent l'enveloppe, sans quota ni tri par sévérité (tri par indice, pour
le déterminisme multi-thread seulement). C'est **exactement** ce que fait le
schéma 3D canonique : Pandolfi & Ortiz vérifient le critère
`t_eff = √(t_n² + β⁻²t_s²) ≥ σ_c` par facette et « les facettes où le critère
est satisfait sont marquées pour traitement », **sans ordonnancement ni
throttling**, la seule garde étant topologique (`CornerToOpen`). Notre
« insertion en tapis » n'est donc pas une déviation d'implémentation : c'est le
comportement du schéma fondateur dans un champ uniformément au seuil.

Deux écarts mineurs relevés, sans lien avec le motif : notre balayage est fait
**à chaque pas** (Pandolfi & Ortiz testent tous les N pas — mais espacer sans
corriger la continuité temporelle aggrave l'overshoot, cf. §3), et notre
pénalité d'activation vaut `4 E/h`, dans le bas de la fourchette de la
littérature (`10 E/h` en ligne de base).

## 2. Deux mesures qui réorientent le diagnostic

**Mesure A — nous comptons bien la même chose qu'eux.** Le compteur `nBroken`
s'incrémente sur `J.D >= 1.0` (FdemSolver.cpp:3468-3470), c'est-à-dire la
**séparation complète**, l'observable que comptent aussi les codes
intrinsèques. Le facteur 3 sur le nombre de fissures est donc réel, et non un
artefact de comptage — objection méthodologique de la revue levée.

**Mesure B — les insertions ne sont pas avortées.** Distribution de
l'endommagement des joints **insérés** (dernière trame) :

| | référence homogène (t = 0,69 s) | Weibull m = 6 (t = 0,96 s) |
|---|---|---|
| insérés | 36 629 | 40 875 |
| D ≥ 0,999 (rompus) | **27 971 — 76,4 %** | **31 022 — 75,9 %** |
| D < 0,05 (quasi intacts) | 3 958 — 10,8 % | 4 658 — 11,4 % |
| D médian | 1,000 | 1,000 |

Trois joints insérés sur quatre vont jusqu'à la rupture totale. L'hypothèse
d'un « tapis d'insertions qui ne s'ouvrent jamais » est **réfutée** : rockim
n'insère pas trop tôt, il insère au bon moment — et ce qu'il insère casse.

Le diagnostic se déplace donc : le problème n'est pas *quand* on insère, c'est
que **28 000 fissures réellement ouvertes refusent de coalescer** en
macro-fractures.

## 3. Ce que la littérature dit du vrai coupable

**Le désordre est l'opérateur de sélection de mode.** Dans un champ lisse à
seuil uniforme, l'ensemble des points au seuil n'est pas un ensemble de points
isolés : c'est une **région entière** qui franchit l'enveloppe dans une fenêtre
de temps étroite. Aucune fissure ne peut alors drainer l'énergie élastique de
son voisinage, donc aucune ne prend le dessus. C'est le constat fondateur de
Tang (1997, RFPA) — l'hétérogénéité de Weibull élément par élément y est la
condition *nécessaire* de la localisation — et la raison pour laquelle
Pandolfi & Ortiz stockent une **limite de traction par facette dès 2002**, puis
Zhou & Molinari (2004) formalisent le **tirage de Weibull par facette** comme
frein physique. Levy & Molinari (2010) montrent que ce sont les **défauts** qui
fixent le nombre de fragments.

**Corollaire inconfortable, et important pour la thèse.** Falk, Needleman &
Rice (2001) établissent que les lois à raideur initiale finie (donc
intrinsèques) **modifient les propriétés élastiques du corps** et que
l'espacement des surfaces cohésives y devient une longueur caractéristique de
la simulation. Kubair & Geubelle (2003) ajoutent que l'extrinsèque est
*numériquement plus stable* que l'intrinsèque. Autrement dit : la belle
localisation en blocs du schéma intrinsèque est **en partie produite par la
complaisance artificielle de ses interfaces** — des bandes molles préexistantes
qui concentrent la déformation. Le schéma de référence n'est pas « plus juste »,
il est **régularisé par un artefact**. Notre objectif ne doit donc pas être
d'imiter l'intrinsèque, mais d'introduire **de façon contrôlée et calibrée**
l'hétérogénéité que l'intrinsèque introduit par accident.

**Le paramètre qui gouverne la LONGUEUR des fissures n'est pas la variance,
c'est la corrélation.** Une pointe de fissure n'avance que si l'arête suivante
est faible *aussi*. Avec un tirage indépendant par arête — un bruit blanc dont
la longueur de corrélation vaut la taille de maille — la probabilité que la
voisine soit sous la moyenne est ~1/2 à chaque pas : la fissure s'arrête après
quelques arêtes. **C'est très exactement le nuage de fissures courtes que nous
observons.** Les travaux récents sur champs aléatoires corrélés en mécanique
des roches (RRFPA 2025 : ACF exponentielle-cosinus, θ = 5 mm ; champ de module
calibré par nanoindentation sur granite : θ = 0,57 mm, CoV 0,18) donnent la
parade : un champ **corrélé spatialement**, avec θ de l'ordre de 1 à 2 tailles
de grain, résolu par `dx ≤ θ/3`.

À noter aussi : l'hétérogénéité n'est **pas monotone**. Trop peu → tapis
simultané ; trop → re-diffusion, chaque élément faible cassant pour son propre
compte. Le motif « blocs » vit dans une fenêtre intermédiaire (m ≈ 3-6 en
convention RFPA, où *m grand = homogène*).

## 4. Ce que notre run Weibull confirme — et ce qu'il ne règle pas

Le run `tunnel_ref_s5_lam1_weib.cfg` (m = 6, tirage **indépendant par joint**,
seule clé ajoutée au deck de référence) donne, à temps simulé égal (t = 0,618 s) :
**25 965 joints rompus contre 26 586** pour l'homogène, soit **−2,3 %**. La
distribution d'endommagement est inchangée (75,9 % contre 76,4 % de rompus).

C'est cohérent avec la littérature, pas contradictoire : un tirage indépendant
par arête **désynchronise les instants d'insertion** (ce qu'on voulait) mais ne
crée **aucun défaut de taille physique** capable de guider une macro-fissure.
La variance seule ne suffit pas ; il manque la corrélation.

## 5. Ce qu'il faut faire ensuite, par ordre de rentabilité

1. **Champ de résistance corrélé** — `strengthCorrLength` existe déjà dans
   rockim (champ gaussien échantillonné via copule, indépendant du maillage,
   `fieldSeed` séparé). Balayer θ ≈ 0,5 / 1 / 2 m à l'échelle du tunnel
   (l'équivalent de 1-2 tailles de « grain » mésoscopique, avec dx = 0,21 m on
   résout θ ≥ 0,6 m). **Prédiction falsifiable** : fissures plus longues et
   moins nombreuses, *à nombre total d'insertions décroissant*.
2. **Normalisation par le volume effectif** — Zhou & Molinari en font la
   condition de l'objectivité au maillage (« le concept de volume effectif
   réduit significativement la dépendance au maillage »). rockim la possède
   déjà, en clé séparée et opt-in : `jointSizeEffect` applique
   `(Zeff/V_J)^{1/m}` joint par joint (FdemSolver.cpp:1834-1841, avec garde
   `jointSizeEffectClamp`). Elle n'était **pas activée** dans le run Weibull :
   à activer dans le balayage suivant, sans quoi raffiner multiplie les
   maillons faibles et gonfle mécaniquement le nombre d'insertions.
3. **Contrôle croisé** : jitter géométrique de maillage à seuils homogènes
   (Molinari et al. 2007), pour séparer le bruit géométrique du bruit matériau.
4. **Orientations** des joints rompus contre orientations d'arêtes, pour
   quantifier ce qui reste de biais de maillage.
5. **Contrôle intrinsèque** sur le même maillage (`insertion = intrinsic`),
   pour mesurer chez nous l'écart de motif entre les deux schémas — et vérifier
   la thèse de Falk-Needleman-Rice sur nos propres données.

Un point de calibration à ne pas oublier : l'effet d'échelle Weibull abaisse la
résistance macroscopique (σ ∝ V^{−1/m}) ; introduire de l'hétérogénéité impose
de **recalibrer la moyenne** sous peine de perdre le calage GBM Red Bohus et la
validation Abaqus.

---

### Références principales

Camacho & Ortiz, *IJSS* 33 (1996) 2899 — l'insertion à seuil.
Ortiz & Pandolfi, *IJNME* 44 (1999) 1267 ; Pandolfi & Ortiz, *Eng. Comput.* 18
(2002) 148 — schéma 3D canonique, limite de traction par facette.
Falk, Needleman & Rice, *J. Phys. IV* 11 (2001) Pr5-43 — intrinsèque vs
extrinsèque, l'espacement cohésif comme longueur caractéristique.
Kubair & Geubelle, *IJSS* 40 (2003) 3853 — stabilité comparée.
Papoulia, Sam & Vavasis, *IJNME* 58 (2003) 679 ; Sam, Papoulia & Vavasis,
*EFM* 72 (2005) 2247 — continuité temporelle à l'activation.
Zhou & Molinari, *IJSS* 41 (2004) 6573 — Weibull par facette, volume effectif.
Levy & Molinari, *JMPS* 58 (2010) 12 — les défauts fixent le nombre de fragments.
Tang, *IJRMMS* 34 (1997) 249 — RFPA, hétérogénéité et localisation.
Yan, Zheng & Wang, *IJRMMS* 169 (2023) 105439 — l'insertion adaptative FDEM.
Wang et al., *Front. Earth Sci.* 12 (2024) 1517816 — le cas tunnel de référence.
