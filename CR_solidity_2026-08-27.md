# Compte rendu — ce qu'on a fait avec Solidity

2026-08-26 / 27. Branche `joint-handoff`, worktree `rockim_p4`.

---

## Le point de départ

On cherchait à reproduire l'impact à insert unique de Yang, Xiang, Naderi,
Wang, Aising, Ugarte et Latham (IJRMMS 191, 2025) avec rockim. Plusieurs
conventions de leur modèle n'étaient pas publiées : instant de gel du DIF,
filtrage du taux de déformation, régularisation tangentielle du contact,
traitement d'un contact né en recouvrement. Elles étaient jusqu'ici
**déduites** — de la thèse de Guo (2014), des articles, du papier ARMA — et
listées comme « écarts insurmontables faute de source ».

Une recherche a montré que ce n'était pas vrai : **leur solveur est public.**

`github.com/ImperialCollegeLondon/solidity-solver-open` — LGPL-3.0, C,
17 000 lignes, format `.Y3D`, activement maintenu (dernier push 2026-03-31).
C'est la lignée Munjiza décrite par Guo et Yang. Le dépôt a été cloné et lu.

---

## Ce que la lecture a donné

### Onze points où rockim était déjà juste

Rien à changer : loi de volume néo-hookéenne avec terme visqueux `η·D` ;
viscosité comme propriété matériau ; courbe z de Munjiza (a = 0,63, b = 1,8,
c = 6) ; couplage elliptique des modes ; branche élastique parabolique ;
raideur double en compression ; trois points d'intégration aux milieux
d'arêtes avec le poids A/6 ; Mohr-Coulomb avec coupure en traction ; DIF
appliqué à la résistance **et** à son énergie de rupture ; frottement de
Coulomb à glissement mémorisé ; **aucun amortissement dans l'intégrateur**.

Ce dernier point mérite d'être souligné : leur seule dissipation hors joints
et frottement est la viscosité de volume. Le « Mass Damping Coefficient » de
leur Table 1 est le η de `T += η·D`, et comme rockim écrit `2μD`, cela fixe
`bulkViscosity = 2000` pour leur η = 4000.

### Six conventions à porter

| ce qu'on a lu | clé rockim ajoutée |
|---|---|
| la plage d'adoucissement vaut `max(2·δp, 3Gf/ft)` et la rupture est à `δp + ot` | `jointDeltaC = solidity` |
| le joint ne meurt qu'au-delà d'**un** point d'intégration rompu (deux sur trois) | `jointFailRule = majority` |
| le DIF est une variable **locale** de la boucle élément, reprise à chaque pas | `strainRateDIFArm = continuous` |
| au pas de naissance d'un contact issu d'un joint mort, la **pénalité de la paire** est recalée pour que la force soit continue | `gcBirth = penalty` |
| le taux de déformation entre **brut**, sans lissage | `strainRateFilter = none` |
| `ktss = 2/7 · penalty` | `potTangentFactor = 1.4286` (config seule) |

Deux d'entre elles ont levé des points réputés bloquants. **L'instant de gel
du DIF** n'a pas à être deviné : il n'y a pas de gel du tout. Et **la
naissance d'un contact en recouvrement** est bien traitée chez eux, mais par
un mécanisme entièrement différent de celui de rockim — ils rééchelonnent la
pénalité au lieu de relâcher la force sur une constante de temps.

### Ce qui est désactivé chez eux

Le DIF lui-même (facteur mis à 1,0) et le modèle d'endommagement diffus du
volume. Ce second point recoupe ce que leur article de 2026 dit explicitement
du calcaire et du grès : reproduits **sans** le modèle d'endommagement.

### Ce que la source a confirmé comme absent

`jointResidualMu` — une piste venue de Y-Geo, pas d'eux. Leur variable
équivalente vaut zéro en dur : un joint rompu ne porte aucun cisaillement, et
tout le frottement vient du contact. La clé reste désactivée dans la config de
réplication.

---

## Comment on l'a porté

**Principe VIII, sans exception.** Les six capacités sont des valeurs
nouvelles ; aucun défaut n'a changé. La preuve n'est pas déclarative : un
binaire a été compilé depuis `HEAD` dans un worktree jetable et passé sur la
même suite de vérification. Sur les **29 tests d'origine, zéro différence de
valeur physique**, et les mêmes quatre échecs de plateforme MSVC aux mêmes
décimales.

**Dix nouveaux repères** (six en 2D, quatre en 3D) verrouillent les capacités
ajoutées. La suite fast passe de 29 à 40 tests, 36 au vert.

**Des gardes d'exclusivité partout.** `gcBirth = penalty` exige
`contact = potential` et refuse `gcBirthTau` ; `strainRateFilter = none`
refuse `strainRateTau` ; `jointFailRule = majority` exige
`jointQuadrature = midedge`. Le but est qu'aucune clé ne puisse être active et
muette, ni écrite et ignorée — un piège qui a mordu quatre fois sur ce projet.

**Une refactorisation** au passage : `setJointLengths()` centralise le calcul
de `dnE`/`dnF`/`slipF`, jusqu'ici dupliqué en trois sites qui avaient déjà
divergé une fois.

---

## Ce que ça a donné sur l'impact

Un run de réplication a été lancé le 26 à 17:01 et arrêté le 27 à 01:30, à
**19,3 % (t = 106 µs sur 550)**, sur décision de coût.

### Ce qui marche

Trois vérifications indépendantes du montage, toutes au microseconde ou mieux :

- l'**instant d'impact** du piston, prédit à 18,8 µs par le jeu et la vitesse,
  observé à 19,2 ;
- l'**arrivée de l'onde** à l'insert, prédite à 66,7 µs par la longueur du bit
  et la célérité de l'acier, observée entre 65,9 et 67,8 ;
- la **masse de l'outil**, déduite des données par conservation de la quantité
  de mouvement : 1,2182 ± 0,0004 kg, soit le bit seul à 0,03 % près.

Et **σ_zz max = 205,7 MPa** à la jauge, dans leur fourchette 200–260. Premier
des sept critères : validé.

### Ce qui ne marche pas

**321 joints rompus, et sur la dernière frame dépouillée : 85 en traction,
zéro en cisaillement.** Même signature que le run de référence adaptatif
(841 / 3). Les six conventions portées corrigent la loi de volume, la loi de
joint, la quadrature, le DIF, la naissance du contact et le rapport
tangent/normal — et **ne font pas apparaître le broyage**.

La chaîne causale du déficit de frottement est alors complète : pas de
cisaillement → pas de surfaces broyées → pas de surfaces qui glissent → pas de
frottement. Chez eux ce poste porte 32 J sur 49,3 ; ici 0,69 J contre 14,4 J
de fissuration.

### La réserve, qui compte autant que le résultat

On s'est arrêté à **0,305 mm d'enfoncement sur les 1,53 attendus**. La pression
sous l'insert n'a jamais atteint le régime où le broyage se produirait, et
c'est cette pression qui active le cisaillement dans un critère de
Mohr-Coulomb. **Ce n'est donc pas une réfutation, c'est un run interrompu.**

### Le coût, qui a décidé de l'arrêt

Temps par frame : 12,5 → 26 → 43,5 → 56 → 84 → 110 → 103 min, plafonnant vers
105. Soit ~58 h pour aller au bout, dont ~20 h auraient suffi pour atteindre
le pic d'enfoncement (vers 50 %) et quatre critères de plus.

Reste une anomalie non expliquée : **128 ms par pas pour 43 000 tétraèdres sur
14 threads**. Suspicion — en schéma intrinsèque les 81 797 joints existent
tous dès t = 0, et la boucle de contact interroge une table de hachage pour
chaque paire candidate à chaque pas uniquement pour découvrir qu'un joint
vivant la porte. À instrumenter : l'optimisation serait sans effet sur la
physique.

---

## Ce qu'on a appris en se trompant

Cinq erreurs, toutes détectées par un contrôle et non par relecture. Elles sont
consignées parce que chacune produisait un chiffre plausible.

1. **`imp_lib.broken()` comptait 206 facettes au lieu de 88.** Le critère
   `damage >= 0.999` a cessé d'être valide dès l'ajout de `jointFailRule =
   majority` : `J.D` est le max sur les points, donc un seul point rompu suffit
   à mettre le champ à 1. Corrigé en `tBreak >= 0`. Ce bug gonflait la fissure
   radiale et le rayon de cratère d'un facteur 2,3.
2. **Le déplacement du bit n'est pas la pénétration de l'insert.** Le bit se
   fait comprimer par le piston à son sommet pendant que l'insert, 24 cm plus
   bas, n'a pas encore bougé : 0,102 mm lu là où il y en avait 0,004.
3. **Une enveloppe de fissuration « en étoile » entièrement artefactuelle.**
   Deux causes cumulées : un centre pris sur sept facettes et faux de 0,67 mm,
   et un rayon *maximum* par secteur corrélé à +0,79 au nombre de facettes du
   secteur. Recentré et passé au 90ᵉ percentile : Rayleigh p = 0,83, isotrope.
   J'étais à un pas d'annoncer des fissures radiales naissantes.
4. **La force-pénétration, trois fois de suite.** Signe inversé (−65 kN) ;
   masse bit+insert appliquée à l'accélération du bit seul, alors que l'insert
   à 9·10⁵ m/s² pèse 58 kN à lui seul ; et l'idée fausse que la jauge et la
   reconstruction par Newton devaient coïncider, alors que le transit de l'onde
   dure 48 µs pour un événement de 100.
5. **Une bannière jamais imprimée**, placée avant la lecture de sa clé.
   Troisième occurrence de ce piège d'ordre sur ce projet.

La leçon qui revient : sur ce genre de dépouillement, **le contrôle croisé
attrape ce que la relecture laisse passer**. Instant prédit contre instant
observé, masse déduite contre masse déclarée, test statistique contre lecture
à l'œil.

---

## Écarts qui restent avec leur code

**Physique : aucun connu et non porté.**

**Discrétisation** : maillage de 42 882 éléments contre leurs 230 788 ; et
leurs tétraèdres de volume sont des T10 quadratiques avec F-bar là où rockim
est en T4 — leur élément *joint*, lui, est bien le TET4JOINT à six nœuds que
rockim reproduit.

---

## Suite proposée

1. **`jointFrictionScaled = 1`.** Le solveur avertit lui-même que la
   combinaison actuelle (`jointShearUnload = origin` avec
   `jointFrictionScaled = 0`) rend le glissement de joint **réversible, sans
   hystérésis** : le joint ne dissipe rien en cisaillement, par construction.
   C'est la forme littérale de leur article, c'est une clé de config, et ça
   vise exactement le déficit.
2. **Instrumenter les 128 ms/pas** avant toute nouvelle réplique fidèle.
3. **Reprendre le run jusqu'à ~50 %** si l'on veut trancher le mode de rupture
   à pleine pression.

---

## Références

- code : `github.com/ImperialCollegeLondon/solidity-solver-open` (LGPL-3.0)
- site : `solidityproject.com`
- détail technique et repères : `DOCUMENTATION_rockim.md` § 5.4 quinquies
- état du run et pièges de dépouillement : `BILAN_replique_solidity_2026-08-27.md`
