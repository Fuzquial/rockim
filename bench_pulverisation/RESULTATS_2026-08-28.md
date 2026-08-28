# Banc pulvérisation — résultats du 2026-08-28

Trois runs de 150 µs sur `box3d_h45.msh` (19 326 tets), sphère R 8,51 mm,
1,28 kg, 9,5 m/s, St Anne, physique complète (coulomb + bulkDamage granite).
Binaire `build_fix2` (WP6 + les cinq réparations + reprises). ~9 min/run.

| | A — complet | B — contrôle | C — A + `jointDeath = damage` |
|---|---|---|---|
| `contactResidualMu` | 0,18 | absent | 0,18 |
| v finale (m/s) | −7,751 | −7,659 | **−6,670** |
| enfoncement (mm) | 1,225 | 1,222 | 1,196 |
| joints rompus | 2 550 | 2 412 | **506** |
| fragments | 522 | 483 | **61** |
| volume détaché (m³) | **1,84e−5** | 1,46e−5 | 1,91e−6 |
| dissipation frottement (J) | **8,10** | 8,90 | 2,26 |
| travail contact général (J) | 6,05 | 6,74 | 2,82 |
| dissipation pulvérisation (J) | 2,76 | 2,93 | **6,29** |
| éléments à D = Dmax | 9 | 9 | **22** |
| évaluations au µ résiduel | 349 071 | — | **596 599** |
| joints morts en compression | 343 (25,6 %) | 368 (30,1 %) | **268 (53,0 %)** |

## A vs B — le mécanisme WP6 agit, mais pas comme je l'avais prédit

**Critère 3 (zéro artefact) : PASS.** Premier contact outil à 24,8 µs,
premier engagement du µ résiduel à 60,8 µs — le mécanisme ne s'engage
jamais avant le contact, ni avant que la pulvérisation ait commencé.

**Signature directe : −9 % de dissipation par frottement** (8,10 contre
8,90 J) et −10 % de travail du contact général. C'est exactement ce que
doit produire un µ qui tombe de 0,6 à 0,18 sur la matière pulvérisée.

**Éjection accrue : +26 % de volume détaché** (1,84e−5 contre 1,46e−5 m³),
+8 % de fragments, +6 % de joints rompus. La mobilité des fragments
augmente : prédiction du deck VÉRIFIÉE.

**Prédiction RÉFUTÉE — la vitesse.** J'avais inscrit « v finale |A| < |B| »
en raisonnant « perte d'appui → l'outil s'enfonce plus → plus d'énergie
dans la roche ». Le run donne l'inverse : A garde PLUS de vitesse
(−7,751 contre −7,659). L'erreur de raisonnement est identifiable :
`contactResidualMu` réduit la résistance **tangentielle**, pas la portance
**normale**. Moins de frottement = moins de dissipation = outil moins
freiné. L'enfoncement n'augmente que de 0,2 % (2,5 µm), marginal.

**Conséquence pour la réplication de Yang** : l'effondrement de portance
qu'ils décrivent ne peut pas venir du seul µ résiduel. Chez eux il naît du
COUPLE dégradation de raideur (D → 0,9) + friction réduite ; la raideur est
déjà dans `bulkDamage` (actif dans les trois runs) et ne suffit visiblement
pas. Le canal manquant est probablement la **raideur de pénalité non
dégradée** que la revue adverse avait signalée comme non traitée (point 4) :
un élément à 10 % de raideur reçoit encore une pénalité de contact pleine.
À instruire avant de conclure sur le rebroussement non monotone.

## C vs A — `jointDeath = damage` ouvre bien le relais, mais change TOUT

**Critère 2 : PASS, et largement.** Les évaluations au µ résiduel passent de
349 071 à 596 599 (**+71 %**), la proportion de joints morts en compression
de 25,6 % à 53,0 %, la charge lâchée au relais de 118 à 196 kN (+66 %). Le
canal roche/roche s'ouvre exactement comme la chaîne causale le prévoyait.

**Mais l'effet dépasse de loin le mécanisme visé** : 5× moins de joints
rompus, 8,5× moins de fragments, 9,6× moins de volume détaché, 3,6× moins
de frottement, 2,3× plus de pulvérisation volumique, et une force de pointe
qui monte à 49 kN contre 34 kN. L'outil décélère beaucoup plus (2,83 m/s
cédés contre 1,75).

Lecture : en tuant les joints dès l'endommagement complet, le milieu devient
discontinu bien plus tôt ; la charge passe au contact, la roche travaille en
assemblage granulaire plutôt qu'en continuum fissurant. Moins de joints
cassent parce que le chemin de charge a changé, et les éléments se
déforment davantage (d'où la pulvérisation accrue).

## DÉCISION sur le défaut `jointDeath`

**Le défaut reste `separation`.** Un changement qui déplace les comptages
d'un facteur 5 à 10 n'est pas un correctif : c'est un choix constitutif
distinct, qui invaliderait silencieusement tout résultat antérieur et pour
lequel nous n'avons AUCUNE donnée expérimentale départageant les deux.
Le levier est implémenté, documenté et désormais chiffré — il s'active
explicitement au deck quand on veut le canal roche/roche.

## Ce que le banc ne dit pas

Fenêtre de 150 µs, donc phase de charge uniquement : ni décharge, ni
rebond, ni écaillage. Maillage 1,9 mm (grossier au regard de ℓ_ch = 14 mm,
mais conforme à h ≤ ℓ/3). Les tendances A/B/C sont robustes ; les valeurs
absolues ne sont pas des prédictions de fidélité.
