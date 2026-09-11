# L'insertion adaptative est-elle utilisable pour un impact ? — le point du 2026-09-11 (nuit)

**Question de Fernando** : « ma plus-value était l'adaptatif, mais j'ai l'impression qu'il n'est
pas adapté pour un impact. Est-ce possible, et si oui, quelles corrections ? »

**Réponse courte** : oui, un impact tourne en adaptatif — sept bancs cette nuit, tous stables,
bilan fermé, pas de temps ×2,2 et calcul 5 à 15 fois moins cher qu'en intrinsèque. Mais
l'impression est fondée : à physique identique, l'adaptatif forme une **peau de 2 mm** sous
l'insert là où l'intrinsèque forme un **cône broyé de 24 mm** (§5), et ce n'est pas un défaut
corrigé depuis le 22/08 qui l'explique. La cause est mesurée (§7) : le critère d'insertion lit la
**moyenne** des deux tétraèdres et dilue de moitié l'élément qui porte l'anneau de traction
hertzien. **La correction est `facetAverage = max`** (nouvelle valeur, `rockim_g1y7.exe`) : ×15
de joints rompus, cône à 8 mm, premier joint à 79 µs au lieu de 90. Une loi de volume sur la
roche (`lawPhase`, nouvelle clé, avec `dpr` et cap) ne change rien : ce n'est pas une
plasticité qui manquait. Reste un facteur 6 avec l'intrinsèque, qui est le ratcheting diffus des
joints de pénalité, propriété de la discrétisation intrinsèque sur laquelle Yang a calibré.

---

## 1. Ce que le schéma fait dans un impact, mécanisme par mécanisme

| mécanisme | intrinsèque (Yang / Solidity) | adaptatif (Yan 2023, rockim) | différence réelle |
|---|---|---|---|
| roche intacte | joints de pénalité partout : complaisance E/(1 + 1/pf), ratcheting diffus (×3,24 de joints endommagés mesuré le 30/08) | continuum EF exact, nœuds liés par groupes | l'adaptatif est **plus fidèle** à la roche intacte |
| amorçage d'une fissure | le joint quitte sa branche élastique quand SA traction atteint l'enveloppe | la facette s'insère quand la contrainte moyenne des deux tétraèdres atteint l'enveloppe, naissance **au pic** avec continuité de traction (dn0 = σ/pj) | même seuil ; l'intrinsèque endommage par point d'intégration (majority), l'adaptatif insère la facette entière |
| zone broyée sous l'insert | ruine en cisaillement et en traction de milliers de joints à contrainte **modérée** (von Mises médian 100 MPa, σ₁ médian 43 MPa à 80 µs), cascade, cône de 24 mm ; bulkDamage par-dessus | facettes insérées seulement là où la **moyenne** des deux tétraèdres franchit l'enveloppe : peau de 2 mm, continuum élastique dessous (1,36 GPa sur des tets isolés) | **la différence principale** (§5, §7) — réduite ×15 par `facetAverage = max` |
| radiales, latérales (traction) | rupture des joints en mode I, z-curve | insertion en traction puis même z-curve | même loi après insertion |
| pas de temps | ressorts de tous les joints à pf = 20 | facettes liées budgétées à la pénalité des joints INSÉRÉS (pf = 4 chez Yan) | **×2,4 mesuré sur s = 1** (2,43 ns contre 1,00) |
| gros blocs / chips | joints disponibles partout : une fissure peut suivre n'importe quelle facette adoucie | la fissure ne progresse que là où le critère tire : « n'isole pas de gros blocs » (étude du 25/08), remède `insertionTipFactor = 1,6` | à mesurer sur le cratère |

**Ce qui n'était PAS le schéma.** Le run du 22/08 (St Anne, `impact_stanne_fidele_s15`, rebond
0,99, 1 754 joints) tournait **sans** `bulkDamage` (livré la nuit suivante), **sans**
`jointShearRange = coulomb` (le ×54 sur la plage de mode II, levé le 27/08) et à `contactMu = 0,6`
(3,3× trop haut, bilan du 03/09). Le bilan du 22/08 le disait déjà : « ce qui manque est une
dissipation progressive de la zone broyée — précisément ce que le modèle de pulvérisation
fournit ». L'explosion 3D du 07/08 a touché **les deux schémas** (CFL sur le diamètre inscrit,
corrigée). L'explosion du 11/09 matin était celle de la loi de Camacho (TSL rigide), pas du
schéma d'insertion.

## 2. Ce qui interdit ou fausse un impact adaptatif, et comment on le contourne

1. **`strainRateDIFArm = envelope` est refusé en adaptatif** (il n'y a pas d'enveloppe à
   franchir sans joint) → `insertion` (gel à l'insertion) ou `continuous` (Solidity). Le deck
   pose `continuous`, comme l'intrinsèque : même DIF.
2. **`adaptive + jointDeath = damage + gcBirth = penalty` est incohérent** (HANDOFF contact du
   02/09 §8.3 : facteur de naissance clampé à 0,01 à vie). On retire `gcBirth = penalty` (la rampe
   par défaut s'applique) ou on passe en `jointDeath = separation`.
3. **`jointElastic = parabolic` + adaptatif** : la continuité de traction à l'insertion est
   écrite pour la branche linéaire (dn0 = σ/pj) ; sur la parabole de Guo la traction à dn0 vaut
   2σ − σ²/ft, soit un saut de +25 % de ft à σ = ft/2 (insertion pilotée par le cisaillement).
   Exact au pic (σ = ft), faux en dessous. Retiré du deck adaptatif ; **corrigé dans le code**
   (`rockim_g1y4`, `activateJoint` : r = 1 − √(1 − σ/ft), linéaire inchangé) — non encore
   mesuré sur le banc.
4. **`groupContinuum`** (acier et carbure sans joints) fonctionne en adaptatif : les facettes
   permanentes sont exclues du balayage d'insertion (`Joint::perm`, 11/09 soir).
5. **`insertionPenaltyFactor`** : 4 (Yan) donne le gain sur dt ; 20 rendrait le joint inséré
   identique à l'intrinsèque et annulerait ce gain (les facettes liées sont budgétées à cette
   pénalité, puisqu'elles peuvent s'insérer à tout pas).

## 3. Historique utile (pour ne pas le redécouvrir)

| date | fait | source |
|---|---|---|
| 06/08 | adaptatif validé sur l'UCS de Yan : UCS 47,8 MPa, E 99,1 %, dt ×2, bande 65° | [[rockim-insertion-adaptative]] |
| 07/08 | explosion 3D en phase débris, **les deux schémas** ; CFL sur le diamètre inscrit réel | [[rockim-instabilite-3d-debris]] |
| 22/08 | premier impact adaptatif (St Anne) : rebond 0,99, 1 754 joints, pas d'étoile radiale ; diagnostic = dissipation de la zone broyée manquante | `bench_impact/BILAN_fissures_radiales.md` |
| 25/08 | le critère « n'insère pas trop tôt » (98,8 % des insertions servent) ; il n'isole pas de gros blocs ; `insertionTipFactor = 1,6` | `BILAN_insertion_adaptative.md` |
| 27/08 | plage de mode II ×54 (`jointShearRange = coulomb`), six conventions Solidity, replica **intrinsèque** | `rockim_p4/BILAN_replique_solidity` |
| 30/08 | decks Kuru : adaptatif retenu « l'unique écart revendiqué » parce que l'intrinsèque multi-corps « se désassemble » (résolu le 11/09 : ft = 1e12 + groupContinuum) | `bench_impact/configs/impact_kuru9*.cfg` |
| 02/09 | triplet interdit adaptive + damage + gcBirth = penalty | HANDOFF contact §8.3 |
| 03/09 | hypothèse « l'écart vient du schéma » posée, banc préparé, **jamais lancé** | `BILAN_impact3d_yang_2026-09-03` §9 |
| 11/09 | deck v2 intrinsèque, contact parallèle, et cette comparaison | ce document |

## 4. Les decks

- `configs/yang2026_bench_s25_adaptive.cfg` : le banc s = 2,5 en adaptatif — quatre lignes
  changent par rapport à `yang2026_bench_s25.cfg` (`insertion = adaptive`,
  `insertionPenaltyFactor = 4`, `jointElastic = parabolic` retiré, `gcBirth = penalty` retiré).
- `configs/yang2026_impact_adaptive.cfg` : la même variante à l'échelle 1 (dt 2,43 ns mesuré).
- `configs/yang2026_bench_s25_sep.cfg`, `_adaptive_sep.cfg` : variantes `jointDeath =
  separation` + `jointResidualMu = 0.18` (un joint rompu comprimé reste porté par la loi de
  joint avec le frottement de la roche, au lieu de devenir une paire de contact) — c'est le
  levier de **coût** en phase de fracture, voir §5.

## 5. Mesures (banc s = 2,5, 10 563 tets, 14 fils, `rockim_g1y3.exe`)

Lecture à **100 µs** (`python tools/bench_compare.py --t 100e-6 …`), même deck de base, même
piston à 9 m/s (42,8 J) :

| | intrinsèque, `damage` (référence v2) | intrinsèque, `separation` + µ_res 0,18 | **adaptatif** (`damage`, rampe) | adaptatif, `separation` | adaptatif + `jointFrictionMobilised` |
|---|---|---|---|---|---|
| dt (ns) | 2,83 | 2,83 | **6,32** | 6,32 | 6,32 |
| 110 µs en | ~50 min (tué à 104) | 952 s | **209 s** | 183 s | 227 s |
| σ_zz min à mi-bit (MPa) | 173 | 173 | 174 | 174 | 174 |
| premier joint rompu (µs) | 64 | 64 | **90** | 90 | 64 |
| joints rompus | **6 090** | 6 760 | **27** | 25 | 79 |
| facettes insérées (110 µs) | (tous existent) | — | 728 (3,7 %) | 677 | 607 |
| pénétration (mm) | 0,299 | 0,297 | 0,284 | 0,284 | 0,294 |
| v_z bit (m/s) | −6,02 | −5,95 | −5,68 | −5,68 | −5,88 |
| éléments (J, history) | −27,3 | −89,3 | −7,7 | −7,7 | −5,3 |
| contact (J, history) | −3,5 | **+137 (pompe)** | −0,9 | −0,9 | −0,4 |
| bilan final (KE) | — | **31 → 255 J, création** | 31 → 25 J, sain | 31 → 25 J, sain | 31 → 27 J, sain |

**Où sont les joints rompus** (dernier VTU, roche seule) :

| | n | r médian / p90 / max (mm) | profondeur médiane / p90 / max (mm) |
|---|---|---|---|
| intrinsèque `damage` à 80 µs | 2 275 | 8,7 / 14,4 / 26 | 3,3 / 10,1 / 23,6 — **un cône broyé** |
| adaptatif à 110 µs | 89 | 4,9 / 8,6 / 10 | 0,7 / 1,7 / 2,8 — **une peau** |
| intrinsèque `separation` à 110 µs | 7 021 | 20,7 / 39 / 67 | 9,3 / 30 / 61 — la pompe |

**L'état des éléments sous l'insert** (zone r < 15 mm, profondeur < 15 mm, ~1 250 tets) :
intrinsèque à 80 µs, von Mises p50 100 MPa, max **347** ; adaptatif à 110 µs, p50 50 MPa,
max **1 359 MPa**, 27 éléments au-dessus de 500 MPa — et leurs facettes sont **toutes rompues** :
ce sont des tétraèdres isolés entre l'insert et la roche, élastiques, qui portent le GPa sans
qu'aucune loi ne les fasse céder. Le critère d'insertion, lui, fait son travail (rien de lié
autour de ces éléments) ; il ne peut pas produire un cône broyé dans un continuum élastique.

**Deux résultats annexes.** (1) `jointFrictionMobilised = damage` (insertion dès |τ| ≥ c sous
compression) ne change pas le compte (607 facettes) : ce n'est pas le terme frottant qui bride
l'adaptatif. (2) La variante intrinsèque `jointDeath = separation` + `gcBirth = penalty`
**crée de l'énergie** (contact −203 J au bilan, KE ×8, 7 000 joints rompus jusqu'à 67 mm de
l'axe) et le garde-fou `budgetAbortPct` ne l'a pas vue : le poste « intégration » (+68,9 J)
absorbe la création dans l'identité comptable. Une borne physique (KE ≤ KE₀ + travail des
sources) a été ajoutée au garde-fou (`rockim_g1y5`) ; la variante est disqualifiée en l'état.

## 6. Verdict et suite

1. **Oui, l'adaptatif tient un impact** : stable, bilan fermé, dt ×2,2, calcul ×5 à ×15 moins
   cher que l'intrinsèque en phase de fracture, aucune clé qui bloque une fois les quatre lignes
   du §4 posées.
2. **Mais il ne reproduit pas le cône broyé de Yang à paramètres égaux.** L'intrinsèque le
   forme par la ruine en cisaillement de milliers de joints comprimés (83-93 % de ruptures en
   cisaillement) — une propriété de la *discrétisation* à joints de pénalité, sur laquelle la
   calibration de Yang (G_I, G_II, δ₀, δ_f) repose. L'adaptatif garde un continuum élastique sous
   l'insert. C'est exactement le diagnostic du 22/08, que les correctifs du 27/08 et du 03/09
   n'ont pas déplacé : l'écart est **le schéma**, et il est **physique** (un continuum sans loi
   de volume ne s'écrase pas).
3. **La correction n'est pas une clé de joint mais une loi de volume sur la roche** — celle de
   votre note (matrice viscoplastique + ω_c + joints extrinsèques). Ce soir : `lawPhase = rock`
   (nouvelle clé, la loi ne porte que sur la roche, le train de frappe reste élastique) et le
   banc adaptatif + `law = dpr` (cône DP en compression, c et φ de la Table 1, UCS_MC = 236 MPa
   = Kuru, Rankine du volume éteint, traction aux joints) : résultat au §7.
4. **Reste à faire** : les deux schémas à 300 µs et plus (rebond), la calibration de la loi de
   volume sur Kuru (leur UCS 235 MPa, pas de triaxial publié), `insertionTipFactor = 1,6` sur le
   cratère, et la correction de `jointElastic = parabolic` en adaptatif (faite, `rockim_g1y4`,
   non encore mesurée).

## 7. Trois corrections essayées, une qui marche

Même banc, même lecture à 100 µs (`tools/bench_compare.py`) :

| | adaptatif (référence §5) | + `law = dpr` sur la roche (`lawPhase`) | + `dpr` + cap de compaction (`capP0` 235 MPa) | **+ `facetAverage = max`** |
|---|---|---|---|---|
| binaire (ancre bit-identité) | g1y3 (8/8) | g1y6 (8/8) | g1y6 (8/8) | **g1y7 (8/8)** |
| premier joint (µs) | 90 | 100 | 101 | **79** |
| joints rompus à 100 µs | 27 | 2 | 0 | **412** |
| facettes insérées à 110 µs | 728 | 197 | 195 | **2 781 (14 %)** |
| joints rompus à 110 µs | 91 | 2 | 2 | **983** (486 traction / 497 cisaillement) |
| zone rompue (r max / prof. max, mm) | 9,7 / 2,3 | — | — | **14,7 / 8,3** (insérées jusqu'à 22 mm) |
| éléments (J) : prélevé / stocké | 6,7 / 3,2 | 9,2 / 7,7 | 8,8 / 6,7 | 6,6 / 3,3 |
| pénétration (mm), v_z bit (m/s) | 0,284, −5,68 | 0,276, −5,55 | 0,277, −5,58 | 0,291, −5,79 |
| 110 µs en | 209 s | 185 s | 187 s | 243 s |

**La loi de volume ne fait rien** : un cône de Drucker-Prager à φ = 61,6° ne cède pas sous
l'insert (2 J de travail plastique), le cap de compaction à 235 MPa non plus — la zone broyée de
l'intrinsèque ne naît pas d'une plasticité à haute pression mais de la **ruine de joints à
contrainte modérée** (éléments à von Mises médian 100 MPa, σ₁ médian 43 MPa) : ce sont les
facettes bien orientées qui cassent, en traction ou en cisaillement à faible σ_n, et la cascade
fait le reste.

**Ce qui bride l'adaptatif est la moyenne des deux tétraèdres dans le critère d'insertion.**
Sous l'insert, l'élément qui porte l'anneau de traction hertzien a pour voisin un élément
comprimé : la moyenne divise la traction de facette par deux et le critère ne tire pas.
`facetAverage = max` (nouvelle valeur, opt-in) retient le plus chargé des deux : premier joint
à 79 µs, **×15 de joints rompus, ×4 de facettes insérées**, un cône qui descend à 8 mm (insérées
à 22 mm) au lieu d'une peau de 2 mm, et 50 % de ruptures en cisaillement contre 93 % en
intrinsèque. Encore 6× moins de ruptures que l'intrinsèque à 100 µs — le reste est le
ratcheting diffus des joints de pénalité (×3,24 mesuré le 30/08), qui n'existe pas dans un
continuum. C'est la première fois qu'un réglage rapproche vraiment les deux zones broyées.

**Suite** : porter `facetAverage = max` au 2D (`FdemSolver`) ; les deux schémas à 300 µs sur le
banc ; le critère de Camacho-Ortiz par partition nodale des forces (la traction réellement
transmise par la facette liée) comme forme définitive — `max` en est l'approximation par
excès ; puis la loi de volume de la note, qui reste nécessaire pour la compaction à long terme
mais n'est pas ce qui manquait ici.
