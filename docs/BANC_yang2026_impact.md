# Reproduction du banc d'impact de Yang et al. 2026 dans rockim

**Date** : 2026-09-11 · **Version** : `rockim_g1fix.exe` · **Statut** : deck monté et
instrumenté, **aucun run long lancé**, quatre corrections en attente de validation.

---

## 1. Objet

Reproduire **à l'identique** la loi constitutive et le montage d'impact à insert unique
de :

> X. Yang, J. Xiang, S. Naderi, Y. Wang, J. Aising, I. Ugarte, J.-P. Latham,
> *High-fidelity modelling of fragmentation and pulverisation in hard granite under
> percussion loading: a FDEM-based approach*,
> **Int. J. Rock Mech. Min. Sci. 206 (2026) 106660**. Code : **Solidity** (Imperial
> College). Données expérimentales : Aising et al., Mines Paris-PSL / Drillco.

Trois articles du même groupe complètent le dossier et ont été dépouillés :

| réf. | apport |
|---|---|
| Yang et al., **IJRMMS 191 (2025) 106125** | validation **multi-critères** (7 critères), algorithme de retrait des fragments |
| Yang et al., **JRMGE 17 (2025) 6095** | taxonomie des fissures, **inserts multiples**, effet du coefficient de glissement |
| Naderi et al., **JRMGE 17 (2025) 6868** | substitut ANN, géométrie d'insert paramétrée, **seule source sur la pénalité** |

*(Les fichiers `ICL.pdf` et `1-s2.0-S1365160925001029-main.pdf` sont des doublons
respectifs du Yang 2026 et du Yang 2025 IJRMMS.)*

---

## 2. Correspondance loi ↔ clés rockim

Tout vient de leur **Table 1** et de leurs équations, sans substitution.

| grandeur (article) | valeur | clé rockim |
|---|---|---|
| densité, module, Poisson | 2626 kg/m³, 60 GPa, 0,24 | `rho`, `E`, `nu` |
| résistance en traction | 10,98 MPa | `ft` |
| cohésion | 29,84 MPa | `cohesion` |
| coefficient de frottement interne **1,85** | φ = atan(1,85) | `frictionDeg = 61.61` |
| G_I = 50, G_II = 1000 J/m² | rapport **20** | `Gf = 50`, `gfShearFactor = 20` |
| joints six-nœuds **pré-insérés** | intrinsèque | clé `insertion` **absente** (défaut) |
| adoucissement | z-curve Munjiza (a=0,63, b=1,8, c=6) | `jointSoftening = munjiza` |
| Mohr-Coulomb à coupure de traction (éq. 1) | terme frottant plafonné à f_t | `jointShearEnvelope = yang` |
| DIF compression (éq. 1) | 0,77 + 0,56 ε̇^**0,07** | `strainRateDIF = yang`, `difExpS = 0.07` |
| DIF traction (éq. 2) | 0,95 + 0,41 ε̇^**0,17** | `difExpT = 0.17` |
| pulvérisation (éq. 3-4) | δ₀ = 0,014 mm, δ_f = 0,4 mm, D_max = 0,9 | `bulkDamage = yang` + 3 clés |
| glissement 0,6 intact / 0,18 pulvérisé | | `contactMu = 0.6`, `contactResidualMu = 0.18` |
| carbure | 15250 kg/m³, 600 GPa, ν = 0,2 | `phase.insert.*` |
| acier | 7850 kg/m³, 200 GPa, ν = 0,29 | `phase.bit.*`, `phase.piston.*`, … |

### Deux clés imposées par rockim, sans équivalent dans l'article

- **`strainRateDIFArm = envelope`** — en schéma intrinsèque il n'existe aucun instant
  d'insertion où figer le DIF ; le solveur **refuse le deck** sans cette clé. `envelope`
  gèle le DIF au franchissement de l'enveloppe de rupture.
- **`jointDeath = damage`** — avec le défaut `separation`, un joint écroui **en
  compression** sous l'indenteur ne meurt jamais, donc le relais de contact roche/roche
  ne s'engage pas et **`contactResidualMu` n'est jamais atteint** — précisément là où
  Yang le fait agir. Avertissement levé par le banc court.

---

## 3. Deux trouvailles sur les articles eux-mêmes

### 3.1 L'exposant du DIF en traction : 0,07 (2025) contre 0,17 (2026)

L'article de **2025** (IJRMMS 191, éq. 3) imprime `0,95 + 0,41 ε̇^0,07` pour la
traction ; celui de **2026** (éq. 2) imprime `0,95 + 0,41 ε̇^0,17`.

rockim portait déjà une option `strainRateDIF = yang-fig2` à **0,1707**, déduite en
lisant leur figure 2b parce que le 0,07 de 2025 rendait la loi discontinue. **Les
auteurs ont corrigé eux-mêmes en 2026** : la lecture de la figure était la bonne. Le
deck pose désormais l'exposant publié, `difExpT = 0.17`.

### 3.2 La raideur de pénalité n'est pas publiée

Recherche plein texte sur les quatre PDF (`penalt`, `stiffness`) :

| article | pénalité de joint |
|---|---|
| Yang 2025, JRMGE | **aucune mention** |
| Yang 2025, IJRMMS 191 | **aucune mention** |
| Yang 2026, IJRMMS 206 | **aucune mention** |
| Naderi 2025 | seule source — et elle donne **trois valeurs incompatibles** |

Naderi :

1. **Table 1** : `k = 900 GPa`, grès `E = 37 GPa` → **k/E = 24,3**
2. **§2.2** : « *set to ten times the elastic modulus* (Mahabadi 2012 ; Li et al. 2019) » → **10**
3. **§2.3.2** : « *the values of E and C in the dynamic model are increased tenfold* » →
   si E = 370 GPa dans le calcul, alors **k/E = 2,4**

*(La troisième est en soi douteuse : 370 GPa pour un grès dépasse le carbure de
tungstène et multiplie la célérité par √10. Rapportée parce qu'elle est écrite.)*

Grandeur distincte, également chez Naderi : la raideur de **contact** vaut 110 GPa en
normal et tangentiel, soit 3 E.

**Conséquence** : reproduire la pénalité de Yang est impossible — il ne la publie pas.
C'est une incertitude irréductible de plusieurs points sur le module effectif, à
assumer et à écrire.

---

## 4. Le maillage

Produit par l'outil **déjà existant** `tools/make_impact_mesh.py` (spec 005, WP4), qui
reprend leur Fig. 5 (géométrie) et leur Fig. 6 (tailles) :

```
python tools/make_impact_mesh.py meshes/impact_yang_s1.msh 1.0
```

**Six corps** : `rock` (cylindre R 125 × 150 mm), `insert` (carbure, hémisphère
R 8,51 mm, 23,2 mm de haut), `bit` (acier, Ø30 × 242 mm), `piston` (Ø26,5 × 260 mm),
`circlip`, `plate`.

| corps | tétraèdres | part | masse |
|---|---|---|---|
| roche | 98 546 | **82,2 %** | 19,320 kg |
| insert | 9 698 | 8,1 % | 0,065 kg |
| bit | 5 699 | 4,8 % | 1,290 kg |
| **piston** | **1 033** | **0,9 %** | 1,057 kg |
| circlip + plaque | 4 925 | 4,1 % | 0,198 kg |
| **total** | **119 901** | | |

Gradation 1 mm dans la boule R 12,5 mm → 2 mm à R 25 → ~13 mm au bord.
**Qualité** (insphère/circonsphère, 1 = régulier) : médiane 0,88, minimum 0,158,
**aucun élément sous 0,10**. 231 471 joints, 45 832 nœuds.

Outils de visualisation écrits pour l'occasion :

- `tools/fig_mesh3d.py` — coupe **exacte** (les tétraèdres sont réellement découpés par
  le plan, pas approximés par une tranche), zoom, **profil de la taille réalisée** et
  qualité. Le profil est le contrôle falsifiant : un champ de taille peut être juste et
  mal réalisé.
- `tools/fig_montage_impact.py` — assemblage coloré par corps, jeux initiaux,
  chronologie chiffrée en pas de temps et en heures.

Variante **roche seule** à la gradation littérale de leur Fig. 6 :
`meshes/impact3d_yang.msh`, **249 175 tétraèdres** contre leurs 230 788 (écart 8 %),
via `tools/make_impact3d_mesh.py --cylinder --mid` (deux options opt-in ajoutées, le
chemin par défaut reste **bit-identique**, vérifié par MD5).

---

## 5. Mesures

### 5.1 Coût, mesuré et non estimé

Trois runs courts (22, 531 et 2 123 pas) :

| grandeur | valeur |
|---|---|
| pas de temps | **dt = 9,424·10⁻¹⁰ s** |
| coût par pas | **0,125 s** (croît avec la durée : 0,100 sur 531 pas) |
| initialisation | ~9,5 s |

| durée simulée | pas | temps |
|---|---|---|
| 1,2·10⁻⁴ s | 127 000 | 4,4 h |
| 2·10⁻⁴ s | 212 000 | 7,4 h |
| 8·10⁻⁴ s (leurs figures) | 849 000 | **29,4 h** |

rockim est **2,7× plus prudent que Yang** sur dt (il annonce 2,5·10⁻⁹ s).

### 5.2 Chronologie — pourquoi un run court ne montre rien

| événement | temps | pas | calcul |
|---|---|---|---|
| le piston part (jeu 0,20 mm) | 0 | 0 | — |
| il touche le bit | 22 µs | 23 579 | 0,8 h |
| **l'onde sort du bit et charge la roche** | **70 µs** | **74 410** | **2,6 h** |
| pic de force (leur Fig. 9a) | ~100 µs | 106 242 | 3,7 h |
| rebond du bit (leur Fig. 9d) | ~470 µs | 498 842 | 17,3 h |

22 µs pour franchir le jeu à 9 m/s, puis 48 µs pour traverser 242 mm d'acier à
5 048 m/s. **Le run court s'est arrêté à 2 µs** : 0 joint cassé, 0 DIF gelé,
`contactResidualMu` jamais engagé, énergie cohésive 10⁻¹⁸ J. Attendu.

Les **22 µs de vol du piston sont un artefact du mailleur** (jeu laissé pour
l'initialisation du contact général) : 0,8 h de calcul pour zéro physique. Les 48 µs de
traversée du bit, eux, sont physiques.

### 5.3 Banc de pénalité — la pénalité **change la physique**

Barreau 20×20×40 mm, 5 534 tétraèdres, **joints rendus incassables**
(`ft = cohesion = 1e12`) pour isoler la seule complaisance élastique. `E = 60` GPa.
Decks dans `configs_bench/pen_f*.cfg`.

| f | E apparent | vs f = 80 | prédit 1D | dt (ns) | gain dt | coût |
|---|---|---|---|---|---|---|
| 5 | 48,24 GPa | −13,87 % | −15,62 % | 5,256 | ×1,80 | 89 s |
| 10 | 51,49 | −8,06 % | −7,95 % | 3,984 | ×1,36 | 123 s |
| **20** | **53,87** | −3,82 % | −3,57 % | 2,928 | ×1,00 | 232 s |
| 40 | 55,29 | −1,27 % | −1,22 % | 2,113 | ×0,72 | 233 s |
| 80 | 56,01 | 0 | 0 | 1,510 | ×0,52 | 339 s |

En intrinsèque chaque interface porte un ressort `p = f·E/h` ; élément et joint sont en
série :

```
E_eff = E / (1 + 1/f)
```

La formule est vérifiée à **0,1–0,3 point** sur quatre décades. Seul f = 5 décroche.

Trois conséquences :

1. **f = 20 → 10 amollit l'éprouvette de 4,4 %.** Ce n'est pas un réglage neutre.
2. **Le gain sur dt suit f^−0,44, pas f^−0,50** : le CFL du continuum met un plancher.
   20 → 10 ne gagne que **×1,36**.
3. **L'asymptote mesurée est E_∞ = 56,4 GPa, pas 60** — il reste −6 % qui ne viennent
   *pas* de la pénalité (inertie à 1,25 /s, amortissement, mesure aux mors). Seuls les
   **écarts relatifs entre lignes** sont exploitables.

Compensation exacte : `E_entrée = E_cible × (1 + 1/f)` → 63,0 GPa à f = 20.
**Non appliquée** : elle ferait sortir le deck de « exactement leur Table 1 ».

**Décision retenue** : garder `jointPenaltyFactor = 20`, le plus proche de la seule
valeur tabulée publiée (24,3 chez Naderi). À noter : f = 20 et f = 40 coûtent le même
temps sur ce petit banc (232 vs 233 s) — à revérifier sur le maillage d'impact, où la
détection de contact pèse bien plus lourd.

---

## 6. Erreurs à corriger

### 6.1 La jauge de contrainte est au mauvais endroit — **la plus sérieuse**

Le solveur **décale le maillage de +0,149 m** (il ramène le bas de la roche à z = 0).
Non pris en compte à l'écriture du deck.

| | repère maillage | repère solveur |
|---|---|---|
| bit | [23,4 ; 265,2] mm | **[173,4 ; 415,2] mm** |
| milieu du bit | 144,3 mm | **294,3 mm** |
| clé posée `gauge.bit` | — | **[150 ; 200] mm** |

La fenêtre attrape les **26 mm du bas du bit**, juste au-dessus de l'insert. Or Yang
place ses jauges **au milieu** du bit (leur Fig. 9a). La contrainte lue serait celle du
**champ proche du contact**, pas celle de l'onde incidente — et c'est le **critère n° 1
de leurs sept**.

> **Correction : `gauge.bit = 0.28 0.31`**

### 6.2 Un outil fantôme

Le scénario `percussion` crée **par défaut** une sphère de 0,5 kg lancée à 8 m/s
(16 J), en plus des six corps :

```
# toolShape = sphere    (defaut)
# toolMass  = 0.5       (defaut)
# impactSpeed = 8       (defaut)
```

`toolShape = none` a été oublié — le deck de référence `fdem3d_bench1_insert.cfg` le
pose systématiquement dès qu'un corps est lancé par `groupVel`. La sphère vole 15 mm
au-dessus du piston, donc elle ne le rattraperait qu'à ~1,9 ms (hors du run), mais elle
fausse le bilan d'énergie et ajoute une masse non voulue.

> **Correction : `toolShape = none`**

### 6.3 `gravity = 0`

Yang écrit « *gravity, which is also applied to the model* ». Le défaut de rockim est 0.

Effet chiffré : sur 800 µs, g ajoute 7,8·10⁻³ m/s à un impact de 9 m/s, soit **0,09 %**
— physiquement négligeable, mais c'est un écart à « exactement leur loi ». Le solveur
prévient en outre que le travail de la pesanteur n'entre pas dans le bilan sans
`energyBodyForces = on`.

> **Correction : `gravity = 9.81`**

### 6.4 Cosmétique : `phase.<corps>.fraction = 1`

Sans objet avec des groupes physiques nommés. Produit l'avertissement trompeur
« rock 95,7 % (target 16,67 %) ».

> **Correction : retirer les six lignes**

---

## 7. Vérifié et correct — ne rien changer

- **La base de la roche est bloquée automatiquement.** Dès que `absorbing ≠ all`, le
  solveur pose `flag_ = FIXED` sur la face z = 0 (`Fdem3dSolver.cpp`, commentaire
  « *percussion AND shear: the block needs its support* »). Conforme à Yang
  (« *the lower surface of the rock is fixed* »). Aucune clé à ajouter.
- `insertion` absente → **intrinsèque**, comme Yang. Ne **pas** mettre `adaptive` :
  l'insertion adaptative est Yan, Zheng & Wang 2023, un autre schéma.
- `difExpT = 0.17` bien appliqué (confirmé par `config_effective.cfg` ; le
  « yang-fig2 » du journal est une étiquette d'affichage, pas le mode réel).
- `jointDeath = damage`, `jointShearEnvelope = yang`, `strainRateDIFArm = envelope`.

---

## 8. Choix non contraints par l'article

| clé | valeur | justification |
|---|---|---|
| `jointPenaltyFactor` | 20 | Yang ne la publie pas ; 24,3 chez Naderi (Table 1). Sensibilité mesurée §5.3. |
| `gcRestitution` | 0,2 (défaut) | Le solveur avertit qu'il faut la figer dès que le rebond est une métrique — il l'est. À balayer. |
| `absorbing` | none | Yang : base bloquée, flancs **libres**, justifié par le rapport 1:7 entre extension de fissure (21 mm) et diamètre (250 mm). |
| amortissement | défaut | Yang n'en mentionne aucun (impact franc). |
| `dtFactor` | 0,15 | dt résultant 2,7× sous celui de Yang. |

---

## 9. Le modèle est HOMOGÈNE — pas de GBM

Journal du run :

```
voronoi: 6 grains, 6 phase(s), 0 grain-boundary joints
joints intra/homo/hetero: 231471/0/0
```

Les « 6 grains » sont les six **corps**, pas des grains minéraux. **Zéro joint de
frontière de grain.**

C'est **fidèle à Yang** : leur Table 1 n'a qu'une colonne pour le Kuru Grey, et ils
assument leur hétérogénéité autrement — « *the unstructured mesh in the rock also
introduces the rock's inherent heterogeneity into the simulation* » (§3.1).

À noter : maille fine **1 mm** contre granulométrie **0,27–1,5 mm** — le maillage est
à l'échelle du grain, donc GBM implicite, mais **sans contraste de propriétés ni
frontières affaiblies**.

Un vrai GBM demanderait 35 à 90 éléments par grain (règle de la thèse), soit ~0,2 mm
dans la zone d'impact : ×100 d'éléments **et** dt divisé. À garder comme **variante
A/B** après le run homogène de référence, jamais avant.

---

## 10. Ce que rockim n'a pas, face aux quatre articles

Le noyau physique **est déjà dans rockim** (DIF des deux signes aux formules exactes,
`bulkDamage = yang` avec les constantes de leur Table 1 en défaut, `contactMu.<phase>` +
`contactResidualMu`, `gauge.<corps>`, `trackGroups`, métriques de cratère). Manquent :

1. **L'algorithme de retrait des fragments dynamique** (Yang 2025 §2.3) — test
   mécanique par anti-gravité (v = 2,5 mm/s, a = 98,1 m/s², 0,5 ms, β = 0,8).
   `crater_metrics.py` n'a qu'un proxy **géométrique**, qui ne voit pas les fragments
   mécaniquement coincés. Affecte la **masse de fragments**, la métrique la plus
   utilisée. *(Écarté par l'utilisateur comme secondaire.)*
2. **La taxonomie des fissures** (médiane / latérale / radiale / inclinée) — empêche de
   tester leur résultat central : les latérales naissent du déchirement par recouvrance
   élastique, pas de la fermeture des médianes.
3. **Inserts multiples et coalescence de cratères** — `make_impact3d_mesh.py` est
   « insert unique » ; Yang 2025 JRMGE fait trois inserts Ø9 mm à 13,86 mm d'entraxe.
4. **Géométrie d'insert paramétrée** (courbure `c` de Naderi : plat / hémisphérique /
   ogival).
5. **Rotation combinée à la percussion** — trou **de la discipline**, pas de rockim.
6. **Boucle d'optimisation de design + UQ** (Naderi). rockim a en revanche
   `bayes_bench.py` (émulateur GP + MCMC), plus fort qu'un MLP pour la *calibration*.

**Écart méthodologique de fond** : les quatre articles calibrent G_I/G_II sur des
**essais d'impact** contre sept critères ; notre calibration Red Bohus porte sur un
**triaxial**. Les deux chemins ne se croisent pas — le premier contraint l'énergie de
fissuration, le second l'enveloppe de rupture.

---

## 11. Décisions en attente

1. **Appliquer les quatre corrections** du §6 ?
2. **Réduire le jeu piston-bit** de 0,2 à 0,02 mm (0,8 h gagnées, zéro effet physique) ?
3. **Durée** : 1,2·10⁻⁴ s (4,4 h, onde + pic + début de fissuration, **pas** le rebond)
   ou 8·10⁻⁴ s (29,4 h, les sept critères) ?
4. **Vitesse du piston** : un seul point à 9 m/s, ou la série 5,68 / 9 / 11 / 13 ?
5. **Doublon `insertion = adaptive`** pour un A/B, une fois la référence intrinsèque
   acquise ? Ce serait un résultat en soi : mêmes critères pour la moitié du temps.

---

## 12. Fichiers produits

| chemin | rôle |
|---|---|
| `configs/yang2026_impact.cfg` | le deck, entièrement commenté ligne à ligne |
| `meshes/impact_yang_s1.msh` | maillage six corps, 119 901 tets |
| `meshes/impact3d_yang.msh` | roche seule, gradation littérale, 249 175 tets |
| `tools/fig_mesh3d.py` | vue d'un maillage 3D : coupe exacte, zoom, profil, qualité |
| `tools/fig_montage_impact.py` | assemblage, jeux, chronologie chiffrée |
| `tools/make_impact3d_mesh.py` | `--cylinder` et `--mid` ajoutés, défaut bit-identique |
| `configs_bench/pen_f*.cfg` | banc de pénalité, 5 points |
| `results/mesh_yang4corps_vue.png` | figure du maillage |
| `results/montage_impact.png` | figure du montage |
| `docs/BANC_yang2026_impact.md` | ce document |

Runs de mesure seulement (aucun run long) : `out_yang_smoke`, `out_yang_smoke2`,
`out_yang_court_f20`, `out_pen_f{5,10,20,40,80}`.
