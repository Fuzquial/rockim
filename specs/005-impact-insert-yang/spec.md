# Feature Specification : Impact à insert unique (Yang et al. 2025-2026)

**Feature Branch** : `005-impact-insert-yang`

**Created** : 2026-08-21, après lecture des deux articles transmis par
F. Uzquiano.

**Status** : Draft — spécification à valider, rien n'est implémenté.

**Input** : Donner à rockim la capacité de reproduire les essais d'impact à
insert unique d'Imperial College / Mines Paris-PSL :

* Yang, Xiang, Naderi, Wang, **Aising**, Ugarte, Latham — *Multi-criteria
  validation of hi-fidelity numerical model of impact breakage*, IJRMMS 191
  (2025) 106125 — calcaire St Anne et grès Rhune, 7 critères de validation ;
* mêmes auteurs — *High-fidelity modelling of fragmentation and pulverisation
  in hard granite under percussion loading*, IJRMMS 206 (2026) 106660 —
  granite Kuru Grey, modèle de pulvérisation, rebond non linéaire.

Les données expérimentales des deux papiers viennent de **Mines Paris-PSL**
(Aising et al.) : essais de frappe instrumentés, masse de fragments, CT,
caméra rapide. Contrainte posée par F. Uzquiano : **conserver l'insertion
adaptative** (schéma extrinsèque Yan 2023) — leur code Solidity est
intrinsèque, rockim gardera son schéma.

---

## 1. Ce que rockim possède DÉJÀ (vérifié dans le source le 2026-08-21)

L'inventaire change radicalement le périmètre : **plus de la moitié du papier 1
est déjà implémentée**, souvent à l'identique.

| brique de l'article | état dans rockim | où |
|---|---|---|
| FDEM 3D tets + joints 6 nœuds, critère Mohr-Coulomb à cut-off (leur éq. 1) | fait | `Fdem3dSolver.cpp` |
| insertion adaptative 3D (liaison de nœuds) | fait, testée (`fdem3d_tension_adaptive`) | `buildBindingTables` / `insertionSweep` |
| DIF en taux de déformation (leurs éq. 2-3) | fait, opt-in `strainRateDIF = yang \| yang-fig2` — **exige déjà `insertion = adaptive`** | `Fdem3dSolver.cpp:63-91` |
| exposant traction 0,07 vs 0,17 | les DEUX : `yang` = 0,07 littéral (papier 1), `yang-fig2` = 0,1707 | `:72` |
| **tri des fragments** (leur §2.3 : anti-gravité, v0, particule de référence, β) | **fait, intégral** — clés `fragBrushStart/V/Accel/Dir/Beta/ZeroV` | `armBrush` / `bodyForces:2975+` |
| fragments = composantes connexes | fait | `computeFragments` / `fragId_` |
| multi-corps par volumes physiques gmsh (aucun joint inter-corps, matériau par corps `groupPhase.<nom>`) | fait | `buildMeshFile:613-676` |
| vitesse initiale par corps `groupVel.<nom> = "vx vy vz"` (le piston !) | fait | `:679+` |
| gravité (`gravity`) | fait | `:1248` |
| contact frottant régularisé (`contactMu`, tanh) | fait | `:2318, :2634` |
| scénario PERCUSSION avec outil maillé (`toolShape = none`) | fait | `Fdem3dSolver.hpp:549` |

Note de chronologie : le papier 2 (2026) imprime l'exposant traction **0,17**
dans son éq. 2 là où le papier 1 (2025) imprimait 0,07 — la lecture « coquille
d'exposant » consignée dans rockim le 2026-08-18 est donc confirmée par les
auteurs eux-mêmes. Les campagnes utiliseront `yang-fig2`.

Piège de source : la Table 5 du papier 1 a les colonnes acier/carbure
**inversées** (acier à 15 250 kg/m³ et 600 GPa). La Table 1 du papier 2 fait
foi : carbure ρ 15 250, E 600, ν 0,20 ; acier ρ 7 850, E 200, ν 0,29.

---

## 2. Ce qui MANQUE — cinq lots

### WP1 — Pulvérisation : endommagement volumique des tétraèdres (papier 2, le cœur)

Leur éq. 3-4 : au-delà d'une déformation d'amorçage ε_d, la contrainte devient
σ = C_d (1−D) σ̄ avec D linéaire en « déplacement effectif » δm entre δm⁰ et
δm^f, plafonné à D_max (0,9 pour Kuru Grey). Ce n'est PAS un modèle de joints :
il dégrade la raideur des **éléments volumiques** sous l'insert et représente le
sur-broyage/la poudre que les joints méso ne peuvent pas capturer. C'est ce qui
produit le rebond **non linéaire** du bit (fig. 1c) et l'éjection des fragments
pulvérisés — les deux signatures du granite absentes des roches sédimentaires.

* clés : `bulkDamage = off | yang` (défaut off = bit-identique),
  `bulkDamageDelta0`, `bulkDamageDeltaF`, `bulkDamageDmax`, `bulkDamageCd` ;
* leur δm est un déplacement effectif d'élément (max/min des composantes) ;
  l'approximation rockim (h_e × déformation équivalente) sera documentée et
  confrontée sur l'UCS avant tout impact ;
* l'énergie dissipée par la dégradation = **poste séparé** du bilan (même
  discipline que brushWork_, jamais dans sumW) ;
* interaction avec `crushCap` (déjà présent en 3D) à trancher : le cap
  volumique et la dégradation de raideur ne doivent pas se cumuler en double
  comptage.

### WP2 — Liaison entre corps (insert carbure ↔ bit acier)

Aujourd'hui `buildFromTets` ne crée AUCUN joint entre volumes physiques : les
corps n'interagissent que par contact. Or l'insert est brasé au bit, et le
circlip/la plaque tiennent l'ensemble. Clé `groupBond.<A>.<B> = joints`
créant les joints conformes à l'interface entre deux volumes physiques
maillés conformément (résistance = phase dédiée, très haute pour un brasage).
Défaut absent = comportement actuel, bit-identique.

### WP3 — Jauge de contrainte et histoire par corps

Deux des sept critères sont des vitesses du bit (indentation, rebond), un
troisième est la contrainte maximale à mi-bit (leur jauge, fig. 8).
`history.csv` ne connaît aujourd'hui que l'outil analytique (toolVx…).
Ajouter : vitesse moyenne par groupe (`vz_<nom>` par corps), profondeur
d'indentation (z min de l'insert), et une sonde σzz moyennée dans une tranche
du bit (`gauge.<nom> = z0 z1`). Colonnes ajoutées seulement si les clés sont
posées.

### WP4 — Générateur de maillage et decks

`tools/make_impact_mesh.py` (gmsh, MSH 2.2, volumes physiques nommés) :

* piston Φ26,5 × 260 mm (acier) ; bit Φ30 × 265 mm (acier) + insert
  hémisphérique R 8,51 mm et circlip (carbure, WP2) ; plaque de charge
  119,9 × 40 × 6 mm percée Φ31 (acier) ; roche cylindre Φ250 × 150 mm ;
* raffinement : 1 mm dans l'hémisphère Φ25 sous l'impact → 2 mm jusqu'à
  Φ50 → 10 mm au bord ; insert 0,7 mm (leur fig. 6 ; ~230 k éléments,
  dt ≈ 2,5e-9 s — leur Table 2) ;
* decks : `impact_stanne.cfg`, `impact_rhune.cfg` (Table 4 papier 1 : E 57/37
  GPa, ν 0,31/0,20, ft 7,0/7,5 MPa, c 18,8/33,6, GI 12/40 J/m², GII 800/1400,
  tanφ 1,0, μ 0,6), `impact_kuru.cfg` (Table 1 papier 2 : E 60, ν 0,24,
  ft 10,98, c 29,84, GI 50, GII 1000, tanφ 1,85, **μ 0,18**, δm⁰ 0,014,
  δm^f 0,4, D_max 0,9) ;
* toutes les campagnes : `insertion = adaptive`, `strainRateDIF = yang-fig2`,
  `gravity = 9.81`, `groupVel.piston = "0 0 -V"`, tri des fragments
  `fragBrushV = 2.5e-3`, `fragBrushAccel = 98.1`, `fragBrushBeta = 0.8`
  (leurs valeurs §2.3).

### WP5 — Dépouillement multi-critères

Les 7 métriques de leur fig. 8 : contrainte max au bit, vitesse d'indentation,
vitesse de rebond, profondeur d'indentation, masse de fragments, longueur de
fissure radiale, rayon de cratère — plus la classification des fissures
(radiales / médianes / latérales / side) par position et orientation des joints
rompus, et les planches façon fig. 9-16 (vue de dessus, coupe verticale,
rouge = traction, jaune = cisaillement). Outils Python sur les VTU, style
`figlib` maison (PDF vectoriel, Computer Modern).

---

## 3. Trajectoire — la même que les auteurs

Le papier 1 (calcaire, grès) n'utilise PAS le modèle de pulvérisation : il
n'est nécessaire qu'au granite. Donc les « beaux impacts » arrivent AVANT le
gros lot :

1. **WP2 + WP3** (petits, code) → build dédié, suite fast 19/19, chemins par
   défaut bit-identiques ;
2. **WP4** (outillage) → maillage + blanc court pour mesurer le coût réel ;
3. **campagne papier 1** : St Anne à 2,475 / 6,235 / 10,66 m/s + Rhune à
   8,52 / 14 m/s — validation contre leurs fig. 9-12 (cratère en étoile,
   fissures radiales 20-40 mm, masses 0-1,6 g) ;
4. **WP1** (pulvérisation) → validation d'abord sur UCS statique, puis
   **campagne granite Kuru Grey** : 5,68 / 8,04 / 9 / 11 / 13 m/s — cible :
   le rebond NON monotone de leur fig. 1c et la masse de poudre ;
5. WP5 en continu dès la première campagne.

Coût de calcul, à mesurer au blanc : leur run = 230 788 éléments, dt 2,5e-9 s,
5 h chez eux (Table 2) ; estimation rockim 14 cœurs : 6-12 h par impact.
Environ 12 impacts pour les deux campagnes → planifier par lots nocturnes,
un job à la fois.

## 4. Garde-fous

* chaque capacité neuve est **opt-in, défaut bit-identique**, suite fast 19/19
  après chaque lot — la discipline des essais 0-3 du banc AbuAisha ;
* l'insertion adaptative est conservée partout ; le DIF l'exige déjà
  structurellement, la contrainte est donc gratuite ;
* binaire dédié par lot (`rockim_i1.exe`, …), jamais d'écrasement d'un exe en
  cours de run ;
* sauvegarde + REVERT.md avant chaque patch ;
* figures : PDF vectoriel Computer Modern, fissures en rendu éléments.
