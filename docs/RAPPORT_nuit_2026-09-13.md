# RAPPORT DE NUIT — 12 → 13 septembre 2026 (rockim_g1, impact Yang 2026)

Fernando, 2 h : « je te laisse la main toute la nuit, corrige tout, demain matin tout doit être en
ordre, un bel essai 3D impact très raffiné ; trois tentatives par problème, puis le suivant ; et
refais un audit ». Ce rapport est le point d'entrée du matin. Il se lit seul ; les détails sont dans
`docs/ETAT_yang2026_2026-09-11.md` §13, `docs/PLAN_loi_joint_2026-09-12.md`, `docs/AUDIT_2026-09-13.md`,
le CHANGELOG et le journal du vault (`JOURNAL/2026-09-12.md`).

## 1. Ce qui tourne au réveil

**Le run s = 1** : `out_yang2026_v3`, journal `results/yang2026_v3.log`, binaire `rockim_g1y15.exe`,
deck `configs/yang2026_impact_v3_plastic.cfg` (120 185 tets, 232 408 facettes, dt 1,0 ns, T 400 µs,
40 trames, 14 fils). Lancé à 2 h 42. Rythme mesuré avant les gains de l'audit : 0,46 µs/min, soit
10 à 14 h pour 400 µs — **au réveil il sera vers 150-250 µs** : l'insert a touché la roche (~60 µs),
le cône et le cratère se forment, le point de retour du bit (~255-270 µs) est proche ou passé, le
rebond n'est pas encore lu. Il continue tout seul ; les figures se font à tout instant sur les trames
écrites :

```bash
python tools/fig_joints_cuts.py out_yang2026_v3 --title "s = 1, v3-P" --stem results/fig/coupes_s1_v3
```
```bash
python tools/fig_joints_only.py out_yang2026_v3 --title "s = 1, v3-P" --stem results/fig/joints_s1_v3
```
```bash
python tools/yang_report.py out_yang2026_v3 results/yang2026_v3.log
```

Garde-fous armés : borne physique KE ≤ KE₀ + sources (`budgetAbortPct = 2`) — s'il aborte, le journal
le dit en clair ; le moniteur de session relève la progression toutes les 15 min.

**État à 7 h 15 (113 µs, 12 trames, 0,25 µs/min en phase de fracture → fin vers 20 h)** :

| critère | s = 1 à 113 µs | jumeau s = 2,5 (200 µs) | Yang 9 m/s |
|---|---|---|---|
| σ_zz max au bit | 176 MPa (à 34 µs) | 173 | ~160 |
| vitesse d'indentation max | **7,37 m/s** (106 µs) | 6,27 | 5,62 |
| enfoncement | 0,44 mm (en cours) | 0,83 | ~1,0 |
| joints rompus / fragments | 1 023 / 362 | 294 / 124 | — |
| traction / cisaillement | 36 % / 64 % | 45 / 55 | — |
| rayon de la zone rompue | 5,2 mm | — | cratère ~7 |
| profondeur de la zone rompue | 5,5 mm | — | — |
| pulvérisés | 0 | 0 | ~360 |
| bilan | joints −1,1 J en cours de charge (stocké), contact −2,2, frottement −1,7 ; pas d'abort | joints +1,95, leapfrog +0,85 | — |

Lecture : le cratère se forme sous la pointe sphérique (zone de 5 mm de rayon, 5,5 mm de
profondeur, deux tiers cisaillement, coupes `results/fig/coupes_s1_v3_110us*.png`), aucune radiale
encore, aucune pulvérisation (δ_m = h·ε : à 1,4 mm de maille il faut 1 % de distorsion déviatorique).
La vitesse d'indentation dépasse Yang de 31 % (et le jumeau s = 2,5 de 18 %) : le bit maillé est 9 %
plus léger que le leur (audit D), le contact roche est plus souple par phase (`potStiffnessByPhase`),
et le maillage fin résiste moins — à trancher avec la suite du run (le point de retour dira ce que la
roche a absorbé).

## 2. La loi retenue, et pourquoi (le conseil)

La référence intrinsèque du 11/09 **créait de l'énergie** : relancée sous la borne KE elle aborte à
81 µs (4,3 J créés, poste joints −16 J). Bissection 3D en 11 variantes : la clé est
`jointShearUnload = origin` (l'éq. 18 de Yan transcrite), la même qu'en 2D ; les conventions
Solidity (`midedge` surtout) n'en modulent que l'amplitude. Mathématiquement (plan §1.3, vérifié par
le conseil) : sous `origin` la sécante τ_lim(σ_n)/s_max suit la compression courante, ∂τ/∂δ_n ≠ ∂σ/∂δ_s,
aucun potentiel ; un cycle à glissement fixé crée ½(k₂−k₁)s². `plastic` (retour radial, Φ = τ·ṡ_p ≥ 0)
est conforme ; Solidity elle-même est non conservative au second ordre dans le même régime (leur code
n'a aucun bilan d'énergie).

Le conseil (4 critiques + synthèse, REVISE) a ajouté : (M7) le DIF continu réécrit dnE à chaque pas,
donc la sécante de mode I est une raideur qui suit ft(t) — corrigé par `jointSecantRatchet = on`
(sécantes non croissantes, Φ ≥ 0) ; (D1) plastic + coulomb n'avait jamais tourné ; le test dt/2 du
terme leapfrog passe par construction ; le « miroir 2D/3D exact » est faux.

**Matrice fracture × conservation (banc s = 2,5)** :

| loi | t | rompus | joints (J) | leapfrog (J) | verdict |
|---|---|---|---|---|---|
| `origin` + conventions (référence) | 81 µs | 2 273 | −16,3 | +7,0 | abort |
| `origin` sans conventions (B4) | 90 µs | 1 375 | −0,10 | +1,7 | sain, casse |
| `plastic` seul | 200 µs | 22 | +2,9 | +1,3 | sain, ne casse pas (e = 0,22 non convergé) |
| **`plastic` + coulomb** | 200 µs | **226** (65 % traction) | **+2,0** | **+0,9** | **sain, casse — retenu** |
| `origin` + coulomb + ratchet | 82 µs (arrêté) | 41 | — | — | candidat D3, non conclu |

Loi du run s = 1 = **v3-P** : `plastic` + `jointShearRange = coulomb` (garde levée) +
`jointSecantRatchet = on` + rampe de naissance + les trois clés de l'audit.

**Jumeau s = 2,5 du deck v3-P** (`out_yang_bench_s25_v3P`, g1y15, 200 µs, 5 498 s à 4 fils) — les
trois clés de l'audit contre plastic + coulomb seul : 294 rompus (226), 124 fragments (93), 0
pulvérisé (2), joints +1,95 J (+2,02), leapfrog +0,85 (+0,94), résidu 2e-9 %, σ_zz 173 MPa, indentation
6,27 m/s à 108 µs, enfoncement 0,83 mm. Sain, +30 % de fracture, mêmes cinématiques. C'est la référence
s = 2,5 à mettre en face du run s = 1.

## 3. Ce qui a été corrigé cette nuit (binaires g1y11 → g1y15, tout opt-in, ancres 8/8 + 9ᵉ deck)

| binaire | clé / geste | défaut corrigé | preuve |
|---|---|---|---|
| g1y11 | `jointSecantRatchet = on` (2D+3D) | sécantes des éq. 17/18 qui remontaient avec σ_n (M1) ou ft(t) (M7) | conseil, PLAN §1.3 |
| g1y11 | garde coulomb ⇒ origin levée | la plage 3G_II/f_s est une longueur de référence, pas une raideur (M4) | plastic+coulomb : 226 rompus, +2 J |
| g1y11 | avertissement `origin` (formule M1) | rien n'est refusé (40 decks archivés rejouables) | — |
| g1y13 | `gcBirth = relay` | `penalty` injectait 1 J au premier contact ; `ramp` coupe le chemin d'effort au relais | codé, mesure A/B en cours (2 fils) |
| g1y13 | `facetAverage = max` en 2D | le critère moyenné diluait l'anneau hertzien (item du 11/09) | port du 3D |
| g1y14 | `jointNormalProxy = law` (2D+3D) | σ_n = pj·dn au lieu de 2pj·dn sous parabolic dans s_E/coulomb/DIF (audit A #1) | Solidity lit σ_tmp = pe·o/el |
| g1y14 | ratchet armé au cap | verrou possible à la naissance d'un joint inséré sous traction (A #4) | — |
| g1y14 | `bulkDamagePhase` | acier et carbure endommageables (audit D) | marge 1,04-1,6 à 175 MPa |
| g1y14 | `potStiffnessByPhase = min` | contact roche/roche 10× trop raide (E carbure pour toutes les paires, B #10) | — |
| g1y14 | fusion parallèle des forces de joint | 95 % du poste joints était série (30,7/32,4 ms à s = 1, audit C #1) | 8/8 + 9ᵉ deck IDENTIQUE |
| g1y14 | ping-pong des fragments, f_ en parallèle | recopies 96×256×72 o par face ; 11,5 Mo de zéros en série | 8/8 + 9ᵉ deck IDENTIQUE |
| g1y15 | `bulkDamagePhase` résolu après les phases | le deck de santé refusait la clé (1ʳᵉ tentative) | run de santé 20 µs OK |

Ancres : `results/bitid_g1y11..15.log` (8/8) et `bitid_jointlaw_g1y1*.log` (9ᵉ deck, les cinq clés du
deck v2 sous fracture). Outils nouveaux : `tools/fig_joints_only.py`, `tools/fig_joints_cuts.py`,
`tools/law_matrix.py`.

## 4. Ce qui reste, classé (pour la décision de Fernando)

1. **Le run s = 1 lui-même** : le laisser aller au bout (T 400 µs, ~10-14 h) ; lire les sept critères
   sur `history.csv` et les figures de coupes ; comparer à Yang (σ ~160 MPa, 5,62 m/s, 4,65 m/s de
   rebond, ~1 mm, cratère 7 mm, radiales 9-10 mm, ~360 pulvérisés). Le rebond exige T ≥ 600-750 µs
   (Yang lit jusqu'à 750) : un second run, ou une reprise, à décider.
2. **Maillage** (audit D) : le s = 1 est à 1,37 mm au cœur (Yang 1,0) et 1,39 mm sur l'insert (Yang
   0,7) — la pulvérisation δ_m = h·ε est ∝ h. Le générateur vise 1,0 mm au cœur (`make_impact_mesh.py`,
   « s = 1 reproduit l'article ») mais gmsh rend une médiane de 1,37 mm : un maillage conforme demande
   s ≈ 0,75 (~290 k tets, dt ×0,75, coût ×3 : ~35 h à ce rythme).
3. **Masses** : bit + insert + circlip 1,367 kg (Yang 1,51), piston 1,057 (facettisation) : corriger par
   ρ ou par le maillage avant la comparaison des vitesses.
4. **`gcBirth = relay`** : A/B fait (banc s = 2,5, plastic + coulomb, 200 µs) : neutre — 231 rompus
   contre 226, joints +2,10 J contre +2,02, leapfrog +0,88 contre +0,94, 211 paires calées (facteur
   moyen 0,53). Rien à changer sur cette loi ; à remesurer sur une configuration qui tue des milliers
   de joints en compression.
5. **Mesures à coder** (audit B) : `dtContactAudit` (raideur de contact réellement vue par nœud — le
   budget de dt suppose 2 paires/nœud), `gcBirthWork` (injection de la rampe), `gcFricSplit`.
6. **Forme de D** (éq. 4 vs Fig. 4a linéaire : Dmax à 0,107 mm ou 0,4 mm), Cd = 1, δ_m = h·√(2/3)‖dev ε‖
   : trancher avant de comparer la masse de fragments ; `fragBrush*` à réarmer pour le critère 5.
7. **2D** : `relay` absent ; deux règles `midedge` ; campagne du 03/09 à relire (elle tournait sur
   `origin` + amortissements).
8. **Consolidation à froid** : `JointLaw.hpp` (le miroir 2D/3D n'est pas exact), pilote de cycle au
   niveau facette, instrumentation ∮τ·dδ_s par facette, CSR de grpsOfVert_, AVX2 à re-certifier.

## 5. Où sont les choses

- Code : `simulations/FDEM/rockim_g1`, branche `g1`, poussée sur `origin/g1` (Fuzquial/rockim).
- Docs : `docs/ETAT_yang2026_2026-09-11.md` (§13 = la nuit), `docs/PLAN_loi_joint_2026-09-12.md`,
  `docs/AUDIT_2026-09-13.md`, `CHANGELOG.md`, `DOCUMENTATION_rockim.md` (nouvelles clés).
- Vault : `JOURNAL/2026-09-12.md`, `ETUDES/rockim-impact-yang-insert-unique.md`.
- Figures : `results/fig/joints_*`, `joints_aretes_*`, `coupes_*` (bancs s = 2,5), `coupes_s1_v3*`
  (run s = 1, à générer).
