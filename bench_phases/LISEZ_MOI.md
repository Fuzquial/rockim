# `bench_phases` — banc court du matériau PAR PHASE en éléments finis (`mode = fem3d`)

Écrit le 2026-09-06 avec la capacité elle-même. **Aucun run long** : les huit decks tournent en
moins d'une seconde chacun, sur des maillages de quelques centaines à quelques milliers de tets.
Le banc doit être vert **avant** de dimensionner la moindre campagne.

## Lancer

```powershell
$env:OMP_NUM_THREADS = "1"
python bench_phases\make_bar2.py bench_phases\bar2.msh                 # barre à 2 couches nommées
python bench_phases\make_bar2.py bench_phases\bar1.msh --sans-noms     # la même, sans $PhysicalNames
foreach ($d in "f3a_homogene_dur","f3b_homogene_mou","f2_deux_phases_identiques",
               "f3c_serie_reuss","f3d_serie_permutee","f5a_voronoi_1phase",
               "f5b_voronoi_3phases_egales","f5c_voronoi_3phases_contrastees") {
  & build\rockim.exe bench_phases\$d.cfg bench_phases\out_$d > bench_phases\log_$d.txt 2>&1
}
python bench_phases\depouille.py bench_phases      # 9 contrôles
python bench_phases\erreurs.py                     # 16 cas fautifs
```

## Ce que chaque contrôle falsifie

**F2 / F5 — neutralité, octet à octet.** Des phases aux propriétés **strictement identiques** à la
fiche globale doivent rendre `history.csv` **octet à octet** identique au deck sans `phases` — sur le
chemin fichier (F2) **et** sur le chemin voronoï (F5). L'indice de phase traverse la masse condensée,
la CFL, la pénalité de contact, l'impédance de Lysmer, la viscosité de volume et la loi : si un seul
produit avait été réassocié en chemin, l'égalité tomberait. La **variante qui doit échouer** est le
même deck avec un contraste : elle doit différer.

**F3 — borne de Reuss, en forme fermée. C'est le contrôle décisif.** Le mode de défaillance le plus
grave de ce chantier est silencieux : câbler `Elem.phase`, les champs `.vtu` et les bilans de
fractions, et **oublier** `laws_[e.phase]->stress`. La figure serait parfaite, la carte de grains
colorée, et tous les chiffres seraient ceux d'un bloc homogène — F2 et F5 **passeraient quand même**.

La barre à deux couches en série a une solution fermée : `1/E_app = f₁/E₁ + f₂/E₂`. Les mors
(bas encastré, haut à x,y bloqués) raidissent les extrémités ; on ne compare donc pas à la théorie
nue mais aux **deux runs homogènes du même maillage**, ce qui élimine l'effet de mors au premier
ordre. Mesuré : **44,158 GPa contre 43,678 attendu, soit 1,10 %**. Et
`E_app(série)/E_app(dur) = 0,5253`, qui vaudrait **1 exactement** si la loi n'était pas indexée.

**F3d — commutativité.** Les deux couches permutées par `groupPhase` : la série est commutative,
`E_app` doit être inchangé (−1,48 % mesuré, l'écart résiduel est l'asymétrie des mors). Vérifie du
même coup que `groupPhase.<groupe> = <phase>` est réellement pris en compte.

**F4 — CFL sur `c_P` MAXIMALE des phases.** En passant d'un bloc uniforme à 60 GPa au contraste
Red Bohus (quartz 83,1 GPa), `dt` doit **chuter** de √(60/83,1) = 0,8497. Mesuré 0,8497. Laissée sur
la fiche globale, la CFL serait 18 % trop grande : le schéma centré n'explose pas forcément, il
devient **bruyant**, et ce bruit se lit comme de la fissuration.

**Gardes imprimées par le solveur lui-même** (elles lèvent, elles n'avertissent pas) : audit de masse
`sum(m)` contre `sum_p rho_p V_p` recalculé indépendamment (3,7·10⁻¹⁶ mesuré) ; aire des faces
extérieures = aire de la boîte en voronoï (6·10⁻¹⁶) ; aucune face vue plus de deux fois ; au moins une
face partagée entre grains distincts ; en `mesh = file` multi-groupes, au moins une face conforme
entre groupes.

## `erreurs.py` — 16 cas qui DOIVENT échouer, avec le BON message

Une clé lue, validée et **sans effet** est le motif interdit numéro 1 du dépôt : le deck a l'air de
dire quelque chose et le calcul ne l'entend pas. Le test vérifie le **message**, pas seulement le code
de retour — un refus qui envoie chercher au mauvais endroit coûte une demi-journée.

`phases` sur `mesh = grid` · `phase.<nom>.law` · clé de loi écrite par phase · nom de phase inconnu ·
`gbAlphaTen` · `gb.<a>.<b>.ft` · `groupBond.<A>.<B>` · groupe physique sans phase · `phases` sans
groupes nommés · Weibull × phases sans `phaseWeibull` · `E < 0` · `nu ≥ 0,5` · `fraction` manquante ·
**fiche de phase sans clé `phases`** (sans elle `PhaseSet` nomme sa fiche unique « rock », donc
`phase.rock.E` serait acceptée, consommée et **inerte** — le trou le plus fin du lot) ·
**`groupPhase` hors `mesh = file`** · **`phaseWeibull` posée sans objet**.

## Ce que le banc NE prouve PAS

- **Ce n'est pas un GBM cohésif.** Nœuds partagés = aucun joint = aucune frontière de grain
  intrinsèquement faible. La frontière concentre la contrainte, elle ne s'ouvre pas.
- **Verrouillage volumique.** Les tets linéaires sans B-bar sur-raidissent d'autant plus que `nu` est
  grand : à `nu` variable entre phases l'artefact est corrélé à la phase et **sous-estime** le
  contraste. F3 tourne à `nu` commun, il ne mesure donc pas ce biais.
- **Objectivité de maillage.** F3 est élastique : ni bande de fissuration, ni érosion. La qualité des
  tets de la tessellation (rapport diamètre inscrit / lc médian 0,734, min 0,418) est imprimée mais
  n'a pas été confrontée au `lc_c` dépendant du taux.
