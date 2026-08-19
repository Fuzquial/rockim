# Historique des versions rockim — chantier « coupe PDC » (aout 2026)

Chaque correctif est **opt-in** : cle absente = comportement historique, suite de
non-regression `fast` bit-identique. Trois binaires coexistent, aucun n'ecrase
les precedents.

## Les binaires

| exe | contient | runs qui l'ont utilise |
|---|---|---|
| `rockim_tun.exe` | etude tunnel EDZ (in situ + excavation) + `cutterThick` | `out_cut_heilman`, `out_cut_v2`, `out_cut_v3`, `out_cut_nogc`, tous les runs tunnel |
| `rockim_a1.exe` | + `gcSurfaceRefresh` (A') + `jointContactPenalty` (EPFL) | `out_cut_a1`, `out_cut_epfl` |
| `rockim_a2.exe` | + `toolImpulseCap` (A) | `out_cut_a2` |

**Aveu a consigner.** `rockim_a1.exe` a d'abord ete compile avec A' seul, puis
**ecrase** en y ajoutant EPFL. Le binaire A'-seul n'existe plus. Comme les deux
cles sont opt-in et independantes, `out_cut_a1` reste reproductible bit pour bit
avec l'exe actuel (`gcSurfaceRefresh = eager`, `jointContactPenalty` absent donc
`fixed`) — mais l'instantane de source correspondant, lui, est perdu. D'ou les
instantanes systematiques ci-dessous a partir de la.

## Les instantanes de source

| dossier | etat |
|---|---|
| `v0_tunnel/` | avant A' — ce que contient `rockim_tun.exe` |
| `v2_epfl/` | A' + EPFL — ce que contient `rockim_a1.exe` |
| `v3_impulse/` | A' + EPFL + A — ce que contient `rockim_a2.exe` |

(`v1` manquant : c'est l'etat A'-seul, perdu par l'ecrasement decrit plus haut.)

## Les trois correctifs

### A' — `gcSurfaceRefresh = legacy | eager`

Rafraichir le cache des faces liberees des qu'un joint **se separe** (`nDead_`)
et non quand un autre **casse** (`nBroken_`), et supprimer la grille `% 8`.
Vise la latence non bornee documentee par le commentaire du code lui-meme.

**Resultat : AGGRAVE.** Vitesse nodale 2 544 -> 4 777 m/s a la trame 10,
fragments 368 -> 442, injection outil 77 -> 85 kJ/m. L'injection suit le NOMBRE
de basculements cohesif/contact, pas leur retard : declarer plus tot en cree
davantage. Correctif conserve comme instrument de mesure, defaut `legacy`.

### EPFL — `jointContactPenalty = fixed | adaptive`

k- = k+(D) = (1-D)·pj : la penalite de compression suit la secante
d'endommagement, supprimant le saut de raideur a dn = 0.

Source : T. Ghesquiere-Dierickx, J.-F. Molinari & G. Anciaux, *Stability of
Extrinsic Cohesive-Zone Model with Penalty-Based Contact in Explicit Dynamic
Fragmentation Simulations*, arXiv:2511.14323v1, 18 novembre 2025, EPFL,
section 4 « Adaptive Contact Penalty ».

**Resultat : DE LOIN LE MEILLEUR.** Injection outil 77 286 -> 5 618 J/m (÷14),
vitesse max 2 544 -> 229 m/s (÷11), fragments 368 -> 106, zero element enfui sur
15 trames sur 16, et un PLATEAU de fissuration au lieu d'une inflation. Bat meme
la suppression totale du contact general, parce qu'il traite les 18 108 joints
VIVANTS et non les seules 285 faces exterieures.

Reserve des auteurs, a citer avec le resultat : *« it cannot be viewed as a
viable long-term remedy »* — l'interpenetration croit avec D. Mesure : profondeur
mediane doublee (0,010 -> 0,023 mm), mais nombre de noeuds concernes divise par
deux et maximum divise par 2,2.

### A — `toolImpulseCap = kappa`

|Fc| <= kappa · 2 · |v_outil| · m_i / dt. L'ecretage historique bornait la
PENETRATION (0,6 h) ; celui-ci borne l'IMPULSION, la grandeur qui lance les
noeuds. Borne physique du choc elastique contre une masse infinie, pas un
reglage. kappa = 1 est la borne stricte.

En cours d'evaluation (`out_cut_a2`).

## Les mesures qui ont guide le chantier

| run | injection outil | / travail corps rigide | v max | fragments |
|---|---|---|---|---|
| `out_cut_v3` (reference) | 77 286 J/m | **408** | 2 544 m/s | 368 |
| `out_cut_nogc` (contact general coupe) | 11 944 | 219 | — | 45 |
| `out_cut_a1` (A') | 84 745 | 416 | 7 153 | 442 |
| `out_cut_epfl` (EPFL) | 5 618 | 110 | 229 | 106 |
| `out_cut_a2` (EPFL + A) | en cours | | | |

Le residu du bilan B4 lit `[OK]` a 1e-10 % dans TOUS ces runs : la pompe siege
dans un canal comptabilise, donc le moniteur ne peut pas la voir. Les deux
observables qui l'exposent sont le rapport injection/travail-de-corps-rigide et
la vitesse nodale rapportee a 2 v_outil.

## Deux impasses ecartees par la mesure, a ne pas refaire

**Le tunneling n'existe pas ici.** Meme au pire cas jamais mesure (24 132 m/s),
un noeud parcourt h/8 par pas ; a 229 m/s, h/812. Toute la famille detection
continue de collision (CCD, volumes balayes, hierarchies englobantes) est sans
objet.

**L'ensemble actif n'a pas de trou.** Les 50/285 faces activees du run EPFL
inquietaient ; le rapport faces actives / fragments est en fait stable d'un run
a l'autre (v3 : 202/368 ; EPFL : 50/106). La couverture suit l'activite, elle ne
la limite pas. Et l'interpenetration reelle, mesuree comme un noeud strictement
interieur a un triangle d'un AUTRE fragment, est bornee a 0,6 h partout — donc
il n'y a pas de traversee, seulement de l'ejection.
