# Résultats — groupe `diag_contact` (origine de l'injection d'énergie par le contact, tunnel réduit)

Calculs du 2026-10-06, 22 h 54 à 23 h 20 (UTC), session cloud Linux, binaire `f638b0f`
(code = `e9ba9a6`), **4 fils OpenMP**, environ 5 min chacun. Sorties brutes dans
`sim_out/diag_contact/` (hors git). Chiffres extraits du bloc de bilan de fin de `.log` (V2/B4) ;
vitesse nodale max lue sur les 14 trames VTU (champ `velocity`), le solveur ne l'imprime pas.

## Maillage : régénéré, donc nouvelle base

`meshes/tunnel_hs_red.msh` n'est pas sous git (`*.msh` ignoré). La recette exacte du maillage de la
session de rédaction (22 730 triangles, h = 0,40 m) est introuvable ; il a été régénéré par
`tools/make_unstructured_mesh.py tunnelhs 100 100 0.40 12 2.0 … 1` : **23 004 triangles**, h inscrit
min/méd/max 0,129/0,235/1,498 m, pire angle 31,9°, l_cz/2 = 0,278 m < hFine (TROP GROSSIER, comme
l'original). Gmsh y laisse 4 nœuds orphelins (centres des arcs du fer à cheval) que rockim refuse ;
ils ont été retirés (nœuds non référencés par un triangle, renumérotation). `meshes/drop_orphans.py`
d'`etude_lois_fem` ne traite que les tétraèdres.

Pour que la comparaison reste valable, la **base `red_tip16` a été rejouée sur ce maillage**
(`base_tip16` ci-dessous). Sur le maillage de la session de rédaction : contact 1,94e6 J/m,
cohésif 9,57e5 J/m ; ici 2,61e6 et 9,98e5 J/m. Même ordre de grandeur et même signature : le travail
net du contact vaut 2,6 fois l'énergie cohésive.

Chaque variante ne change qu'une clé par rapport à la base (vérifié dans `config_effective.cfg`) ;
l'avertissement « jointShearUnload = origin est NON CONSERVATIF » disparaît bien avec `plastic` et
avec `ratchet` (« secantes de decharge NON CROISSANTES … Phi >= 0 »).

## Bilan (J/m, t = 0,25 s ; écart à la base entre parenthèses)

| run | clé changée | travail net du contact | travail des joints (cohésif + amortisseur) | cohésif | frottement (dans le contact) | Cundall | KE finale | résidu B4 | joints rompus / insérés | v nodale max (trames) |
|---|---|---|---|---|---|---|---|---|---|---|
| `base_tip16` | — | 2,611e6 | 1,075e6 | 9,983e5 | 2,708e6 | 9,127e6 | 1,906e5 | −1,328e6 (6,62 %) | 11 536 / 13 117 | 25,85 m/s |
| `tip16_gcpen` | `gcBirth = penalty` | 2,058e6 (**−21,2 %**) | 1,041e6 | 9,679e5 (−3,0 %) | 2,615e6 | 8,550e6 (−6,3 %) | 1,080e5 | −1,340e6 (7,15 %) | 11 462 / 13 193 | 11,88 m/s |
| `tip16_mu0` | `contactMu = 0` | 1,235e7 (**+373 %**) | 1,508e6 | 1,379e6 (+38,2 %) | 0 | 1,957e7 (+114 %) | 2,206e5 | −1,333e6 (3,37 %) | 12 202 / 13 825 | 23,80 m/s |
| `tip16_plastic` | `jointShearUnload = plastic` | 5,582e6 (**+114 %**) | 1,257e6 | 1,178e6 (+18,0 %) | 4,059e6 | 1,235e7 (+35,3 %) | 1,503e5 | −1,304e6 (5,02 %) | 11 535 / 13 403 | 11,27 m/s |
| `tip16_ratchet` | `jointSecantRatchet = true` | 2,187e6 (**−16,2 %**) | 1,011e6 | 9,461e5 (−5,2 %) | 2,497e6 | 8,515e6 (−6,7 %) | 1,663e5 | −1,325e6 (7,02 %) | 11 784 / 13 406 | 11,18 m/s |
| `tip16_mu0_gcpen` | `contactMu = 0` + `gcBirth = penalty` | 9,556e6 (**+266 %** ; −22,6 % sur `mu0`) | 1,351e6 | 1,237e6 (+23,9 %) | 0 | 1,625e7 (+78,0 %) | 1,864e5 | −1,328e6 (4,00 %) | 11 679 / 13 286 | 30,15 m/s |

Contact en potentiel de Munjiza (p = kt = 1e10 N/m), insertion adaptative (pénalité 4 E/h),
amortissement local `dampingLocal = 0.15` (deck ligne 99). Signe : le bilan imprime le contact en −travail injecté ; « travail net du
contact » est ici le travail INJECTÉ (positif = source).

## Lecture

1. **Aucune des quatre hypothèses n'explique l'injection seule.** La meilleure variante (`gcpen`)
   retire 21 % ; le contact injecte encore 2,1 fois l'énergie cohésive.
2. **Le frottement n'est pas la source, il la masque.** Sans frottement, le travail net injecté
   est multiplié par 4,7 (1,24e7 J/m) et Cundall dissipe deux fois plus. Le frottement dissipe
   environ 2,7e6 J/m dans la base, autant que l'injection nette : la partie NORMALE du contact
   injecterait donc de l'ordre de 5e6 J/m.
3. **La décharge `origin` non conservative n'est pas la source principale.** `plastic` double
   l'injection (+114 %) ; le cliquet (`ratchet`, piste n° 1 du chapitre coupe) ne la réduit que de
   16 %. Le cliquet supprime bien la création d'énergie de la décharge (Φ ≥ 0), et ce qui reste ne
   vient donc pas d'elle.
4. **La naissance du contact sur un joint mort contribue (−21 % avec `gcBirth = penalty`)** et
   divise par deux la vitesse nodale max (25,9 → 11,9 m/s). `plastic` et `ratchet` divisent aussi
   la vitesse max par 2,3 : le pic de 25,9 m/s de la base est un événement local, pas l'injection
   globale.
5. [Lecture corrigée par le diagnostic 3 : le résidu varie avec dt.] **Le résidu B4 est le même dans les cinq calculs (−1,30 à −1,34e6 J/m)**, alors que les autres
   termes varient d'un facteur 2 à 5. C'est un terme fixe non compté par le bilan (piste : travail
   de la pré-contrainte in situ ou du relâchement d'excavation), indépendant du contact. À
   départager avant de lire le bilan comme un diagnostic fin.
6. Piste suivante suggérée : la partie normale du potentiel au relais joint → contact (65 % des
   joints meurent en compression ; charge normale lâchée au relais 1,69e6 kN/m cumulés dans la
   base). Un calcul `contactMu = 0` + `gcBirth = penalty` séparerait naissance et potentiel normal.

## Diagnostic 2 : `tip16_mu0_gcpen` (2026-10-06, 23 h 35, même maillage de 23 004 triangles)

Deck `configs/diag_contact/tip16_mu0_gcpen.cfg`, seul `meshFile` repointé vers le maillage
régénéré de la base (la version versionnée depuis, `tunnel_hs_red.msh` d'origine à 22 730 triangles,
n'est pas celle de `base_tip16`). Durée 263 s.

| comparaison | travail net du contact | effet de `gcBirth = penalty` |
|---|---|---|
| avec frottement : `base_tip16` → `tip16_gcpen` | 2,611e6 → 2,058e6 | −21,2 % (−5,5e5 J/m) |
| sans frottement : `tip16_mu0` → `tip16_mu0_gcpen` | 1,235e7 → 9,556e6 | −22,6 % (−2,8e6 J/m) |

La naissance du contact sur un joint mort pèse la même PART (≈ 22 %) avec ou sans frottement. Sans
frottement et avec naissance par pénalité, le contact injecte encore 9,6e6 J/m, soit **7,7 fois
l'énergie cohésive** (1,24e6 J/m) : **la source dominante est la réponse normale du potentiel de
contact elle-même**, pas la naissance ni la décharge des joints. La vitesse nodale max monte à
30,1 m/s (trame 7). Le résidu B4 reste à −1,328e6 J/m, toujours identique.

Piste suivante : la raideur de pénalité normale (p = 1e10 N/m) contre le pas de temps et la taille
des éléments (biais O(dt) du potentiel de Munjiza, que le solveur dit « petit » : il ne l'est pas
ici), par exemple un balayage de `potPenaltyFactor` (× 0,1 / × 10) ou `dtFactor` / 2 sur `tip16_mu0`.

## Diagnostic 3 : raideur normale du potentiel et pas de temps (2026-10-07, 00 h 38 à 01 h 35)

Base : `tip16_mu0` (sans frottement, pour isoler la partie normale). Même maillage de 23 004
triangles (`meshFile` seul repointé). Les trois calculs sont allés au bout (T = 0,25 s), aucun
arrêt à 30 min. Le premier essai de `pot01` a été perdu dans un redémarrage du conteneur et relancé.

| run | clé changée | p (N/m) | dt (s) | travail net du contact | cohésif | Cundall | KE finale | résidu B4 | joints rompus / insérés | v nodale max (trames) | durée |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `tip16_mu0` | — | 1e10 | 5,586e-6 | 1,235e7 | 1,379e6 | 1,957e7 | 2,206e5 | −1,333e6 (3,37 %) | 12 202 / 13 825 | 23,80 m/s | 311 s |
| `tip16_mu0_pot01` | `potPenaltyFactor = 0.1` | 1e9 | 5,586e-6 | 1,435e6 (**−88,4 %**) | 8,878e5 (−35,6 %) | 9,743e6 (−50,2 %) | 2,071e5 | −1,296e6 (5,98 %) | 11 068 / 13 044 | 12,59 m/s | 690 s |
| `tip16_mu0_dt05` | `dtFactor = 0.1` | 1e10 | 2,793e-6 | 7,266e6 (**−41,2 %**) | 8,630e5 (−37,4 %) | 1,367e7 (−30,1 %) | 1,601e5 | **−6,543e5** (2,24 %) | 11 755 / 13 270 | 22,83 m/s | 1 472 s |
| `tip16_mu0_pot10` | `potPenaltyFactor = 10` | 1e11 | 3,484e-6 | 5,574e7 (**+351 %**) | 2,186e6 (+58,5 %) | 6,948e7 (+255 %) | 4,341e5 | −8,464e5 (0,66 %) | 13 193 / 14 562 | 31,42 m/s | 1 272 s |

Lecture :

1. **L'injection est un artefact de la raideur normale du potentiel.** Elle suit p de près :
   ×8,6 de 1e9 à 1e10 N/m, ×4,5 de 1e10 à 1e11 N/m (pente log-log 0,93 puis 0,65). À
   p = 1e9 N/m, le contact n'injecte plus que 1,6 fois l'énergie cohésive (1,44e6 contre 8,9e5 J/m),
   contre 9,0 fois à 1e10 et 25 fois à 1e11.
2. **Elle dépend aussi du pas de temps** : dt divisé par 2 retire 41 % (un biais d'ordre 1 en dt
   retirerait 50 %). C'est la signature d'une erreur d'intégration du potentiel (raideur de contact
   mal résolue par le pas), pas d'un terme physique. Avec `pot10`, le solveur réduit lui-même dt
   (3,48e-6 au lieu de 5,59e-6 s), ce qui n'empêche pas l'injection de quadrupler.
3. **Le résidu B4 « constant » des diagnostics 1 et 2 ne l'est pas** : il est divisé par 2 avec dt
   (−1,33e6 → −6,5e5 J/m) et tombe à −8,5e5 avec `pot10` (dt réduit). C'est aussi une erreur
   d'intégration en O(dt), et pas un terme de pré-contrainte ou d'excavation non compté comme le
   suggérait la lecture 5 plus haut. Cette lecture est donc retirée.
4. L'énergie cohésive et le nombre de joints rompus bougent avec p (−36 % et −9 % à p/10 ; +59 %
   et +8 % à 10 p) : le faciès lui-même dépend de la raideur de contact sur ce banc. Le chapitre ne
   doit pas comparer de faciès de tunnel obtenus à des p différents.
5. Recommandation : p = 1e9 N/m (`potPenaltyFactor = 0.1`) ramène l'injection à l'ordre de
   l'énergie cohésive sans changer le pas de temps, pour un coût ×2,2 en temps de calcul (690 s
   contre 311 s, plus de joints restent en contact longtemps). Vérifier sur la base AVEC frottement
   (`red_tip16`) et sur le maillage de production avant de l'adopter ; vérifier aussi
   l'interpénétration maximale au contact à p = 1e9 (non imprimée ici).

## Diagnostic 4 : `potPenaltyFactor = 0.1` AVEC frottement (2026-10-07, 02 h 35 à 03 h 13)

Decks `configs/diag_contact/red_tip16_pot01.cfg` et `red_adapt_pot01.cfg` (une seule clé ajoutée
par rapport à `tunnel/red_tip16.cfg` et `tunnel/red_adapt.cfg`, vérifié par diff) ; base
`red_adapt` rejouée sur le même maillage (`red_adapt_base`, elle n'existait pas). Même maillage de
23 004 triangles pour les quatre calculs (`meshFile` seul repointé), même dt (5,586e-6 s), T =
0,25 s, aucun arrêt. `red_adapt` diffère de `red_tip16` par l'absence de facteur de pointe
(`insertionTipFactor = 1.6`, `insertionTipDamage = 0.5` retirés). Conteneur environ 2,3 fois plus
lent qu'en début de nuit : les durées ne se comparent pas à celles des diagnostics 1 et 2.

R_EDZ, demi-axes, déplacements de couronne, reins et radier : `tunnel_edz/tools/edz_metrics.py`
(dernière trame, joints rompus de `fdem_final_joints.csv`). U paroi moyenne : moyenne de |U| sur
les 86 nœuds des arêtes libres de la cavité (hors boîte extérieure), dernière trame.
**L'interpénétration maximale au contact n'est pas mesurable** : le solveur ne l'imprime pas, et
`rn`, `rs` de `fdem_final_joints.csv` sont les rapports de rupture, pas des ouvertures. Le
déplacement maximal de paroi en tient lieu, faute de mieux.

| run | travail net du contact | frottement | cohésif | Cundall | résidu B4 | joints rompus / insérés | R_EDZ p95 / max (m) | demi-axes h / v (m) | U max (couronne / reins / radier) (m) | U paroi moyenne (max) (m) | v nodale max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `base_tip16` (p = 1e10) | **+2,611e6** | 2,708e6 | 9,983e5 | 9,127e6 | −1,328e6 (6,6 %) | 11 536 / 13 117 | 19,56 / 29,77 | 17,27 / 16,16 | 0,160 / 0,222 / 0,160 | 0,130 (0,178) | 25,85 m/s |
| `red_tip16_pot01` (p = 1e9) | **−8,017e5** | 1,469e6 | 8,602e5 (−13,8 %) | 7,329e6 (−19,7 %) | −1,305e6 (7,1 %) | 11 302 / 13 455 (−2,0 %) | 20,58 / 30,23 | 17,84 / 16,49 | 0,207 / **0,538** / 0,322 | 0,153 (0,331) | 5,64 m/s |
| `red_adapt_base` (p = 1e10) | **−5,003e5** | 6,916e5 | 3,707e5 | 2,256e6 | −1,413e6 (21,2 %) | 5 700 / 6 667 | 14,46 / 16,89 | 12,14 / 13,40 | 0,147 / 0,145 / 0,174 | 0,086 (0,154) | 3,84 m/s |
| `red_adapt_pot01` (p = 1e9) | **−8,654e5** | 8,157e5 | 5,103e5 (+37,6 %) | 3,239e6 (+43,6 %) | −1,382e6 (14,3 %) | 8 154 / 9 251 (+43,1 %) | 17,13 / 19,31 | 14,93 / 14,89 | 0,156 / 0,180 / 0,166 | 0,110 (0,159) | 4,79 m/s |

Signe : travail net du contact positif = injecté (source), négatif = dissipé.

Lecture :

1. **Avec frottement, p = 1e9 N/m supprime l'injection** : le contact de `red_tip16` passe de
   +2,6e6 J/m injectés à −8,0e5 J/m dissipés, et la vitesse nodale max tombe de 25,9 à 5,6 m/s. Le
   facteur de pointe était donc la condition qui exposait le biais (beaucoup de joints morts en
   compression), pas sa cause.
2. **`red_adapt` ne présentait déjà pas d'injection nette à p = 1e10** (−5,0e5 J/m) : son bilan est
   dominé par le frottement. Pour lui, p/10 rend le contact un peu plus dissipatif encore.
3. **Mais le faciès bouge avec p, et c'est le point à retenir pour le chapitre.** Sur `red_tip16`,
   l'EDZ et le nombre de joints rompus sont presque inchangés (R_EDZ p95 +5 %, joints −2 %), mais le
   déplacement maximal aux reins passe de 0,22 à 0,54 m et la paroi converge davantage (U moyen
   +18 %, max ×1,9) : la paroi s'enfonce plus dans le massif fragmenté, ce qui est l'effet attendu
   d'un contact plus mou (interpénétration plus grande, non mesurée). Sur `red_adapt`, p/10 fait
   rompre 43 % de joints en plus, élargit l'EDZ de 18 % (14,5 → 17,1 m) et augmente U paroi de 28 %.
   Le contact raide retenait des blocs que le contact mou laisse se déplacer.
4. **Le résidu B4 ne s'améliore pas** (6,6 → 7,1 % ; 21,2 → 14,3 %) : il reste dominé par l'erreur
   d'intégration en O(dt) vue au diagnostic 3, que p ne corrige pas à dt fixé. Le résidu de 21 % de
   `red_adapt_base` est le plus fort de la nuit et mérite d'être signalé.
5. **Conclusion pratique** : p = 1e9 N/m règle l'injection mais ne laisse pas le faciès invariant ;
   on ne peut pas l'adopter sans une étude de convergence en p (et en dt) des observables du chapitre
   (R_EDZ, U paroi, joints rompus), au moins à 3e9 et 1e9 N/m, et une mesure de l'interpénétration.
   Tant que cette étude manque, les chiffres de tunnel du chapitre sont à donner avec p et dt.

## Réserves

Maillage réduit trop grossier (hFine = 0,40 m > l_cz/2), 4 fils (non bit-identique), un seul tirage,
t = 0,25 s seulement ; vitesse max échantillonnée sur 14 trames, pas sur tous les pas.
