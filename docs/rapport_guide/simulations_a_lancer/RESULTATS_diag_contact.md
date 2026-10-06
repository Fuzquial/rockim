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
5. **Le résidu B4 est le même dans les cinq calculs (−1,30 à −1,34e6 J/m)**, alors que les autres
   termes varient d'un facteur 2 à 5. C'est un terme fixe non compté par le bilan (piste : travail
   de la pré-contrainte in situ ou du relâchement d'excavation), indépendant du contact. À
   départager avant de lire le bilan comme un diagnostic fin.
6. Piste suivante suggérée : la partie normale du potentiel au relais joint → contact (65 % des
   joints meurent en compression ; charge normale lâchée au relais 1,69e6 kN/m cumulés dans la
   base). Un calcul `contactMu = 0` + `gcBirth = penalty` séparerait naissance et potentiel normal.

## Réserves

Maillage réduit trop grossier (hFine = 0,40 m > l_cz/2), 4 fils (non bit-identique), un seul tirage,
t = 0,25 s seulement ; vitesse max échantillonnée sur 14 trames, pas sur tous les pas.
