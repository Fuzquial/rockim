# Port de la branche `joint-handoff` dans `g0` — verdict ligne par ligne

**Date** : 2026-09-05 · **Branche source** : `joint-handoff` (e0af1a4) ·
**Cible** : `g0` (nee de `f2-2026-09-02`, commit de naissance `1a2b2a9`) ·
**Perimetre** : `src/FdemSolver.cpp`, `src/Fdem3dSolver.cpp`,
`include/rockim/FdemSolver.hpp`, `include/rockim/Fdem3dSolver.hpp`.

**VERDICT GLOBAL — AUCUNE LIGNE DE `src/` NI DE `include/` N'A ETE PORTEE.**
Les trois « pertes probables » signalees a la revue (le poste d'energie separe
`brushWork_`, le court-circuit `if (muCRes_ < 0.0) return muC_;` de
`contactResidualMu`, le budget de pas de temps Signorini A1 qui sort `kp_` du
CFL) sont **toutes les trois deja dans `g0`** : `f2` les a reimplementees, et
dans les trois cas en les ETENDANT. Le reste de l'ecart est du refactor : soit
du code de `main` que `f2` a lui-meme reecrit — et que `joint-handoff`, parti de
`main`, a simplement conserve —, soit la meme intention exprimee autrement.

---

## 1. Comment l'ecart a ete mesure

`joint-handoff` est parti de `main` (139df12) et n'a jamais ete fusionne ; `f2`
est parti du meme `main` et a diverge de ~5 400 lignes. Comparer les deux
branches directement produit un diff illisible. On a donc mesure, **fichier par
fichier**, l'ensemble des lignes presentes dans la version `joint-handoff` et
absentes de la version `g0` (comparaison sur la ligne DEPOUILLEE de son
indentation, donc insensible au reformatage), puis on a marque chaque ligne
selon qu'elle existait deja dans `main` ou qu'elle est **NEUVE** (ecrite par
`joint-handoff`) :

| fichier | lignes de `joint-handoff` absentes de `g0` | dont NEUVES |
|---|---:|---:|
| `src/FdemSolver.cpp` | 71 | 3 |
| `src/Fdem3dSolver.cpp` | 28 | 0 |
| `include/rockim/FdemSolver.hpp` | 3 | 2 |
| `include/rockim/Fdem3dSolver.hpp` | 4 | 3 |
| **total** | **106** | **8** |

Le chiffre de « 87 lignes » de la revue est de cette famille : la valeur exacte
depend de la normalisation (106 lignes brutes ; 8 lignes reellement ecrites par
la branche ; 5 si l'on retire les commentaires). **Les 98 lignes marquees
`main` ne sont pas du contenu de `joint-handoff`** : ce sont des lignes de
`main` que `f2` a reecrites pour son propre compte et que la branche, plus
ancienne, a conservees intactes. Elles ne peuvent pas etre « perdues » par
`g0` : `g0` a retenu la version plus recente.

Ci-dessous, le verdict est rendu **par site** (bloc contigu), avec la ligne
concernee dans `joint-handoff`, ce qu'elle fait, ce qu'il y a a la place dans
`g0`, et la decision : **REFACTOR** (ne rien porter) ou **PERTE REELLE**
(porter). Les numeros de ligne de `g0` sont ceux de l'etat au moment de la
lecture ; ils peuvent glisser, les ancres textuelles restent valables.

---

## 2. Les trois « pertes probables » de la revue — les trois sont infirmees

### 2.1 `brushWork_`, le poste d'energie separe du balai de tri

| | |
|---|---|
| `joint-handoff` | `src/FdemSolver.cpp:3559` `brushWork_ += bw;  // POSTE SEPARE, jamais dans sumW` ; `src/Fdem3dSolver.cpp:3890` idem ; `include/rockim/FdemSolver.hpp:970` et `Fdem3dSolver.hpp:425` `double brushWork_ = 0.0;` |
| `g0` | `src/FdemSolver.cpp:4715` `brushWork_ += bw;  // poste SEPARE (cf. energyBodyForces)` ; `src/Fdem3dSolver.cpp:4370` idem ; membres declares en `FdemSolver.hpp` (~1118) et `Fdem3dSolver.hpp:453` |
| **verdict** | **REFACTOR.** Le membre, l'accumulation et l'affichage separe existent a l'identique. `f2` va plus loin : la cle `energyBodyForces` (`eBody_`) permet de FAIRE ENTRER `gravWork_` et `brushWork_` dans `sumW` quand on le demande (`FdemSolver.cpp:6627-6628`, `8206`, `8218` ; `Fdem3dSolver.cpp:2584-2585`, `4779`, `4790`), et le message de fin le dit (`(DANS le bilan B4 : energyBodyForces = on)` / `(hors bilan B4)`). Seul le libelle du commentaire differe. **Rien a porter.** |

### 2.2 Le court-circuit `if (muCRes_ < 0.0) return muC_;` de `contactResidualMu`

| | |
|---|---|
| `joint-handoff` | `include/rockim/FdemSolver.hpp:1051` et `1054`, `include/rockim/Fdem3dSolver.hpp:498` et `501` : `if (muCRes_ < 0.0) return muC_;` puis `if (!p) return muC_;` — 4 des 8 lignes NEUVES du diff |
| `g0` | memes deux tests, mais rendant `mu` et non `muC_` : `include/rockim/FdemSolver.hpp` (~1292-1296) et le corps jumeau de `include/rockim/Fdem3dSolver.hpp` |
| **verdict** | **REFACTOR, et le changement est OBLIGATOIRE.** `f2` a ajoute au meme `ctcMu()` le frottement PAR PHASE (`muPerPhase_`, Table 1 de Yang 2026) et le couplage 1-D (`cplMode_`, Solidity `Y3Did.c` l. 995 / 1044 / 1263-1265 / 1292). La variable locale `mu` part de `muC_` puis est eventuellement remplacee par le `muPhase_` le plus FAIBLE de la paire. Rendre `muC_` au lieu de `mu` **annulerait le frottement par phase** dans le chemin non-pulverise : ce serait une regression, pas un port. Quand aucune des deux capacites n'est armee, `mu == muC_` et les deux ecritures coincident. **Rien a porter.** |

### 2.3 Le budget de pas de temps Signorini A1 (sortir `kp_` du CFL)

| | |
|---|---|
| `joint-handoff` | `src/Fdem3dSolver.cpp:2033-2034` : `dtMin = std::min(dtMin, 2.0*sqrt(m_[i] / (K[i] + (toolSig_ ? 0.0 : nExtra*kp_))));` ; en 2D, le meme choix ecrit en ligne |
| `g0` | `src/Fdem3dSolver.cpp:2348` `double kContact = toolSig_ ? 0.0 : kp_;` sous le commentaire « A1 : en Signorini l'outil n'a plus de raideur (condition sur la VITESSE) — kp_ sort du budget » ; `src/FdemSolver.cpp:4151` `double kContact = toolSig_ ? kpPlaten_ : std::max(kp_, kpPlaten_);` sous le meme commentaire A1 |
| **verdict** | **REFACTOR, et `f2` est PLUS COMPLET.** La sortie de `kp_` du budget sous `toolSignorini` est presente dans les DEUX solveurs. `f2` a de plus : la raideur des PLATINES (`kpPlaten_`) qui reste au budget en 2D, celle du contact general en SHPB, celle du potentiel (`potP_`, `potKt_`), la cle `dtBudgetTangential` avec son avertissement quand `potKt_` depasse le budget sans y entrer (Xiang, Munjiza, Latham & Guises, *Eng. Comput.* 26(6) 2009, p. 677), la borne diffusive de `bulkViscosity` et celle de l'amortisseur nodal. **Rien a porter.** |

---

## 3. Les 4 autres lignes NEUVES de `joint-handoff`

| ligne (`joint-handoff`) | contenu | ce qu'il y a dans `g0` | verdict |
|---|---|---|---|
| `src/FdemSolver.cpp:5175` | `double cap = ctcMu(elemOf_[i]) * rn;   // WP6` — le cap de Coulomb sur l'IMPULSION du contact de Signorini outil-roche, evalue au mu residuel | l'ALGEBRE impulsion / saut de vitesse a ete sortie vers l'en-tete partage `include/rockim/ToolSignorini.hpp` (`toolsig::impulse`, qui applique `cap = mu * rn`) ; la lambda `nodeSig` elle-meme reste dans chaque solveur (2D `src/FdemSolver.cpp:6462`, 3D `src/Fdem3dSolver.cpp:4139`) et lui passe deja `ctcMu(elemOf_[i])` : `src/FdemSolver.cpp:6400` et `6464` (`m_[i], dt_, ctcMu(elemOf_[i]), toolSigRelax_);  // WP6 sur mu`) | **REFACTOR** — meme physique, code partage 2D/3D |
| `src/FdemSolver.cpp:5263-5264` | commentaire `// WP6 : volontairement PAS de ctcMu ici — le plateau est une frontiere de machine, pas un support de fragments.` | `src/FdemSolver.cpp:6645-6648`, meme commentaire ETENDU : `// WP6/WP7 : … Il garde donc contactMu GLOBAL, sans frottement par phase ni couplage (1-D).` | **REFACTOR** — la decision de perimetre est conservee et precisee |
| `include/rockim/Fdem3dSolver.hpp:496` | commentaire `// evaluation residuelle (atomic : appele depuis les boucles OMP).` | le `#pragma omp atomic` et les compteurs `nCtcPulv_` / `tCtcPulv0_` sont en place dans le `ctcMu()` de `g0` ; seul le libelle du commentaire a change | **REFACTOR** |

---

## 4. Les 98 lignes de `main` que `f2` a reecrites (aucune n'appartient a `joint-handoff`)

Regroupees par site. Toutes : **REFACTOR, rien a porter.**

| site (`joint-handoff`) | ce que la ligne faisait | ce que `g0` fait a la place |
|---|---|---|
| `FdemSolver.cpp:1149-1152` | `mesh = file` refuse pour `disc` ET `shpb` | `FdemSolver.cpp:1421-1425` : seul `shpb` est refuse ; `disc + mesh = file` est desormais SUPPORTE (message dedie) |
| `FdemSolver.cpp:1808` | `"meshFile: element references unknown node id"` | `FdemSolver.cpp:2085-2087` (3D : `1218`) : erreur NOMMEE C3 — numero de l'element ET id du noeud fautif |
| `FdemSolver.cpp:1869` | `Tessellation::build(..., gm == "delaunay", gh)` | `FdemSolver.cpp:2173-2177` : trois arguments de plus (`gSpread`, `pSize`, `gRandom` — polydispersite Laguerre-Newton et `grainMeshRandom`) |
| `FdemSolver.cpp:1925` / `Fdem3dSolver.cpp:1286` | `throw ... "inverted element in mesh gen"` / `"degenerate tet"` | `FdemSolver.cpp:2270` et `Fdem3dSolver.cpp:1499` : `guards::degenerateError(...)` — erreur NOMMEE C3 (mode, id, sommets, positions, aire ou volume) |
| `FdemSolver.cpp:2711` | `} else if (scen_ == Scenario::PERCUSSION) {` (fond encastre) | `FdemSolver.cpp:3808` : `} else if (scen_ == Scenario::PERCUSSION \|\| fixBottomShear_) {` |
| `FdemSolver.cpp:3056` | `K[i] = 2.0 * phases_.mat[...].E * thk_;` | `FdemSolver.cpp:4131-4133` : `(tiOn_ ? tiEmax_ : phases_.mat[...].E)` — isotropie transverse |
| `FdemSolver.cpp:3099` | `double cfl = hCfl / phases_.maxCp();` | `FdemSolver.cpp:4250` : `hCfl / (tiOn_ ? tiCp_ : phases_.maxCp())` — celerite de Christoffel |
| `FdemSolver.cpp:3119` | `dt_ = dtFactor * min(min(dtMin, cfl), dtVis);` | `FdemSolver.cpp:4283 s.` : le minimum porte en plus sur la borne de l'amortisseur nodal et les autres bornes ajoutees depuis |
| `FdemSolver.cpp:3196` | verrou du pic a 30 % (`nBroken_ > 0 && …`) | `FdemSolver.cpp:4403-4406` : meme verrou, condition reecrite |
| `FdemSolver.cpp:3276-3293` (18 l.) et `Fdem3dSolver.cpp:2128-2141` (14 l.) | garde E5 : ~256 noeuds echantillonnes tous les 1024 pas, decalage tournant | **remplacee par C4 (w20)** : `checkFinite()` teste TOUTES les composantes de `u`, `v`, `f` de TOUS les noeuds tous les `nanCheckEvery` pas (defaut 256) et NOMME le premier fautif — `FdemSolver.cpp:4412-4416`, `Fdem3dSolver.cpp:2511`, `include/rockim/Guards.hpp`. Le detecteur de `f2` est strictement plus fort |
| `FdemSolver.cpp:3364` | `double szz = nuP_[e.phase] * (s(0) + s(1));` | `FdemSolver.cpp:4579` : `tiOn_ ? tiZz_.dot(eps) : nuP_[…] * (s(0)+s(1))` |
| `FdemSolver.cpp:3535`, `3540` ; `Fdem3dSolver.cpp:3868`, `3873` | boucle OpenMP et report nodal de `bodyForces()` | `FdemSolver.cpp:4681 s.`, `Fdem3dSolver.cpp:4393 s.` : meme boucle, comptabilite V2/B4 du travail de la pesanteur ajoutee |
| `FdemSolver.cpp:3683` ; `Fdem3dSolver.cpp:3982` | `" J/m (hors bilan B4)\n"` | `FdemSolver.cpp:4913-4915`, `Fdem3dSolver.cpp:4518 s.` : message conditionnel `eBody_` (cf. §2.1) |
| `FdemSolver.cpp:4661` ; `Fdem3dSolver.cpp:3190` | `if (muC_ > 0.0 && potKt_ > 0.0)` dans le contact par potentiel | meme test, `muC_` remplace par `ctcMu(...)` (frottement par phase / couplage) |
| `FdemSolver.cpp:5165-5186` (corps de `nodeSig`) | contact de Signorini outil-roche ecrit dans le solveur 2D | son ALGEBRE est **sortie** vers `include/rockim/ToolSignorini.hpp`, partagee 2D/3D ; la GEOMETRIE et la lambda restent dans chaque solveur, a dessein (cf. §3) |
| `FdemSolver.cpp:5673` | `double vol = -0.5 * A2 * thk_;` (fermeture de cavite hydro) | `FdemSolver.cpp:7147` : `-0.5 * (A2 + aClose) * thk_` — terme de fermeture ajoute |
| `FdemSolver.cpp:5990-5995`, `6091-6095` | amortissement local de Cundall (branche groupe et branche noeud) | meme arithmetique, pilotee par `dampNow` (bascule `dampingSwitchT`) au lieu de `damping_` |
| `FdemSolver.cpp:6202-6220` (writeFrame) | deux appels `writeTriMesh` dupliques, `bdOn_` ou non | `FdemSolver.cpp:7815 s.` : une carte NOMMEE `vtk::ScalarField ef{...}` que les champs optionnels enrichissent — byte-identique quand aucun optionnel n'est arme |
| `FdemSolver.cpp:6224`, `6239` | champs du `.vtu` des joints | `g0` ecrit en plus `tInsert` (catalogue AE de nucleation) et `weakPlane` |
| `FdemSolver.cpp:6328-6329` | corps de `countInserted()` | `FdemSolver.cpp:8069 s.` : meme fonction, exclusion des joints pre-rompus (`J.pre`) ajoutee |
| `Fdem3dSolver.cpp:1862-1863` | `toolShape must be sphere \| flat \| none (3D)` | `Fdem3dSolver.cpp:2137` : `sphere \| flat \| pdc \| none` — l'outil PDC 3D a ete ajoute |
| `Fdem3dSolver.cpp:2122` | `peakF_ = std::max(peakF_, tool_.F.norm());` | meme mesure, reecrite avec la metrologie d'outil de `f2` |
| `Fdem3dSolver.cpp:3565`, `3630` | `if (tool_.flat) {` | `Fdem3dSolver.cpp:4076`, `4153` : `} else if (tool_.flat) {` — une branche `pdc` a ete inseree avant |

---

## 4bis. Contre-verification adverse (2026-09-06)

La mesure de l'ecart a ete refaite independamment (lignes depouillees de leur
indentation, dedupliquees, marquees contre `main`) : **120 lignes** de
`joint-handoff` absentes de `g0`, dont **9 neuves**. Le tableau du paragraphe 1 en
annonce 106 et 8 ; l'ecart tient a la normalisation, mais l'ordre de grandeur et
surtout la CONCLUSION sont reproduits. La neuvieme ligne est
`double k = kPara * J.pj * J.A0 / 3.0;` (`src/Fdem3dSolver.cpp`) : elle n'apparait
« absente » que parce que le port de `insertion = none` l'a elle-meme reecrite en
`double k = noJoints_ ? 0.0 : kPara * J.pj * J.A0 / 3.0;`. Rien n'est perdu.

Dix verdicts ont ete recontroles dans le code de `g0`, tous confirmes : `brushWork_`
et `energyBodyForces` ; `ctcMu()` rendant `mu` sous `muPerPhase_` / `cplMode_` ;
`kContact` (budget A1) dans les DEUX solveurs ; `disc + mesh = file` desormais
supporte, seul `shpb` etant refuse ; `guards::degenerateError` ; `fixBottomShear_`
dans la branche percussion ; `checkFinite` / `nanCheckEvery` (C4) a la place de
l'echantillon E5 ; `toolShape = sphere | flat | pdc | none` en 3D ;
`countInserted()` excluant `J.pre` ; l'extraction Signorini (avec la reserve de
redaction corrigee aux paragraphes 3 et 4).

## 5. Consequences

- **Aucune modification de `src/` ni de `include/` n'est faite au titre de
  `joint-handoff`.** L'ancre de bit-identite n'est donc pas concernee par ce
  lot ; la mesure `8/8 IDENTIQUE` qui l'encadre vaut avant comme apres.
- Les seuls fichiers repris de la LIGNEE `joint-handoff` sont les deux scripts
  Python que sa descendante `dif-intrinseque` a ajoutes
  (`bench_impact/tools/fig_bilan.py`, et la coquille de chaine brute corrigee
  dans `bench_impact/tools/fig_fp.py`) — voir le commit dedie.
- `joint-handoff` reste au depot (branche + worktree `rockim_p4`) : rien n'est
  supprime (decision 2 de Fernando). Ce document est la trace de la lecture,
  pour qu'elle n'ait pas a etre refaite.
