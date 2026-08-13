# rockim — feuille de route au 2026-08-13

*Mise à jour après le chantier « détection du contact ». Trois objectifs
directeurs, fixés par Fernando : (1) parité structurelle avec le code de
l'article (MultiFracS / Munjiza), (2) capacité à reproduire un impact
validable contre le banc de Mines Paris, (3) robustesse et vitesse.*

---

## Fait (2026-08-11 → 13)

- [x] **Audit complet** du solveur (architecture, physique, robustesse, perf) — `AUDIT_rockim_2026-08-11.md`
- [x] **Loi de joint 3D corrigée** : dashpot bilatéral sur joint intact, écrêtage de la résultante réservé au joint rompu, borne MOOSE `cd ≤ m_eff/dt`, compteur `dampWork_` + verdict
- [x] **f(D) de Yan (`jointSoftening = yan`) + breakMode/failMode en 3D**
- [x] **Contact général 3D : détection parallélisée** (bit-identique au sériel)
- [x] **Profileur `ROCKIM_PROF` en 3D**
- [x] **`mesh = file`** — import Gmsh MSH 2.2 (triangles 2D, tets 3D) + `tools/make_unstructured_mesh.py` + `meshes/box3d_h45.msh` + `configs/fdem3d_percussion_base.cfg`
- [x] **Avertissement au démarrage** si maillage structuré + scénario de fissuration
- [x] **`tools/verify_suite.py`** — runner de non-régression (fast/full/all, charge nulle, dampWork ≤ 0), **23/23 PASS**
- [x] **`DOCUMENTATION_rockim.md`** — référence exhaustive des clés et des sorties, extraite du code
- [x] **`bulkViscosity`** (terme 2µD, éq. 6 de Yan) + preset `configs_yan/article_exact_base.cfg` + `PLAN_article_exact.md`
- [x] **Dépôt git** (branches `main` / `article-exact`) livré en bundle
- [x] **Optimisation de la détection du contact 2D+3D** — géométrie hoistée, sphère de rejet exacte, dédoublonnage O(1) : **×2,45 sur le contact, ×1,94 sur le mur, bit-identique**
- [x] **A2 — décharge en cisaillement sur la sécante à l'origine (éq. 18)**, 2D **et** 3D : clé `jointShearUnload = plastic | origin`, défaut bit-identique (17/17), s_p de Munjiza, `J.slip` réinterprété en origine figée (continuité à l'insertion conservée), avertissement au démarrage sur la réversibilité du frottement, **6 nouveaux repères** dont le premier repère du parc piloté par le cisaillement (UCS fig. 17)

---

## A. Parité structurelle avec Munjiza / l'article — branche `article-exact`

| # | chantier | effort | note |
|---|---|---|---|
| **A1** | **Activation adaptative du contact** (Fukuda) — `act_` contient tout l'extérieur et le balaie à chaque pas ; il reste ~40 % du pas sans débris | court | ⚠️ pas de règle « aucun joint mort → rien » : le SHPB multi-corps a besoin du contact dès t = 0 |
| ~~A2~~ | ~~Décharge en cisaillement sur sécante à l'origine (éq. 18)~~ | **FAIT 08-13** | `jointShearUnload = origin`, 2D+3D. UCS 51,07 → 50,37 MPa (−1,4 %), part de cisaillement 47 → 51 % |
| **A3** | **Contact par POTENTIEL de Munjiza + détection NBS/MR** (éq. 2-5) | gros (1-2 j) | LE chantier structurel. Opt-in `contact = penalty \| potential`. Test décisif : collision élastique sans frottement → `gcWork ≈ 0` machine |
| A4 | Quadrature de Gauss des joints (vs 2-3 points nœud-à-nœud) | moyen | écart déclaré dans les en-têtes |
| A5 | z(D) de Munjiza | court | optionnel : f(D) de Yan couvre le besoin |

## B. Reproduire un impact validable (banc Yang / Aising, IJRMMS 2025)

| # | chantier | effort | note |
|---|---|---|---|
| **B1** | **Physical groups dans `mesh = file`** | moyen (~200 l.) | **LE verrou** : débloque d'un coup le multi-corps piston-taillant-roche, la géométrie réelle (cylindre, insert hémisphérique) et l'affectation correcte des CL |
| **B2** | Jauge de contrainte dans le taillant | court | copie des moniteurs SHPB 2D (`monEl1_`, moyenne surfacique) |
| **B3** | Métriques de cratère : rayon, longueur des fissures radiales | court | post-traitement Python sur les VTU de joints |
| **B4** | **Bilan d'énergie par sous-système** | court (~80 l.) | fissuration, frottement, Lysmer, Cundall, plasticité. Sert la validation *et* la chasse aux pathologies. Zéro risque (instrumentation pure) |
| B5 | Algorithme de « brossage » des fragments (anti-gravité, β = 0,8) | court | pour comparer la masse simulée à la masse collectée |
| B6 | DIF (effet de vitesse) sur les joints | court | le taux de déformation est déjà calculé depuis `bulkViscosity` |
| **B7** | **Récupérer les données expérimentales d'Aising** | — | pas du code, mais bloquant pour toute validation |

## C. Robustesse et dette technique (issues de l'audit)

| # | chantier | effort |
|---|---|---|
| C1 | Validation d'entrées : `PhaseSet::validate` partout, nx/ny/nz > 0, dt fini, clés inconnues signalées, `getb` strict | court |
| C2 | `st.x0` non renseigné hors `fem3d` → dpdfh/Weibull **muets** sous fdem/fdem3d | court |
| C3 | `ftScale` ignoré par `dpr`/`saksala` (champ écrit dans les VTU mais sans effet) | court |
| C4 | `pairKey` non trié en `Dem3d` → perte silencieuse de l'historique tangentiel | très court |
| C5 | Gardes NaN testant `u_[0]`, qui peut être un nœud FIXED (détecteur aveugle) | très court |
| C6 | Bornes MOOSE / plafond d'impulsion absents hors FDEM (fem, fem3d, dem, dem3d) | court |
| C7 | Contraction de **faces minces** dans `Tessellation3` (c'est le dt 3D qui paie) | moyen |

## D. Performance restante

| # | chantier | effort | gain attendu |
|---|---|---|---|
| D1 | VTK **binaire** + `frames.csv` tenu ouvert | court | dominant en E/S sur gros maillages |
| D2 | `MatState` par loi (≈ 700 o/élément quelle que soit la loi) | court | ×4-5 mémoire en élastique/dpr |
| D3 | OpenMP dans `fem`, `dem`, `dem3d` (zéro pragma aujourd'hui) | moyen | ces modules sont 100 % sériels |
| D4 | Grille DEM reconstruite tous les k pas ; plancher 1e-12 en 2D | court | |
| D5 | `toolContact` en O(N nœuds) → requête de grille locale | court | |
| D6 | `RandomField` anisotrope : rayon de troncature isotrope | court | ×100 en 3D à 10:1 |
| D7 | **Restart / checkpoint** | moyen | indispensable dès les runs de plusieurs dizaines d'heures |
| D8 | Mass / density scaling quasi-statique | moyen | aucune recette transposable trouvée (cf. veille) |

## E. Physique et capacités

| # | chantier | effort |
|---|---|---|
| **E1** | SHPB : extrémité absorbante/libre après l'impulsion (rebroyage sans fin aujourd'hui) | court |
| **E2** | **Recalibration Red Bohus de zéro** (jeu invalidé par la correction d'amortissement) | gros (campagne) |
| E3 | Validation « rebond complet » T = 4e-4 sur maillage non structuré + ajout au runner | court |
| E4 | Re-mesure des démos percussion 3D (les références ont changé) | court |
| E5 | Parité 3D : platines / UCS / brésilien en 3D | moyen |
| E6 | Parité GBM 2D→3D en adaptatif (~15-18 h de calcul) | calcul |
| E7 | Impacts répétés (pour une vitesse de pénétration) | moyen |
| E8 | Tets composites (verrouillage volumique de la zone broyée) | gros |

## F. Logistique

| # | chantier |
|---|---|
| F1 | Cloner le bundle **hors OneDrive**, rebuild MSVC, re-baseliner `verify_suite` côté Windows |
| F2 | Sous-dossier dédié aux livrables (ne plus encombrer la racine `rockim/`) |
| F3 | Remote git (URL + token) pour un vrai `push` |

---

## Prochains pas recommandés

1. **A1 — activation adaptative du contact** : le plus court des chantiers restants à fort rendement (~40 % du pas passé à balayer un extérieur sans débris).
2. **B4 — bilan d'énergie par sous-système** : zéro risque, sert deux objectifs à la fois — et c'est lui qui rendra *mesurable* la réversibilité du frottement signalée par A2.
3. **B1 — physical groups** : le verrou qui ouvre la reproduction du banc de ton labo.
4. **A3 — contact par potentiel + NBS** : le gros morceau, à faire quand une séance entière peut y être consacrée.
