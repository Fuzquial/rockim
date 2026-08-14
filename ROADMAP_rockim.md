# rockim — feuille de route au 2026-08-13

*Mise à jour après le chantier « détection du contact ». Trois objectifs
directeurs, fixés par Fernando : (1) parité structurelle avec le code de
l'article (MultiFracS / Munjiza), (2) capacité à reproduire un impact
validable contre le banc de Mines Paris, (3) robustesse et vitesse.*

---

## Décision de périmètre (Fernando, 2026-08-14) : rockim = FDEM

**rockim se concentre sur `fdem` / `fdem3d`.** La modélisation continue
(lois fem3d : dpr, saksala, saksala2011, dpdfh) est **déléguée à
Abaqus + VUMAT**, qui la gère bien — c'est la production de la thèse.
Concrètement :

- `fem`/`fem3d` et `dem`/`dem3d` passent en **mode gelé** : aucune capacité
  nouvelle, mais on CONSERVE les selftests (références Fortran 8e-14 /
  4,7e-12), les repères de la suite et le pipeline de validation croisée
  rockim ↔ Abaqus (`export_abaqus.py`) — c'est de la vérification, pas du
  développement ;
- la **spec 002 (OpenMP étendu) se réduit à `FdemSolver`/`Fdem3dSolver`** —
  l'efficacité mesurée (~55 % sur 18 threads, 2026-08-14) devient la cible ;
- le **checkpoint/restart (T030-T036)** n'inventorie plus que les deux
  solveurs FDEM — le point critique du chantier rétrécit d'autant ;
- **D2 (mémoire)** se cible sur le MatState du bulk FDEM (élastique/DP —
  maigre) plutôt que sur les grosses lois continues ;
- le bulk FDEM reste conforme à l'acquis « grains élastiques (ou DP),
  tout le softening dans les cohésifs ».

---

## Fait (2026-08-11 → 13)

- [x] **Audit complet** du solveur (architecture, physique, robustesse, perf) — `AUDIT_rockim_2026-08-11.md`
- [x] **Loi de joint 3D corrigée** : dashpot bilatéral sur joint intact, écrêtage de la résultante réservé au joint rompu, borne MOOSE `cd ≤ m_eff/dt`, compteur `dampWork_` + verdict
- [x] **f(D) de Yan (`jointSoftening = yan`) + breakMode/failMode en 3D**
- [x] **Contact général 3D : détection parallélisée** (bit-identique au sériel)
- [x] **Profileur `ROCKIM_PROF` en 3D**
- [x] **`mesh = file`** — import Gmsh MSH 2.2 (triangles 2D, tets 3D) + `tools/make_unstructured_mesh.py` + `meshes/box3d_h45.msh` + `configs/fdem3d_percussion_base.cfg`
- [x] **Avertissement au démarrage** si maillage structuré + scénario de fissuration
- [x] **`tools/verify_suite.py`** — runner de non-régression (fast/full/all, charge nulle, dampWork ≤ 0), **30 repères** (dont percussion 2D/3D, SHPB multi-corps et UCS cisaillante)
- [x] **`DOCUMENTATION_rockim.md`** — référence exhaustive des clés et des sorties, extraite du code
- [x] **`bulkViscosity`** (terme 2µD, éq. 6 de Yan) + preset `configs_yan/article_exact_base.cfg` + `PLAN_article_exact.md`
- [x] **Dépôt git** (branches `main` / `article-exact`) livré en bundle
- [x] **Optimisation de la détection du contact 2D+3D** — géométrie hoistée, sphère de rejet exacte, dédoublonnage O(1) : **×2,45 sur le contact, ×1,94 sur le mur, bit-identique**
- [x] **A2 — décharge en cisaillement sur la sécante à l'origine (éq. 18)**, 2D **et** 3D : clé `jointShearUnload = plastic | origin`, défaut bit-identique (17/17), s_p de Munjiza, `J.slip` réinterprété en origine figée (continuité à l'insertion conservée), avertissement au démarrage sur la réversibilité du frottement, **6 nouveaux repères** dont le premier repère du parc piloté par le cisaillement (UCS fig. 17)
- [x] **A1 — activation adaptative du contact (Fukuda)**, 2D **et** 3D : clé `gcActivation = full | adaptive`, défaut bit-identique. Règles C/A/B (peau endommagée + anneau, autre corps par union-find, voisinage d'un contact porté), cadence par v_max, cache des faces mortes au timing du mode full. **Percussion 3D ×2,32 bit-identique** (1130 → 488 s, 4 % des faces activées), percussion 2D bit-identique phase débris comprise, UCS −15 % aux mêmes chiffres, SHPB multi-corps armé dès t = 0 et identique sur 83 % du run (au-delà : enveloppe chaotique, le full sous OMP=2 diverge 8× plus tôt). **7 nouveaux repères**, zéro trou d'activation (instrumentation RKM_GCLOG : aucune face touchée avant activation)
- [x] **A3 — contact par POTENTIEL de Munjiza (éq. 2-5), 2D ET 3D** : clé `contact = penalty | potential`, cœur géométrique pur (`PotentialContact.hpp` : triangle-triangle en 2D, **polyèdre tet-tet** en 3D), intégrale de bord **exacte**, 3e loi machine, **conservatif** (selftests : ΔKE/KE₀ = 3,7e-12 en 2D, 2,0e-8 en 3D, transferts élastiques exacts ; SHPB incassable : gcWork −2,6 J/m pour 766 en jeu), frottement incrémental éq. 4-5 (vectoriel en 3D), détection élément-élément O(N), relève de naissance par **aire/volume** (le pen0_ du potentiel — la rampe temporelle injecte, mesuré). Gardes 3D nées du contrôle zeroload : plancher de volume + fermeture du polyèdre (les tets tangents faisaient des kN parasites au repos). SHPB : onde incidente = penalty à 3e-6 près. Zone broyée conservative (percussion 2D : e 0,55 → 0,71). Combiner avec gcActivation (potentiel+adaptatif 3D : 112 s vs 682 s). **8 nouveaux repères.** LE chantier structurel du groupe A est clos — le groupe A est **terminé**

---

## A. Parité structurelle avec Munjiza / l'article — branche `article-exact`

| # | chantier | effort | note |
|---|---|---|---|
| ~~A1~~ | ~~Activation adaptative du contact (Fukuda)~~ | **FAIT 08-13** | `gcActivation = adaptive`, 2D+3D. Percussion 3D **×2,32 bit-identique** (1130 → 488 s). Le SHPB multi-corps est armé dès t = 0 par la règle « autre corps » |
| ~~A2~~ | ~~Décharge en cisaillement sur sécante à l'origine (éq. 18)~~ | **FAIT 08-13** | `jointShearUnload = origin`, 2D+3D. UCS 51,07 → 50,37 MPa (−1,4 %), part de cisaillement 47 → 51 % |
| ~~A3~~ | ~~Contact par POTENTIEL de Munjiza (éq. 2-5)~~ | **FAIT 08-13 (2D + 3D)** | `contact = potential`, conservatif (ΔKE/KE₀ = 3,7e-12 / 2,0e-8), détection O(N), relève par aire/volume, gardes anti-slivers 3D. À combiner avec `gcActivation = adaptive` |
| A4 | Quadrature de Gauss des joints (vs 2-3 points nœud-à-nœud) | moyen | écart déclaré dans les en-têtes |
| A5 | z(D) de Munjiza | court | optionnel : f(D) de Yan couvre le besoin |

## B. Reproduire un impact validable (banc Yang / Aising, IJRMMS 2025)

| # | chantier | effort | note |
|---|---|---|---|
| **B1** | ~~Physical groups dans `mesh = file`~~ **FAIT (V1, 2026-08-14)** | — | groupes = corps (zéro joint inter-groupes), phase homonyme / `groupPhase.<nom>`, `groupVel.<nom>`, `trackGroup`, `toolShape = none` (impacteur maillé), résumé par corps, générateur `bench1` (bloc + insert sphérique), contrôle `zeroload_bench1_3d` (0 casse, gcWork = 0 exact) + impact insert→roche |
| **B2** | ~~Jauge de contrainte dans le taillant~~ **FAIT (V2, 2026-08-14)** | — | `grpFx/y/z` (force de contact nette sur le corps suivi, les deux lois) + `grpSzz` (σzz moyen volumique du corps) dans history ; la version « tranche à position donnée » attendra le banc V3 |
| **B3** | ~~Métriques de cratère~~ **FAIT (V2, 2026-08-14)** | — | `tools/crater_metrics.py` : R_crater/R_max/profondeur/aire cassée/volumes endommagé-détaché/fissures radiales + bras endommagés par secteur, multi-corps géré ; validé sur bench1 (cône hertzien R 4,6 mm) et percussion longue (R 8,5 mm, étoile 12 bras à 18 mm) |
| **B4** | ~~Bilan d'énergie par sous-système~~ **FAIT (V2, 2026-08-14)** | — | travaux par famille (éléments, joints, dashpot, contact+frottement, Cundall, Lysmer, outil, platines) + correction leapfrog exacte f²dt²/2m ; résidu : 0,017 % (percussion 2D), 0,005 % (3D), zéro machine au repos ; bloc `energy budget` au résumé, 6 colonnes history, extracteur `budget` dans la suite ; 2D et 3D, bit-neutre (12/12) |
| B5 | ~~Brossage des fragments~~ **FAIT (V2, 2026-08-14)** | — | `crater_metrics.py --brush beta` : volume brossable = β × détachés dans le bol (percussion longue : 96,4 mm³ sur 120,5 détachés) ; masse = ρ × volume |
| B6 | DIF (effet de vitesse) sur les joints | court | le taux de déformation est déjà calculé depuis `bulkViscosity` |
| **B7** | **Récupérer les données expérimentales d'Aising** | — | pas du code, mais bloquant pour toute validation |

## D. Hydro (OUVERT le 2026-08-14, à lancer APRÈS la clôture de la partie
## technique — décision Fernando)

Couplage hydro-mécanique 2D d'abord : réseau d'écoulement sur le graphe des
joints (loi cubique, conductivité ∝ ouverture³ — l'ouverture est déjà
mesurable), pression appliquée aux lèvres, couplage explicite décalé.
Cible de validation : initiation de fracture depuis un forage pressurisé
(`confineFaces = bore` = le cas limite à fluide infiniment mobile).
Prérequis ABSOLUS : fiabilité 3D close (E0-E3), B4 avec travail du
confinement, checkpoint/restart. Référence marché : module hydro d'Irazu.

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
| **D0** | **Potentiel 3D en régime débris** (N1, diagnostic aux compteurs 2026-08-14 ; fait : seaux réutilisés, ordre canonique, SAT complet caché, clip sans copie → 682→448 s à T = 5e-5, tout bit-neutre ; longue 3 474 s vs 488 pénalité). Le poste dominant de la longue est l'**intégration exacte des 33 M de clips avec force** (~70 µs pièce, ≈ 2/3 du run : la subdivision aux 12 plans sur les contacts persistants des débris) ; clips vides ~12 % ; les axes d'arêtes du SAT n'y séparent plus rien | moyen-gros | pistes : intégration allégée en contact persistant (quadrature vs exactitude — à arbitrer), warm-start du polyèdre, cadence des contacts stationnaires, chapeau sans atan2/tri (ordre-équivalence à prouver), plancher PHYSIQUE `potFloor` opt-in |
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

*Le groupe A (parité structurelle avec l'article) est terminé : A1 + A2 + A3 faits le 08-13, A4/A5 restent optionnels.*

1. **B4 — bilan d'énergie par sous-système** : zéro risque, et il est devenu URGENT — A3 a montré que la loi de contact pilote un tiers du bilan d'énergie de l'impact ; il faut pouvoir attribuer chaque joule (fissuration, frottement, Lysmer, Cundall, contact) pour arbitrer penalty/potential contre l'expérience.
2. ~~B1 — physical groups~~ **FAIT (V1, 2026-08-14)** : l'insert maillé impacte l'éprouvette (deux corps Gmsh, zeroload propre) — la suite est V2/V3 du plan v2 (bilan d'énergie, cratère, piston-taillant-roche).
3. **E2 — recalibration** : contrepartie obligatoire d'A3 (les calibrations historiques compensaient le puits d'énergie du contact quasi-plastique).
4. **D0 — chemin clip-vide du potentiel 3D** : le diagnostic est posé (compteurs au résumé) ; fait : seaux, ordre canonique, SAT complet, clip sans copie (682→448 s à T = 5e-5, bit-neutre) ; reste : la percussion longue à ~7× la pénalité — pistes au tableau D et au plan v2.
