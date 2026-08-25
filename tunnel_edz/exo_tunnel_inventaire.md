# Inventaire « exo tunnel / exo hole » (banc 6) — DP-DFH sur cavité pressurisée

Racine : `C:\Users\fuzquianoalricabi\simulations\CONTINUUM\exo_hole_plate\`
Sous-dossier campagne confinée : `...\exo_hole_plate\confine_lc5\`
Copie archivée (drive/base) : `C:\Users\fuzquianoalricabi\OneDrive - Université Paris Sciences et Lettres\Documents\phd_geothermie\CONTINUUM\OBJECTIVITE_MAILLAGE\03_exo_hole_2d\` + `SYNTHESE.md`

Deux synthèses écrites font foi :
- `C:\Users\fuzquianoalricabi\simulations\CONTINUUM\exo_hole_plate\RESULTATS_bench6.md` (banc 6)
- `C:\Users\fuzquianoalricabi\simulations\CONTINUUM\exo_hole_plate\confine_lc5\RESULTATS_exo.md` (campagne confinée)
- Fiche projet : `C:\Users\fuzquianoalricabi\simulations\phd\CONTINUUM.md` §4, §5bis, §7 ; `PHD.md` §41

⚠️ Distinction déjà documentée (CONTINUUM.md §5bis) : **« l'exo hole » par défaut = le banc 6** = decks `TUN_*` + `vumat_hole.f`. Les decks `EXC_*` + `vumat_perc.f` (dans `confine_lc5\`) = **uniquement** la campagne confinée. Les deux VUMAT diffèrent d'**une seule ligne** (vérifié par `diff` : ligne 189, `DELD=0.98` vs `DELD=1.0d9`).

---

## 1. Les decks Abaqus `TUN_*`

### 1.1 Géométrie et modèle (identiques pour tous les decks)

Générateurs : `tunnel_sweep.py`, `tunnel_mesh.py`, `tunnel_lc.py`, `tunnel_lc2.py`, `tunnel_tri.py`, `tunnel_nc3.py`, `tunnel_nc5tri.py`, `tunnel_nc10tri.py` (tous dans la racine ci-dessus).

| Paramètre | Valeur | Source |
|---|---|---|
| Domaine | plaque carrée 200 × 200 mm (`HALF = 100`) | `tunnel_sweep.py` l.14 |
| Cavité | trou circulaire centré, **rayon a = 10 mm** (`RHOLE`) | idem |
| Épaisseur | 1 mm (`HomogeneousSolidSection … thickness=1.0`) | idem |
| Formulation | 2D **déformation plane**, `TWO_D_PLANAR`, `nlgeom=YES` | idem + `TUN_t003_lys.inp` |
| Éléments | **CPE4R** (+ CPE3 de remplissage) pour la série quad ; **CPE3 seuls, `elemShape=TRI`** pour la série non structurée | `tunnel_sweep.py` / `tunnel_tri.py` |
| Densité | 2,62e-9 t/mm³ | idem |
| `*Depvar` | 16, `delete=1` | idem |
| Pas de temps | `*Dynamic, Explicit`, **pas de mass scaling** | `TUN_t003_lys.inp` l.98582 |
| **Bulk viscosity** | **`0., 0.` — DÉSACTIVÉE** (pour ne pas polluer la signature de vitesse) | patch de `add_lysmer_any.py` |
| Sorties | `*Output, field, number interval=30` ; `S, SDV, LE, STATUS` ; `U, V` ; history `PRESELECT` (énergies) | patch de `add_lysmer_any.py` |

### 1.2 Chargement

- Pression interne **uniforme sur la paroi du trou**, `PMAX = 250 MPa`, `*Dsload, amplitude=Ramp / HOLE, P, 250.`
- Rampe **SmoothStep** de 0 → 250 MPa sur `TSTEP`, puis **maintien** à 250 MPa.
- Durée du step = `max(TSTEP, THOLD)` avec `THOLD = 40 µs` (temps de propagation pour laisser les fissures radiales se former jusqu'à Rd).
- Taux balayés (`tunnel_sweep.py`, `TSTEPS = [300, 30, 10, 3, 1] µs`) :

| deck | TSTEP | ṗ = 250/TSTEP (MPa/µs) | ε̇ paroi ≈ ṗ/E (/s) |
|---|---|---|---|
| `TUN_t300_lys` | 300 µs | 0,833 | ~16 |
| `TUN_t030_lys` | 30 µs | 8,33 | ~160 |
| `TUN_t010_lys` | 10 µs | 25 | ~480 |
| `TUN_t003_lys` | 3 µs | 83,3 | ~1600 |
| `TUN_t001_lys` | 1 µs | 250 | ~4800 |

(ε̇ = valeurs annoncées dans `RESULTATS_bench6.md` §1 ; le générateur imprime `pdot*1e6/52000`.)

### 1.3 Conditions aux limites — **frontières absorbantes de Lysmer**

Script : `C:\Users\fuzquianoalricabi\simulations\CONTINUUM\exo_hole_plate\add_lysmer_any.py`
(usage : `abaqus python add_lysmer_any.py TUN_t003` → `TUN_t003_lys.inp` ; suffixe `_lys` = deck effectivement calculé).

- Chaque nœud des 4 bords est relié à un nœud « sol » coïncident par un élément **`CONN2D2`** (Abaqus 2D Explicit n'a **pas** de `DASHPOTA` — cf. mémoire « Limites Abaqus 2D Explicit »).
- Amortisseurs `*Connector Damping` : `c_p = ρ·c_P·L_trib·t`, `c_s = ρ·c_S·L_trib·t`, avec `L_trib = 8·HALF / n_bord` **adaptatif au maillage**.
- Trois familles : `CB_X` (bords ±x : normal = c_p, tangent = c_s), `CB_Y` (bords ±y : inversé), `CB_C` (coins : c_p sur les deux axes).
- Valeurs pour `TUN_t003_lys.inp` : c_p = 0,0127862, c_s = 0,00738214.
- `*Boundary` : `LYS_GND, 1, 6` (sols encastrés) + `LYS_BND, 6, 6`.
- **Aucun encastrement mécanique du domaine** — pas de symétrie, plaque entière (360°).

### 1.4 Le catalogue exact des decks (avec comptes d'éléments recomptés dans les `.inp`)

**(a) Balayage de vitesse** — maillage quad fixe dx = 1 mm, générateur `tunnel_sweep.py` :
`TUN_t300 / t030 / t010 / t003 / t001` → **48 561 éléments** chacun (47 102 CPE4R + 1 459 CPE3).

**(b) Objectivité au maillage à ṗ = 83 MPa/µs** — générateur `tunnel_mesh.py` (TSTEP = 3 µs figé) :

| deck | dx visé | éléments (recomptés) | dx mesuré (`span_of.py`) |
|---|---|---|---|
| `TUN_t003_h200` | 2,0 mm | 12 093 | 1,839 mm |
| `TUN_t003_h100` (= `TUN_t003`) | 1,0 mm | 48 561 | 0,918 mm |
| `TUN_t003_h050` | 0,5 mm | 196 107 | 0,456 mm |
| `TUN_t003_h025` | 0,25 mm | **0 (fichier vide, 1 184 octets)** | **ÉCHEC** |

Le cas 0,25 mm (≈640 k éléments) est **explicitement documenté comme un échec reproductible** : `writeInput` de CAE produit un `.inp` vide (limite mémoire) — `RESULTATS_bench6.md` §4.

**(c) Objectivité liée à ℓc, quad medial-axis, ṗ = 25 MPa/µs (ℓc = 1,926 mm)** — `tunnel_lc.py` + `tunnel_lc2.py` :

| deck | dx | éléments | dx mesuré |
|---|---|---|---|
| `TUN_nc1` | ℓc = 1,926 | 13 047 | 1,777 |
| `TUN_nc2` | ℓc/2 = 0,963 | 52 635 | 0,882 |
| `TUN_nc3` | ℓc/3 = 0,642 | 107 159 | 0,616 |
| `TUN_nc5` | ℓc/5 = 0,385 (uniforme) | 277 327 | 0,384 |
| `TUN_nc10` | ℓc/10 = 0,193 **gradé** (2 mm au bord) | 359 141 | 0,199 |

⚠️ `tunnel_lc2.py` force `setMeshControls(elemShape=QUAD_DOMINATED, algorithm=MEDIAL_AXIS)` — c'est **ce maillage-là qui produit l'artefact** décrit au §3.4.

**(d) Objectivité en maillage NON structuré (triangles CPE3 libres, advancing-front), 2 vitesses** — `tunnel_tri.py`, `tunnel_nc5tri.py`, `tunnel_nc10tri.py` :

| deck | ṗ (MPa/µs) | ℓc (mm) | dx | éléments | dx mesuré |
|---|---|---|---|---|---|
| `TUN_v08tri_nc1/2/3` | 8,333 (TSTEP 30 µs) | 5,116 | ℓc, ℓc/2, ℓc/3 | 3 196 / 12 517 / 26 715 | 5,335 / 2,699 / 1,853 |
| `TUN_v25tri_nc1/2/3` | 25 (TSTEP 10 µs) | 1,926 | ℓc, ℓc/2, ℓc/3 | 22 309 / 88 197 / 197 884 | 2,023 / 1,017 / 0,679 |
| `TUN_nc5tri` | 25 | 1,926 | ℓc/5 = 0,385 uniforme | **540 153** | 0,410 |
| `TUN_v25tri_nc10` | 25 | 1,926 | ℓc/10 gradé | 60 232 | 1,236 (moy. globale) |

**(e) `TUN_TEST`** (juillet 2026) — run de re-validation `vumat_hole.f` vs `vumat_perc.f`, comparé à `TUN_t030_lys` par `_cmp_bench6.py`. Résultat recalculé : **1 253 / 49 361 éléments supprimés (2,54 %)**, D≥0,98 sur 2,58 % ; le run de référence `TUN_t030_lys_field.npz` a `D moyen = 0,0255` — accord exact. (C'est le test qui a fixé la règle « ne pas intervertir les deux VUMAT ».)

### 1.5 Runner

`C:\Users\fuzquianoalricabi\simulations\CONTINUUM\exo_hole_plate\run_bench6_sweep.ps1` — 5 jobs **séquentiels** :
`abaqus.bat job=$j user=vumat_hole.f cpus=14 double=both interactive`, précédé de `. load_env.ps1` (vcvars64 2022 + oneAPI 2024.2).

---

## 2. La VUMAT : `vumat_hole.f`

Chemin : `C:\Users\fuzquianoalricabi\simulations\CONTINUUM\exo_hole_plate\vumat_hole.f` (28 941 octets, 743 lignes).
En-tête interne : *VUMAT_KSTDFH_PSIVAR_PHICAP.F — DP-DFH Bohus, psi(p) plafonnée à phi=51,7*.
Références citées dans le fichier : Shariati, Saadati, Hild arXiv:2201.01870 (DP-DFH) ; Forquin & Hild 2010, *Adv. Appl. Mech.* 44 (obscuration).

### 2.1 Carte `*User Material, constants=10` (identique dans TOUS les decks TUN_* et EXC_*)

```
*User Material, constants=10
52000.,    0.25,    51.7,   153.3,     15.,     24.,    23.5,      1.
    0.38, 4.18879
```

| # | nom | valeur | sens |
|---|---|---|---|
| 1 | E | 52 000 MPa | Young |
| 2 | ν | 0,25 | Poisson |
| 3 | β | 51,7° | friction Drucker-Prager |
| 4 | d | 153,3 MPa | cohésion DP |
| 5 | ψ | 15° — **IGNORÉ** | remplacé par ψ(p̄) variable, voir §2.2 |
| 6 | **m** | **24** | module de Weibull |
| 7 | **σ_w** | **23,5 MPa** | échelle de Weibull à Z_eff (carte *quasi-statique* Bohus cavité ; la carte percussion utilise 120) |
| 8 | **Z_eff** | **1 mm³** | volume de référence Weibull |
| 9 | **k** | **0,38** | vitesse d'obscuration (relative à c) |
| 10 | **S** | **4,18879 = 4π/3** | facteur de forme du volume obscurci |

⚠️ Gotcha documenté (`RESULTATS_bench6.md` en-tête) : la carte doit être **scindée 8 + 2** (Abaqus ne lit que 8 constantes par ligne) — sinon k et S sont mal lus. Le `.inp` généré le fait correctement.

`V_el` (volume élémentaire Weibull) est **figé en dur à 1** dans `vumat_hole.f` (lignes `V_el = ONE`), donc σ_k = σ_w·(Z_eff/V_el)^(1/m) = 23,5 MPa constant. C'est le réglage « banc » ; la variante structurale correcte (`rc99`, V_el = charLength³) est **une autre VUMAT** — cf. CONTINUUM.md §3 « Pièges ».

### 2.2 Ce que fait la loi

**(1) Plasticité Drucker-Prager, sur contraintes EFFECTIVES.**
f = q − p̄·tanβ − d, potentiel non associé g = q − p̄·tanψ, retour radial `dlam = f/(3G + K·tanβ·tanψ)`, apex traité.
**Option « fidélité Saadati 2022 »** : le retour DP n'est appliqué **que si p̄ > 0** (compression). En traction, la contrainte principale monte élastiquement jusqu'au seuil DFH — c'est le DFH seul qui gère la rupture en traction.
**Dilatance ψ(p̄) variable** : `psivar = clamp(160.345 − 0.213793·p̄, 0, 51.7)` — 62° à 460 MPa, 0° à 750 MPa, plafonnée à φ = 51,7° pour l'admissibilité.

**(2) Endommagement DFH — anisotrope, 3 scalaires dans un repère FIGÉ.** C'est le point clé pour la transposition.

- **Tirage Weibull par élément, déterministe** (`kst_seed`) : hash spatial 64 bits (xorshift64) des coordonnées **initiales** `coordMp` → 3 uniformes → inversion Weibull `sc(i) = σ_k·(−ln(1−u))^(1/m)`, puis **tri croissant**. Pas de `RANDOM_NUMBER` → reproductible d'un run à l'autre, mais **dépendant du maillage** (le tirage change si les centroïdes changent).
- **Gel de la direction (§4a du code)** : au **premier amorçage** (σ_I effective ≥ sc(1), la plus faible des 3 résistances), on diagonalise la contrainte effective, on construit R = [v_max | v_mid | v_min] avec det(R) = +1, et on stocke **les angles d'Euler ZYX dans SDV 7-9**. **Cette base ne tourne plus jamais** (repère corotationnel Abaqus mis à part). `ti(1) > 0` ⇔ repère figé.
- **Amorçage par direction** : dans ce repère figé, chaque direction i s'amorce indépendamment quand σ_ii ≥ sc(i), horodatée dans `ti(i)` (SDV 13-15).
- **Croissance = obscuration locale**, forme fermée :
  `D_i = 1 − exp(−x_i³)`, avec l'incrément `dx_i = (S·λ_t)^(1/3) · k · c · dt`, `c = √(E/ρ)`,
  `λ_t = max( (σ_ii/σ_w)^m / Z_eff , 1/V_el )` (densité de fissures Weibull, avec **plancher d'une fissure par élément**).
  Le code reconstruit `x_i = (−ln(1−D_i))^(1/3)` à chaque pas → **aucun SDV supplémentaire**. Exact pour λ_t constant : D = 1 − exp(−λ_t·S·(k·c·(t−t_ini))³). C'est explicitement présenté comme un **remplacement de la chaîne à 3 intégrateurs de Forquin-Hild 2010** (3 intégrateurs × 3 directions ne tiennent pas dans 16 SDV).
  → **le « volume d'ombre » n'est jamais construit géométriquement** : il est replié analytiquement dans Z_o(t) = S·(k·c·t)³, c'est-à-dire, avec S = 4π/3, **une sphère de rayon k·c·t centrée sur la fissure amorcée, croissant à k·c = 0,38 × 4455 = 1 693 m/s**.
- **Contrainte nominale** : dans le repère figé, `σ_ii ← (1−D_i)·σ_ii` **si σ_ii > 0 seulement** (traction) ; cisaillements `σ_ij ← min(f_i,f_j)·σ_ij`. **Fissure fermée (compression) = plein transfert** — pas de dégradation.
  L'effectif est reconstruit exactement à l'incrément suivant en divisant par les mêmes facteurs (inverse exact).
- **Plafonds** : `DCAP = 0.99` (raideur résiduelle 1 %, empêche la distorsion des éléments fissurés) ; **`DELD = 0.98`** = seuil de suppression d'élément sur DMAX. Dans `vumat_perc.f` (confiné/percussion) la même ligne vaut **`1.0d9`** = aucune suppression.

**Résumé pour la question posée** : l'endommagement est **anisotrope mais pas tensoriel** — 3 scalaires D₁, D₂, D₃ portés par une **triade orthonormée figée à l'instant du premier amorçage**, stockée en angles d'Euler. Pas de rotation ultérieure, pas de tenseur d'endommagement d'ordre 2 ou 4, pas de suivi de fissure discrète.

### 2.3 Carte SDV (`*Depvar` = 16)

| SDV | contenu |
|---|---|
| 1 | STATUS (1 actif / 0 supprimé) |
| 2 | **DMAX = max(D₁,D₂,D₃)** ← c'est ce que lisent tous les extracteurs |
| 3 | PEEQ (plastique équivalent DP) |
| 4-6 | D₁, D₂, D₃ |
| 7-9 | Euler ZYX du repère figé (rad) |
| 10-12 | résistances Weibull tirées sc₁≤sc₂≤sc₃ (MPa) |
| 13-15 | temps d'amorçage t_ini par direction (s) |
| 16 | σ_max historique effective (diagnostic) |

### 2.4 Les autres VUMAT présentes dans le dossier (contexte, non utilisées par le banc 6)

`vumat_kstdfh.f`, `vumat_kstdfh_psivar.f`, `vumat_kstdfh_psivar_rc99.f` (**V_el = charLength³**), `..._rescap.f`, `..._rc99_visc.f` (viscosité Ḋ ≤ 1/τ), `vumat_psivar_rc99_veff1.f`, `vumat_veff1_visc.f`, `vumat_dfhcap.f`, `vumat_cyclic.f`, et `confine_lc5\vumat_perc.f`.

---

## 3. Les résultats sur disque

### 3.1 Données brutes

- **`bench6_sweep.csv`** (20 lignes) : colonnes `base,tstep_us,pdot_MPaus,N_cracks,r_count_mm,core_r_mm,Rd_mm,n_cracked,area_cr_mm2,Ediss,ALLIE,ALLPD,ALLWK,ALLKE,ALLAE`.
  ⚠️ **Piège de lecture** : pour toutes les lignes `TUN_nc*` et `TUN_v*tri_*`, les colonnes `tstep_us` et `pdot_MPaus` valent **0** — c'est un artefact du regex `_t(\d+)` de `crack_count.py` sur les noms `nc`/`tri`, **pas** une valeur physique. Les vrais taux sont : `nc*` et `v25tri_*` et `nc5tri` → ṗ = 25 MPa/µs ; `v08tri_*` → ṗ = 8,33 MPa/µs (lus dans `tunnel_lc.py`, `tunnel_tri.py`).
- **20 fichiers `TUN_*_lys_field.npz`** : champ D par élément (`X4`, `Y4` polygones, `D` = SDV2, `labels`, `ft`), produits par `dump_field.py` (`abaqus python dump_field.py TUN_t300_lys`). Élément supprimé → absent de l'odb → `D = NaN` (tracé « fissure complète »).
- **`crack_count_cc_result.npz`, `arms_spacing_result.npz`, `arms_v2_result.npz`** : résultats des méthodes de comptage indépendantes (contenus dumpés au §3.3).
- odb conservés localement : `TUN_TEST.odb` (188 Mo), `confine_lc5\EXO_TEST.odb` (183 Mo). **Aucun JSON** dans le dossier.
- Livrables zippés : `bench6_livrable.zip` (20 Mo), `partage_bench6.zip` (20 Mo), dossier `partage_bench6\{data,figures,scripts}`.

### 3.2 Scripts de dépouillement (tous dans la racine)

| script | rôle |
|---|---|
| `crack_count.py` | `abaqus python crack_count.py TUN_t001` → N **angulaire** (clusters d'angle sur l'anneau [max(30, r_cœur+8), 80] mm, seuil D>0,5, 360 bins 1°, CMIN=3, fusion de trous ≤1 bin) + rayon du cœur plein (fraction ≥ 0,85) + Rd + aire fissurée + énergies. Append dans `bench6_sweep.csv`. |
| `span_of.py` | **N_span** = clusters **connexes 2D** (union-find, lien = 1,8·dx, adaptatif au maillage) dont l'extension radiale ≥ **12 mm**. C'est l'observable retenue comme robuste. |
| `crack_count_cc.py` | méthode connexe de validation (N_shell + N_span) |
| `count_arms_spacing.py`, `count_arms_v2.py` | comptage par espacement inter-bras |
| `angular_fft_crackcount.py` | FFT angulaire (→ `angular_fft_result.png`) |
| `raycast_cracks.py` | ray-casting (→ `raycast_diagnostic.png`) |
| `radiality_analysis.py`, `count_profile.py`, `analyze_sweep.py` | contrôles annexes |
| `dump_field.py`, `extract_tunnel.py`, `extract_fields_sweep.py` | extraction odb → npz |

### 3.3 Ce que montrent les résultats — **N, motif, effet du taux**

**(a) Le motif est RADIAL et DISCRET.** σ_θθ = +p(a/r)² en traction ⇒ fissures radiales émanant du trou (mode attendu). Trois régimes nets (`RESULTATS_bench6.md` §1, figure `fig_bench6_montage.pdf`) :

| ṗ (MPa/µs) | régime | motif observé |
|---|---|---|
| 0,83 | quasi-statique | **peu de fissures radiales LONGUES**, jusqu'au bord (maillon-faible) |
| 8,3 | dynamique | bras radiaux compacts, cœur dense |
| 25 | transition | plus de bras |
| 83 | dynamique | multiplication des bras, cœur ≈ R_d |
| 250 | **comminution** | **nuage aréal diffus, plus de bras discrets** |

Rd de Lamé QS = a·√(p/σc) ≈ 33 mm — cohérent avec la taille du cœur dense.

**(b) N croît avec ṗ — la signature d'obscuration EST présente.** Valeurs `N_span` **recalculées par mes soins** avec `span_of.py` sur les npz du disque, identiques à celles publiées :

| ṗ (MPa/µs) | 0,83 | 8,33 | 25 | 83,3 | 250 |
|---|---|---|---|---|---|
| **N_span** | 4 | 6 | 7 | 10 | 15 |
| N_shell (crack_count_cc) | 19 | 10,4 | 12,4 | 17 | 60,4 |
| N angulaire (`bench6_sweep.csv`) | 9 | 9 | 9 | 18 | 1 (*) |
| **fraction fissurée [15;80] mm** | 0,086 | 0,036 | 0,050 | 0,125 | **0,302** |
| rayon du cœur plein r_c (mm) | 17 | 17 | 19 | 27 | 28 |
| n éléments fissurés | 2 788 | 1 253 | 1 606 | 3 429 | **8 267** |
| E_diss (mJ) | 206,3 | 18,45 | 37,26 | 131,7 | 180,1 |

(*) le N angulaire tombe à 1 à ṗ = 250 : l'anneau est intégralement endommagé, le comptage angulaire dégénère — c'est la signature de la comminution.

**(c) L'exposant est PLAFONNÉ à ~0,27, pas 0,89.** Prédiction point-matériel 1D : N ∝ ṗ^(m/(m+3)) = ṗ^0,889. Ajustements **recalculés** :

| fenêtre | b mesuré |
|---|---|
| t030→t003 (dynamique non saturé) | **0,223** |
| t030→t001 | **0,272** |
| 5 points | 0,224 |
| `arms_spacing_result.npz` | b_dyn = 0,241 / b_full = 0,185 |
| `arms_v2_result.npz` | b_dyn = 0,282 / b_full = 0,105 |

**5 méthodes indépendantes** (histogramme angulaire, espacement inter-bras, FFT angulaire, ray-casting, clusters connexes 2D) rejettent 0,89 de façon décisive et convergent sur **b ≈ 0,22–0,31**. Aucune sous-fenêtre n'atteint 0,5. Figure : `fig_bench6_Nexp.pdf/.png`.

**Mécanisme du plafonnement, tel qu'écrit** (`RESULTATS_bench6.md` §2) : au-delà de ṗ ≈ 80 MPa/µs, le mécanisme **bascule de la multiplication de bras radiaux discrets vers la comminution aréale**. À ṗ = 250 : 239 clusters dont **15 seulement** s'étendent ≥ 12 mm ; fraction fissurée 0,30. Cause invoquée : **p = 250 MPa ≫ σc dynamique (~30-37 MPa)** ⇒ le **front d'onde lui-même** endommage sur son passage (endommagement réparti), en plus des fissures radiales.

**(d) Objectivité au maillage — ÉNERGIE et MOTIF objectifs, tracé exact non.**

Balayage à ṗ = 83 MPa/µs (`fig_bench6_mesh_fields.pdf`, `fig_bench6_mesh_conv.pdf`) :

| dx (mm) | éléments | N_span | N angulaire | fraction | E_diss (mJ) |
|---|---|---|---|---|---|
| 1,84 | 12 093 | 11 | 20 | 0,149 | 125 |
| 0,92 | 48 561 | 10 | 18 | 0,125 | 132 |
| 0,46 | 196 107 | 16 | 23 | 0,138 | 144 |

Balayage lié à ℓc à ṗ = 25 MPa/µs, **quad medial-axis** (`fig_lc.py`, valeurs recalculées `span_of.py` en accord) :

| dx (mm) | ℓc/dx | N_span | fraction | E_diss (mJ) |
|---|---|---|---|---|
| 1,777 (=ℓc) | 1,08 | 6 | 0,049 | 35,48 |
| 0,882 (ℓc/2) | 2,18 | 8 | 0,047 | 36,29 |
| 0,616 (ℓc/3) | 3,13 | 9 | 0,051 | 35,68 |
| 0,384 (ℓc/5) | 5,02 | 7 | 0,058 | 35,80 |
| 0,199 (ℓc/10) | 9,68 | 9 | **0,105** | 35,06 |

→ **E_diss constante à ±2 % jusqu'à ℓc/10** ; la fraction fissurée, elle, « explose » à ℓc/10.

**(e) L'explosion à ℓc/10 était un ARTEFACT DU MAILLAGE STRUCTURÉ.** Refait en **triangles CPE3 libres** (`fig_tri.py`, `fig_cracks.pdf`) :

| ṗ = 25, ℓc = 1,926 | ℓc/dx | E_diss (mJ) | fraction |
|---|---|---|---|
| `v25tri_nc1` | 0,95 | 28,05 | 0,045 |
| `v25tri_nc2` | 1,89 | 28,27 | 0,040 |
| `v25tri_nc3` | 2,84 | 30,32 | 0,043 |
| `nc5tri` | 4,70 | 31,34 | 0,048 |

| ṗ = 8,33, ℓc = 5,116 | ℓc/dx | E_diss (mJ) | fraction |
|---|---|---|---|
| `v08tri_nc1` | 0,96 | 6,327 | 0,014 |
| `v08tri_nc2` | 1,90 | 5,69 | 0,035 |
| `v08tri_nc3` | 2,76 | 6,434 | 0,026 |

Conclusion écrite (`SYNTHESE.md` §3) : **fraction STABLE (~0,044) en triangles ; en quad medial-axis elle passait de 0,049 à 0,105 avec des rayons parasites alignés sur les axes (N_span 7 vs 2) et ~+15 % d'énergie.** ⇒ **mailler NON structuré est obligatoire pour le faciès.** L'énergie converge à ±5 % dès **dx < ℓc/2**, aux deux vitesses ; descendre à ℓc/5 ou ℓc/10 n'apporte plus rien.

**Limite explicitement assumée** : les *statistiques* convergent, mais le **tracé exact des fissures diffère d'un maillage à l'autre** (le tirage Weibull est par élément). Pour un faciès reproductible il faudrait un **champ aléatoire de σ_w à longueur de corrélation fixe** (grains ~1 mm, `*Initial Conditions, type=SOLUTION`) — **NON FAIT à ce jour**.

**(f) Contrôles honnêtes** (`RESULTATS_bench6.md` §4) :
- Lysmer OK : pas d'anneau d'endommagement circonférentiel au bord, pas de fissures circonférentielles fantômes. **Mais** à p = 250 ≫ σc le front d'onde produit un endommagement réparti jusqu'aux coins (r > 90 mm), et les cas lents montrent des fissures radiales atteignant le bord (domaine fini).
- **La suppression d'élément (DELD = 0,98) n'a jamais été déclenchée dans le sweep**, D plafonnant à DCAP = 0,99 ; la raideur résiduelle 1 % a suffi (aucun crash de distorsion en 2D). *(Nuance : le run de contrôle `TUN_TEST` montre 2,54 % d'éléments supprimés — donc la suppression peut mordre selon le cas.)*
- Le point QS (t300) est **hors loi de puissance** (plateau + fissures longues type maillon-faible) : **exclu du fit dynamique**.

### 3.4 Figures présentes (PDF vectoriel + PNG + script `.py` pour chacune)

`fig_bench6_Nexp` (N vs ṗ, log-log, fit + référence 0,89) · `fig_bench6_montage` (5 champs D du sweep) · `fig_bench6_mesh_fields` (champs des 3 maillages) · `fig_bench6_mesh_conv` (N et E vs dx) · `fig_lc` (champs ℓc, ℓc/3, ℓc/5 + convergence) · `fig_tri` (E normalisée et fraction vs ℓc/dx, tri vs quad) · `fig_cracks` (faciès 2 vitesses × 3 raffinements, triangles) · `fig_cap`, `fig_hs`, `fig_rate_m`, `fig_triax`, `fig_cyclic`, `fig_cyclic_cycles` (études annexes DP/cap/triaxial/cyclique du même dossier).
Style : `font.family serif`, CMU Serif, `mathtext.fontset cm`, sortie PDF vectoriel — conforme aux règles du projet.
Champs bruts en PNG : `field_t300_t030_t010_t003_t001.png`, `field_nc5_nc5tri.png`, `sweep_fields.png`, `hole_meshes.png`.

### 3.5 Campagne CONFINÉE (`confine_lc5\`, 15 runs `EXC_*`) — pour mémoire

Méthode : pression isotrope sur les 4 bords extérieurs, rampée en Step-1 (150 µs, quasi-statique) et **maintenue** pendant la pressurisation de la cavité en Step-2. Triangles CPE3 non structurés gradués, dx = ℓc/5 au trou. VUMAT `vumat_perc.f` (**DELD = 1e9, aucune suppression**). 3 taux (8,3 / 25 / 83) × 5 confinements (0/20/50/75/100 MPa). Coût mesuré 13–74 s/run.

Fraction d'aire fissurée (D > 0,5) :

| ṗ \ P (MPa) | 0 | 20 | 50 | 75 | 100 |
|---|---|---|---|---|---|
| 8,3 | 0,736 | 0,365 | 0,111 | 0,042 | 0,011 |
| 25 | 0,662 | 0,324 | 0,117 | 0,046 | 0,018 |
| 83 | 0,734 | 0,432 | 0,120 | 0,048 | 0,022 |

→ **le confinement supprime la fissuration d'un facteur ~35** ; ALLPD ÷24 à ÷153 selon ṗ. Effet de taux **secondaire** devant l'effet de confinement. ⚠️ à P ≤ 20 les fissures touchent le bord (Rd ≳ 100 mm = demi-largeur) → absolus plafonnés ; tendance propre dès P ≥ 50. Figure `confined_exo.pdf/png`.

---

## 4. `compute_lc.py` — la formule exacte et les valeurs

Chemin : `C:\Users\fuzquianoalricabi\simulations\CONTINUUM\exo_hole_plate\compute_lc.py` (1 232 octets).
Copie archivée : `...\OBJECTIVITE_MAILLAGE\01_theorie_lc\compute_lc.py`.

**Formule** (rapport DP-DFH éq. 45-46) : ℓc = **k · c · t_c**, avec c = √(E/ρ) et t_c donné par la saturation de la probabilité d'obscuration P_obs(t_c) = 1 :

> t_c = [ (m+1)(m+2)(m+3)·Z_eff / (6·S·(kc)³) · (σ_w/Σ̇)^m ]^(1/(m+3))

Implémentation en **log** pour éviter l'overflow de (σ_w/Σ̇)^m :
`ln t_c = [ ln((m+1)(m+2)(m+3)Z_eff) − ln(6·S·(kc)³) + m·ln(σ_w/Σ̇) ] / (m+3)`

Le taux de contrainte pariétal est pris **égal à ṗ** (`Sdot = pdot_MPa_us * 1e6` MPa/s). Constantes du script : E = 52 000, ρ = 2,62e-9, σ_w = 23,5, m = 24, k = 0,38, S = 4,18879, Z_eff = 1.

**Sorties exactes (script exécuté)** :

```
c = 4.4550e+06 mm/s = 4455 m/s | kc = 1.6929e+06 mm/s
pdot         tc(us)     sigc   lc(mm)     lc/2     lc/3
0.83        23.4741     19.5   39.740   19.870   13.247
8.33         3.0221     25.2    5.116    2.558    1.705
25.00        1.1377     28.4    1.926    0.963    0.642
83.33        0.3902     32.5    0.661    0.330    0.220
250.00       0.1469     36.7    0.249    0.124    0.083
```

Lecture : ℓc ∝ Σ̇^(−m/(m+3)) — **la longueur caractéristique dépend du taux et diverge en quasi-statique** (39,7 mm à ṗ = 0,83 ; 0,25 mm à ṗ = 250). σ_c est la résistance dynamique induite (19,5 → 36,7 MPa). Pour la carte percussion (σ_w = 120), la synthèse donne ℓc ≈ 2,31 mm à ε̇ ~ 2000 s⁻¹ (V = 11 m/s).

**Règle opérationnelle qui en découle** (checklist obligatoire, `SYNTHESE.md` §5 / `phd/CONTINUUM.md` §4) :
1. calculer ℓc au taux de l'étude et mailler **dx < ℓc/2** dans la zone utile (ℓc/5 et au-delà : gain nul) ;
2. maillage **NON structuré** (CPE3 libres en 2D, tétras gmsh en 3D) ;
3. **pas de mass scaling global** (gonfle ρ → baisse c → **fausse ℓc**) ; sélectif à plancher constant toléré pour des tendances ;
4. V_el = charLength³ pour toute étude structurale (V_el = 1 réservé aux bancs point-matériel) ;
5. DELD = 1e9 pour le continuum ; DELD ≈ 0,98 réservé aux cas où l'évacuation de matière est le sujet (banc 6) ;
6. métriques d'objectivité = **énergie / volume endommagé / pénétration** (convergent). **L'étendue discrète des fissures ne converge pas** sans champ σ_w corrélé.

---

## 5. Ce qui est ABSENT (à ne pas inventer)

- **ABSENT** : tout fichier JSON de résultats (aucun `.json` dans le dossier).
- **ABSENT** : le maillage dx = 0,25 mm à ṗ = 83 (`TUN_t003_h025.inp` fait 1 184 octets et **0 élément** — échec `writeInput` documenté, reproductible).
- **ABSENT** : un run `TUN_t003_h025_lys` — le script `fig_bench6_mesh_conv.py` le prévoit en argument optionnel mais aucun `.npz` correspondant n'existe.
- **ABSENT** : le champ aléatoire de σ_w à longueur de corrélation fixe (perspective identifiée, jamais implémentée) — c'est la raison écrite pour laquelle le **tracé** des fissures n'est pas reproductible d'un maillage à l'autre.
- **ABSENT** : toute mesure d'orientation/statistique des directions gelées (SDV 7-9 sont écrites par la VUMAT mais **aucun script du dossier ne les dépouille** — tous lisent SDV2 = DMAX).
- **ABSENT** : une étude d'objectivité avec V_el = charLength³ **sur ce cas 2D** (les decks TUN utilisent tous V_el = 1 figé). La variante rc99 n'a été balayée qu'en percussion 3D.
- **ABSENT** : valeur d'exposant approchant 0,5 ou 0,89 dans une quelconque sous-fenêtre (explicitement testé et rejeté).
- **Valeurs à ne pas lire au premier degré** : `tstep_us = 0` et `pdot_MPaus = 0` dans `bench6_sweep.csv` pour les lignes `nc*`/`tri*` (artefact de regex, cf. §3.1).

---

## 6. Lecture pour la transposition FDEM (ce que le banc 6 établit factuellement)

Ce que la logique DP-DFH **fait réellement** dans cette étude, et qui doit être confronté au solveur FDEM :

1. **Il n'y a pas de « zone d'ombre » géométrique.** L'obscuration est repliée analytiquement : Z_o(t) = S·(k·c·t)³ avec S = 4π/3, soit une **sphère de rayon k·c·t** (k·c = 1 693 m/s, soit 38 % de c) autour de chaque site amorcé, et D = 1 − exp(−λ_t·Z_o). Aucun voisinage n'est parcouru, aucun test de visibilité n'est fait. C'est une **statistique de Poisson locale**, pas une interaction entre éléments.
2. **La direction est gelée une fois pour toutes** au premier amorçage (repère principal de la contrainte **effective** à cet instant, stocké en Euler ZYX). Ce n'est pas un tenseur d'endommagement : ce sont **3 scalaires portés par une triade figée**, avec dégradation **en traction seulement** et plein transfert en compression.
3. **La régularisation vient de ℓc = k·c·t_c**, rate-dépendante, et **elle marche** : énergie dissipée finie et stable (±2 à ±5 %) sur un facteur 4 à 10 de raffinement, motif (cœur ≈ Rd + bras radiaux aux mêmes positions) stable, **pas d'effondrement sur un élément** — contrairement à l'adoucissement local du banc 4.
4. **Ce qui NE se transpose pas** : (a) l'exposant m/(m+3) = 0,89 du fragment 1D **ne se retrouve pas** sur un comptage de fissures 2D structural (0,22–0,27 mesuré, 5 méthodes) ; (b) le **tracé** des fissures dépend du maillage via le tirage Weibull par élément ; (c) le maillage structuré **fabrique** des rayons alignés sur les axes et gonfle la fraction fissurée de 0,049 à 0,105 — un solveur FDEM avec un maillage régulier hériterait du même biais.
5. **Point de bascule physique identifié** : dès que p ≫ σc dynamique, le front d'onde endommage sur son passage et le mécanisme passe de « multiplier des fissures propres » à « pulvériser » — la loi ne sature pas, elle **change de mode**.