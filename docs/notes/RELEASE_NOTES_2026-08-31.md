# rockim — état du 2026-08-31 (pour mise à jour GitHub)
*Base : bundle du 2026-08-30 (chantier A11/B8/B10 + 107 contrôles de suite —
EN AVANCE sur le GitHub arrêté au 25/08) + les correctifs et ajouts de la
session d'audit/réplique des 30-31/08. AUCUNE modification du C++ du solveur
dans cette session : src/ et include/ sont ceux du bundle, byte pour byte.*

## Correctifs (fichiers modifiés)
* **CMakeLists.txt** — accepte Eigen >= 3.3 SANS borne haute (Eigen 5.x rejeté
  avant : compatibilité même-version-majeure de son fichier de version) ;
  repli téléchargement en forme directe non dépréciée (CMP0169, CMake >= 4) ;
  messages explicites hors-ligne. Vérifié : configure + build OK avec Eigen
  3.4.0 et 5.0.1, code identique 0 warning.
* **tools/verify_suite.py** — `--update-refs <f.json>` / `--refs <f.json>`
  implémentés (l'en-tête les promettait) : baselines de plateforme pour les
  références CHIFFRÉES seulement, invariants toujours en dur ; _meta de
  provenance dans le JSON. Vérifié : 44/44 sur macOS ARM contre
  refs_macos_arm.json (fourni), qui échouait 34/44 contre la baseline Linux.
* **LISEZ_MOI.md** — l'avertissement Apple Silicon (« un ou deux échecs de
  justesse ») remplacé par la mesure (10/44, bit-identiques entre Eigen 3.4 et
  5.0, invariants tous verts) et la procédure de baseline par plateforme.
* **tools/make_unstructured_mesh.py** — Mesh.Algorithm = 5 (Delaunay) partout :
  le frontal (6) pavait en quasi-équilatéraux alignés (fabric d'orientation
  7,9 mesuré contre 1,4) — trois directions de fissure imposées.
* **tools/make_impact_mesh.py** — sous champ de fond gradué, frontal ET
  Delaunay pavent la zone fine en hexagones (fabric 10,8) : MeshAdapt (1) +
  OptimizeNetgen mesurés seuls corrects (fabric 1,9, h_min 0,43 mm à l'échelle
  test). Matrice A/B/C/D/B2 du 2026-08-30.
* **bench_impact/tools/fig_impact.py** — titre paramétrable `--title` (il était
  codé en dur « St Anne 10,66 m/s »). ⚠️ Défaut connu à corriger : la ligne
  « fissure radiale max » duplique le rayon de cratère ; la vraie mesure de
  radiales est celle de crater_metrics.py (portée par secteur au-delà du cratère).

## Ajouts
* **bench_impact/tools/fig_pulv.py** — équivalent de la fig. 18 de Yang et al.
  2026 (nPulv(t) + bdWork), le seul manque de post-traitement.
* **bench_impact/configs/impact_kuru9*.cfg** — les decks de réplique du papier 2
  (IJRMMS 206:106660), Table 1 vérifiée mot à mot : pilote (s15), fidèle (s10),
  intrinsèque, et article-exact (intrinsic + strainRateDIFArm = envelope +
  contactMu 0,18 + jointSoftening munjiza — toutes divergences de modèle fermées).
* **meshes/impact_kuru_s15.msh / _s10.msh** — maillages gradués 1/2/10 mm
  (mailleur corrigé, graine 2 pour le s10 : meilleur pire-élément sur 3 graines).
* **configs/fdem3d_percussion_{pulv,weibull,etoile,etoile_v2}.cfg** — decks de
  démonstration (pulvérisation armée, champ Weibull corrélé, campagne étoile).
* **refs_macos_arm.json** — baseline de plateforme macOS ARM (tier fast, 1 fil,
  Apple Clang 21/libc++, 2026-08-30).

## Résultats de réplique à connaître (run fidèle 9 m/s, 2026-08-30)
Résidu d'énergie 8,1e-5 % ; cratère 9,6 mm (cible 10-12) ✓ ; MAIS vitesse
d'indentation +25 % (transfert piston à caler), perte d'énergie du bit 68 %
(cible 31,6 %), radiales quasi absentes (2 mm contre 13-16). Voir
RAPPORT_replique_kuru9_2026-08-30.md (artifact de session). Suspects par
ordre : transfert piston -> contactMu. Le run article-exact (en cours au
moment de cet export) teste le second.

## Non inclus (régénérables)
* meshes/box3d_star_*.msh (15 Mo, campagne étoile) :
  `python tools/make_unstructured_mesh.py box3d 0.12 0.12 0.08 0.0025 meshes/box3d_star_h25_delaunay.msh 1`
* Les sorties de runs (out_*) et figures de session — dans les artifacts.
