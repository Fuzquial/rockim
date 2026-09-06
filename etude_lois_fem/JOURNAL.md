# Journal du chantier « briques constitutives en FEM continu » (rockim fem3d)

*Proposition : [PROPOSITION_etude_lois_fem.md](PROPOSITION_etude_lois_fem.md). Dossier vivant du code :
`FDEM/rockim_f2/` (sources `src/`, `include/rockim/`), pas le worktree. Règles : un job à la fois,
14 threads ; validation de la liste de clés + coût avant chaque lancement de MATRICE ; croissance par
addition (clés opt-in, défaut bit-identique) ; banc court qui doit échouer avant tout run long.*

## Binaires (tous dans `FDEM/rockim_f2/`, un nom par build)

| exe | sources | usage |
|---|---|---|
| `rockim_f2w0.exe` (03/09 23:29) | vivantes INTACTES (identiques à f2p pour fem3d) | référence de bit-identité (`bitid_w0/`), phase A |
| `rockim_f2w1*.exe` | intermédiaires (liens PARTIELS après modif d'en-tête : INVALIDES pour fem3d) | ne pas utiliser |
| `rockim_f2w2.exe` (03/09 23:56) | items 1-5 sans `eroCode`, dpr lisant `capP0` sans garde | phase B (bancs courts), `selftest-triax` OK |
| `rockim_f2w3.exe` (04/09 00:04) | + `eroCode` (canal d'érosion), `dprCap` (garde du cap sur dpr) | **bit-identique w0 sur les 10 configs** (`bitid_w3/`) |
| `rockim_f2w4.exe` (00:21) | + lame : corps = derrière les DEUX faces (B4 : 15 MN avant) ; `erodeStrainMax` (plafond d'étirement, B3 : cascade sans lui) | bancs B2-B4 rejoués |
| `rockim_f2w5.exe` (00:33) | + `erodeWfrac` (spall sur l'énergie dissipée) | A2b : wPlas −2·10⁴ J → défaut d'apex trouvé |
| `rockim_f2w6.exe` (00:53) | + `dpApex` (pas de retour DP au-delà de l'apex), lame : corps borné | phase B : T3 OK, lame µ=0 exacte, mini-bloc stable sans érosion, E1 bande 1,9 lc (h) / 3,0 lc (h/2) ; A2b sain (wPlas +10,7 J, 38 érodés) |
| `rockim_f2w7.exe` (01:31) | + `dpTension = off`, lame par bissectrice, jauge σ_xx en tranche, `bottomPressure`, wPlas nominal | E1 objectif (1,05), A2b sain (wPlas +8,4 J), bit-identique w0 ; la bissectrice fait perdre le contact de la lame |
| **`rockim_f2w8.exe` (01:57)** | lame : règle par hauteur rétablie | **binaire de la matrice** : lame exacte (tan 20° à 10° et 60° de dépouille), témoin du noyau d'origine (wPlas −2,4·10⁴ J), **bit-identique w0 (10 configs)** |
| `rockim_f2w11.exe` (04/09 13:39, rebuild de revue 14:54) | + loi `cdp` opt-in (Concrete Damaged Plasticity), `selftest-cdp`, `matpoint`, fem3d `pullDelay` / `gripSection` / `triaxStats` — chantier `CONTINUUM/calib_bohus_triax/cdp_rockim/`, doc §5.19 ; revue : ftScale honoré par cdp, `cdpViscosity`, clés `cdp*` inconnues refusées, pilote `Hold1D` réécrit (crochet + Illinois), colonne `residu` / rc 2 de matpoint, exception de la loi relancée hors région OpenMP (`elementForces`) | **bit-identique w10 sur les 10 configs + `fem3d_shear_short` à OMP 14 (`history`, `frames`, `.vtu`)** et sur 7 configs à OMP 4 (`bitid_w11/`, `comparaison*.txt`) ; ne remplace PAS w10 pour la matrice ; chaîne OMP 14 finie : 11/11 identiques. **Rebuild 15:16 (soir)** : + `cdpCompLength` (crack band en compression, opt-in : branche post-pic des tables de compression remise à l'échelle L_ref/lc, gardes G2c/G3/G4), banc (m) de `selftest-cdp` (79 [OK]), colonne `eps_in_c` de matpoint — bit-identique w10 sur 3 configs à OMP 4 (`bitid_w11/w11/`, `comparaison.txt`) |
| `rockim_f2w15.exe` (05/09 05:44, build complet en place, `Remove-Item *.obj` + `build_f2.cmd`) | + `bulkViscosity = b1 b2` (fem3d, viscosité de volume à la Abaqus/Explicit, opt-in : contrainte visqueuse q I sur les forces internes seulement, colonne `wBulk` en fin de `history.csv` si active, dt CFL × (√(1+b1²) − b1)), `Elem::Jbv` (J du pas précédent), doc §5.20 ¶ bulkViscosity | **bit-identique w14 sur `PQ_cdpI_P020_court` à OMP 2** (`bitid_w15/comparaison_w15.txt`, history/frames/4 vtu, cmp) ; banc `bitid_w15/SELFTEST_bv.md` : `0 0` ≡ sans clé, `0.06 1.2` sonnerie de tête −23 % / wBulk > 0 / dt × 0,94180, `0.3 1.2` −62 % / dt × 0,74403, porte compression du terme quadratique vérifiée ; pas de scénario à volume constant |
| `rockim_f2w16.exe` (05/09 06:13, build complet en place, `Remove-Item *.obj` + `build_f2.cmd`) | + `kinematics = biot \| hencky` (fem3d, opt-in, défaut biot = chemin historique) : hencky = ε = ln U par décomposition spectrale de Fᵀ F (`include/rockim/Fem3dKinematics.hpp`, R = F U⁻¹ exacte), contrainte de la loi lue comme Cauchy co-rotationnelle, forces internes sur la configuration courante P = J R σ U⁻¹ ; doc §5.20 ¶ « Cinématique hencky » | **bit-identique w15 sur `PQ_cdpI_P020_court` à OMP 2** et `kinematics = biot` ≡ sans clé (`bitid_w16/comparaison_w16.txt`, cmp 6/6) ; banc `bitid_w16/SELFTEST_hencky.md` : formes fermées à −0,01 % (biot E(λ−1), hencky λ^(−2ν) E ln λ), mauvaise forme rejetée (+27 % à λ = 0,8), rotation rigide P = 0 (test unitaire 23/23), petites déformations −0,215 % |
| `rockim_f2w19.exe` (05/09 17:15, build complet en place, `Remove-Item *.obj` + `build_f2.cmd` ; après w17 `toolLockXY` et w18 `probes`, non listés ici) | + `tensionDamage = scalar \| fixed` + `tensionShearRetention = β` (dpr/saksala, opt-in, défaut scalar = chemin actuel) : endommagement de traction à direction figée (fixed crack model Rashid 1968 / Rots 1989), jusqu'à trois directions figées dans le complément orthogonal, cinétique de bande de fissuration (f_t, G_f, l_c) par direction, S_ii (1 − d_i) si ouverte, S_ij (1 − β max d^ouv), `MatState::Fcm`, champs VTU/probes `dFix1..3`, `rockim selftest-fixed` ; doc §5.20 ¶ « Endommagement de traction à direction figée » | **bit-identique w18 sur `PQ_cdpI_P020_court` (cdp) ET sur le deck court dpr `bitid_w19/C_T1_R_P000_court_grid` à OMP 2** (`bitid_w19/comparaison_w19.txt`, cmp 6/6 + 6/6, `scalar` explicite ≡ sans clé 6/6, contrôle `fixed` 6/6 différent) ; banc `bitid_w19/selftest_fixed/SELFTEST_fixed.md` 33/33 : raideur y après fissure x = E (fixed) contre 0,1 E (scalar, rejeté), unilatéral, 45°, stress locking 0,242/0,924, G_f à 0,1 % |
| `rockim_f2w20.exe` (05/09 18:08, build complet en place, `Remove-Item *.obj` + `build_f2.cmd`) | gardes des entrées C3/C4/C6 (`include/rockim/Guards.hpp`, registre des clés par mode `tools/keys_by_mode.json` / `KeysByMode.hpp`) : nœud orphelin / masse nulle / tet dégénéré = erreur, NaN → code 3 + ERROR.txt, `hydro*` et clé d'un autre mode refusées ; doc §8.10 | **8/8 IDENTIQUE** contre l'ancre bitid w18 (`tools/bitid.py`, 4 fils) ; banc `bitid_w20/selftest_gardes/` 30/30 ; suite fast 48/48 OMP 1 |
| `rockim_f2w21.exe` (05/09 20:12, build complet en place) | **clé inconnue = erreur** (décision Fernando 20:00, C1) : `Config` suit les clés consommées, `include/rockim/KeyGuard.hpp` (obsolète / autre mode / famille dynamique / suggestion Levenshtein ≤ 2 / inconnue, `unknownKeys = error \| warn`), C7 `config_effective.cfg`, `tools/obsolete_keys.json`, `tools/scan_decks.py` ; doc §8.11 | **8/8 IDENTIQUE** (ancre w18, 4 fils, `bitid_w21/bitid_w21.log`) ; banc `bitid_w21/selftest_cles/` 37/37 ; relecture adverse : D1 (clés à ≥ 2 lecteurs hors mode acceptées) et D2 (`config_effective.cfg` non rejouable) → **remplacé par w22** |
| `rockim_f2w22.exe` (05/09 20:51, build complet en place, `build_w22.log`, sha256 7e683603c53b2f90…) | corrections D1/D2 : règle des **lecteurs** (`kReaders` : clé non consommée légitime ssi le mode courant ou le code partagé la lit, sinon « sans effet en mode Y : cles des modes Z1, Z2 »), `config_effective.cfg` rejouable (défauts commentés, clés du deck non lues gardées actives), `Config::seal()` (plus de verrou ni de chaîne de défaut dans les getters de `step()`) ; `scan_decks.py` : gardes pré-init + fichiers sans clé de solveur | **8/8 IDENTIQUE** (ancre w18, 4 fils, `bitid_w22/bitid_w22.log`) ; banc `selftest_cles` 47/47 (A 21, C 6 rejeux bit-identiques, B 20) ; 0 deck existant sur 957 touché par la règle des lecteurs ; suite fast : `bitid_w22/suite_fast_w22.log` |

## Chronologie

- **03/09 23:29** — build w0, jeu de bit-identité (10 configs fem3d, hachages `history.csv`/`frames.csv`,
  reproductibles à 14 threads). Débit : 4-6·10⁶ tet-pas/s à 14 threads (dpr 6, dpdfh/sk2011 3).
- **23:36 → 23:53 — Phase A** (9 runs, `phaseA/`, w0) : voir `phaseA/RESULTATS_phaseA.md`. Verdicts à
  P = 0 : le cut-off tractif est la brique dominante (÷10 d'érodés en passant à l'apex), le seuil et non
  l'adoucissement ; cap dimensionné inerte ; sans érosion fem3d diverge (A8 NaN à 50 µs) ; l'érosion de
  la référence à ft = 9 MPa n'est PAS un comptage (11 % du bloc). Route R1′ confirmée par les coûts.
- **23:40 → 00:05 — Code items 0-6** (par addition) : `MatLaw` (meridian = power, compDamage = crackband,
  dprCap/capP0, erodeDc, compteurs wPlas/wDamT/wDamC, eroCode, `selftest-triax`) ; `Fem3dSolver`
  (confinement suiveur + topPressure + jauge, toolDelay, symmetryY, activeNodes, erodeDetMin + compteurs
  + eRemoved, toolShape = blade, fieldStats) ; `main.cpp` (selftest-triax) ; doc §5.18 ;
  outillage `make_cfgs.py` (liste blanche), `extract_fem3d.py`, `matrice/make_matrice.py` (42 decks).
- **`selftest-triax` (w2) : OK** — linéaire / mc / corde-par-puissance / pente 3,82 exacts à 10⁻⁷ % ;
  puissance exacte de 20 à 500 MPa (+1,5 % au seul point σ₃ = 0 : tangente verticale) ; contrôles
  « doivent rater les données » OK (−22,8 % à 50, −49,9 % à 20) ; crack band en compression : aire
  adoucie × h = A_c G_IIc à 1,6 %, wDamC = A_c fc0²/2E à 1 %.
- **Piège trouvé** : un lien partiel après modification de `MatState` garde le constructeur inline d'un
  vieux `.obj` (COMDAT) → ω_c lu à 1 avant le premier pas, sans erreur. Règle : rebuild complet après tout
  `.hpp` (mémoire `reference-rockim-rebuild-complet-apres-hpp`).
- **Bit-identité w2 vs w0** : 8/9 configs terminées identiques ; `fem3d_percussion_dpr` différait parce
  que ce deck porte `capP0 = 250e6` (héritage saksala) que dpr lisait désormais → garde `dprCap` ajoutée,
  à re-prouver avec w3.
- **Phase B** (15 bancs courts, `phaseB/`, w2) : en cours à 00:05 ; verdicts par `verdict_phaseB.py`.
- **00:16 → 00:35 — Phase B, première passe (w2 puis w4)** : T3 confiné OK (jauge −99,98 MPa pour −100,
  KE 10⁻¹¹ J, état vierge en dpr puissance ; ressorts = 1 → −75,3 MPa comme prédit) ; sans soupape le
  mini-bloc diverge (NaN) ; le canal ω_c normalisé érode. Trois défauts trouvés et corrigés : (1) la lame
  comptait comme pénétrant tout nœud derrière le plan de coupe (15 MN) → corps = intersection des deux
  demi-espaces ; (2) la soupape det F seule laisse les orphelins étirer leurs voisins (60 % du mini-bloc,
  eRemoved 10¹⁵⁵ J) → `erodeStrainMax` ; (3) **le seuil hérité `erodeD = 0,98` supprime l'élément à 78 %
  de ft** (D = perte de raideur, pas chute de contrainte) : 19 % de Gf dissipés sur la barre E1 →
  `erodeWfrac` (spall sur l'énergie de la bande) ; la référence R passe à erodeD off + erodeWfrac 0,98.
  Critères de script corrigés : KE au maximum (pas à la fin), dt à ny impair réduit de 18 % (pas ×1,5).

- **00:45 → 01:00 — Défaut du noyau d'origine trouvé** : le retour DP radial à p fixe n'a pas d'apex ; au-delà
  de k − 3αp = 0 (traction hydrostatique 18,8 MPa, atteinte dès 0,03 % de déformation effective) dλ > √J₂/G
  et le déviateur change de signe ; la part compressive ainsi fabriquée échappe au split unilatéral →
  injection d'énergie (mini-bloc : 10¹⁵⁵ J en 16 µs, quel que soit le plafond de déformation ; A2b :
  wPlas = −2·10⁴ J). Le seuil hérité `erodeD = 0,98` masquait tout en supprimant les éléments tendus
  à 0,5 % de déformation. Correctif opt-in `dpApex = true` (aucun retour DP au-delà de l'apex, la traction
  revient à Rankine — l'OPTION-3 de la VUMAT DP-DFH) ; **promotion en défaut = décision de Fernando**
  (elle changerait les résultats des démos dpr du dépôt). Les E1 de w5 sont pollués par ce défaut et
  rejoués sur w6. Lame : « face la plus proche » envoyait sous l'outil les nœuds près de la pointe
  (dépouille 60° : l'outil enjambait la matière) → face de coupe au-dessus de la pointe, dépouille en
  dessous, corps borné par sr, sc ≥ 0. Plafond d'étirement passé à 0,3 (le 0,15 du fem 2D, ×2).

- **01:05 → 01:45 — Revue adverse du code (16 agents, 4 lentilles + contre-expertise ; 11 trouvailles
  confirmées sur 12)** : bit-identité à clés absentes prouvée par lecture sur tout le diff (chaque branche
  nouvelle gardée, arithmétique conservée, sorties inchangées) ; corrigé dans w7 : (1) **bloquant** — entre
  p = 0 et l'apex, la jambe tractive du cône DP coulait dans l'espace effectif pendant que Rankine
  endommageait (E1 : 3,3 Gf·A de « plasticité tractive ») → `dpTension = off` (DP en compression seule,
  la vraie OPTION-3) ; (2) lame : la règle « face de coupe au-dessus de la pointe » rendait la dépouille
  inatteignable et créait une falaise de pénalité → face la plus proche DANS le coin (bissectrice) ;
  (3) jauge de confinement en tranche `symmetryY` (lisait (1+ν)/2 σ_xx) → ⟨σ_xx⟩ ; (4) deck de boue
  D_T1_mud_P100 dérivait en corps rigide (pas d'appui) → `bottomPressure` ; (5) soupapes det F et |ε|
  confondues dans la matrice (soup05 / soupOff ne changeaient que det F) → les deux clés ensemble ; étiquette
  « soupape géométrique (det F / |ε|) » ; (6) wPlas compté en contrainte NOMINALE, cap compris.
  **Connu, non corrigé (documenté)** : Y_t = ½σ⁺:C⁻¹:σ⁺ sans terme croisé −ν/E tr σ⁺ tr σ⁻ (la partition
  wDamT/wDamC est approchée en état mixte : ne pas boucler un bilan d'énergie avec) ; la loi puissance a
  une tangente verticale en σ₃ = 0 (+45 % de q à 1 MPa : propriété du modèle Hoek-Brown/MH, sans donnée) ;
  à ω_c → A_c = 0,98 la raideur compressive résiduelle est 2 % (un élément écrasé sous confinement peut
  finir dans la soupape : c'est compté) ; `gfShearFactor` n'est pas lu en fem3d ; `merFc0` n'a pas de
  garde de continuité en σ₃ = 0 (à laisser à sa valeur par défaut).

- **06:08 → 06:35 — Matrice lancée (GO de Fernando) puis ARRÊTÉE après 5 runs** : `C_T1_R_P100` montrait
  100 % du bloc à D ≥ 0,5 AVANT l'impact. **Défaut n° 4 du noyau d'origine** : le cut-off de Rankine est
  piloté par la déformation principale maximale (l'en-tête disait « contrainte effective ») ; la dilatation
  de Poisson d'une zone comprimée (2νP/E = 7,5·10⁻⁴ sous 100 MPa, νσ/E sous l'UCS) dépasse k0 = ft/E =
  1,2·10⁻⁴ et met D ≈ 0,85 partout sans traction (fig. 3 : D ≈ 0,9 dans tout le bloc à P = 0). Le banc
  T3 ne l'avait pas vu (il testait D ≥ 0,9). Clé `rankineDrive = stress` (w9) ; T3 refait avec V(D ≥ 0,5) :
  0 % en contrainte, 100 % en déformation (témoin qui DOIT échouer). E1 en contrainte : bande 1,52 lc.
- **06:40 → 07:10 — Question de Fernando « jamais de maillage structuré ? »** : la grille de Kuhn de fem3d
  EST structurée (miroir + jitter 0,45). Réponse : port du lecteur Gmsh `mesh = file` de fdem3d dans fem3d
  (w10), maillages Delaunay 3D + Netgen : bloc T1 gradué 1,5 mm sous l'insert → 3 mm au loin (37,8 k tets,
  h_in min 0,39 mm), 1,0 mm (88 k), deux réalisations (centre du champ décalé de 0,3 mm) ; tranche T2
  0,5 mm (17,8 k tets, h_in min 0,08-0,10 mm quel que soit l'algorithme : slivers d'une tranche de
  4 couches), 0,35 mm (50 k), deux réalisations (0,52 / 0,48 mm). `dtFactor = 0,7` sur le diamètre inscrit
  (= 0,3 sur l'arête de la grille). Matrice régénérée (`mesh = file`, `rankineDrive = stress`).

- **07:00 — Matrice RELANCÉE sur `rockim_f2w10.exe`** (bit-identique w0 sur 10 configs), maillages Delaunay
  Gmsh, référence R = dpr + dpApex + dpTension off + rankineDrive stress + meridian power + ω_c crack-band +
  spall sur l'énergie + soupape géométrique comptée + masque. Fumées : élastique 20 s, T3 confiné −99,98 MPa,
  R à P = 0 : 4,1 J plastiques / 0,21 J traction / 0,34 J ω_c, pic 55 kN, 0 érodé, 116 s (la grille en
  Rankine-déformation donnait 8,4 J et 39 kN : l'artefact de Poisson pesait ×2 sur la référence) ; tranche
  élastique 174 s pour 0,5 ms (≈ 20 min par coupe de 3,5 ms). Ordre C → D → E, ≈ 8 h.

- **07:20 — Bande de bruit T1 sur Delaunay (cœur 1,5 vs 1,0 mm)** : restitution −30 % (P = 0) / −16 %
  (P = 100), travail absorbé +40 % / +47 % ; dispersion inter-réalisations ±2 %. La règle pré-écrite
  (bande > 15 % → référence plus fine) s'applique : **référence T1 = cœur 1,0 mm** (88 k tets, 4-5 min/run),
  bande = cœur 0,75 mm (173 k tets) et cœur 1,5 mm (demi-poids). Matrice régénérée (44 decks), les
  4 runs déjà faits (R à P = 0/100 en 1,0 et 1,5 mm) sont recyclés sous leurs nouveaux noms ; le runner
  saute les runs déjà terminés. Premiers faits : à P = 100 la référence est presque élastique (e_r 0,75,
  4 J absorbés, aucune plasticité) alors que G_IIc = 1 N/mm (b_c ×10) donne e_r 0,12 et 14 J : **G_IIc est
  LE paramètre sensible sous confinement**, comme la critique adverse l'anticipait. V(D ≥ 0,9) rétrogradé
  au rang 2 : avec kf = 174 k0, D = 0,9 correspond encore à 95 % de ft portés (D mesure la raideur).
  Relance 07:2x sur w10.

- **07:55 → 08:20 — Convergence de la référence R à P = 0 (sondes `convergence/`)** : cœur 1,5 / 1,0 /
  0,75 / 0,5 mm → travail absorbé 6,7 / 9,5 / 12,0 / 12,9 J, restitution 0,58 / 0,41 / 0,25 / ≈0,19,
  pic 55 / 50 / 47 / 46 kN. La zone plastique sous l'insert (rayon de contact ≈ 3 mm) exige ≈ 6 éléments
  dans le rayon : la référence T1 doit passer au **cœur 0,5 mm** (maillage `T1_c05.msh`, 98 k tets,
  ≈ 10 min/run seul) ; sonde 0,35 mm en cours pour borner la bande. Coût matrice T1 ≈ 4 h.

- **08:25 — Référence T1 = cœur 0,5 mm** (`T1_c05.msh`, 98 k tets ; réalisations `T1_c05_s2/s3`, point fin
  `T1_c035` 124 k tets, point grossier `T1_h075`). Les runs à cœur 1,0 mm (référence, s2, s3, GIIc1) sont
  archivés sous `out_*_h10*` comme points de bruit documentés, hors verdict. À P = 100 la référence est déjà
  convergée à 1,0 mm (e_r 0,626 → 0,613 à 0,75 mm ; ±2 % inter-réalisations) ; à P = 0 il faut 0,5 mm.
  La partie COUPE de la phase C est lancée seule (`run_T2C.sh`, 08:33) pendant la sonde 0,35 mm ; la
  matrice complète repart ensuite avec le runner qui saute les runs faits.

- **08:59 — Sonde 0,35 mm rendue** : W_abs 14,0 J (0,5 mm : 12,9, +8 %), plastique 9,9 J (+10 %), pic 44,7 kN
  (−2 %), restitution ≈ 0,125 (0,5 mm : ≈ 0,19). Le travail absorbé et la partition convergent à < 10 % à
  0,5 mm ; la restitution (1 − W/E₀, petit nombre) reste sensible et sera exclue par la règle des 15 %.
  **Matrice lancée 09:00 sur w10, ordre C1 (percussion bruit) → D (percussion ablations) → C2 (coupe bruit)
  → E (coupe ablations)** ; coût : percussion ≈ 12 min/run (24 runs ≈ 5 h), coupe ≈ 40-60 min/run (20 runs,
  nuit). La tranche extrudée (Delaunay 2D × 4 couches) n'apporte rien (h_in min 0,11 contre 0,08 mm, ×1,5 de
  tets) : abandonnée, maillages `T2x_*` conservés pour mémoire.

## Matrice (GO de Fernando le 04/09 ; relancée 07:00 sur w10)

`matrice/MATRICE.md` : 42 decks (C bruit 17, D percussion 13, E coupe 12), runner `run_matrice.sh <exe> <C|D|E|all>`.
Coût estimé au débit mesuré : T1 ≈ 2-2,5 min (dx 1,5), ≈ 8 min (dx 1,0) ; T2 ≈ 10 min (146 k pas) ;
total ≈ 5 h de machine. À valider phase par phase (liste exacte des clés dans chaque `.cfg`).

- **01:57 → 02:12 — w8, clôture des bancs** : `phaseB/RESULTATS_phaseB.md`. Tout ce qui devait passer passe,
  tout ce qui devait échouer échoue ; deux lignes d'information (masque inerte sans orphelin, µ = 0,4 en
  élastique). La matrice (42 decks, `matrice/`) est prête pour `rockim_f2w8.exe`, **en attente du GO**.

## 2026-09-04 11:45 — figures 4-6 de la phase C (fig_matrice.py)
- fig4 convergence R (P 0/100, 1,5 → 0,35 mm + 2 germes à 0,5) ; fig5 F(t)/enfoncement/partition R P0, R P100, GIIc1 P100 ; fig6 coupes |y−24| < 0,4 mm (ω_c, ε̄ᵖ, noir = supprimés). Les polygones d'éléments supprimés dont les nœuds se sont envolés (arête > 2 mm) sont filtrés.
- FAIT : dans C_T1_R_P100_GIIc1, 2 656 érodés = 2 639 par la soupape géométrique (nEroGeo) + 17 par la loi (spall). La loi ne retire pas elle-même le matériau broyé (ω_c = 0,98, résiduel 2 %) : c'est la soupape detF/strain qui fait le cratère (65 mm³). À discuter avec `erodeDc` (existe, opt-in) lors des verdicts.
- Sens physique : à P = 100 la référence est quasi élastique (e_r 0,48, 8,3 J, F 64 kN) ; G_IIc 10 → 1 N/mm fait basculer en poinçonnement (F plafond ~8 kN, 1,6 mm et continue, e_r 0,07, wDamC 6,1 J).

## 2026-09-04 11:53 — PAUSE de la matrice (demande de Fernando : calcul important à lancer)
- Lanceur run_matrice.sh (PID 32396/1508) et run D_T1_noWc_P000 (rockim PID 35520, démarré 11:43) arrêtés par Stop-Process. 19 runs terminés (C1 : 15, D : 4 — R_P050, lin_P000, apex_P000). Le dossier out_D_T1_noWc_P000 incomplet sera écrasé à la reprise (run_matrice.sh ne saute que les runs avec « wall time » dans le log).
- Reprise = relancer run_matrice.sh DÉTACHÉ (même commande qu'au lancement, Start-Process bash), après accord de Fernando.

## 2026-09-04 12:20 — fig7 force–pénétration (fig_FP.py)
- 4 panneaux : R à P 0/50/100 ; ablations P = 0 (R, lin, apex) ; P = 100 R vs G_IIc 1 ; convergence R P = 0 (1,5 → 0,35 mm). Force = toolFz filtrée 10 µs (compression positive), pénétration comptée depuis le premier contact (|F| > 200 N).
- Lecture : à P = 0 la charge est concave (plastification + Rankine sous l'insert dès le début) et la boucle est large ; à P ≥ 50 la charge est quasi linéaire et la décharge presque parallèle (boucle mince, restitution 0,45-0,48). lin ≡ R sur toute la charge (courbure inerte) ; apex diffère seulement après 0,45 mm (pic +7 %, décharge plus élastique). G_IIc = 1 : pic à 20 kN puis plateau 7-12 kN sur 1,3 mm (poinçonnement). Maillage : la charge converge dès 0,75 mm, la décharge dès 0,5 mm.
- Le dossier out_D_T1_noWc_P000 interrompu a été renommé out_D_T1_noWc_P000_interrompu (rien supprimé) pour ne pas polluer resultats.csv/VERDICTS.md.

## 2026-09-05 00:05 — DÉFAUT TROUVÉ : broche fantôme sous l'insert (noeud orphelin des maillages T1)
- Tous les maillages T1 (h15/h10/h075/c05/c035 et germes s2/s3) contiennent UN noeud orphelin (jamais référencé par un
  tétraèdre) : le point géométrique du champ de taille Gmsh (24, 24, 32), EXACTEMENT au point d'impact. fem3d lui donne
  une masse nulle et l'épingle (FIXED) — mais toolContact lui appliquait la force de pénalité kp × pénétration : une
  broche rigide fixe sous le pôle de la sphère. Découvert le 05/09 via le port Signorini (NaN par division par sa masse).
- Ordre de grandeur (C_T1_R, kp = E h_min = 9,3e6 N/m) : broche max 5,9 kN à P = 0 (pic 45,6 kN, 13 %), 4,9 kN à
  P = 100 (pic 63,9 kN, 8 %), énergie élastique de broche ~1,9 J restituée (sur 16 J : pèse sur e_r). Les phases C1 et D
  (percussion) sont BIAISÉES par cette broche ; la coupe (T2, sans point de champ) ne l'est pas.
- Correctifs : garde `if (m_[i] <= 0) continue;` dans Fem3dSolver::toolContact (change les résultats des maillages à
  orphelin seulement) ; maillages nettoyés `*_clean.msh` (meshes/drop_orphans.py, check_orphans.py) ; clé
  `toolContact = signorini` portée dans fem3d (w13, banc de Hertz en cours). Décision Fernando : refaire C1 + D (28 runs,
  ~5 h) avec les maillages propres et le contact retenu.

## 2026-09-05 05:35 → 06:10 — w15 : viscosité de volume opt-in `bulkViscosity = b1 b2` (agent codeur, tour 1)
- Pourquoi : `audit_crush/SYNTHESE_audit_broyage.md` — sur le maillage Abaqus identique rockim porte 17-20 % de moins qu'Abaqus
  au pic et son bulbe broyé continue de croître après le pic (+43 % contre +9 %) ; suspect principal du post-pic = la viscosité
  de volume d'Abaqus/Explicit (b1 0,06, b2 1,2), absente de rockim.
- Code (`src/Fem3dSolver.cpp`, `include/rockim/Fem3dSolver.hpp`, originaux dans `bitid_w15/orig/`) : par élément
  ε̇_vol = (J_{n+1} − J_n)/(dt J_{n+1}) (J = det F déjà calculé, mémoire `Elem::Jbv`), L_e = lc = V0^{1/3}, c_d = cP() non
  endommagé ; p₁ = b1 ρ c_d L_e ε̇_vol (tout signe), p₂ = ρ (b2 L_e ε̇_vol)² en compression seule (ε̇_vol < 0) ; σ_bv = q I,
  q = p₁ − [ε̇<0] p₂, du signe de ε̇_vol (s'oppose au taux : en compression rapide ajoute de la compression), ajoutée à la
  contrainte de la LOI pour les forces internes seulement (svm/pm/szz/VTU/jauges = matériau, comme Abaqus) ;
  wBulk = Σ V0 q ε̇_vol dt ≥ 0 (résumé + dernière colonne de history.csv si actif, réduction OpenMP en ordre de threads) ;
  `computeStableDt` : CFL × (√(1 + b1²) − b1) (part linéaire de ξ ; la part quadratique −b2² L_e ε̇/c_d n'est pas portée,
  la borne du ressort de contact n'est pas touchée). Clé absente OU `0 0` = chemin inchangé (aucune colonne). Deux nombres
  exigés, ≥ 0, finis. Doc : §5.20, ¶ « Viscosité de volume » (formules, choix de L_e, dt, non porté).
- Bit-identité (`bitid_w15/comparaison_w15.txt`, `bitid_w15.sh`) : `PQ_cdpI_P020_court` sans clé, w14 vs w15, OMP 2 :
  history e45a7cbc0609c9f8, frames 5ddf25b84b9c15f1, vtu[4] df85a16d5b944875 pour les deux, cmp 6/6 identiques.
- Banc falsifiant (`bitid_w15/selftest_bv/`, `SELFTEST_bv.md`) : bloc élastique 6×6×24 mm, échelon −2 m/s, damping 0,
  kp/10 (la borne 2√(m/kp) est héritée en scénario tension et masquait le facteur sur dt). (i) `0 0` ≡ sans clé octet pour
  octet ; (ii) `0.06 1.2` : sonnerie de tête −23,2 % (tiers central −9,7 %, jauge moyennée), wBulk 8,1e-5 J, dt × 0,94180 ;
  (iii) `0.3 1.2` : −61,8 % / −46,3 %, wBulk 2,2e-4 J, dt × 0,74403 ; témoin sans clé au même dt : −6 % seulement (c'est la
  viscosité, pas le dt) ; (iv) facteurs de dt exacts à 5 décimales ; signe VÉRIFIÉ (wBulk > 0 dans les 8 runs actifs, sonnerie
  réduite et non amplifiée ; quadratique seul en rampe lente : traction 1,8e-38 J = 0, compression 5,4e-10 J > 0) ;
  (v) AUCUN scénario fem3d à volume constant (cisaillement pur / rotation rigide) : non couvert par un run, seulement
  pullV = 0 (8,9e-24 J) et la construction (q ne dépend que de det F). Percussion élastique : toolFz −9,5 %, dt borné par kp.
- Deck écrit, NON lancé : `perc3d/PQ2_cdpH_pulseZ2_abqmesh_bv.cfg` (= abqmesh + `bulkViscosity = 0.06 1.2` + `frames = 40`) ;
  marqueur `W15_READY` laissé au relecteur (la file de nuit `queue_nuit.sh` attend les deux). Rien touché dans perc3d/ sinon,
  w14 et la file de nuit intacts. Limites à relire : L_e = V0^{1/3} (Abaqus : longueur caractéristique du C3D4), part
  quadratique du dt non portée, pas d'équivalent ALLVD complet.

## 2026-09-05 06:05 → 06:40 — w16 : cinématique de Hencky opt-in `kinematics = hencky` (agent codeur, tour 1)
- Pourquoi : `cdp_rockim/audit_crush/SYNTHESE_audit_broyage.md`, constat 1 — sur le même maillage et la même carte CDP rockim
  porte 20 % de moins qu'Abaqus au pic ; la couche de contact est à J = 0,84-0,94 (20-30 % de déformation) où la mesure de
  Biot du fem3d (ε = sym(Rᵀ F) − I, P = R σ_c, référence) et la mesure logarithmique d'Abaqus/Explicit (Cauchy, forces sur la
  configuration courante) diffèrent de 5-19 % sur la force transmise.
- Code (`src/Fem3dSolver.cpp`, `include/rockim/Fem3dSolver.hpp`, nouveau `include/rockim/Fem3dKinematics.hpp`, originaux dans
  `bitid_w16/orig/`) : clé `kinematics = biot | hencky` (autre valeur refusée). hencky, par élément dans `elementForces` :
  décomposition spectrale de C = Fᵀ F (`SelfAdjointEigenSolver` 3×3) → U⁻¹ = Σ n nᵀ/λ_i, R = F U⁻¹ (polaire exacte, plus
  d'itération), ε = ln U = Σ ln λ_i n nᵀ passée à la loi (interface `MatLaw::stress` inchangée ; epsP devient logarithmique) ;
  σ_c rendue = Cauchy co-rotationnelle ; P = J R σ_c U⁻¹ (= J σ F⁻ᵀ, F⁻ᵀ = R U⁻¹), f_a = −V0 P ∂N_a/∂X avec les gradients de
  référence déjà stockés (≡ −V σ ∂N_a/∂x). Garde-fou det F ≤ 1e-9 ou λ_min ≤ 1e-9 → chemin biot de l'élément (comme aujourd'hui).
  `bulkViscosity` : q I ajouté à σ_c avant P dans les deux cinématiques. OpenMP : tout local à l'élément. Inchangés et
  documentés : masse lumpée V0, Lysmer, contact, lc = V0^{1/3}, CFL sur la géométrie initiale (non recalculée), `erodeDetMin`
  (det F), `erodeStrainMax` lit ‖ε‖ = ‖ln U‖ en hencky (0,3 ⇔ λ ≈ 0,74 au lieu de 0,70), sorties svm/pm/szz/slat/VTU/jauges
  `sigZZmid` etc. = contrainte de la loi donc Cauchy en hencky, `epsAxMid`/`epsVolMid` = ln U (ε_vol = ln J), wEl = ½ σ_c : (ε − ε_p)
  par volume de référence (nominal). Approximation documentée : ln U TOTAL (hyperélastique) contre le taux intégré d'Abaqus —
  identiques pour un trajet coaxial (uniaxial, pôle de l'insert), différents en cisaillement fini.
- Bit-identité (`bitid_w16/comparaison_w16.txt`, `bitid_w16.sh`) : `PQ_cdpI_P020_court` sans clé, w15 vs w16, OMP 2 :
  history e45a7cbc0609c9f8, frames 5ddf25b84b9c15f1, vtu[4] 180a5a99311bbe44 pour w15, w16 ET w16 + `kinematics = biot`
  (= les hachages history/frames de comparaison_w15.txt : chaîne w14 → w15 → w16 bit-identique), cmp 6/6 identiques pour
  les deux paires, peak tool force 7791,04 N dans les trois logs (123-132 s chacun).
- Banc à formes fermées (`bitid_w16/selftest_hencky/`, `SELFTEST_hencky.md`) : bloc élastique 4×4×8 mm (768 tets), E 77,66 GPa,
  ν 0,29, scénario tension, `gripLateralFree = true` (vérifié dans `integrate()` : mors z = 0 et z = H libres en x, y ; ni Lysmer ni
  confinement en tension → faces libres, état uniaxial homogène : sigZZmid = σ_mors à 1,0000 en biot), rampe 200 µs puis ±1,8 m/s,
  Cundall 0,7, quasi statique (ρcv = 29 MPa contre 15,5 GPa), λ = 1 + epsAxGrip, σ = F_mors/A0 lissée 21 lignes et interpolée à
  la cible. (a) λ = 0,8 et 1,2 : biot à −0,01 % de E(λ−1), hencky à −0,01 % / −0,00 % de λ^(−2ν) E ln λ (RMS 0,02 % sur le
  trajet) ; mauvaise forme rejetée : hencky contre la forme biot +26,98 % (0,8) / −17,99 % (1,2), biot contre la forme hencky
  −21,3 % / +21,9 % ; rapport hencky/biot mesuré +26,99 % = forme fermée. Jauge du tiers central en hencky : |sigZZmid| = E ln λ
  à −0,01 % (Cauchy) et σ_mors/|sigZZmid| = 1,1382 = λ^(−2ν) (aire courante) : deux mesures indépendantes de P = J R σ U⁻¹.
  (b) rotation rigide : aucun scénario fem3d ne la permet → test unitaire `test_kinematics.cpp` compilé à part (23/23) :
  F = Q → ln U = 0, R = Q, P = 0 à 1e-16 ; F' = Q F → ε inchangée, R' = Q R, P' = Q P à 1e-15 ; formes fermées à 2e-16 ;
  det F ≤ 0 → false. Constat annexe : le R itéré du chemin biot historique n'est exact qu'à ~1e-7 à 20 % de déformation.
  (c) λ = 1,002 : hencky/biot −0,215 %, chacun à −0,005 % de sa forme. (i) clé absente ≡ `kinematics = biot` octet pour
  octet (history, frames, 3 vtu) et w15 ≡ w16 sans clé sur le banc ; `kinematics = log` refusé (rc 1). Coût ×1,7 en élastique
  (spectrale), dt inchangé.
- Deck écrit, NON lancé : `perc3d/PQ2_cdpH_pulseZ2_abqmesh_hencky.cfg` (= abqmesh + `kinematics = hencky` + `frames = 40` + en-tête) ;
  marqueur `W16_READY` laissé au relecteur. Rien touché dans perc3d/ sinon, w14/w15 et la file de nuit intacts (w14 tournait
  sur 14 threads pendant tout le tour ; tests à OMP 2). Limites à relire : ln U total vs taux intégré (non coaxial), lois plastiques
  en hencky non testées sur banc (epsP logarithmique), CFL non recalculée sous forte compression, erodeStrainMax sur ‖ln U‖.

## 2026-09-05 10:31 — MATRICE v2 LANCÉE (percussion d'abord), décision Fernando « fais tout de 1 à 4 »
- Enseignements de la nuit (cdp_rockim/JOURNAL.md) reportés sur l'étude : broche fantôme (maillages `*_clean`, Signorini), viscosité
  de volume + cinématique de Hencky (25 % du pic, tout le post-pic), suppression = mécanisme du verdict poinçonnement/rebond.
- `matrice_v2/` : 49 decks (44 v1 + 5 témoins `_noero`), liste exacte des clés dans `MATRICE_v2.md` (validée par Fernando ce matin),
  vérifiée deck par deck (diff v1 → v2, 0 orphelin, chargement w16). Coût : percussion C1 + D 10,2 h mesurées ; coupe C2 + E 30-38 h
  ESTIMÉES (aucun T2 n'a jamais tourné).
- File `matrice_v2/queue_matrice.sh` (détachée, w16, 14 threads, log `queue_matrice.log`) : C1 (11) → D (16) → `C_T2_R_P000` seul
  pour mesurer le coût réel des T2 ; la phase E sera dimensionnée sur ce coût (si > 1 h/run : ablations à P = 0 et 100 seulement,
  doublon `E_T2_R_P100_noero` ≡ `C_T2_R_P100_soupOff` retiré).
- Dépouillement v2 en cours d'écriture (agents) : `depouille_v2.py`, `verdicts_v2.py`, `fig_matrice_v2.py` — verdicts sur l'énergie
  absorbée, le cratère et le bulbe avec témoins de suppression, plus sur la force maximale seule.

## 2026-09-05 10:50 → 11:20 — w17 : outil bloqué en x, y `toolLockXY = true` (agent codeur, tour 1)
- Pourquoi : contre-expertise `perc3d/DECHARGE_abaqus_vs_rockim.md` § 6 (anomalie relevée) — dans le quart de bloc l'insert
  posé sur l'arête (0, 0) dérive vers les plans de symétrie (toolX −0,055 mm au pic, −0,18 mm à 320 µs, −0,81 mm à 800 µs) :
  le quart de sphère ne reçoit qu'un quart de la pression de contact, rien n'équilibre F_x, F_y et la vitesse x, y de l'outil
  rigide est libre ; Abaqus bloque le nœud de référence du bouton (`*Boundary RPB1 1,2`).
- Code (`include/rockim/Fem3dSolver.hpp`, `src/Fem3dSolver.cpp`, originaux dans `bitid_w17/orig/`, patch exact
  `bitid_w17/patch_w17.py`) : `Tool3::lockXY` (false) lu dans `placeTool()` par `toolLockXY` ; `Tool3::integrate` : après
  `v += dt/m F`, `v.x = v.y = 0` si actif, donc `x += dt·0` — x, y EXACTEMENT fixes ; la force de contact latérale reste
  calculée par `toolContact()` et comptée (`toolFx`, `peakF_`) ; `work_` = −F·v dt n'en reçoit rien (RP bloqué). Scénario
  shear : clé IGNORÉE avec `[FEM3D] WARNING` (outil piloté en déplacement, v_x = cutSpeed). Message `[FEM3D] toolLockXY : …`
  au démarrage. Ajout console seulement (aucun fichier) : ligne de résumé « tool lateral » (x, y de fin, max |dx|, |dy|,
  |Fx|, |Fy| ; membres `toolX0_`, `toolDxMax_`…), toujours imprimée en percussion/shear — l'observable de la dérive.
- Build COMPLET en place (Remove-Item *.obj, 13 .obj recompilés, `build_f2.cmd rockim_f2w17.exe` = vcvars64 + cl /O2 /openmp,
  `build_w17.log`, 11:03) ; `rockim_f2w16.exe` intact (06:13, sha 569d6b9a122207dc), matrice v2 non touchée (14 threads pendant
  tout le tour ; tests à OMP 2).
- Bit-identité (`bitid_w17/comparaison_w17.txt`, `bitid_w17.sh`) : `PQ_cdpI_P020_court` sans clé, w17 vs sorties conservées de
  w16 (même deck, même OMP 2, 06:18 — w16 non relancé) : history e45a7cbc0609c9f8, frames 5ddf25b84b9c15f1, vtu[4]
  180a5a99311bbe44, cmp 6/6 IDENTIQUE, pic 7791,04 N ; chaîne w14 → w15 → w16 → w17 bit-identique sans clé.
- Banc falsifiant (`bitid_w17/selftest_lock/`, `SELFTEST_lock.md`, `analyse_lock.py`, PASS 5/5) : même deck court + `frames = 13`,
  sans clé / avec clé (128 / 141 s). Sans clé : toolX ≠ 0 sur 910/2102 pas dès le premier contact (73,8 µs), −1,150 µm en x et
  −1,143 µm en y à 130 µs (57 µs de contact, croissance ~t³ : cohérent avec −55 µm au pic du run complet). Avec la clé :
  toolX = `0` sur 2102/2102 pas, toolY = `0` sur 15/15 images (test exact sur la chaîne), max |Fx| = 1476 N / |Fy| = 1410 N
  toujours calculés et comptés (1344 / 1016 à outil libre), message présent seulement avec la clé. Contrôle shear
  (`shear_lock_warn.cfg`, mini grille 6×6×4 élastique, lame, 0,15 s) : WARNING imprimé et toolX final = −1e-4 + 4 × 5,0542e-6
  exactement (clé ignorée).
- **Constat à relire** : la clé n'est PAS neutre sur le deck court — pic |Fz| brut −11 %, filtré 10 µs −14 % (7361 → 6352 N),
  travail −6 % (1,506 → 1,417 J) à enfoncement final égal (0,06 µm d'écart), alors que la dérive de position n'y est que de
  1 µm. L'effet vient de la VITESSE latérale de l'outil libre (F_lat ≈ 1 kN ≈ 0,17 F_z sur 0,125 kg : 10⁻² m/s en oscillation),
  qui entre dans `vrel = v_i − tool.v` de chaque nœud de contact et fixe la direction du frottement au pôle
  (`contactVreg` = 10⁻³ m/s) ; avec la clé elle est nulle, comme pour le RP d'Abaqus. Sur 800 µs (−0,81 mm) la géométrie
  compte aussi. Mécanisme à confirmer par le relecteur (par ex. `contactMu = 0` : l'écart doit tomber).
- Deck écrit, NON lancé : `perc3d/PQ2_cdpH_pulseZ2_abqmesh_bvh_nrb.cfg` = `_hencky` (bv + hencky) + `absorbing = all` + `toolLockXY
  = true`, `bottomFree` retiré (sans objet : avec `all` le fond est absorbé, la branche bottomFree n'est lue que sinon), en-tête
  explicatif. Vérifié dans `Fem3dSolver::setupBoundaries` (branche box) : la face supérieure z = H n'est classée sur aucun axe
  (nAxis = −1 → aucun Lysmer sur le dessus) et `quarterModel` saute les faces x = 0 et y = 0 → Lysmer sur le fond z = 0 et les
  deux faces extérieures x = W, y = D seulement, ρ c_p A/3 par nœud et par face en normal (ρ c_s A/3 en tangent ; Abaqus : une seule
  valeur 1,76 N s/mm par nœud ≈ ρ c_p A, donc vraisemblablement normaux seuls — direction à vérifier dans le .inp), `absorbSpringFactor = 0` → pas de ressort — l'emprise des 492 dashpots NRB du deck Abaqus. VÉRIFIÉ sur le
  maillage (`bitid_w17/check_nrb_mesh.py` → `.txt`, coordonnées exactes, tol 1e-9 = 1e-6) : fond 157 nœuds, x = W 186, y = D 186,
  soit **492 nœuds distincts = les 492 dashpots** ; dessus (2120 nœuds) et plans de symétrie (2134 / 2135) non absorbés.
  Rien touché dans perc3d/ sinon, ni dans matrice_v2/. Doc §5.20 ¶ « Outil bloqué latéralement ».
- 11:20 — dépouillement v2 livré et vérifié (recalcul indépendant de C_T1_R_P000 v1/v2 à < 1e-3) : `matrice_v2/depouille_v2.py`
  (112 colonnes : énergie, partition avec wBulk, cratère loi/soupape, bulbe d_c/D depuis history + vtu, coupe T2 sur la fenêtre
  x ∈ [4, 12] mm), `verdicts_v2.py` (règle proposée : inerte si |Δ| ≤ max(bruit phase C, 10 %) ; « sous la soupape » si |Δ| ≤ écart
  R vs témoin `_noero` ; COMPTE au-delà ; 0 → x = « apparition, à juger » — à valider par Fernando), `fig_matrice_v2.py` (5 figures,
  Computer Modern via usetex). T1 = bloc COMPLET (pas de ×4), T2 = tranche symmetryY. VERDICTS_v2.md se remplit avec la phase D.

## 2026-09-05 13:35 → 14:00 — w18 : sondes de point matériel `probes = x,y,z ; …` (agent codeur, tour 1)
- Pourquoi : Fernando veut suivre des POINTS MATÉRIELS du bloc (sous l'insert, à côté, en profondeur, au bord) et tracer leurs
  trajets contrainte-déformation pour voir quels mécanismes s'activent selon la loi (cône, cap, Rankine, ω_c, tables cdp).
- Code (`include/rockim/Fem3dSolver.hpp`, `src/Fem3dSolver.cpp`, originaux dans `bitid_w18/orig/`) : clé opt-in fem3d
  `probes = x,y,z ; x,y,z ; …` (mètres, repère du bloc recalé [0,W]×[0,D]×[0,H], dessus z = H). `setupProbes()` (appelé dans
  `init()` après `computeStableDt`) : analyse stricte (virgule décimale = erreur nommée), localisation barycentrique par les
  gradients dN déjà stockés (tol 1e-9, premier tétraèdre trouvé), sinon centroïde le plus proche + `[FEM3D] WARNING: probe k
  lies in no tetrahedron (OUTSIDE the block …)` ; plusieurs sondes par élément (un CRENEAU par élément distinct, `probeSlot_`,
  message « share element »). `elementForces` : `prb = !probes_.empty()` — sans clé UN booléen par appel, aucun test par élément ;
  avec, `probeRec` dépose σ de la LOI (avant viscosité de volume, = celle des forces) et ε passée à la loi (Biot / ln U) dans le
  créneau (un thread par élément : pas de course) ; érosion : σ = 0, ε = dernière valeur. `historyRow` (const) appelle `probesRow`
  via `std::unique_ptr<std::ofstream>` : `<outputDir>/probes.csv` à la cadence de history.csv, précision 9, flush par ligne.
  En-tête `t` + par sonde `p<k>_{elem, sxx syy szz sxy syz sxz, exx … exz (tensoriels), p (compression > 0), q, s1 ≥ s2 ≥ s3,
  trEpl, eplEq = √(2/3 dev:dev), epvEq (loi), dt = st.D, dc = st.Dc, pc (0 sans cap), detF, eroded}` + `epsTpl epsCpl dcdp` si cdp,
  `Dv1..3` si dpdfh ; champs sans objet = 0. Repère co-rotationnel (documenté : p, q, s_i invariants).
- Build COMPLET en place (Remove-Item *.obj, 13 .obj, `build_f2.cmd rockim_f2w18.exe`, `build_w18.log`, 13:45) ; w16 (569d6b9a…)
  et w17 (81334bf4…) intacts, matrice v2 (14 threads) et perc3d (4 threads) non touchés ; tous les tests à OMP 2.
- Bit-identité (`bitid_w18/comparaison_w18.txt`, `bitid_w18.sh`) : (a) `PQ_cdpI_P020_court` sans clé, w18 vs sorties conservées de
  w17 : history e45a7cbc0609c9f8, frames 5ddf25b84b9c15f1, vtu[4] 180a5a99311bbe44, cmp 6/6 IDENTIQUE, pic 7791,04 N — chaîne
  w14 → w18 ; (b) même deck + `probes` (5 sondes, `PQ_cdpI_P020_court_probes.cfg`) : mêmes hachages, cmp 6/6, un seul fichier en
  plus (probes.csv, 2103 lignes, 146 colonnes) ; 136 s / 136 s.
- Banc falsifiant (`bitid_w18/selftest_probes/`, `SELFTEST_probes.md`, `analyse_probes.py` → `resultats.txt`, **PASS 14/14**) :
  (i) bloc élastique 8³ mm sous confinement isotrope 30 MPa (tension `pullV = 0`, `pullDelay > T`, `topPressure = confiningPressure`,
  `gripLateralFree`) : p = 29,990 MPa (−0,03 %), q = 5·10⁻¹¹ MPa, ε_ii = −P/3K à 0,03 % ; contrôle QUI DOIT ÉCHOUER sans
  `topPressure` : p = 20,0 = 2P/3, q = 30,0 = P, rejeté ; (ii) compression élastique 4×4×8 `triaxStats` : `p1_szz` = `sigZZmid` à
  0,21 % max (0,01 % moyen) sur 2415 lignes (P2 0,09 %) ; (iii) sonde à (4, 4, 12) mm hors du bloc : WARNING, élément le plus
  proche dans la couche du dessus (z_c = 7,75 mm) ; (iv) trois sondes à 10⁻⁷ m : même élément 1752, 26 colonnes identiques sur
  3335 lignes.
- Lecture du deck court cdp avec sondes (quart de bloc, pôle = arête (0, 0), 57 µs de contact ; figures
  `bitid_w18/figures_probes/fig_probes_P1..5.pdf`) : P1 (pôle −1 mm) p max 383 / q max 698 MPa, s3 −840 MPa, ε_c^pl 1,04 %,
  d_c 0,92 = plafond de la table, det F 0,926, chute post-pic visible sur σ_zz-ε_zz (Biot, deck sans hencky) ; P2 (−5 mm) q 391,
  d_c 0,37 ; P3 (lèvre) q 40 MPa, excursions en traction σ_zz +7 MPa, d_t 0,008 naissant ; P4 (bord du bulbe) et P5 (champ
  lointain) élastiques (q 70 / 24 MPa). Le trajet NOMINAL passe sous le méridien EFFECTIF cdp (trait fin) : attendu, (1 − d) σ̄.
- Livrables matrice : `matrice_v2/probes_D.txt` = la ligne `probes = 24e-3,24e-3,31e-3 ; 24e-3,24e-3,27e-3 ; 28e-3,24e-3,31.5e-3 ;
  32e-3,24e-3,29e-3 ; 42e-3,24e-3,26e-3` (P1 pôle −1 mm, P2 pôle −5 mm, P3 lèvre r 4 / −0,5, P4 bulbe r 8 / −3, P5 lointain
  r 18 / −6 ; bloc complet 48×48×32, insert (24, 24, 32)) — decks D NON modifiés (relecteur après validation) ;
  `matrice_v2/fig_probes.py` (une figure par point, 3 panneaux : p-q + méridien du deck — dpr linéaire / puissance, apex, cap
  `capP0` ; cdp effectif ; σ_zz-ε_zz ; d_t, d_c, ε_pl depuis le contact toolFz > 200 N — 3 runs max, couleurs C0/C3/C2, PDF + PNG
  Computer Modern usetex), testé sur deux probes.csv synthétiques (dpr puissance + cap vs linéaire, érosion) et sur le banc cdp ;
  usage `python fig_probes.py D_T1_R_P050 D_T1_lin_P050 D_T1_apex_P050` (runs AVANT `--cfg`). Doc §5.20 ¶ « Sondes de point
  matériel ».
- Limites à relire : composantes dans le repère co-rotationnel (pas celui du bloc) ; ε_pl équivalente = norme de von Mises de
  ε_pl (la cdp n'a pas de scalaire unique : ε_t^pl et ε_c^pl sont donnés à part) ; méridien cdp tracé en contrainte effective
  face à un trajet nominal ; localisation O(n_el × n_sondes) à l'init (négligeable, 98 k tets × 5) ; rien lancé dans matrice_v2
  ni perc3d.

## 2026-09-05 16:55 → 17:30 — w19 : endommagement de traction à direction figée `tensionDamage = fixed` (agent codeur, tour 1)
- Pourquoi (demande de Fernando, 05/09) : le d_t de Rankine est un scalaire qui dégrade toute la partie tendue du tenseur — un élément
  fissuré en x perd aussi sa raideur en y. Méthode de la direction figée à l'amorçage (fixed crack model, Rashid 1968 ; Rots &
  Blaauwendraad 1989) appliquée au noyau continu : endommagement PAR DIRECTION, repère figé au premier dépassement de f_t, raideur perdue
  normalement au plan de fissure et conservée parallèlement ; cinétique = la bande de fissuration de Rankine existante (f_t, G_f, l_c) par
  direction, pas l'obscuration ; géométrie copiée de `dpdfh` (D_i dans un repère `eul` figé).
- Code (originaux dans `bitid_w19/orig/`) : clé opt-in `tensionDamage = scalar` (défaut, bit-identique) `| fixed` + `tensionShearRetention = β`
  (défaut 1) pour dpr/saksala ; refusée pour saksala2011 (port fidèle, endommagement piloté par la déformation viscoplastique, pas de bande),
  cdp, dpdfh. `MatState::Fcm {nAct, d[3], kap[3], eul[3]}` ; `PlasticDamageLaw::fixedCrackUpdate` : tenseur pilote T = σ_eff/E (stress) ou
  ε − ε_p (strain), amorçage quand la valeur propre max de T dans le COMPLÉMENT ORTHOGONAL des directions figées dépasse k_0 (1re : repère
  principal complet, main droite, Euler ZYX comme `dfhk::eulr` ; 2e : rotation de (n2, n3) autour de n1, θ = ½ atan2(2 T_23, T_22 − T_33) ;
  3e : n3), croissance κ_i = max(κ_i, n_iᵀ T n_i), d_i = 1 − k_0/κ_i exp(−(κ_i − k_0)/k_f) = la formule scalaire ; nominal S = Rᵀ σ R,
  S_ii (1 − d_i) si S_ii > 0 (unilatéral), S_ij (1 − β max(d_i^ouv, d_j^ouv)) (β = 1 = min(f_i, f_j) de la VUMAT DP-DFH), σ = R S_n Rᵀ ;
  ω_c sur la partie spectrale négative soustraite après ; wDamT += ½⟨S_ii⟩²/E dd_i (trois directions) ; s.D = max d_i, s.kappa = max κ_i.
  Fem3dSolver : `fixedDam_`, champs VTU `dFix1..3` (branche ajoutée, stats/vtkCap honorés), colonnes probes `dFix1..3`. main :
  `rockim selftest-fixed [out.csv]`. Doc §5.20 ¶ « Endommagement de traction à direction figée » ; fiche `MatLaw.hpp`.
- Build COMPLET en place (Remove-Item *.obj, 13 .obj, `build_f2.cmd rockim_f2w19.exe`, `build_w19.log`, 17:15) ; w18 (13:45) intact,
  matrice v2 (14 threads) non touchée, rien lancé dans matrice_v2/ ni perc3d/ ; tous les tests à OMP 2.
- Bit-identité (`bitid_w19/comparaison_w19.txt`, `bitid_w19.sh`, `bitid_w19_run.log`) : (a) `PQ_cdpI_P020_court` (cdp) sans clé, w19 vs
  w18 conservé : e45a7cbc0609c9f8 / 5ddf25b84b9c15f1 / 180a5a99311bbe44, cmp 6/6, pic 7791,04 N — chaîne w14 → w19 ; (b) deck COURT dpr
  `bitid_w19/C_T1_R_P000_court_grid.cfg` (= C_T1_R_P000 à T 100 µs, 2 images, maillage interne 24 × 24 × 16 mm, 55 296 tets — le noyau
  touché) sans clé, w18 vs w19 : e23cbec3a88c2e6f / 08d05ce828a8e2af / bd067aa893e231ac, cmp 6/6, pic 24 623,9 N, 106 / 102 s ;
  (c) `tensionDamage = scalar` explicite ≡ sans clé, cmp 6/6 ; (d) contrôle `fixed` : 6/6 DIFFÉRENT, pic 23 348,9 N (−5,2 %), wDamT 0,166 vs
  0,122 J, 2686 érodés (canal erodeWfrac, wDamT sommé sur trois directions) contre 0, damage ≥ 0,9 sur 17 187 vs 23 234 tets, dFix1 > 0 sur
  50 429, deux directions 33 094, trois 8288, `damage` = max(dFix) à 0,0 ; 76,6 s (rotation R S Rᵀ moins chère que le spectral).
- Banc falsifiant (`bitid_w19/selftest_fixed/`, `SELFTEST_fixed.md`, `fixed.csv`, **33 PASS / 0 FAIL**) au point matériel, carte Bohus dpr
  (dpApex, dpTension off, rankineDrive stress), pilotage en déformation ε = e n nᵀ − ν e (I − n nᵀ) : (0) 5 clés invalides refusées,
  sans clé ≡ `scalar` explicite (memcmp 304 pas, 0 composante différente) ; (a) traction x à d = 0,9, décharge, traction y : pente
  σ_yy/ε_yy = **1,000 E** (fixed, PASS) contre **0,0997 E = (1 − D) E** (scalar, REJETÉE par le critère E à 1 %) — le test qui sépare ;
  même cinétique uniaxiale (5·10⁻¹⁶), d1 inchangé bit pour bit ; y s'amorce à E ε_yy = f_t (1,02 k_0), n2 = y ; (b) compression x
  après fissure x : E dans les deux (−27,000 MPa) ; (c) 45° : eul (45,000, 0, 0)°, n1 = (0,7071, 0,7071, 0), objectivité 2,6·10⁻¹⁵ ;
  (d1) tel qu'écrit (x maintenu, ε_yy croissant) : σ_xy = 0 par symétrie, y s'amorce à S_22 = f_t, d1 0,5 → 0,763 (Poisson) ;
  (d2) DONNÉE rotation (x maintenu à d 0,5 + traction selon 45°) : σ_xy^nom/σ_xy^eff = 0,242 (β = 1 = 1 − max(d1 0,758, d2 0,505)) et
  0,924 (β = 0,1), désalignement de la principale nominale 13,9° / 13,4° vs effective (stress locking de Rots ; scalar 0,183) ;
  (e) ∫σ dε × l_c/G_f = 0,9999, wDamT × l_c/G_f = 1,0008, les deux modèles.
- Livrable matrice : `matrice_v2/D_T1_fixed_P000.cfg` = C_T1_R_P000 + `tensionDamage = fixed` + la ligne `probes` de `probes_D.txt`,
  NON lancé (binaire w19 obligatoire : un binaire antérieur ignore la clé en silence).
- Points à relire : (1) wDamT somme les trois directions → `erodeWfrac` se déclenche plus tôt (jusqu'à 3 G_f/l_c par élément) — garder ou
  passer sur max_i w_i ; (2) saksala2011 exclu (la spec le nommait) ; (3) le (d) « tel qu'écrit » donne σ_xy = 0 par construction, la
  donnée de stress locking vient de la variante (d2) ; (4) eul2/eul3 d'une fissure unique non contraints (plan propre dégénéré) ;
  (5) pilotage `rankineDrive = strain` de la branche figée non couvert par le banc.

## 2026-09-05 18:25 — matrice v2 : PERCUSSION COMPLÈTE (C1 11 + D 16 runs, w16/w18, sondes sur D)
- Coût mesuré : ×1,5 à ×1,7 de la v1 (bv + hencky + partage des cœurs avec les agents) ; 27 runs en 7 h 50.
- Résultats (résumé, détail `matrice_v2/resultats_v2.md`, `VERDICTS_v2.md`, figures/) : à P = 0 seule la coupure de traction compte
  (apex 23 MPa : +22 % de pic, e_r ×2, V(d_t > 0,9) ÷2,7 ; lin, noWc, cap inertes à 3 %) ; dès 50 MPa la traction s'éteint (V_D09
  11 863 → 140 mm³) et le méridien devient la seule brique active (lin : +13 % à 50, +15 % à 100, W −26 %, e_r 0,52 vs 0,34) ;
  cap de compaction −5 % de pic / +8 % de W à 100 (second ordre) ; boue inerte ; noWc +4 % ; G_IIc = 1 = poinçonnement par suppression
  (2 766 él., 14 kN). d_c < 0,5 partout sauf G_IIc 1 : avec G_IIc = 10 N/mm le continu ne broie pas. Témoins sans suppression = R
  à 0,1 % aux trois P (R ne supprime rien). Sondes : sous le pôle 2-3,5 GPa axial, 8-12 % de déformation, d_c ≤ 0,3 ; à P = 0 d_t
  précoce sous le pôle (traction de Hertz sous le contact) ; lèvre rompue en traction à q < 40 MPa.
- Suite de file : C_T2_R_P000 (coût T2) 18:24 → références sondées P0/P100 (w18) → D_T1_fixed_P000 (w19).

## 2026-09-05 19:00 — dépouillement v2 de la percussion (agent dépouilleur) : VERDICTS_v2.md livré, deux règles côte à côte
- `depouille_v2.py` relancé avec VTU sur les 27 runs T1 (+ `C_T2_R_P000` en cours, non jugé) ; `verdicts_v2.py` réécrit : règle **A** a priori de la
  PROPOSITION (20 % et > bruit, inerte < 10 %, exclue si bruit > 15 %, rang 1 = e_r, δ_max, δ_res, partition, volumes) ET règle **B** proposée
  (seuil max(bruit, 10 %), « sous la soupape », rang 1 = W/e_r, cratère, V_ωc05) → `VERDICTS_v2_tables.md` (annexe) ; `VERDICTS_v2.md` = livrable
  style article, règles présentées comme NON validées.
- **Témoins `_noero` bit-identiques à R** aux trois P (SHA-256 des history.csv : b846d029… / 3a9dba07… / 5896c123…) : R ne supprime rien, la classe
  « sous la soupape » est vide ; ils servent de référence SONDÉE à P = 0 et 100 (les `_probes` doivent leur être identiques : table « Répliques »).
- **Résiduel négatif corrigé** : c'était le VOL LIBRE de l'outil après séparation (e_r 0,3 → 4,5-4,7 m/s × 110-190 µs = 0,5-0,7 mm), pas le
  soulèvement (Poisson +0,024 mm à 100 MPa, mesuré dans les VTU). Nouvelle convention δ_res = d_sep + [dz_surf(t_sep) − dz_surf(t_c)] (surface
  libre à r > 18 mm dans les VTU : dérive rigide du bloc sur les Lysmer −0,057/−0,067/−0,072 mm retirée) : R = 0,409 / 0,240 / 0,210 mm ;
  colonnes `d_sep_mm`, `d_fin_mm`, `vol_libre_mm`, `derive_bloc_mm`, `soulevement_mm`, `history_sha256` ajoutées.
- Verdicts (A / B) : lin inerte à 0, COMPTE à 50 et 100 (W −19/−26 %, e_r +41/+50 %, δ_res −44/−58 % ; sondes : q 5,3-6,2 GPa au pôle contre
  2,5-2,7, ε_pl 2-3 % contre 7,8 — c'est l'extrapolation hors données qui compte) ; apex COMPTE à 0 (F +22 %, e_r ×2, wDamT −45 %, V_D09 ÷2,7),
  inerte à 50/100 (A dit COMPTE à 50 sur wDamT +31 % = 2 mJ : artefact d'une composante sans poids) ; noWc inerte partout (d_c ≤ 0,29, wDamC 2-4 %
  de W : la prédiction « dominant à 100 MPa » est falsifiée en percussion) ; cap inerte à 0, COMPTE à 100 (δ_res +28 %, wPlas +21 %, e_r −16 %) ;
  boue inerte (A : wDamT −76 % = 6 mJ) ; G_IIc 1 = poinçonnement, run INVALIDE (eRemoved 22 % de W, contact non fini à T).
- Bruit phase C : W 2,5 %, F 1,6 %, δ_res 7,6 %, e_r 22,7 % à P = 0 (maillage grossier), wDamC 47-52 %, V_D09 17 % à 100.
- Figures : `fig_v2_6_pic_energie_briques` (synthèse barres F/W/e_r/δ_res × brique × P), partition vérifiée (Σ = 16 J, confWork constant pendant
  le contact), sondes par brique `figures/probes_P<P>_<groupe>/fig_probes_P1..3` (`make_probes_figs.py`, R + 2 briques, légendes sous les panneaux).
- Cohérence : 3 runs au hasard recalculés hors depouille (écart 0), hachages égaux.
- En attente : `D_T1_R_P000_probes`, `_P100_probes`, `D_T1_fixed_P000` (file après C_T2, ≈ 1 h 50/run de coupe) — les scripts les intègrent à la
  relance ; coupe non jugée. Rien lancé, decks et out_* intacts.

## 2026-09-05 17:40 → 19:00 — w20 : gardes des entrées C3 / C4 / C6 (agent codeur, tour 1)
- Pourquoi (plan de robustesse du 05/09, chantier C ; constats de l'état des lieux) : import Gmsh fem3d gardait un nœud orphelin en
  silence (masse 0, épinglé FIXED = la broche fantôme sous le pôle d'impact) ; « degenerate tet » sans identifiant ; détecteur de
  NaN AVEUGLE sur `u_[0]` en fem3d, échantillon 1/256 en fdem 2D/3D, rien en fem/dem/dem3d ; `hydro` accepté en 3D sans un mot ;
  clés d'un autre mode ignorées en silence. Aucune décision à prendre : aucun comportement légitime à préserver.
- Code (originaux `bitid_w20/orig/`, patch rejouable `bitid_w20/patch_w20.py`, 67 remplacements) : nouveau `include/rockim/Guards.hpp`
  (checkOrphans, checkMasses, checkDegenerate/degenerateError, checkFinite, exception `NanError`) ; `Config::keys()` ;
  `include/rockim/KeysByMode.hpp` GÉNÉRÉ par `tools/gen_keys_by_mode.py` avec le registre `tools/keys_by_mode.json` (méthode dans
  la docstring et DOC §8.10) ; `main.cpp` : garde `hydro*` hors fdem, `keysbymode::check`, `catch (NanError)` → `ERROR.txt` + code 3 ;
  six solveurs : `nanEvery_ = geti("nanCheckEvery", 256)`, `checkFinite()` à cadence fixe + dans `finalize()`, contrôles de maillage
  après lumping. Grille fem3d `geometry = cylinder` : nœuds hors cylindre comptés/imprimés, toujours épinglés (seule exception).
- Build COMPLET en place (Remove-Item *.obj, `.\build_f2.cmd rockim_f2w20.exe`, `build_w20.log`, 18:08 ; `NoDefaultCurrentDirectoryInExePath`
  oblige le `.\`) ; w18/w19 intacts, matrice v2 non touchée, rien lancé dans matrice_v2/ ni perc3d/, tests à OMP 2, suite à OMP 1.
- Bit-identité w19 = w20 (`bitid_w20/bitid_w20.sh`, `comparaison_w20.txt`, OMP 2, cmp octet par octet) : dpr court
  `C_T1_R_P000_court_grid` e23cbec3a88c2e6f / 08d05ce828a8e2af / bd067aa893e231ac, 6/6 IDENTIQUE, pic 24 623,9 N ; fdem 2D
  `fdem2d_court` (= verify_fdem_voronoi_tension) 0c3f2f1ec938759d / c06c08e9385b7766, 28/28 ; fdem3d `fdem3d_court` (= verify_fdem3d_tension
  à T 100 µs) 4994fe26a0164531 / 7e6ec512835aac19, 9/9 ; **cdp `PQ_cdpI_P020_court` REFUSÉ par w20** : `perc3d/Q1_c05.msh` porte
  l'orphelin (nœud 9, (0, 0, 32) mm = l'arête d'impact du quart de bloc) — la chaîne w14 → w19 (e45a7cbc…) tournait AVEC la broche ;
  rebasée sur `PQ_cdpI_P020_court_clean.cfg` (`Q1_c05_clean.msh`) : w19 = w20, 6/6 IDENTIQUE, history e45a7cbc0609c9f8 / frames 5ddf25b84b9c15f1 / vtu a1e108561139734a, pic 7 791,04 N — history et frames sont ceux de la chaîne w14 → w19 sur le maillage à orphelin (la garde `m ≤ 0` du contact, w13, excluait déjà ce nœud), seul le hachage VTU change (un nœud de moins). `check_orphans.py` sur les 30 + 4 maillages : orphelin
  dans TOUS les `T1_*.msh` non nettoyés (11) et `Q1_c05.msh` ; aucun dans les `T2_*`, `*_clean`, `P_cdpQ_v11_block`, `Q2_abq` → les decks
  de matrice_v2 (T1_*_clean, T2_*) passent la garde ; 16 decks perc3d pointent encore `Q1_c05.msh`.
- Coût de la garde C4 : deck dpr court `C_T1_R_P000_court_grid` (55 296 tets, 4 209 nœuds, OMP 2, matrice v2 en parallèle), w20 : `nanCheckEvery = 256` 103,7 / 104,7 s contre `nanCheckEvery = 0` 104,7 / 102,3 s (+0,7 %, dans le bruit) ; `nanCheckEvery = 1` (balayage à CHAQUE pas) 107,6 s (+4 %) ; history.csv identique dans les trois cas (`bitid_w20/cost_*.log`).
- Banc falsifiant `bitid_w20/selftest_gardes/` (`SELFTEST_gardes.md`, `run_gardes.py`, **30 / 30 PASS**) : T1_c05.msh original → erreur C3
  (fem3d et fdem3d), T1_c05_clean → passe ; tet plat / sliver / triangle plat à la main → erreurs nommées ; NaN → code 3 + ERROR.txt
  dans les 6 solveurs ; hydro en fem3d/fdem3d → refus ; clé fdem3d en fem3d, fem3d en fdem3d et en fdem → refus ; 9 témoins passent.
  Deux attentes du banc corrigées (pas le code) : l'orphelin est le nœud 9 (pas 18587) ; `dtFactor = 50` en traction 2D / DEM ne produit
  pas de NaN nodal (joints tous cassés, énergies NaN, champs finis, rc 0) → decks à dt × 1e200.
- Suite `verify_suite.py --tier fast` OMP 1 avec w20 : **TOUT PASSE (48/48)**, OMP 1, ~40 min sous charge (matrice v2 + bitid en parallele ; 22 min a vide), `bitid_w20/suite_fast_w20.log` + `.json` — aucun faux positif des gardes C3/C4/C6 sur les 48 tests.
- Points à relire : (1) code de retour NaN = **3** comme demandé (le plan C4 disait 2) ; (2) énergies internes non finies hors garde
  (fdem `joints : nan J/m` à rc 0) — ajouter les scalaires d'énergie au balayage est trivial mais leurs noms diffèrent par solveur ;
  (3) le registre est un instantané : nouvelle clé = relancer `gen_keys_by_mode.py` avant le build (à mettre dans la porte B7) ;
  (4) `meshOrphans = drop`, `meshMinDihedralDeg` et `fpTrap` du plan non faits (pas demandés au tour 1) ; (5) le point d'entrée
  `bitid_w12/PQ_cdpI_P020_court.cfg` du protocole bitid doit être remplacé par `_clean` dans `tools/bitid.py` (fichier de l'autre agent,
  non touché).
- 19:28 — **w20 (gardes des entrées, vague 0 du plan de robustesse) livré et relu** : nœud orphelin / masse nulle / tet dégénéré = erreur
  nommée dans les 6 solveurs ; `nanCheckEvery` (défaut 256, +0,7 %) = balayage réel de u, v, f → ERROR.txt + rc 3 ; `hydro*` refusé hors
  fdem ; registre des clés par mode (`tools/keys_by_mode.json`, 388 clés, 150 propres) → « clé X sans effet en mode Y » ; 30/30 bancs qui
  doivent échouer, 0 refus sur 899 decks du dépôt, suite fast 48/48 à OMP 1, bit-identité w19 = w20 sur 4 decks (49/49 fichiers).
  Conséquence voulue : `bitid_w12/PQ_cdpI_P020_court.cfg` (Q1_c05.msh, orphelin id 9) est REFUSÉ par w20 → ancre cdp rebasée sur
  `bitid_w20/PQ_cdpI_P020_court_clean.cfg` (mêmes history/frames que la chaîne w14→w19) ; 16 decks de perc3d pointent encore Q1_c05.msh.
  À harmoniser : code de retour NaN 3 (code) vs 2 (plan C4) ; les scalaires d'énergie NaN (joints fdem) ne sont pas balayés.
- 19:55 — **vague 0 du plan de robustesse terminée** : B5/B6 (`tools/bitid.py`, 8 decks `tests_f2/bitid/`, `tools/bitid_refs.json` ancré sur w18 à 4 threads, `tools/exe_manifest.json` : 47 exe dont 30 non documentés, `CHANGELOG.md` créé), C3/C4/C6 (w20). Passe complète w20 contre l'ancre : **8/8 IDENTIQUE**. Règle mesurée : les hachages fem3d dépendent du nombre de threads (pic 7 916,75 N à 4 threads contre 7 791,04 à 2) → comparer à threads égaux.
- 20:10 — **références sondées** D_T1_R_P000_probes / P100_probes finies : bit-identiques à C_T1_R_P000/P100 (43,43 kN, 0,615 mm) → les
  sondes n'altèrent rien, la comparaison brique/R aux cinq points est complète. **Démonstration direction figée** D_T1_fixed_P000 (w19,
  660 s) : F 42,8 kN (R 43,4), W 14,3 J, e_r 0,106, V_D09 5 314 mm³ (R 11 863) MAIS 20 157 éléments supprimés par le canal `erodeWfrac`
  (cratère 5 520 mm³) : l'énergie de traction SOMMÉE sur les trois directions atteint le seuil w ≥ frac·G_f/l_c beaucoup plus tôt qu'en
  scalaire → run non exploitable en l'état ; décision : critère d'érosion sur max_i w_i pour le modèle figé (w22, après w21), puis rejeu.
- Fil de coupe ARRÊTÉ à 19:58 : le GIF de C_T2_R_P000 (`matrice_v2/figures/gif_coupe_C_T2_R_P000.gif`, `make_gif_coupe.py`) montre un
  FENDAGE en coin du bloc 22 × 10 mm, pas une rainure (bande supprimée pointe → coin opposé, fragment supérieur poussé, F_x 55 N dans la
  fenêtre = poussée d'un fragment libre) ; 75 % des suppressions avant la première image. Le protocole T2 est à redessiner (bloc plus
  grand, appui du fond, passe, images toutes les 50 µs) avant toute phase C2/E.

## 2026-09-05 20:25 — HÉTÉROGÉNÉITÉ marche 2 lancée (go Fernando) ; décisions 1-3 du plan de robustesse
- `heterogeneite/` : 9 decks W_m{3,6,12}_s{1,2,3}_P000 = C_T1_R_P000 (R, bloc complet, T1_c05_clean) + `matWeibullM`, `fieldSeed`,
  `strengthCorrLength = 1e-3` (champ corrélé, un grain) + sondes des 5 points ; runner `run_hetero.sh` (w18, 14 threads, ~13 min/run).
- Décisions Fernando : (1) clé inconnue = ERREUR avec message expliquant la raison, `unknownKeys = warn` pour la transition (w21 en cours) ;
  (2) pas de suppression : nouvel arbre `FDEM/rockim_g0/` né sous git (branche g0, état f2 à la validation de w21), `rockim_f2` gelé à la
  fin des campagnes, worktrees p2-p4 archivés en bundle et laissés en place ; (3) « fusionne tout dans g0 sans rien toucher dans les
  dossiers respectifs » : les trois branches non fusionnées (joint-handoff, insertion-pointe, dif-intrinseque) sont portées dans g0 en
  clés opt-in avec bit-identité, f2 et les worktrees intacts.

## 2026-09-05 20:00 → 20:25 — w21 : la clé inconnue est une ERREUR (C0 / C1 / C7, agent codeur, tour 1)
- Décision Fernando 20:00 : clé inconnue = erreur avec la raison (obsolète → nouveau nom ; autre mode → le mode ; faute de
  frappe → suggestion ; sinon « inconnue de rockim »), `unknownKeys = warn` pour les vieux decks. Code : `Config` suit les
  clés consommées (stockage partagé entre copies, `shared_ptr`), `include/rockim/KeyGuard.hpp` (audit après `init()`,
  Levenshtein ≤ 2, politique, `config_effective.cfg`), registre étendu (`kKnown`, `kDynamicPrefix`, `kObsolete`,
  `tools/obsolete_keys.json`), `tools/scan_decks.py` + rapport, banc `bitid_w21/selftest_cles/` (17 decks, 17/17 + 20
  suite = 37/37), `rockim_f2w21.exe` (20:12), bitid 8/8 IDENTIQUE (4 fils, 20:35). C0 : `fragBrushV` → `fragBrushV0`
  dans 7 decks `bench_impact/configs/` + les 7 de `rockim_f2_wt` (les runs Kuru / pulvérisation / St Anne antérieurs
  tournaient SANS brossage). 8 clés `dfhPsi*` inconnues dans `impact3d_dpdfh*.cfg` : à trancher.
- Relecture adverse (20:25 → 20:45, `bitid_w21/relecture/RELECTURE_w21.txt`) : **W21_READY non écrit**. D1 : la règle
  « registre = légitime » laissait passer sans message toute clé lue par ≥ 2 solveurs dont aucun du mode (adverse a1 :
  `jointSoftening`, `insertion`, `bulkDamage`, `fragBrushV0` en fem3d → rc 0, « 4 du deck non lues »). D2 :
  `config_effective.cfg` n'était PAS rejouable en fem3d ni fdem (défauts écrits en clés actives → gardes « satellite
  orpheline »). Sa chaîne : bitid relecteur 7 ECHEC rc 0xC0000142 (lancement de processus, non probant ; celui du codeur
  sur le même exe = 8/8), run_cles 37/37, suite fast jamais exécutée (`--exe` relatif).

## 2026-09-05 20:45 → 21:00 — w22 : corrections D1 / D2 de la relecture (agent codeur, tour 2)
- Originaux w21 : `bitid_w22/orig/`. Registre : `tools/gen_keys_by_mode.py` exporte les **lecteurs** de chaque clé
  (`kReaders` dans `KeysByMode.hpp`, `readers` du json ; `shared` pour le code partagé et `ALWAYS_COMMON`).
- **D1** — `KeyGuard::audit` : clé non consommée du registre → légitime ssi le mode courant ou `shared` est lecteur (lecture
  conditionnelle : `capP0` si `dprCap`, `hydroStart` si `hydro`, `gravity` après init) ; sinon `cle 'X' sans effet en mode
  Y : cle du mode Z seulement` (un lecteur, sous-chaîne w20 conservée) ou `… : cles des modes Z1, Z2`. `tools/scan_decks.py`
  applique la même règle : **0 deck existant touché sur 957** (les seuls « sans effet » sont les decks adverses du
  relecteur, désormais exclus du balayage comme les `selftest_*`).
- **D2** — `config_effective.cfg` : lignes actives = clés du deck moins les fautives (consommées « # deck ligne n » ; non
  consommées légitimes « ; non lue a l'initialisation ») ; défauts du code **commentés** `# cle = valeur (defaut)`. Rejeu des
  six témoins `c11_ok_*` : rc 0 et `history.csv` bit-identique (partie C du banc ; fem3d 1799dac8f9b79eda, fdem
  0de2add74fea515d, fdem3d 5406f6365629bb4d, fem d81c00f7c1b7c601, dem 0d0f2ffce319d33b, dem3d 12efc6934ee61cc1).
- **Coût** (remarque du relecteur) — `Config::seal()` appelé par `main.cpp` après l'audit : les getters lisent la table
  sans verrou ni chaîne de défaut (`confineGaugeTime` à chaque pas). `Config::deck()` pour le fichier C7.
- `tools/scan_decks.py` : gardes pré-init de `main.cpp` modélisées (`tg_3d.cfg` thermal en fdem3d désormais vu), fichiers
  sans clé de solveur (matpoint `mp*`, cartes `law = cdp`, temporaires de selftest : 76) mis à part ; les gardes de
  maillage C3 restent hors de portée du juge statique (documenté). Rapport `tools/scan_decks_2026-09-05.md` régénéré :
  881 decks de solveur, 3 avec faute (tg_3d thermal ; `dfhPsi*` ×2 decks), 0 obsolète, 0 faute de frappe.
- `src/main.cpp` repassé en CRLF (diff réel 52 lignes contre l'original w20). Build complet `rockim_f2w22.exe` (20:51,
  `build_w22.log`, 0 avertissement, sha256 7e683603c53b2f90…). Banc `run_cles.py` A (21 decks : 17 w21 + c12 lecteurs
  multiples fem3d, c13 clés dem en fem3d, c14 `kpFactor` en fdem, c15 lecteurs multiples + `warn`) + C (6 rejeux) :
  **27/27 PASS** (20:54). Chaîne `bitid_w22/chaine_w22.sh` (run_cles A+C+B, bitid 4 fils, suite fast OMP 1) : résultats
  dans `bitid_w22/chaine.log`, `run_cles_ABC_w22.log`, `bitid_w22.log`, `suite_fast_w22.log`.
- 20:56 — `run_cles.py --suite 20` : **47/47 PASS** (A 21 + C 6 + B 20). 21:15 — **bitid w22 : 8/8 IDENTIQUE** (ancre w18, 4 fils,
  mêmes pics que w20/w21 : 7 916,75 / 6 775,14 / 30 390,2 / 19 541,3 / 205 788 N). Decks adverses du relecteur rejoués sur w22
  (`bitid_w22/adverse_w21_sur_w22.txt`) : a1 rc 1 (D1 fermé), a2-a10/r2 inchangés.

## 2026-09-06 — matériau PAR PHASE en éléments finis (`mode = fem3d`) — chantier « tessellation FEMDEM réutilisée »

**Demande de Fernando** : « utiliser le matériel qu'on a déjà côté femdem pour la tessellation Voronoï,
puis corriger et améliorer rockim pour qu'il puisse associer en fonction des phases un matériau plutôt
qu'un autre » — c'est-à-dire un essai **en éléments finis** (milieu continu, pas de joints discrets)
sur une microstructure de grains où chaque grain porte les propriétés de sa phase minérale.

**Pourquoi ce n'est pas un confort.** Le champ de Weibull des campagnes (`matWeibullM`) fait varier la
**résistance** sur un bloc de raideur parfaitement **uniforme** : il n'y a donc aucun contraste
élastique, donc aucune raison mécanique qu'une contrainte se concentre quelque part, donc le calcul ne
**pouvait pas** localiser. Un contraste de raideur — quartz 83,1 / feldspath 70 / biotite 29,3 GPa —
est ce qui concentre les contraintes aux frontières de grain. La capacité rend l'expérience possible.

### Ce qui a été écrit

- `include/rockim/Fem3dSolver.hpp` — `Elem` gagne `int phase, grain` ; `PhaseSet phases_` ;
  `std::vector<std::unique_ptr<MatLaw>> laws_` (une instance **par phase**, même clé `law`) avec
  `law_` conservé en **pointeur non propriétaire** sur `laws_[0]` (les sites qui n'appellent que
  `name()` / `sigmaCdp()` / `viscousOverstress()` sont inchangés) ; tables chaudes `rhoP_`, `rhoCdP_` ;
  `tetPhase_` / `tetGrain_` ; `phaseWeibull_`.
- `src/Fem3dSolver.cpp` — `phaseKeyGuards()` (clés de **joint** refusées, avec la raison),
  `auditPhaseKeys()` (famille `phase.<nom>.<propriété>`, deux messages distincts : nom de phase
  inconnu **vs** clé de loi écrite par phase), `buildMeshVoronoi()` (Tessellation3 **réutilisée**,
  compactage déterministe des sommets orphelins), lecture des `$PhysicalNames` + `groupPhase` dans
  `buildMeshFile()` (le lecteur **lisait les ntags et les jetait**), `checkMassAudit()`, masse
  condensée / CFL / pénalité / Lysmer / viscosité de volume **par phase**,
  `laws_[e.phase]->stress(...)`, champs `.vtu` `phase`/`grain`/`matE`/`matRho`/`matFt`, bilan par
  phase, qualité des tets, avertissement de verrouillage volumique, verdict `tension` supprimé en
  multiphase.
- `src/main.cpp` — deux barrières élargies à `fem3d` (`mesh = voronoi`, `phases`).
- `include/rockim/KeysByMode.hpp` — **régénéré** par `tools/gen_keys_by_mode.py` (le registre est une
  **sortie** du scan des sources, pas une entrée à éditer).
- `bench_phases/` — `make_bar2.py` (barre à deux couches nommées, avec ou sans `$PhysicalNames`),
  8 decks, `depouille.py` (9 contrôles), `erreurs.py` (13 cas fautifs).

### Le contrôle qui compte

Le mode de défaillance le plus grave est **silencieux** : câbler `Elem.phase`, les champs `.vtu`, les
bilans de fractions… et oublier `laws_[e.phase]->stress`. La figure est parfaite, la carte de grains
est colorée, et tous les chiffres sont ceux d'un bloc homogène. Aucun banc « trois phases identiques =
bit-identique » ne le détecte (le contraste changerait quand même la masse et la CFL).

**Le banc de Reuss le détecte.** Barre à deux couches en série (`dur` E₁ = 83,1 GPa en bas, `mou`
E₂ = 29,3 GPa en haut, interface à H/2, nœuds partagés), `law = elastic`, 384 tets, 0,6 s de calcul.
La solution est fermée : `1/E_app = f₁/E₁ + f₂/E₂`. Les mors raidissent les extrémités, on calibre
donc sur les **deux runs homogènes du même maillage**, ce qui élimine l'effet de mors au premier
ordre. **Mesuré 44,158 GPa contre 43,678 attendu : 1,10 %.** `E_app(série)/E_app(dur) = 0,5253` — il
vaudrait **1 exactement** si la loi n'était pas indexée par la phase.

### Résultats du banc (2026-09-06)

| contrôle | verdict |
|---|---|
| F2 neutralité, `mesh = file` : 2 phases identiques = deck sans `phases` | `history.csv` **octet à octet** |
| F2 variante qui **doit** échouer (contraste) | diffère ✔ |
| F3 borne de **Reuss** | 44,158 contre 43,678 GPa, **1,10 %** |
| F3 le contraste est vu par la **loi** | rapport 0,5253 (vaudrait 1 si non branché) |
| F3d commutativité par `groupPhase` | −1,48 % |
| F5 voronoï : 3 phases identiques = deck sans `phases` | **octet à octet** |
| F5 variante qui **doit** échouer (Red Bohus) | diffère ✔ |
| F4 CFL sur `c_P` **max** | dt ×0,8497 = √(60/83,1), exact |
| `erreurs.py` | **16/16** refusés avec le bon message |

Chemin voronoï en fem3d, cube 6 mm, `grainSize` 2 mm : 40 grains, 1514 tets, **391 nœuds partagés**,
aire extérieure = aire de la boîte à **6·10⁻¹⁶**, 609 faces intérieures entre grains distincts, audit
de masse à **3,7·10⁻¹⁶**, fractions réalisées 34,2 / 45,8 / 20,0 % pour 35 / 45 / 20 % visées, qualité
des tets : rapport diamètre inscrit / lc médian **0,734** (min 0,418 ; tétraèdre régulier 0,833).

### Preuves de non-régression (2026-09-06)

- **`tools/bitid.py` : 8/8 IDENTIQUE** contre l'ancre w18 (`rockim_f2w18.exe`, 655a87329c68b0d8…),
  `OMP_NUM_THREADS = 4`, exe `4fa8e13a647c5492…` issu d'une **reconstruction complète** (objets
  supprimés, `.hpp` modifié). Pics retrouvés **exactement** : 7 916,75 / 6 775,14 / 30 390,2 /
  19 541,3 / 205 788 N. Les trois decks `fem3d` exercent tous les chemins touchés — `cdp_PQ` (mesh =
  file, masse du chemin fichier, pénalité, quarterModel, confinement), `dpr_T1` (98 342 tets,
  **viscosité de volume** + hencky + signorini : c'est lui qui mesure la descente de `bvRhoC` dans
  `processElem`, et il porte une colonne `wBulk` dans `history.csv`), `sk2011_cyl` (chemin **grille**
  + **Lysmer cylindre**, la seconde boucle, celle qu'on oublie). `fdem3d_kuru9` (6 corps, **3 phases**,
  `groupBond`) prouve que le FEMDEM n'a pas bougé. Rapport : `bench_phases/VERDICT_bitid.txt`,
  `bench_phases/bitid_phases_fem3d.json`.
- **`tools/verify_suite.py` (tier fast, `OMP_NUM_THREADS = 1`) : 48/48 PASS avant et 48/48 après**,
  test par test et **valeur mesurée par valeur mesurée** (seules les durées diffèrent). La référence
  « avant » a été prise sur un build propre de `HEAD` (`build_base/`, `f2e90c78090d2812…`) : rien
  n'était en échec avant, rien ne l'est après — **aucun échec de plateforme préexistant à cacher**.
  `bench_phases/suite_avant_HEAD.txt` et `bench_phases/suite_apres_phases.txt`.
- **Piège de plateforme rencontré** : `verify_suite.py --exe <chemin RELATIF>` échoue
  (`FileNotFoundError` au premier test — c'est le « `--exe` relatif » déjà consigné le 05/09), et
  Apex One retient quelques secondes un exe fraîchement lié. Passer un chemin **absolu** et relancer.

### Limites, à écrire dans tout livrable qui s'appuie dessus

1. **Ce n'est pas un GBM cohésif.** Nœuds partagés = aucun joint = aucune frontière de grain
   intrinsèquement faible. La frontière **concentre** la contrainte, elle ne **s'ouvre pas** en tant
   que surface. Les clés de joint sont refusées, pas ignorées.
2. **Verrouillage volumique.** Tétraèdres linéaires, un point d'intégration, ni B-bar ni intégration
   sélective (vérifié). Ils sur-raidissent d'autant plus que `nu` est grand : si `nu` varie d'une phase
   à l'autre, l'artefact est **corrélé à la phase** et **sous-estime** le contraste — c'est-à-dire
   qu'il pousse vers la conclusion « le contraste fait moins que prévu ». Le solveur avertit.
3. **`lc = V0^(1/3)` sur des tets en cône.** La tessellation produit des éventails plats, surtout aux
   frontières de grain, là où l'étude regarde. Tant que ce point n'est pas mesuré contre le `lc_c`
   dépendant du taux (`OUTILS/compute_lc.py`), faire l'essai physique sur `mesh = file` (Gmsh) et ne
   garder `mesh = voronoi` que pour la géométrie.
4. **Contrainte de pointe à l'interface.** Le désaccord élastique à un coin de grain est singulier :
   site et instant d'amorçage restent dépendants du maillage. « La fissure suit les frontières de
   grain » n'est pas une preuve en soi.
5. **Clés d'option de loi globales.** Une seule loi pour toutes les phases ; `erodeD`, `dfh*`, `cdp*`,
   `meridian`… valent pour tout le bloc. Écrites par phase, elles sont refusées avec ce message.

### Trois trous fermés à la relecture de mon propre code (le même jour)

1. **`phase.rock.E` sans clé `phases`.** Sans `phases`, `PhaseSet::from` rend **une** fiche qu'il
   nomme « rock » et ne lit **aucune** clé `phase.*`. Une fiche `phase.rock.E = 70e9` passait donc
   l'audit (le nom « rock » existe !), était **consommée** par `keysWithPrefix`, jugée légitime — et
   restait parfaitement **inerte**. C'est exactement le motif interdit. Refusée en nommant la cause.
2. **`groupPhase.<groupe>` hors `mesh = file`.** Sans groupes physiques il n'y a personne à associer :
   la tessellation nomme ses grains par leur phase, la grille n'a aucun groupe. La clé était lue pour
   personne. Refusée.
3. **`phaseWeibull` posée sans objet** (pas de `matWeibullM`, ou une seule phase) : la permission
   n'autorisait rien. Refusée, comme les autres « satellites orphelins » du dépôt.

### Défaut PRÉEXISTANT relevé au passage (non corrigé ici, à trancher)

Sur le chemin `mesh = grid` de fem3d, `hmin_` est le pas **nominal** `min(dx, dy, dz)`
(`Fem3dSolver.cpp` l. 356) alors que le chemin fichier prend le **diamètre inscrit** 6V/A. Or
`Fdem3dSolver.cpp` l. 2457-2464 documente exactement ce piège et son prix : le découpage de Kuhn donne
un 6V/A médian ≈ 0,4 × arête (pire ≈ 0,2), le `hmin` nominal **surestime donc le dt stable jusqu'à
5 fois**, et c'est ce qui a produit 2,4 MJ d'énergie de bloc pour 16 J incidents le 2026-08-07 côté
FEMDEM. En fem3d seul `dtFactor = 0,3` masque l'écart, et à 0,2 d'aspect le rapport monte à 1,5 > 1.
**Ce chantier n'y touche pas** (aucune phase n'est admise sur le chemin grille), mais le défaut est
réel et indépendant.
