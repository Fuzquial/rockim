# P4 : essais de flexion trois points et brésilien de Guo 2014 en fdem3d

Bancs P4.1 et P4.2 de la campagne V&V (docs/VV_campagne.md), avec leur volet B1 : reproduction des
essais du chapitre 3 de la thèse de Guo (2014, Solidity 3D) et confrontation à des références
physiques indépendantes du code. Script : `p4_essais_guo.py`. Valeurs de Guo extraites :
`guo_reference.json`. Préparé le 2026-10-04 ; seul un essai à blanc a tourné (voir la fin).

## 1. Objectif

P4.1 vérifie que rockim donne la raideur élastique d'une poutre en flexion trois points et une
charge de rupture compatible avec la mécanique de la rupture cohésive (effet d'échelle du module de
rupture). P4.2 vérifie le champ de contrainte élastique d'un disque comprimé diamétralement
(solution de Hondros) et la résistance brésilienne f_bt = 2F/(πDt) au regard de ft. Les deux
bancs mesurent aussi la sensibilité à la vitesse de chargement et au maillage, puis comparent les
pics à ceux de Solidity (B1).

## 2. Montages

### P4.1, flexion trois points (Guo §3.3, tableau 3.1, fig. 3.1 à 3.7)

Poutre de 320 × 40 × 20 mm (x, y, z), ρ = 2 340 kg/m³, E = 26 GPa, ν = 0,2, ft = 3 MPa,
c = 10,5 MPa, φ = 30°, Gf = 30 J/m². Trois platines d'acier (ρ = 7 850 kg/m³, E = 200 GPa,
ν = 0,28) de 5 × 5 × 20 mm : deux appuis sous les extrémités (axes à x = 2,5 et 317,5 mm,
portée L = 315 mm ; Guo ne cote pas la portée, environ 314 mm sur sa figure 3.1a) et un poinçon au
milieu. Comme chez Guo, les trois platines bougent : poinçon vers le bas, appuis vers le haut,
chacun à Vy, avec une rampe linéaire de 1 ms. Frottement 0,6 entre lèvres de fissure, 0,1 entre
acier et roche.

Maillage structuré à 24 tétraèdres par cube (sommets, centres des faces, centre du cube), le
découpage de Guo : 49 152 tets à 5 mm et 393 216 à 2,5 mm, exactement ses nombres. Les platines
sont maillées à 5 mm (96 tets chacune ; Guo : 87 pour les trois).

### P4.2, essai brésilien (Guo §3.4, tableau 3.2, fig. 3.8 à 3.12)

Disque D = 40 mm, t = 15 mm, axe z, mêmes ρ, E, ν, ft que P4.1, c = 15 MPa, Gf = 50 J/m². Deux
platines planes d'acier de 24 × 5 × 15 mm en contact tangent avec le disque, mobiles en sens
opposés à Vy, rampe linéaire de 0,2 ms. Guo prend des platines rigides ; ici elles sont en acier
continu (E = 200 GPa, huit fois la roche), ce qui élargit le contact de Hertz de 6 % (E* de 24,1 GPa au lieu de
27,1 GPa) et ne change pas la charge de rupture au premier ordre.

Maillage non structuré Gmsh : 49 856 tets à 1,2 mm pour le disque (Guo : 51 690), 11 101 à 2 mm ;
platines à 2 mm (3 277 tets ; Guo : 2 854). Le cercle est construit en quatre arcs pour qu'un
sommet du maillage tombe exactement sur les génératrices de contact : pas de jeu initial.

## 3. Références physiques

### P4.1

Raideur élastique entre la flèche relative de l'axe neutre au milieu et au droit des appuis, poutre
d'Euler-Bernoulli corrigée du cisaillement de Timoshenko (κ = 5/6, G = E/(2(1 + ν))) :

δ/F = L³/(48 E I) + L/(4 κ G A), I = b h³/12, A = b h.

Avec L = 315 mm : 48EI/L³ = 4,259 MN/m, correction de cisaillement 4,6 %, K = 4,070 MN/m. Le
script évalue la formule aux positions réelles des sommets retenus par `point.` (porte-à-faux en
rotation rigide si le sommet est hors portée).

Module de rupture : σ_r = 3 F L/(2 b h²). Pour une fissure cohésive, σ_r/ft dépend de h/l_ch avec
l_ch = E Gf/ft² = 86,7 mm, soit h/l_ch = 0,46. Le rapport tend vers 1 pour une poutre très
grande devant l_ch et croît quand h/l_ch diminue (Petersson 1981, Hillerborg 1983) ; à
h/l_ch ≈ 0,5 il vaut environ 1,6 à 1,8 pour un adoucissement linéaire. Bande retenue : 1,3 à 2,1.
La charge qui porte la contrainte de flexion élastique à ft vaut 203 N.

Première fréquence de flexion (appuis simples) : f₁ = (π/2L²)√(EI/(ρA)) = 609 Hz, période
T₁ = 1,64 ms. Le rapport t_pic/T₁ qualifie le caractère quasi statique d'un run.

### P4.2

Contraintes de Hondros (1959) sur le diamètre chargé, charge F répartie sur un arc 2α, ρ = r/R,
p = F/(2αRt), traction positive :

σ_xx = (2p/π) [ (1 − ρ²) sin 2α / (1 − 2ρ² cos 2α + ρ⁴) − arctan( (1 + ρ²)/(1 − ρ²) tan α ) ]

σ_yy = −(2p/π) [ (1 − ρ²) sin 2α / (1 − 2ρ² cos 2α + ρ⁴) + arctan( (1 + ρ²)/(1 − ρ²) tan α ) ]

Au centre, σ_xx tend vers 2F/(πDt) quand α tend vers 0. Le demi-angle α = b/R vient de la
largeur de contact de Hertz cylindre-plan, b = √(4 (F/t) R/(π E*)), 1/E* = (1 − ν²)/E +
(1 − ν_a²)/E_a ; à F = 2 827 N, b = 0,45 mm, α = 0,022 rad. Sur le diamètre chargé σ_xx est la
contrainte principale majeure, comparée à `sigma1` des frames VTU (moyenne des éléments du disque
dont le centroïde est à moins d'une maille de l'axe, par tranches de R/8).

Résistance brésilienne : f_bt = 2F/(πDt). Pour un disque sans défaut et un critère de rupture en
traction, la contrainte du centre atteint ft quand F = πDt ft/2 = 2 827 N. Le rapprochement
élastique des platines à cette charge, estimé par Johnson (1985), vaut 0,047 mm ; il sert à choisir
la durée T, jamais de critère. Temps de traversée du disque : D/c = 12 µs.

## 4. Mise en œuvre

fdem3d n'a ni scénario `brazilian` (il n'existe qu'en 2D, `FdemSolver`) ni flexion. Les deux
essais sont donc montés en `scenario = loads` (§5.21 de la documentation) avec des corps nommés :

- éprouvette `poutre` ou `disque`, phase `roche`, insertion adaptative (variante principale) ou
  intrinsèque avec `jointPenaltyFactor = 20` (variante proche de Solidity) ;
- platines `appui_g`, `appui_d`, `poincon` ou `plateau_h`, `plateau_b`, phase `acier`
  (`phase.acier.*`, `groupPhase.<corps>`), `groupContinuum.<corps> = true` ;
- `velocity.<platine> = 0 ∓Vy 0` et `amplitude.<platine> = 0 0 t_r 1` (rampe linéaire de Guo) ;
- `contact = potential`, `contactMu = 0.6`, `contactMu.roche = 0.6`, `contactMu.acier = 0.1`
  (règle de paire : minimum, donc 0,1 entre acier et roche) ;
- `contactForcePairs` : force de contact de chaque platine sur l'éprouvette, colonnes
  `Fc_<platine>_<éprouvette>_y`. C'est la mesure de F : la réaction `RF_<platine>_y` de la
  liaison contient en plus l'inertie de la platine pendant la rampe (mesuré à l'essai à blanc) et
  ne sert que de contrôle ;
- P4.1 : points `mil`, `sg`, `sd` sur l'axe neutre (y = h/2, z = b/2) au milieu et au droit des
  appuis, et `fib` en fibre inférieure au milieu, avec `force.<g> = 0 0 0` pour écrire leur
  déplacement ; flèche δ = −(U_mil − (U_sg + U_sd)/2) ;
- P4.2 : F = (F_haut + F_bas)/2 (éq. 3.1 de Guo), ε_yy = Δd/D avec Δd le rapprochement des
  platines (éq. 3.2) ;
- `writeRuptureFields = true` ; amorçage lu dans la dernière frame `fdem3d_joints_*.vtu` :
  barycentre des cinq premières facettes rompues (plus petits `tBreak`).

Aucune clé nouvelle, aucune modification de `src/`. Pas d'amortissement (`dampingLocal = 0`, défaut
de `loads`) : Guo applique un modèle viscoélastique aux platines seulement.

### Variantes et durées

T = t_r/2 + d/(2 Vy), avec d = deux fois la flèche élastique au pic attendu (σ_r = 1,7 ft) pour
P4.1 et 2,5 fois le rapprochement de Johnson à 2 827 N pour P4.2, majoré de 50 % en intrinsèque.
Le pic doit être suivi d'une chute de 30 % pour être déclaré atteint.

| Run | maille | Vy (m/s) | insertion | T (ms) | dt (ns) | tets |
|---|---|---|---|---|---|---|
| p41_ad_h5_v0p05 | 5 mm | 0,05 | adaptative | 2,20 | 15,3 | 49 440 |
| p41_ad_h5_v0p02 | 5 mm | 0,02 | adaptative | 4,74 | 15,3 | 49 440 |
| p41_ad_h5_v0p01 | 5 mm | 0,01 | adaptative | 8,99 | 15,3 | 49 440 |
| p41_ad_h2p5_v0p02 | 2,5 mm | 0,02 | adaptative | 4,74 | 7,63 | 393 504 |
| p41_in20_h5_v0p02 | 5 mm | 0,02 | intrinsèque, pf 20 | 6,86 | 9,19 | 49 440 |
| p42_ad_h1p2_v0p1 | 1,2 mm | 0,1 | adaptative | 0,68 | 3,31 | 53 133 |
| p42_ad_h1p2_v0p05 | 1,2 mm | 0,05 | adaptative | 1,27 | 3,31 | 53 133 |
| p42_ad_h1p2_v0p02 | 1,2 mm | 0,02 | adaptative | 3,02 | 3,31 | 53 133 |
| p42_ad_h2_v0p02 | 2 mm | 0,02 | adaptative | 3,02 | 5,69 | 14 359 |
| p42_in20_h2_v0p02 | 2 mm | 0,02 | intrinsèque, pf 20 | 4,48 | 3,39 | 14 359 |

Le pas de temps a été mesuré par un lancement de quelques pas sur chaque maillage complet.

Choix de la vitesse. Guo trouve en flexion un pic à 0,02 m/s supérieur de 1,5 % seulement au pic
convergé (0,01, 0,005 et 0,002 m/s confondus, fig. 3.3), et au brésilien des pics identiques à
0,02 et 0,01 m/s (fig. 3.9). La vitesse principale est donc 0,02 m/s, la plus rapide de son
domaine convergé : elle divise le coût par 4 par rapport à 0,005 m/s, la vitesse de sa figure 3.7.
Une réserve : la poutre de rockim est environ deux fois plus raide que celle de Guo (sa raideur
sécante est 2,0 MN/m contre 4,07 MN/m en théorie), elle atteint donc son pic deux fois plus tôt,
vers 1,6 T₁ à 0,02 m/s contre environ 4,6 T₁ chez Guo. La série 0,05 / 0,02 / 0,01 m/s à 5 mm
mesure la convergence propre à rockim au lieu de la supposer. Au brésilien, t_pic vaut environ
100 fois le temps de traversée à 0,02 m/s : le régime est quasi statique, et la série
0,1 / 0,05 / 0,02 m/s le vérifie.

### Coût estimé

Débit mesuré à l'essai à blanc, un fil : 1,45 × 10⁶ tet-pas/s en flexion, 0,9 × 10⁶ au brésilien
(contact par potentiel plus chargé). À deux fils, environ 1,7 fois plus. Estimations à deux fils :

| Run | pas | durée à 2 fils |
|---|---|---|
| p41_ad_h2p5_v0p02 | 6,2 × 10⁵ | 28 h |
| p41_in20_h5_v0p02 | 7,5 × 10⁵ | 5 h |
| p41_ad_h5_v0p01 | 5,9 × 10⁵ | 3,3 h |
| p41_ad_h5_v0p02 | 3,1 × 10⁵ | 1,7 h |
| p41_ad_h5_v0p05 | 1,4 × 10⁵ | 0,8 h |
| p42_ad_h1p2_v0p02 | 9,1 × 10⁵ | 9 h |
| p42_ad_h1p2_v0p05 | 3,8 × 10⁵ | 3,8 h |
| p42_ad_h1p2_v0p1 | 2,1 × 10⁵ | 2 h |
| p42_in20_h2_v0p02 | 1,3 × 10⁶ | 3,5 h |
| p42_ad_h2_v0p02 | 5,3 × 10⁵ | 1,4 h |

Total environ 58 h de créneau à deux fils ; sur quatre créneaux, la durée est fixée par le run
de 2,5 mm (28 h). Guo à 0,005 m/s en 2,5 mm demanderait environ 31 ms simulées, soit 4 × 10⁶ pas
et une douzaine de jours : hors de portée, d'où le choix de 0,02 m/s. Les frames VTU du run de
2,5 mm pèsent environ 190 Mo chacune : `frames = 4` en flexion, 10 au brésilien.

## 5. Métriques

| Banc | Métrique | Définition |
|---|---|---|
| P4.1 | K_sim | pente de F en fonction de δ entre 10 et 40 % du pic, avant le pic |
| P4.1 | K_err | K_sim/K_th − 1 |
| P4.1 | F_pic, δ_pic, t_pic/T₁ | maximum de F, flèche et temps réduit associés |
| P4.1 | σ_r/ft | 3 F_pic L/(2 b h²)/ft |
| P4.1 | amorçage | premières facettes rompues à moins de 2h du milieu et à moins de h de la fibre inférieure |
| P4.2 | err_centre | σ1 moyen à |y| ≤ 0,15 R contre σ_xx de Hondros au centre, à la frame où F ≈ F_pic/2 |
| P4.2 | L2 | écart L2 relatif du profil σ1(y) à Hondros sur |y| ≤ 0,8 R |
| P4.2 | f_bt/ft | 2F_pic/(πDt)/ft |
| P4.2 | amorçage | premières facettes rompues à |y − y_c| ≤ R/2 |
| tous | contrôles | F des appuis contre F du poinçon (équilibre), RF contre Fc, résidu du bilan B4 |

## 6. Critères d'acceptation (fixés le 2026-10-04, avant tout calcul complet)

| Critère | Seuil |
|---|---|
| P4.1 raideur, 2,5 mm, adaptative | ≤ 3 % |
| P4.1 raideur, 5 mm, adaptative | ≤ 6 % (rigidité de flexion des tets linéaires) |
| P4.1 σ_r/ft, runs adaptatifs | entre 1,3 et 2,1 |
| P4.1 vitesse, pics à 0,01 et 0,02 m/s (5 mm) | ≤ 3 % |
| P4.1 maillage, pics 5 mm et 2,5 mm à 0,02 m/s | ≤ 10 % (Guo : −7 %) |
| P4.1 amorçage en fibre inférieure au milieu | oui |
| P4.1 B1, pic contre Guo (298,8 N, 2,5 mm) | ≤ 15 % |
| P4.2 Hondros au centre | ≤ 5 % |
| P4.2 Hondros, profil sur 0,8 R | ≤ 10 % |
| P4.2 f_bt/ft, runs adaptatifs | entre 0,85 et 1,15 |
| P4.2 vitesse, pics à 0,02 et 0,05 m/s (1,2 mm) | ≤ 3 % |
| P4.2 maillage, pics 2 mm et 1,2 mm à 0,02 m/s | ≤ 10 % |
| P4.2 amorçage central | oui |
| P4.2 B1, f_bt contre Guo (1,96 MPa) | ≤ 15 % |
| bilan B4 | résidu ≤ 1 % de l'échelle |

Prédictions écrites avant le calcul. En flexion, le pic de rockim devrait être proche de celui de
Guo (σ_r/ft = 1,47 chez lui, dans la bande). Au brésilien, un échec du critère B1 est attendu :
Guo obtient f_bt/ft = 0,65 et l'attribue à un amorçage par cisaillement sous les platines ; un
modèle qui respecte Hondros et amorce au centre doit donner f_bt proche de ft, à 50 % au-dessus de
Guo. Les critères f_bt/ft et B1 ne peuvent pas passer ensemble ; le premier est la référence
physique, le second documente l'écart à Solidity. La raideur de Guo (2,0 MN/m en flexion, environ
3 MN/m au brésilien, avec ε_yy ≈ 1,8 % au pic, 15 à 20 fois le rapprochement élastique) n'est pas
une cible : elle traduit la souplesse de ses joints intrinsèques et de ses pénalités.

## 7. Lancer

Depuis la racine du dépôt :

```
python3 vv/P4_essais_guo/p4_essais_guo.py ref                  # références et durées T
python3 vv/P4_essais_guo/p4_essais_guo.py prepare              # maillages et decks
python3 vv/run_queue.py P4_essais_guo --slots 4 --threads 2    # campagne
python3 vv/P4_essais_guo/p4_essais_guo.py analyse              # resultats.json, figures
python3 vv/P4_essais_guo/p4_essais_guo.py all --smoke --threads 1   # essai à blanc
python3 vv/P4_essais_guo/p4_essais_guo.py run --only p42_ad_h2     # un sous-ensemble
```

Sorties : `out/<run>/` (deck, journal, `history.csv`, frames), `meshes/`, `resultats.json`,
`fig_flexion.pdf/png` (F en fonction de δ, droite de Timoshenko, pic de Guo, charge à σ_r = ft),
`fig_bresilien.pdf/png` (F en fonction de ε_yy ; profil σ1 contre Hondros).

## 8. Essai à blanc (2026-10-04, un fil)

Maillages grossiers (poutre à 10 mm, 6 432 tets ; disque à 4 mm, 5 097 tets), Vy = 1 m/s, rampe
de 10 et 5 µs, T = 60 et 40 µs : 11 s et 26 s de calcul, résidu B4 de 0,0076 % et 0,0002 %. La
chaîne fonctionne de bout en bout : corps et phases reconnus, frottement par phase (0,6 et 0,1),
`contactForcePairs`, points de mesure, rupture en flexion (14 facettes), lecture des frames et
profil de Hondros. Les valeurs ne signifient rien : à 1 m/s et 40 µs la réponse est dynamique
(pic du disque à 2,1 ft, flexion dominée par l'inertie). Elles sont dans `resultats_smoke.json`.
L'essai a montré que la réaction de liaison contient l'inertie de la platine (2 900 N à 3 µs au
brésilien, contre zéro de contact) : la mesure de F passe par `contactForcePairs`.

## 9. Doutes et limites

- Portée de Guo non cotée (314 ± 4 mm) ; ici 315 mm.
- Guo ne dit pas si δ_y est absolu ou relatif aux appuis ; avec trois platines mobiles la
  différence est d'un facteur voisin de 2 sur la raideur.
- Le texte de Guo donne 1 844,76 N au brésilien à 0,01 m/s, sa figure 3.9 plafonne à 1 809 N.
- Gf_II : `gfShearFactor` reste à son défaut (10) ; Guo ne donne pas l'énergie en mode II au
  chapitre 3.
- Platines en acier continu au brésilien, rigides chez Guo.
- Pas d'arrêt automatique après le pic en `scenario = loads` (`stopPeakDrop` ne vaut qu'en
  `tension`) : T est fixé à l'avance et le script signale un pic non atteint.
