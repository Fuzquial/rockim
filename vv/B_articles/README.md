# Partie III : reproductions d'articles (B4 à B7) et faisabilité de B3

Ce dossier reprend les reproductions d'articles déjà faites avec rockim avant la campagne V&V
(`bench_abuaisha/`, `tunnel_edz/`, `configs_yan/`, `configs/stanne2025_*`, `configs/yang2026_*`),
les rend rejouables avec le binaire courant (`build_nofma/rockim`, macOS arm64) et fixe, avant tout
nouveau calcul, ce qu'un re-run doit retrouver. B1 (Guo) et B2 (Lisjak) sont traités ailleurs.

Ces bancs sont anciens : les écarts aux articles ont été constatés avant qu'aucun critère ne soit
écrit. Ici, le critère principal d'un re-run est la reproduction du résultat antérieur de rockim.
L'écart à l'article n'est contrôlé qu'a posteriori : il doit rester celui qui a été annoncé (même
signe, dérive bornée en points). Aucun de ces bancs n'est une validation au sens de la partie IV :
tous comparent rockim à la simulation calée d'un autre code, ou à un critère analytique.

## Organisation

| Élément | Contenu |
|---|---|
| `b_articles.py` | `JOBS(args)`, actions `prepare`, `smoke`, `run`, `analyse`, `all` ; références, critères, maillages |
| `decks/` | copies des decks d'origine (fins de ligne et encodage normalisés, contenu inchangé) |
| `meshes/` | maillages régénérés par les outils d'origine (ignorés par git) |
| `out/<cas>/` | sorties brutes ; `out/_fumee/` : fumée du 2026-10-04 et `smoke.json` |
| `resultats.json` | par cas : valeur publiée, rockim antérieur, rockim courant, critères et verdicts |
| `fig_b4_pression.pdf` | pression du fluide p(t) des runs AbuAisha antérieurs (historiques archivés) et courants |
| `fig_b_ecarts.pdf` | écart à la valeur publiée, rockim antérieur et courant, par grandeur |

Le script réécrit chaque deck dans `out/<cas>/deck.cfg` en ne changeant que `meshFile`, `outputDir`,
et `T` ou `frames` quand le tableau `CASES` le dit. Il n'ajoute aucune clé.

## Critères (fixés le 2026-10-04, avant tout re-run)

| Nature de la grandeur | Tolérance sur l'écart au rockim antérieur |
|---|---|
| avant rupture : onde élastique, pic de pression qui précède la propagation | 2 % |
| cinématique d'impact : enfoncement, pentes, instant de retournement | 5 % |
| post-pic ou fissuration : rayon d'EDZ, comptages, cratère, convergence | 10 % (15 % pour la masse de fragments et les facettes rompues de St Anne) |
| UCS de Yan | 0,15 MPa, tolérance du repère `ucs_yan_adaptive` de la suite |

Justification. Les maillages sont régénérés par Gmsh 4.15.2 sous macOS et non par la version
Windows d'origine ; leur réalisation diffère (voir chaque banc), et l'effet d'une graine de maillage
n'a jamais été mesuré sur ces bancs. Le chaos de plateforme est établi sur les grandeurs post-pic :
l'UCS de Yan donne 327 joints rompus et 1 288 insérés sous Linux, 313 et 1 308 sous MSVC.

Écart à l'article, a posteriori : même signe que l'écart antérieur et dérive au plus de 2 points
(B4, B6) ou 5 points (B5, B7). Seule exception a priori : l'ouverture de Parker (B4), jamais mesurée
après le correctif de signe du 2026-08-20, doit tenir à 10 % de la solution fermée (critère du
PLAN_ROBUSTESSE du 2026-09-05, G5).

## Fumée du 2026-10-04 (1 fil, T réduit, une trame)

Tous les decks démarrent avec le binaire courant. Le juge statique (`tools/scan_decks.py`) ne relève
aucune clé inconnue ni obsolète, et le binaire lit toutes les clés : les seules « non lues à
l'initialisation » sont des lectures conditionnelles légitimes (`absorbing`, `verifyFt`,
`budgetAbortPct`, `budgetAbortMin`).

| Cas | Éléments | dt (s) | Pas du cas complet | Coût estimé à 2 fils |
|---|---|---|---|---|
| B4_aniso, B4_iso | 189 403 triangles, 283 998 joints | 2,235e-8 | 147 600 et 156 600 | 2,2 h et 2,3 h |
| B4_e3_aniso, B4_e3_iso12 | idem | 2,235e-8 | 165 500 et 156 600 | 2,5 h et 2,3 h |
| B4_parker_c | 22 444 triangles | 8,34e-8 | 240 000 | 0,35 h |
| B4_parker | 192 666 triangles | 1,87e-8 | 1 068 500 | 11 h |
| B5_s5 | 106 222 triangles, 159 155 joints | 3,09e-6 | 80 800 | 0,6 h |
| B6_ucs | 2 410 triangles, 3 543 arêtes | 1,43e-8 | arrêt au pic | 33 s mesurés à 1 fil |
| B7_kuru_s25 | 10 597 tétraèdres | 2,48e-9 | 121 000 | 1 h |
| B7_stanne_s25 (jusqu'à 134 µs) | 24 148 tétraèdres | 4,10e-9 | 32 700 | 0,4 h |
| B7_stanne_137 | 109 160 tétraèdres | 1,65e-9 | 181 700 | 9 h (6 à 15 h) |

Le coût par pas est corrigé du temps de mise en place (runs de trois pas : 3,3 s pour le forage,
2,7 s pour rock137), divisé par 1,6 pour le second fil, et ne contient pas la croissance du contact
en phase de fracture. Sur le run St Anne rock137 de Windows, cette croissance a multiplié le coût
moyen par 2,9 ; elle est incluse dans l'estimation de B7_stanne_137 et de B7_kuru_s25. La machine
était chargée par une autre campagne pendant la fumée : les temps varient de 30 % d'un run à l'autre.

## B4 : AbuAisha et al. 2017, fracturation hydraulique en paroi de forage

Référence. AbuAisha, Eaton, Priest et Wong, J. Petrol. Sci. Eng. 154 (2017) 100-113, code Y-Geo
(FDEM 2D), section 3 et annexe A. Étude d'origine : `bench_abuaisha/` (README, VALIDATION_hydro.md)
et `ETUDES/rockim-bench-abuaisha-hydro.md`.

Cas. Forage R = 50 mm dans un bloc 8 × 8 m, maille de 3 mm en paroi (105 éléments sur le pourtour),
contraintes effectives in situ, pompe à débit imposé (`hydro = on`, `hydroInjection = rate`),
matériau de leur Table 1. Cas rejoués : `B4_aniso` (σ'_H = 6,8, σ'_h = 4,6 MPa) et `B4_iso`
(4,6 et 4,6 MPa), pompe dès t = 0 ; en complément lourd, `B4_e3_aniso` et `B4_e3_iso12`, protocole
de leur §3.2 (`hydroStart = 4e-4`). `B4_parker_c` et `B4_parker` : ouverture d'une fissure sous
pression uniforme, solution fermée de Parker (leur annexe A).

Grandeurs comparées. Pression de rupture (pic de `hydroP`), pression à la première insertion,
instant du pic ; ouverture au centre de la fissure de Parker.

Résultats antérieurs (binaires `rockim_e1.exe` à `rockim_e3.exe`, Windows/MSVC, 20 au 22/08/2026,
historiques archivés dans `FDEM/rockim/bench_abuaisha/historiques/` du dépôt de thèse) :

| Cas | Référence publiée | rockim antérieur | Écart |
|---|---|---|---|
| anisotrope, pompe dès t = 0 | 12,0 MPa (leur éq. 10, Hubbert-Willis) | 14,993 MPa à 2,928 ms | +24,9 % |
| isotrope 4,6 MPa | 14,2 MPa (2σ + ft) | 16,078 MPa à 3,177 ms | +13,2 % |
| anisotrope, protocole de l'article | 12,0 MPa | 14,999 MPa (insertion 13,324) | +25,0 % |
| isotrope 3,5 MPa, protocole de l'article | 12,0 MPa (état non publié) | 13,755 MPa (insertion 12,596) | +14,6 % |
| essai 2 (mouillage dès l'insertion), iso et aniso | 12,0 MPa | 12,806 et 14,145 MPa | +6,7 % et +17,9 % |
| Parker, w(0) | 0,0640 mm | aucune valeur après le correctif de signe | |

La valeur numérique de l'article est elle-même incertaine : environ 12,5 MPa dans le texte, 11,69 MPa
lus sur leur figure 11b. Le contrôle de signe, le croisement de volume (+0,5 à +0,7 %) et l'écart de
Lamé (−2,0 %) relèvent de la vérification et ne sont pas repris ici.

Indépendance du calage. Les paramètres sont ceux de leur Table 1, sans ajustement sur la pression de
rupture. La référence principale est un critère analytique (Kirsch et Hubbert-Willis), que leur
propre FDEM dépasse aussi ; il s'agit d'une reproduction du modèle de l'article, pas d'une validation.
Trois choix de rockim ne viennent pas de l'article : `fluidBulk = 2.2e9` (non publié, il fixe
l'horloge), l'essai 1 qui porte l'état isotrope à 3,5 MPa pour égaler la cible anisotrope, et l'essai
2 (`hydroWetDamage = 0`), choisi après avoir constaté le dépassement et sans fondement physique
(VALIDATION_hydro §4.3). Les cas rejoués ici excluent l'essai 2.

Statut de reproductibilité. Decks repris tels quels ; toutes les clés sont lues. Le maillage
`hf_bore.msh` est régénéré par `tunnel_edz/tools/make_circle_mesh.py` : 189 403 triangles et
283 998 joints, contre 284 124 joints actifs cités pour le maillage d'origine, donc une autre
réalisation. Le deck isotrope d'origine (`T = 3.0e-3`) s'arrête avant le pic (3,18 ms) ; l'historique
archivé va à 4 ms, donc le run de référence ne venait pas de cette version du deck : `T` est porté à
3,5e-3 (3,3e-3 pour l'anisotrope). Les runs ne sont pas lancés.

Coût. 2,2 h et 2,3 h à 2 fils pour les deux cas principaux, 0,35 h pour Parker grossier ; 2,3 à
2,5 h par cas du protocole de l'article et 11 h pour Parker au maillage de production.

Commande.

```bash
python3 vv/B_articles/b_articles.py run --cas B4_aniso B4_iso B4_parker_c --threads 2
python3 vv/B_articles/b_articles.py analyse
```

## B5 : Wang et al. 2024, zone endommagée d'un tunnel profond

Référence. Wang, Qiao, Zheng, He, Hu et Yan, Front. Earth Sci. 12:1517816 (2024), tunnel de Hutou
Beishan, code MultiFracS (FDEM 2D). Étude d'origine : `tunnel_edz/` et
`FDEM/rockim/tunnel_edz/RESULTATS_2026-08-17.md` du dépôt de thèse.

Cas. Tunnel en fer à cheval dans un carré de 100 m, maille fine 0,22 m, in situ hydrostatique
σ0 = 5 MPa, excavation par relâchement de la traction de paroi (`excavRelease`), insertion
adaptative, matériau de leur Table 1, T = 0,25 s. Cas rejoué : `B5_s5` (deck
`tunnel_ref_s5_lam1.cfg`, run d'origine `out_tun_ref_iso`).

Grandeurs comparées. Rayon maximal de l'EDZ, déplacement maximal, nombre de fissures par mode ;
convergence moyenne de paroi (non publiée, robuste aux blocs détachés).

Résultats antérieurs (`rockim_tun.exe`, Windows/MSVC, 2026-08-17) :

| Grandeur | Wang et al. | rockim antérieur | Écart |
|---|---|---|---|
| rayon maximal de l'EDZ | 19 m | 17,19 m (p95 14,49 m) | −9,5 % |
| déplacement maximal | 0,347 m | 0,389 m, porté par un bloc détaché ; paroi moyenne 0,069 m | +12 % |
| joints rompus | non publié | 14 935 (traction 4 817, mixte 2 371, cisaillement 7 747) | |
| hiérarchie des modes | cisaillement > mixte > traction | cisaillement > traction > mixte (seuil q = 0,5) | en partie |

Balayage σ0 = 3 à 7 MPa (non rejoué) : EDZ 12,8 / 15,4 / 17,2 / 17,4 / 18,6 m contre 11,1 / 16 /
19 / 22 / 22,8 m publiés, soit +15 % à 3 MPa puis −4 à −21 %, avec une saturation au-delà de
5 MPa que l'article ne montre pas.

Indépendance du calage. Paramètres de leur Table 1, non ajustés sur l'EDZ. La comparaison oppose
rockim à la sortie d'un autre code (MultiFracS) : reproduction code à code. Plusieurs choix
numériques diffèrent de l'article et pèsent sur le résultat : excavation par convergence-confinement
au lieu de la réduction du module du noyau, insertion adaptative au lieu d'intrinsèque, pénalité
20 E au lieu de 100 E, et `dampingLocal = 0.15`, dont le passage à 0,7 divise l'endommagement par 4
sur le maillage de fumée. Le résidu d'énergie du contact (416 kJ/m injectés contre 545 kJ/m de
fissuration) restait ouvert.

Statut de reproductibilité. Bloquant trouvé : le maillage régénéré par
`tools/make_unstructured_mesh.py tunnelhs` contient 4 nœuds de construction (centres des arcs)
qu'aucun triangle ne référence, et la garde C3 du binaire courant le refuse (« noeuds orphelins »),
ce que `rockim_tun.exe` ne faisait pas. `b_articles.py` retire ces nœuds de sa copie
(`drop_orphans_2d`, triangles inchangés). Le maillage obtenu a 159 155 joints contre 159 269 à
l'origine : autre réalisation. Fumée correcte ; run non lancé.

Coût. 0,6 h à 2 fils (44 min sous Windows à l'origine).

Commande.

```bash
python3 vv/B_articles/b_articles.py run --cas B5_s5 --threads 2
```

## B6 : Yan, Zheng et Wang 2023, UCS avec insertion adaptative

Référence. Yan, Zheng et Wang, IJRMMS 169 (2023), fig. 17 à 19. Études d'origine :
`ETUDES/rockim-insertion-adaptative.md` et `FDEM/rockim/BILAN_insertion_adaptative.md` ; repère de
suite `ucs_yan_adaptive` (`configs_yan/ucs_adap.cfg`).

Cas. Éprouvette 36 × 72 mm, maillage Delaunay interne d'environ 1,5 mm (graine 4211), plateaux
frottants (μ = 0,1) à 0,1 m/s, matériau de leur Table 1 (calé par eux sur Tatone), loi `yan`.

Grandeurs comparées. UCS, joints rompus et insérés ; module apparent à titre d'information.

| Grandeur | Yan et al. | rockim antérieur | rockim courant (2026-10-04) | Écart courant / antérieur |
|---|---|---|---|---|
| UCS | environ 51 MPa (figure) | 51,0735 MPa (suite, Linux) | 51,430 MPa | +0,70 % (critère 0,29 %) |
| joints rompus | | 327 (313 sous MSVC) | 438 | +34 % |
| joints insérés | | 1 288 (1 308 sous MSVC) | 1 511 | +17 % |
| module apparent | 15 GPa en entrée | 99,1 % du théorique (BILAN) | 15,86 GPa, 99,6 % de E/(1 − ν²) | |

Verdict. La reproduction du résultat antérieur échoue sur les trois grandeurs ; l'écart à l'article
(+0,8 % contre +0,1 %) reste dans la dérive admise, a posteriori. Le run est déterministe : le deck
d'origine et sa copie donnent des historiques identiques octet pour octet. Le repère
`ucs_yan_adaptive` de `tools/verify_suite.py` échouerait donc avec ce binaire. Deux causes sont
possibles et non départagées : une évolution du code depuis la prise du repère, ou un maillage
Voronoi différent (il est tiré par `std::uniform_real_distribution`, dont la suite de valeurs
n'est pas garantie identique d'une bibliothèque standard à l'autre ; le nombre d'éléments
d'origine n'a pas été consigné).

Indépendance du calage. Les paramètres de Yan sont calés par eux sur les essais de Tatone ; l'accord
à 51 MPa reproduit leur simulation calée et ne valide rien d'indépendant.

Coût. 33 s à 1 fil.

Commande.

```bash
python3 vv/B_articles/b_articles.py run --cas B6_ucs --threads 1
```

## B7 : Yang et al. 2025 et 2026, impact d'un insert unique

Références. Yang, Xiang, Naderi, Wang, Aising, Ugarte et Latham, IJRMMS 191 (2025) 106125,
calcaire de St Anne à 10,66 m/s ; et IJRMMS 206 (2026) 106660, granite de Kuru à 9 m/s ; code
Solidity (FDEM 3D). Études d'origine : `docs/COMPARAISON_yang2025_stanne_2026-09-14.md`,
`BANC_yang2026_impact.md`, `ETAT_yang2026_2026-09-11.md`, `mesures_ICL_2026-09-12.json`,
`ETUDES/rockim-impact-yang-insert-unique.md`.

Cas. Train de frappe complet (piston, bit, insert brasé), roche cylindrique de 125 × 150 mm. Trois
cas : `B7_stanne_137` (deck `stanne2025_rock137_visc0.cfg`, 300 µs, le run de comparaison) ;
`B7_stanne_s25` (même deck sur la roche grossière rock25, arrêté à 134 µs, contrôle de
reproductibilité bon marché) ; `B7_kuru_s25` (banc court Kuru à s = 2,5, loi `plastic`, 300 µs).

Grandeurs comparées. Pic de la première onde à mi-bit, vitesse d'indentation (pente entre 10 et
90 % de l'enfoncement), enfoncement maximal, fin de charge ; rayon de cratère et masse de fragments
(St Anne) ; force de la roche sur l'insert et son impulsion (rock25).

Résultats antérieurs, St Anne rock137 (`rockim_g1y19.exe`, Windows/MSVC 14 fils, 105 690 s,
2026-09-14), valeurs publiées lues sur les figures 9 et 10 (±10 %) :

| Grandeur | Yang FDEM | Yang essai | rockim antérieur | Écart à Yang FDEM |
|---|---|---|---|---|
| contrainte au bit, première onde | 200 MPa | 200 MPa | 212,3 MPa | +6 % |
| vitesse d'indentation | 8,0 m/s | 7,9 m/s | 7,34 m/s | −8 % |
| enfoncement maximal | 1,10 mm | 1,05 mm | 1,292 mm | +17 % |
| fin de charge | 254 µs (texte) | | 284,3 µs | +12 % |
| rayon de cratère | 13,5 mm | 13 mm | 12,0 mm | −11 % |
| longueur radiale | 32 mm | 31 mm | 29,6 mm, ou 17,7 mm selon la définition | −8 % ou −45 % |
| masse de fragments | 2,5 g | 1,2 g | 4,70 g | +88 % (+290 % à l'essai) |
| vitesse de rebond | 5,6 m/s | 4,2 m/s | non mesurable à 300 µs | |

Kuru (Yang 2026). Aucun run complet avec la loi finale. Le seul jeu comparé est le banc court
s = 2,5 à loi `plastic` (`rockim_g1y10` et `g1y11`, Windows 10 fils, 4 090 s, 2026-09-13) : retour
du bit à 270 µs (Yang 255), vitesse d'indentation maximale instantanée 6,07 m/s (5,62),
enfoncement 0,87 mm (environ 1,0), contrainte au bit 174 MPa (environ 160), 32 joints rompus et un
rebond de 0,22 contre 0,83 : « sain et faux » selon l'étude d'origine. Le run s = 1 de la nuit du
13/09 a été arrêté vers 113 µs avec un piston 36 % plus énergétique que celui de l'article.

Indépendance du calage. Les paramètres matériau sont ceux de leurs Tables 1 et 4, calés par les
auteurs. La loi de joint retenue (`plastic` + plage `coulomb` + `jointSecantRatchet`) et la pénalité
(26,316 = p0/(2E)) ont été choisies pendant la même campagne, sur les bancs Kuru et sur des critères
de conservation de l'énergie, en regardant les critères de Yang. La comparaison St Anne n'est donc
pas indépendante du chemin de mise au point, et les valeurs publiées sont des lectures de figures.

Statut de reproductibilité. Toutes les clés sont lues. Maillages régénérés : rock25 24 148
tétraèdres contre 24 010, rock137 109 160 contre 108 667 (autres réalisations). La commande de
`REPRODUIRE_stanne_radiales_2026-09-14.md` (rapport de raffinement 1,37) donne 54 185 tétraèdres :
le maillage rock137 de la série correspond à un rapport 1,0 (`MAILLAGE_serie_2026-09-13.md`), que
`b_articles.py` utilise. Le maillage Kuru s = 2,5 est copié de `meshes/` (10 597 tétraèdres, identique
à sa régénération). Le binaire courant avertit que la raideur tangentielle du contact dépasse de
1,43 fois celle que budgète le pas de temps (decks St Anne) ; à surveiller. Runs non lancés.

Coût. B7_stanne_s25 0,4 h et B7_kuru_s25 1 h à 2 fils ; B7_stanne_137 environ 9 h à 2 fils
(entre 6 et 15 h), réservé aux runs lourds.

Commande.

```bash
python3 vv/B_articles/b_articles.py run --cas B7_stanne_s25 B7_kuru_s25 --threads 2
VV_B_LOURDS=1 python3 vv/run_queue.py B_articles --slots 1 --threads 2   # avec St Anne rock137
```

## B3 : Fukuda et al. 2020, faisabilité

Référence. Fukuda, Mohammadnejad, Liu, Zhang, Zhao, Dehkhoda, Chan, Kodama et Fujii, RMRE 53 (2020)
1079-1112, code 3D Y-HFDEM IDE sur GPU.

Ce qui est publié. Les vérifications de contact (sphère sur un plan, énergie normalisée en fonction
de l'amortissement ; bloc glissant de Xiang avec μ = 0,5) ; un UCS et un brésilien quasi statiques
sur un calcaire (D = 51,7 mm, épaisseur 25,95 mm, maille de 1,5 mm, 187 852 tétraèdres pour le
disque, plateaux rigides courbes, 0,01 m/s, échelle de masse 5, Table 1 complète) ; un brésilien
dynamique au SHPB sur marbre. Le pic brésilien n'est publié que sous forme de courbe (fig. 13d) et la
comparaison à Hondros que sous forme de profils (fig. 13e). Les paramètres sont calés par essais et
erreurs sur l'UCS de laboratoire, et le SHPB sur l'essai lui-même.

Verdict. Le bloc glissant est déjà couvert par P2.1 et la sphère sur un plan par P2.2. Le brésilien
3D est faisable avec les clés existantes (fdem3d intrinsèque, plateaux maillés comme en P4.2) mais
non trivial : il faut mailler des plateaux courbes, rockim n'a pas d'échelle de masse, et le coût
est de plusieurs jours (pas de 4,5 ns chez eux avec l'échelle de masse, de l'ordre de 2 ns sans elle,
soit 5 à 10 millions de pas sur 188 000 tétraèdres). La
seule grandeur chiffrable serait un pic lu sur figure, issu d'un jeu de paramètres calé, et le profil
de Hondros relève de la vérification élastique déjà prévue en P4.2. B3 n'est pas implémenté ; son
apport propre se réduit à la sensibilité du mode de rupture au maillage structuré (leur fig. 14),
qualitative.

## Points ouverts

1. Lancer les cas par défaut (environ 7 h de calcul à 2 fils au total) quand la machine est libre,
   puis `analyse`.
2. B6 : trancher entre évolution du code et maillage Voronoi dépendant de la plateforme, avant de
   rebaser le repère `ucs_yan_adaptive` sous macOS.
3. B4 : mesurer l'effet d'une graine de maillage (proposé à l'origine, jamais fait), sans quoi la
   tolérance de 2 % sur le pic reste une hypothèse.
4. B7 : la vitesse de rebond, critère le plus sensible au contact, demande T = 450 à 500 µs sur
   rock137 (deck `stanne2025_rock137_visc0_T500.cfg`, non repris ici).
