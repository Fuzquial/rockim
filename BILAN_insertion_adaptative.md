# L'insertion adaptative face aux macro-fissures — bilan de l'étude

Campagne du 2026-08-22 au 2026-08-25, close à la demande de Fernando (le
volet DP-DFH est GELÉ en l'état, section 8). Point de départ : sur le tunnel
profond de Wang et al. (Front. Earth Sci. 12:1517816, 2024, code MultiFracS)
comme sur l'impact de Yang et al. (IJRMMS 191, 2025), rockim reproduit les
grandeurs globales mais produit un **nuage diffus de fissures courtes** là où
les articles montrent des **macro-fissures localisées découpant des blocs**.
Hypothèse de Fernando : « le schéma adaptatif est la raison principale ; le
critère d'insertion et la manière d'insérer sont incorrects ». L'étude a
consisté à instrumenter cette hypothèse, la confronter à la littérature, et
la tester par contrôles croisés.

Tout est committé (rockim_p1 `main`, rockim_p2 branche `insertion-pointe`) ;
figures et métriques sur le drive (`FDEM\rockim\tunnel_edz`, `...\bench_impact`) ;
volets bibliographiques dans `biblio_insertion/`.

---

## 0. Deux résultats préalables qui cadrent tout le reste

**(a) Le déficit de déformation était un artefact de durée, pas de schéma.**
À T = 0,25 s la paroi convergeait encore linéairement (~6 mm/21 ms, le
dernier cinquième du temps portait 20 % du déplacement). Wang et al.
intègrent jusqu'à stabilisation (300 000 pas). Prolongé à T = 1,0 s, notre
U de paroi passe de 0,148 à 0,574 m hors débris (leur publié : 0,347 —
notre run n'était TOUJOURS pas stabilisé à 1 s, gain de 15 % sur le dernier
cinquième). Toute comparaison de motif doit donc se faire à **temps depuis
le relâchement égal**, ce que tous les outils de l'étude font.

**(b) Le balayage λ valide la physique globale.** Nombre de fissures
décroissant avec λ (24 297 → 9 069 pour λ = 0,5 → 1,5 ; leur fig. 16 :
~14 000 → ~5 000, même rapport ~2,7) ; l'ellipse d'EDZ bascule d'axe
horizontal à vertical en traversant λ = 1, comme chez eux. Le désaccord
n'est PAS global : il est dans la **texture** du réseau de fissures.

## 1. Le critère n'insère pas « trop tôt » — il insère trop uniformément

Audit du code ([src/FdemSolver.cpp:1925](FDEM/rockim/rockim_p1/../rockim_p2/src/FdemSolver.cpp))
et compteurs du run de référence : sur ~37 000 joints insérés, 76 % vont à
la rupture complète et **98,8 %** sont endommagés ou rompus. Presque aucune
insertion « gaspillée » : le critère (traction moyenne des deux éléments
adjacents ≥ enveloppe, éq. 7-8 de Yan 2023) tire exactement quand il doit.

Le mécanisme du nuage est ailleurs : autour d'un tunnel, la solution
élastique met **tout un anneau au-dessus de l'enveloppe simultanément** (le
rayon plastique classique). Sur un continuum lisse et homogène, un critère à
seuil par facette répond par une insertion en tapis — c'est sa réponse
*correcte* à ce champ.

## 2. Ce que dit la littérature (revue en 6 volets, `biblio_insertion/`)

- **Le critère de rockim est conforme aux fondateurs.** Camacho & Ortiz 1996
  évaluent la traction par partition nodale des forces (équivalent mécanique
  de notre moyenne des deux éléments) ; Pandolfi & Ortiz 2002 vérifient par
  facette, tous les N pas, **sans ordonnancement ni plafond** — toutes les
  facettes au seuil sont fracturées au même passage. Notre « tapis » n'est
  pas une déviation du schéma, c'est le schéma.
- **Mais les fondateurs n'ont jamais tourné en homogène** : la structure
  `Facet` de Pandolfi-Ortiz stocke une limite de traction PAR FACETTE depuis
  2002, et Zhou & Molinari 2004 ont formalisé le Weibull par facette comme
  « throttle physique ». Leurs blocs en compression viennent de résistances
  aléatoires + maillages non structurés + avance en pointe.
- **Nuances propres à rockim** relevées par l'audit : critère en OU
  (traction OU cisaillement) là où la littérature récente utilise une norme
  effective plus sélective ; pénalité d'activation 4 E/h en bas de la
  fourchette publiée (10–10⁴ E/h) ; et un écart de fidélité DORMANT — le
  défaut `shearEnvelope` tronque la branche de traction (fs = c au lieu de
  c − σn·tanφ pour 0 < σn < ft) ; la forme de l'article est `shearEnvelope =
  yang`, non posée dans les decks tunnel. Non testé (voir §9).

## 3. L'hétérogénéité ne change pas le régime — résultat négatif propre

Le remède canonique de la littérature, testé trois fois sur le cas tunnel
(mêmes decks, T = 0,7–1,0 s, comparaison à temps égal) :

| variante | fissures (t = 0,41 s) | blocs | mono-élément | côté moyen |
|---|---|---|---|---|
| homogène (référence) | 21 378 | 8 079 | 75,7 % | 0,27 m |
| Weibull m = 6, bruit blanc | 20 294 (−5 %) | — | — | — |
| corrélé θ = 1 m | 22 995 (+7 %) | 9 517 | 80,5 % | 0,24 m |
| corrélé θ = 4 m (t = 0,29 s) | −5 % vs réf. | −1,4 % | 75,0 % | 0,25 m |

Ni le tirage par joint, ni le champ corrélé court, ni les zones faibles de
19 éléments (θ = 4 m) ne délogent la granulation. θ = 1 m est même
contre-productif (nucléation facilitée sans coalescence accrue). Seul effet
robuste de θ = 4 : la brisure de symétrie du motif à λ = 1 (physique). Le
balayage θ = 2 a été sauté (encadré par deux cas défavorables).

## 4. Quatre instruments, quatre faits (outils dans `tunnel_edz/tools/`)

1. **`crack_clusters.py`** — le réseau est déjà PERCOLANT partout : une
   composante connexe contient 96-99 % des joints rompus, dans TOUS les
   runs. L'hypothèse « ça ne coalesce pas » est fausse ; c'est elle qui
   avait motivé le balayage de corrélation.
2. **`block_sizes.py`** — la bonne observable est le dual : la taille des
   blocs intacts. Verdict : **l'anneau est granulé à l'échelle de
   l'élément** (3 blocs sur 4 = un seul triangle, côté moyen 1,3 dx). La
   taille de bloc est fixée par le maillage, pas par la physique — défaut
   d'objectivité du MOTIF (l'énergie, elle, reste finie).
3. **`nucleation_vs_propagation.py`** — l'adaptatif nuclée plus qu'il ne
   propage : 43,7 % de propagation contre 58,9 % pour l'intrinsèque (trames
   tardives : 55-69 % contre 77-81 %). C'est la définition mécanique du
   nuage.
4. **`crack_coherence.py`** — surprise : les DEUX schémas sont au niveau du
   hasard directionnel (angle médian entre arêtes rompues voisines 53,1° /
   54,2° contre 55,5° pour un tirage aléatoire ; alignées < 30° : 23,2 /
   22,2 / 20,7 %). Les macro-fissures intrinsèques ne sont PAS localement
   colinéaires — elles zigzaguent d'arête en arête. Leur continuité est une
   organisation à grande échelle, pas une sélection d'orientation. Corollaire :
   « guider la direction » ajouterait une propriété qu'AUCUN des deux
   schémas ne possède.

## 5. Le contrôle intrinsèque — la vraie nature de la différence

Même maillage (tunnel_hs_iso, 106 298 tris), mêmes lois, mêmes CL, même
in-situ ; seul le schéma change (deck `tunnel_intrinseque.cfg`, excavation
décalée à 0,15 s pour absorber le tassement initial des joints). À temps
depuis le relâchement égal :

| | adaptatif | intrinsèque |
|---|---|---|
| joints rompus | 21 259 | 36 527 (+72 %) |
| blocs | 8 079 | 16 430 (×2) |
| mono-élément | 75,7 % | 84,9 % |
| plus gros blocs LIBRES | 37 m² | **318, 120, 114, 104, 68 m²** |
| propagation | 43,7 % | 58,9 % |
| prix | — | dt ÷ 1,9 ; E_eff = 0,95 E |

**Le schéma d'insertion ne détermine pas la finesse de la granulation — il
détermine la capacité à isoler de gros blocs.** L'intrinsèque granule PLUS
près de la paroi, mais ses 159 269 interfaces préexistantes offrent des
chemins continus qui percolent loin et bouclent des contours fermés autour
de grands volumes ; l'adaptatif ne crée de surface que là où le critère
tire, et ses fissures essaiment sans jamais fermer un contour. La fig. 11
intrinsèque montre les macro-fractures radiales continues jusqu'à 20-25 m
(U paroi 0,386 m contre 0,347 publié) ; `fig_blocs_compare.png` juxtapose
les trois régimes.

Le diagnostic initial de Fernando est donc validé dans sa conclusion (c'est
bien le schéma) mais pas dans son mécanisme (ce n'est ni l'instant ni la
direction d'insertion — c'est l'absence de chemins préexistants).

## 6. Premier remède implémenté : `insertionTipFactor` (rockim_p2)

Capacité opt-in (défaut 1,0 = bit-identique) : une facette adjacente à un
joint de D ≥ `insertionTipDamage` (défaut 0,5) voit son enveloppe divisée
par le facteur — correction assumée du fait que la moyenne d'éléments CST
écrase la singularité de pointe (2,7 éléments par zone cohésive en mode I).
Balayage complet sur le tunnel :

| facteur | 1,0 | 1,3 | **1,6** | 2,0 | cible intrinsèque |
|---|---|---|---|---|---|
| propagation | 43,7 % | 52,3 % | **57,4 %** | 55,3 % | 58,9 % |

**Optimum net à 1,6**, à 1,5 point de la cible ; au-delà, le relâchement
re-nucléé des sites parasites. Réserve honnête : le facteur augmente aussi
le NOMBRE total de rompus (+69 % à temps égal pour 1,6) ; le motif
(`fig_blocs_compare`, panneau central) montre des blocs polygonaux nettement
plus grands que la référence, mais la métrique de blocs n'a pas été
consolidée en chiffres pour tip16 — à faire avant d'en faire un réglage
recommandé.

## 7. Idées évaluées et écartées (avec la raison)

- **Zone d'ombre DFH géométrique dans le FDEM tunnel** : hors régime —
  kc = 738 m/s donne un rayon d'obscuration de 406 m sur la durée du run
  (4× le domaine) ; l'obscuration décrit un chargement en un coup, pas une
  cavité qui recharge en permanence. Reste pertinente pour la PERCUSSION
  (50 µs → 40 mm). Et en FDEM les fissures déchargent déjà mécaniquement —
  la superposer compterait le relâchement deux fois.
- **Mémoire de plan (direction figée à l'amorçage, à la DP-DFH)** : design
  posé (§4 l'a requalifiée : ce serait un AJOUT, aucun des deux schémas ne
  sélectionne de direction) ; non codée, l'étude ayant été close avant.
- **Course de vitesse insertion/onde** : réfutée par la mesure — 5 joints
  insérés par temps de traversée d'élément ; le déchargement a le temps
  d'agir, le tapis vient du rechargement permanent et de l'uniformité.
- **« Le maillage est trop grossier »** : réfuté trois fois — leur maillage
  est PLUS grossier (83 458 contre 106 298 éléments), le mode dominant
  (cisaillement, ℓ_ch = 1,56 m = 7,5 éléments) est résolu, et un bloc de
  1,3 dx signale un motif fixé par dx : raffiner ferait des blocs plus
  petits, pas plus grands.

## 8. Le détour DP-DFH (GELÉ le 2026-08-25, à la demande de Fernando)

Question : rockim sait-il refaire les cas de la thèse SANS cohésifs, en
éléments finis purs ? Acquis avant le gel :

- **`law = dpdfh` existait déjà** (portage fidèle de vumat_kstdfh.f, hash
  spatial bit-compatible, repère figé à l'amorçage, intégrateur
  d'obscuration en forme fermée — jamais géométrique).
- **Exo tunnel (banc 6) refait en FDEM continuum** : fissures RADIALES
  discrètes (~20 armes contre N_span = 10 chez Abaqus au même taux — même
  ordre, comptages non identiques), cœur dense, motif qualitativement
  conforme. C'est la démonstration que la loi de volume localise là où les
  cohésifs homogènes granulent.
- **`dfhPsiVar` ajouté** (opt-in) : découverte au passage — la carte dit
  ψ = 15° mais `vumat_hole.f` l'IGNORE et calcule ψ(p̄) = clamp(160,345 −
  0,213793·p̄ ; 0 ; 51,7) ; sous p̄ = 509 MPa l'écoulement est donc ASSOCIÉ.
  La constante #5 des cartes est morte. Porté à l'identique dans rockim.
- **`insertion = none` ajouté** (2D + 3D) : le continuum pur véritable.
  Motivé par un post-mortem complet : neutraliser les joints par ft = 10¹²
  est un DÉTONATEUR (un élément pulvérisé à D → DCAP se distord, sa
  contrainte aberrante franchit même 10¹², et les joints activés portent
  dnE = ft/pj = 8 cm : 53 J → −89 GJ en 10 µs). Deux bugs débusqués en
  testant la clé (liaison des nœuds et intégration des groupes gardées par
  `adaptive_`) — corrigés, suite 19/19.
- **Impact 3D sur bloc 200×200×120** : énergies saines, cratère de surface
  PROPRE (~12 mm à 40 µs, l'ordre du 10-12 mm de Yang) — mais l'endommagement
  ne s'arrête jamais (disque uniforme de 45 mm en fin de run) et l'insert ne
  rebondit pas (9,5 → 5,4 m/s). Cause identifiée : AUCUNE suppression
  d'éléments (dfhDeld = 1e9, convention percussion), la zone broyée reste
  une mousse porteuse. Le banc 6 supprimait à DELD = 0,98 — premier réglage
  à essayer à la reprise.
- **Tunnel EDZ en DP-DFH : échec instructif** — 96,8 % du domaine endommagé.
  Erreur de transposition : le c = 0,8 MPa de leur Table 1 est une cohésion
  DE JOINT ; en continu il implique UCS = 3,0 MPa quand leur massif fait
  25 MPa (Rc/σ₀ = 5). Il fallait c_macro ≈ 6,65 MPa. Correction d'une ligne,
  non relancée (gel).

## 9. Suites classées, à la reprise

1. Consolider la métrique de blocs sur tip16/tip13 (les chiffres manquent) ;
   si les blocs suivent la propagation, faire de `insertionTipFactor = 1.6`
   un réglage documenté du cas tunnel.
2. Tester `shearEnvelope = yang` sur le tunnel (écart de fidélité dormant,
   §2) — un run, comparaison à temps égal.
3. Mémoire de plan à l'amorçage (design §7) — la capacité que ni
   l'adaptatif ni l'intrinsèque ne possède.
4. DP-DFH à la reprise du gel : (a) impact avec dfhDeld = 0,98 (le cratère
   par suppression, comme le banc 6) ; (b) tunnel avec la carte macro
   corrigée ; (c) l'horloge d'obscuration comme plancher de D des joints
   FDEM en percussion.
5. Le portage 3D de `insertionTipFactor` est fait (rockim_p2) mais n'a
   jamais tourné sur l'impact — le banc d'essai naturel est la pulvérisation
   A/B de spec 005.

## 10. Inventaire

**Runs** (tous archivés, VTU locaux, métriques + figures sur le drive) :
référence 0,25 s / stabilisation 1,0 s / Weibull 1,0 s / corrélé θ = 1
(tué à 60 %) / θ = 4 complet / intrinsèque (tué à 96 %) / tip13-16-20
complets / balayage λ 4 runs / exo tunnel dpdfh ×2 / impact 3D ×4 / tunnel
dpdfh (échec carte). Coût total ≈ 60 h machine.

**Capacités ajoutées à rockim** (toutes opt-in, défauts bit-identiques,
suite 19/19 à chaque étape — principe VIII respecté) : `insertionTipFactor`
+ `insertionTipDamage` (2D + 3D), `insertion = none` (2D + 3D, + 2 bugs de
liaison corrigés), `dfhPsiVar` + `dfhPsi0/KPsi/PsiMax`.

**Outils de mesure** (`tunnel_edz/tools/`) : `crack_clusters.py`,
`block_sizes.py`, `nucleation_vs_propagation.py`, `crack_coherence.py`,
`fig_blocs_compare.py`, `fig_wang11.py`, `fig_wang11_blocs.py`,
`fig_lambda.py`, `wall_history.py` ; (`bench_impact/tools/`)
`fig_surface_i3d.py`, `fig_fp_i3d.py`, `fig_vitesse_i3d.py`.

**Biblio** : `biblio_insertion/` (6 volets, sources marquées
[VERIFIE]/[MEMOIRE]) + `tunnel_edz/exo_tunnel_inventaire.md`.
