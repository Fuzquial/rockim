# Bancs d'insertion des joints (rockim 2D, `mode = fdem`) — critères d'acceptation

Rédigé le 2026-10-07, avant les calculs des bancs 2, 3 et 4. Binaire
`/home/user/rockim/build_contact/rockim` (commit b869dae). Unités SI dans les decks.

> Écart de procédure à déclarer. Les 24 calculs du banc 1 (0,2 à 0,4 s chacun) ont été
> lancés une première fois avant la rédaction de ce fichier, par erreur d'ordre. Les critères
> du banc 1 ci-dessous sont ceux fixés à la conception du mini-maillage (tolérances de 2 %
> sur le pic et l'énergie, 3 % sur le pic en mode mixte, 5 % d'écart adaptatif/intrinsèque) ;
> ils n'ont pas été retouchés après lecture des chiffres. Le banc 1 est relancé après
> l'écriture de ce fichier ; les chiffres retenus sont ceux de la relance.

## Notations et références

- f_t = 10 MPa, c = 25 MPa, φ = 0° (banc 1) ou 40° (bancs 2-4), G_I = Gf, G_II = gfShearFactor · Gf.
- Longueur de zone cohésive retenue : ℓ_cz = ℓ_ch = E·G_I/f_t² (longueur de Hillerborg ; celle
  de Rice, 9π/32 ℓ_ch = 0,88 ℓ_ch, est du même ordre). Bancs 2-4 : E = 50 GPa, G_I = 40 J/m²,
  donc ℓ_cz = 20 mm et h = ℓ_cz/2, /4, /8 = 10 ; 5 ; 2,5 mm.
- Pénalité : p = facteur·E/h_loc, h_loc = moyenne des diamètres inscrits des deux triangles
  (`mesh = file` passe par la branche « voronoi » de `assignJointProps`). En adaptatif le
  facteur est `insertionPenaltyFactor` (raideur des joints insérés seulement), en intrinsèque
  `jointPenaltyFactor`.
- Référence de la loi adaptative : Yan, Zheng & Wang 2023 (IJRMMS 169:105439) — critère
  d'insertion éq. 7-8 (σ ≥ f_t ou τ ≥ f_s, contrainte moyenne des deux éléments, fig. 3),
  f(D) éq. 11, D en mode I/II/mixte éq. 12, 14, 16, énergies éq. 13 et 15 (l'aire sous σ(o)
  vaut G_I, sous τ(s) vaut G_II). Leur banc brésilien (fig. 12-13 : pénalité 1 E à 100 E ;
  fig. 15-16 : h = 1,5 / 1,0 / 0,75 mm) est la référence qualitative de la convergence en
  maille de l'adaptatif.

## Banc 1 — joint isolé (mini-maillage de 8 CST, `scripts/mesh_b1.py`)

Bande 10 mm de large, joint traversant ; ν = 0 et mors latéralement libres : le champ est
uniaxial et uniforme, exact pour des CST. Les six arêtes parasites sont presque verticales :
à la rupture du joint cible elles restent sous 30 % de f_t et 55 % de c (rapports imprimés
par `mesh_b1.py`). Trois trajets : traction pure (joint horizontal) ; cisaillement pur
(joint à 45° sous compression uniaxiale, φ = 0 donc f_s = c, aucun frottement, σ_n = −τ) ;
mode mixte (joint à 45° en traction, σ_n = τ). Adaptatif contre intrinsèque, pénalité 10 et
100 E/h, adoucissement `linear` et `yan`.

| critère | mesure | seuil |
|---|---|---|
| B1-a pic mode I | max σ_n au joint / f_t | 1 ± 0,02 |
| B1-b pic mode II | max τ au joint / c | 1 ± 0,02 |
| B1-c pic mixte | σ_n au joint au pic / f_t (critère « ou », branche de traction car τ = σ_n < c) | 1 ± 0,03 |
| B1-d énergie mode I | travail cohésif final / (G_I·L·t), L = 10 mm | 1 ± 0,02 |
| B1-e énergie mode II | travail cohésif final / (G_II·L·t), L = 14,14 mm | 1 ± 0,02 |
| B1-f énergie mixte | travail cohésif final, aucune référence exacte (dépend de la loi codée) | entre G_I·L et G_II·L ; écart adaptatif/intrinsèque ≤ 5 % à pénalité et loi égales ; écart 10 E/100 E ≤ 5 % |
| B1-g courbe | σ_n(δ_n) mode I : adoucissement jusqu'à δ_c ; δ_c attendu = 2 G_I/f_t = 14 µm (linear) | rupture (σ < 1 % f_t) à δ_c ± 10 % |
| B1-h pénalité | pic et énergie entre 10 E et 100 E | écart ≤ 2 % (mode I et II) |

## Banc 2 — barre en traction quasi statique, maillage non structuré

Barre 20 × 40 mm (gmsh, `tools/make_unstructured_mesh.py box2d`), mors latéralement libres
(champ uniaxial uniforme : aucun site d'amorçage imposé), ν = 0,25, φ = 40°, G_II = 10 G_I,
`jointSoftening = yan`, pullV = 0,02 m/s, rampe 20 µs, Cundall 0,7, jointXi = 0. Pas de
retour élastique brutal attendu (énergie élastique au pic 0,4 J/m < G_I·W = 0,8 J/m).
Trois tailles h = 10 ; 5 ; 2,5 mm (plus 1,25 mm si le coût le permet). Adaptatif
(`insertionPenaltyFactor = 100`) contre intrinsèque (`jointPenaltyFactor = 100`), contact
`potential` sans correctif (`contactCandidates = active`, `gcBirth = ramp`,
`potForceExact = false`) contre contact corrigé (`vertex`, `offset`, `true`).

Mesures : charge de rupture σ_pic ; joints insérés (adaptatif) ou endommagés D > 0
(intrinsèque) ; joints rompus (D ≥ 1) ; N_fiss = composantes connexes des joints rompus
(partage d'un sommet de la configuration initiale) ; fissure principale = plus grande
composante ; énergie cohésive finale W.

| critère | seuil |
|---|---|
| B2-a pic | σ_pic / f_t ∈ [0,95 ; 1,10] pour chaque calcul |
| B2-b convergence du pic | écart entre les deux maillages les plus fins ≤ 3 % |
| B2-c énergie | W / (G_I·W_barre·t) ∈ [1,0 ; 1,5] (1 = fissure droite unique, la tortuosité majore) |
| B2-d convergence de l'énergie | écart entre les deux maillages les plus fins ≤ 10 % |
| B2-e une fissure | la fissure principale traverse la barre (étendue en x ≥ 0,95 W) ; hors d'elle, ≤ 10 % des joints rompus |
| B2-f nucléation diffuse | part des joints insérés (adaptatif) ou endommagés (intrinsèque) dont le milieu est à plus de 2h de la fissure principale : ≤ 20 % = « localise », > 50 % ou N_fiss ≥ 5 = « nuclée partout », entre les deux = « mixte » |
| B2-g contact corrigé | en traction quasi statique les lèvres se touchent peu : σ_pic et W changent de ≤ 1 % ; travail du contact ≤ 0 avec le correctif |

## Banc 3 — plaque à entaille latérale (SENT), traction

Plaque 40 × 80 mm, entaille rectangulaire depuis le bord gauche à mi-hauteur, longueur
a = 10 mm, ouverture 1 mm (ligament 30 mm). Même matériau et même pilotage que le banc 2
(contact corrigé). h = 2,5 et 1,25 mm. Adaptatif (100 E) contre intrinsèque (100 E). Pas de
scénario de flexion en 2D : la poutre entaillée n'est pas faite.

Part de propagation (adaptation de `tunnel_edz/tools/nucleation_vs_propagation.py`, avec
les instants `tBreak` de la trame finale au lieu des différences de trames) : un joint qui
rompt est « propagation » s'il partage un sommet initial avec un joint rompu strictement
avant lui, « nucléation » sinon ; même chose à l'insertion avec `tInsert`.

| critère | seuil |
|---|---|
| B3-a amorçage | le premier joint rompu est à moins de 2h de la pointe d'entaille |
| B3-b traversée | la fissure principale part de la pointe (≤ 2h) et atteint le bord droit (x ≥ W − h) |
| B3-c trajectoire | écart moyen \|y − y_entaille\| ≤ 0,1 W = 4 mm, écart max ≤ 0,2 W |
| B3-d longueur | longueur de la fissure principale / ligament ∈ [1,0 ; 1,3] |
| B3-e propagation | part de propagation à la rupture ≥ 80 % ; joints rompus hors fissure principale ≤ 10 % |

## Banc 4 — sensibilité de l'adaptatif au critère (banc 2, h = ℓ_cz/4 = 5 mm, contact corrigé)

Base : `insertionCriterion = or`, `facetAverage = arith`, `insertionHoldSteps = 1`,
`insertionTipFactor = 1`. Une clé changée à la fois : `elliptic` ; `volume` ; `max` ;
`insertionHoldSteps = 5` ; `insertionTipFactor = 1,6`. Répété à h = 2,5 mm si le coût le
permet. Attendus écrits avant calcul (ce sont des prédictions, pas des seuils de validité) :

| variante | prédiction |
|---|---|
| elliptic | insère plus tôt (traction et cisaillement se cumulent) : pic ≤ base, ≥ 0,9 base ; plus de joints insérés |
| volume | quasi neutre sur un maillage quasi uniforme : pic à ± 1 %, nombre de rompus à ± 20 % |
| max | insère plus tôt que la moyenne : pic ≤ base ; plus de joints insérés |
| hold 5 | en quasi statique sans oscillation : pic à ± 1 % ; motif inchangé |
| tip 1,6 | favorise la propagation : N_fiss ≤ base, part de propagation ≥ base, pic à ± 2 % (l'amorçage est inchangé) |

Verdict « sensible » si le pic varie de plus de 2 % ou le nombre de rompus, N_fiss ou W de
plus de 10 %.

## Règles d'exécution

OMP_NUM_THREADS = 2, nice -n 10, `timeout 1200` (20 min), au plus deux calculs en parallèle
(un seul tant que le rejeu tunnel tourne). Pour chaque calcul : durée, commande, chiffres
extraits. Les `.vtu` sont supprimés après dépouillement.

## Addendum 1 (même jour, après le banc 1, AVANT la variante stable du banc 1 et avant les bancs 2-4)

1. Le banc 1 tel que conçu est instable après le pic : la bande de 80 mm stocke au pic
   f_t²·H/(2E) = 80 J/m² > G_I = 70 J/m² (condition 1D de retour brutal H > 2 ℓ_ch, ici
   ℓ_ch = 70 mm) ; le joint s'ouvre de 0 à 18 µm entre deux trames (66 µs) pour 2,7 µm de
   course des mors. B1-g (ouverture de rupture) n'y est pas mesurable ; les critères d'énergie,
   qui ne dépendent que de la loi, restent jugeables. Variante stable ajoutée, `b1s` :
   mêmes maillages, E = 200 GPa (stockage 20 J/m² en mode I, 0,57 G_II·L en mode II),
   400 trames, mêmes seuils B1-a à B1-h. Les ratios des arêtes parasites ne dépendent pas
   de E.
2. Erreur de calcul dans la description du banc 2 : l'énergie élastique de la barre au pic
   vaut σ²/(2E)·W·H = 0,8 J/m (et non 0,4) pour H = 40 mm, soit exactement G_I·W : limite
   de stabilité. Correctif de géométrie avant tout calcul du banc 2 : barre 20 × 30 mm
   (0,6 J/m < 0,8 J/m, H = 1,5 ℓ_ch). Les seuils B2-a à B2-g sont inchangés.

## Addendum 2 (après `b1s`, avant les calculs `b1s_mixc` et les bancs 2-4)

Le trajet « mixte » avec mors latéralement libres ne pilote pas le rapport ouverture/glissement :
le bloc supérieur choisit lui-même sa direction (ouverture quasi pure sous `yan`, glissement
quasi pur sous `linear`). Variante `b1s_mixc` : même deck que `b1s_mix`, mors bloqués
latéralement (`gripLateralFree = false`, le champ reste uniforme avant le pic car ν = 0) ;
le déplacement relatif est alors vertical, donc δ_n = δ_s (trajet proportionnel à 45°).
Critères B1-c et B1-f inchangés, appliqués à cette variante.

## Addendum 3 (pendant le banc 2, avant les calculs b2L et b4)

À T = 1,25 ms (course 25 µm), les barres ADAPTATIVES h = 10 et 5 mm ne sont pas séparées :
contrainte résiduelle 0,70 et 0,18 f_t, et la contrainte remonte après une première chute
(plusieurs fissures partielles à des hauteurs différentes). Les barres intrinsèques sont
séparées (≤ 0,06 f_t). Pour juger B2-c/B2-d/B2-e sur une rupture complète, les calculs
adaptatifs du banc 2 sont refaits à T = 3 ms (course 60 µm), série `b2L` ; le banc 4 (tout
adaptatif) passe aussi à T = 3 ms. Les calculs à 1,25 ms sont conservés et rapportés. Une
barre est dite « séparée » si σ_fin < 0,05 f_t ; sinon B2-c à B2-e sont déclarés
« non jugeables » pour ce calcul. Seuils inchangés.

## Addendum 4 (pendant le banc 3)

Coût du contact corrigé : `b2_h1.25_adapt_corr` (124 s sans correctif) et `b3_h2.5_adapt`
ont atteint la limite de 20 min (arrêt à 89 % et 91 % de T, après la rupture complète : la
dernière trame est dépouillée, le calcul est déclaré « interrompu »). La plaque SENT h = 2,5 mm
est séparée à t ≈ 0,75 ms ; les deux calculs h = 1,25 mm du banc 3 sont donc réduits à
T = 1,0 ms (course 20 µm). `b2L_h1.25_adapt_corr` (T = 3 ms) n'est pas lancé : il dépasserait
20 min de façon certaine. Seuils inchangés.

## Addendum 5 (banc 3, h = 1,25 mm)

Avec le contact corrigé, `b3_h1.25_adapt` n'avait parcouru que 24 % de T = 1 ms après 13 min :
arrêté. Le banc 2 ayant montré que le correctif ne change ni le pic ni le motif en adaptatif
(et l'énergie de moins de 0,2 %), les deux calculs h = 1,25 mm du banc 3 sont relancés avec
le contact `potential` sans correctif (gain de coût ×5 à ×10 mesuré au banc 2). Seuils inchangés.
