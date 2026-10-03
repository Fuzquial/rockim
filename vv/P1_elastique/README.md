# P1.2 et P1.3 : élasticité linéaire de l'élément en 3D

## Objectif

Vérifier que l'élément tétraédrique linéaire de rockim, en `fem3d` et en `fdem3d` continu
(insertion adaptative sans rupture), donne la solution statique de l'élasticité linéaire sur deux
cas à solution exacte : un champ homogène (patch test, P1.3) et le champ de Kirsch autour d'un trou
circulaire (P1.2). Le premier doit être reproduit à la précision machine près, le second converger
en maillage avec un ordre mesuré. Les deux cas reprennent ceux de Lei et Rougier (2016), cas 1 du
cube et plaque trouée de rayon a = 0,1 m.

## P1.3, patch test : cube sous chargement biaxial

### Problème

Cube de côté L = 1 m, matériau ρ = 2 650 kg/m³, E = 60 GPa, ν = 0,25. Traction morte de
compression PX = 1 MPa sur la face `xmax` et PY = 0,5 MPa sur la face `ymax`, rouleaux sur `xmin`
(u_x = 0), `ymin` (u_y = 0) et `bottom` (u_z = 0), face `top` libre. Lei et Rougier chargent les
quatre faces latérales à p = 1 000 Pa ; nous retenons PX différent de PY pour rompre la symétrie,
et une charge mille fois plus forte pour que le déplacement (environ 15 µm) soit lu dans les VTU
avec une résolution relative de 3·10⁻⁸ (coordonnées écrites à 12 chiffres).

### Solution exacte

Le champ est homogène : σ_xx = −PX, σ_yy = −PY, σ_zz = 0, et

ε_xx = (σ_xx − ν σ_yy)/E, ε_yy = (σ_yy − ν σ_xx)/E, ε_zz = −ν (σ_xx + σ_yy)/E,

u(x, y, z) = (ε_xx x, ε_yy y, ε_zz z), origine au coin `c000`.

Un champ de déplacement linéaire appartient à l'espace des tétraèdres linéaires : l'erreur de
discrétisation est nulle sur tout maillage. La formulation co-rotationnelle (déformation de Biot,
P = R σ) est exactement linéaire pour un étirement pur sans rotation, ce qui est le cas ici.

### Métriques et critères

- Erreur nodale : max |u − u_ex| / max |u_ex| sur tous les nœuds (toutes les copies en `fdem3d`),
  u lu comme différence des positions du dernier VTU et du premier. Critère : 1·10⁻⁶ sur chaque
  maillage et chaque mode.
- Déplacements moyens de groupes U_xmax_x, U_ymax_y, U_top_z de `history.csv`, contre ε_xx L,
  ε_yy L, ε_zz L. Critère : 2·10⁻⁵ (le fichier n'a que 6 chiffres significatifs).
- Trois maillages non structurés : h = 0,25, 0,125 et 0,0625 m (L/h = 4, 8, 16).

## P1.2, Kirsch en 3D

### Problème

Quart de plaque [0, W] × [0, W] × [0, t], W = 2 m, percé d'un quart de trou de rayon a = 0,1 m
centré à l'origine (W/a = 20). Symétries par rouleaux sur `xmin` et `ymin`, déformation plane par
rouleaux en z sur `bottom` et `top`, traction morte σ_x^∞ = SX = −1 MPa sur `xmax` et λ SX sur
`ymax`, avec λ = 0 et λ = 0,5. Le trou (`hole`) est libre.

Maillage : maillage 2D graduel extrudé en deux couches de tétraèdres, taille
h(r) = h_p (1 + 1,5 (r − a)/a), plafonnée, avec h_p = s a/10 à la paroi et t = 2 h_p. Trois
facteurs d'échelle s = 1, 0,5 et 0,25 (h_p = a/10, a/20, a/40 ; 4 620, 18 036 et 73 098 tets).
Les nœuds du plan médian sont libres en z, le problème reste donc tridimensionnel.

### Solution exacte

Kirsch (Jaeger et Cook), pour une traction uniaxiale s selon x, θ mesuré depuis x :

σ_rr = s/2 (1 − a²/r²) + s/2 (1 − 4a²/r² + 3a⁴/r⁴) cos 2θ

σ_θθ = s/2 (1 + a²/r²) − s/2 (1 + 3a⁴/r⁴) cos 2θ

σ_rθ = −s/2 (1 + 2a²/r² − 3a⁴/r⁴) sin 2θ

et en déformation plane, avec G = E/(2(1 + ν)) et κ = 3 − 4ν :

u_r = s/(4G) [ r ((κ − 1)/2 + cos 2θ) + (a²/r)(1 + (1 + κ) cos 2θ) − (a⁴/r³) cos 2θ ]

u_θ = −s/(4G) [ r + (1 − κ) a²/r + a⁴/r³ ] sin 2θ.

La charge selon y s'obtient en remplaçant θ par θ − π/2, et les deux se superposent. À la paroi en
(0, a), σ_xx = (3 − λ) SX : facteur de concentration 3 pour λ = 0, 2,5 pour λ = 0,5. La plaque
finie (W/a = 20) écarte la solution de celle de la plaque infinie d'un terme d'ordre (a/W)²,
soit quelques dixièmes de pour cent.

### Mise en œuvre des mesures

- `fem3d` : sondes `probes` (contrainte constante par tétraèdre, traction positive) le long de
  l'axe θ = 90° (25 points, r de a à 4a, resserrés à la paroi) et de l'axe θ = 0 (10 points), à
  z = t/4. Le journal imprime le centroïde de chaque tétraèdre sondé ; la solution exacte est
  évaluée au point sondé et au centroïde.
- `fem3d` et `fdem3d` : sommets de paroi `point.wA` en (a, 0, t/2) et `point.wB` en (0, a, t/2),
  et à 2a (`point.mA`, `point.mB`), chargés d'une force nulle pour que `history.csv` écrive leur
  déplacement. La position réelle du sommet retenu est relue dans le journal.
- `fdem3d` n'a pas de sondes ; il est jugé sur les déplacements, sur les deux maillages les plus
  grossiers (son coût par pas est environ 25 fois celui de `fem3d` et son pas 5 fois plus petit).

### Métriques et critères (maillage le plus fin de chaque variante)

| métrique | définition | critère |
|---|---|---|
| kt_centroid | σ_xx du tétraèdre de paroi en (0, a) contre l'exact à son centroïde | 2 % |
| kt_raw | même contrainte contre (3 − λ) SX | 5 % |
| line_centroid | erreur L2 relative de σ_xx le long de θ = 90°, r de a à 4a, exact aux centroïdes | 3 % |
| u_wall | déplacement des sommets wA (u_x) et wB (u_y) | 1 % |
| ordre | pente log-log de l'erreur L2 brute de σ_xx contre h_p (théorique 1, contrainte constante par élément) | 0,8 au moins |
| stationnaire | variation relative de U_wB_y sur les 10 derniers pour cent de T | 1·10⁻⁴ |

## Régime quasi statique

`scenario = loads`, charges en `amplitude.<g> = ramp` (montée en cosinus sur 5 traversées L/c de
l'onde P, c = 5 212 m/s), `dampingLocal = 0,5`. Durée : 60 traversées du cube, 40 traversées de la
plaque (L = W). Au smoke test, U atteint sa valeur finale à 1·10⁻⁵ près après 40 traversées (cube)
et 25 (plaque).

## Clés rockim utilisées

`mode`, `law = elastic` (fem3d), `insertion = adaptive`, `ft`, `cohesion`, `frictionDeg`, `Gf`
(fdem3d), `scenario = loads`, `mesh = file`, `meshFile`, `rho`, `E`, `nu`, `dampingLocal`, `T`,
`frames`, `historyFlush`, `fix.<g>`, `traction.<g>`, `amplitude.<g> = ramp`, `point.<g>`,
`force.<g>`, `probes` (fem3d), `outputDir`.

## Lancer

Depuis la racine rockim :

    python3 vv/P1_elastique/p1_elastique.py prepare        # maillages et decks
    python3 vv/run_queue.py P1_elastique --slots 4 --threads 2
    python3 vv/P1_elastique/p1_elastique.py analyse        # resultats.json, fig_*.pdf et .png

ou `p1_elastique.py all --threads 2` en série. Option `--only <variantes>` pour restreindre,
`--smoke` pour le seul maillage grossier (sorties dans `out_smoke/`, analyse dans
`resultats_smoke.json`).

Figures : `fig_patch` (erreur nodale contre h), `fig_kirsch_profil` (σ_xx le long de θ = 90°),
`fig_kirsch_convergence` (erreurs contre h_p/a).

## Smoke test (2026-10-04, OMP_NUM_THREADS = 1, maillage le plus grossier)

| run | résultat | durée |
|---|---|---|
| patch_fem3d, h = 0,25 m | erreur nodale 7,5·10⁻⁸ ; U_xmax_x −1,45833·10⁻⁵ m (exact −1,458333·10⁻⁵) ; résidu du bilan 1·10⁻¹¹ % | 0,16 s |
| patch_fdem3d, h = 0,25 m | erreur nodale 7,5·10⁻⁸, mêmes U | 4,6 s |
| kirsch_fem3d λ = 0, s = 1 | Kt = 2,988 (−0,40 %) ; contre le centroïde +5,1 % ; L2 brute 8,6 %, aux centroïdes 5,6 % ; u_wA −0,29 %, u_wB +0,19 % | 45 s (60 traversées) |
| kirsch_fem3d λ = 0,5, s = 1 | Kt = 2,500 (0,00 %) ; contre le centroïde +4,2 % ; u_wB −0,98 % | 45 s |
| kirsch_fdem3d λ = 0, s = 1 | le deck tourne (T réduit à 2·10⁻⁴ s, non significatif), 0 joint inséré | 38 s |

Sur le maillage grossier, la contrainte du tétraèdre de paroi est plus proche de la valeur à la
paroi que de la valeur exacte à son centroïde : le critère kt_centroid ne se jugera qu'au
maillage le plus fin.
