# tools/depouiller — figures d'un run a OUTIL MAILLE, depuis un .npz

Les outils historiques (`screen3d/rapport.py` et cie) lisent `toolFz`,
`toolVz`, `toolZ` : les colonnes de l'outil ANALYTIQUE. Dans les cas a outil
maille (`toolShape = none` + `groupVel`), elles sont NULLES — c'est aussi
pourquoi le resume affiche `peak tool force : 0 N`. Ces cinq scripts lisent
l'instrumentation WP3 a la place, depuis le `.npz` produit par
`tools/pack_run.py`.

    python tools/pack_run.py out_pulv_coulomb bench_impact/donnees/P1.npz
    python tools/depouiller/rapport_npz.py   bench_impact/donnees/P1.npz P1
    python tools/depouiller/coupe_npz.py     bench_impact/donnees/P1.npz P1
    python tools/depouiller/surface_npz.py   bench_impact/donnees/P1.npz P1
    python tools/depouiller/film_npz.py      bench_impact/donnees/P1.npz P1
    python tools/depouiller/filmhaut_npz.py  bench_impact/donnees/P1.npz P1

| script | figure |
|---|---|
| `rapport_npz` | F-p par DEUX routes (jauge du bit / quantite de mouvement), vitesses et rebond, pulverisation datee, canaux d'energie |
| `coupe_npz` | coupe verticale exacte : roche coloree par `bulkD`, fissures en TRAITS |
| `surface_npz` | facies de surface : aretes DEBOUCHANTES en traits + rose azimutale, pulverisation, cratere |
| `film_npz` | film de la coupe |
| `filmhaut_npz` | film du facies de surface |

## Conventions etablies (et pieges deja payes)

* **Force** : deux routes independantes. La jauge `szz_bit x A_bit`
  (A = pi R^2, R = 15 mm) est l'analogue de la jauge de deformation du banc ;
  la route inertielle est d/dt de la quantite de mouvement des DEUX corps
  brases. Elles different du transit d'onde dans le bit (~0,13 m d'acier).
* **Penetration** : sous la SURFACE LIBRE, jamais depuis la position
  initiale de l'outil — le maillage laisse un jeu (0,200 mm sur `impact_pulv`,
  et non 0,02 mm comme le dit un commentaire perime du deck). Le negliger
  surestime l'enfoncement d'autant.
* **Surface libre** : max des NOEUDS, pas des barycentres d'elements (0,43 mm
  d'ecart sur un maillage a 2 mm).
* **Lissage** : padding MIROIR, jamais zero. Attention : le miroir force la
  derivee a s'annuler au bord — pour une pente en fin de serie, ajuster sur
  une fenetre interieure.
* **Fissures** : rendu ELEMENTS. En coupe, des traits (intersection exacte du
  triangle et du plan) ; en vue de dessus, la seule ARETE DEBOUCHANTE, c'est-
  a-dire les deux sommets qui sont dans le plan de la surface libre en
  configuration de reference. Projeter les facettes entieres donne un tapis
  illisible.
* **Les fissures suivent la matiere** : le `.npz` ne garde les joints qu'a la
  derniere frame. Les dessiner figes fait flotter des traits dans le vide des
  que la matiere est ejectee. Le maillage etant DISCONTINU, les sommets de
  joint sont EXACTEMENT des noeuds d'elements : un hachage spatial les apparie
  (8 004 / 8 004, erreur 0,0000 mm sur P1) et la geometrie se retrouve a
  n'importe quelle frame.
* **Cratere** : separer le relief EN PLACE (|dz| < 3 mm) des ECLATS EN VOL,
  sinon un fragment a +10 mm ecrase l'echelle et la cuvette (~1 mm) disparait.
* **Datation des frames** : uniforme (verifie par appariement de la position
  d'outil ; l'appariement devient ambigu apres le rebroussement, la meme
  position etant atteinte a l'aller et au retour).
