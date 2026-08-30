# Guo 2014, chapitre 6 — Dolosse et Core-Loc : la validation dynamique du modèle 3D

*Fiche du 2026-08-28. Source : thèse Guo (Imperial College), p. 212-260, lues
intégralement sur PDF fourni par Fernando. Publié : Guo, Latham & Xiang,
Computers & Structures 146 (2015) 117-142.* **[V]** = lu dans ce PDF ;
**[D]** = analyse/proposition personnelle.

## 1. Les deux essais Dolosse (§6.3) — validés contre Burcharth 1981b **[V]**

Dolos H = 1 m (réduit du 2,32 m de Burcharth, rapports a/H = 0,36, b/H = 0,20,
c/H = 0,056 conservés ; 0,43 t), béton non armé : ρ 2340, E 26 GPa, ν 0,2,
**ft 3,3 MPa, c 16,5 MPa, φ 30°, G_f 50 J/m²** (Table 6.2). Socle
viscoélastique SANS rupture (« estimation conservative »), μ = 0,6
béton-béton, 0,1 béton-acier. Maillages NON STRUCTURÉS ~2,9 cm : 170 282 tets
(chute), 187 294 (pendule), dt = 5e-8 s. Une seule frappe (pas de fatigue au
modèle → fissure partielle, pas de scission complète ; Burcharth frappait
6-8 fois).

- **Chute (marteau)** : l'unité pivotée à 26° sur l'arête d'un fluke (angle
  choisi pour retrouver les contraintes du 2,32 m à ~8°), lâchée sous
  gravité, impact plat-plat du fluke droit. Séquence numérique : amorçage au
  **coin haut** tige-fluke (~0,5 ms) → milieu de tige → coin bas ; tout se
  fige au rebond (~2,4 ms).
- **Pendule (missile)** : cylindre acier (masse = Dolos/5), h = 0,2 m →
  v0 = √(2gh) = 2 m/s, impact à H/12 du bout du fluke. Séquence : écrasement
  local d'UN élément au point d'impact (cisaillement), puis première fissure
  de traction au **coin bas** tige-fluke (~0,8 ms).

**La signature discriminante** : chute et pendule s'amorcent en sens
OPPOSÉS (haut→bas vs bas→haut), et les deux sens sont ceux des essais
physiques de Burcharth ET des observations de site (Crescent City, Myrick &
Melby 2005 : fractures fluke-tige et mi-tige). C'est une cible qualitative
robuste — invariante d'échelle, de maillage fin, de calage précis.

**Pré-fissures de surface** (1 taille d'élément aux congés, comme Burcharth
en observait au démoulage) : les nouvelles fissures partent des POINTES des
pré-fissures (localisation déplacée), mais motif global et étendue quasi
identiques à l'intact → « surface cracks had negligible influence on the
strength » — la conclusion expérimentale de Burcharth, retrouvée.

## 2. Core-Loc multi-corps (§6.4) — le patron d'audit énergétique **[V]**

5 unités Core-Loc pleine échelle (3,31 m, 18,72 t) posées à 16 mm au-dessus
d'une pente, lâchées sous gravité, T = 0,1 s ; 98 734 tets (0,15 m par
unité), dt = 2e-7 s. Chaque unité est dépouillée en énergie (E_p éq. 6.2 +
E_k éq. 6.3) **contre un run jumeau purement viscoélastique** :

- unité 4 (vol libre, aucune collision) : courbes rupture et élastique
  **exactement superposées**, E_p + E_k conservée — le zéro de l'instrument ;
- unité 5 (collisions sans rupture) : deux crans d'énergie = deux collisions
  (l'amortissement visqueux ∝ gradient de vitesse) ; petit RETARD de la
  courbe avec joints — la complaisance ajoutée par les éléments joints,
  assumée et expliquée ;
- unités 1-3 (rompues) : chaque événement de fissuration = perte PERMANENTE
  visible sur E_p + E_k là où l'élastique restitue ; toutes les fissures
  sont déclenchées par des collisions.

Diversité des modes en un seul run : grandes fractures de traction nettes
(l'unité 3 quasi coupée en deux), écrasement local aux points de contact
(« down to element size » — un maillage plus fin l'améliorerait, aveu),
cisaillement qui décapite l'unité 1. Limites avouées : pas de fatigue/cyclage
(pas de rocking), chargement de houle futur (proxy cyclique ou couplage CFD
ch. 5), CPU.

## 3. Ce que rockim peut en faire **[D]**

1. **bench_dolos** (2 configs, tier full) : géométrie fig. 6.2 entièrement
   paramétrée (H, a, b, c + congés — gmsh CSG), matériau Table 6.2, chute
   26° ET pendule 2 m/s. Critères PASS qualitatifs et robustes :
   (a) chute : premier amorçage au coin HAUT tige-fluke, puis mi-tige ;
   (b) pendule : premier amorçage au coin BAS — l'INVERSION est le critère ;
   (c) une frappe = fissure partielle, pas de scission ;
   (d) variantes pré-fissurées : amorçage déplacé aux pointes, motif global
   conservé. Coût estimé à l'échelle rockim : ~50k steps × ~0,3 s ≈ 4-5 h
   par essai à maillage comparable — tier full ; une version réduite
   (H 0,3 m, ~30k tets) pour un contrôle qualitatif rapide.
2. **L'audit jumeau rupture/élastique par corps** : rockim a déjà les corps
   (WP2/WP3 : groupBond, jauges par corps) et les bilans d'énergie ; le
   patron « même deck, insertion coupée (elastic twin), E_p + E_k par corps
   superposées » est un repère de non-régression puissant et peu coûteux —
   l'unité 4 de Guo est exactement notre zeroload multi-corps, en gravité.
3. La complaisance des joints (retard unité 5) est l'argument SOURCE à citer
   quand on compare des chronologies intrinsèque vs adaptatif : Guo
   documente le biais de raideur des joints intrinsèques sur un cas propre.

## 4. Bilan de la moisson de thèse

Avec ce chapitre, l'extraction utile de Guo 2014 est COMPLÈTE : ch. 2 (loi
de joint — le correctif coulomb en vient), §2.3.4-2.3.5 (contact, intégration,
dt), §2.4-2.5 (sensibilité au maillage — la règle h ≤ ℓ/3), §3.5 (polyaxial —
bench_polyaxial), ch. 5 (couplage fluide — registre des faces mouillées pour
la spec 004), ch. 6 (Dolosse/Core-Loc — bench_dolos et l'audit jumeau).
Restent volontairement non extraits : ch. 4 (remaillage adaptatif, intérêt
tunnel secondaire) et le ch. 3 hors §3.5 (UCS/BTS, couverts par nos propres
calibrations).
