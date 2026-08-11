# Branche article-exact — supprimer les différences structurelles avec Yan et al. 2023

*Directive Fernando (2026-08-11) : zéro différence structurelle entre rockim et
le code de l'article (MultiFracS). État des lieux, ce qui est fait, plan du reste.*

## Déjà identique (rien à faire)

Critère d'insertion (éq. 7-8) ; f(D) (éq. 11) et moteur mixte (éq. 12-16) ;
irréversibilité (Fukuda) ; calibration énergétique des ouvertures (éq. 13/15) ;
compression σ = p·o (éq. 9) ; décharge mode I sécante origine (éq. 17) ;
continuité de contrainte à l'insertion (§2.4) ; scission des nœuds (fig. 7 —
équivalence liaison/scission démontrée). En mode adaptatif, la naissance au pic
(dn0) rend la convention d'ouverture équivalente à leur extrinsèque pur (l'écart
résiduel = la largeur élastique ft/p ≈ nm).

## Fait sur cette branche (commit de ce jour)

- **éq. 6, terme 2µD** : clé `bulkViscosity` [Pa·s] — taux de déformation
  D = sym(Ḟ F⁻¹) tourné co-roté, ajouté à la contrainte avant rotation retour.
  Dissipatif par construction (2µ D:D ≥ 0). Défaut 0 = branche non exécutée,
  bit-identique (suite fast 7/7 revalidée). Valeur article : 7,6e3 (Table 1).
- **Preset `configs_yan/article_exact_base.cfg`** : toutes les formes littérales
  en un bloc — jointFrictionScaled = 1 (éq. 10), absorbFactor = 2 (éq. 21),
  bulkViscosity = 7,6e3, dampingLocal = 0 (leur amortissement est le visqueux,
  pas Cundall), jointXi = 0 (pas de dashpot de joint chez eux), adaptatif 4E/h,
  matériau Table 1 en commentaire.

## Chantier A2 — éq. 18, décharge en cisaillement (PETIT, prochaine séance)

Aujourd'hui : plasticité à retour (décharge sur la sécante de pénalité,
glissement plastique conservé) — coïncide avec l'éq. 18 en charge monotone.
À faire : `jointShearUnload = origin` — état s_max/τ_max par point, décharge et
recharge sur la sécante à l'origine, cap f(D)·fs inchangé. ~40 lignes + états.
Validation : charge nulle (0 joint), cas de cisaillement pur monotone
bit-identique à l'actuel, cycle charge-décharge-recharge contre tracé manuel
de l'éq. 18. Risque faible, isolé dans processJoint.

## Chantier B — contact par POTENTIEL de Munjiza + NBS (LE morceau, éq. 2-5)

Aujourd'hui : pénalité nœud-arête (grille de cellules, birth-gap, contact
débris quasi-plastique, frottement tanh). L'article : force normale distribuée
= p_n × intégrale de bord du gradient du potentiel sur la zone de RECOUVREMENT
des deux triangles (éq. 2-3, conservatif par construction), tangentiel
incrémental à ressort avec cap de Coulomb (éq. 4-5), détection NBS.

Plan proposé (une séance dédiée, branche article-exact) :
1. **Détection** : NBS de Munjiza 1998 (mappage triangles→cellules, listes
   chaînées par ligne) — remplace la grille actuelle pour les paires
   triangle-triangle candidates. La détection parallèle par listes de threads
   est conservée (elle est indépendante de la loi de force).
2. **Force normale** : clipping triangle-triangle (Sutherland-Hodgman, 2D),
   potentiel φ = min(3·A_i/A) par triangle, force = p_n ∮ n_Γ (φ_c − φ_t) dΓ
   sur le bord de l'intersection, appliquée aux nœuds par les poids
   barycentriques des segments de bord. ~400-600 lignes.
3. **Tangentiel** : ressort incrémental par PAIRE (éq. 4) avec cap de Coulomb
   (éq. 5) — remplace le tanh régularisé pour ces paires.
4. **Bascule** : `contact = penalty (défaut, inchangé) | potential` — tout
   l'existant reste bit-identique ; le potentiel est opt-in.
5. **Validation, dans l'ordre** : (a) charge nulle (0 joint, 0 force de
   contact parasite) ; (b) conservativité : deux blocs en collision élastique
   sans frottement → gcWork ≈ 0 à la précision machine (LE point fort du
   potentiel, notre pénalité ne l'a pas) ; (c) suite complète (les verify ne
   touchent pas le contact général : inchangés attendus) ; (d) A/B percussion
   maillage Yan penalty vs potential ; (e) SHPB (le contact EST le chemin de
   charge — comparaison fig. 23 refaite).
   Effort : 1-2 jours de travail concentré. C'est le seul endroit où la veille
   avait tranché « à ne pas reprendre » pour un bac à sable — la directive
   « zéro différence structurelle » le remet au programme, avec le runner pour
   filet.

## Non-différences (documentées, pas de chantier)

- **Cinématique d'élément** : leur éq. 6 (Green + split E/(1+ν), E/(1−ν),
  1/√|det F|) vs notre co-rotationnel + Biot + D-matrice : deux formulations
  objectives équivalentes aux petites déformations élastiques (les deux tournent
  sans déformation parasite) ; les modules apparents coïncident (E/(1−ν²)
  vérifié fig. 19b). Porter l'algèbre exacte n'apporterait aucun observable et
  casserait tous les bit-repères — à ne faire que si une comparaison terme à
  terme au point matériel devient nécessaire.
- **Intégration temporelle** : différences centrées explicites des deux côtés.
