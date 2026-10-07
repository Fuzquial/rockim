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

## Chantier A2 — éq. 18, décharge en cisaillement — **FAIT (2026-08-13)**

Clé `jointShearUnload = plastic | origin`, 2D **et** 3D, défaut inchangé.

- `plastic` (défaut) : plasticité à retour radial — décharge sur la sécante de
  **pénalité**, glissement plastique conservé.
- `origin` : l'éq. 18 — décharge **et** recharge sur la sécante à l'**origine**
  passant par (s_max, τ_env(s_max)), avec τ_env = min(p·s_max, cap). Symétrique
  exact de l'éq. 17 déjà en place pour le mode I.

Trois points de conception qui n'étaient pas dans le plan initial et qui se sont
imposés à l'écriture :

1. **Le glissement au pic est celui de Munjiza**, s_p = (c + tanφ·|σ_n|)/p,
   évalué sur l'enveloppe de Mohr-Coulomb **non endommagée** et sur la part
   *géométrique* p·dn de la contrainte normale (règle maison : un terme visqueux
   ne fixe jamais une résistance). C'est ce qui rend le calcul **non circulaire**
   — s_p ne dépend pas de D, alors que le cap, lui, en dépend — et c'est ce qui
   fait démarrer l'endommagement de mode II au **même instant** que le retour
   radial : les deux modes coïncident en charge monotone, comme annoncé.
2. **s_max est mis à jour AVANT la traction normale**, parce qu'en `origin` le
   moteur rs de l'éq. 16 en dépend, et que ce moteur alimente D, donc f(D), donc
   l'enveloppe normale du même pas.
3. **`J.slip` change de sens** en `origin` : ce n'est plus un glissement
   plastique mais l'**origine figée** de la sécante, celle qu'`activateJoint()`
   estampe à −τ₀/p à l'insertion adaptative. La continuité de contrainte en
   cisaillement au moment de l'insertion — le pendant de `dn0` en mode I — est
   donc conservée telle quelle, gratuitement.

⚠️ **Piège signalé au démarrage** : l'éq. 18 place *tout* le cap dans la
sécante, frottement de Coulomb compris. Avec `jointFrictionScaled = 0` (défaut
rockim) le glissement frottant devient **réversible** — plus d'hystérésis. La
forme littérale de l'article est donc `origin` **+** `jointFrictionScaled = 1`,
et c'est ce que pose désormais `configs_yan/article_exact_base.cfg`.

## Chantier A3 — contact par POTENTIEL de Munjiza (éq. 2-5) — **2D FAIT (2026-08-13)**

Clé `contact = penalty | potential`, défaut inchangé (bit-identique). Le cœur
géométrique est PUR et isolé dans `include/rockim/PotentialContact.hpp` :
clip de Sutherland-Hodgman du recouvrement triangle-triangle, potentiel
φ = 3·min(λ) (la « tente » de Munjiza), intégrale de bord **exacte** (subdivision
aux traversées des six médianes, trapèze exact par morceau), répartition nodale
consistante aux mêmes points spatiaux des deux côtés (3e loi machine).

- **Le test décisif** (`rockim selftest-potential2d`) : collision élastique
  sans frottement, frontale puis oblique — transfert exact (vA → 0 machine,
  vB → v0), **ΔKE/KE₀ = 3,7e-12**, quantité de mouvement 5e-15. Le compteur de
  travail de convention solveur porte un biais O(dt) documenté (~8e-4), la
  conservation se juge sur ΔKE.
- **Conservation au niveau solveur** : SHPB incassable sans frottement,
  gcWork = −2,6 J/m pour 766 J/m en jeu (0,3 % = biais du compteur), là où le
  penalty par défaut (quasi-plastique, gcRest = 0,2) dissiperait ~80 % de
  chaque rebond par construction.
- **SHPB réel** : onde incidente identique au penalty à 3e-6 près (pics
  −0,924477e-3 vs −0,924474e-3), transmise −0,219 vs −0,210e-3 ; le disque
  casse moins (595 vs 812) — la loi d'interface diffère, écart physique.
- **La leçon architecturale du chantier** : rockim cache le recouvrement des
  paires voisines derrière le joint vivant (sa pénalité porte la compression) ;
  quand un joint meurt COMPRIMÉ, la paire naît déjà en recouvrement et son
  énergie potentielle p∫φ surgirait du néant (+936 J/m mesurés sur la
  percussion 2D sans relève, +665 sans frottement). Une rampe temporelle ne
  suffit PAS (approche bradée / sortie plein tarif → +179 J/m injectés sur les
  premiers cycles) : la relève doit être un décalage d'ÉTAT — référence sur
  l'**aire** de recouvrement, décroissante en gcBirthTau, sortie sous la
  référence LIBRE (signe absorbant garanti, le pen0_ exact du potentiel).
  Résidu final +27 J/m, borné par un repère de suite. La solution de fond —
  l'architecture Munjiza où le contact porte la compression des paires
  voisines dès t = 0 et où le joint ne porte que la cohésion — est notée
  comme option de phase 3.
- **Physique** : la zone broyée devient CONSERVATIVE — percussion 2D : rebond
  e 0,55 → 0,71, 5 casses contre 174 à conditions égales. Le contact
  quasi-plastique du penalty était un canal de dissipation majeur des runs
  historiques ; le potentiel est la forme de l'article.
- **Phase 2 — FAIT (2026-08-13, même séance)** : portage 3D tet-tet. Clip du
  polyèdre de recouvrement par les 4 demi-espaces (face de coupe reconstruite
  par tri angulaire, valide par convexité), φ = 4·min(λ), intégration exacte
  par subdivision aux 12 plans de médiane, lumping nodal consistant. Selftest
  pointe-contre-face : transfert exact, ΔKE/KE₀ = 2,0e-8, quantité de
  mouvement machine (le tip-contre-tip de deux tets a été essayé et écarté :
  répulsion quasi nulle aux sommets — φ → 0 — les cônes broutent en position
  dégénérée). **Deux gardes nées du contrôle zeroload** : plancher de volume
  relatif et contrôle de fermeture du polyèdre — les tets exactement TANGENTS
  (voisins par arête/sommet d'un maillage qui pave l'espace) produisaient des
  slivers à volume quasi nul mais à grandes faces mal refermées, soit des kN
  parasites au repos (5 joints cassés à charge nulle, zéro après gardes, et
  le selftest frontal y a gagné un ordre de grandeur : 1,5e-7 → 2,0e-8).
  Percussion 3D Gmsh : mêmes 6 joints que la pénalité à T = 5e-5, gcWork
  −5e-5 J. Coût par paire supérieur au nœud-face (détection à optimiser —
  réutilisation des seaux) : combiner avec gcActivation = adaptive.

## Chantier B — détection NBS d'origine (note)

La détection actuelle (binning AABB en grille de hachage, paires uniques par
stamp) est O(N) comme le NBS de Munjiza & Andrews 1998 — le NBS « liste
chaînée » historique est une optimisation de constante, pas de complexité.
L'écart structurel restant avec l'article est donc porté par la phase 3D
ci-dessus, pas par la détection.

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
