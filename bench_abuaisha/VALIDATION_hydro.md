# Statut de validation du couplage hydro-mécanique de rockim

Benchmark AbuAisha, Eaton, Priest et Wong, *Hydro-mechanically coupled FDEM
framework to investigate near-wellbore hydraulic fracturing in homogeneous and
fractured rock formations*, *J. Petrol. Sci. Eng.* **154** (2017) 100–113,
section 3 et annexe A.

État au 2026-08-21, 04 h 10. Source : `simulations/FDEM/rockim/rockim_p1/`.

---

## 1. Objet et portée du document

Ce document sépare trois catégories que la lecture des figures seules confond
facilement : ce qui est vérifié par une mesure, ce qui relève d'une hypothèse
assumée du modèle, et ce qui reste ouvert. Il répond à une question précise,
posée le 2026-08-21 : les équations sont-elles correctement implémentées et la
physique est-elle correcte ?

La réponse courte est que le signe et la mécanique du couplage sont vérifiés par
plusieurs contrôles indépendants, mais que le modèle de fluide comporte des
simplifications majeures, que le critère de mouillage n'a aucun fondement
physique, et que la suite de non-régression ne couvre pas le module.

---

## 2. Le correctif de signe du 2026-08-20

`hydroForces()` appliquait la pression suivant la normale sortante du solide au
lieu de son opposé. Le fluide serrait la cavité au lieu de l'ouvrir, et le
forage rompait à 6,6 MPa en produisant un breakout aligné sur sigma'_h, au lieu
des 12 MPa attendus en traction.

La forme corrigée, `src/FdemSolver.cpp:5017` :

```cpp
Eigen::Vector2d half = -0.5 * hydroP_ * L * thk_ * n;  // vers l'INTERIEUR du solide
```

La racine du bug est une mauvaise lecture de leur équation 7, qui s'écrit
`F = -(p/2)[y2-y1 ; x2-x1]`. Ce vecteur n'est pas orthogonal au segment, son
produit scalaire avec `(dx, dy)` valant `2 dx dy`. La coquille est donc dans la
seconde composante seule, et non dans le signe de tête. C'est pourtant ce signe
de tête, le bon, qui avait été supprimé en croyant corriger la coquille.

Deux éléments de l'article tranchent. Leur section 3.1 pose « negative sign to
compressive stresses », donc une pression p > 0 impose sigma = -p I et la
traction t = sigma . n = -p n. Et Y-Geo triangule en sens direct (Munjiza 2004,
Mahabadi 2012, Lisjak 2014a), donc leur normale est bien sortante, comme celle
de rockim.

Le contrôle qui existait alors, `parker_compare.py`, mesurait `max(y) - min(y)`,
une valeur absolue. Il a validé une interpénétration comme une ouverture. La
leçon est générale et vaut d'être retenue : un contrôle de signe qui passe par
une norme ne contrôle pas le signe.

Retour arrière disponible dans `_backup_avant_correctif_signe_2026-08-20/`.

---

## 3. Contrôles quantitatifs passés

| contrôle | référence | mesure | écart |
|---|---|---|---|
| croisement de volume, cas anisotrope | contrainte moyenne 5,70 MPa | 5,73 | +0,5 % |
| croisement de volume, cas isotrope | 4,60 MPa | 4,63 | +0,6 % |
| croisement de volume, essai 1 | 3,50 MPa | 3,52 | +0,7 % |
| déplacement de paroi, deux chemins de chargement | Lamé, 21,771 µm | 21,336 | −2,0 % |
| écart entre les deux chemins de chargement | 0 | 0,000e+00 | exact |
| ouverture de fissure mûre, cas anisotrope | Sneddon, 92,7 µm | 89,0 | −4 % |
| volume de cavité, somme w·L contre théorème de Green | — | 5,34 contre 5,65e−5 m²/m | 6 % |
| bilan énergétique | — | — | 0,0125 % |
| suite de non-régression, tier fast | valeurs de référence Linux | identiques | 19/19 |

Le croisement de volume est le plus probant des trois premiers. La cavité doit
repasser au-dessus de son volume initial exactement quand la pression de fluide
égale la contrainte in situ moyenne : en dessous le confinement l'emporte et le
trou se referme, au-dessus le fluide l'emporte et le trou s'ouvre. Trois
confinements différents donnent trois accords à moins de 0,7 %, avec un biais
systématiquement positif et minuscule, le fluide devant dépasser légèrement le
confinement pour le vaincre. Un signe inversé rendrait ce croisement impossible :
le volume ne remonterait jamais.

Le contrôle de signe dédié est `tools/hydro_sign_check.py`, sur les configs
`signe_conf.cfg` et `signe_hydro.cfg`, dix secondes chacune.

---

## 4. Hypothèses du modèle, et ce qu'elles excluent

### 4.1 Absence d'écoulement

La pression est un scalaire unique pour toute la cavité mouillée, recalculé à
chaque pas depuis la masse injectée et le volume total (`FdemSolver.cpp:4982`) :

```cpp
hydroP_ = (r > 1e-300) ? hydroP0_ + fluidK_ * std::log(r) : hydroP0_;
```

Il n'existe donc ni viscosité, ni gradient de pression, ni délai de propagation.
La pression en pointe d'une fissure de 400 mm égale instantanément celle du
forage. C'est le modèle de leur annexe A, ce qui rend la reproduction fidèle au
benchmark, mais ce n'est pas la physique de la fracturation hydraulique.

### 4.2 Absence de fuite dans la matrice

La recherche des termes `Darcy`, `leak`, `poro` et `permeab` dans le solveur ne
rend aucune occurrence. Il n'y a ni fuite, ni poroélasticité, ni couplage avec
une pression de pore. Les contraintes manipulées sont effectives par convention,
la pression de formation de 28 MPa étant soustraite en amont dans les decks.

### 4.3 Critère de mouillage

C'est le point le plus faible du modèle. La clé `hydroWetDamage`
(`FdemSolver.cpp:4671`) fixe le seuil d'endommagement à partir duquel une
interface conduit le fluide :

```cpp
wetDmin_ = cfg_.getd("hydroWetDamage", 1.0);
```

Son défaut, 1.0, reproduit le comportement historique où seule une interface
entièrement rompue conduit. La valeur 0.0 fait conduire toute interface insérée
dès sa naissance. Aucune de ces deux valeurs n'est physique : la première laisse
une zone cohésive en train de céder rester étanche, la seconde fait pénétrer le
fluide dans une interface intacte qui porte encore toute sa cohésion. Ce sont
deux bornes, et rien dans le code ni dans la littérature consultée ne dit où se
situe la réalité entre elles.

L'enjeu est chiffré : sur le cas isotrope à sigma' = 3,5 MPa, le pic vaut
13,750 MPa à la borne haute et 12,806 à la borne basse. L'écart de 0,94 MPa
représente la moitié du dépassement au seuil analytique.

### 4.4 Module de compressibilité du fluide

`fluidBulk = 2.2e9` est une hypothèse assumée. L'article ne donne pas cette
valeur. Retrouver leur instant de rupture, t = 1,08 ms d'après la légende de
leur figure 9, en demanderait environ 4,4 GPa. Ce paramètre ne change pas la
pression de rupture, que la mécanique fixe, mais il fixe toute l'horloge du
calcul.

### 4.5 Frontières

Les calculs tournent avec `absorbing = none` et des appuis bloqués, conformément
à leur figure 6. Le rayonnement de la rupture reste donc piégé dans le domaine,
ce qui se voit sous forme de cercles concentriques dans le champ de contrainte
principale. Sans effet sur le pic, qui est atteint avant, mais à prendre en
compte pour la phase de propagation.

---

## 5. Points ouverts

Le pic dépasse encore le seuil analytique de 6,7 à 17,9 % selon le cas. Le
budget qui décompose ce dépassement a été vérifié deux fois par prédiction avant
lancement, à 0,3 % et 1,6 %, mais sa troisième prédiction a manqué de 6,0 % et
l'une de ses hypothèses est tombée : l'incubation n'explique pas la moitié de
l'écart entre états de contrainte, puisque supprimer l'incubation laisse cet
écart intact. Le détail est en section 7.2.

Deux hypothèses restent non mesurées. La contribution du nombre de sites
d'amorçage s'appuie sur une unique réalisation de maillage, et seul un essai à
plusieurs graines trancherait. La moitié statique de l'écart, désormais chiffrée
à 0,677 MPa, n'est attribuée à la profondeur d'échantillonnage que par
raisonnement, sans essai dédié.

Les vitesses de particule atteignent 1,72 m/s contre une échelle plafonnée à
0,12 m/s dans leurs figures 12 et 13, soit quatorze fois plus, ce que 25 % de
pression supplémentaire n'explique pas. Le suspect est leur amortissement
visqueux, dont le coefficient mu = 5,6e5 kg/(m·s) vaut exactement 2h racine(E rho)
pour h = 30 mm, soit dix fois le critique de Munjiza pour leurs éléments de 3 mm.
Coquille d'exposant dans leur Table 1 ou sur-amortissement délibéré, la question
n'est pas tranchée.

Enfin, la valeur de référence elle-même est incertaine. Leur texte annonce
environ 12,5 MPa quand leur propre figure 11b se lit à 11,69 MPa.

---

## 6. Couverture de la suite de vérification

`tools/verify_suite.py` compte une soixantaine de tests distincts, tous tiers
confondus. Le mot `hydro` n'y apparaît pas une seule fois. Le module dont le
signe était inversé est donc précisément le seul que la suite ne teste pas, et
il ne le teste toujours pas.

Le contrôle `tools/hydro_sign_check.py` existe et passe, mais il n'est câblé
dans aucun tier. Rien n'empêche le bug de revenir au prochain remaniement.

C'est la lacune la plus concrète et la plus vite comblée de ce document : les
deux configs de contrôle tournent en dix secondes chacune.

---

## 7. Campagne d'essais

### 7.1 Terminés

| essai | question | résultat |
|---|---|---|
| runs de référence | reproduire leur figure 11 | anisotrope 14,993 MPa contre cible 12,00 (+24,9 %) ; isotrope 16,078 contre 14,20 (+13,2 %) |
| essai 1 | neutraliser l'artefact de normalisation en portant la cible isotrope à 12,00 | pic 13,750 contre 13,79 prédits, soit 0,3 % ; l'écart entre les deux états passe de 11,7 à 10,3 points, donc 88 % en est physique |
| essai 2, isotrope | supprimer l'incubation à zone cohésive sèche | pic 12,806 contre 12,61 prédits, soit 1,6 % |

L'essai 1 a permis, grâce aux colonnes `nInserted` et `nDamaging` ajoutées à
`history.csv`, de séparer pour la première fois deux contributions que le seul
compte de joints rompus mélangeait. De la cible analytique à la première
insertion, la pression gagne 0,613 MPa : c'est le dépassement statique, l'écart
entre la contrainte que l'équation 10 suppose en paroi et celle que le champ
éléments finis y produit. De la première insertion au pic, elle gagne encore
1,137 MPa : c'est l'incubation.

L'essai 2 attaque cette seconde part. Son résultat sur le cas isotrope :

| | essai 1 | essai 2 |
|---|---|---|
| première insertion | 12,613 MPa | 12,635 |
| première rupture | 13,671 | 12,685 |
| pic | 13,750 | 12,806 |
| incubation, pic moins insertion | 1,137 | 0,171 |
| dépassement sur la cible de 12,00 | +14,6 % | +6,7 % |

La pression d'insertion ne bouge pas, 0,17 % d'écart. Le patch ne touche donc
pas au champ statique, seulement à ce qui suit la naissance de l'interface : la
part statique du budget reste intacte et seule l'incubation s'effondre, de 85 %.
Le résidu de 0,171 MPa est le temps qu'il faut au fluide pour pénétrer et à la
fissure pour s'ouvrir.

Deux effets secondaires méritent d'être notés. L'écart entre première insertion
et première rupture passe de 1,058 MPa à 0,051 : avec le fluide présent dans la
zone cohésive dès sa naissance, une interface qui s'endommage rompt presque
immédiatement. Et l'invariant des sept joints rompus au pic, observé dans les
trois runs à incubation sèche, tombe à 22 : il était propre à ce régime.

Le correctif du timbre de cache a été nécessaire pour que le patch morde. Le
front mouillé n'était recalculé que si `nBroken_` avait changé
(`FdemSolver.cpp:4774`), ce qui suffit tant que seule une interface rompue
conduit, puisque la topologie mouillée ne peut alors changer qu'à une rupture.
Dès que le seuil descend, l'équivalence tombe et le front serait resté figé entre
deux ruptures. Le timbre porte désormais sur le nombre d'interfaces conductrices.
La signature le confirme : en fin de run, 3571 faces mouillées pour 1864
interfaces insérées et 729 rompues, soit un rapport de 1,86 sur les insérées une
fois retirées les 105 faces initiales du forage.

### 7.2 L'essai 2 anisotrope, et la falsification d'une hypothèse du budget

Terminé le 2026-08-21 à 04 h 03, en 6398 s. Pic mesuré 14,145 MPa contre
13,34 prédits dans le bandeau du deck, soit **6,0 % d'erreur**. C'est la première
prédiction manquée de la campagne, après 0,3 % sur l'essai 1 et 1,6 % sur l'essai
2 isotrope.

La cause de l'erreur est identifiée. La prédiction supposait que l'incubation
anisotrope, budgétée à 1,653 MPa, s'effondrerait comme celle du cas isotrope,
tombée de 1,137 à 0,171. Elle n'est descendue qu'à 0,833.

Conséquence plus lourde que l'erreur elle-même : **l'écart entre les deux états de
contrainte n'a pas diminué.** Il passe de 1,243 MPa, soit 10,36 points de la
cible, à 1,339 MPa, soit 11,15 points. L'hypothèse selon laquelle l'incubation en
expliquait la moitié est donc fausse. L'incubation contribue au dépassement du
seuil analytique dans les deux cas, mais de façon assez semblable pour
s'annuler dans la différence.

En revanche, l'essai livre une mesure neuve. La pression d'insertion est pour la
première fois instrumentée dans les deux états de contrainte, ce qui permet de
décomposer l'écart :

| | anisotrope | isotrope | écart |
|---|---|---|---|
| insertion | 13,312 MPa | 12,635 | +0,677 |
| incubation, pic moins insertion | 0,833 | 0,171 | +0,662 |
| pic | 14,145 | 12,806 | +1,339 |

L'écart se partage donc presque exactement en deux moitiés : **51 % de champ
statique**, c'est-à-dire la contrainte que le champ éléments finis produit en
paroi avant toute fissuration, et **49 % d'incubation résiduelle**, celle qui
résiste au mouillage direct dans le cas anisotrope.

Une hypothèse, non établie, pour cette résistance : la concentration de
contrainte anisotrope étant plus localisée, moins de sites s'ouvrent
simultanément et le fluide dispose de moins de chemin connexe à exploiter tôt.
Un essai à plusieurs graines de maillage la testerait.

Conséquence pour la suite : l'essai de raffinement de paroi devient prioritaire,
puisqu'il attaque directement la moitié statique, désormais mesurée à 0,677 MPa
et non plus inférée.

| | essai 1 | essai 2 |
|---|---|---|
| pic anisotrope | 14,993 MPa (+24,9 %) | 14,145 (+17,9 %) |
| pic isotrope | 13,750 (+14,6 %) | 12,806 (+6,7 %) |
| écart entre les deux | 10,36 points | 11,15 points |

### 7.3 En cours

Aucun calcul en cours.

### 7.4 Proposés, non lancés

| essai | ce qu'il tranche | clé | coût |
|---|---|---|---|
| raffinement de paroi | la moitié statique de l'écart entre états de contrainte, mesurée à 0,677 MPa, et la profondeur d'échantillonnage. **Devenu prioritaire par le résultat de l'essai 2 anisotrope** | `hFine` du générateur de maillage | élevé, le pas de temps tombe avec la maille |
| frontières absorbantes | les ondes réfléchies piégées | `absorbing = sides` | environ 2 h 15 |
| graines de maillage | le nombre de sites d'amorçage, seule hypothèse non mesurée du budget | argument de graine de `make_circle_mesh.py`, actuellement 1 | 3 graines sur 2 cas, environ 13 h 30 |
| amortissement | l'écart de quatorze sur les vitesses | `dampingLocal`, actuellement 0,15 | 3 valeurs, environ 7 h |
| insertion intrinsèque | le schéma d'insertion, structurellement le leur | `insertion = intrinsic` | très élevé, 284 124 joints actifs dès t = 0 |

---

## 8. Provenance et emplacements

| élément | emplacement |
|---|---|
| solveur | `src/FdemSolver.cpp`, en-tête `include/rockim/FdemSolver.hpp` |
| decks | `bench_abuaisha/configs/` |
| outils de dépouillement | `bench_abuaisha/tools/` |
| figures, rangées par campagne | `bench_abuaisha/figs/` sur le drive |
| historiques des runs | `bench_abuaisha/historiques/` sur le drive |
| journal de non-régression | `bench_abuaisha/verify_e2_2026-08-20.log` |
| retours arrière | `_backup_avant_correctif_signe_2026-08-20/`, `_backup_avant_essai2_2026-08-20/` |

Binaires : `rockim_e1.exe` pour l'essai 1, `rockim_e2.exe` pour l'essai 2. Le
second ajoute la clé `hydroWetDamage` dont le défaut rend le chemin
bit-identique au premier, ce que la suite confirme en 19/19 aux valeurs de
référence inchangées.
