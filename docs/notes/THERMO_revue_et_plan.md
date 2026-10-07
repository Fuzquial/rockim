# Choc thermique pour pré-fissurer avant percussion — revue et plan

**2026-09-01. F. Uzquiano, Mines Paris – PSL.**
Revue en six angles, **69 références**, chacune soumise à un vérificateur chargé
de la réfuter : **9 réfutées (13 %)**. Les références marquées « vérifiée »
ci-dessous ont vu leur notice ouverte ; celles dont seul le résumé a été lu le
disent.

---

## 0. Les trois choses à savoir avant tout le reste

### Saksala l'a déjà fait, et c'est notre auteur

- **T. Saksala (2020)**, *Thermal shock assisted percussive drilling: a numerical
  study on the single-bit axisymmetric case*, **IJRMMS**.
- **T. Saksala (2020)**, *3D numerical modelling of thermal shock assisted
  percussive drilling*, **Computers and Geotechnics 128, 103849**.
- **T. Saksala & A. Ibrahimbegovic (2020)**, *Thermal shock weakening of granite
  rock under dynamic loading: 3D numerical modeling based on embedded
  discontinuity finite elements*.

C'est **exactement** la question posée — choc thermique pour affaiblir avant
percussion — traitée numériquement par l'auteur dont la thèse porte déjà la loi
(`law = saksala`, `saksala2011` porté et vérifié à 8·10⁻¹⁴, `rapport_saksala.tex`,
`figs_saksala/`). Ce n'est donc pas une exploration : **c'est une cible de
reproduction**, avec un modèle de matériau que nous avons déjà.

### Mon verrou en ℓ_cz était MAL POSÉ — et il disparaît

J'avais annoncé que $\ell_{cz} = E G_{Ic}/f_t^2 = 14$–25 mm serait plus grand que
la couche fissurée, donc bloquant. **C'est faux.** Tarasovs & Ghassemi montrent
que la longueur qui gouverne un réseau de fissures de choc thermique n'est pas
$\ell_{cz}$ — évaluée à $f_t$ — mais

$$\xi = \left[\frac{K_{Ic}(1-\nu)}{E\alpha\,\Delta T}\right]^2 = \ell_{cz}\left(\frac{f_t}{\sigma_{th}}\right)^2$$

Comme $\sigma_{th}$ vaut 8 fois $f_t$ à $\Delta T = 100$ K, **$\xi$ est 60 à 240
fois plus petit** : 1,5 mm à 100 K, 0,38 mm à 200 K, 0,10 mm à 400 K. La zone de
process tient donc largement dans la couche fissurée dès 100 K. Le verrou
s'évanouit dès qu'on l'évalue **à la contrainte réellement appliquée**.

### Ma contrainte de paroi était une borne supérieure, pas une valeur

$E\alpha\Delta T/(1-\nu)$ suppose un contact thermique **parfait** (Dirichlet).
Cai et al. (2023) simulent **31,4 MPa** après 100 s d'injection d'azote liquide
là où la formule idéale donne 85 MPa ; Cha et al. (2017) mesurent que la paroi
n'atteint pas le point d'ébullition de l'azote avant plus de 20 minutes —
l'ébullition en film plafonne le flux. **Facteur de perte 2,5 à 3.**

Donc à $\Delta T = 100$ K : ~13–15 MPa réels, non 38. Toujours 2,6 à 3 fois
$f_t$, donc ça fissure — mais la marge est bien plus mince, et Lu & Fleck (1998)
donnent le cadre pour un nombre de Biot fini.

---

## 1. Les prédictions falsifiables, avec leurs préfacteurs

Tout repose sur **deux longueurs** : la longueur thermique transitoire
$L = \sqrt{4\kappa t}$ et la longueur matérielle $\xi$ ci-dessus.

| grandeur | loi | source |
|---|---|---|
| profil | $T(z) = T_0 - \Delta T\,\mathrm{erfc}(z/L)$ | Tarasovs & Ghassemi éq. 1 |
| contrainte | $\sigma_{th} = E\alpha(T_0-T)/(1-\nu)$ | éq. 2 — **confirme la mienne mot pour mot, $1/(1-\nu)$ compris** |
| **profondeur** | $a/\xi = (L/\xi)^{1{,}09}$, soit $a \approx 1{,}3$–$1{,}5\,L$ | éq. 7 |
| **espacement** | $d/\xi = 4{,}4\,(z/\xi)^{0{,}77}$ | éq. 8 |
| amorçage | $\Delta T_c = f_t(1-\nu)/(E\alpha) = 13{,}0$ K | Kingery 1955 (voir §4) |

**Correction n° 3 à mon estimation : les fissures ne s'arrêtent pas au front
thermique, elles le dépassent de ~40 %.** Ma règle $2z^*\sqrt{\alpha t}$
sous-estimait la profondeur.

Sur Red Bohus, $\Delta T = 100$ K, $t = 40$ s : couche fissurée **19–22 mm**
(quasi indépendante de $\Delta T$), espacement **8–11 mm** près de la surface et
**38–47 mm** au fond. Bahr et al. (2010) recoupent par $a\,d = 1{,}74\,L^2$ à un
facteur 2 près — **toute valeur dans 18–40 mm à 40 s est compatible avec les deux
approches indépendantes**.

### Le test le plus discriminant : deux échelles, pas une

Bourdin, Marigo, Maurini & Sicsic (PRL 2014) : l'espacement **à la naissance** du
motif vaut $\lambda^* \approx \sqrt{l_0\,l} \approx 2$ mm à 100 K, alors que le
motif **développé** est centimétrique (20–40 mm). **Le facteur ~10 entre les deux
EST le doublement de période** de Bažant & Ohtsubo (1979) : une fissure sur deux
s'arrête, l'espacement passe de $s$ à $2s$ à $4s$.

Un calcul correct doit montrer **les deux échelles successivement**. Signature
mesurable : le nombre de fissures actives chute par facteurs 2, et l'histogramme
des longueurs est **multimodal**, pas unimodal.

---

## 2. Le vrai obstacle n'est pas la mécanique, c'est le bilan d'énergie

C'est le résultat le plus important de la revue, et il est sobre.

**Le gain mécanique existe et il est grand, mesuré :**

- **Satish et al. (2006)** — basalte micro-ondé puis **foré en percussion** :
  **+42 % de vitesse de pénétration**. C'est notre expérience exacte.
- **Liu et al. (2025)** — jet d'**azote liquide** (notre mécanisme, le froid) plus
  indenteur sur granite chaud : **+140 % de volume broyé, −22 % de force normale
  crête**.
- **Rossi, Saar & Rudolf von Rohr (2020)**, accès libre, lu intégralement —
  granite du Grimsel : énergie spécifique de coupe **54,64 → 2,62 J/mm³ (÷21)**,
  taux d'enlèvement **×30**.

**Mais le bilan n'est presque jamais bouclé.** Rossi et al. écrivent noir sur
blanc que leur calcul « ne tient pas compte de l'énergie thermique nécessaire ».
Quand on le boucle :

| source | dépensé | gagné |
|---|---|---|
| Satish 2006 | 100 kWh/t | +42 % pénétration |
| Hassani 2016 | jusqu'à 740 kWh/t | −30 % UCS |
| Hartlieb & Grafe 2017 | — | 4,7 kWh/t d'énergie de coupe économisée |

**On dépense 20 à 150 fois ce qu'on économise.** La seule échappatoire publiée est
Kingman et al. (2004) : monomode, **très forte densité de puissance, impulsion
très courte** — 0,4 à 0,83 kWh/t suffisent alors pour −30 à −50 % de résistance.

> **La leçon transposable, et elle doit piloter tout le reste : ce qui décide du
> bilan n'est pas la DOSE totale mais la DENSITÉ DE PUISSANCE et la BRIÈVETÉ du
> dépôt.** Un refroidissement lent et diffusif est presque sûrement perdant. Si
> l'idée doit marcher, c'est par un choc bref et intense.

---

## 3. Comment la littérature résout la séparation des échelles

**Ce n'est pas du sous-cyclage : c'est un schéma étagé.** La documentation Itasca
(PFC 7.0 et 3DEC 9.0, ouverte et lue) donne la même formule dans les deux codes :

$$\frac{\Delta t_{th}}{\Delta t_{mec}} = \frac{L_c}{\kappa}\sqrt{\frac{K + 4G/3}{\rho}}$$

Sur nos chiffres : **2,7·10⁶** pour une maille de 1 mm, **1,1·10⁷** pour 4 mm —
notre facteur 10⁷ est confirmé et localisé.

La stratégie : on avance d'**un pas thermique**, puis on prend des **sous-pas
mécaniques dont le temps n'est pas compté** — la mécanique n'est plus une
dynamique, c'est une **relaxation quasi-statique**. Chiffrage pour notre cas
(maille 0,3 mm, $\Delta t_{th} = 0{,}06$ s, $\Delta t_{mec} = 1{,}2\cdot10^{-8}$ s) :

- explicite pur : 40 s coûtent **3,3·10⁹ pas** — infaisable ;
- étagé : ~670 pas thermiques × ~2000 sous-pas = **1,3·10⁶ pas**, **gain ×2500**.

Convergence indépendante : Bourdin *et al.* traitent le même problème de trempe
**en quasi-statique**. Le quasi-statique n'est donc pas un pis-aller numérique,
**c'est la physique du problème**.

**La contrainte thermique entre de deux façons, toutes deux additives :**

- **Lignée Yan** (Han et al. 2022, texte intégral lu) : force nodale équivalente,
  $\mathrm{d}\sigma = -\delta_{ij}\,3K^*\alpha\,\mathrm{d}T$ (éq. 5) et
  $f_n = -\tfrac12 3K^*\alpha\,\mathrm{d}T\,n_j L$ (éq. 6), injectée comme charge
  de volume sur le triangle. **Elle ne touche pas à la loi de joint** : elle
  charge le triangle, qui charge le joint.
- **Lignée Latham-Xiang** (Joulin et al. 2020), celle de rockim : décomposition
  multiplicative $F = F_e F_\theta$ avec $F_\theta = (1+\alpha\,\mathrm{d}T)I$.

Les deux coïncident en petites déformations à $(\alpha\Delta T)^2 = 6\cdot10^{-7}$
près à 100 K — **un test de non-régression gratuit.**

Et **Yan & Jiao (2020)** donnent la rétroaction inverse : transfert de chaleur
discret tenant compte de la **résistance thermique des fractures**. C'est le
couplage à deux sens, si on le veut un jour.

---

## 4. Les neuf réfutations, et les deux qui changent quelque chose

13 % des références proposées n'ont pas survécu au contre-examen. Deux méritent
d'être retenues parce qu'elles auraient été citées à tort :

- **Le paramètre $R = \sigma_f(1-\nu)/(E\alpha)$ n'est PAS de Hasselman 1969.**
  Il est de **W. D. Kingery (1955)**, *Factors Affecting Thermal Stress Resistance
  of Ceramic Materials*, J. Am. Ceram. Soc. 38(1), 3-15,
  doi:10.1111/j.1151-2916.1955.tb14545.x. Hasselman 1969 apporte autre chose : la
  théorie unifiée amorçage/propagation.
- **Le banc de validation céramique n'est pas Jiang et al. 2012** (qui existe et
  est réel, mais ne porte pas la figure voulue) : c'est **Shao et al. (2011)**,
  J. Am. Ceram. Soc. 94(9), 2804-2807, doi:10.1111/j.1551-2916.2011.04728.x,
  dont la fig. 5(d) est celle que Bourdin et al. reproduisent.

---

## 5. Plan d'implémentation pour rockim

Issu de l'audit de code sur `rockim_f2`, avec les points de greffe vérifiés.

### Étape 1 — la pré-contrainte thermique (petit, ~1 journée)

| élément | où |
|---|---|
| stockage | `std::vector<double> Tel_;` à côté de `hEl_` (`FdemSolver.hpp:503`) + `k3aP_` de $3K\alpha_T$ par phase, rempli dans la boucle `FdemSolver.cpp:668-676` (`Material::K()` existe déjà). **1,52 Mo** armé pour 190 000 éléments, **0 octet** désarmé |
| setup | une `setupThermal()` juste après `setupExcavation()` (`:1072`) |
| clés | `thermTemp`, `thermTref`, `thermAlpha`, `thermStart`, `thermRamp` — **surtout pas `T`**, déjà pris pour la durée du run |
| assemblage | **une ligne** à `:3996`, juste après le terme in situ : `if (thermOn_) P += sTh * R;` avec une rampe calquée sur `excavRelief()` (`:6472-6479`) |
| sortie | carte `vtk::ScalarField` nommée dans `writeFrame` (`:6856-6876`), colonnes d'historique **en queue** |

> **Piège identifié par l'audit, et il est sérieux :** ne **pas** fondre le
> thermique dans `insituS_`. `excavationForces()` le relit (`:6492`) et
> relâcherait la part thermique de la traction de paroi — ce qui changerait **en
> silence** tous les runs tunnel.

### Étape 2 — le phasage vers la percussion (petit)

**`toolStart` n'existe pas** : `tool_.v` est posé à $t=0$ (`:3000`) et l'outil est
intégré à chaque pas. Il faut l'ajouter, avec une garde d'ordre sur le modèle
`toolStop`/`fragBrushStart` (`:468-474`). Et dimensionner `toolGap` (défaut
10⁻⁴ m) contre le **soulèvement thermique** $\alpha_T\Delta T H \approx 0{,}16$ mm
pour 100 K sur 0,2 m — sinon contact prématuré.

### Étape 3 — le schéma étagé (moyen, et c'est le vrai chantier)

Sans lui, une rampe thermique de durée $t_{th}$ coûte $t_{th}/\Delta t$ pas **au
pas de la percussion** : aucun mass scaling n'existe et `dt_` est figé (`:1139`).
C'est l'étape qui décide si le sujet est finançable.

### Le raccourci, disponible aujourd'hui

**`preBrokenJoints`**, ajouté ce matin, permet de poser une couche pré-fissurée
**sans calculer le choc thermique** : une famille de segments à l'espacement
donné par la loi de Tarasovs & Ghassemi ($d = 8$–11 mm en surface, 38–47 mm au
fond, profondeur 19–22 mm). On mesure alors directement **le gain en percussion**
— la seule grandeur qui décide de l'intérêt du procédé — sans écrire une ligne
de thermique.

**C'est de loin le meilleur rapport valeur/effort du dossier**, et je le
recommande en premier.

### Les maillages existent

`make_cut_mesh.py W H depth notchLen hFine bandH hFar out.msh [seed]` avec
`notchLen ≤ 0` donne un rectangle avec une **bande fine sous la face
supérieure** — exactement la géométrie d'une couche fissurée sous un fond de
trou. `make_circle_mesh.py` couvre le forage. Aucun générateur 3D de forage
n'existe.

---

## 6. Ce qui reste à trancher

1. **Pré-contrainte simple ou thermo-élasticité complète ?** La première est une
   ligne, mais elle est aveugle au `crushCap`, au `bulkDamage` et à la loi de
   volume. La seconde touche trois branches constitutives.
2. **Reproduire Saksala 2020 d'abord ?** Son cas axisymétrique à un insert est la
   cible naturelle, et nous avons sa loi.
3. **Le bilan d'énergie fait-il partie du livrable ?** Si oui, il faut instrumenter
   l'énergie thermique déposée, pas seulement l'énergie mécanique économisée —
   sinon on publiera un ÷21 qui ne veut rien dire, comme la littérature le fait
   régulièrement.
