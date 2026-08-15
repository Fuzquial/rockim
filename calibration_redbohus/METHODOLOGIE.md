# Méthodologie de calibration Red Bohus — proposition à valider

*Rédigé le 2026-08-15 après lecture de Ye et al. 2025 (IJRMMS 194, 106233),
Bu et al. 2026 (IJRMMS 199, 106400) et Jiang et al. 2025 (Sci. Rep. 15:34923).
Complète le README (cibles et phases) en fixant le PROTOCOLE exact.*

---

## 1. Objectif et critères de succès

Obtenir un jeu de paramètres FDEM **avec son incertitude** qui reproduise le
comportement macroscopique du Red Bohus, et **prédise** des états non vus.

| Critère | Seuil de réussite | Justification |
|---|---|---|
| UCS | ±15 % | la dispersion expérimentale est de 17 % (4 essais : 112 à 158 MPa) — viser mieux n'a pas de sens |
| BTS | ±20 % | Jiang mesure un COV > 8 % sur le BTS numérique ; Bu plafonne à R² = 0,73 ; Ye renonce à l'ajuster |
| q(σ₃ = 20 et 50) | ±10 % | essais très reproductibles (±0,5 %), ce sont les cibles exigeantes |
| E macro | ±15 % | trois valeurs expérimentales concurrentes (52 / 57,3 / 77,7 GPa) |
| Faciès | qualitatif | brésilien diamétral, UCS en cisaillement/fissures multiples, triaxial en bande unique |
| **Prédiction** σ₃ = 75 et 100 | **±20 %, non vus au calage** | le test que les trois articles ne font pas |

---

## 2. Les cinq leçons de la littérature, et ce qu'on en fait

1. **Le conflit traction/compression est STRUCTUREL, pas une maladresse.**
   Ye privilégie l'UCS et laisse filer le BTS ; Bu le démontre (BTS R² 0,74 →
   0,90 dès qu'on retire c et φ des cibles) ; Jiang obtient son pire R² (0,65)
   sur la résistance en traction. → **On l'anticipe** : optimisation
   multi-objectif (front de Pareto, pas de somme pondérée), et le compromis
   retenu est un CHOIX documenté, pas un échec.
2. **Réduire la dimension avant d'apprendre** (Jiang : Pearson, exclut
   ν_bloc et k_n/k_s dont |r| < 0,2). → **Phase 1 = criblage + corrélations**,
   les paramètres non influents sont FIGÉS et déclarés.
3. **Les essais confinés sont indispensables** pour c et φ macro (Bu en fait
   1 944 ; Ye, qui s'en passe, n'a pas d'enveloppe). → σ₃ = 20 et 50 MPa
   entrent dans la calibration.
4. **Bon émulateur direct ≠ bon inverseur** (Bu : la forêt aléatoire gagne en
   R² mais échoue en inversion, ils basculent sur un réseau profond). → on
   compare les émulateurs **sur la tâche inverse**, pas sur le R² direct.
5. **Le BTS numérique est très dispersé** (Jiang : COV 8-20 %). → **3 graines
   de maillage** pour tout point brésilien, la cible est la moyenne.

---

## 3. Phase A — sélection de l'architecture (le préalable)

rockim offre trois architectures crédibles, une par école. Aucune n'est
évidente a priori : **on tranche par la mesure**, sur un protocole identique.

| # | Architecture | Bulk | Paramètres à calibrer | École |
|---|---|---|---|---|
| **A1** | GBM Voronoï, bulk élastique | `elastic` + crushCap | 6 : ft, c, φ des joints, Gf, gfShearFactor, crushCap | Bu / Jiang |
| **A2** | bulk élasto-plastique Mohr-Coulomb | **`law = mc`** (implémenté et vérifié 1e-9 %) | 9 : les 5 précédents + mcCohesion, mcFrictionDeg, mcDilationDeg | **Ye** |
| **A3** | bulk DP-DFH **déjà calibré** sur Red Bohus | `law = dpdfh` (carte VUMAT validée 4,7e-12) | **5 : les joints seuls** | thèse |

**Protocole A** : chaque architecture est lancée sur les 4 essais (UCS, BTS,
triaxial 20 et 50), à son jeu de départ, 3 graines pour le BTS → **36 runs**,
~1 h. Critères de sélection : (a) les quatre essais cassent-ils vraiment
(le taux de rejet de 72 % de Ye vient de là) ; (b) l'enveloppe simulée est-elle
**concave** ; (c) le ratio UCS/BTS est-il dans le bon ordre de grandeur (cible
12,3) ; (d) coût par run.

> **Attendu** : A3 est le mieux posé (bulk figé par les campagnes Abaqus,
> seuls les joints restent inconnus) et c'est la seule loi qui sépare
> vraiment traction et compression — ce que Ye appelle explicitement de ses
> vœux en conclusion. Mais la mesure décide.

---

## 4. Phase B — criblage et corrélations (méthode Jiang)

Sur l'architecture retenue : plan **mono-variable**, chaque paramètre pris à
3 valeurs autour du centre, les autres au centre → ~20 jeux × 4 essais =
**80 runs**, ~1 h 30.

Sorties enregistrées : E, UCS, BTS, q(20), q(50), plus les **descripteurs de
faciès** (fraction de cisaillement, diamétralité du brésilien, angle de la
bande, nombre de fissures). Puis **matrice de corrélation de Pearson**
paramètres × sorties.

**Livrable** : la liste des paramètres retenus (|r| > 0,3 sur au moins une
cible) et des paramètres figés, avec leur valeur et sa justification.

---

## 5. Phase C — base de données (plan d'expériences)

**Hypercube latin** sur les paramètres retenus (4 à 6 selon la phase B),
bornes élargies à la Bu/Jiang, **200 jeux × 4 essais = 800 runs**. Coût
unitaire 40-120 s selon la loi → **6 à 12 h en série, 2 à 3 h en parallèle
(6 runs × 3 threads)**. Une nuit suffit largement.

Garde-fous : `stopPeakDrop = 0.3` (on n'achète pas le post-pic profond),
`budgetAbortPct = 5` + `budgetAbortMin = 0.05` (moniteur E2 armé), 3 graines
pour le BTS. Les jeux dont un essai échoue sont **conservés et étiquetés**
(ils informent l'émulateur sur les frontières du domaine valide — c'est
l'information que le pré-classement de Ye jette).

---

## 6. Phase D — émulateur et inversion

1. **Émulateurs comparés** : processus gaussien (krigeage) et forêt aléatoire,
   split 75/25, graine fixée — **et le juge est la tâche INVERSE**, pas le R²
   direct (leçon Bu) : on retire 10 jeux de la base, on inverse leurs sorties,
   on mesure l'erreur sur les paramètres retrouvés ;
2. **Inversion en deux passes** (protocole Bu) : grille grossière sur tout le
   domaine, puis fine autour du meilleur ;
3. **Front de Pareto multi-objectif** (protocole Ye) sur trois objectifs :
   F1 = compression (UCS + enveloppe), F2 = traction (BTS), F3 = faciès —
   chacun en **erreur relative sur la courbe entière** quand la courbe
   expérimentale existe (UCS et triaxiaux), sur le pic sinon ;
4. **Postérieur bayésien** sur l'émulateur GP, vraisemblance pondérée par les
   **écarts-types expérimentaux réels** (UCS ±21,4 ; q ±2,8 et ±2,6 ; BTS
   ±0,98), avec un terme de discrépance de modèle. → intervalles de crédibilité
   et **corrélations entre paramètres** : c'est notre apport, aucun des trois
   articles ne le fournit.

---

## 7. Phase E — validation

- **3 graines** au jeu retenu, sur les 4 essais de calibration → vérification
  des seuils du §1 ;
- **Prédiction pure** : σ₃ = 75 et 100 MPa, jamais vus pendant le calage ;
- **Faciès** comparés aux modes de Jiang et aux photos de Dumoulin ;
- **Contrôle 3D** : UCS + triaxial 20 sur maillage grossier — l'étude
  d'objectivité du 14/08 (pics 75,0-76,9 MPa sur un facteur 6 en éléments)
  autorise le maillage grossier ;
- **Bilan d'énergie** < 1 % sur tous les runs retenus (il ferme à 1e-12 %
  depuis les corrections du 14-15/08).

---

## 8. Budget total

| Phase | Runs | Temps |
|---|---|---|
| A — architecture | 36 | ~1 h |
| B — criblage | 80 | ~1 h 30 |
| C — base | 800 | ~3 h (parallélisé) |
| D — émulateur | 0 | minutes |
| E — validation | ~30 | ~2 h (dont 3D) |
| **Total** | **~950** | **~8 h de machine**, étalées sur 2-3 nuits |

---

## 8bis. L'algorithme, brique par brique — justification sourcée

```
1.  LHS maximin de 60 jeux sur les paramètres retenus
2.  simuler les 4 essais par jeu (3 graines pour le brésilien)
3.  répéter jusqu'à convergence :
        a. entraîner un GP (Matérn 5/2, ARD) par sortie
        b. choisir 10 nouveaux jeux par acquisition EHVI multi-objectif
        c. les simuler, les ajouter à la base
4.  NSGA-II sur l'émulateur → front de Pareto
5.  MCMC sur l'émulateur → postérieur, intervalles, corrélations
6.  validation par runs directs + prédiction σ₃ = 75 et 100
```

### (a) Le cadre : émulateur GP + calibration bayésienne + discrépance

**Kennedy & O'Hagan (2001)**, *Bayesian calibration of computer models*,
JRSS-B **63**(3), 425-464 — le cadre de référence du domaine. Deux processus
gaussiens : l'**émulateur** η qui interpole le code coûteux avec son incertitude
d'interpolation, et la **discrépance** δ qui représente l'écart structurel
modèle/réalité. Le postérieur des paramètres intègre les trois sources
d'incertitude : interpolation, discrépance, erreur d'observation.

→ **Pourquoi c'est le bon cadre ici** : notre enveloppe expérimentale est
concave alors que les joints sont Mohr-Coulomb — il y a un biais structurel
connu. Sans terme de discrépance, le postérieur serait artificiellement
resserré autour d'un compromis biaisé (Brynjarsdóttir & O'Hagan, 2014).
Les fondations du krigeage pour codes coûteux : **Sacks, Welch, Mitchell &
Wynn (1989)**, *Design and Analysis of Computer Experiments*, Statistical
Science 4(4), 409-435.

### (b) Le GP plutôt qu'un réseau de neurones

Régime « petites données » : 5-9 paramètres, quelques centaines de points.
Le GP donne une **variance de prédiction** — sans elle, ni enrichissement
adaptatif ni postérieur honnête. **Bu et al. (2026)** obtiennent d'ailleurs
R² = 0,96 avec leur GP à l = 4 mm ; ils ne l'écartent que dans un régime
(99 roches, 3 456 points, 7 paramètres) qui n'est pas le nôtre — et leur
propre résultat montre que le meilleur R² **direct** (forêt aléatoire) est
mauvais en **inverse**, ce qui disqualifie le critère qu'ils auraient utilisé
pour éliminer le GP.

### (c) ARD comme criblage — remplace l'analyse de Pearson

Le noyau à **longueurs de corrélation par dimension** (Automatic Relevance
Determination) apprend une échelle par paramètre ; l'inverse de la longueur
mesure la pertinence de la dimension, et sur des entrées standardisées les
longueurs se lisent directement comme des mesures d'importance — procédure
de *screening* documentée comme telle dans le **toolkit MUCM / mogp-emulator**
(`ProcAutomaticRelevanceDetermination`), origine **Neal (1996)** puis
**Williams & Rasmussen**, *Gaussian Processes for Machine Learning* (2006),
§5.1.

→ **Avantage sur le Pearson de Jiang et al. (2025)** : l'ARD capture les
effets **non linéaires et les interactions**, là où un coefficient de Pearson
ne voit que la corrélation linéaire — et il sort gratuitement de l'émulateur,
sans les 80 runs d'un criblage mono-variable.
→ **Réserve à garder** : l'ARD naïf a des limites documentées pour la
sélection de variables (**Paananen et al., arXiv:1712.08048**) — on croisera
donc les longueurs ARD avec une sensibilité de la prédictive avant de figer
un paramètre.

### (d) L'enrichissement adaptatif plutôt qu'un plan unique

Origine : **Jones, Schonlau & Welch (1998)**, *Efficient Global Optimization
of Expensive Black-Box Functions*, J. Global Optim. 13, 455-492 (critère
Expected Improvement). Pour la **calibration** spécifiquement :
- **Teixeira et al. (2025)**, *Surrogate-aided Bayesian calibration with
  adaptive learning strategies*, Mech. Syst. Signal Process. **237**, 113014 ;
- **Adaptive GP surrogates for Bayesian inference** (arXiv:1809.10784) et
  **Posterior sampling with adaptive GP in Bayesian parameter identification**
  (arXiv:2411.17858) : le GP est raffiné par apprentissage actif, les points
  d'entraînement sont choisis pour **représenter au mieux le postérieur à
  échantillonner** — « réduction significative de l'effort de calcul par
  rapport aux plans statiques » ;
- **Adaptive multi-output GP for large-scale parameter estimation**,
  Engineering Computations **41**(6), 2024 ;
- côté DEM : la littérature récente identifie explicitement le **nombre de
  simulations d'entraînement** comme le verrou de la calibration par
  surrogate (Sci. Rep. **14**, 2024, transfer learning pour DEM).

→ **Ce que ça change** : les trois articles tirent tout d'un coup (3 456,
328, 231 runs) ; l'enrichissement place les points là où l'émulateur est
incertain **dans la zone plausible** — d'où l'estimation ~120 runs au lieu
de 800.

### (e) NSGA-II sur l'émulateur + acquisition EHVI

- **Deb et al. (2002)**, *A fast and elitist multiobjective genetic algorithm:
  NSGA-II*, IEEE Trans. Evol. Comput. 6(2), 182-197 — l'algorithme employé
  par **Ye et al. (2025)** ;
- **Emmerich, Deutz & Klinkenberg (2011)** : *Hypervolume-based expected
  improvement — monotonicity properties and exact computation* (EHVI) ;
  version parallèle différentiable : **Daulton, Balandat & Bakshy (2020)**,
  NeurIPS (qEHVI).

→ **Pourquoi les deux et pas l'un ou l'autre** : NSGA-II est robuste mais
**exige un très grand nombre d'évaluations** — inacceptable sur le simulateur,
gratuit sur l'émulateur. Les comparaisons publiées (TSEMO vs ParEGO vs EHVI
vs NSGA-II, 9 problèmes à budget 150 évaluations) montrent que l'optimisation
bayésienne est **bien plus économe en échantillons** : on l'utilise donc pour
CHOISIR les runs (EHVI), et NSGA-II pour explorer le front une fois
l'émulateur en place.

### (f) Le front de Pareto plutôt qu'une somme pondérée

Justification empirique interne aux trois articles : le conflit
traction/compression est mesuré indépendamment par **Ye** (il doit privilégier
l'UCS), **Bu** (BTS R² 0,74 → 0,90 en retirant c et φ des cibles) et **Jiang**
(R² 0,65 sur la résistance en traction, sa pire sortie). Une somme pondérée
cacherait ce compromis dans un choix de poids arbitraire ; le front le rend
**visible et chiffré**.

### (g) Option multi-fidélité (phase E)

**Adaptive sampling of multi-fidelity Gaussian process** (arXiv:1907.11739)
et EHVI multi-fidélité (Emmerich et al., 2021) : apprendre la relation
2D (40-120 s) → 3D (30-60 min) sur quelques dizaines de paires pour prédire
le 3D au prix du 2D. Pertinent pour le contrôle 3D final ; à trancher après
la phase A.

---

## 9. Ce qui reste à trancher (validation demandée)

1. **Architecture** : lancer la phase A comparative (recommandé), ou aller
   directement sur A3 (bulk DP-DFH figé, joints seuls) ?
2. **Élasticité cible** : 77,7 GPa / 0,29 (fit des 12 branches triaxiales,
   recommandé) — ou 52 / 0,25 pour rester cohérent avec la carte DP-DFH de
   la thèse, sachant qu'en A3 le bulk EST cette carte ?
3. **Traction** : cible BTS 10,3 MPa (recommandé, on simule un brésilien) ou
   σt 18,3 MPa de Saadati/Shariati (traction directe équivalente Weibull) ?
4. **Domaine** : calibrer sur σ₃ ≤ 50 et prédire 75/100 (recommandé), ou
   calibrer sur les quatre confinements ?
5. **Pondération du conflit traction/compression** : front de Pareto et choix
   documenté a posteriori (recommandé), ou priorité explicite donnée d'emblée
   à la compression (choix de Ye) ?
