# LOT 4 — Bilan de rockim contre l'état de l'art Imperial
# et LE point qui commande tout le reste

*Fiche du 2026-08-29. Livrable 4 du [brief](../MISSION_etat_de_l_art_2026-08-29.md) §5.*

---

## 0. HONNÊTETÉ SUR LA COUVERTURE DE CET AUDIT

L'audit systématique prévoyait onze éléments de formulation, chacun audité puis
**contre-audité**. **Quatre ont abouti** — loi de joint, pénalité, contact et
détection, rampe de naissance — avant qu'une limite de session ne coupe les
autres. Les sept restants (frottement, DIF, pulvérisation, fragments et énergie,
insertion et périmètre, maillage et pas de temps, avantages de rockim) ont été
traités **par lecture directe du rédacteur**, plus ciblée et donc moins
exhaustive. **Les contre-audits n'ont pas tourné du tout.**

Conséquence de méthode : les statuts des quatre éléments audités portent des
numéros de ligne vérifiés et sont contradictoires ; **les autres sont des
constats de première lecture, à reprendre**. Ils sont marqués **[1re LECTURE]**.
Le §9 dit ce qui reste à faire.

---

## 1. LE RÉSULTAT PRINCIPAL — et ce n'était pas la question posée

Le brief demandait ce qui manque à rockim. **La réponse la plus importante est
que rockim porte, dans son propre code source, 117 attributions à un code qui
n'est pas celui d'Imperial.**

    grep -rc "Y3D[a-z]*\.c|[Ss]olidity" src include tools bench_impact/configs

| fichier | occurrences |
|---|---|
| `src/Fdem3dSolver.cpp` | 46 |
| `src/FdemSolver.cpp` | 43 |
| `include/rockim/Fdem3dSolver.hpp` | 16 |
| `include/rockim/FdemSolver.hpp` | 12 |
| `tools/verify_suite.py` | 14 |
| decks `bench_impact/configs/*.cfg` | 41 |

La CORRECTION 2 du 29/08 avait purgé les **bilans**. **Le code, les en-têtes, la
suite de vérification et les decks n'ont pas été touchés.** Or c'est là que les
attributions comptent : ce sont elles qui finiront dans le manuscrit, dans les
figures et dans les réponses aux relecteurs.

### 1.1 Table de rachat — ce qui peut désormais être re-sourcé

Les lots 2 et 3 ont établi la formulation publiée. Chaque attribution peut donc
être rejugée :

| clé de rockim | attribution actuelle | verdict |
|---|---|---|
| **`jointFailRule = majority`** | « la règle de Solidity (`nfail>1`, Y3Dfd.c l. 1175) » | **RACHETÉE.** Manuscrit UCL **p. 14** : « A joint element is labelled as failed when **at least two integration points have zero stress components** ». Même règle, vraie source. **Remplacer la citation.** |
| **`jointDeltaC = guo`** | Guo 2014 éq. 2.30 | **DÉJÀ BONNE.** Confirmée par le manuscrit UCL éq. 10 p. 12 (`Gf ≈ ⅓ f δc`, aire de la parabole approchant l'exponentielle). |
| **plancher `st = max(2 sp, 3 GfII/dpefs)`** | « Y3Dfd.c l. 1110-1126 » + Guo éq. 2.24/2.30 | **PARTIELLEMENT RACHETÉE.** La partie Guo tient ; le plancher `2 sp` n'a **aucune source publiée** — à déclarer comme choix de rockim. |
| **`jointDeltaC = solidity`** | « Y3Dfd.c l. 1099 » | **NON RACHETABLE.** Aucun équivalent publié. À renommer, et à documenter comme *convention relevée dans une transposition d'équipe, non attribuable à Imperial*. |
| **taux de déformation `none`** | « ce que fait LEUR code : Y3Dfd.c l. 1448 » | **NON RACHETABLE.** Le mot « LEUR » est faux. Aucune source publiée ne décrit la mesure du taux. |
| **`gcBirth = penalty`** | « Y3Did.c l. 915-964 » (13 occurrences) | **NON RACHETABLE.** Et c'est la voie de naissance de contact que le §7 désigne comme correctif candidat : elle ne peut pas s'appuyer sur cette source. |
| **`contactMu.<phase>`, règle de paire = MINIMUM** | « Solidity Y3Did.c l. 1292 » | **NON RACHETABLE — et c'est un cas d'école.** Le [lot 2c](2026-08-29_lot2c_frottement_tangentiel.md) §4 a établi qu'**Imperial ne publie AUCUNE règle de paire**, sur cinq sources. Le « minimum » est donc un choix de rockim. Il peut être bon ; il ne peut pas se réclamer d'eux. |
| **`contactDamageCoupling = solidity`** | « raideur normale ET frottement multipliés par `d_fact = min(1−D_i, 1−D_j)`, effondrement /1000 sous 0,041 (Y3Did.c l. 995, 1044, 1263-1265) » | **NON RACHETABLE, ET C'EST LE PLUS GRAVE.** Le [lot 2b](2026-08-29_lot2b_couplage_endommagement_contact.md) a établi que le mot « penalty » **n'apparaît pas une seule fois** dans l'article de pulvérisation, et que **le seul couplage (1−D) publié porte sur la contrainte d'élément**. rockim a donc implémenté, sous un nom qui promet une réplication, un mécanisme **qu'aucune publication d'Imperial ne décrit**. |

### 1.2 L'action

**Effort FAIBLE, valeur ÉLEVÉE, et c'est la seule action de ce lot qui touche à
l'intégrité scientifique du dépôt.** Trois gestes :

1. remplacer les citations **rachetables** par leur vraie source (article, page,
   équation) ;
2. renommer les clés qui promettent une réplication qu'elles ne font pas —
   `contactDamageCoupling = solidity` devient p. ex.
   `contactDamageCoupling = hypothese_rockim`, `jointDeltaC = solidity` devient
   `jointDeltaC = transposition` ;
3. porter en tête des deux solveurs l'avertissement de la CORRECTION 2, pour que
   la prochaine session ne recommence pas.

**Tant que ce n'est pas fait, toute phrase du manuscrit qui dit « comme
Imperial » en s'appuyant sur ces clés est indéfendable.**

---

## 2. LE TABLEAU DE BILAN

Statuts : **PRÉSENT** conforme · **PARTIEL** incomplet ou autrement ·
**ABSENT** à écrire · **DIVERGENT** autre choix, assumé et documenté.
**Blocage** = sans cela un résultat publié d'Imperial n'est pas reproductible.

### 2.1 Loi de joint cohésive — **audité et vérifié**

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| fonction d'adoucissement z(D), a=0,63 b=1,8 c=6,0 | **PRÉSENT** | `YanSoftening.hpp:53-68`, défauts `:46-50` ; `FdemSolver.cpp:645-654` accepte `munjiza` comme **alias** de `yan` — le dépôt avait déjà tranché que Yan 2023 reprend la z-curve de Munjiza | aucun | non |
| σ trois branches (compression 2·pj, durcissement parabolique, adoucissement z·ft) | **PRÉSENT** (opt-in) | `FdemSolver.cpp:4003` (`pjC = 2.0*J.pj`), `:3987-3990` (parabole), `:3991` (min) ; défaut `linear`, clé `jointElastic = parabolic` | aucun | non |
| D mixte (éq. 13) | **PRÉSENT** (opt-in) | `FdemSolver.cpp:3958, 3917, 3968` ; forme littérale sous `jointShearUnload = origin` | aucun | non |
| irréversibilité D = Dmax et plafond à 1 | **PRÉSENT** | `FdemSolver.cpp:3970` | aucun | non |
| rupture : ≥ 2 points d'intégration sur 3 | **PRÉSENT** (opt-in) | `Fdem3dSolver.cpp:2777-2783`, endommagement **par point** `:2552` ; clé `jointFailRule = majority`, qui **exige** `jointQuadrature = midedge` | aucun | non |
| `min(branche élastique, z·ft)` : l'endommagement de cisaillement coupe aussi le pic de traction | **DIVERGENT** | `FdemSolver.cpp:3991` et son commentaire `:3973-3975`. L'éq. 11 d'Imperial porte **ft plein** en pré-pic | faible | non |
| « rompu » = D ≥ 1, contre « zero stress components » | **PARTIEL** | à D = 1 un point rockim transmet encore σ = pjC·δn en compression et tanφ·|σn| en cisaillement | faible | **à trancher** |

### 2.2 Pénalité de joint — **audité et vérifié**

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| dimension et entrée dans la loi | **PRÉSENT** | `J.pj = pf*E/h` [Pa/m], `FdemSolver.cpp:2122`, `Fdem3dSolver.cpp:1521` ; `pf·E` est homogène au p₀ d'Imperial | aucun | non |
| niveau : 20 par défaut, **26,32 au deck de réplique** | **DIVERGENT documenté** | `impact_imperial.cfg:209, 228` ; justification `FdemSolver.cpp:622-627` ; dérivation `DOCUMENTATION_rockim.md:321-332` | aucun | non |
| **h = diamètre inscrit (6V/A) contre « mean length of the EDGES »** | **PARTIEL** | `Fdem3dSolver.cpp:1392` (`6.0*e.V0/Atot`) contre manuscrit UCL p. 12 éq. 6 | faible | **OUI** |
| impression de la pénalité sous `intrinsic`, **3D** | **PRÉSENT** | `Fdem3dSolver.cpp:544-562`, avec avertissements croisés dans les deux sens | aucun | non |
| impression de la pénalité sous `intrinsic`, **2D** | **ABSENT** | `git show --stat d4be57b` : **1 file changed**, 3D seulement. En 2D, sous le défaut `intrinsic`, la raideur appliquée n'apparaît **nulle part** | faible | non |
| impression de la raideur `pj` réellement portée (Pa/m) | **ABSENT** | aucun `std::cout` de `J.pj`. Sur un maillage gradué 1→10 mm, `pj` varie d'un facteur 10 et la bannière « 20 E/h » ne dit pas laquelle porte le joint le plus fin | faible | non |
| validation de `jointPenaltyFactor` | **ABSENT** | `pf = 0` donne `dnE = ft/0 = +inf` puis NaN, **en silence** — alors que dix lignes plus bas ft, coh, E et φ sont validés | faible | non |

### 2.3 Contact — **audité et vérifié**

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| potentiel de Munjiza, force distribuée (éq. 17) | **PRÉSENT** | `PotentialContact.hpp:116-208` (2D), `:427-610` (3D) ; conservation vérifiée à 3,7e-12 (2D) et 2,0e-8 (3D) | aucun | non |
| exclusion des couples à joint **vivant** | **PRÉSENT** | `FdemSolver.cpp:4589-4591`, `Fdem3dSolver.cpp:3173-3178` | aucun | non |
| détection **paresseuse** : rien à l'intérieur du continu | **PARTIEL** | jeu actif = faces extérieures + faces des joints **morts** (`FdemSolver.cpp:4287-4309`) — le point capital est acquis | — | non |
| naissance des couples par les **six groupes nodaux** (éq. 15-16) | **ABSENT** | aucune construction topologique ; la liste est **reconstruite géométriquement à chaque pas** | moyen | non |
| structure d'accélération « NBS » | **PARTIEL** | grille **dense** sur boîte fixe, tous les seaux vidés à chaque pas : coût O(cellules), pas O(N) — ~1,2e6 seaux pour un cube de 50 mm à h = 0,5 mm | faible | non |

> **Réserve de sobriété sur R4.** Le dépôt a **déjà mesuré** le plafond du gain :
> sur la percussion longue le poste dominant est l'intégration exacte des 33 M de
> clips **avec force** (≈ 2/3 du run) ; les clips vides ne pèsent que ~12 %
> (`DOCUMENTATION_rockim.md:183`). **R4 attaquerait 10-15 % du mur, pas un
> facteur 7.** Ma recommandation R4 du lot 3 doit être relue avec ce chiffre.

### 2.4 Rampe de naissance du contact — **audité et vérifié**

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| objectif de l'éq. 18 (pas de force sortie du néant) | **PRÉSENT** | traité depuis le chantier A3 ; clé `gcBirth` | aucun | non |
| forme : facteur **temporel** linéaire sur ~10 pas | **DIVERGENT mesuré** | rockim retranche un **décalage d'état**, pas un facteur de temps. Mesures du dépôt : **+936 J/m** sans relevé, **+179 J/m** avec rampe temporelle, **+27 J/m** avec le décalage d'état (`FdemSolver.hpp:706-712`) | — | non |
| **déclenchement conditionnel** (cisaillement sous compression) | **ABSENT** | rockim relève **toute** naissance, y compris en traction franche où il n'y a rien à relayer. Les deux ingrédients existent pourtant : `J.bmode` (`FdemSolver.cpp:4170`) et le signe de `J.fDeath` — mais `bmode` est déclaré **sortie seule** (`FdemSolver.hpp:203`) | faible | non |
| durée comptée en pas | **PARTIEL** | `gcBirthTau = 1e-6 s` (`FdemSolver.cpp:1020`) ; au dt de 1,93e-9 s cela couvre **≈ 518 pas**, ~50× le ntotal ≈ 10 d'Imperial — et **aucune ligne du dépôt ne justifie ni ne balaie cette valeur** | faible | non |

### 2.5 Les sept autres éléments — ~~[1re LECTURE]~~ **SUPERSÉDÉ par le §2.6**

*Ce tableau était une première lecture. Le §2.6, ajouté le 2026-08-29 au soir, le
remplace : les sept éléments y sont repris avec numéros de ligne vérifiés.*


| élément | constat | à reprendre |
|---|---|---|
| **DIF** | **PRÉSENT et au-delà** : `YangDif.hpp` code les deux lois, exposant de traction **paramétré**, parité 2D/3D garantie par construction, quatre contrôles de non-régression. Le dépôt avait diagnostiqué la coquille **le 2026-08-18**, avant moi et mieux (cf. lot 2b, encadré de crédit) | vérifier que le DIF s'applique bien à **G_I et G_II**, et pas seulement aux résistances |
| **pulvérisation** | à confronter à l'éq. 3-4 de 2026 | la définition de δ_m (`Fdem3dSolver.cpp:2274`) et son homogénéité |
| **couplage D → contact** | **`contactDamageCoupling = solidity` existe** et fait exactement ce qu'aucune publication ne décrit (§1.1) | requalifier, puis décider si on le garde comme **extension assumée** |
| **règle de paire du frottement** | **`contactMu.<phase>` avec règle du MINIMUM**, sourcée sur le code non-Imperial | requalifier ; le choix reste défendable en propre |
| **périmètre des joints** | 34 507 joints sur les trois corps contre la roche seule chez Imperial (lot 3 §3) | **vérifier s'il existe déjà un filtre par phase** — si oui, R1 coûte une clé de deck |
| **bilan d'énergie** | rockim mesure et **imprime** le résidu ; Imperial l'obtient **par soustraction** et le dit | chiffrer précisément, c'est une avance de premier ordre |
| **fragments** | à vérifier | existence d'une identification par connectivité et d'un post-traitement de retrait |

### 2.6 LES SEPT ÉLÉMENTS REPRIS À LA MAIN — vérifiés ligne à ligne

*Le workflow d'audit ayant échoué deux fois sur une limite de session, ces sept
éléments ont été lus directement. Chaque ligne a été ouverte et relue.*

#### Frottement tangentiel

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| ressort tangentiel + plafond de Coulomb | **PRÉSENT** | `Fdem3dSolver.cpp:3310-3315` : `Ft -= potKt_ * dt_ * vt;` puis `cap = ctcMu(eLo,eHi)*Fn; if (Ftn > cap) Ft *= cap/Ftn;`. C'est la forme intégrale de `f_t = −k_t δ_t` avec écrêtage de Coulomb — l'éq. 8-9 d'Imperial | aucun | non |
| **terme visqueux tangentiel `− η v_t`** | **ABSENT** | rockim n'a que le ressort. `potXi_` (`FdemSolver.cpp:836`) est une fraction d'amortissement critique du contact **NORMAL**, et **2D seulement**. Rien de tangentiel | faible | non |
| **convention de signe** | **PARTIEL** | rockim prend `−k_t δ_t` — la convention du **chapitre de 2017**, pas celle de l'article de **2009** qui écrit `+k_t δ_t` (lot 2c §3ter). Le choix n'est **documenté nulle part** | faible | non |
| **règle de paire** | **DIVERGENT non sourçable** | `Fdem3dSolver.hpp:524-529` : `if (muPerPhase_) { m = min(muPhase_[A], muPhase_[B]); }` — **le plus faible gouverne**, commenté « Solidity Y3Did.c l. 1292 ». **Imperial ne publie aucune règle** (six sources) | aucun | non |

> **JE CORRIGE UNE ERREUR QUE J'AI RÉPÉTÉE TOUTE LA SESSION.** J'ai écrit à
> plusieurs reprises « le `k_t = 2/7` de rockim ». **C'est inexact.** La raideur
> tangentielle de rockim vaut `potKt_ = potTangentFactor × E × h`, avec
> **`potTangentFactor` = 1,0 par défaut** (`Fdem3dSolver.cpp:701`,
> `FdemSolver.cpp:835`). Le 2/7 n'est pas dans le code : c'est le **rapport
> k_t/k_n imposé par les decks d'impact** — `potTangentFactor = 1.4286 = 5,0 × 2/7`
> (`impact_imperial.cfg:377-379`) — et sa seule source est le commentaire
> `ktss = 2.0/(7.0)*d1pepe[icoup]`, qui cite le code non-Imperial.
> **À ajouter à la table de rachat du §1.1, colonne NON RACHETABLE** : ni
> l'article de 2009 ni le chapitre de 2017 ne donnent de valeur de k_t.

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| **banc analytique du frottement** (rectangle glissant, `L = v²/2µg`) | **ABSENT** | `tools/verify_suite.py` a bien des contrôles de frottement (`jointdeath_friction_2d`, `jointResidualMu`, `jointFrictionScaled`) mais **tous portent sur le frottement de JOINT**, aucun sur le chemin **tangentiel de contact**. Le banc publié (lot 2c §3ter) est à solution analytique fermée | faible | non |

#### DIF — **le point le plus propre du dépôt**

| point | statut | preuve |
|---|---|---|
| lois et bornes | **PRÉSENT** | `YangDif.hpp:45-58`, exposant de traction **paramétré**, parité 2D/3D par construction |
| **points d'application** | **PRÉSENT, strictement conforme** | `Fdem3dSolver.cpp:1779-1789`, `stampDif` : `J.ft *= dT; J.Gf *= dT; J.coh *= dC; J.GfII *= dC;` — **exactement** la règle d'Imperial (traction → f_t et G_I ; compression → cohésion et G_II), et le **frottement interne n'est pas touché** |
| rafraîchissement continu | **PRÉSENT** | `refreshDif` (`:1796-1806`) repart des valeurs de base (`snapBase`) : jamais de composition, le facteur peut redescendre |

**Rien à faire. C'est conforme, testé et documenté au-delà de la source.**

#### Pulvérisation et couplage

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| δ_m dimensionnellement correct | **PRÉSENT** | `Fdem3dSolver.cpp:2362` : `dm = hEl_[eI] * sqrt(2.0/3.0) * ed.norm()` — une **longueur**, conforme à δ_m d'Imperial | aucun | non |
| **couplage D → frottement** | **DIVERGENT non sourçable** | `Fdem3dSolver.hpp:530-545`, sous `cplMode_` : `dr = min(1−D_A, 1−D_B)` **et `if (dr < 0.041) dr *= 1e-3`**. rockim reproduit le nombre magique 0,041 **et** l'effondrement par 1000 du code non-Imperial. **Aucune publication d'Imperial ne décrit ce mécanisme** (lot 2b) | — | non |
| couplage D → pénalité de contact | **présent en opt-in**, même statut | `jointContactPenalty = adaptive` | — | non |

> Le couplage est **continu**, pas en tout-ou-rien — ma première lecture du lot 4
> le supposait binaire. Il est plus proche de la forme du code non-Imperial que
> je ne le pensais, ce qui **aggrave** le problème d'attribution du §1.1 au lieu
> de l'atténuer : ce n'est pas une inspiration lointaine, c'est une
> transcription, constante magique comprise.

#### Bilan d'énergie — **rockim est très au-dessus**

| point | statut | preuve |
|---|---|---|
| postes ventilés et imprimés | **PRÉSENT** | `Fdem3dSolver.cpp:4411-4429` : KE initiale → KE bloc, poste éléments, dont **visqueux (2µD)** et dont **pulvérisation (bulkDamage)**, avec le nombre d'éléments pulvérisés |
| **résidu mesuré ET run interrompu** | **PRÉSENT, sans équivalent** | `Fdem3dSolver.cpp:2262-2276` : `resid = (ke − keInit_) − sumW` ; si le résidu dépasse `budgetAbortPct` de l'échelle, le run **s'arrête** en imprimant le **nœud le plus rapide et sa position** |

**Imperial obtient amortissement et erreur PAR SOUSTRACTION (A952 p. 3). rockim
mesure le résidu et refuse de continuer s'il dérive.** C'est l'écart le plus net
de tout le bilan, et il est à l'avantage du dépôt.

#### Insertion et périmètre

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| deux schémas (intrinsèque, adaptatif) | **PRÉSENT + capacité en plus** | Imperial n'a que l'intrinsèque | aucun | non |
| **filtre de joints par phase ou par corps** | **ABSENT** | recherche : `jointPhase`, `jointsIn`, `noJoint`, `jointBodies`, `jointMaterial`, `skipJoint` → **aucune occurrence** dans `src`, `include` ni les decks | **MOYEN** | non |

> **CORRECTION DE MON PLAN.** L'étape 5 du lot 5 (« restreindre les joints à la
> roche ») supposait « 0,5 à 1 j, faible **si un filtre par phase existe** ».
> **Il n'existe pas.** L'effort est **MOYEN** : il faut ajouter le filtre à la
> pose des joints, dans les deux solveurs, avec sa clé, sa bannière et son
> contrôle de non-régression.

#### Maillage, pas de temps, garde-fous

| point | statut | preuve | effort | blocage |
|---|---|---|---|---|
| garde-fou crack-band | **PRÉSENT, sans équivalent publié** | `MatLaw.cpp:1300-1313` : lève si `Gf/(lcMax·ft) − 0,5·ft/E ≤ 0,05·ft/E`, en nommant `E·Gf/ft²` dans le message | aucun | non |
| borne CFL sur le **vrai** diamètre inscrit | **PRÉSENT** | `Fdem3dSolver.cpp:2131-2133`, avec l'histoire du bug du 2026-08-07 en commentaire (2,4 MJ d'énergie de bloc pour 16 J incidents) | aucun | non |
| borne diffusive du terme visqueux, **élément par élément** | **PRÉSENT** | `Fdem3dSolver.cpp:2134-2146` | aucun | non |
| **k_t dans le budget de pas de temps — 2D** | **PRÉSENT** | `FdemSolver.cpp:3134` : `if (contactPot_) kContact = max(kContact, max(potP_, potKt_));` | aucun | non |
| **k_t dans le budget de pas de temps — 3D** | **ABSENT** | `Fdem3dSolver.cpp:2115-2122` : `dtMin = min(dtMin, 2·sqrt(m/(K + nExtra·kp_)))`. **`potKt_` n'y entre pas**, et `kContact` n'existe pas dans ce fichier | faible | non |

> **C'est la question que la source primaire du frottement a fait naître, et la
> réponse est mauvaise.** Xiang, Munjiza, Latham & Guises (2009) p. 677
> avertissent que « in order to reduce the numerical error for calculation of
> **tangential forces, the smaller time step is required** », alors que le même
> calcul **sans frottement** est stable au pas plus grand. Le 2D compte `potKt_`
> dans son budget ; **le 3D ne le compte pas** — et c'est le 3D qui porte
> l'impact. Conséquence : dès que le frottement travaille sous l'insert, la
> marge de stabilité n'est pas celle qu'on croit. **Le correctif est une ligne**,
> le miroir exacte de `FdemSolver.cpp:3134`.

---

---

## 3. CE QUE ROCKIM A ET QU'IMPERIAL N'A PAS

Le brief l'exige, et ce n'est pas une politesse — la liste est longue et plusieurs
entrées sont scientifiques, pas cosmétiques.

1. **Une suite de non-régression de 98 contrôles**, pas 42 comme le brief le
   croyait : 42 *fast*, 45 *full*, 8 *all* (`tools/verify_suite.py`). **Imperial
   ne publie aucune suite.**
2. **Un bilan d'énergie fermé et imprimé**, là où Imperial obtient amortissement
   et erreur numérique **par soustraction** et l'écrit (ARMA 24-0952 p. 3).
3. **Le refus plutôt que le silence.** `jointElastic = parabolic` lève une
   exception si `jointSoftening` n'est ni `yan` ni `munjiza`
   (`FdemSolver.cpp:196-205`), avec le commentaire « LECON DU 2026-08-25, deux
   fois dans la meme seance : une capacite active et muette est indiscernable
   d une capacite inerte ». **Aucun code publié ne fait ça.**
4. **Des avertissements croisés** quand un deck pose la clé de l'autre schéma
   d'insertion (`Fdem3dSolver.cpp:552-561`).
5. **Le pas de temps budgète explicitement la pénalité**, facteur 2 de la
   parabole compris (`FdemSolver.cpp:3105-3112`).
6. **Une vérification analytique au point matériel** : `tools/yan_point.cpp`
   intègre σ(δn) par trapèzes et compare l'aire à G_fI — la loi expédiée est
   testée, pas une ré-implémentation.
7. **Un garde-fou crack-band** (`MatLaw.cpp:1304-1314`) qui lève si le plus gros
   élément dépasse E·G_f/f_t². **Imperial n'a aucun critère d'objectivité
   publié** — leur Table 2 de balayage de maillage ne publie que le nombre
   d'éléments et le temps CPU (lot 3 §8.3).
8. **La parité 2D/3D** posée en principe et instrumentée (`YangDif.hpp` la nomme
   comme raison d'être).
9. **Un diagnostic instrumenté** de la pathologie propre à la pénalité
   intrinsèque : part des joints au-dessus de D = 0,01 au pic, D moyen,
   « diffuse ratcheting » (`FdemSolver.cpp:6840-6843`).
10. **Un plafond d'impulsion dur** sur le contact (20 m/s par pas et par nœud,
    `FdemSolver.cpp:5019-5026`).
11. **Un levier mesuré et séparé** : pénalité contre schéma d'insertion, à
    pénalité égale +1,5 point, le reste étant de la pénalité
    (`BILAN_interference_2026-08-29.md:207-213`). **Imperial ne fait aucune étude
    de ce genre.**
12. **La provenance du 3000 GPa retrouvée** — règle de Turon, Dávila, Camanho &
    Costa (2007), K = α·E/t avec α ≈ 50 — **et une contradiction interne
    d'Imperial mise au jour** : Guo recommande E ≤ p₀ ≤ 10E deux phrases après
    avoir cité Turon, et les auteurs de l'article ont suivi Turon
    (`DOCUMENTATION_rockim.md:321-332`).
13. **Une prédiction confirmée** : l'exposant 0,1707 dérivé de la figure 2(b) le
    2026-08-18, imprimé 0,17 par l'article de 2026 (lot 2b).

**Les entrées 7, 11, 12 et 13 sont des résultats scientifiques**, pas de
l'ingénierie logicielle. Elles ont leur place dans le manuscrit.

---

## 4. JE CORRIGE MA PROPRE RECOMMANDATION R2 — elle était fausse

Le [lot 3](2026-08-29_lot3_insertion_maillage.md) §10 recommandait de « porter la
pénalité de joint à ≈ 50 E » au motif que rockim serait « 2,6 fois trop bas ».
**Deux erreurs.**

**Erreur 1 — le facteur 2 de la convention.** Le dépôt avait **déjà** établi
(`DOCUMENTATION_rockim.md:321-332`) que sous la convention parabolique de Guo la
raideur vaut p₀/(2h), donc l'équivalent rockim de leur p₀ = 3000 GPa est
**p₀/(2E) = 26,32**, et **le deck de réplique est déjà à 26,32**
(`impact_imperial.cfg:209, 228`). Le « 20 » que je comparais est le **défaut**,
pas la valeur de réplique. Appliquer 50 aurait été **deux fois trop raide**.

**Erreur 2 — et elle est neuve.** L'équivalence 26,32 fait s'annuler h des deux
côtés. **Elle n'en a pas le droit** : ce ne sont pas les mêmes longueurs.

| | définition | valeur pour un tétraèdre régulier d'arête a |
|---|---|---|
| **rockim** | diamètre de la sphère inscrite, 6V/A (`Fdem3dSolver.cpp:1392`) | **0,4082 a** |
| **Imperial** | « the **mean length of the edges** of the joint element » (manuscrit UCL p. 12, éq. 6) | **a** |

Rapport **2,4495**. À `pf` égal, rockim divise par une longueur 2,45 fois plus
petite : **sa pénalité est 2,45 fois plus raide**. L'équivalence correcte est
donc

    pf = (p0 / 2E) x (h_inscrit / h_arete) = 26,32 x 0,4082 = 10,74

**Le deck de réplique est à 26,32, soit ≈ 2,45 fois trop raide.**

> **R2 corrigée.** Ne pas monter la pénalité : **la descendre**, ou mieux,
> **changer la mesure de h**. Deux voies :
> * **voie deck**, immédiate : `jointPenaltyFactor = 10.74` sur le banc de
>   réplique, en documentant que c'est 26,32 corrigé du rapport de longueurs ;
> * **voie code**, propre : calculer h comme la **longueur moyenne des arêtes de
>   la facette de joint**, ce qui rend `jointPenaltyFactor` directement
>   comparable au p₀ publié et supprime le facteur pour toujours. Effort faible,
>   mais **change les résultats** — donc opt-in et bannière, selon la règle du
>   dépôt.
>
> **Réserve** : le rapport 0,4082 vaut pour un tétraèdre **régulier**. Sur un
> maillage réel il varie. C'est un argument de plus pour la voie code.

**C'est le seul point BLOQUANT trouvé sur la formulation** : tant que la
longueur de référence diffère, l'ouverture au pic et la complaisance de joint ne
peuvent pas coïncider avec les leurs, quel que soit le facteur.

---

## 5. LE VRAI BLOCAGE DE LA RÉPLICATION N'EST PAS UNE ÉQUATION MANQUANTE

Le dépôt a mesuré que **la branche normale du contact INJECTE de l'énergie** :
**+3,66 J sur l'impact 3D P1 (6,9 % de KE₀), et 11,1 J (20 %) en insertion
intrinsèque** (`BILAN_interference_2026-08-29.md:118-131`,
`HANDOFF_2026-08-29.md:134-136`), l'explication mécanique étant le facteur de
naissance du contact.

**Tant qu'un canal injecte 20 % de l'énergie d'entrée, le partage publié — 2,6 %
à la fissuration, 64,9 % au frottement (ARMA 24-0952) — n'est pas
reproductible.** Aucune des recommandations du lot 3 n'attaque ça.

Le correctif candidat existe et n'a **jamais été testé sur l'impact 3D** :
`gcBirth = penalty`, qui recale la pénalité de la paire sur la force du joint
mourant (`FdemSolver.cpp:4606-4632`). Il exige `contact = potential`
(`FdemSolver.cpp:1033-1041`) — ce que les decks d'impact posent déjà.

**C'est l'essai numéro un.** Et il faut le re-sourcer au passage (§1.1).

---

## 6. PLAN D'ACTION ORDONNÉ

| # | action | effort | gain mesurable | bloque la réplication ? |
|---|---|---|---|---|
| **A1** | Tester `gcBirth = penalty` sur l'impact 3D | faible | l'injection de 11,1 J par la branche normale doit tomber ; sans ça le partage 2,6/64,9 % est hors d'atteinte | **OUI** |
| **A2** | Aligner la longueur h (deck à 10,74, ou mieux : h = arête moyenne en code) | faible | ouverture au pic et complaisance de joint coïncidant avec les leurs | **OUI** |
| **A3** | Re-sourcer les 117 attributions, renommer les clés non rachetables | faible | intégrité du manuscrit | non, mais **indéfendable en l'état** |
| **A4** | Imprimer la pénalité en 2D et la raideur `pj` réelle (min/moy/max) | faible | rend visible ce qui gouverne dt et la complaisance | non |
| **A5** | Valider `jointPenaltyFactor` (> 0) | faible | supprime un NaN silencieux | non |
| **A6** | Restreindre les joints à la roche (R1 du lot 3) | **moyen** — aucun filtre par phase n'existe (§2.6) | raideur d'outil, pas de temps | non |
| **A7** | Armer la naissance de contact sur `bmode` et le signe de `fDeath` | faible | supprime le relevé en traction franche, où il n'est qu'une perte de portance | non |
| **A8** | Balayer `gcBirthTau` (518 pas contre ~10 chez eux) | faible | valeur aujourd'hui non justifiée | non |
| **A9** | Détection de contact par événement (R4 du lot 3) | moyen | **10-15 % du CPU seulement** — le poste dominant est ailleurs | non |
| **A10** | Rendre optionnel le `min(élastique, z·ft)` | faible | partage traction/cisaillement sous l'insert | non |
| **A11** | **Porter `potKt_` dans le budget de pas de temps du 3D** (une ligne, miroir de `FdemSolver.cpp:3134`) | **faible** | supprime une marge de stabilité illusoire dès que le frottement travaille — les auteurs de la source qualifient eux-mêmes le point d'« alarmant » | non |
| **A12** | Ajouter le banc analytique du rectangle glissant à `verify_suite.py` | **faible** | premier contrôle du chemin tangentiel de CONTACT, à solution fermée | non |

**A1 et A2 d'abord, et seulement elles, avant tout autre chantier.** Ce sont les
deux seuls blocages ; le reste est du raffinement ou de l'hygiène.

---

## 7. CE QUI RESTE À AUDITER

Les sept éléments non audités formellement (§0), et **surtout les onze
contre-audits**, qui n'ont pas tourné. Trois points précis à reprendre :

* le **couplage D → frottement** de rockim : continu ou en tout-ou-rien ? le WP7
  parlait d'une « porte par phase » ;
* le **périmètre des joints** : un filtre par phase existe-t-il déjà ? de sa
  réponse dépend si A6 coûte une clé ou une journée ;
* le **bilan d'énergie** : chiffrer précisément en quoi il est plus fermé que le
  leur, c'est une entrée de manuscrit.
