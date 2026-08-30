# CONTRE-AUDIT DU LOT 4 — les corrections, et le plan refait
*2026-08-30. Dix vérificateurs indépendants, un par section, chacun renvoyé au code
aux lignes citées. Ce document CORRIGE
[`biblio_insertion/2026-08-29_lot4_bilan_rockim.md`](../biblio_insertion/2026-08-29_lot4_bilan_rockim.md)
et le §6 de [`SYNTHESE_etat_de_l_art_2026-08-29.md`](../SYNTHESE_etat_de_l_art_2026-08-29.md).
**En cas de contradiction, c'est ce document qui fait foi.** Les originaux sont
conservés — on ne réécrit pas l'historique.*

---

## 0. LE RÉSULTAT, ET IL EST MAUVAIS

**129 verdicts. 52 confirmés — 40 %.**

| verdict | n |
|---|---|
| CONFIRMÉ | 52 |
| statut trop **flatteur** | 23 |
| **preuve fausse** | 21 |
| statut trop **sévère** | 17 |
| blocage mal jugé | 8 |
| effort irréaliste | 8 |

**6 critiques, 32 hautes, 34 moyennes.** Les critiques et les hautes sont traitées
aux §1-§8 ci-dessous ; **les 34 moyennes ont été dépouillées une par une le
2026-08-30 et sont au [§9](#9-les-34-verdicts--moyenne--dépouillés-un-par-un-le-2026-08-30)** —
onze d'entre elles changent une affirmation que l'on réutiliserait au manuscrit.
Le dépouillement des sources, lui, était à **91 %**
(614/673). L'écart n'est pas du bruit, et son motif est net :

| section | confirmés | comment elle a été produite |
|---|---|---|
| `loi-joint-penalite` | **13/17** | audit par agent **+ contre-audit** |
| `attributions-117` | 7/20 | à la main, en fin de session |
| `avance-rockim` | 5/15, dont **9 trop flatteurs** | à la main |
| `plan-action` | **2/17** | à la main |

**Ce qui est passé par un cycle complet tient ; ce qui a été fait à la main en
fin de session est le pire.** Et les neuf « trop flatteurs » sur les avantages de
rockim disent que la consigne « ne pas présenter rockim comme systématiquement en
retard » a été sur-corrigée en complaisance.

---

## 1. LES DEUX « BLOCAGES » TOMBENT

### 1.1 A2 — la longueur h n'est pas bloquante, et son critère était vide

**Trois défauts cumulés, tous établis :**

1. **Ce n'est pas un blocage.** Un blocage désigne une **capacité absente**.
   rockim possède la capacité : la longueur de référence diffère, mais **le
   facteur qui la multiplie est libre au deck**, et le document lui-même qualifie
   la voie deck d'« immédiate ». Un deck mal réglé n'est pas un blocage.
2. **L'ordre de grandeur ne tient pas.** Un écart de ~2 % sur le module effectif
   ne peut pas empêcher de retrouver un partage 2,6 / 64,9 %, quand le même
   document invoque par ailleurs une injection de 6,9 à 20 % de KE₀. **Deux
   ordres de grandeur entre les deux « blocages ».** Et rockim étant plus RAIDE,
   il est **plus proche** de ce qu'Imperial dit rechercher — ils prennent p₀ aussi
   haut que le pas de temps le permet.
3. **Le critère de réussite était tautologique.** δ_nE = f_t/pj = f_t·h/(pf·E), et
   avec pf = p₀/(2E) cela vaut **2 f_t h / p₀ quelle que soit la définition de h**.
   Le critère ne pouvait pas échouer. Et l'outil désigné pour le vérifier
   (`tools/yan_point.cpp`) **ne voit pas le maillage**.

**Et le chiffre corrigé est encore faux dans le sens flatteur.** Le rapport 2,4495
vaut pour un tétraèdre **régulier**, qui est l'**optimum** : tout maillage réel
s'en écarte **d'un seul côté**, donc le facteur réel est **toujours ≥ 2,4495**.
« Il varie » laissait croire à une dispersion à deux queues. La valeur d'action
est plutôt **≈ 9,6** que 10,74.

**Pire, et non relevé par moi** : `h = 0,5·(hEl_[eA] + hEl_[eB])` n'est même pas
« le diamètre inscrit » — c'est la **moyenne des diamètres inscrits des deux
tétraèdres adjacents**, une longueur d'**élément**, là où Imperial prend une
longueur de **facette** (« of the joint element »). L'écart est donc **plus grand**
que ce que j'écrivais, pas plus petit.

### 1.2 A1 — j'ai mis en titre un chiffre rétracté

Les **11,1 J / 20 %** viennent de `BILAN_interference_2026-08-29.md:118-131`,
c'est-à-dire **l'avant-correction**. Mon propre §3, avantage n° 11, cite les
lignes **207-213** — la mesure corrigée. **Le chiffre applicable à la
configuration adaptative est 3,66 J / 6,9 %.**

Et le `blocage = OUI` n'est **pas démontré** : rien n'établit que le partage
publié soit hors d'atteinte sans le correctif. La formulation honnête est : *une
injection résiduelle existe et vaut d'être supprimée ; qu'elle interdise le
résultat publié n'est pas montré.*

**Enfin, l'issue d'A1 est à deux faces**, et je n'en présentais qu'une. La seule
mesure comparative du dépôt (`DOCUMENTATION_rockim.md:516-519`, percussion **2D
sans débris**) est favorable de peu — résidu −0,9174 J/m en `ramp` contre
−0,9949 en `penalty`. Rien ne garantit le signe en 3D avec débris.

---

## 2. LE §1 — « MON RÉSULTAT PRINCIPAL » — EST FAUX CINQ FOIS

1. **Le chiffre 117 est la somme de quatre lignes sur six** de mon propre tableau,
   qui totalise **172**. Et le paragraphe suivant insiste précisément que les deux
   lignes omises « n'ont pas été touchées ». **Je me contredis à trois lignes
   d'intervalle.**
2. **Le bon chiffre est ~78.** Ce qui engage l'intégrité, ce sont les **78 lignes**
   citant `Y3D*.c l. NNNN` à l'appui d'un choix de modèle. Le simple mot
   « Solidity » (118 occurrences) **nomme un projet public** — ce n'est pas une
   faute. 78 est plus petit **et** plus précis.
3. **Ma commande de preuve ne renvoie rien.** `grep -rc "…|…"` en BRE cherche le
   caractère `|` littéral. Il faut **`grep -rEc`**. Un rapporteur qui la recopie
   conclut que le résultat n'existe pas.
4. **« La CORRECTION 2 avait purgé les bilans » est FAUX.** Aucun n'a été purgé —
   et le plus lu de tous, **`DOCUMENTATION_rockim.md` §5.4 quinquies**, le guide
   utilisateur, expose toujours la table `jointDeltaC = solidity` /
   `jointFailRule = majority` / `strainRateFilter = none` avec ses citations
   `Y3Dfd.c l. NNNN`, **sans le moindre avertissement**. C'est la pièce qui
   atteindra le plus sûrement un relecteur, et A3 l'aurait manquée.
5. **Le périmètre est incomplet** : `configs/` (distinct de `bench_impact/configs/`),
   `build_sol.cmd`, et **les specs WP6/WP7** — c'est-à-dire la spécification écrite
   de `contactDamageCoupling`, que je désignais comme « le plus grave » — sont
   **entièrement hors du filet**.

Et la prémisse elle-même était fausse : `solidity-solver-open` **est** le dépôt
public d'Imperial (LGPL-3.0). Déjà corrigé dans
[`A03_resourcer_attributions.md`](A03_resourcer_attributions.md).

---

## 3. TROIS ACTIONS DU PLAN ÉTAIENT DANGEREUSES OU IMPOSSIBLES

### 3.1 A6 « restreindre les joints à la roche » — **DÉSASSEMBLE L'OUTIL**

`src/Fdem3dSolver.cpp:1687`, verbatim : « **Nodes are already duplicated per
tet** ». La duplication des nœuds est **globale et inconditionnelle**, et les
joints sont créés pour **toute** face interne (`jt_.push_back(J)`, l. 1464).

**Donc supprimer les joints de l'acier et du carbure sous `insertion = intrinsic`
— le schéma du deck de réplique — laisse des tétraèdres à nœuds dédoublés que plus
rien ne relie. Le taillant tombe en poussière de tets tenus par le seul contact.**

**La bonne implémentation n'est pas un filtre.** Imperial fait autrement, et
rockim a déjà ce qu'il faut :

* **outil rigide** — c'est ce que fait D1 pour les plateaux brésiliens
  (« assumed to be rigid, material properties are not needed »), et rockim a
  `toolSig_` (contact de Signorini) qui **sort déjà `kp_` du budget de pas de
  temps** ;
* ou **souder les nœuds** de la phase outil (pas de dédoublement), ce qui suppose
  de toucher la construction du maillage.

**Effort réel : inconnu, et supérieur à 2 j.** L'action est **retirée du plan**
tant qu'elle n'est pas re-spécifiée.

### 3.2 A8 « balayer `gcBirthTau` » — **impossible en l'état**

Le deck de réplique **ne peut littéralement pas accepter la clé** : y écrire
`gcBirthTau` fait lever le solveur. Le budget qu'A8 mobilise est **nul en valeur**.
Ce qu'il faudrait balayer, ce sont `gcBirthPenMin` / `gcBirthPenMax`, dont le deck
dit lui-même que les bornes sont **arbitraires**.

### 3.3 A11 « une ligne, miroir de `FdemSolver.cpp:3134` » — **instruction dangereuse**

**A11 est FAITE** (cf. [`A11_dt_tangentiel.md`](A11_dt_tangentiel.md)), et le
document le sait — il porte l'encadré « CORRIGÉ ET MESURÉ ». Mais la ligne A11
**subsiste dans le tableau du plan** avec la consigne « miroir de
`FdemSolver.cpp:3134` », **alors que le même document explique quinze lignes plus
haut que ce miroir divise le pas de temps par ≈ 32**.

**Un exécutant qui ne lit que le tableau applique littéralement l'instruction.**
C'est la forme la plus dangereuse d'erreur de rédaction : une consigne fausse
survivant à sa propre réfutation, dans le même fichier.

Le patch effectif fait **~80 lignes sur deux fichiers**, pas une. Effort :
**FAIBLE mais pas trivial**.

---

## 4. QUATRE STATUTS TROP SÉVÈRES — des accusations à retirer

1. **« Le terme visqueux tangentiel −η v_t est ABSENT » — FAUX.**
   `src/Fdem3dSolver.cpp:3683, 3752, 3769` : `ftv = −ctcMu(...)·fn·tanh(vtn/vReg_)·v̂ₜ`.
   Pour |v_t| ≪ v_reg, tanh(x) ≈ x, donc **ftv ≈ −(µ f_n / v_reg)·v_t**,
   c'est-à-dire **exactement −η v_t avec η_eff = µ f_n / v_reg**
   (`contactVreg`, défaut 1e-3).
   **L'énoncé juste est plus intéressant que mon accusation** : les deux branches
   de contact de rockim ont chacune **une moitié** de l'équation 8 d'Imperial —
   la branche **pénalité** a le terme visqueux et **pas de ressort**, la branche
   **potentiel** a le ressort et **pas de terme visqueux**. **Aucune des deux ne
   réunit les deux termes.** Ça reste un écart, mais ce n'est pas une absence.
2. **`jointDeltaC = solidity` « non rachetable » — FAUX.** Deux des trois
   ingrédients sont **sur la même page 12** du manuscrit UCL dont je me sers deux
   lignes plus haut pour racheter `jointDeltaC = guo`. La différence entre les
   deux n'est pas publié/non-publié : ce sont **deux lectures de la borne
   d'intégration ambiguë de l'éq. 10**.
3. **« Transcription, constante magique comprise »** — la charge d'intégrité
   tombe avec la prémisse. Reprendre une forme d'un **dépôt public sous LGPL-3.0**
   n'est pas une faute. Ce qui reste vrai et suffit : **aucune publication
   d'Imperial ne décrit ce couplage** (lot 2b), donc le statut **DIVERGENT** est
   le bon — sans le vocabulaire d'accusation.
4. **Le critère de réussite d'A3 commandait une destruction.** « `grep … rend
   zéro` » **effacerait la mention `ImperialCollegeLondon/solidity-solver-open,
   LGPL-3.0, lu le 2026-08-26` de `tools/verify_suite.py`** — c'est-à-dire la
   provenance la mieux écrite du dépôt.

---

## 5. LES AVANTAGES DE ROCKIM — le compte honnête

Neuf des treize entrées sont **trop flatteuses**. Les trois plus nettes :

* **Entrée 7, garde-fou crack-band** : il ne s'applique qu'aux lois `dpr` et
  `saksala` (`MatLaw.cpp:1305`), or le deck de réplique pose
  **`bulkModel = neohookean`** — **le garde-fou ne se déclenche jamais sur le banc
  d'impact**. Et il porte sur la loi de **volume**, alors qu'en FDEM l'énergie de
  rupture est portée par les **joints**. À **retirer** de la catégorie.
* **Entrée 10, plafond d'impulsion 20 m/s** : il vit dans la branche **pénalité**
  ; le banc de réplique tourne en **`contact = potential`**. **Inerte sur la
  totalité du banc.** C'est exactement la faute que le dépôt s'interdit et
  signale ailleurs de sa propre initiative (`FdemSolver.cpp:829-831`).
* **Entrées 1 et 3, « Imperial ne publie aucune suite » / « aucun code publié ne
  fait ça »** : indéfendable maintenant qu'on sait que leur code **est** public.
  L'énoncé correct est « aucune suite de non-régression n'est **décrite dans les
  publications** d'Imperial ».

**Le compte honnête des « résultats scientifiques » :** **un** de plein droit
(n° 11, la séparation mesurée des leviers pénalité/schéma), **un** publiable après
reformulation (n° 13, l'exposant 0,1707 prédit et confirmé), **une** note de bas
de page (n° 12, la règle de Turon — une coïncidence numérique plausible, non
démontrée), et **une à retirer** (n° 7).

**Un résultat solide et une reformulation publiable, ce n'est pas maigre.** C'est
simplement quatre fois moins que ce que j'annonçais.

---

## 6. LE PLAN REFAIT

*L'ancien plan est à **2/17**. Celui-ci le remplace intégralement.*

**Il n'y a plus de « blocage ».** Il y a des corrections utiles, ordonnées par
rapport gain/risque.

| # | action | effort | critère de réussite, **falsifiable** |
|---|---|---|---|
| **B1** | **Poser `jointPenaltyFactor ≈ 9,6`** au deck de réplique (voie deck d'A2), en documentant que c'est 26,32 corrigé du rapport de longueurs **et que le rapport réel est ≥ 2,4495** | **1 ligne de deck** | le partage d'énergie et la courbe force-pénétration doivent **bouger** ; s'ils ne bougent pas, la pénalité n'était pas le levier qu'on croyait et il faut le dire |
| **B2** | **Retirer du plan A6 tel qu'il était écrit**, et re-spécifier : outil **rigide** via `toolSig_` (qui existe), ou soudure des nœuds de la phase outil | re-spéc. d'abord | avec l'outil rigide, la vitesse de rebond doit se rapprocher de la mesure **sans** que le champ de fissuration dans la roche change |
| **B3** ✅ **FAIT le 2026-08-30** | **Corriger le §1 du lot 4** : chiffre **79** (et non 78 — voir §9 M-3, la différence est la citation de `specs/…/WP7`), commande `grep -rEc`, périmètre étendu à `configs/`, `build_sol.cmd`, `specs/` ; et **porter l'avertissement sur `DOCUMENTATION_rockim.md` §5.4 quinquies**, la pièce la plus lue | faible | ✅ `grep -rEc "Y3D[a-z]*\.c"` sur le périmètre déclaré donne **79**, le nombre imprimé au §1 récrit. Bandeau porté, plus deux sources d'article ajoutées (A3.3) |
| **B4** | **A3 redéfini** (A3.1-A3.3 de la fiche A03) : citations **complètes** (dépôt, licence, commit, fichier, ligne), **trois** statuts au lieu de deux, sources d'article ajoutées là où elles existent. **Ne rien renommer, ne rien effacer.** | moyen | aucune citation `Y3D*.c` nue ne subsiste ; **la mention LGPL de `verify_suite.py` est toujours là** |
| **B5** | **Essai `gcBirth = penalty`** sur l'impact 3D — **après B1**, pour ne pas refaire un run de ~40 h | 1 run | l'injection de la branche normale doit passer **sous 1 % de KE₀** (référence : **3,66 J / 6,9 %**, pas 11,1 J). **Une hausse est un résultat valide** et doit être publiée telle quelle |
| **B6** | **Balayer `gcBirthPenMin` / `gcBirthPenMax`**, dont le deck dit les bornes arbitraires — **et non `gcBirthTau`**, que le deck refuse | faible | le facteur moyen de calage ne colle à aucune borne ; s'il y colle, c'est le clamp qui décide à la place de la physique |
| **B7** ✅ **FAIT le 2026-08-30** | **Corriger le §3** : entrées **7 et 10 RETIRÉES** (toutes deux inertes sur le banc de réplique — `bulkModel = neohookean` pour l'une, `contact = potential` pour l'autre), **1, 3, 6, 8, 9, 12 et 13 requalifiées**, « résultats scientifiques » ramenés de **quatre à un** (+ une reformulation publiable) | faible | ✅ chaque entrée survivante porte une preuve revérifiée par moi contre le dépôt, commande donnée |
| **B8** | **Porter le diagnostic *diffuse ratcheting* sur la percussion 3D** (§9, M-8). Il n'existe aujourd'hui qu'en **2D** et **sur l'essai brésilien** — c'est-à-dire nulle part sur le cas comparé à Imperial | faible (le code 2D existe : `gDfrac_`, `gDmean_` à transposer) | la part de joints au-dessus de D = 0,01 au pic **imprimée** sur `impact_imperial.cfg`, et **son écart entre insertion intrinsèque et adaptative**. C'est la jauge qui manque à B1 et B5 |
| **B9** | **Figer le compte des sources muettes sur une liste NOMINATIVE** (§9, M-11). Aujourd'hui **4, 5, 6, 7 et 8** circulent dans trois documents pour deux affirmations voisines | très faible | **un seul nombre par affirmation**, chaque source nommée dans une table du lot 2c, reprise partout ailleurs par renvoi |

**Ce qu'il ne faut PAS faire**, mis à jour :

* ~~monter la pénalité à 50 E~~ (faux d'un facteur ~5) ;
* ~~implémenter un maillage adaptatif~~ (Imperial n'en a pas) ;
* ~~passer `contactMu` à 0,18 sur St Anne~~ (leur valeur y est 0,6) ;
* ~~« restreindre les joints à la roche » par un filtre~~ — **désassemble l'outil** ;
* ~~balayer `gcBirthTau`~~ — **le deck refuse la clé** ;
* ~~appliquer « miroir de `FdemSolver.cpp:3134` »~~ — **diviserait dt par ≈ 32** ;
* ~~faire tomber `grep Y3D*.c` à zéro~~ — **effacerait une provenance légitime**.

**Cinq interdits que la liste d'origine oubliait** (ajoutés le 2026-08-30, §9 M-24) :

* ~~modifier `hEl_` globalement en faisant B1/A2~~ — il gouverne **aussi** le CFL,
  la borne diffusive, le crack-band et le **seuil de pulvérisation**
  (`src/Fdem3dSolver.cpp:2204, 2218, 359, 2434`). Toute correction de longueur doit
  rester **confinée à la pose des joints**, sinon elle change la physique en silence ;
* ~~recopier `max(potP_, potKt_)` du 2D en 3D~~ — **piège d'unités** : en 3D `potP_`
  est en Pa et `potKt_` en N/m ; la copie diviserait le pas de temps par ≈ 32 ;
* ~~effacer les citations au code public d'Imperial~~ — c'est une source **publique,
  identifiable et sous licence** : il faut l'**ancrer**, pas l'effacer ;
* ~~poser `gcBirthTau` avec `gcBirth = penalty`~~ — **le solveur lève** ;
* ~~SUPPRIMER les joints de l'outil~~ — il faut les **figer** ; les supprimer
  **désassemble l'outil** (§3.1).

**Et un interdit à resserrer** : « si la pulvérisation change les résultats, c'est
un artefact » **ne suit pas** de la source citée. Ce qu'elle établit (IJ206 §2.2
p. 4, reproduction « *without introducing the additional damage model* ») est que
le modèle est **superflu sur ce cas**, pas que tout effet serait un artefact.
L'interdit défendable est : **« ne pas l'activer sur la réplique St Anne »**.

---

## 7. CE QUI SURVIT SANS RETOUCHE

* **Toute la formulation des lots 2 et 3.** Le contre-audit n'y touche pas : elle
  a été établie sur sources primaires et contre-vérifiée à 91 %.
* **La section `loi-joint-penalite` du lot 4**, à 13/17 — la seule passée par un
  cycle complet.
* **A11**, faite, mesurée, verrouillée par deux contrôles, avec sa correction
  d'exagération déjà portée.
* **Le fait que la loi de joint de rockim soit complète et fidèle**, et que le DIF
  soit conforme jusqu'aux points d'application.

---

## 8. LA LEÇON, ET ELLE EST CHIFFRÉE

**91 % de confirmés sur ce qui est passé par audit + contre-audit. 40 % sur ce que
j'ai écrit à la main sans relecture adverse. 12 % sur le plan.**

Le plan — la pièce sur laquelle on allait travailler — était la plus fausse. Ce
n'est pas un hasard : c'est la seule qui ne soit adossée à aucune source, ni
article ni ligne de code. **Rien de ce qui engage une décision ne devrait sortir
d'ici sans être passé par une relecture adverse.**

---

## 9. LES 34 VERDICTS « MOYENNE » — dépouillés un par un le 2026-08-30

*Les §1 à §8 ci-dessus ont traité les 6 critiques et l'essentiel des 32 hautes.
**Les 34 moyennes non confirmées n'avaient jamais été ouvertes.** Elles le sont
ici. Le tri n'est pas fait sur l'étiquette de gravité — plusieurs d'entre elles
touchent une affirmation que je réutiliserais au manuscrit, ce qui les rend plus
lourdes que leur label. **J'ai revérifié moi-même chaque preuve chiffrée contre
le dépôt ; les commandes sont données.** Là où mon comptage diffère de celui du
vérificateur, je le dis.*

### 9.0 Résultat du dépouillement

| disposition | n |
|---|---|
| **CONFIRMÉ, et une affirmation réutilisable change** | **11** |
| CONFIRMÉ, sans effet sur une affirmation réutilisable | 13 |
| **PARTIELLEMENT retenu** (le fond tient, la preuve ou la portée est à corriger) | 7 |
| **REJETÉ** (le vérificateur se trompe, ou son reproche ne porte pas) | 3 |

**Aucun des 34 ne renverse une conclusion de formulation des lots 2 et 3.** Ils
portent tous sur le lot 4, le plan, ou la manière de présenter un avantage de
rockim. C'est cohérent avec le motif du §0 : ce qui est adossé à une source
primaire tient ; ce qui juge le dépôt à la main, non.

---

### 9.1 LES ONZE QUI CHANGENT UNE AFFIRMATION RÉUTILISABLE

#### M-1 (entrée 13 du §3) — « une prédiction confirmée : l'exposant 0,1707 » → **c'est une INFÉRENCE, pas une prédiction. Et le chiffre cité en preuve est faux.**

C'est celle que je vous ai vendue comme publiable. Elle l'est — mais pas sous ce
mot, et pas avec ce nombre.

**Deux défauts, tous deux vérifiés par moi.**

**(a) « Prédiction » n'est pas soutenable.** L'article confirmateur est
*IJRMMS* **206** (2026) 106660 — **la même année** que la dérivation du
2026-08-18. Rien dans le dépôt n'établit qu'il n'était pas déjà paru à cette
date ; la fiche dit seulement qu'il a été **lu** le 2026-08-29 sur un PDF fourni.
Une prédiction demande d'établir l'antériorité, et je ne peux pas.
*(Le [lot 2b](../biblio_insertion/2026-08-29_lot2b_couplage_endommagement_contact.md)
écrit « une source publiée un an plus tard » : c'est vrai de l'article de **2025**,
pas de celui de 2026.)*

**(b) Le raccord bas imprimé dans `include/rockim/YangDif.hpp` est FAUX.**
L'en-tête écrit « 0,1707 raccorde EXACTEMENT les deux bornes (**1,0031** en bas,
1,8500 en haut) ». Recalculé :

```
0,95 + 0,41·(5e-6)^0,1707 = 1,001039   →  1,0010, et non 1,0031
0,95 + 0,41·(1e2)^0,1707  = 1,849887   →  1,8500 ✓
```

Le raccord **haut** est exact ; le **bas** vaut 1,0010. La valeur correcte était
d'ailleurs déjà écrite au [lot 2b](../biblio_insertion/2026-08-29_lot2b_couplage_endommagement_contact.md)
(l. 217) : **1,0010**. L'artefact cité comme preuve portait donc un chiffre erroné.
**Corrigé dans l'en-tête le 2026-08-30**, avec la note de correction datée.

> **La formulation qui résiste à un rapporteur** — à utiliser telle quelle :
>
> « La loi de DIF en traction publiée en 2025 avec l'exposant 0,07 ne se raccorde
> à **aucune** de ses deux bornes : elle vaut 1,1245 au lieu de 1 en ε̇ = 5·10⁻⁶ /s,
> et 1,5160 au lieu du plateau 1,85 en ε̇ = 10² /s, soit **22 % de discontinuité**.
> L'exposant qui raccorde **simultanément** les deux vaut **0,1707** (1,0010 en
> bas, 1,8500 en haut) — valeur relevée **indépendamment** sur la courbe tracée de
> leur propre figure 2(b). L'article de 2026 imprime **0,17**. »
>
> C'est une démonstration de **cohérence interne**, testable et falsifiable. Elle
> est plus forte ainsi que sous le mot « prédiction », qui invite à vérifier une
> date que je ne peux pas défendre.

Et **le vrai apport n'est pas l'exposant** : c'est la conséquence physique
mesurée. Le saut de 22 % en ε̇ = 10² /s est un **attracteur** en insertion
extrinsèque — un joint qui franchit le seuil voit sa résistance bondir et cesse
de s'insérer, si bien que la population insérée s'empile juste sous 10². Mesure
du 2026-08-18 : **médiane 99,36 /s** (max 99,9988) avec l'exposant littéral,
contre **40,22 /s** avec 0,1707. Cela, aucune publication ne le décrit, et c'est
verrouillé par `dif_yang_litteral_2d` et `dif_yang_fig2_2d`.

#### M-2 (entrée 1 du §3) — « 98 contrôles, pas 42 » → **97, et les tiers ne s'additionnent pas.**

Recompté par dépouillement AST de la liste `TESTS` de `tools/verify_suite.py` à
HEAD :

```
44 fast + 45 full + 8 all = 97 tests, 217 tuples de checks
```

Trois erreurs dans une seule phrase :

1. **98 n'est le total de rien.** Aucun décompte du fichier ne le donne.
2. **La ventilation citée était périmée ET incohérente** : 42+45+8 = **95**, pas
   98. Le 42/45/8 est l'état d'**avant** le commit `53a4371` (chantier A11), qui a
   ajouté `dtbudget_tangential_defaut_3d` et `dtbudget_tangential_on_3d`.
3. **Les tiers s'EMBOÎTENT, ils ne s'additionnent pas.**
   `tools/verify_suite.py:894` : `TIERS = {"fast": ["fast"], "full": ["fast","full"],
   "all": ["fast","full","all"]}`. Donc `--tier full` rejoue 44+45 = **89**, et
   `--tier all` = **97**. L'invocation par défaut (`:962`, `default="fast"`)
   n'exécute que **44** contrôles.

**Corollaire, et il retourne le reproche fait au brief** : le « 42 » du brief
n'était pas une erreur — c'était **exactement le tier par défaut de l'époque**.

> **Formulation à utiliser** : « 97 contrôles de non-régression, dont **44 au tier
> rapide par défaut**, 89 au tier `full`, 97 au tier `all` ; **217 assertions**
> chiffrées au total. »

Citer « 98 » au manuscrit, c'est publier un nombre que personne ne pourra
reproduire.

#### M-3 (§1 du lot 4) — « 117 attributions » → **le nombre est faux sous TOUTES les conventions de comptage.**

Recompté par moi, cinq façons, sur le périmètre étendu
(`src include tools bench_impact/configs specs configs build_sol.cmd`) :

| ce qu'on compte | commande | n |
|---|---|---|
| lignes citant un fichier `Y3D*.c` | `grep -rEc "Y3D[a-z]*\.c"` | **79** |
| … dont celles portant un numéro de ligne `l. NNNN` | `grep -rEc "Y3D[a-z]*\.c[^ ]* l\. *[0-9]"` | **72** |
| occurrences de `Y3D*.c` | `grep -rhoE "Y3D[a-z]*\.c" \| wc -l` | 79 |
| lignes citant `Y3D*.c` **ou** « solidity » | `grep -rEc "Y3D[a-z]*\.c\|[Ss]olidity"` | **182** |
| occurrences idem | `grep -rhoE ... \| wc -l` | 208 |

**117 n'est aucun de ces nombres.** Et le tableau du §1 du lot 4 se contredisait
lui-même : ses propres lignes somment à 172 (périmètre restreint), pas 117.

**Le vérificateur annonçait 78** pour la première ligne ; **je trouve 79** — la
différence est la citation de `specs/005-impact-insert-yang/WP7_couplage_contact.md`,
hors de son périmètre. Sa ventilation par fichier est exacte pour le reste.
*(C'est le seul écart de comptage entre lui et moi sur les 34.)*

#### M-4 (§1 du lot 4) — la commande de preuve imprimée **ne rend rien**.

`grep -rc "Y3D[a-z]*\.c|[Ss]olidity" …` : **exit 1, aucune sortie**. Sans `-E`
(ni `\|`), l'alternance `|` est un caractère **littéral** en expression régulière
de base : le motif cherche la chaîne « Y3D….c|Solidity », qui n'existe nulle part.

Ce n'est pas un détail de forme. **Le lecteur qui recopie la commande — un
rapporteur, typiquement — conclut que le résultat principal du document n'existe
pas.** Il faut **`grep -rEc`**.

#### M-5 (§2.6 du lot 4) — « ce n'est pas une inspiration lointaine, c'est une **transcription** » → **le mot est faux dans l'autre sens.**

Ce que le dépôt établit (CORRECTION 2, `BILAN_interference_2026-08-29.md:245-260`) :
dans le code public, le facteur d'endommagement est câblé à zéro, donc `d_fact = 1`,
donc le couplage `penalty *= d_fact` et `mu = mud*d_fact` **ne s'exécute jamais**.
Côté rockim, `cplDf`/`ctcMu` sont pilotés par `el_[e].bdD`, l'endommagement
volumique **propre à rockim**, exigé à `src/Fdem3dSolver.cpp:648-652`.

**rockim n'a donc pas transcrit un COMPORTEMENT : il a repris la FORME d'une
branche qui ne s'exécute jamais dans le code lu, et l'a branchée sur un moteur
d'endommagement qui, lui, est celui de l'article de 2026 (éq. 3-4).** Ce n'est ni
une réplication, ni une transcription — c'est **l'activation d'une branche
morte**. C'est un geste défendable et intéressant ; ce n'est pas ce que le mot
« transcription » annonce, et l'« aggravation » revendiquée ne suit pas de la
preuve invoquée.

#### M-6 (§3, entrée 12 du lot 4) — le bilan d'énergie « refuse de continuer s'il dérive » → **l'ARRÊT est opt-in, et le deck de réplique ne l'arme pas.**

Vérifié : `src/Fdem3dSolver.cpp:2316` `eAbortPct_ = cfg_.getd("budgetAbortPct", 0.0);`
puis `:2322` `if (eAbortPct_ <= 0.0 || eAbort_) return;`. Sur les **22 decks** de
`bench_impact/configs`, **quatre** l'arment (`impact_kuru9`, `impact_kuru11`,
`impact_p2_nombres`, `impact_p2_facies`, tous à 2 %). **`impact_imperial.cfg` ne le
pose PAS** — et c'est précisément le deck que le §5 et les actions A1/A2 désignent
comme le deck de réplique St Anne.

> **La distinction juste, et elle reste un avantage** : la **MESURE** du résidu est
> **inconditionnelle** — `src/Fdem3dSolver.cpp:4519-4522` imprime toujours
> « residu : … J (… % de l'echelle) [OK|CHECK] ». C'est l'**INTERRUPTION** qui est
> opt-in. Formulé ainsi, l'argument tient devant un relecteur ; formulé comme
> « refuse de continuer s'il dérive », il tombe dès qu'on ouvre `impact_imperial.cfg`.

**Correctif à coût nul, à faire** : poser `budgetAbortPct` sur les decks de
réplique. *(Non fait : c'est une modification de deck, elle change ce que le run
produit — à décider par vous.)*

#### M-7 (§3, entrée 12 du lot 4) — « postes ventilés » → **le poste GRAVITAIRE n'existe pas, et c'est celui qu'ARMA distingue.**

Six des sept postes attendus sont là et vérifiés. Le septième, non :
`src/Fdem3dSolver.cpp:4038` `bodyForces()` applique la pesanteur **sans aucun
compteur de travail**, et `grep gravWork include/rockim/Fdem3dSolver.hpp` → **0**.
Or `gravity = 9.81` est posé dans `impact_imperial.cfg`, `impact_kuru9.cfg`,
`impact_kuru11.cfg` et huit autres decks. Et ARMA 24-0952 pose explicitement
l'énergie **potentielle gravitaire** dans ses éq. 3-7.

**Magnitude honnête** : sur 600 µs d'impact les déplacements sont en micromètres,
le travail de la pesanteur est de l'ordre de **10⁻⁴ J** contre ~49 J injectés —
**invisible**. Le défaut est **structurel, pas numérique**. Il devient réel dès
qu'un run est long ou quasi statique.

**Second trou, plus gênant** : `src/Fdem3dSolver.cpp:4063`
`brushWork_ += bw; // POSTE SEPARE, jamais dans sumW`, imprimé « hors bilan B4 ».
Dès que le tri anti-gravité des fragments est armé, **son travail tombe
entièrement dans le résidu** — et `budgetAbortPct` peut couper un run **sain** sur
un artefact purement numérique.

#### M-8 (§3, entrée 9 du lot 4) — le diagnostic « diffuse ratcheting » → **il ne s'imprime que sur l'essai BRÉSILIEN, et en 2D.**

Vérifié : `src/FdemSolver.cpp:6840-6843` est enfermé dans
`if (scen_ == Scenario::BRAZILIAN)` (ouvert `:6770`), et
`grep "gDfrac_\|diffuse ratcheting\|sub-critical damage" src/Fdem3dSolver.cpp`
→ **0 occurrence**.

**Le périmètre est l'inverse de ce qu'on croit en lisant l'entrée** : la
pathologie de la pénalité intrinsèque est diagnostiquée sur un essai de
**calibration 2D**, et **pas du tout** sur la percussion 3D — c'est-à-dire pas sur
le cas comparé à Imperial, celui-là même où le §5 mesure une injection d'énergie
en insertion intrinsèque.

> **Ce qu'il faut en tirer est une ACTION que le plan n'avait pas** : porter ce
> diagnostic sur le scénario de percussion 3D. Cela donnerait à A1 et A2 la jauge
> qui leur manque. Ajoutée au plan refait (§6) comme **B8**.

#### M-9 (§3, entrée 6 du lot 4) — `tools/yan_point.cpp` « la loi expédiée est testée » → **à moitié seulement, et rien n'est automatisé.**

Trois écarts, vérifiés :

1. **« La loi expédiée » n'est vraie qu'à moitié.** L'en-tête dit bien que `f(D)`
   et son intégrale sont les fonctions **partagées** (`rockim/YanSoftening.hpp`) —
   mais la ligne suivante dit que le pilote **« reproduces »** / **« mirrors »** la
   mise à jour de traction de `jointForces()`. **Cette moitié-là est une
   ré-implémentation**, susceptible de diverger du solveur sans que rien ne le
   signale.
2. **Il ne « compare » rien** : il imprime `GfI_total`, `GfI_elastic`,
   `GfI_fracture`, `GfI_target` — **aucun ratio, aucun verdict PASS/FAIL**.
3. **Il n'est pas automatisé** : `grep -c yan_point CMakeLists.txt build*.cmd`
   → **0 partout** (31 scripts de build vérifiés). Et le contrôle `yan_integral`
   de la suite (`tools/verify_suite.py:140-144`) **ne l'appelle pas** : il
   verrouille la valeur imprimée par le **SOLVEUR**, avec le commentaire explicite
   « stdout n'imprime que 6 décimales : la vérification à 1e-12 vit dans
   `tools/yan_point.cpp` ».

> **Formulation défendable** : « un pilote au point matériel qui réutilise la
> fonction d'adoucissement **expédiée** et trace σ(δ) et τ(δ) pour les figures ; la
> valeur de son intégrale est verrouillée dans la suite par le solveur lui-même. »
> L'avantage existe — aucune publication d'Imperial ne décrit un tel pilote — mais
> il est **manuel, non compilé par le build, sans verdict**, et la moitié du chemin
> testé est une copie.

#### M-10 (§3, entrée 8 du lot 4) — « la parité 2D/3D posée en principe et **instrumentée** » → **couverte des deux côtés, mais jamais COMPARÉE ; et quatre ruptures sont ouvertes.**

Le principe est bien nommé (`include/rockim/YangDif.hpp:8-13`, en toutes lettres :
une divergence 2D/3D sur les bornes en dur « serait MUETTE et fausserait toute
comparaison dimensionnelle »). Et l'instrumentation existe : **13 paires** de
contrôles homonymes `_2d`/`_3d`.

**Mais aucun contrôle de la suite ne COMPARE une grandeur 2D à son homologue 3D** :
les 13 paires portent des références **indépendantes**. Chacune verrouille son
solveur ; **aucune ne teste l'égalité**. Et la parité n'est pas acquise — le
document l'établit lui-même : impression de la pénalité absente en 2D, `potXi_`
en 2D seulement, `potKt_` hors budget de pas de temps en 3D jusqu'à A11. **J'en
ajoute une quatrième, vérifiée** : le diagnostic de l'entrée 9 (ci-dessus M-8) n'a
**aucun** équivalent 3D. Et une cinquième que je relève ici : **`groupBond`
n'existe qu'en 3D** (`grep -c "groupBond\|gbond_" src/FdemSolver.cpp` → **0**).

> **La formulation honnête est PLUS forte que celle qu'elle remplace** : « la
> parité 2D/3D est une **règle de constitution** ; deux lois (`YanSoftening.hpp`,
> `YangDif.hpp`) sont **partagées** entre les solveurs pour la garantir par
> construction ; 13 capacités sont contrôlées des deux côtés — et **le dépôt tient
> le registre de ses propres ruptures de parité**, dont **cinq** sont ouvertes. »
>
> **C'est le registre des ruptures qui est l'avantage, pas la parité.**

#### M-11 (§2.6 et plan) — « Imperial ne publie aucune règle de paire » → **le décompte de sources n'est stable nulle part.**

Trois nombres circulent pour la **même** question :

| où | ce qui est écrit |
|---|---|
| lot 2c §6 (l. 392) | « **quatre** sources muettes » |
| lot 4 §1.1 (l. 88) | « sur **cinq** sources » |
| lot 4 §2.6 (l. 193) | « (**six** sources) » |

*(Et pour k_t, le lot 2c §4 dit « **Sept** sources muettes », le lot 4 « huit ».)*

**Un décompte qui varie de 4 à 6 pour la même affirmation ne peut pas être cité.**

Et le statut « **non sourçable** » est à retirer : le lot 2c (l. 318) écrit
lui-même que « ici tous les corps sont du **même matériau** : **la question ne se
pose pas**, ce qui explique peut-être qu'elle ne soit jamais traitée ». Ce n'est
pas un silence coupable, c'est un **cas qui ne se présente pas chez eux**.

> **Formulation à utiliser** : « Aucune des sources dépouillées ne traite la
> combinaison du frottement pour une **paire de matériaux différents** — et pour
> cause : dans tous leurs essais publiés les corps en contact sont du **même**
> matériau, la question ne se pose donc pas. Le **minimum** retenu par rockim est
> un choix propre, à assumer comme tel et à justifier physiquement (le plus faible
> gouverne le glissement), sans se réclamer d'eux **ni** leur reprocher un
> silence. » **Le nombre de sources ne doit pas apparaître** tant qu'il n'est pas
> figé sur une liste nominative.

---

### 9.2 LES TREIZE CONFIRMÉS QUI NE CHANGENT AUCUNE AFFIRMATION RÉUTILISABLE

Enregistrés, non repris dans le corps des fiches — ils portent sur des
estimations d'effort et des étiquettes de statut internes au plan, lequel est
déjà **remplacé** par le §6 ci-dessus.

| # | ce qui était écrit | verdict | ce qui est retenu |
|---|---|---|---|
| M-12 | §4 : « quel que soit le facteur » l'ouverture au pic ne peut pas coïncider | trop sévère | **Faux : `pj = pf·E/h` est un produit** (`src/Fdem3dSolver.cpp:1521`, site **unique**) — un facteur scalaire sur `pf` est **strictement** équivalent à un facteur sur `h`. Mesuré : `pf = 26,32 × 0,3648 = 9,60` fait coïncider l'ouverture au pic **en moyenne** sur 36 682 joints. **Le vrai argument est ailleurs, et il est plus fort** : la **dispersion joint-par-joint** (0,240 à 0,452, soit un `pf` équivalent de **6,3 à 11,9** sur le même maillage) — qu'aucun scalaire ne corrige. *(Déjà porté au §1.1.)* |
| M-13 | §2.2 : `pf` noté « DIVERGENT, blocage non » et `h` noté « blocage OUI », **sur le même nombre** | blocage mal jugé | Incohérence interne : la même erreur de raideur est bloquante écrite du côté de la longueur, non bloquante du côté du facteur. **C'est « blocage = non » qui est la bonne des deux.** |
| M-14 | decks : « 41 occurrences » | preuve fausse | Ni 39 (lignes) ni 42 (occurrences). **Seule ligne de la table jamais ventilée par fichier — donc la seule jamais vérifiée.** |
| M-15 | `gcBirth = penalty` « NON RACHETABLE » | trop sévère | **PARTIELLEMENT rachetable** : le problème **et** une rampe sont publiés (manuscrit UCL **p. 17, éq. 18**). Seul le **ré-échelonnement de pénalité** est propre à rockim. *(Porté sur `DOCUMENTATION_rockim.md` §5.4 quinquies le 2026-08-30.)* |
| M-16 | A3 « effort faible, 1 j » | irréaliste | 172 à 182 lignes, 3 statuts à distinguer, 2 sources de 4 900 et 7 100 lignes, plus la table de rachat à refaire. **La fiche A03 le dit dans ses propres mots.** |
| M-17 | A1 « l'essai n° 1, effort faible » | irréaliste | Trois obstacles matériels : la ligne de base A/B a **disparu** (scratchpad volatile) ; A1 n'est **pas** un essai à une variable (il faut poser **aussi** `jointDeath = damage`, sinon `fDeath ≈ 0` et `gcBirth = penalty` dégénère) ; le coût est d'une **nuit de calcul**. |
| M-18 | A2 « effort faible, 1 j » | irréaliste | 1 j tient pour le **seul patch 3D opt-in**. 1,5 à 2 j avec la parité 2D, la fiche et les mesures. **Et le périmètre n'est pas défini** : `hEl_` gouverne aussi CFL, borne diffusive, crack-band et **seuil de pulvérisation** — le toucher globalement change tout cela en silence. |
| M-19 | A8 « balayer `gcBirthTau`, 0,5 j, gain mesurable » | trop flatteur | **`gcBirthTau` est INERTE sous `gcBirth = penalty`** — et le solveur **refuse** de poser les deux clés ensemble. Le deck de réplique pose déjà `gcBirth = penalty`. **A8 ne peut porter que sur la configuration que le plan cherche à abandonner** : soit c'est un repli conditionnel à l'échec d'A1 — et il faut l'écrire — soit il sort de la liste. *(Sorti, §6.)* |
| M-20 | A12 « banc analytique, effort faible » | irréaliste | Aucun scénario existant ne lance un corps **libre** à vitesse horizontale sur un plan fixe ; la gravité existe, **pas** le mécanisme de vitesse initiale de corps ; et `verify_suite.py` n'a **aucun** extracteur de distance d'arrêt. |
| M-21 | A12 : critère « distance d'arrêt à 1 % de v²/2µg » | trop flatteur | **Le critère est écrit pour V1** (rectangle maillé) alors que la voie recommandée est **V2** (point matériel) — où la distance vaut v²/2µg **par construction** : le critère redevient **tautologique**. Et trois chiffres circulent pour la même action (« MOYEN à ÉLEVÉ », « 1-2 j », « 0,5-1 j + refactor »). |
| M-22 | §2.1 : mort de joint « à trancher », porté comme blocage potentiel | blocage mal jugé | `bench_impact/configs/impact_imperial.cfg:187` pose **déjà** `jointDeath = damage`. La capacité existe, elle est documentée, le deck de réplique l'arme. **Aucun résultat publié n'est hors d'atteinte de ce fait → à déclasser en « non ».** |
| M-23 | §2.3 : détection paresseuse notée **PARTIEL**, effort « — » | trop sévère | **Un PARTIEL sans manque nommé est un statut vide.** La preuve citée supporte **PRÉSENT** : le jeu actif est même **en-deçà** de l'extérieur sous `gcActivation = adaptive` (4 % des faces, ×2,32 bit-identique). |
| M-24 | liste « ce qu'il ne faut PAS faire » : six interdits | trop sévère | Quatre solidement étayés ; le cinquième repose sur la prémisse désavouée ; le sixième est **trop large** (la source établit que le modèle est **superflu** sur ce cas, pas que tout effet serait un artefact). **Et la liste est trop COURTE** — cinq interdits manquent, ajoutés au §6. |

---

### 9.3 LES SEPT PARTIELLEMENT RETENUS

| # | le reproche | ce que je retiens, ce que je rejette |
|---|---|---|
| M-25 | « le choix de convention de signe n'est documenté nulle part » → **preuve fausse** | **Retenu.** `src/FdemSolver.cpp:4704-4709` porte bien le commentaire « ressort tangentiel a **HISTOIRE** : Ft += -kt (vrel.t) dt, plafonne au cap de Coulomb ». **Mais je retiens surtout ce que le vérificateur ajoute** : le bloc **3D** (`src/Fdem3dSolver.cpp`) porte l'attribution **sans** cette phrase. **C'est une asymétrie de documentation 2D/3D réelle** — sixième rupture de parité (cf. M-10). |
| M-26 | « rockim prend −k_t δ_t, pas la convention de 2009 » → **trop sévère** | **Retenu.** Le lot 2c écrit lui-même que les deux conventions sont **physiquement équivalentes**, l'une orientant simplement la force contre le glissement. Ce n'est pas un PARTIEL. **Et le vrai écart est ailleurs, non relevé** : la **remise à zéro** de δ_t — `if (H.step < stepCount_ - 1) H.Ft.setZero();` efface tout l'historique de glissement dès qu'une paire saute **un seul pas**, et le lot 2c §4 dit précisément que cette remise à zéro « n'est pas décrite » chez eux. **C'est cela qu'il fallait marquer PARTIEL.** |
| M-27 | « aucun contrôle sur le chemin tangentiel de contact » → **preuve fausse** | **Retenu.** `tools/verify_suite.py:83` définit `"gcfric": r"dont frottement …"`, et **trois** contrôles la verrouillent (`:716`, `:737`, `:742`). **Rejeté en partie** : deux des trois noms que le lot 4 citait (`jointResidualMu`, `jointFrictionScaled`) ne sont pas des tests mais des **clés de deck** — l'erreur est là, pas dans l'absence de couverture. **Ce qui reste vrai** : aucun banc **analytique** (à solution fermée) du frottement de contact — et c'est ce que A12/B-A12 vise. |
| M-28 | A6 « aucun filtre par phase n'existe » → **trop sévère** | **Retenu, et il change l'effort.** `src/Fdem3dSolver.cpp:1159-1180` : la clé `groupBond.<A>.<B> = joints` **existe déjà**, avec sa table `gbond_`, sa validation et sa bannière ; et `:1436-1452` porte littéralement le branchement « groupes différents et non liés → pas de joint, la face part au contact général ». **La recherche par mots-clés avait échoué parce qu'elle cherchait cinq noms qui n'existent pas.** La plomberie d'A6 est une **clé sœur** de `groupBond`, pas une écriture ex nihilo. **Mais la sémantique — FIGER plutôt que SUPPRIMER — reste ce qui coûte** (cf. §3.1 : supprimer désassemble l'outil). **Le plan avait jugé les deux à l'envers.** |
| M-29 | « `potKt_` n'entre pas dans le budget de pas de temps ; `kContact` n'existe pas dans ce fichier » → **preuve fausse** | **Retenu sur la forme, rejeté sur le fond.** Le constat **était vrai à l'écriture** (22:47) ; A11 a réécrit le fichier à 23:10. **Ce qui reste à corriger** : la **ligne du tableau affiche toujours ABSENT** alors que l'encadré juste au-dessous dit le contraire — **le tableau et son encadré se contredisent sur la même page**. Corrigé. *(Le vérificateur note honnêtement qu'il ne peut pas trancher l'état de 22:47, faute de copie antérieure.)* |
| M-30 | « la marge de stabilité n'est pas celle qu'on croit **dès que le frottement travaille** » → **trop sévère** | **Retenu, et l'argument physique est bon.** `src/Fdem3dSolver.cpp:3384-3386` écrête `Ft` au cap de Coulomb : **quand le frottement TRAVAILLE (glissement), le ressort est saturé et n'apporte plus AUCUNE raideur** — il ne raidit qu'en **adhérence**. La phrase désigne donc précisément le régime où l'omission **ne peut pas** mordre. Les mesures d'A11 sont par ailleurs **reproduites à l'identique** par le vérificateur (−2,21 % et −4,66 %, dans une marge de 6,7×). **Le tableau du §6 portait encore la formule exagérée sans chiffre** : corrigé. |
| M-31 | « 518 pas contre ~10 chez eux » → **preuve fausse** | **Retenu.** Le `dt` de 1,93e-9 s appartenait à un **autre run** (gradient St Anne). Sur les seuls `dt` 3D mesurés et publiés par le dépôt (**1,30e-8** et **1,93e-8** s, fiche A11 §5-6) : **77 et 52 pas**, soit **5 à 8×**, pas ~50×. **Et la comparaison elle-même est boiteuse** : `gcBirthTau` est une constante de décroissance **exponentielle**, pas la longueur d'une rampe linéaire. *(Porté sur `DOCUMENTATION_rockim.md` le 2026-08-30.)* |

---

### 9.4 LES TROIS REJETÉS

| # | le reproche | pourquoi je ne le retiens pas |
|---|---|---|
| M-32 | « `gcBirth = penalty` : la ligne mérite PARTIELLEMENT RACHETABLE » — reproche **annexe** : « la ligne renvoie au §7 pour le correctif, or le §7 est *ce qui reste à auditer* » | Le fond est retenu (M-15). **Ce reproche-ci porte sur un renvoi interne** qui disparaît avec la réécriture du §1 : il n'y a rien à corriger séparément. |
| M-33 | « le geste 2 (renommer les clés) casse la traçabilité des repères historiques » | **Le geste 2 est ANNULÉ**, pas re-chiffré : la fiche A03 §5 établit que le nom `solidity` est **exact** et ne doit pas être renommé. Le reproche vise une action qui n'existe plus. |
| M-34 | « le vérificateur ne peut pas exclure qu'un deck SHEAR détourné y arrive à moindres frais » (banc de frottement) | **Ce n'est pas un verdict, c'est une réserve de méthode** — et elle est déjà la conclusion de la fiche A12 : « à trancher par essai plutôt que par lecture ». Rien à corriger. |

---

### 9.5 CE QUE CE DÉPOUILLEMENT AJOUTE AU PLAN

Deux actions **nouvelles**, issues des moyennes et absentes de tous les plans
précédents :

* **B8 — porter le diagnostic « diffuse ratcheting » sur le scénario de percussion
  3D** (M-8). Aujourd'hui il n'existe qu'en 2D et sur l'essai brésilien, c'est-à-dire
  nulle part sur le cas comparé à Imperial. **C'est la jauge qui manque à A1 et A2** :
  sans elle, on mesure l'injection d'énergie sans pouvoir dire quelle fraction des
  joints ratchette. Effort : faible (le code 2D existe, `gDfrac_`/`gDmean_` à
  transposer). Critère : la part de joints au-dessus de D = 0,01 au pic **imprimée**
  sur `impact_imperial.cfg`, et son écart entre insertion intrinsèque et adaptative.

* **B9 — figer le compte des sources muettes sur une liste NOMINATIVE** (M-11).
  Tant que 4, 5, 6, 7 et 8 circulent pour deux affirmations voisines, aucune ne peut
  aller au manuscrit. Effort : très faible (une table dans le lot 2c, reprise
  partout). Critère : **un seul nombre par affirmation**, chaque source nommée.

Et **deux corrections de deck** que je n'applique pas sans votre accord, parce
qu'elles changent ce que les runs produisent :

* poser `budgetAbortPct` sur `impact_imperial.cfg` et `impact_imperial_coulomb.cfg`
  (M-6) — sinon l'avantage revendiqué n'est pas armé là où on le revendique ;
* décider si `brushWork_` doit entrer dans `sumW` (M-7) — en l'état, le tri
  anti-gravité des fragments verse son travail dans le résidu, et un garde-fou
  d'énergie peut couper un run sain.
