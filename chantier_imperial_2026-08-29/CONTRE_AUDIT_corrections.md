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

**6 critiques, 32 hautes.** Le dépouillement des sources, lui, était à **91 %**
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
| **B3** | **Corriger le §1 du lot 4** : chiffre 78, commande `grep -rEc`, périmètre étendu à `configs/`, `build_sol.cmd`, specs WP6/WP7 ; et **porter l'avertissement sur `DOCUMENTATION_rockim.md` §5.4 quinquies**, la pièce la plus lue | faible | `grep -rEc` donne le même nombre que le document, dans le périmètre déclaré |
| **B4** | **A3 redéfini** (A3.1-A3.3 de la fiche A03) : citations **complètes** (dépôt, licence, commit, fichier, ligne), **trois** statuts au lieu de deux, sources d'article ajoutées là où elles existent. **Ne rien renommer, ne rien effacer.** | moyen | aucune citation `Y3D*.c` nue ne subsiste ; **la mention LGPL de `verify_suite.py` est toujours là** |
| **B5** | **Essai `gcBirth = penalty`** sur l'impact 3D — **après B1**, pour ne pas refaire un run de ~40 h | 1 run | l'injection de la branche normale doit passer **sous 1 % de KE₀** (référence : **3,66 J / 6,9 %**, pas 11,1 J). **Une hausse est un résultat valide** et doit être publiée telle quelle |
| **B6** | **Balayer `gcBirthPenMin` / `gcBirthPenMax`**, dont le deck dit les bornes arbitraires — **et non `gcBirthTau`**, que le deck refuse | faible | le facteur moyen de calage ne colle à aucune borne ; s'il y colle, c'est le clamp qui décide à la place de la physique |
| **B7** | **Corriger le §3** : retirer l'entrée 7, requalifier 1, 3 et 10, ramener les « résultats scientifiques » à un plus une reformulation | faible | chaque entrée survivante porte une preuve qui résiste à un rapporteur |

**Ce qu'il ne faut PAS faire**, mis à jour :

* ~~monter la pénalité à 50 E~~ (faux d'un facteur ~5) ;
* ~~implémenter un maillage adaptatif~~ (Imperial n'en a pas) ;
* ~~passer `contactMu` à 0,18 sur St Anne~~ (leur valeur y est 0,6) ;
* ~~« restreindre les joints à la roche » par un filtre~~ — **désassemble l'outil** ;
* ~~balayer `gcBirthTau`~~ — **le deck refuse la clé** ;
* ~~appliquer « miroir de `FdemSolver.cpp:3134` »~~ — **diviserait dt par ≈ 32** ;
* ~~faire tomber `grep Y3D*.c` à zéro~~ — **effacerait une provenance légitime**.

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
