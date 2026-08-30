# LOT 4 — Bilan de rockim contre l'état de l'art Imperial

> # ⚠️ CE DOCUMENT EST CORRIGÉ PAR UN CONTRE-AUDIT
> **129 verdicts, 52 confirmés — 40 %.** Six critiques, trente-deux hautes,
> trente-quatre moyennes. Les corrections, et **le plan refait**, sont dans
> [`chantier_imperial_2026-08-29/CONTRE_AUDIT_corrections.md`](../chantier/CONTRE_AUDIT_corrections.md)
> — dont le **§9** dépouille les 34 moyennes une par une (2026-08-30).
> **En cas de contradiction, c'est ce document-là qui fait foi.**
>
> **État au 2026-08-30 : les §1 et §3 ont été RÉCRITS ici même**, leurs versions
> d'origine conservées intégralement aux §1bis et §3bis. Ce qui a changé :
> * **§1** — la prémisse (`solidity-solver-open` **EST** le code d'Imperial), le
>   nombre (**79 lignes** citant un `Y3D*.c`, dont 72 avec numéro de ligne — jamais
>   117), la commande (il faut **`grep -rEc`** ; sans `-E` elle ne rend **rien**),
>   le périmètre (+`specs/`, `configs/`, `build_sol.cmd`) et **sept des huit
>   verdicts** de la table de rachat ;
> * **§3** — **neuf des treize avantages** étaient trop flatteurs : deux **retirés**
>   (n° 7 et 10, inertes sur le banc de réplique), quatre **affaiblis ou
>   reformulés**, et le compte des « résultats scientifiques » passe de **quatre à
>   un** (+ un publiable après reformulation) ;
> * **§6** — le plan contient **trois actions dangereuses ou impossibles** (A6
>   désassemble l'outil, A8 est refusée par le deck, A11 portait une consigne qui
>   aurait divisé le pas de temps par 32) : il est **remplacé** par le §6 du
>   contre-audit ;
> * **les deux « blocages » tombent** (§1.1 du contre-audit).
>
> Les §0, §2, §4, §5, §6 et §7 ci-dessous sont **conservés tels quels** — on ne
> réécrit pas l'historique — et restent sous l'autorité du contre-audit.

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
Le ~~§9~~ **§7** dit ce qui reste à faire *(ce document n'a que sept sections ;
renvoi corrigé le 2026-08-30)*.

---

## 1. LE RÉSULTAT PRINCIPAL — **récrit le 2026-08-30 après contre-audit**

> **Ce §1 a été entièrement réécrit.** Sa version d'origine du 2026-08-29 est
> conservée **intégralement** au [§1bis](#1bis-la-version-dorigine-du-2026-08-29-conservée)
> — on ne réécrit pas l'historique, mais on ne laisse pas non plus un lecteur
> pressé tomber d'abord sur un texte faux. **Ce qui a changé** : la prémisse
> (§1.0), le nombre, la commande et le périmètre (§1.1), la nature même du défaut
> (§1.2), sept des huit verdicts de la table de rachat (§1.3) et l'estimation
> d'effort (§1.4).

### 1.0 La prémisse d'origine était fausse, et il faut le dire avant tout le reste

La version du 2026-08-29 annonçait comme résultat principal que « rockim porte
117 attributions **à un code qui n'est pas celui d'Imperial** ».

**C'est faux.** `solidity-solver-open` **EST** le code d'Imperial College London :
dépôt public `ImperialCollegeLondon/solidity-solver-open`, **LGPL-3.0**, C,
17 000 lignes, format `.Y3D` — la lignée Munjiza de la thèse de Guo et des
articles de Yang *et al.* Cloné et lu le **2026-08-26**, provenance documentée en
**quatre endroits** du dépôt (`CR_solidity_2026-08-27.md:19`,
`BILAN_replique_solidity_2026-08-27.md:11`, `DOCUMENTATION_rockim.md`
§5.4 quinquies, `tools/verify_suite.py:443`). L'histoire de l'erreur — une
sur-correction du 29/08 au soir, reprise par le brief de mission — est retracée
dans [`../chantier_imperial_2026-08-29/A03_resourcer_attributions.md`](../chantier/A03_resourcer_attributions.md) §3.

**Ce qui reste vrai, et c'est l'avertissement utile** — plus étroit, mais solide :

> Le code public **n'est pas la version qui a produit l'article de 2026**. Son
> facteur d'endommagement d'élément y est câblé à zéro (`Y3Dfd.c` l. 749-751,
> `df = R0`) et son DIF est neutre (`dpeftdif = R1`), alors que l'article publie
> les équations (3)-(4) d'un modèle d'endommagement. **Lire une FORME dans ce code
> et en conclure une implémentation de ce que décrit l'article de 2026 reste une
> faute** — non pas parce que ce serait le code de quelqu'un d'autre (c'est bien le
> leur, même lignée, mêmes auteurs), mais parce que **ce n'est pas la version dont
> l'article parle**. La lecture la plus simple : *version ouverte en retard sur la
> version interne*, ce qui est banal pour un code de recherche.

Donc les noms de valeur `solidity` sont **exacts** et **ne seront pas renommés**.
C'est un geste que la version d'origine proposait, et il est **annulé**.

### 1.1 Le nombre, la commande et le périmètre — les trois étaient faux

**Le nombre.** « 117 » n'est le total de rien, sous aucune convention de comptage,
et il contredisait le tableau qui le supportait trois lignes plus bas (dont les
lignes somment à 172). Recompté le 2026-08-30 sur le périmètre **étendu**
`src include tools bench_impact/configs specs configs build_sol.cmd` :

| ce qu'on compte | commande | n |
|---|---|---|
| **lignes citant un fichier source `Y3D*.c`** | `grep -rEc "Y3D[a-z]*\.c" …` | **79** |
| … **dont** celles portant un numéro de ligne `l. NNNN` | `grep -rEc "Y3D[a-z]*\.c[^ ]* l\. *[0-9]" …` | **72** |
| lignes citant `Y3D*.c` **ou** le mot « solidity » | `grep -rEc "Y3D[a-z]*\.c\|[Ss]olidity" …` | **182** |
| occurrences de `Y3D*.c` \| de `Y3D*.c` ou « solidity » | `grep -rhoE … \| wc -l` | 79 \| 208 |

**La commande.** Celle imprimée dans la version d'origine — `grep -rc "…\|…"`,
**sans `-E`** — ne rend **rien du tout** (exit 1, aucune sortie) : en expression
régulière **de base**, le `|` est un caractère **littéral**, et le motif cherche la
chaîne « Y3D….c|Solidity », qui n'existe nulle part. Un rapporteur qui recopie la
commande conclut que le résultat principal du document n'existe pas. **Il faut
`grep -rEc`.**

**Le périmètre.** La version d'origine s'arrêtait à
`src include tools bench_impact/configs`. Elle manquait **`specs/`** (5 lignes dans
`WP7_couplage_contact.md`, 1 dans `WP6_contact_residuel.md`, 1 dans `spec.md`),
**`configs/`** (`fdem3d_bench1_insert.cfg`, `p1_banc.cfg`) et **`build_sol.cmd`** —
c'est-à-dire les spécifications de lots de travail et le script de build qui porte
le nom dans son nom.

**Ventilation à jour, sur le périmètre étendu** (lignes citant `Y3D*.c`) :

| fichier | lignes |
|---|---|
| `src/Fdem3dSolver.cpp` | 22 |
| `src/FdemSolver.cpp` | 20 |
| `include/rockim/Fdem3dSolver.hpp` | 9 |
| `include/rockim/FdemSolver.hpp` | 8 |
| `bench_impact/configs/impact_imperial_coulomb.cfg` | 7 |
| `bench_impact/configs/impact_imperial.cfg` | 6 |
| `bench_impact/configs/impact_kuru9.cfg`, `impact_kuru11.cfg` | 2 + 2 |
| `tools/verify_suite.py` | 2 |
| `specs/005-impact-insert-yang/WP7_couplage_contact.md` | 1 |
| **total** | **79** |

### 1.2 Le vrai problème n'est pas l'attribution — c'est la REPRODUCTIBILITÉ

La prémisse fausse masquait un défaut réel, et qui subsiste entier :

1. **Ces 72 citations ne sont pas vérifiables.** Elles disent `Y3Dfd.c l. 1099`
   et rien d'autre : ni dépôt, ni licence, ni **commit**, ni date de lecture. Or le
   dépôt est **activement maintenu** (dernier push relevé le 2026-03-31) : **les
   numéros de ligne bougent**. Un rapporteur qui reclone aujourd'hui ne retrouvera
   pas nécessairement ces lignes. *Une référence sans version n'est pas une
   référence.*
2. **Trois statuts sont confondus en deux.** Ce que disent les **articles publiés**,
   ce que fait le **code public**, ce que fait la **version interne** (inconnue,
   non consultable) — c'est cette distinction manquante qui a produit toute la
   confusion de la session, dans les deux sens.
3. **Quelques conventions n'ont réellement aucune source publiée**, et elles
   doivent être assumées comme des choix de rockim (table ci-dessous). C'est un
   sous-ensemble bien plus petit que ce qu'annonçait la version d'origine.

### 1.3 Table de rachat — **récrite**, sept verdicts sur huit ont changé

*Rappel de la convention : « racheté » = une source **d'article** (auteur, page,
équation) peut remplacer ou compléter la citation de code.*

| clé de rockim | attribution actuelle | verdict **révisé le 2026-08-30** |
|---|---|---|
| **`jointFailRule = majority`** | « `nfail>1`, Y3Dfd.c l. 1175 » | **RACHETÉE — et c'est la seule à DEUX sources indépendantes.** Manuscrit UCL **p. 14** : « *A joint element is labelled as failed when **at least two integration points have zero stress components*** ». Le code et le texte **concordent**. **Garder les deux citations, pas remplacer l'une par l'autre.** *(inchangé)* |
| **`jointDeltaC = guo`** | Guo 2014 éq. 2.30 | **DÉJÀ BONNE.** Confirmée par le manuscrit UCL **éq. 10 p. 12** (`Gf ≈ ⅓ f δc`). *(inchangé)* |
| **plancher `st = max(2 sp, 3 GfII/dpefs)`** | « Y3Dfd.c l. 1110-1126 » + Guo éq. 2.24/2.30 | **PARTIELLEMENT RACHETÉE.** La partie Guo tient ; le plancher `2 sp` n'a **aucune source d'article** — mais il est **relevé dans un code publié sous LGPL-3.0**, ce qui n'est pas rien. À présenter comme *convention d'implémentation relevée dans le code public d'Imperial, sans contrepartie publiée*. |
| **`jointDeltaC = solidity`** | « Y3Dfd.c l. 1099 » | ~~NON RACHETABLE, à renommer~~ → **CONVENTION DE CODE, SOURCÉE.** Le nom est **exact** : il désigne un code public identifiable. **Ne pas renommer.** À compléter d'un ancrage de version (§1.4, A3.1). |
| **taux de déformation `strainRateFilter = none`** | « ce que fait LEUR code : Y3Dfd.c l. 1448 » | ~~NON RACHETABLE, le mot « LEUR » est faux~~ → **le mot « LEUR » est JUSTE.** C'est bien leur code. Reste vrai : **aucune source d'article ne décrit la mesure du taux** — donc statut *convention de code, sans contrepartie publiée*. |
| **`gcBirth = penalty`** | « Y3Did.c l. 915-964 » (13 occurrences) | ~~NON RACHETABLE~~ → **PARTIELLEMENT RACHETÉE, et le rachat est substantiel.** Le **problème** et **une rampe** sont publiés : manuscrit UCL **p. 17** décrit exactement la difficulté (« *the overlap between tetrahedral elements due to compression will generate an initial non-zero contact force …, which can cause instability problems* ») et **l'éq. (18)** publie le remède, `f = (n_c/n_total)·f_initial`, « *n_total … (usually 10)* ». **Seul le RÉ-ÉCHELONNEMENT DE PÉNALITÉ est propre à rockim** — et il est mesuré (+936 / +179 / +27 J/m). |
| **`contactMu.<phase>`, règle de paire = MINIMUM** | « Solidity Y3Did.c l. 1292 » | **CHOIX DE ROCKIM, à assumer — mais sans reprocher un silence.** Aucune source dépouillée ne traite la combinaison pour une **paire de matériaux différents**, et pour cause : dans tous leurs essais publiés **les corps en contact sont du même matériau**, la question ne se pose pas. Le minimum se justifie physiquement (le plus faible gouverne le glissement). ⚠️ **Ne pas citer de nombre de sources** : 4, 5 et 6 circulent dans trois documents pour la même affirmation (cf. contre-audit §9, M-11). |
| **`contactDamageCoupling = solidity`** | « raideur normale ET frottement × `d_fact = min(1−D_i, 1−D_j)`, effondrement /1000 sous 0,041 (Y3Did.c l. 995, 1044, 1263-1265) » | ~~NON RACHETABLE, ET C'EST LE PLUS GRAVE~~ → **reste le point le plus délicat, mais pour une autre raison.** Le [lot 2b](2026-08-29_lot2b_couplage_endommagement_contact.md) tient **entièrement** : le mot « penalty » n'apparaît **pas une fois** dans l'article de pulvérisation, et le seul couplage (1−D) publié porte sur la **contrainte d'élément**. **Mais la forme vient bien de leur code public** — où elle est **INERTE** (`df = R0` ⇒ `d_fact = 1`). rockim a donc **activé une branche morte**, en la branchant sur son propre moteur d'endommagement (`el_[e].bdD`, exigé à `src/Fdem3dSolver.cpp:648-652`). **Ce n'est ni une réplication, ni une « transcription » : c'est un geste propre, défendable, qu'il faut décrire comme tel.** |

**Sept des huit lignes ont changé de verdict.** Le décompte des « non rachetables »
passe de **cinq à zéro** : aucune ne l'est au sens où la version d'origine
l'entendait (« sans source »). Elles se répartissent désormais en *rachetées par
un article* (2), *conventions de code publiées mais sans contrepartie d'article*
(4) et *choix propres de rockim, à assumer* (2, dont un partiel).

### 1.4 L'action — **elle change de nature, et son effort n'est pas « faible »**

Les trois gestes de la version d'origine sont remplacés par ceux de la fiche
[`A03_resourcer_attributions.md`](../chantier/A03_resourcer_attributions.md) §5 :

| # | quoi | pourquoi |
|---|---|---|
| **A3.1** | remplacer les citations **nues** `Y3Dfd.c l. 1099` par une citation **complète** : dépôt, licence, **commit** (ou date de lecture), fichier, ligne | 72 citations concernées ; le dépôt est maintenu, donc les lignes bougent, donc rien n'est vérifiable en l'état |
| **A3.2** | distinguer partout **trois** statuts : **articles** / **code public** / **version interne** | c'est la distinction manquante, et c'est elle qui a produit la confusion **dans les deux sens** |
| **A3.3** | ajouter la source **d'article** *en plus* du code là où elle existe | fait le 2026-08-30 pour `jointFailRule` (UCL p. 14) et `gcBirth` (UCL p. 17, éq. 18), dans `DOCUMENTATION_rockim.md` §5.4 quinquies |

~~**Effort FAIBLE**~~ → **effort MOYEN.** Le contre-audit l'établit : ce n'est pas
un geste mécanique mais **79 arbitrages de source**, un par site, répartis sur
10 fichiers dont deux sources de 4 900 et 7 100 lignes, plus la table de rachat à
refaire (fait ici) et les deux bilans qui la citent. **La fiche A03 le dit dans ses
propres mots** : « cela ne doit pas être fait à la va-vite un soir de session ».

**Ce qui reste vrai de la conclusion d'origine, et qu'il faut garder** :

> Tant que les citations ne sont pas ancrées sur une version, **toute phrase du
> manuscrit qui dit « comme Imperial » en s'appuyant sur ces clés est
> invérifiable** — non pas indéfendable, la nuance compte : la source existe et
> elle est publique, mais le lecteur ne peut pas la retrouver.

---

## 1bis. La version d'origine du 2026-08-29, conservée

*Ce qui suit est le §1 tel qu'il a été écrit le 2026-08-29, **inchangé**. Il est
conservé pour l'historique et parce que le raisonnement qui l'a produit est
lui-même instructif. **Il ne doit pas être cité** : sa prémisse, son nombre, sa
commande et sept de ses huit verdicts sont corrigés au §1 ci-dessus.*
*Seul changement de forme : les **niveaux de titre** de ce bloc ont été abaissés
en `#####` pour ne pas créer de doublons d'ancres avec le §1 récrit. **Aucun mot
du texte n'a été touché.***

<details>
<summary><b>Déplier le §1 d'origine (prémisse fausse — ne pas citer)</b></summary>

##### 1. LE RÉSULTAT PRINCIPAL — et ce n'était pas la question posée

> **⚠️ LA PRÉMISSE DE TOUT CE §1 EST FAUSSE. Lire d'abord
> `chantier_imperial_2026-08-29/A03_resourcer_attributions.md`.**
> `solidity-solver-open` **EST** le code d'Imperial : dépôt public
> `ImperialCollegeLondon/solidity-solver-open`, **LGPL-3.0**, C, 17 000 lignes,
> cloné et lu le 2026-08-26 — provenance documentée en **quatre endroits** du
> dépôt (`CR_solidity_2026-08-27.md:19`, `BILAN_replique_solidity_2026-08-27.md:11`,
> `DOCUMENTATION_rockim.md:409`, `tools/verify_suite.py:443`). La CORRECTION 2 du
> 29/08 avait sur-corrigé, et le brief de mission a repris sa formulation.
> **Ce qui reste vrai** : ce code n'est pas la version qui a produit l'article de
> 2026 (son facteur d'endommagement est à zéro), donc y lire une forme et en
> conclure une implémentation de l'article reste une faute — mais pour cette
> raison-là, pas parce que ce serait le code de quelqu'un d'autre.
> **Le comptage des 117 tient ; les verdicts de la table de rachat ci-dessous,
> non.** L'action A3 est suspendue et redéfinie dans la fiche du chantier.

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

##### 1.1 Table de rachat — ce qui peut désormais être re-sourcé

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
| **`contactMu.<phase>`, règle de paire = MINIMUM** | « Solidity Y3Did.c l. 1292 » | **NON RACHETABLE — et c'est un cas d'école.** Le [lot 2c](2026-08-29_lot2c_frottement_tangentiel.md) §4 a établi qu'**aucune des 8 sources dépouillées ne traite la règle de paire** — liste nominative arrêtée le 2026-08-30 (B9). ⚠️ Et ce n'est pas un silence coupable : dans tous leurs essais publiés les corps sont du **même matériau**. Le « minimum » est donc un choix de rockim. Il peut être bon ; il ne peut pas se réclamer d'eux. |
| **`contactDamageCoupling = solidity`** | « raideur normale ET frottement multipliés par `d_fact = min(1−D_i, 1−D_j)`, effondrement /1000 sous 0,041 (Y3Did.c l. 995, 1044, 1263-1265) » | **NON RACHETABLE, ET C'EST LE PLUS GRAVE.** Le [lot 2b](2026-08-29_lot2b_couplage_endommagement_contact.md) a établi que le mot « penalty » **n'apparaît pas une seule fois** dans l'article de pulvérisation, et que **le seul couplage (1−D) publié porte sur la contrainte d'élément**. rockim a donc implémenté, sous un nom qui promet une réplication, un mécanisme **qu'aucune publication d'Imperial ne décrit**. |

##### 1.2 L'action

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


</details>

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
| **règle de paire** | **DIVERGENT non sourçable** | `Fdem3dSolver.hpp:524-529` : `if (muPerPhase_) { m = min(muPhase_[A], muPhase_[B]); }` — **le plus faible gouverne**, commenté « Solidity Y3Did.c l. 1292 ». **aucune des 8 sources dépouillées ne la traite** — liste nominative au [lot 2c §4](2026-08-29_lot2c_frottement_tangentiel.md), arrêtée le 2026-08-30 (B9) | aucun | non |

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

> **⚠️ CORRIGÉ ET MESURÉ LE 2026-08-29 (chantier A11).** Le manque est réel et il
> est corrigé (`chantier_imperial_2026-08-29/A11_dt_tangentiel.md`), mais
> **l'affirmation ci-dessous est exagérée**. Effet mesuré sur le pas de temps à
> `potTangentFactor = 1,4286`, la valeur des decks d'impact : **−2,21 %** en
> insertion intrinsèque, **−4,66 %** en adaptative. Ce n'est pas un précipice :
> le budget nodal est dominé par la raideur des ressorts de joint, devant
> laquelle le terme de contact pèse peu. Le correctif se justifie par la parité
> 2D/3D et par la source, **pas par une urgence** — aucun run n'explosait à cause
> de ça. Un piège d'unités a été trouvé au passage : en 3D `potP_` est en **Pa**
> et `potKt_` en **N/m**, si bien que recopier le `max(potP_, potKt_)` du 2D
> aurait divisé le pas par ≈ 32 en silence.

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

## 3. CE QUE ROCKIM A ET QU'IMPERIAL N'A PAS — **récrit le 2026-08-30 après contre-audit**

> **Neuf des treize entrées d'origine étaient trop flatteuses.** Le contre-audit
> l'a établi entrée par entrée, et j'ai revérifié moi-même chaque preuve chiffrée
> ci-dessous. La version d'origine est conservée **intégralement** au
> [§3bis](#3bis-la-version-dorigine-du-2026-08-29-conservée).
>
> **Le motif de l'erreur, et il faut le nommer** : la consigne du brief — « ne pas
> présenter rockim comme systématiquement en retard » — a été **sur-corrigée en
> complaisance**. Une liste d'avantages qui ne survit pas à l'ouverture d'un
> fichier ne défend pas le dépôt : elle l'expose.

**Deux règles de rédaction, tirées des corrections ci-dessous, et qui valent pour
tout le manuscrit :**

1. **« Imperial ne publie pas X » est presque toujours faux, maintenant qu'on sait
   que leur solveur est public.** L'énoncé correct est **« X n'est décrit dans
   aucune PUBLICATION d'Imperial »** — ce qui reste un argument, et un argument
   vérifiable.
2. **Une capacité qui n'est pas ARMÉE sur le banc dont on parle ne compte pas
   comme un avantage sur ce banc.** Trois des treize entrées échouent sur ce seul
   test. C'est exactement la faute que le dépôt s'interdit ailleurs de sa propre
   initiative (`src/FdemSolver.cpp:829-831`).

### 3.1 La liste corrigée

| # | l'avantage, **requalifié** | statut |
|---|---|---|
| **1** | **101 contrôles de non-régression au 2026-08-30, dont 44 au tier rapide par défaut**, 91 au tier `full`, 101 au tier `all` ; **225 assertions** chiffrées. ⚠️ **Ce nombre bouge à chaque chantier** (95 avant A11, 97 après, 101 depuis B10) : le citer sans sa date, c'est refaire l'erreur du « 98 ». ~~« 98, pas 42 »~~ : 98 n'est le total de rien, la ventilation 42/45/8 était périmée **et** somme à 95, et **les tiers s'emboîtent au lieu de s'additionner** (`tools/verify_suite.py:894`). Le « 42 » du brief était **exactement le tier par défaut de l'époque**. ~~« Imperial ne publie aucune suite »~~ → **« aucune suite de non-régression n'est décrite dans leurs publications »**. | **TENU, requalifié** |
| **2** | **Un bilan d'énergie fermé et imprimé**, là où Imperial obtient amortissement et erreur numérique **par soustraction** et l'écrit (ARMA 24-0952 p. 3). ⚠️ **Deux réserves, §3.2.** | **TENU, avec réserves** |
| **3** | **Le refus plutôt que le silence.** `jointElastic = parabolic` lève si `jointSoftening` n'est ni `yan` ni `munjiza` (`src/FdemSolver.cpp:196-205`), avec en commentaire *« une capacite active et muette est indiscernable d une capacite inerte »*. ~~« Aucun code publié ne fait ça »~~ → **« cette discipline n'est décrite dans aucune publication d'Imperial »** ; et je n'ai pas dépouillé leur code sur ce point, donc je ne peux rien dire de ce qu'il fait. | **TENU, requalifié** |
| **4** | **Des avertissements croisés** quand un deck pose la clé de l'autre schéma d'insertion (`src/Fdem3dSolver.cpp:552-561`). | **TENU** |
| **5** | **Le pas de temps budgète explicitement la pénalité**, facteur 2 de la parabole compris (`src/FdemSolver.cpp:3105-3112`). ⚠️ En **3D**, la raideur **tangentielle** du contact n'y entrait pas jusqu'au chantier **A11** du 2026-08-29 ; elle y entre désormais sous `dtBudgetTangential = on` (opt-in, défaut bit-identique), effet mesuré **−2,2 %** (intrinsèque) et **−4,7 %** (adaptative). | **TENU, daté** |
| **6** | **Un pilote au point matériel** (`tools/yan_point.cpp`) qui **réutilise la fonction d'adoucissement expédiée** (`rockim/YanSoftening.hpp`) et trace σ(δ) et τ(δ). ~~« la loi expédiée est testée, pas une ré-implémentation »~~ : **à moitié seulement** — l'en-tête dit lui-même que le pilote *« reproduces »* / *« mirrors »* la mise à jour de traction de `jointForces()`, donc **cette moitié EST une ré-implémentation**. Il **ne compare rien** (il imprime quatre nombres, aucun ratio, **aucun verdict**) et **n'est appelé ni par le build ni par la suite** (`grep -c yan_point CMakeLists.txt build*.cmd` → **0** sur les 31 scripts). La valeur de l'intégrale est verrouillée dans la suite **par le solveur lui-même** (`yan_integral`, `tools/verify_suite.py:140-144`). | **AFFAIBLI** |
| **7** | ~~Un garde-fou crack-band qui lève si le plus gros élément dépasse E·G_f/f_t²~~ | **RETIRÉ.** `src/MatLaw.cpp:1304` : `if (kind == "dpr" \|\| kind == "saksala")` — **le garde-fou ne s'applique qu'à ces deux lois**. Or le deck de réplique pose `bulkModel = neohookean` (`impact_imperial.cfg:190`) : **il ne se déclenche jamais sur le banc d'impact.** Et il porte sur la loi de **volume**, alors qu'en FDEM l'énergie de rupture est portée par les **joints**. |
| **8** | **La parité 2D/3D est une RÈGLE DE CONSTITUTION** : deux lois (`YanSoftening.hpp`, `YangDif.hpp`) sont **partagées** entre les solveurs pour la garantir *par construction* — `YangDif.hpp:8-13` écrit qu'une divergence sur les bornes en dur « serait MUETTE et fausserait toute comparaison dimensionnelle ». **13 capacités** sont contrôlées des deux côtés. ~~« instrumentée »~~ : **aucun contrôle ne COMPARE une grandeur 2D à son homologue 3D** — les 13 paires portent des références **indépendantes**. Et **cinq ruptures de parité sont ouvertes** (impression de la pénalité absente en 2D ; `potXi_` 2D seulement ; `potKt_` hors budget 3D jusqu'à A11 ; le diagnostic de l'entrée 9 sans équivalent 3D ; `groupBond` **3D seulement** — `grep -c "groupBond\|gbond_" src/FdemSolver.cpp` → **0**). **C'est le REGISTRE des ruptures qui est l'avantage, pas la parité.** | **REFORMULÉ — et plus fort ainsi** |
| **9** | **Un diagnostic instrumenté du *diffuse ratcheting*** de la pénalité intrinsèque (part des joints au-dessus de D = 0,01 au pic, D moyen — `src/FdemSolver.cpp:6840-6843`). ⚠️ **Il est enfermé dans `if (scen_ == Scenario::BRAZILIAN)` et n'existe qu'en 2D** (0 occurrence en 3D). **Le périmètre est l'inverse de ce qu'on croit** : la pathologie est diagnostiquée sur un essai de **calibration 2D**, et **pas du tout sur la percussion 3D** — le cas même comparé à Imperial. | **AFFAIBLI — et il en sort une action (B8)** |
| **10** | ~~Un plafond d'impulsion dur sur le contact (20 m/s par pas et par nœud)~~ | **RETIRÉ du banc de réplique.** `capF = 20·m/dt` vit dans `generalContact()` — la branche **pénalité** (`src/FdemSolver.cpp:5019`, miroir 3D `:3685`). Le banc de réplique tourne en **`contact = potential`** (`impact_imperial.cfg:79`) : **le plafond est inerte sur la totalité du banc.** *(Il reste un vrai garde-fou pour les runs en pénalité — ne pas le supprimer, juste ne pas le revendiquer ici.)* |
| **11** | **Un levier mesuré et séparé** : pénalité contre schéma d'insertion, **à pénalité égale +1,5 point**, le reste étant de la pénalité (`BILAN_interference_2026-08-29.md:207-213` — **la mesure CORRIGÉE**, pas celle de :118-131). Aucune publication ne décrit une telle séparation. | **TENU — le plus solide de la liste** |
| **12** | **La provenance du 3000 GPa retrouvée** — règle de Turon, Dávila, Camanho & Costa (2007), K = α·E/t avec α ≈ 50 — **et une contradiction interne d'Imperial mise au jour** : Guo recommande E ≤ p₀ ≤ 10E deux phrases après avoir cité Turon, et les auteurs de l'article ont suivi Turon. ⚠️ **Coïncidence numérique plausible, non démontrée** : rien n'établit qu'ils ont appliqué Turon consciemment. **Note de bas de page, pas résultat.** | **TENU, déclassé** |
| **13** | **Une INFÉRENCE indépendante, confirmée** : l'exposant littéral 0,07 de la loi de DIF en traction publiée en 2025 ne raccorde **aucune** de ses deux bornes (**1,1245** au lieu de 1 en ε̇ = 5·10⁻⁶ /s ; **1,5160** au lieu du plateau 1,85 en ε̇ = 10² /s, soit **22 % de discontinuité**) ; l'exposant qui raccorde **simultanément** les deux vaut **0,1707** (**1,0010** en bas, **1,8500** en haut), valeur relevée indépendamment sur la courbe tracée de leur figure 2(b) ; **l'article de 2026 imprime 0,17**. ~~« une prédiction »~~ : **non soutenable** — l'article confirmateur est de la **même année** (*IJRMMS* **206**, 2026, 106660) et rien n'établit qu'il n'était pas déjà paru le 2026-08-18. ⚠️ **1,0031** figurait dans l'en-tête de `YangDif.hpp` : **c'était faux**, corrigé en **1,0010** le 2026-08-30. | **TENU, reformulé** |

### 3.2 Les deux réserves sur l'entrée 2 (bilan d'énergie)

Elles ne l'annulent pas, elles la bornent — et la deuxième est opérationnelle.

**(a) L'ARRÊT est opt-in, et le deck de réplique ne l'arme pas.**
`src/Fdem3dSolver.cpp:2316` `eAbortPct_ = cfg_.getd("budgetAbortPct", 0.0)`, puis
`:2322` `if (eAbortPct_ <= 0.0 || eAbort_) return;`. Sur les **22 decks** de
`bench_impact/configs`, **quatre** l'arment (`impact_kuru9`, `impact_kuru11`,
`impact_p2_nombres`, `impact_p2_facies`, à 2 %) — **`impact_imperial.cfg` non**.

> **La distinction juste** : la **MESURE** du résidu est **inconditionnelle**
> (`src/Fdem3dSolver.cpp:4519-4522` imprime toujours « residu : … J (… % de
> l'echelle) [OK|CHECK] ») ; c'est l'**INTERRUPTION** qui est opt-in. Formulé
> ainsi, l'argument tient devant un relecteur. Formulé comme « refuse de continuer
> s'il dérive », il tombe dès qu'on ouvre le deck.

**(b) Le poste GRAVITAIRE n'existe pas — et c'est celui qu'ARMA distingue.**
`src/Fdem3dSolver.cpp:4038` `bodyForces()` applique la pesanteur **sans aucun
compteur de travail**, et il n'existe aucun `gravWork_` dans l'en-tête. Or
`gravity = 9.81` est posé dans **20 des 22 decks**, `impact_imperial.cfg` compris, et
ARMA 24-0952 pose explicitement l'énergie potentielle gravitaire dans ses éq. 3-7.

*Magnitude honnête* : sur 600 µs d'impact les déplacements sont micrométriques,
le travail de la pesanteur est de l'ordre de **10⁻⁴ J** contre ~49 J injectés —
**invisible**. Le défaut est **structurel, pas numérique**. Il devient réel dès
qu'un run est long ou quasi statique.

*Second trou, plus gênant* : `src/Fdem3dSolver.cpp:4063`
`brushWork_ += bw; // POSTE SEPARE, jamais dans sumW`. Dès que le tri anti-gravité
des fragments est armé, **son travail tombe entièrement dans le résidu** — et un
garde-fou `budgetAbortPct` peut alors couper un run **sain** sur un artefact
purement numérique.

> ### ⚠️ MESURÉ le 2026-08-30 — ma magnitude « invisible » était la bonne réponse
> ### à la mauvaise question
>
> Le poste gravitaire **existe désormais** (`energyBodyForces`, fiche
> [`B10_bilan_energie_forces_volumiques.md`](../chantier/B10_bilan_energie_forces_volumiques.md)), et la mesure renverse la
> conclusion ci-dessus. Sur `configs/fdem3d_percussion.cfg`, `gravity = 9.81` :
>
> | T | travail de la pesanteur | résidu B4 au défaut | résidu, poste compté |
> |---|---|---|---|
> | 2·10⁻⁵ s | 1,05738·10⁻⁶ J | **1,05718·10⁻⁶ J** | −2,015·10⁻¹⁰ J |
> | 2·10⁻⁶ s | 1,71228·10⁻¹⁰ J | **1,71064·10⁻¹⁰ J** (116 %, `[CHECK]`) | −1,645·10⁻¹³ J |
>
> **Le résidu ÉTAIT le travail non compté de la pesanteur**, à 0,02 % près. La
> pesanteur est négligeable devant l'énergie **injectée** — c'est ce que je
> mesurais — et **dominante** devant l'énergie **non expliquée**, qui est
> justement ce que le verdict juge.
>
> **Et le garde-fou coupe un run sain à cause de ça, mesuré** : à `T = 9·10⁻⁶ s`
> avec `budgetAbortPct = 2`, le défaut **ABORTE** à t = 8,344·10⁻⁶ s sur un
> résidu de **175 % de l'échelle** qui vaut, à 0,3 % près, le travail de la
> pesanteur (2,37808·10⁻⁹ J) — sur un run à **zéro joint rompu**, hotspot à
> `|v| = 7,95·10⁻⁵ m/s. Clé armée, le run va au bout, résidu 0,25 %, `[OK]`.
>
> Le second trou (`brushWork_`) est traité par la même clé. **La décision d'armer
> `budgetAbortPct` sur les decks de réplique devient sûre**, ce qu'elle n'était
> pas tant que le bilan était incomplet.

> ### ⚠️ MESURÉ le 2026-08-30 — et ma magnitude « invisible » était la bonne
> ### réponse à la mauvaise question
>
> Le poste gravitaire **existe désormais** (`energyBodyForces`, fiche
> [`B10_bilan_energie_forces_volumiques.md`](../chantier/B10_bilan_energie_forces_volumiques.md)), et la mesure renverse la
> conclusion ci-dessus. Sur `configs/fdem3d_percussion.cfg` avec `gravity = 9.81` :
>
> | T | travail de la pesanteur | résidu B4 au défaut | résidu, poste compté |
> |---|---|---|---|
> | 2·10⁻⁵ s | 1,05738·10⁻⁶ J | **1,05718·10⁻⁶ J** | −2,015·10⁻¹⁰ J |
> | 2·10⁻⁶ s | 1,71228·10⁻¹⁰ J | **1,71064·10⁻¹⁰ J** (116 %, `[CHECK]`) | −1,645·10⁻¹³ J |
>
> **Le résidu ÉTAIT le travail non compté de la pesanteur**, à 0,02 % près.
> La pesanteur est négligeable devant l'énergie **injectée** — c'est ce que je
> mesurais — et **dominante** devant l'énergie **non expliquée**, qui est
> justement ce que le verdict juge.
>
> **Et le garde-fou coupe un run sain à cause de ça, mesuré** : à `T = 9·10⁻⁶ s`
> avec `budgetAbortPct = 2`, le défaut **ABORTE** à t = 8,344·10⁻⁶ s sur un
> résidu de **175 % de l'échelle** qui vaut, à 0,3 % près, le travail de la
> pesanteur (2,37808·10⁻⁹ J) — sur un run à **zéro joint rompu**. Clé armée, le
> run va au bout, résidu 0,25 %, `[OK]`.
>
> Le second trou (`brushWork_`) est traité par la même clé. **La décision
> d'armer `budgetAbortPct` sur les decks de réplique devient sûre**, ce qu'elle
> n'était pas tant que le bilan était incomplet.

### 3.3 Le compte honnête des « résultats scientifiques »

La version d'origine en annonçait **quatre** (entrées 7, 11, 12, 13). Le compte
qui résiste :

| | |
|---|---|
| **de plein droit** | **entrée 11** — la séparation **mesurée** des deux leviers, pénalité contre schéma d'insertion, à pénalité égale. Aucune publication ne fait cette étude. |
| **publiable après reformulation** | **entrée 13** — l'exposant 0,1707, présenté comme **inférence indépendante** et non comme prédiction. Et **le vrai apport n'est pas l'exposant** : c'est la conséquence physique mesurée — le saut de 22 % en ε̇ = 10² /s est un **attracteur** en insertion extrinsèque, la population insérée s'empilant juste sous le seuil (médiane **99,36 /s** avec l'exposant littéral contre **40,22 /s** avec 0,1707, mesure du 2026-08-18). Cela, aucune publication ne le décrit. |
| **note de bas de page** | **entrée 12** — la règle de Turon : coïncidence numérique plausible, non démontrée. |
| **retirée** | **entrée 7** — le garde-fou crack-band ne se déclenche pas sur le banc. |

**Un résultat solide et une reformulation publiable, ce n'est pas maigre — c'est
simplement quatre fois moins que ce que la version d'origine annonçait.** Et les
deux qui restent sont, eux, défendables ligne de code à l'appui.

---

## 3bis. La version d'origine du 2026-08-29, conservée

*Ce qui suit est le §3 tel qu'il a été écrit le 2026-08-29, **inchangé** — seuls
les niveaux de titre sont abaissés pour ne pas créer de doublons d'ancres.
**Il ne doit pas être cité** : neuf de ses treize entrées sont corrigées ci-dessus.*

<details>
<summary><b>Déplier le §3 d'origine (neuf entrées trop flatteuses — ne pas citer)</b></summary>

##### 3. CE QUE ROCKIM A ET QU'IMPERIAL N'A PAS

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


</details>

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

> # ⚠️ CE PLAN EST REMPLACÉ — NE PAS L'EXÉCUTER
> **Le contre-audit le note 2 sur 17** — la pire section du document, et la seule
> qui n'était adossée à aucune source, ni article ni ligne de code. **Le plan
> applicable est le §6 de**
> [`../chantier_imperial_2026-08-29/CONTRE_AUDIT_corrections.md`](../chantier/CONTRE_AUDIT_corrections.md)
> (actions **B1 à B9**). Ce qui suit est conservé pour l'historique.
>
> **Trois actions sont dangereuses ou impossibles telles qu'écrites** :
> * **A6** « restreindre les joints à la roche » — **désassemble l'outil** ; et son
>   effort est mal jugé dans l'autre sens : `groupBond` et son branchement existent
>   déjà (`src/Fdem3dSolver.cpp:1159-1180` et `:1436-1452`), la plomberie est une
>   clé sœur, c'est la **sémantique** (figer plutôt que supprimer) qui coûte ;
> * **A8** « balayer `gcBirthTau` » — **le solveur refuse la clé** en même temps que
>   `gcBirth = penalty`, que le deck de réplique pose déjà : `τ` y est **inerte**.
>   Et le « 518 pas contre ~10 » était faux (le `dt` venait d'un autre run) : c'est
>   **77 et 52 pas**, soit **5 à 8×** ;
> * **A11** « miroir de `FdemSolver.cpp:3134` » — appliquée littéralement, cette
>   consigne **diviserait le pas de temps par ≈ 32** (piège d'unités : en 3D `potP_`
>   est en Pa, `potKt_` en N/m). L'action a été faite correctement le 2026-08-29,
>   fiche [`A11_dt_tangentiel.md`](../chantier/A11_dt_tangentiel.md) ;
>   son gain réel est **−2,2 %** et **−4,7 %**, pas la suppression d'un précipice —
>   et le frottement, **quand il travaille**, sature au cap de Coulomb et n'apporte
>   plus aucune raideur, donc c'est précisément le régime où l'omission ne mord pas.
>
> **Et les deux « blocages » (A1, A2) n'en sont pas** : la colonne « bloque la
> réplication ? » est fausse deux fois (contre-audit §1). Les estimations d'effort
> d'A1, A2, A3 et A12 sont toutes irréalistes (§9, M-16 à M-21).

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
| **A12** | Ajouter le banc analytique du rectangle glissant | ~~faible~~ **MOYEN à ÉLEVÉ** — rockim n'a aucun scénario capable de l'exprimer, voir `chantier_imperial_2026-08-29/A12_banc_frottement.md` | premier contrôle du chemin tangentiel de CONTACT, à solution fermée | non |

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
