# A03 — Re-sourcer les attributions : **SUSPENDU**
# La prémisse sur laquelle je l'avais fondée est fausse

*Chantier du 2026-08-29. Cette fiche arrête une action et dit pourquoi.*

---

## 1. Ce que j'allais faire, et pourquoi je ne le fais pas

Le [lot 4](../biblio_insertion/2026-08-29_lot4_bilan_rockim.md) §1 annonçait, comme
résultat principal du bilan :

> « **117 attributions à un code qui n'est pas celui d'Imperial** subsistent dans
> le code, les en-têtes, la suite de vérification et les decks. »

et proposait de les « re-sourcer », c'est-à-dire de remplacer les citations
rachetables et de **renommer les clés qui promettent une réplication qu'elles ne
font pas**.

**Cette phrase est fausse. C'est bien le code d'Imperial.**

## 2. La preuve, et elle était dans le dépôt depuis trois jours

`tools/verify_suite.py:442-443` :

> « ---- les trois conventions relevees dans le CODE de Solidity ----
> (**ImperialCollegeLondon/solidity-solver-open, LGPL-3.0, lu le 2026-08-26**). »

`CR_solidity_2026-08-27.md:19-21` :

> « Une recherche a montré que ce n'était pas vrai : **leur solveur est public.**
> `github.com/ImperialCollegeLondon/solidity-solver-open` — **LGPL-3.0, C,
> 17 000 lignes, format `.Y3D`, activement maintenu (dernier push 2026-03-31).
> C'est la lignée Munjiza décrite par Guo et Yang. Le dépôt a été cloné et lu.** »

Et encore : `BILAN_replique_solidity_2026-08-27.md:11`,
`DOCUMENTATION_rockim.md:409-410` — quatre endroits, concordants, avec licence,
volume, format et date de lecture.

**[MÉTA]** Recoupement externe : le dépôt `ImperialCollegeLondon/solidity-solver-open`
existe bien sous l'organisation GitHub d'Imperial College London — C, LGPL-3.0,
activité récente. Je n'ai pas pu l'ouvrir (pas d'accès sortant dans ce
conteneur), mais quatre mentions internes concordantes plus la métadonnée
publique suffisent.

## 3. D'où vient l'erreur

Elle a une histoire, et elle est instructive.

1. **26-27/08** — le dépôt trouve le code public, le clone, le lit, et documente
   la provenance en quatre endroits. **Correct.**
2. **29/08 (soir)** — la CORRECTION 1 du `BILAN_interference` s'appuie sur cette
   lecture pour réfuter la dérivation de la pénalité à 52,6 E. **Erreur de fond**
   (elle lisait `mat.txt`, un deck d'exemple sans rapport), mais la source était
   la bonne.
3. **29/08 (nuit)** — la CORRECTION 2 constate que le facteur d'endommagement y
   est câblé à zéro, donc que ce code **ne peut pas être celui qui a produit
   l'article de 2026**. **Constat juste.** Mais elle en tire une conclusion trop
   large : « ce dossier n'est PAS le code d'Imperial ».
4. **Le brief de mission** reprend cette formulation et l'aggrave : « une
   transposition partielle que mon équipe a faite à partir de ce qui avait pu
   être lu dans les articles ». **C'est ce que j'ai cru toute la session.**

**La CORRECTION 2 a sur-corrigé.** Elle avait raison sur le fait — le facteur
d'endommagement est à zéro — et tort sur ce qu'il fallait en conclure.

## 4. Ce qui reste vrai, et il faut le dire précisément

**Le code public d'Imperial n'est pas la version qui a produit l'article de
2026.** Son facteur d'endommagement est câblé à zéro (`Y3Dfd.c` l. 749-751,
constat de la CORRECTION 2), alors que l'article publie les équations (3) et (4)
d'un modèle d'endommagement d'élément. La lecture la plus simple est celle d'une
**version ouverte en retard sur la version interne**, ce qui est banal pour un
code de recherche.

Donc l'avertissement de méthode **tient**, mais pour une raison plus étroite et
plus juste :

> Lire une **forme** dans `solidity-solver-open` et en conclure une
> **implémentation de ce que décrit l'article de 2026** reste une faute. Ce n'est
> pas parce que ce serait « le code de quelqu'un d'autre » — c'est bien le leur,
> même lignée, mêmes auteurs — mais parce que **ce n'est pas la version dont
> l'article parle**.

Et la conclusion du [lot 2b](../biblio_insertion/2026-08-29_lot2b_couplage_endommagement_contact.md)
est **inchangée** : l'article de 2026 ne mentionne pas une seule fois le mot
`penalty`, et le seul couplage (1−D) qu'il publie porte sur la contrainte
d'élément. Ce qui change, c'est le statut de la source alternative : ce n'est pas
« un code sans autorité », c'est **le code publié du même groupe, dans une
version antérieure**.

## 5. Ce que devient l'action A3

**Elle change de nature, et elle rétrécit beaucoup.**

Ce qu'il ne faut **PAS** faire — et que j'allais faire :

* ~~renommer `contactDamageCoupling = solidity` en `hypothese_rockim`~~ : le nom
  `solidity` est **exact**, il désigne un code public identifiable ;
* ~~renommer `jointDeltaC = solidity`~~ : idem ;
* ~~déclarer « non rachetable » la règle de paire du frottement~~ : elle est
  relevée dans un code publié sous LGPL-3.0. Ce n'est pas la même chose qu'une
  invention sans source.

Ce qu'il reste à faire, et qui garde du sens :

| # | quoi | pourquoi |
|---|---|---|
| **A3.1** | remplacer les citations **bare** `Y3Dfd.c l. 1099` par une citation **complète** : dépôt, licence, commit ou date de lecture, fichier, ligne | une référence sans dépôt ni version n'est pas vérifiable par un rapporteur — et le code est « activement maintenu », donc les lignes bougent |
| **A3.2** | distinguer partout **trois** statuts, et non deux : ce que disent les **articles**, ce que fait le **code public**, ce que fait la **version interne** (inconnue) | c'est la distinction qui manquait, et c'est elle qui a produit toute cette confusion |
| **A3.3** | ajouter, là où c'est vrai, la source d'article **en plus** du code — par ex. `jointFailRule = majority` a désormais **deux** sources concordantes : `nfail>1` dans le code, et « at least two integration points have zero stress components » au manuscrit UCL p. 14 | deux sources indépendantes valent mieux qu'une |

**A3.1 à A3.3 ne changent aucun comportement.** Ce sont des commentaires et des
messages. Mais ils demandent de reprendre la table de rachat du lot 4 §1.1 en
entier, sur une prémisse neuve — et cela ne doit pas être fait à la va-vite un
soir de session.

## 6. Ce que cette fiche déclenche ailleurs

* **Lot 4 §1** : le titre « attributions à un code qui n'est pas celui
  d'Imperial » et la table de rachat sont **faux dans leur prémisse**. Un
  avertissement y est porté ; la réécriture reste à faire.
* **Synthèse §5.3** : idem.
* **Le brief de mission §0** : sa prémisse est fausse. Elle est conservée telle
  quelle — on ne réécrit pas l'historique — mais elle ne doit plus être citée
  sans cette fiche.

## 7. Ce qu'il faudrait vérifier ensuite

1. **Le dépôt public a-t-il un tag ou un commit correspondant aux articles ?**
   Le CR du 27/08 note un dernier push au 2026-03-31 ; l'article de 2026 est
   postérieur. Un `git log` du clone répondrait.
2. **Le facteur d'endommagement est-il absent de TOUTE l'histoire du dépôt, ou
   seulement de la tête ?** S'il apparaît dans une branche ou un commit, la
   CORRECTION 2 tombe entièrement.
3. **Les auteurs disent-ils quelque part que le code public diffère de celui des
   articles ?** Un README ou une note de version trancherait.

Ces trois questions se règlent en dix minutes **avec le clone sous la main**.
Elles ne peuvent pas se régler depuis ce conteneur, qui n'a pas d'accès sortant.
