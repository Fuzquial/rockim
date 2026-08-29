# MISSION — état de l'art Imperial / Solidity, et bilan des manques de rockim
*Brief rédigé le 2026-08-29 pour une session dédiée. Commanditaire :
F. Uzquiano. À exécuter dans l'ordre des lots ; le livrable final est le lot 5.*

---

## 0. AVERTISSEMENT — à lire avant toute chose

**`/home/user/solidity` N'EST PAS le code d'Imperial.** C'est une transposition
partielle, faite par l'équipe, à partir de ce qui avait pu être lu dans les
articles. Ne l'utilise **pas** comme source primaire, et ne cite jamais une de
ses lignes comme « ce que fait Imperial ».

La preuve, si tu as besoin de t'en convaincre — `src/Y3Dfd.c` l. 749-751,
dans `CauchyTet4` :

    /*calculate damage factor */
    df=R0;
    *deldam=df;

Le facteur d'endommagement est **câblé à zéro** : les paramètres `d1pem0`,
`d1pemf`, `d1pedm` sont lus au deck et transportés jusqu'ici, mais le calcul
de D à partir de δ_m — l'équation (4) de Yang et al. 2026 — est absent. Donc
`d1df[] = 0`, donc `d_fact = 1`, donc le couplage `penalty *= d_fact` et
`mu = mud*d_fact` de `Y3Did.c` est **inerte**. La tuyauterie existe, la
physique n'y est pas.

La session précédente s'est fait piéger : elle a lu une **forme** et en a
conclu une **implémentation**, puis a écrit dans un bilan « je l'ai lu dans
leur source ». C'est corrigé (CORRECTION 2 de `BILAN_interference_2026-08-29.md`),
mais l'erreur est instructive : **ce dossier peut servir de piste
d'architecture, jamais de preuve.**

---

## 1. Le mandat, en une phrase

Établir, **à partir de sources publiées et vérifiables**, comment fonctionne
réellement le FDEM d'Imperial (Solidity / Y3D / VGeST) sur les problèmes
d'impact et de forage percussif — puis dresser le bilan honnête de ce qui
manque à rockim pour faire la même chose.

---

## 2. LOT 1 — cartographier les sources sur Solidity

**But** : savoir ce qui est publiquement connu du solveur, et où sont les trous.

À chercher : articles de revue, actes de conférence, thèses de doctorat,
rapports techniques, dépôts de code, documentation, pages de projet (VGeST,
Virtual Geoscience Workbench), manuels de Munjiza.

Noms à couvrir : **A. Munjiza** (l'auteur du FDEM et de ses livres),
**J.-P. Latham**, **J. Xiang**, **X. Yang**, **N. Guo** (ou l'auteur de la
thèse 2014 en notre possession), et les co-auteurs récurrents du groupe
*Applied Modelling and Computation Group* d'Imperial.

**Livrable 1** : une bibliographie annotée. Pour chaque source — référence
complète, ce qu'elle apporte (formulation ? algorithme ? validation ?
paramètres ?), et son **statut d'accès** : libre, dépôt institutionnel
(Spiral d'Imperial héberge souvent les post-prints), ou **paywall**
(Elsevier, Springer, Taylor & Francis…).

> **Les paywalls ne sont pas un obstacle final** : Fernando a un accès
> institutionnel Mines Paris-PSL. Dresse la liste précise de ce qu'il faut
> qu'il télécharge, par ordre de valeur, avec la raison en une ligne. Ne
> tente pas de contourner un paywall — demande la pièce.

---

## 3. LOT 2 — les articles d'impact, équations et algorithmes

**But** : reconstituer la formulation exacte, équation par équation.

Deux articles sont déjà en notre possession et **partiellement dépouillés** :
* Yang et al., *IJRMMS* **191** (2025) 106125 — calcaire St Anne et grès
  Rhune, 7 critères de validation ;
* Yang et al., *IJRMMS* **206** (2026) 106660 — granite Kuru Grey,
  **pulvérisation**. Fiche déjà écrite :
  [`biblio_insertion/yang2026_pulverisation.md`](biblio_insertion/yang2026_pulverisation.md)
  — la relire d'abord, ne pas refaire ce travail.

À reconstituer et à consigner, avec numéro d'équation et de page :
* la **loi de joint cohésif** (modes I, II, mixte), y compris la forme de la
  courbe d'adoucissement et la définition exacte des plages ;
* l'**algorithme de pulvérisation** : définition de δ_m, seuils, ce que
  devient la raideur, et surtout **comment l'endommagement se couple au
  contact** (raideur normale ? frottement ? les deux ?) — c'est le point
  qui a coûté le plus cher à la session précédente ;
* le **DIF** (facteur d'accroissement dynamique) et son armement ;
* le **contact** : potentiel de Munjiza, détection, pénalité, frottement, et
  la règle de choix des propriétés pour une **paire de matériaux différents** ;
* le **retrait des fragments** et la mesure de la masse détachée ;
* les **calibrations** publiées, avec leurs unités, et ce que les auteurs
  disent de leur transférabilité d'une roche à l'autre.

**Livrable 2** : un mémo « formulation Imperial », équations numérotées,
chaque ligne portant sa source (référence, page, numéro d'équation).

---

## 4. LOT 3 — les schémas d'insertion et le maillage adaptatif

**But** : trancher une question qui a fait perdre une journée.

À établir :
* ce que fait exactement l'insertion **intrinsèque** (tous les joints présents
  dès l'origine) : quelle pénalité élastique, comment la souplesse artificielle
  qu'elle introduit est traitée, et si les auteurs corrigent le module pour la
  compenser ;
* ce que fait l'insertion **adaptative** / extrinsèque (joints créés à la
  demande) : critère, seuil, et qui l'utilise ;
* le **maillage adaptatif** de la thèse de Guo (raffinement en cours de calcul,
  critères, coût) — Fernando peut fournir la thèse **par extraits** (elle est
  lourde) : demande-lui les sections précises dont tu as besoin, ne demande pas
  le fichier entier ;
* le **périmètre des joints** : sont-ils posés dans tous les corps, ou dans la
  roche seule ? Avec quelle conséquence sur le pas de temps et sur la raideur
  de l'outil ?

**Livrable 3** : un mémo comparatif des schémas, avec leurs coûts et leurs
effets connus, et une recommandation motivée pour rockim.

---

## 5. LOT 4 — le bilan rockim contre l'état de l'art

**But** : la question de Fernando, « ce qu'il manque à rockim ».

Méthode : pour chaque élément de formulation établi aux lots 2 et 3, statuer
sur rockim par lecture du code (`src/Fdem3dSolver.cpp`, `src/FdemSolver.cpp`,
`include/rockim/`) :

| statut | signification |
|---|---|
| **PRÉSENT** | implémenté et conforme — donner la ligne |
| **PARTIEL** | implémenté autrement, ou incomplet — dire en quoi |
| **ABSENT** | à écrire — estimer l'effort |
| **DIVERGENT** | rockim fait un autre choix, délibérément — donner la raison et le bilan qui la porte |

Trois exigences :
1. **Ne pas présenter rockim comme systématiquement en retard.** Il a des
   choses qu'Imperial n'a pas — une suite de non-régression à 42 contrôles,
   des capacités opt-in à défaut bit-identique, un bilan d'énergie fermé,
   des bannières obligatoires. Le dire.
2. **Chiffrer l'effort** de chaque manque, et le **gain attendu**.
3. **Distinguer** ce qui bloque une réplication de ce qui est un raffinement.

**Livrable 4** : le tableau de bilan, une ligne par élément, trié par valeur.

---

## 6. LOT 5 — le plan

**Livrable final** : un document unique qui présente les quatre lots, se
termine par un **plan de travail ordonné** — quoi implémenter, dans quel
ordre, à quel coût, avec quel critère de réussite mesurable pour chaque étape
— et distingue clairement ce qui est établi par source, ce qui est inféré, et
ce qui reste ouvert.

---

## 7. Ce qui est DÉJÀ acquis — ne pas le refaire

* `biblio_insertion/yang2026_pulverisation.md` — l'article 2026 dépouillé
  (δ_m est une **longueur** ; 0,18 est le frottement de la **roche** dans
  leur Table 1 ; le mécanisme **dépend de la roche** : rupture au rebond pour
  le calcaire St Anne, à la pénétration pour le granite).
* `biblio_insertion/guo2014_*.md` — quatre fiches sur la thèse Guo (maillage
  §2.4, contact et intégration §2.3.4-5, couplage fluide ch. 5, Dolosse ch. 6).
* `BILAN_interference_2026-08-29.md` — l'expérience adaptatif/intrinsèque et
  **ses deux corrections** ; lire les corrections, elles annulent une partie
  du corps du texte.
* `specs/005-impact-insert-yang/` — les lots WP1 à WP7 déjà implémentés.
* `HANDOFF_2026-08-29.md` — l'état des runs et les questions ouvertes.

---

## 8. Exigences de méthode (non négociables)

* **Traçabilité** : toute affirmation sur Imperial porte sa source — auteur,
  année, page, équation. Sans source, c'est une hypothèse, et elle est
  étiquetée comme telle.
* **Distinguer lu / inféré / supposé.** C'est l'erreur qui a coûté cher cette
  fois-ci.
* **Ne rien affirmer sur le code d'Imperial** à partir de `/home/user/solidity`.
* Les fiches vont dans `biblio_insertion/`, datées, en français, au format des
  fiches existantes.
* Commit et push sur `joint-handoff` **et** `claude/rockim-recent-mpka9c`.
