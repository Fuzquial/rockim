# CHANTIER « réplique Imperial » — ouvert le 2026-08-29

Ce dossier porte **toute la documentation** des changements apportés à rockim à
la suite de la mission d'état de l'art
([`MISSION_etat_de_l_art_2026-08-29.md`](../MISSION_etat_de_l_art_2026-08-29.md)
et sa [synthèse](../SYNTHESE_etat_de_l_art_2026-08-29.md)).

**Règle du dépôt appliquée ici sans exception** : le code source est modifié
**en place**, mais toute capacité nouvelle est **opt-in**, le **défaut reste
bit-identique**, elle **s'annonce** au journal, et elle est **verrouillée par un
contrôle de non-régression**. Une capacité active et muette est indiscernable
d'une capacité inerte — et un **risque** silencieux est pire encore.

## Ce que chaque fiche doit contenir

1. **La source**, avec auteur, année, page et numéro d'équation. Pas de source =
   pas de changement.
2. **Ce que le code faisait avant**, avec ses numéros de ligne.
3. **Ce qu'il fait maintenant**, avec la clé qui l'arme et son défaut.
4. **La preuve de bit-identité au défaut**, mesurée.
5. **La mesure de l'effet** quand la clé est armée — un chiffre, pas une
   intention.
6. **Les contrôles de non-régression** ajoutés, par leur nom.
7. **Ce que le changement ne fait pas**, et les réserves honnêtes.

## Index

**[`CONTRE_AUDIT_corrections.md`](CONTRE_AUDIT_corrections.md) — à lire en premier.**
Dix vérificateurs indépendants ont noté le bilan du lot 4 à **52/129**. Ce document
porte les corrections et **le plan applicable**. Il prime sur le lot 4 et sur le
§6 de la synthèse.

| fiche | action | état |
|---|---|---|
| [`A11_dt_tangentiel.md`](A11_dt_tangentiel.md) | porter la raideur tangentielle du contact dans le budget de pas de temps du 3D | **FAIT** |
| [`../SOURCES_SOLIDITY.md`](../SOURCES_SOLIDITY.md) | **la provenance des citations `Y3D*.c`, une fois pour toutes** — dépôt, licence, date de lecture, trois statuts, et la table des **13** références auxquelles se réduisent les 72 citations | **FAIT (B4a)** — reste **B4b** : vérifier les 13 contre le clone et relever le commit, impossible depuis ce conteneur |
| [`B08_ratcheting_percussion_3d.md`](B08_ratcheting_percussion_3d.md) | **porter le diagnostic de *diffuse ratcheting* sur la percussion 3D** — il n'existait qu'en 2D et sur l'essai brésilien | **FAIT** — × 3,24 de joints ratchetant en insertion intrinsèque, pic d'effort inférieur de 8,1 %, **0 joint rompu** de part et d'autre. Et un avertissement mesuré : la jauge ne se compare pas d'une pénalité à l'autre |
| [`B10_bilan_energie_forces_volumiques.md`](B10_bilan_energie_forces_volumiques.md) | **fermer le bilan d'énergie** : le septième poste d'ARMA (l'énergie gravitaire) n'existait pas, et le travail du tri était hors bilan — clé `energyBodyForces`, **2D et 3D** | **FAIT** — le résidu B4 **était** le travail non compté de la pesanteur, à 0,02 % près, et le garde-fou `budgetAbortPct` **coupait un run sain** à cause de lui |
| [`A12_banc_frottement.md`](A12_banc_frottement.md) | banc analytique du rectangle glissant | **SPÉCIFIÉ, non implémenté** — estimation d'effort corrigée, trois voies chiffrées, décision au commanditaire |
| [`A03_resourcer_attributions.md`](A03_resourcer_attributions.md) | re-sourcer les attributions | **SUSPENDU** — la prémisse était fausse : `solidity-solver-open` EST le code public d'Imperial (LGPL-3.0). Action redéfinie en A3.1-A3.3 |
| `A01_gcbirth_penalty.md` | essai `gcBirth = penalty` sur l'impact 3D | à faire |
| `A02_longueur_h.md` | aligner la longueur de référence h sur celle d'Imperial | **suspendu** — attend la confirmation du contre-audit sur le facteur 2,4495 |

## ⚠️ Incident du 2026-08-30 — une branche miroir laissée non compilable

*Consigné ici parce que le mode de panne est générique et se reproduira.*

**Ce qui s'est passé.** Le dépôt porte trois branches qui doivent rester
synchronisées. Deux avancent ensemble (`claude/fdem-imperial-research-9xfhly`
et `joint-handoff`) ; la troisième, `claude/rockim-recent-mpka9c`, est un
**miroir plus étroit** et reçoit les mêmes changements par `git cherry-pick`.
Trois conflits sont survenus le 2026-08-30 et ont été résolus par
`git checkout --theirs`, après avoir **vérifié à chaque fois que la version
retenue était un sur-ensemble du FICHIER en conflit**.

**La vérification était juste pour le fichier, et sans valeur pour l'arbre.**
Elle a posé un `src/Fdem3dSolver.cpp` de 4 995 lignes à côté d'un
`include/rockim/Fdem3dSolver.hpp` resté 160 lignes plus court, qui ne déclare
pas les membres que ce `.cpp` utilise. **La branche n'a plus compilé pendant
trois commits** — neuf erreurs `was not declared in this scope` — sans que rien
ne le signale, parce que **aucun `cherry-pick` n'avait été suivi d'une
compilation**.

**La règle qui en découle, et elle est courte :**

> **Un `cherry-pick` entre branches divergentes se vérifie en COMPILANT, pas en
> comparant des lignes triées.** Comparer les fichiers dit ce qu'un fichier
> contient ; seule la compilation dit si l'arbre tient debout. Et quand la
> résolution consiste à imposer une version, la question n'est pas « ce fichier
> est-il complet ? » mais « **quels AUTRES fichiers ce fichier suppose-t-il ?** »

**La réparation, et comment elle a été contrôlée** (commit `ed31047`) :
`src/`, `include/` et `tools/verify_suite.py` alignés sur la branche de travail,
après avoir établi qu'aucune capacité du miroir n'était perdue — **0 clé de deck
et 0 fonction membre** présentes avant et absentes après, sur **les deux**
solveurs (175 → 190 clés en 2D, 126 → 141 en 3D ; 57 → 61 et 37 → 41 fonctions).
Puis, depuis un **checkout propre** de la branche réparée : `cmake` + `ninja`
sans erreur, les 13 decks référencés par la suite tous présents et **identiques**
à ceux de la branche de travail, maillages identiques, et la suite au tier par
défaut rejouée sur ce binaire-là.

## Ce que ce chantier ne touche pas

Aucun défaut existant n'est modifié. Aucun deck existant n'est modifié sans que
la fiche correspondante le dise. La suite des 104 contrôles (au 2026-08-30 — ce nombre bouge à chaque chantier) reste valable telle
quelle : les contrôles ajoutés ici s'y ajoutent, ils n'en remplacent aucun.
