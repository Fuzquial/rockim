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
| [`A12_banc_frottement.md`](A12_banc_frottement.md) | banc analytique du rectangle glissant | **SPÉCIFIÉ, non implémenté** — estimation d'effort corrigée, trois voies chiffrées, décision au commanditaire |
| [`A03_resourcer_attributions.md`](A03_resourcer_attributions.md) | re-sourcer les attributions | **SUSPENDU** — la prémisse était fausse : `solidity-solver-open` EST le code public d'Imperial (LGPL-3.0). Action redéfinie en A3.1-A3.3 |
| `A01_gcbirth_penalty.md` | essai `gcBirth = penalty` sur l'impact 3D | à faire |
| `A02_longueur_h.md` | aligner la longueur de référence h sur celle d'Imperial | **suspendu** — attend la confirmation du contre-audit sur le facteur 2,4495 |

## Ce que ce chantier ne touche pas

Aucun défaut existant n'est modifié. Aucun deck existant n'est modifié sans que
la fiche correspondante le dise. La suite des 98 contrôles reste valable telle
quelle : les contrôles ajoutés ici s'y ajoutent, ils n'en remplacent aucun.
