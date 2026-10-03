# Patchs d'optimisation du 2026-10-03 — NON VERIFIES, NON APPLIQUES

Trois commits produits par l'agent d'optimisation de la session cloud (base `245e7cb`), sauvegardes
ici pour ne pas les perdre. Annonces bit-identiques par l'agent, mais **pas encore controles** :
ni fusionnes, ni compares bit a bit par une autre execution.

- `0001` fdem3d : contact par potentiel et joints plus rapides
- `0002` fdem 2D et fem3d : fusions et grilles de contact
- `0003` VtkWriter : formatage ASCII parallele

Appliquer sur la branche (depuis la racine du depot) :

    git am patches/optimisation_2026-10-03/000*.patch

Puis, AVANT de les garder : compiler, rejouer `tools/verify_suite.py --tier fast`, et comparer les
sorties octet pour octet avec le binaire d'avant (`git stash` / binaire de la branche sans les patchs)
sur les decks de `tests_f2/bitid/` a nombre de fils egal (`cmp` sur history.csv, frames.csv et les
fichiers finaux). Les ancres de `tools/bitid_refs.json` sont des empreintes Windows/MSVC : elles ne se
comparent pas a un binaire Linux ou macOS, il faut comparer avant/apres sur la meme machine.
Si un deck differe : ne pas garder le patch concerne.
