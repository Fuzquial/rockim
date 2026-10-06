# Journal des initiatives de la nuit du 6 au 7 octobre 2026

Autonomie accordée par le doctorant jusqu'à 7-8 h, initiatives non validées autorisées.
Chaque initiative est inscrite ici AVANT d'être lancée, puis mise à jour avec son résultat.

## Règles que je m'applique

| règle | pourquoi |
|---|---|
| Une initiative = une ligne ici (motif, périmètre, statut, commit) | traçabilité, relecture rapide au réveil |
| Tout se fait sur la branche `claude/adoring-bardeen-qc4ovq`, en commits séparés | chaque initiative se retire par un `git revert` |
| Aucune suppression ; un fichier remplacé garde son original (git ou `_orig/`) | règle 6 du CLAUDE.md |
| Code du solveur : aucune modification poussée. Au plus un correctif PROPOSÉ (fichier .patch + note), jamais intégré | le code fait foi pour toute la thèse ; décision du doctorant |
| Calculs : courts seulement (< 30 min, cloud) ; les longs sont listés pour le poste | consigne du 6/10, 22 h 20 |
| Tout chiffre nouveau dans le rapport : source `fichier:ligne`, vérifié par un second agent | même règle que la relecture |
| Une initiative qui change une conclusion du rapport est signalée en tête du résumé du matin | le doctorant doit la voir avant ses encadrants |
| Pas de PR, pas de message externe, rien hors des deux dépôts | consignes de session |

## Initiatives

| heure (Paris) | initiative | motif | statut |
|---|---|---|---|
| 00 h 20 | Chapitre « coupe 2D au cutter PDC » + 3-4 rejeux courts | demande du doctorant ; même famille de défaut que l'injection du tunnel | en cours |
| 00 h 25 | Recherche et préparation de nouveaux bancs (biblio + web), decks prêts sans calcul complet | demande du doctorant | en cours |
| 00 h 30 | Diagnostic de l'injection d'énergie par le contact (3 calculs, une clé chacun) confié à la session de calcul | réserve n° 1 du rapport | envoyé |
| 01 h 30 | Intégration du rejeu de Yan (session de calcul, RESULTATS_yan.md) : nouvelle section sec-yan-rejeu, résumé, introduction, discussion, annexe C | **change une conclusion** : la campagne d'août tient avec le binaire courant (UCS +0,09 %) ; l'écart macOS du 4/10 vient du tirage de Voronoï | fait |
| 01 h 30 | Relecture : tours 1 et 2 appliqués (delta 30,7 % puis 8,0 %) | méthode validée par le doctorant | fait |
| 01 h 45 | Nouveaux bancs : 7 retenus (écaillage Bohus, brésilien dynamique Kuru, Lamb, anneau, Petersson, taillant Saksala 2014, Bobet-Einstein), decks prêts pour 3 | demande du doctorant | fait (NOUVEAUX_BANCS.md) |
| 01 h 50 | File de la session de calcul complétée : écaillage (2) et SHPB Kuru (8) après le diagnostic du contact, avant l'impact | bancs de VALIDATION contre l'essai, courts ; plus utiles que l'impact prolongé | envoyé |
| 02 h 10 | Chapitre « coupe 2D au cutter PDC » (r05d) intégré, 4 rejeux courts | demande du doctorant | fait. **Constat nouveau** : Signorini ferme la pompe d'outil (R_inj 0,9 contre 408 en août), mais la décharge `origin` sans cliquet crée de l'énergie dans les joints (115 J/m d'un coup) ; avec `jointSecantRatchet = on` le calcul va au bout. Même décharge que tout le tunnel : piste n° 1 pour l'injection du tunnel |
