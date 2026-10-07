# Lire en premier (pour un encadrant)

État au 7 octobre 2026.

## En quinze minutes

1. Le résumé et la discussion du rapport-guide :
   [docs/rapport_guide/rapport_guide_rockim.pdf](docs/rapport_guide/rapport_guide_rockim.pdf),
   pages 1 à 3 et dernier chapitre. Le rapport décrit ce que fait rockim, comment il est vérifié et
   ce qui reste ouvert.
2. Le diagnostic et la correction du contact (section ci-dessous), qui changent la lecture des
   calculs de tunnel et d'impact : [docs/rapport_guide/ENQUETE_CONTACT.md](docs/rapport_guide/ENQUETE_CONTACT.md).
3. La note sur l'insertion des joints, point de méthode encore ouvert :
   [docs/rapport_guide/NOTE_INSERTION.md](docs/rapport_guide/NOTE_INSERTION.md).

## Résultats établis

- Vérification : loi de joint (énergie de rupture restituée à mieux que 0,3 % en traction et en
  cisaillement), contact d'outil de Signorini, bilan d'énergie, onde dans une barre, solutions de
  Kirsch et de Lamé. Suite de non-régression : 65 contrôles sur 66 (un échec connu).
- Reproduction de Yan et al. 2023 (insertion adaptative) : les neuf calculs principaux, rejoués le
  6 octobre avec le binaire actuel, retrouvent la campagne d'août (résistance uniaxiale à +0,03 %,
  angle de frottement 22,92° pour 22,87° chez Yan).
- Couplage hydraulique (AbuAisha et al. 2017) : contrôles analytiques tenus (Lamé à −2,0 %,
  croisement de volume à moins de 0,5 %) ; la pression de rupture est surestimée de 7 à 25 %.
- Première comparaison à un essai : écaillage du granite de Bohus (Saadati et al. 2016), trois
  critères sur cinq tenus.

## Limites à connaître

- Correction du contact (7 octobre 2026). Le contact par potentiel créait de l'énergie : deux
  triangles ne partageant qu'un sommet n'étaient pas surveillés, se recouvraient, puis la paire
  « naissait » avec ce recouvrement. Dans le tunnel, l'énergie ainsi créée dépassait l'énergie de
  fissuration. Le correctif est disponible par options (`contactCandidates = vertex`,
  `gcBirth = offset`, `potForceExact = true`), vérifié par un test indépendant
  ([docs/rapport_guide/correction_contact/](docs/rapport_guide/correction_contact/)) : le contact ne
  crée plus d'énergie. Il n'est pas encore le comportement par défaut. Les résultats de tunnel et
  d'impact antérieurs sont à relire avec cette réserve.
- Insertion adaptative : même avec le contact corrigé, elle produit une fissuration diffuse
  (deux tiers des blocs sont un seul triangle) là où l'article de référence découpe des blocs. Le
  critère d'insertion, fondé sur la contrainte moyenne de deux éléments, est le premier suspect
  ([docs/rapport_guide/NOTE_INSERTION.md](docs/rapport_guide/NOTE_INSERTION.md)).
- Le jeu de calibration Red Bohus de juillet n'est plus reproductible (compression simple 33 MPa
  contre 139 MPa), à la suite de deux correctifs de l'amortisseur de joint du 5 août.
- Les répliques d'impact (St Anne, Kuru) sont des comparaisons de code à code, sur un seul maillage,
  et la réplique de Kuru échoue.

## Branches

La branche de référence est `main`. Les branches `claude/*`, `correctif/*`, `g0`, `g1`, `f2-*` et
`article-exact` sont des branches de travail ou d'archive.
