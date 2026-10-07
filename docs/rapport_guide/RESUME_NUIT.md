# Résumé de la nuit du 6 au 7 octobre 2026 (rapport-guide rockim)

Branche `claude/adoring-bardeen-qc4ovq` du dépôt rockim. Tout est poussé.

## À lire en premier

| fichier | contenu |
|---|---|
| `docs/rapport_guide/rapport_guide_rockim.pdf` | le rapport (120 pages, 7,6 Mo, compilation propre) |
| `docs/rapport_guide/rapport_guide_changements.pdf` | le même, avec les changements depuis la version d'hier soir surlignés (latexdiff ; les tableaux ne sont pas surlignés) |
| `docs/rapport_guide/JOURNAL_NUIT.md` | chaque initiative prise, avec l'heure, le motif et le statut |
| `docs/rapport_guide/simulations_a_lancer/RESULTATS_*.md` | les résultats bruts de la session de calcul |
| `docs/rapport_guide/simulations_a_lancer/NOUVEAUX_BANCS.md` | 7 nouveaux bancs trouvés en bibliographie et en ligne |
| `docs/rapport_guide/FIGURES_A_REFAIRE.md` | figures d'archive à régénérer sur le poste |

## Conclusions qui ont changé cette nuit

1. **Injection d'énergie par le contact dans le tunnel (ancienne réserve n° 1) : cause identifiée.** C'est une erreur d'intégration du potentiel de contact normal : le travail injecté suit la raideur de pénalité (−88 % à p/10, +351 % à 10 p) et le pas de temps (−41 % à dt/2). Ce n'est ni la décharge `origin` (cliquet : −16 %), ni la naissance du contact (≈ −22 %), et le frottement la masque au lieu de la créer. Le résidu du bilan qui semblait « constant » est lui aussi une erreur d'intégration. Conséquence : le faciès du tunnel dépend de la raideur de contact ; aucune comparaison de faciès entre raideurs différentes n'est permise. Diagnostic 4 (avec frottement, même maillage) : p = 1e9 N/m supprime l'injection (+2,6e6 → −8,0e5 J/m dissipés ; v max 25,9 → 5,6 m/s) ; le facteur de pointe exposait le biais sans en être la cause ; la référence adaptative n'avait déjà pas d'injection nette (−5,0e5 J/m). Mais le faciès bouge avec p : en adaptatif +43 % de joints rompus, EDZ p95 14,46 → 17,13 m, U paroi +28 % ; avec pointe, U max aux reins 0,22 → 0,54 m. Le résidu B4 ne s'améliore pas (erreur en O(dt)) et atteint 21,2 % sur la base adaptative ; l'interpénétration n'est pas mesurable. **Conclusion pratique** : p = 1e9 N/m règle l'injection mais ne laisse pas le faciès invariant ; avant d'adopter une valeur, il faut une étude de convergence en p et en dt des observables du chapitre (au moins 3e9 et 1e9 N/m). D'ici là, les chiffres de tunnel se donnent avec p et dt.
2. **Yan et al. 2023 : la campagne d'août tient avec le binaire actuel.** Les neuf calculs principaux rejoués retrouvent août (UCS +0,03 %, brésilien +0,06 %, φ = 22,92°). L'écart macOS du 4 octobre n'était pas le code ; l'hypothèse restante est le tirage du maillage de Voronoï.
3. **Calibration Red Bohus : le jeu de juillet n'est plus reproductible.** Traction 16,67 MPa (tient), compression 33 MPa contre 139 MPa. Identique avec les binaires du 19 et du 25 août : le changement est antérieur à git. Ta fiche du 5 août donne les candidats : deux correctifs de l'amortisseur de joint (141,5 → 85,4 → 32,8 MPa).
4. **Première comparaison à un essai : écaillage du granite de Bohus (Saadati 2016).** 3 critères sur 5 passent nettement (pic de vitesse −2,2 %, plan d'écaillage 55 contre 57 mm, bilan). Lecture au spot (Ø 3 mm impossible à h = 2,5 mm ; plus petit spot 6 × 6 mm) : pic 7,593 m/s (−2,0 %, passe), montée 26,2 µs (échoue, +1,2 µs au-delà de la bande), σ Novikov 22,94 MPa (échoue de 0,24 MPa). Le spot suit le point arrière (écart < 1 %), pas la moyenne de face : verdict 3 sur 5, deux critères hors bande de peu. Détachement de l'écaille avec DIF non tranché (nodif se détache entre 145 et 150 µs ; dif coupé à 130 µs). Géométrie de l'éprouvette supposée.
5. **Coupe 2D au cutter PDC (nouveau chapitre).** Le contact de Signorini ferme la pompe d'énergie de l'outil (rapport d'injection 1,09 contre 408 en août) ; la force de coupe reste 10 fois sous la cible de Heilman, cause non isolée (cinq candidats listés).
6. **Comparaison des schémas d'insertion sur Kuru (6 090 contre 27 joints) :** obtenue avec `origin` sans cliquet et avant une correction ; dans la configuration de référence l'intrinsèque ne rompt que 57 joints. La conclusion « le schéma d'insertion change la physique » est devenue conditionnelle.

## Ce qui a été rejoué cette nuit (binaire actuel)

Yan (9 calculs) ; suite de non-régression rapide (50/51, échec connu) ; calibration de juillet (4) ; diagnostic du contact du tunnel (12 calculs en 4 séries) ; coupe 2D (4) ; écaillage de Bohus (2, puis 2 pour la lecture au spot) ; brésilien dynamique de Kuru (1, arrêté : trop long) ; AbuAisha (contrôle de signe, Parker grossier, pics sur maillages de 12 et 6 mm) ; tunnel réduit (4).

## Le rapport

- Chapitres ajoutés : AbuAisha (couplage hydraulique), tunnel EDZ (Wang 2024), coupe 2D (Heilman 2024), écaillage de Bohus, annexe C (simulations à lancer).
- Titre recentré : « formulation, vérification et comparaisons » (la validation se limite à l'écaillage).
- Relecture : 4 tours de 4 à 7 relecteurs spécialisés (chiffres, fond, clarté, cohérence, français, style), éditeurs séparés, puis contrôles ciblés après chaque intégration. Delta entre versions : 30,7 % → 8,0 % → 6,2 % → 0,6 %, aucune remarque majeure restante au tour 4.
- Figures : 15 repassées en Computer Modern et virgule décimale ; PDF ramené de 20 à 7,6 Mo.

## Restent ouverts

- Choix de la raideur de contact du tunnel : le faciès dépend de p ; étude de convergence en p et dt à faire (le diagnostic 4 est intégré).
- Écaillage : détachement de l'écaille avec DIF (calcul à T = 300 µs) ; spot de Ø 3 mm seulement avec un raffinement local de la face arrière ; géométrie de l'éprouvette à confirmer ; indépendance du DIF de Yang vis-à-vis de cet essai.
- Brésilien de Yan : la valeur d'août est 5,29 (figure) ou 5,30 MPa (fiche) ; le pic absolu de Yan n'est pas chiffré dans l'article.
- Parker et AbuAisha de production, tunnel de production, non rejoués.

## À lancer sur ton poste (14 fils), par ordre d'utilité

1. Tunnel au maillage de production avec le binaire actuel : étude de convergence en p (au moins 1e10, 3e9 et 1e9 N/m) et en dt des observables du chapitre (R_EDZ, U paroi, joints rompus) (≈ 45 min par calcul).
2. Écaillage de Bohus, variante dif, à T = 300 µs (≈ 30 min à 4 fils), pour trancher le détachement de l'écaille.
3. Bissection de la régression de la calibration avec `rockim_f1` / `rockim_f2` (correctifs du 5 août).
4. Jumeau grossier de Kuru en insertion intrinsèque et adaptative dans la configuration de référence.
5. Les 8 brésiliens dynamiques de Kuru à durée complète (45 à 60 min chacun).
6. St Anne prolongé à 500 µs, puis sur le maillage fin `rock073` ; Kuru fin jusqu'au retournement.
7. AbuAisha et Parker de production ; Yan à la maille de 0,75 mm ; bancs V&V longs.

Liste complète et commandes : annexe C du rapport et `simulations_a_lancer/STRATEGIE.md`.

## Questions pour toi

- Quel nom et quelle géométrie exacte pour l'éprouvette d'écaillage de Saadati 2016 ?
- Faut-il garder la borne ℓ_cz < a dans le guide (violée par St Anne et Kuru) ?
- Valeur du brésilien de Yan en août : 5,29 ou 5,30 MPa ?

## Ajout de 6 h 30 : enquête dans le code du contact (lecture seule, rien modifié)

Rapport complet : `docs/rapport_guide/ENQUETE_CONTACT.md` (preuves par une copie instrumentée du code compilée hors dépôt, programmes de test dans `enquete_contact/`). Le rapport-guide n'a PAS encore été mis à jour avec ces conclusions : à valider d'abord.

1. **Cause principale de l'injection** : la recherche de contact ne regarde que les triangles ayant une arête active (`FdemSolver.cpp:7134-7139`). Deux triangles qui ne partagent qu'un sommet ne sont reliés par aucun joint (`:7237`) et peuvent se recouvrir sans être vus. Quand un joint voisin meurt, la paire apparaît d'un coup avec ce recouvrement accumulé, et l'énergie p∫(φA+φB) est créée sans travail. La rampe `gcBirth` ne l'amortit pas : `gcBirthTau` = 1 µs < dt = 5,8 µs (`:1428`, `:7281-7283`). Mesuré sur `tip16_mu0` : 9,09e6 J/m créés à la naissance pour 8,20e6 J/m de travail net du contact ; 100 % viennent de paires à sommet commun. Explique la dépendance en p (énergie ∝ p à recouvrement fixé).
2. **Résidu B4 : sans lien avec le contact.** La correction d'intégration est calculée avant que les rouleaux n'annulent la vitesse (`:9148` puis `:9163`) : l'énergie retirée par les rouleaux égale exactement le résidu, ∝ dt.
3. **Défaut secondaire** : les forces nodales de Munjiza ne dérivent pas d'une énergie quand les triangles se déforment (terme en trop p·∫φ·∇λ), 0,3 % de l'injection ici.
4. **Compteur du travail du contact** : biais O(dt) négatif ; un travail positif est donc une vraie injection, plutôt sous-estimée. Deux commentaires du code ont le signe faux (`:9914-9916`, `:10620-10630`).
5. **Écarté** : le pas de temps inclut bien la raideur du contact (`:5222`).

Correctifs proposés (non appliqués, avec leur test dans le rapport) : inclure les voisins par sommet des éléments actifs dans la recherche (ou neutraliser l'énergie d'une naissance à recouvrement profond et l'imprimer au bilan) ; avertir si `gcBirthTau` < dt ; appliquer la contrainte des rouleaux avant la correction d'intégration ; corriger les forces de Munjiza ; corriger les deux commentaires.

Conséquence pour le rapport : la conclusion « erreur d'intégration de la raideur » (diagnostics 3 et 4) est à remplacer par « énergie de recouvrement créée à la naissance des paires à sommet commun, proportionnelle à p » ; le résidu B4 des tunnels vient des rouleaux. À faire après ta validation.
