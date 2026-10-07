# biblio_insertion : revue bibliographique de l'insertion cohésive extrinsèque et adaptative

## 1. Objet

Six volets de revue bibliographique produits les 24 et 25 août 2026 pour instruire une question du doctorant : le critère d'insertion des joints et la manière de les insérer sont-ils incorrects dans rockim ? Le déclencheur est le constat du tunnel de Hutou Beishan et du banc d'impact : l'insertion adaptative rompt beaucoup de joints courts en nuage là où le schéma intrinsèque découpe des blocs.

## 2. Statut

Archive bibliographique, terminée le 25 août 2026. Chaque source est marquée lue ou de mémoire (`[VERIFIE]`, `[V]`, `[C]` contre `[MEMOIRE]`, `[M]`, `[NV]` selon le volet) ; les sources de mémoire sont à revérifier avant toute citation dans le manuscrit. La première exécution de la revue a été interrompue par une limite de session ; seul le volet 2 en a réchappé, les cinq autres ont été relancés avec une lecture des PDF par tranches bornées, d'où le nombre de sources marquées de mémoire.

Ce qui a changé depuis : la correction du contact du 7 octobre 2026 (`docs/rapport_guide/ENQUETE_CONTACT.md`) montre qu'une partie de la fragmentation observée venait d'une injection d'énergie par le contact, et non du seul schéma d'insertion ; le rejeu du tunnel (`docs/rapport_guide/rejeu_tunnel_fix/RESULTATS.md`) retire environ 24 % des joints rompus en adaptatif mais l'anneau reste granulé. Les conclusions de cette revue sur le schéma d'insertion restent à lire avec ce correctif en tête. La note de fond la plus récente est `docs/rapport_guide/NOTE_INSERTION.md`.

## 3. Résultat principal

- La littérature distingue trois régimes de déclenchement, pas deux : intrinsèque, extrinsèque à seuil, et extrinsèque appliqué à un champ lisse dans une zone à la limite d'écoulement, qui dégénère en insertion « en tapis » (`biblio_fdem.md` §1, `biblio_roches.md` §1).
- Trois pathologies documentées s'additionnent dans la configuration de rockim : discontinuité temporelle à l'activation (Papoulia et al.), absence de mécanisme qui limite l'insertion dans un champ lisse à seuil unique (Zhou et Molinari, Molinari et al.), dépendance au maillage (`biblio_patho.md` §1).
- Un remède publié récurrent est l'hétérogénéité des résistances par facette (Weibull, Zhou et Molinari 2004 ; Tang, RFPA). Chez Paulino et Celes, le critère n'est qu'un des trois piliers, avec la topologie (quels nœuds dupliquer) et le maillage (quelles directions sont offertes) ; corriger le seul critère ne suffit pas (`biblio_paulino.md` §1).
- Conclusions opérationnelles : [docs/notes/BILAN_insertion_adaptative.md](../docs/notes/BILAN_insertion_adaptative.md), §2 (ce que dit la littérature) et §7 (idées écartées avec la raison).

## 4. Figures

Sans objet : dossier de texte.

## 5. Contenu du dossier

| fichier | contenu |
|---|---|
| `biblio_volet2_fondations.md` | Camacho et Ortiz 1996, Ortiz et Pandolfi 1999, Pandolfi et Ortiz 2002, Ruiz et al. : généalogie du schéma extrinsèque, critère effectif, évaluation nodale, absence d'ordonnancement et de plafond, résistances par facette dès 2002 |
| `biblio_impl.md` | Yan, Zheng et Wang, IJRMMS 169 (2023) 105439, confronté ligne à ligne à l'implémentation de rockim (`insertionSweep`, `activateJoint`, gardes de continuité, pénalité) |
| `biblio_patho.md` | pathologies reconnues : discontinuité temporelle à l'activation (Papoulia et al.), dépendance au maillage du nombre de fragments (Molinari), remèdes publiés (Weibull par facette, Zhou et Molinari 2004) |
| `biblio_paulino.md` | école Paulino et Celes : TopS, PPR, isotropie des chemins (pinwheel, Rimoli et Rojas), biais de maillage |
| `biblio_fdem.md` | codes FDEM : Munjiza 2004 (intrinsèque), MultiFracS, Y-Geo et Irazu, Fukuda, HOSS ; comparaisons publiées intrinsèque contre extrinsèque |
| `biblio_roches.md` | applications aux roches : Tang et RFPA (hétérogénéité de Weibull nécessaire à la localisation), Lisjak et l'EDZ, BPM, GBM, tunnels profonds |

Une copie identique des six volets se trouve dans `tunnel_edz/biblio/`.

## 6. Rejouer

Sans objet. La revue alimente le rapport-guide : chapitre tunnel (`docs/rapport_guide/sections/r05c_tunnel.tex`, sous-section sur le schéma d'insertion) et discussion (`r08_discussion.tex`), dans [rapport_guide_rockim.pdf](../docs/rapport_guide/rapport_guide_rockim.pdf).
