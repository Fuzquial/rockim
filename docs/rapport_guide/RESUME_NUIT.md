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

1. **Injection d'énergie par le contact dans le tunnel (ancienne réserve n° 1) : cause identifiée.** C'est une erreur d'intégration du potentiel de contact normal : le travail injecté suit la raideur de pénalité (−88 % à p/10, +351 % à 10 p) et le pas de temps (−41 % à dt/2). Ce n'est ni la décharge `origin` (cliquet : −16 %), ni la naissance du contact (≈ −22 %), et le frottement la masque au lieu de la créer. Le résidu du bilan qui semblait « constant » est lui aussi une erreur d'intégration. Conséquence : le faciès du tunnel dépend de la raideur de contact ; aucune comparaison de faciès entre raideurs différentes n'est permise. La valeur p = 1e9 N/m n'est **pas encore** validée (diagnostic 4, voir plus bas).
2. **Yan et al. 2023 : la campagne d'août tient avec le binaire actuel.** Les neuf calculs principaux rejoués retrouvent août (UCS +0,03 %, brésilien +0,06 %, φ = 22,92°). L'écart macOS du 4 octobre n'était pas le code ; l'hypothèse restante est le tirage du maillage de Voronoï.
3. **Calibration Red Bohus : le jeu de juillet n'est plus reproductible.** Traction 16,67 MPa (tient), compression 33 MPa contre 139 MPa. Identique avec les binaires du 19 et du 25 août : le changement est antérieur à git. Ta fiche du 5 août donne les candidats : deux correctifs de l'amortisseur de joint (141,5 → 85,4 → 32,8 MPa).
4. **Première comparaison à un essai : écaillage du granite de Bohus (Saadati 2016).** 3 critères sur 5 passent nettement (pic de vitesse −2,2 %, plan d'écaillage 55 contre 57 mm, bilan), 2 sont à la limite selon le point de lecture. Verdict partagé ; géométrie de l'éprouvette supposée.
5. **Coupe 2D au cutter PDC (nouveau chapitre).** Le contact de Signorini ferme la pompe d'énergie de l'outil (rapport d'injection 1,09 contre 408 en août) ; la force de coupe reste 10 fois sous la cible de Heilman, cause non isolée (cinq candidats listés).
6. **Comparaison des schémas d'insertion sur Kuru (6 090 contre 27 joints) :** obtenue avec `origin` sans cliquet et avant une correction ; dans la configuration de référence l'intrinsèque ne rompt que 57 joints. La conclusion « le schéma d'insertion change la physique » est devenue conditionnelle.

## Ce qui a été rejoué cette nuit (binaire actuel)

Yan (9 calculs) ; suite de non-régression rapide (50/51, échec connu) ; calibration de juillet (4) ; diagnostic du contact du tunnel (11 calculs en 4 séries) ; coupe 2D (4) ; écaillage de Bohus (2) ; brésilien dynamique de Kuru (1, arrêté : trop long) ; AbuAisha (contrôle de signe, Parker grossier, pics sur maillages de 12 et 6 mm) ; tunnel réduit (4).

## Le rapport

- Chapitres ajoutés : AbuAisha (couplage hydraulique), tunnel EDZ (Wang 2024), coupe 2D (Heilman 2024), écaillage de Bohus, annexe C (simulations à lancer).
- Titre recentré : « formulation, vérification et comparaisons » (la validation se limite à l'écaillage).
- Relecture : 4 tours de 4 à 7 relecteurs spécialisés (chiffres, fond, clarté, cohérence, français, style), éditeurs séparés, puis contrôles ciblés après chaque intégration. Delta entre versions : 30,7 % → 8,0 % → 6,2 % → 0,6 %, aucune remarque majeure restante au tour 4.
- Figures : 15 repassées en Computer Modern et virgule décimale ; PDF ramené de 20 à 7,6 Mo.

## Restent ouverts

- Diagnostic 4 (p/10 avec frottement, sur la référence et la variante pointe) : lancé, résultat à intégrer (voir l'état ci-dessous).
- Écaillage : lecture moyennée sur le spot laser (Ø 3 mm) ; géométrie de l'éprouvette à confirmer ; indépendance du DIF de Yang vis-à-vis de cet essai.
- Brésilien de Yan : la valeur d'août est 5,29 (figure) ou 5,30 MPa (fiche) ; le pic absolu de Yan n'est pas chiffré dans l'article.
- Parker et AbuAisha de production, tunnel de production, non rejoués.

## À lancer sur ton poste (14 fils), par ordre d'utilité

1. Tunnel de production avec le binaire actuel, à la raideur de référence et à p/10 (≈ 45 min chacun).
2. Bissection de la régression de la calibration avec `rockim_f1` / `rockim_f2` (correctifs du 5 août).
3. Jumeau grossier de Kuru en insertion intrinsèque et adaptative dans la configuration de référence.
4. Les 8 brésiliens dynamiques de Kuru à durée complète (45 à 60 min chacun).
5. St Anne prolongé à 500 µs, puis sur le maillage fin `rock073` ; Kuru fin jusqu'au retournement.
6. AbuAisha et Parker de production ; Yan à la maille de 0,75 mm ; bancs V&V longs.

Liste complète et commandes : annexe C du rapport et `simulations_a_lancer/STRATEGIE.md`.

## Questions pour toi

- Quel nom et quelle géométrie exacte pour l'éprouvette d'écaillage de Saadati 2016 ?
- Faut-il garder la borne ℓ_cz < a dans le guide (violée par St Anne et Kuru) ?
- Valeur du brésilien de Yan en août : 5,29 ou 5,30 MPa ?
