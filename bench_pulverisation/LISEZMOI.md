# bench_pulverisation — le banc du WP6 (μ de contact résiduel)

Créé le 2026-08-28 (spec 005, plan `WP6_contact_residuel.md` §4). Trois runs
A/B/C qui isolent l'effet de la **perte d'appui** (contactResidualMu) et le
levier du relais roche/roche (jointDeath), à physique par ailleurs complète
(coulomb + bulkDamage granite + fricScaled), sur le maillage du dépôt
`meshes/box3d_h45.msh` (80×80×60 mm, 19 326 tets, non structuré — aucune
génération requise). Sphère R 8,51 mm, 1,28 kg, 9,5 m/s, St Anne, T = 150 µs.

| run | deck | ce qui change | rôle |
|---|---|---|---|
| A | `configs/pulv_a_full.cfg` | `contactResidualMu = 0.18` | physique complète |
| B | `configs/pulv_b_ctrl.cfg` | — (clé absente) | contrôle : isole WP6 |
| C | `configs/pulv_c_death.cfg` | A + `jointDeath = damage` | relais roche/roche actif sous l'insert |

`gcRestitution = 0.2` est FIGÉ dans les trois decks (levier 2 de la revue :
c'était un défaut silencieux qui contamine toute mesure de restitution). À
balayer (0,1/0,2/0,5) seulement si e devient la métrique de calibration.

## Lancement (conteneur ou machine, ~25-40 min par run)

    build/rockim bench_pulverisation/configs/pulv_a_full.cfg out_pulv_A
    build/rockim bench_pulverisation/configs/pulv_b_ctrl.cfg out_pulv_B
    build/rockim bench_pulverisation/configs/pulv_c_death.cfg out_pulv_C

Vérifier au log de A et C : la bannière `[FDEM3D] contactResidualMu = 0.18`
ET, au résumé final, `contact residuel : N evaluations` avec **N > 0** (la
ligne « JAMAIS engagé » sur A/C = mécanisme mort = FAIL immédiat).

## Critères PASS (prédictions inscrites avant tout lancement)

1. **A vs B — perte d'appui** : enfoncement final A > B ; v finale |A| < |B|
   (plus d'énergie retenue par la roche) ; eFric(A) ≠ eFric(B) net.
   Le mécanisme doit s'engager APRÈS le début du broyage (t premier
   engagement > t premier joint cisaillé), jamais avant le contact.
2. **C vs A — relais roche/roche** : nCtcPulv(C) > nCtcPulv(A) (les paires
   roche/roche s'ajoutent aux paires outil/roche) ; morts de joints en
   compression au log de C (compteur jointDeath).
3. **Zéro artefact** : aucun engagement avant le premier contact outil ;
   budget d'énergie clos aux tolérances habituelles ; B strictement
   identique à un run pré-WP6 du même deck (bit-identique, déjà prouvé
   sur decks jumeaux le 28/08 — refaire ici au moindre doute).

Mesures : `history.csv` (toolVz, enfoncement, eFric, eGc, nPulv) + le
résumé stdout (nCtcPulv, premier engagement). Le ratio Ft/Fn à l'interface
outil s'approxime par |toolFx,toolFy|/toolFz (agrégat sphère) — une vraie
colonne dédiée reste à ajouter si le banc devient quantitatif.
