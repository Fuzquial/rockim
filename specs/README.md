# specs/ — spécifications des développements

Spécifications écrites selon la méthode Spec Kit : un dossier par développement, avec `spec.md`
(le besoin et les critères d'acceptation), parfois `plan.md`, `tasks.md` et une liste de contrôle
`checklists/requirements.md`. Le statut déclaré est celui de l'en-tête de `spec.md` ; l'état
constaté est celui du code au 7 octobre 2026.

| Spécification | Créée le | Statut déclaré | État constaté |
|---|---|---|---|
| [001 — pilier performance, architecture du banc percussif](001-pilier-performance-architecture/spec.md) | 14/08/2026 | clarifiée, prête pour le plan ; `plan.md` et `tasks.md` écrits | les 21 tâches de `tasks.md` restent non cochées ; les performances ont ensuite été traitées hors de ce cadre (section performances du CHANGELOG, 03/10) |
| [002 — registre de lois de joint cohésif](002-registre-lois-joint/spec.md) | 17/08/2026 | brouillon | non réalisé sous cette forme : pas de clé `jointLaw` dans le code ; la loi de joint évolue par clés séparées |
| [003 — cutter PDC en 3D](003-cutter-pdc-3d/spec.md) | 18/08/2026 | brouillon | réalisé en partie : `src/ToolPdc3d.cpp`, `rockim selftest-pdc3d`, deck `configs/cut3d_heilman.cfg` |
| [004 — couplage hydro-mécanique 2D](004-couplage-hydro-mecanique/spec.md) | 19/08/2026 | brouillon, rien n'est implémenté | réalisé en partie : clé `hydro = on` du mode `fdem`, validée par le banc [`bench_abuaisha/`](../bench_abuaisha/README.md) (`VALIDATION_hydro.md`) |
| [005 — impact à insert unique (Yang et al. 2025-2026)](005-impact-insert-yang/spec.md) | 21/08/2026 | brouillon, rien n'est implémenté | réalisé hors de ce cadre : decks `configs/yang2026_*` et `configs/stanne2025_*`, run St Anne du 14/09 |

Les spécifications de lois actuellement suivies sont à la racine de `docs/` :
[`SPEC_loi_note_2026.md`](../docs/SPEC_loi_note_2026.md).
