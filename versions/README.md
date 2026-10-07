# versions/ — historique des binaires du chantier coupe PDC

Le dossier ne contient qu'un fichier, [`HISTORIQUE.md`](HISTORIQUE.md), écrit en août 2026 pendant
le chantier de coupe au cutter PDC 2D. Il trace trois binaires Windows successifs et les trois
correctifs de contact qu'ils ont introduits, chacun derrière une clé optionnelle dont le défaut
garde le comportement historique.

| Binaire | Correctif ajouté | Clé | Verdict mesuré |
|---|---|---|---|
| `rockim_tun.exe` | état de départ (étude tunnel, `cutterThick`) | — | référence `out_cut_v3` |
| `rockim_a1.exe` | rafraîchissement anticipé des faces libérées, puis pénalité de contact adaptative (EPFL, Ghesquière-Dierickx et al. 2025) | `gcSurfaceRefresh`, `jointContactPenalty` | le premier aggrave, le second divise l'énergie injectée par l'outil par 14 |
| `rockim_a2.exe` | plafond d'impulsion outil | `toolImpulseCap` | en cours d'évaluation à la date du document |

Le document mentionne des instantanés de source `v0_tunnel/`, `v2_epfl/` et `v3_impulse/`. Ils ne
sont pas dans ce dépôt, pas plus que les trois binaires.

## Où est la suite de l'historique

| Période | Source |
|---|---|
| août 2026, coupe PDC | `versions/HISTORIQUE.md` |
| depuis le 5 septembre 2026 (tag `g0-0.1.0`) | [`CHANGELOG.md`](../CHANGELOG.md), une section par série de binaires `rockim_g1*` |
| empreintes de bit-identité de chaque binaire | `results/bitid_*.json`, voir [`results/README.md`](../results/README.md) |
| binaire Windows versionné | `rockim_g1y19.exe` à la racine, celui du run St Anne |
