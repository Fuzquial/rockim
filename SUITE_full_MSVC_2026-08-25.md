# Suite `--tier full` sous MSVC : 42/50, et les 8 échecs sont de plateforme

*Constat du 2026-08-25, établi pendant le chantier « DIF intrinsèque ».
À lire avant de s'alarmer d'un tier full rouge sous Windows.*

## Le fait

`python tools/verify_suite.py --exe rockim_d1.exe --tier full` donne **42/50**
sur cette machine (MSVC 2022, Windows). Les huit échecs sont :

| test | mesuré | attendu (baseline Linux) |
|---|---|---|
| `ucs_yan_adaptive` | broken 313, inserted 1308 | 327, 1288 |
| `ucs_yan_origin` | inserted 1087 | 1131 |
| `percussion_2d` | broken 50 | 174 |
| `percussion_2d_gcadaptive` | broken 50, gcact 19 | 174, 21 |
| `percussion_2d_potential` | broken 6 | 5 |
| `shpb_mini` | broken 902 | 812 |
| `shpb_mini_gcadaptive` | broken 903 | 812 |
| `shpb_mini_potential` | broken 572 | 578 |

## Ce ne sont PAS des régressions — c'est mesuré

Les huit tests ont été relancés sur `rockim_d0.exe`, **compilé depuis `main`
sans aucune modification**. Ils échouent avec des valeurs **identiques au
chiffre près** : 313/1308, 1087, 50, 50/19, 6, 902, 903/122, 572.

Autrement dit, la cause est antérieure à tout chantier en cours, et l'identité
des valeurs entre les deux binaires constitue par ailleurs une **preuve de
bit-neutralité sur le tier full** plus forte que celle du tier fast.

## L'explication était déjà écrite

`DOCUMENTATION_rockim.md` §3.3 : « Références en dur (**baseline Linux
2026-08-11 — re-baseliner une fois sous MSVC**) ». Et son §8, règle 9 :
« Reproductibilité : garantie par `seed` PAR binaire ; **MSVC et libstdc++
tirent des nombres différents à graine égale** (Voronoï, phases) — re-baseliner
par plateforme. » C'est aussi le chantier **F1** de la ROADMAP, toujours ouvert.

La signature du symptôme le confirme : les grandeurs **physiques** passent
(`ucs_mpa` 51,04 contre 51,07 ; `dampWork` ≤ 0 partout ; `gcwork`, `gcact`
conformes) et seuls les **comptages entiers** sur des runs chaotiques
(percussion en phase débris, SHPB, UCS post-pic) divergent — exactement ce
qu'un tirage aléatoire différent et une associativité flottante différente
produisent.

## Ce qu'il faut en faire

Re-baseliner le tier full sous MSVC (chantier F1), en gardant les références
Linux à côté plutôt qu'en les écrasant : les deux plateformes sont légitimes et
la suite gagnerait un jeu de références par plateforme. **Tant que ce n'est pas
fait, le critère de non-régression sous Windows est la COMPARAISON avant/après
sur la même machine, pas le compteur PASS/FAIL.** C'est le protocole qui a été
suivi ici.

Les tiers `fast` (21/21) et les contrôles à charge nulle passent, eux, sans
réserve sous MSVC.
