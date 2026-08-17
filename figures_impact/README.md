# Dépouillement des impacts 3D rockim — banc P1 (2026-08-14)

Figures, GIF et scripts des simulations d'impact du 14 août 2026, **rapatriés le
2026-08-17** depuis le dossier temporaire de session où ils avaient été produits
(`AppData\Local\Temp\claude\...\scratchpad`) — un dossier que Windows peut vider
sans avertissement. Les données sources, elles, sont dans `..\out_*`.

## Les trois calculs

| Run | Config | Ce qu'il est |
|---|---|---|
| `..\out_smoke` | `configs\smoke_impact.cfg` | dégrossissage multithread, 6 frames, 54 Mo |
| `..\out_banc_mid` | `configs\p1_banc_mid.cfg` | **le banc de référence** : bloc 120³ mm, insert R = 11 mm à 8 m/s, ~82k tets, T = 1,2e-4 s (contact + charge-décharge Hertz + rebond). ~1 h 30 multithread. 243 Mo |
| `..\out_p1_1t` | `configs\p1_banc.cfg` | mesure de performance P1, 842k tets, 1 thread, 1 frame seulement, 345 Mo |

## Les figures

| Fichier | Contenu |
|---|---|
| `banc_mid_fdelta.png` | **force-pénétration** du banc moyen — la courbe propre |
| `compare_fdelta_abaqus_rockim.png` | la même, superposée aux essais FDEM Abaqus (`simulations\FDEM\percussion_2d_rod`) |
| `banc_mid_progress.png` | avancement du run (pénétration, force, joints rompus) |
| `banc_mid_impact.gif` / `impact_f5.png` | l'impact en coupe, 6 frames |
| `banc_mid_topview.gif` / `topview_f3.png` | **vue de dessus avec le maillage** |
| `banc_mid_impacted.gif`, `..._env.gif` (+ `_f5.png`) | éléments touchés, variante enveloppe |
| `smoke_impact.png` | courbes du smoke test |
| `p1_progress.png` | avancement partiel du run P1 mono-thread |

## Les scripts

Tous se lancent sans argument depuis **ce dossier** (`python plot_fdelta.py`) et
trouvent seuls leurs entrées. Vérifiés le 2026-08-17 : les 8 tournent.

| Script | Lit | Écrit |
|---|---|---|
| `plot_fdelta.py`, `plot_mid.py` | `history_mid.csv` (local) | `banc_mid_fdelta.png`, `banc_mid_progress.png` |
| `plot_p1_progress.py` | `history_snapshot.csv` (local) | `p1_progress.png` |
| `plot_smoke.py` | `..\out_smoke\history.csv` | `smoke_impact.png` |
| `gif_banc_mid.py`, `gif_topview.py`, `gif_impacted.py` | les VTU de `..\out_banc_mid` | les GIF |
| `compare_abaqus_rockim.py` | le npz Abaqus + `history_mid.csv` | la comparaison |

`gif_impacted.py` accepte `--envelope`.

## Deux corrections apportées au rapatriement

1. **Chemins morts.** `plot_smoke.py` pointait sur `Downloads\rockim_p1\...`
   (disparu depuis le déplacement du clone) et les trois `gif_*.py` cherchaient un
   dossier `out_banc_mid_snap` qui n'existait que dans le temp. Repointés en
   relatif sur `..\out_smoke` et `..\out_banc_mid`.

2. **`out_banc_mid\history.csv` a sa dernière ligne tronquée** — 26 colonnes au
   lieu de 28, terminée par `,-` : le run a été interrompu sans que le tampon
   d'écriture soit vidé. Une seule ligne sur 2016, la dernière. Les `gif_*.py`
   plantaient dessus ; ils lisent maintenant avec `invalid_raise=False` et
   l'ignorent (l'avertissement numpy est normal). La copie locale
   `history_mid.csv`, prise pendant le run, est intacte — c'est pourquoi les
   `plot_*.py` n'ont jamais eu le problème.

   À noter comme comportement de rockim : **pas de garantie de vidage de
   `history.csv` à l'arrêt**. Rejoint le chantier « précision d'écriture
   history.csv » déjà identifié.
