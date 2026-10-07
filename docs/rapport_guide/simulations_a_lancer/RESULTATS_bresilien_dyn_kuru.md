# Résultats — banc B11, brésilien dynamique de Kuru (SHPB) : NON CONCLUANT, à relancer sur le poste

Calcul du 2026-10-06/07 (UTC) : `shpb_v05_nodif` de 23 h 57 à 00 h 27, **arrêté par la règle des
30 min** (`timeout`, code 124) à t = 502,5 µs sur T = 747 µs (67 %). Session cloud Linux, binaire
`f638b0f`, 4 fils OpenMP. Sortie partielle : `sim_out/nouveaux_bancs/bresilien_dyn_kuru/shpb_v05_nodif/`
(`history.csv` jusqu'à 502,5 µs ; pas de VTU de joints, le résumé de fin n'est pas écrit).

**Les sept autres decks SHPB n'ont pas été lancés**, et `bd_qs` non plus, comme demandé.

## Pourquoi arrêter la file

Le pic de σ_t se lit sur l'onde transmise (jauge M2, `shpbMonitor2` = 2,4408 m, 0,8 m après le
disque). Dans ce calcul :

| t (µs) | ε_M1 | ε_M2 | sxxC (MPa) | joints rompus |
|---|---|---|---|---|
| 140 | −1,4e-8 (arrivée de l'onde incidente) | 0 | 0 | 0 |
| 301 | −5,50e-4 (palier incident) | 0 | 0 | 0 |
| 322 | −5,18e-4 | 0 | −3,0 (l'onde atteint le disque) | 0 |
| 382 | −2,33e-4 | 0 | −42,8 | 33 (premiers joints rompus) |
| 442 | +3,7e-7 | 0 | −76,8 | 2 347 |
| 462 | +2,5e-5 | −6,0e-9 (arrivée de l'onde transmise) | −268,0 | 8 163 |
| 502,5 | +1,40e-4 | −7,35e-5, **encore en montée** | −379,2 | 17 738 |

L'onde transmise n'arrive à M2 qu'à 462 µs et monte encore à 502,5 µs : **le pic de σ_t n'est pas
dans la fenêtre calculée** (σ_t = 2 E_b ε_M2 × 0,022/(π D) vaut 5,0 MPa au dernier pas, contre
24,5 MPa attendus par Saksala). Le calcul ralentit en outre fortement après la rupture : 60 % des
pas en 15,5 min, puis 7 % seulement dans les 14,5 min suivantes, à mesure que le nombre de joints
insérés et rompus explose (17 738 joints rompus à 502 µs : le disque est broyé aux appuis, pas
seulement fendu).

Chaque deck SHPB demande donc environ 45 min à 1 h à 4 fils sur ce conteneur, au lieu des 20 à
25 min estimées dans NOUVEAUX_BANCS.md §5.3 (l'estimation ne comptait pas le surcoût des joints).
Avec la règle des 30 min, les huit calculs s'arrêteraient tous avant le pic transmis : 4 h de
calcul sans résultat exploitable. La file a été arrêtée à 00 h 27 et le temps reporté sur `impact`.

## Ce qu'on peut déjà lire (informatif, `nodif`, v = 5 m/s)

- Arrivée de l'onde incidente au disque vers 320 µs (1,6 m à 5 048 m/s = 317 µs) : cohérent.
- Premiers joints rompus à 380 µs, environ 60 µs après l'arrivée de l'onde.
- `sxxC` est NÉGATIF pendant toute la montée (−95 MPa à 402 µs) : le centre du disque est en
  compression selon x, alors que le brésilien doit y mettre de la traction. À vérifier avant toute
  exploitation : convention de signe de `sxxC`, orientation de l'axe de chargement (x = axe des
  barres ?), ou écrasement local aux appuis qui domine la mesure au centre.

## À faire sur le poste local (14 fils)

Les huit `shpb_v*_{nodif,dif}` à T complet, puis `bd_qs`. Ou en cloud, avec une limite d'au moins
1 h par calcul. Avant cela, vérifier le signe de `sxxC` sur un deck court.
