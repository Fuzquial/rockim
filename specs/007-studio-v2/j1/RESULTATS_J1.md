# Jalon J1 : deux prototypes de l'écran Résultats, mesurés

**Date** : 2026-10-02
**Cas de référence** : `etude_triax_hetero/out/F7_disc_gbm_P020` (12 995 triangles, 19 339 joints, 26 frames)
**Poste** : Windows 11, GPU Intel Arc intégré, écran 60 Hz, fenêtre 1600 × 1000
**Code** : `studio2/proto/` (jetable). Lancement :

- C, web local : `python studio2/proto/web/serveur.py` puis http://localhost:8765 (`?mesure=1` pour les mesures)
- B, PySide6 + QML : `studio2\proto\qml\.venv\Scripts\python studio2\proto\qml\main.py [--mesure]`

Les deux prototypes lisent le même cache binaire (`commun/cache_run.py`). La conversion d'un run de
165 Mo de VTU ASCII prend 4,9 s et produit 15 Mo, ce qui tient N2 (< 15 s) pour les deux piles.

## Mesures

| Critère | Seuil | C, web local | B, PySide6 + QML |
|---|---|---|---|
| N1 démarrage jusqu'à l'image avec données | < 2 s | **0,45 s** | 2,7 s (7,9 s après modification du QML) |
| N3 ouverture d'un run depuis le cache | < 1 s | 0,10 s | 0,06 à 0,10 s |
| N4 changement de frame, jusqu'à l'image | < 100 ms | 16 ms moyen, 17 ms max | 13 ms moyen, 28 ms max |
| N5 zoom et déplacement | ≥ 60 i/s | 60 i/s, pire écart 17 ms | 60 i/s (178 sans vsync) |
| N6 blocage du fil d'interface | < 50 ms | 150 ms une fois, compilation des shaders au démarrage | 12 ms (démarrage non mesuré) |
| N8 installation | une commande, < 10 min | Python standard + navigateur, aucune dépendance | environnement dédié de 715 Mo ; le PySide6 global du poste est incomplet ; chemin des DLL à déclarer à la main |
| Taille du code | — | 650 lignes | 1 400 lignes |

En C, un changement de frame ne recopie rien : toutes les frames sont en mémoire GPU et la frame
courante est un entier passé au shader. Le rendu seul coûte 0,03 ms. Les 16 ms mesurées sont
l'attente de l'image suivante de l'écran.

En B, chaque changement de frame recalcule le tampon de couleurs en numpy puis le renvoie au GPU.
C'est rapide sur ce cas mais croîtra avec la taille du maillage, alors que C reste constant.

## Aspect

Captures : [capture_C_web.png](capture_C_web.png) et [capture_B_qml.png](capture_B_qml.png).

Les deux tiennent le thème « granite & teal ». B a des détails mieux aboutis, qui ne tiennent pas à
la pile et se transposent tels quels : tuiles de chiffres clés, compteurs dans les puces de
fissures, échelle graphique de 10 mm, arc-en-ciel saturé plus proche d'Abaqus. B affiche en
revanche « 8.33 ms » au lieu de « 8,33 ms » (N10).

## Fragilités observées

- **B** : `grabWindow()` fait planter le processus sans message (rendu Direct3D sur fil séparé) ;
  les erreurs QML sont muettes sous Windows sans variable d'environnement ; la vue 3D corrige les
  couleurs par défaut et délave les champs ; l'API QtGraphs bouge entre versions mineures.
- **C** : la compilation des shaders bloque 150 ms au démarrage (corrigible par compilation
  asynchrone) ; le volet du navigateur intégré à l'application Claude bride l'affichage à 5 i/s,
  les mesures ont donc été faites dans Edge.

## Recommandation

**C, application web locale.** Elle tient tous les seuils sauf un dépassement de N6 au démarrage,
qui se corrige ; B rate N1 et N8 sans réglage simple pour les rattraper. C demande deux fois moins
de code et s'installe sans dépendance. Les détails d'aspect de B sont repris dans C.

Décision à prendre par Fernando, sur ces chiffres et sur les deux captures.
