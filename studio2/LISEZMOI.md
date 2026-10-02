# Rockim, interface v2 (spec 007)

Interface web locale du code maison FDEM, lignée g1, tranche 1 : l'essai triaxial 2D.
Spécification, mesures et journal : [specs/007-studio-v2/spec.md](../specs/007-studio-v2/spec.md).

## Lancer

```
python studio2/app/serveur.py
```

puis ouvrir http://localhost:8770 dans Edge ou Chrome. Aucune dépendance à installer en dehors
de Python 3.10+ avec numpy et vtk (lecture des VTU) ; gmsh seulement pour le maillage homogène.

Options : `--espace DOSSIER` (défaut `studio2/espace`), `--port N`. `--faux-solveur` remplace g1
par `tests/faux_rockim.py` : démonstrations et tests, aucun calcul réel.

L'espace de travail contient la file (`file.json`, `decks/`, `out/`, `logs/`), les calculs
courts (`courts/`), les caches d'affichage (`caches/`) et les réglages (`reglages.json` :
exécutable, jobs, fils, dossiers de campagne à afficher).

## Écrans

| Écran | Rôle |
|---|---|
| Essai | choix physiques -> deck g1 (noyau/formulaire.py), schéma, vérifications, estimation, ajout ou variation de σ₃ |
| Maillage | aperçu du maillage RÉEL par un run g1 d'une frame (quelques secondes), estimé contre réalisé, dessin de fissures |
| Loi des joints | loi cohésive au cliquer-glisser (formules du solveur), essai éclair sur 10 grains (1 à 17 s) |
| File | calculs en cours, actions, courbe et journal en direct ; les calculs survivent à la fermeture |
| Résultats | comparaison de runs (pic filtré), détail : champs de contrainte et de déformation, fissures, temps |

Un clic sur « Générer l'aperçu », « Essai éclair » ou « Ajouter à la file » lance g1 : c'est
l'utilisateur qui décide de chaque calcul.

## Tester

```
python -m pytest studio2/tests
```

380 tests, une quarantaine de secondes. Ils ne lancent jamais g1 (faux solveur), mais lisent les
runs réels de `etude_triax_hetero/out` et `studio2/_travail` quand ils sont présents.

## Organisation

- `noyau/` : toute la logique, sans interface (essai, formulaire, validation, maillage, file,
  calculs courts, résultats, dépouillement, loi des joints) ;
- `app/serveur.py` : API JSON sur le noyau, file pilotée en tâche de fond ;
- `app/static/` : l'interface (HTML, CSS, modules JS sans compilation) ;
- `proto/` : prototypes jetables du jalon J1, gardés pour mémoire.
