# GUIDE rockim-studio — l'interface FDEM (spec 006)

État : M0 livré (modèle/validation/lancement/suivi), M1 en cours (scène 3D
opérationnelle). Ce guide couvre ce qui MARCHE aujourd'hui, validé en
filière complète (gabarit → validation → run → résultats) sur tunnel EDZ,
impact 3D et hydro-frac Abu-Aisha.

## 1. Installation (une fois)

```bat
pip install -r tools\studio\requirements.txt
```

(PySide6, pyvista, pyvistaqt, matplotlib, numpy ; gmsh viendra avec M2.)

## 2. Lancer

```bat
cd tools\studio
python -m rockim_studio                    :: vide (nouveau cas fdem)
python -m rockim_studio ..\..\configs\tunnel_bore_fast.cfg
```

Premier démarrage : menu **Calcul → Exécutable rockim…** pour pointer ton
`rockim.exe` (mémorisé ensuite). La spinbox de la barre = OMP_NUM_THREADS
(0 = laisser l'environnement).

## 3. L'écran

* **Gauche — arbre du modèle** : les groupes (Général, Maillage, Matériau,
  Joints, Contact, Outil, CL, Hydro, Sorties…) filtrés par le mode courant ;
  les clés explicitement posées apparaissent en enfants avec leur valeur.
  Cliquer un groupe ouvre son formulaire à droite.
* **Droite — propriétés** : un widget par clé (énumérations en menu, bornes
  contrôlées, notation `50e9` acceptée). Le **défaut du solveur s'affiche en
  placeholder grisé** ; une valeur posée passe en gras avec un bouton ↺
  pour revenir au défaut. Infobulle = rôle + bornes (DOCUMENTATION §5).
* **Centre — Résultats 3D / Courbes** : scène VTK (rotation/zoom souris,
  slider temporel, champ au choix, case « fissures » = joints rompus en
  rouge sur bulk translucide) ; courbes de history.csv (live pendant le
  run, complet après).
* **Bas — console** : journal du solveur + onglet Validation (erreurs
  bloquantes en rouge : virgule décimale, bornes, énumérations, meshFile
  introuvable, clés inconnues).

Ctrl+Z / Ctrl+Shift+Z : annuler/rétablir. La disposition des panneaux est
mémorisée.

## 4. Les trois filières de production

Menu **Fichier → Nouveau depuis un modèle FDEM** — chaque entrée ouvre une
COPIE du deck de référence (jamais d'écrasement : l'enregistrement demande
un nouveau nom).

### Tunnel EDZ
Gabarits : rapide / production / Weibull (`configs/tunnel_bore*.cfg`).
Maillage requis : `python tools/make_unstructured_mesh.py tunnel 0.1 0.1
0.01 0.002 meshes/tunnel_2d.msh 1`. Si le .msh manque, la validation le dit
AVANT le lancement (avec la commande à copier depuis l'en-tête du deck).

### Impact 3D
Gabarits : smoke (bench1 réduit, quelques minutes — à passer AVANT de payer
le banc), banc St Anne s1,5 et sa variante FIDÈLE (spec 005). Maillages :
`make_unstructured_mesh.py bench1 …` / `tools/make_impact_mesh.py`.

### Hydro-frac Abu-Aisha (spec 004)
Gabarits : ISO grossier (`hf_iso_hydro_c.cfg`, hydro=on, pompe à débit),
ISO/ANISO production. Le groupe **Hydro** de l'arbre expose les clés du
couplage (hydroRate — 20 l/s = 0.02 —, fluidBulk avec son avertissement
d'hypothèse, boreCX/CY/SelectR…). Maillage :
`python tunnel_edz/tools/make_circle_mesh.py 8.0 8.0 0.05 0.012 0.4 0.3
meshes/hf_bore_c.msh 1` (grossier ; production : hFine 0.003).

## 5. Lancer et dépouiller

**F5** : validation (bloquante si erreurs) → choix du `out_dir` → la config
est COPIÉE dans le dossier de sortie (`studio.cfg`, meshFile absolutisé —
traçabilité : chaque out_* garde sa config exacte) → run avec journal et
courbe live (`historyFlush`). **Shift+F5** : arrêt propre (terminate).
En fin de run OK, la scène 3D charge automatiquement les frames ; sinon
**Ctrl+R** ouvre n'importe quel dossier `out_*` existant.

## 5 bis. La sonde nodale (« XY data from ODB »)

Sur n'importe quel run terminé (Ctrl+R pour l'ouvrir) : onglet **Courbes**,
panneau **Sonde nodale**. Saisir x/y (ou cocher « sonde au clic » dans
Résultats 3D et cliquer un nœud), cocher les variables, **Tracer au nœud** —
un sous-graphe par variable, axe temps partagé. Variables servies :
déplacements u_x/u_y/(u_z)/u_mag reconstruits de x(t) − x(0), les champs
aux points (velocity_*) et les champs de la maille adjacente (damage,
sigmaXX, vonMises…). ⚠️ Résolution temporelle = les FRAMES du run (clé
`frames`) — pour un historique dense type capteur, il faudra la sonde côté
solveur (spec 006 §7 S4).

## 6. Tests

```bash
python3 tools/studio/tests/test_registry_roundtrip.py   # registre + 104 cfg
python3 tools/studio/tests/test_model.py                # modèle + fumée GUI
python3 tools/studio/tests/test_results.py              # série VTU
python3 tools/studio/tests/test_end_to_end.py           # run réel complet
```

Après toute modification des clés dans le C++ :
`python3 tools/studio/dev/extract_keys.py` puis committer le JSON (le mode
`--check` est le garde anti-dérive à brancher en suite).

## 7. Limites connues (état M1 partiel)

* Pas encore de géométrie/maillage INTÉGRÉS (M2) : les .msh se génèrent par
  les scripts documentés dans l'en-tête de chaque deck.
* Pas de picking de faces pour poser les CL (M3) : elles s'éditent par les
  clés, groupe « Conditions aux limites ».
* Un seul run affiché à la fois dans les courbes (superposition : fin M1).
* Le rendu DEM (glyphes) n'existe pas et n'existera pas — studio = FDEM.
