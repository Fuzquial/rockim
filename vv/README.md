# Campagne V&V de rockim : bancs, decks et entrées

Ce dossier contient les bancs de la campagne décrite dans `docs/VV_campagne.md`. Chaque banc est
un dossier autonome qui génère ses maillages, écrit ses decks, lance rockim et compare les sorties
à sa référence avec des critères écrits avant le calcul.

## Organisation

| Dossier | Banc | Référence | Solveurs |
|---|---|---|---|
| `V1_onde/` | P1.1 onde de compression dans une barre ; balayage de la pénalité | d'Alembert | fem3d, fdem3d |
| `P1_elastique/` | P1.2 Kirsch en 3D ; P1.3 cube biaxial (patch test) | élasticité linéaire | fem3d, fdem3d |
| `P2_bloc/` | P2.1 bloc glissant jusqu'à l'arrêt | L = v²/(2µg) (Xiang 2009, Fukuda 2020) | fdem3d |
| `P2_sphere/` | P2.2 sphère sur un plan | Hertz (Johnson 1985) | fdem3d |
| `P3_fissure_guo/` | P3.1 fissure pressurisée, sensibilité au maillage ; B1 | Guo 2014 §2.4 ; LEFM | fdem3d |
| `P3_aile_lisjak/` | P3.3 fissure en aile ; B2 | Lisjak 2013 §3.4.3 ; quatre modèles de LEFM | fdem 2D |
| `P4_essais_guo/` | P4.1 flexion trois points ; P4.2 brésilien ; B1 | Guo 2014 ch. 3 ; théorie des poutres, Hondros | fdem3d |
| `B_articles/` | B3 à B7 : AbuAisha 2017, tunnel de Wang 2024, UCS de Yan 2023, impact de Yang 2025 et 2026 ; faisabilité de Fukuda 2020 | valeurs publiées et résultats rockim antérieurs | fdem, fdem3d |
| `B_impact_coupe/` | benchmarks d'impact et de coupe tirés de la bibliographie (Aising, Saksala, études PDC, Labra et Rojek, LeBaron) | valeurs publiées | fdem3d, fem3d, fdem |

## Conventions

Chaque banc contient :

- un script `<id>.py` avec les actions `prepare`, `run`, `analyse` et `all`, et une fonction
  `JOBS(args)` qui renvoie la liste des runs (objets `Job` de `vvcommon.py`) ;
- les critères d'acceptation, en constantes du script, écrits avant le premier calcul ;
- un `README.md` : objectif, problème, référence et formules, mise en œuvre, métriques, critères,
  commande pour reproduire ;
- `resultats.json` et les figures `fig_*.pdf` (versionnées) et `fig_*.png` (aperçus, ignorés) ;
- `meshes/` (maillages générés, ignorés par git) et `out/` (sorties brutes, ignorées).

Les bancs n'utilisent que des clés existantes de rockim et ne modifient pas `src/`. Un run est
considéré comme fini quand son journal `rockim.log` contient la ligne `wall time` ; la file saute
les runs finis, ce qui permet de reprendre une campagne interrompue.

## Lancer

Un banc seul :

```bash
python3 vv/P2_bloc/p2_bloc.py all
```

Plusieurs bancs en parallèle, en respectant slots × fils ≤ nombre de cœurs :

```bash
python3 vv/run_queue.py P1_elastique P2_sphere P3_fissure_guo --slots 4 --threads 2
```

puis l'analyse de chaque banc (`python3 vv/<banc>/<id>.py analyse`). Le binaire par défaut est
`build_nofma/rockim` (Release, OpenMP de Homebrew, `-ffp-contract=off`).
