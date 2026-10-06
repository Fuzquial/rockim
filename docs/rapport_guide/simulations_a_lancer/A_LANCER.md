# Simulations à lancer pour le rapport-guide

Préparé le 2026-10-07 par la session de rédaction du rapport (`docs/rapport_guide/`).
Toutes les configurations sont prêtes dans `configs/<groupe>/`, les maillages spéciaux dans
`meshes/` (chemins relatifs à la racine du dépôt). Les autres maillages sont générés par
rockim (Voronoï, grilles) ou déjà dans `meshes/` du dépôt.

## Avant de lancer

```bash
mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j && cd ..
python3 tools/verify_suite.py --tier fast      # repère : 50/51 sous Linux (t1_toolcontact_penalty échoue, connu)
```

Lancement par groupe, un calcul à la fois (reprise automatique des calculs finis) :

```bash
bash docs/rapport_guide/simulations_a_lancer/lancer.sh yan 2
```

Sorties dans `sim_out/` (hors git). À commiter : un `RESULTATS.md` par groupe dans ce dossier,
avec les chiffres extraits des journaux et des `history.csv`, et les petites courbes utiles en CSV.

## Priorité 1 — rejouer la reproduction de Yan et al. 2023 (`configs/yan/`)

Motif : la campagne d'août n'a jamais été rejouée depuis les correctifs du 19-20/08 ; la seule
reprise (UCS, 04/10) ne reproduit plus le résultat (51,43 contre 51,07 MPa, +34 % de joints rompus).
Chapitre `r05a_yan.tex`.

| config | essai | à extraire | référence | durée estimée (1 fil) |
|---|---|---|---|---|
| `strip_adap`, `strip_i1`, `strip_i100` | bande pesante | tassement du sommet | −6,0087e-4 m (éq. 20) ; Yan 0,23 %, 208 % à 1E | secondes à minutes |
| `bd_adap` | brésilien adaptatif | pic σt = 2P/(πDt) | août : 5,30 MPa ; Yan : rapport 1,21 au 500E | ~1-5 min |
| `ucs_adap` | compression uniaxiale | UCS, E apparent, joints rompus/insérés | août : 51,1 MPa, 15,5 GPa, 327/1 288 | ~1-5 min |
| `tx_adap_0/10/20/40` | triaxial 2D | σ1 au pic (sur la courbe, pas le max glissant), ajustement φ, c | août : φ = 22,8° ; Yan 22,87° | ~5-20 min chacun |

À faire en plus si le temps le permet : la même compression à maille 0,75 mm (l'article ; les
decks sont à 1,5 mm), et retracer le profil diamétral du brésilien avec la bonne référence
σyy = −(2P/πDt)[4/(1−r²) − 1] (la figure archivée a un facteur 2).

## Priorité 2 — situer le jeu Red Bohus de juillet avec le binaire actuel (`configs/calib/`)

Chapitre `r07a_calib.tex`. `R_ucs`, `R_tens`, `R_tx20`, `R_tx50` : UCS, traction, σ1 au pic à 20 et
50 MPa. Références : essais 126,6 / 10,27 (BTS) / 424,8 / 649,0 MPa ; archive de juillet 141,5 MPa
(UCS), 186 (20), 247 (50). Vérifier au passage le biais ν/(1−ν)σ3 des mors bloqués pendant
`pullDelay` (soustraire σyy au moment où l'axial démarre). ~5-10 min chacun.

## Priorité 3 — impact court contre Hertz (`configs/impact/`)

Déjà rejoués une fois (50 µs, jeu ramené à 0,02 mm) : pic 32,2 kN en élastique contre 26,05 kN
de Hertz ; 0,53-0,66 × Hertz en fissurable. À prolonger à T = 1,2e-4 s pour avoir la décharge, la
durée de contact (Hertz 89,4 µs) et la restitution. ~1-2 h chacun à 1 fil, moins à 2-4 fils.

## Priorité 4 — compléments AbuAisha et tunnel (`configs/abuaisha/`, `configs/tunnel/`)

Déjà rejoués dans la session de rédaction (voir les dossiers et le chapitre correspondant) ;
relancer seulement si un chiffre manque au chapitre, ou pour les variantes à maillage fin (`_m6`).

## Bancs de la campagne V&V (dossier `vv/`), si Gmsh Python est disponible

Les maillages passent par l'API Python de Gmsh (absente du conteneur de la session de rédaction).
Si `python3 -c "import gmsh"` fonctionne :

| banc | commande | durée | apport |
|---|---|---|---|
| P2.2 sphère contre Hertz | `python3 vv/P2_sphere/p2_sphere.py all` | ~2 h (4 × 2 fils) | contact contre Hertz avec critères écrits |
| P3.1 fissure pressurisée de Guo | `python3 vv/P3_fissure_guo/p3_fissure_guo.py all` (sans la maille 1,25 mm) | plusieurs heures | objectivité de la rupture en maillage |
| P1.2-P1.3 élasticité | `python3 vv/P1_elastique/p1_elastique.py all` | < 1 h | patch test et Kirsch |
| P4 flexion et brésilien de Guo | `python3 vv/P4_essais_guo/p4_essais_guo.py all` | ~58 h à 2 fils | trop long ici : poste local |

## Hors de portée d'une session cloud (poste local, 14 fils)

- St Anne prolongé à 500 µs pour le rebond : `configs/stanne2025_rock137_visc0_T500.cfg`.
- St Anne sur le maillage fin `rock073` (251 460 tétraèdres) : convergence.
- Kuru maillage fin jusqu'au retournement (> 255 µs).
