# Tunnel EDZ — reproduction de Wang et al. (2024) dans rockim

*Dossier d'étude préparé le 2026-08-17. **Rien ici n'est appliqué au dépôt** :
les deux patchs C++ sont des fichiers texte à coller quand la machine sera
libre. Seul le maillage (fait le même jour) touche `tools/` et `meshes/`.*

Référence : Y. Wang, J. Qiao, S. Zheng, Z. He, Y. Hu, C. Yan, *Application of
FDEM in the study of large deformation mechanisms in deep-buried soft rock
tunnels: a case study*, **Front. Earth Sci. 12:1517816 (2024)**. Code utilisé
par les auteurs : **MultiFracS** (C. Yan) — la même lignée que la loi
`jointSoftening = yan` et l'insertion adaptative déjà portées ici.

---

## 1. Ce qui existe déjà, ce qu'il faut ajouter

| Besoin de l'article | Dans rockim |
|---|---|
| FDEM 2D triangles + joints cohésifs mode I/II/mixte | natif |
| Adoucissement, décharge sécante (leur fig. 2) | `jointSoftening = yan`, `jointShearUnload = origin` |
| Maillage Gmsh gradué, tunnel en fer à cheval | `meshes/tunnel_hs.msh` (fait) |
| **Frontières bloquées en direction normale** | **déjà là** : `lateralRollers` (flancs, `ROLLERX`) + `gripLateralFree` avec `pullV = 0` (haut/bas : `uy = 0`, `ux` libre) |
| Classement traction / cisaillement / mixte | `bmode`, `rnB`, `rsB` déjà écrits dans `fdem_final_joints.csv` → post-traitement Python, **aucun C++** |
| Déplacement max, EDZ, longueur de fissures | post-traitement des VTU + du CSV joints |
| **Contrainte in situ (σ_h, σ_v), λ ≠ 1** | **PATCH 1** |
| **Excavation** | **PATCH 2** |

Deux patchs, une trentaine de lignes de physique en tout.

## 2. Pourquoi ce choix de méthode d'excavation

Wang et al. donnent au noyau une résistance artificielle pendant la mise en
place de l'in situ, puis **réduisent progressivement son module** (*core
modulus reduction*, Farrokh et al. 2006). Transposer ça exigerait les groupes
physiques en 2D (aujourd'hui 3D seulement) **et** un module variable en cours
de calcul.

Le patch 2 fait l'équivalent par le chemin standard des codes de tunnel — la
**méthode convergence-confinement** : le massif est maillé *avec* sa cavité et
pré-contraint ; les faces de la paroi portent donc à t = 0 une traction
déséquilibrée, exactement celle qu'exerçait la roche excavée. On la rétablit
(facteur de relâchement = 1, état initial **rigoureusement en équilibre**),
puis on la fait décroître jusqu'à zéro. Même état initial, même état final,
sans faire varier un module.

La traction appliquée est **σ₀·n** face par face, pas une pression scalaire :
c'est exact pour λ ≠ 1, ce qu'une pression uniforme ne saurait pas faire — et
c'est la moitié de l'article (§5.2).

## 3. Ordre des opérations

```
1.  coller PATCH 1 puis PATCH 2          (voir les deux fichiers .md)
2.  build_tun.cmd                        -> rockim_tun.exe  (~2 min)
3.  python tools/verify_suite.py --exe rockim_tun.exe --tier fast
        => 15/15 attendu, AUCUN repère ne doit bouger : les deux patchs sont
           inertes tant que insituSh/insituSv valent 0 (constitution, I)
4.  V1  contrôle à charge nulle in situ  (~2 min)     configs/verif_zeroload_insitu.cfg
5.  V2  contrôle de Kirsch, élastique    (~10 min)    configs/verif_kirsch_*.cfg
6.  V3  smoke du cas de référence        (~10 min)    configs/tunnel_smoke.cfg
7.  V4  cas de référence sigma0 = 5 MPa  (~30-60 min) configs/tunnel_ref_s5_lam1.cfg
8.  les deux balayages (8 runs)          python tunnel_edz/tools/make_configs.py
```

Compilation sous un **nouveau nom** (`/Fe:rockim_tun.exe`) : l'exe d'un run en
cours est verrouillé par Windows, et on garde l'ancien binaire pour comparer.

## 4. Les quatre contrôles, avec leur critère de réfutation

**V1 — charge nulle in situ** (`verif_zeroload_insitu.cfg`). In situ 5/5 MPa,
`excavStart` placé *après* la fin du run : la cavité n'est jamais relâchée.
Le massif est alors en équilibre exact.
> Attendu : **0 joint cassé**, énergie cinétique résiduelle ~1e-10 J,
> `dampWork ≤ 0`, résidu B4 < 0,01 %, et la jauge `achievedConfinement` doit
> lire **−5,000 MPa**. Tout écart signale une erreur de signe ou de convention
> dans le patch 1 — c'est le test le plus discriminant du dépôt, appliqué ici.

**V2 — solution de Kirsch** (`verif_kirsch_lam1.cfg`, `..._lam05.cfg`). Trou
CIRCULAIRE R = 5 m, matériau rendu incassable (`ft = 1e12`), donc élastique
pur. La contrainte orthoradiale en paroi a une forme fermée :
σ_θ(θ) = (σ_H + σ_V) − 2(σ_H − σ_V)·cos2θ.

| cas | couronne (θ = 90°) | reins (θ = 0°) |
|---|---|---|
| λ = 1 (5/5) | **10 MPa** | **10 MPa** |
| λ = 0,5 (σ_H = 5, σ_V = 10) | 3σ_H − σ_V = **5 MPa** | 3σ_V − σ_H = **25 MPa** |

> Tolérance ±5 % (Kirsch est pour une plaque infinie ; ici a/b = 1/10).
> C'est LE contrôle qui valide d'un coup la pré-contrainte, le relâchement et
> l'anisotropie. `python tunnel_edz/tools/kirsch_check.py out_kirsch_lam1`

**V3/V4 — le cas de l'article.** σ₀ = 5 MPa hydrostatique, matériau de leur
Table 1. Chiffres publiés à retrouver :

| observable | leur valeur |
|---|---|
| rayon d'EDZ | **19 m** (σ₀ = 5 MPa) |
| déplacement max | **0,347 m** |
| hiérarchie des fissures | cisaillement > mixte > traction |
| faciès | X conjugué, spirales logarithmiques |

Balayages : σ₀ = 3/4/5/6/7 MPa → EDZ 11,1 / 16 / 19 / 22 / 22,8 m et
déplacement 0,244 / 0,278 / 0,347 / 0,393 / 0,451 m ; λ = 0,5…1,5 (σ_H = 5 fixé,
σ_V = 5/λ) → l'EDZ passe d'elliptique horizontale à elliptique verticale, la
casse totale décroissant de façon monotone.

## 5. Deux pièges déjà identifiés, traités dans les configs

1. **`crushCap`**. Leurs éléments sont ÉLASTIQUES (leur fig. 1). Le cap
   déviatorique de rockim vaut par défaut 8·cohésion = 6,4 MPa ici, alors que
   le seul champ in situ à λ = 0,5 donne déjà σ_vm ≈ 5,6 MPa et que la
   concentration en paroi le triple : sans lever le cap, tout le pourtour
   plastifierait pour une raison purement numérique. Les configs posent
   `crushCap = 1e12`. (Le patch 1 n'ajoute pas σ₀ à `e.svm`, donc le cap ne
   *voit* pas l'in situ — raison de plus pour le neutraliser explicitement.)
2. **`insertion = adaptive` est requis**, pas seulement souhaitable : les
   joints sont alors *liés* (nœuds partagés exacts) tant que le critère n'est
   pas atteint, si bien que la pré-contrainte d'élément suffit à décrire
   l'état initial. En intrinsèque, chaque joint devrait d'abord se fermer de
   σ₀/p_j pour transmettre l'in situ : transitoire parasite et contrôle à
   charge nulle non exact. Le patch émet un avertissement dans ce cas.

## 6. Budget

dt ≈ dtFactor·h_min/c_p avec h_min = 0,086 m et c_p = 2 149 m/s, divisé par ~2
par la pénalité d'insertion ⇒ **dt ~ 4·10⁻⁶ s**, soit ~60 000 pas pour
T = 0,25 s. À ~0,14 µs/élément/pas (mesure 2D maison), 95 k éléments donnent
**~20-40 min multithread par cas**, ~9 cas = une demi-journée machine. Leur
code est GPU ; à cette taille, en 2D, ça ne change rien.

Repère de cohérence : avec LEUR pénalité (1000 GPa = 100·E, insertion
intrinsèque), le même calcul donne dt ≈ 8·10⁻⁷ s et ~250 000 pas — ils en
annoncent 300 000.

## 7. Contenu du dossier

```
tunnel_edz/
  README.md                     ce fichier
  PATCH_1_insitu.md             pré-contrainte in situ      (3 hunks)
  PATCH_2_excavation.md         relâchement de la cavité    (5 hunks)
  configs/tunnel_ref_s5_lam1.cfg    le cas de référence de l'article
  configs/tunnel_smoke.cfg          le smoke obligatoire
  configs/verif_zeroload_insitu.cfg V1
  configs/verif_kirsch_lam1.cfg     V2a
  configs/verif_kirsch_lam05.cfg    V2b
  configs/sweep/*.cfg (8) + run_sweep.cmd   générés
  tools/make_configs.py         régénère les 8 configs des balayages
  tools/make_circle_mesh.py     maillage à trou circulaire (Kirsch)
  tools/edz_metrics.py          EDZ, classement des fissures, convergence
  tools/kirsch_check.py         sigma_theta en paroi contre la forme fermée
  plot_tunnel_mesh.py           figure de contrôle du maillage
  tunnel_mesh_check.png
```

Maillages déjà produits (dans `meshes/`, hors de ce dossier) :
`tunnel_hs.msh` (94 960 triangles), `tunnel_hs_smoke.msh` (15 674),
`kirsch_r5.msh` (28 824, ~126 éléments sur le pourtour du trou).

## 8. État de vérification — APPLIQUÉ ET VALIDÉ le 2026-08-17

Les deux patchs sont **collés dans le dépôt** et compilés en `rockim_tun.exe`
(`build_tun.cmd`, 58 s). Résultats mesurés :

| contrôle | résultat |
|---|---|
| suite `--tier fast` | **15/15**, valeurs identiques aux références → bit-neutralité des défauts confirmée |
| **V1** charge nulle in situ | **KE = 1,97·10⁻¹⁷ J/m**, 0 arête insérée, 0 joint cassé, dashpot nul, jauge de mors = **5,000 MPa** exactement |
| **V2a** Kirsch λ = 1 | **écart 1,7 %** (PASS à 5 %) ; σ_rr du premier anneau −0,630 mesuré / −0,626 théorique ; paroi extrapolée −9,84 MPa pour −10 |
| **V2b** Kirsch λ = 0,5 | **écart 2,1 %** ; profil angulaire complet reproduit ; paroi extrapolée **−5,51 MPa à la couronne** (théorie −5,49) et **−24,05 aux reins** (théorie −24,57) |

Sélection des faces de cavité, contrôle gratuit : périmètre lu **31,4127 m**
pour 2πR = 31,4159 (Kirsch), et **33,21 m** sur le fer à cheval — cohérent
avec la reconstruction géométrique.

### La leçon payée par V1 : les quatre coins

Le premier V1 donnait **5 185 J/m** d'énergie cinétique parasite. Cause :
`lateralRollers` ne pose `ROLLERX` que sur les nœuds encore `FREE`, or les
rangées haut/bas sont déjà des mors — les **4 nœuds de coin** restaient donc
libres en x tout en portant σ_xx, chacun avec ~7,5 MN de force déséquilibrée.
Correctif retenu : `gripLateralFree = false` (mors encastrés). L'écart à la
lettre de l'article (« normal-fixed ») ne porte alors que sur le degré de
liberté TANGENTIEL des faces haut/bas, à 45 m du tunnel — et **V2 arbitre
quantitativement que c'est sans effet en paroi** (1,7 % et 2,1 %). Le
raffinement propre serait un troisième patch fixant les 4 coins en x ;
non nécessaire au vu de V2.

### Un piège de lecture du bilan B4 en quasi-statique

Quand rien ne bouge, la ligne `integration` (correction leapfrog f²dt²/2m sur
les nœuds contraints, qui portent de très grandes forces) domine le budget et
le `residu` s'affiche en pourcentage d'une échelle quasi nulle — V1 affiche
`[CHECK]` avec KE = 2·10⁻¹⁷ J/m. **En statique, juger sur l'énergie cinétique
et la jauge de contrainte, pas sur le % du résidu B4.**

Restent non exécutés : `edz_metrics.py` (les noms de colonnes viennent de la
lecture du code ; à confirmer au premier run qui casse) et le cas de
référence lui-même.
