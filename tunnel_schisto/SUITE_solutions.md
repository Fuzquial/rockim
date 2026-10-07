# Suite « traversée du litage » — decks, ordre, coût, verdicts attendus

**2026-09-02.** Réponse à la remarque « tout est situé dans un plan de schistosité,
rien n'a réussi à traverser » sur `out_lisjak45_4h`. Base théorique et
bibliographique : `REVUE_traversee_litage.md`. Toutes les solutions sont **en
ajout** : le run de base `tunnel_lisjak45_4h.cfg` reste inchangé, chaque deck
n'en diffère que par ses lignes marquées `<<<` plus `writeJointState = true`
(S8). Binaire : `rockim_f2j.exe` (S8 + S5, bit-identique clés absentes, suite
`fast` en cours). Coût : **≈ 4 h 11 par run**, 14 threads, un job à la fois.

## Ce que le run de base a déjà dit (S6, S8 partiel)

| contrôle | résultat | conséquence |
|---|---|---|
| S6 — ℓ_cz = EΓ/σ² contre l'arête (200 mm) | mode I **le long du litage : 179 mm → non résolue** ; en travers : 556 mm (2,8 éléments) ; mode II : résolues | le rapport 17,5 n'agit pas comme un rapport d'énergies au sens de He–Hutchinson mais comme un rapport résistance × ouverture critique. Le maillage fin (t = 0,35, h = 0,12, ~17 h) est le remède : 179/120 = 1,5 élément |
| théorie (revue §2) | Γ_i/Γ_b = 0,057 est 4,4× sous le seuil 1/4 de He–Hutchinson à 90°, 15× à 30° | l'arrêt sur plan est **prédit**, à toute incidence |
| mesures (revue §3.1) | rapports mesurés 0,10–0,30 en G ; **jamais 0,057** | 0,057 est une constante de calibration de Lisjak, dépendante du maillage qui l'a produite |

## Ordre recommandé (arbre de décision de la revue, §5.4)

| # | deck | variable unique | coût | ce que ça tranche |
|---|---|---|---|---|
| 1 | `S4_ratios1.cfg` | loi γ neutre (ratios = 1), TI conservée | 4 h | le losange vient-il de la **loi cohésive** ? |
| 2 | `S4_isotrope_maillage_lite.cfg` | aucune clé `bedding`, maillage lité seul | 4 h | le losange vient-il du **maillage** ? |
| 3 | S8 sur le run de base | post-traitement (`joint_state_stats.py`) — exige un rerun de la base avec `writeJointState` | 4 h (ou lu sur les runs 1–2) | frontière **ouverte/glissante** (arrêt mécanique) ou **fermée-bloquée** (arrêt énergétique) |
| 4 | `S1_G030.cfg` | G ratio 0,30 (K ≈ 3 converti avec E∥/E⊥ = 3) | 4 h | le losange survit-il à un rapport **mesuré** ? |
| 5 | `S3_lambda130.cfg` | λ = 1,3 (σ_v = 3,85 MPa) | 4 h | **seule comparaison à une donnée in situ chiffrée** : Table 2 d'Armand (Bure) ; réoriente-t-il la zone ? |
| 6 | `S1_G010`, `S1_G025`, `S1_G050` | balayage G | 12 h | où se dissout le losange ; 0,25 doit être la bascule à 90° |
| 7 | `S2_ft0246_G030`, `S2_ft045_G0057`, `S2_ft031_G0057` | résistance et ténacité découplées | 12 h | lequel des deux commande ; ft = 0,31 traverse-t-il malgré G = 0,057 ? (régime « résistance », cf. S6) |
| 8 | `S3_lambda077.cfg` | λ = 1/1,3 | 4 h | symétrie de S3a |
| 9 | `S7_ramp024.cfg` | relâchement 3× plus lent | 4 h | piégeage dynamique (`tip_velocity.py` avant/après) |
| 10 | `S9a_muRes040.cfg` | frottement résiduel 0,40 (+ `jointFrictionScaled = 0`, exigé) | 4 h | sensibilité au frottement des plans délaminés |

**Le résultat qui tranche « physique / artefact »** : runs 1 + 2 + 4. Si le
losange survit aux trois, c'est la physique des critères ; s'il tombe à l'un,
c'est l'ingrédient correspondant. **Le résultat qui tranche « run juste / EDZ
juste »** : run 5 contre la Table 2 d'Armand.

## Les outils (`tunnel_schisto/tools/`)

| outil | mesure | exige |
|---|---|---|
| `edz_sectors.py --dip 45` | profondeur d'enveloppe par secteur, corrigée de la paroi ; rapport le long / en travers (edz_metrics est aveugle à 45°) | rien |
| `joint_state_stats.py --dip 45` | S8 : état de contact des plans-frontières, σ_n et τ par état, **événements de traversée** | `writeJointState = true` |
| `tip_velocity.py` | S7 : vitesse de pointe vs c_R | rien |
| `../rockim/rockim_p1/tunnel_edz/tools/block_sizes.py` | tailles de blocs | rien |
| `../rockim/rockim_p1/tunnel_edz/tools/crack_orientation.py` | sélection conjuguée de Mohr-Coulomb (depuis la racine de p1, chemin absolu) | rien |

Grandeurs à extraire systématiquement (revue §5) : rapport le long / en
travers, fraction de joints rompus à > 60° du litage, nombre d'événements de
traversée, état des plans-frontières, vitesse de pointe maximale.

## S5 — lits faibles et lits forts (nouvelles clés, non encore en deck)

`weakPlaneFactor2` + `weakPlaneFrac2` (+ `weakPlaneSeed`) : une fraction des
plans reçoit un second facteur, tirée par index de plan (les trois pendages
partagent la séquence). Chandler 2016 : 5/7 faibles, 2/7 forts. À combiner avec
`weakPlanes` sur le maillage lité **sans** loi γ : c'est la variante « lits
discrets bimodaux, matrice isotrope » de la revue. Deck à écrire après les
runs 1–2, qui disent d'abord si le maillage seul compte.

## Ce que cette suite ne peut pas faire

- **3D** (S11) : le front d'excavation, source des chevrons de Bure, et le
  contournement d'une délamination finie sont hors de portée du 2D.
- **Le temps** : le flambage des plaques à Mont Terri prend des jours et des
  pressions interstitielles ; un FDEM mécanique sec n'y va pas.
- **Les ténacités de l'argilite cible** : aucune mesure par orientation trouvée
  pour l'Opalinus ni le Callovo-Oxfordien ; la fourchette 0,10–0,30 vient des
  shales.

## Lancement

Depuis la racine de `rockim_f2`, un run à la fois, dans **votre** terminal :

```
powershell -ExecutionPolicy Bypass -File tunnel_schisto\run_solutions.ps1 S4_ratios1 S4_isotrope_maillage_lite S1_G030 S3_lambda130
```

(le script refuse de démarrer si un `rockim` tourne déjà, journalise dans
`out_<deck>\run.log`, et enchaîne les decks donnés dans l'ordre).
