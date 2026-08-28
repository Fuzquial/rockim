# BILAN — les cinq réparations du 2026-08-28

Décision explicite de F. Uzquiano (« je veux que tout soit corrigé donc
corrige tout ») après l'inventaire des défauts relevés pendant le chantier
WP6 (audit du code du 28/08 + revue adverse). Chaque réparation change un
comportement précis ; tout le reste est prouvé bit-identique (4 decks
témoins 2D/3D, avec et sans bulkDamage, `cmp` sur toutes les sorties).

## 1. `toolStop` arrête VRAIMENT l'outil (2D + 3D)

Avant : au franchissement de `toolStop`, `toolContact()` faisait un return
sec — le message « OUTIL ARRETE » mentait, `tool_.v` restait intacte et
l'outil poursuivait sa course à vitesse constante à travers la roche, sans
contact. Après : `tool_.v.setZero()` au déclenchement ; `integrate()` (F = 0
faute de contact) le maintient immobile. Message : « OUTIL ARRETE (v = 0) ».
**Décks affectés** : tout deck de tri (protocole Yang étape 1).
**Preuve** : run court, `toolVz = 0` dans history dès t > toolStop.

## 2. `toolImpulseCap` lisible en percussion 2D

Avant : la clé n'était lue que dans la branche PDC du scénario shear —
en percussion 2D elle était ignorée EN SILENCE. `imp2d_panoplie.cfg` posait
1.0 en croyant l'armer : son contrôle (`imp2d_nocap`) ne contrôlait rien,
les deux runs étaient bit-identiques sur cette clé. L'écrêtage lui-même
(site nodeFc) existait déjà ; seul le READ manquait. Après : lecture dans
la branche percussion + bannière. La branche shear non-PDC (flat/disc de
coupe) reste sans read — aucun deck concerné, résidu assumé et documenté.
**Décks affectés** : `tunnel_edz/configs/imp2d_panoplie.cfg` (le cap
devient réellement actif — c'était le but de ce deck).
**Preuve** : cap = 0.02 vs cap = 0 divergent désormais en percussion 2D.

## 3. `meanTensionCapFactor` gardé `!law_` en 3D (symétrie 2D)

Avant : asymétrie — le 2D gardait le cap par `!law_`, le 3D l'appliquait
PAR-DESSUS la contrainte d'une loi MatLaw (double comptage : une loi
possède sa contrainte). Tranché : la sémantique 2D est la bonne, le 3D
s'aligne. `neoHooke_` et `bulkDamage` restent capés (ce ne sont pas des
lois MatLaw ; comportement inchangé).
**Décks affectés** : tout deck fdem3d avec `law = ...` (impact3d_dpdfh).
**Preuve** : `law = dpr` + cap abaissé → PRE clippait (toolFz 10,4 kN),
POST ne clippe plus (12,9 kN) ; sans `law`, bit-identique.

## 4. `gcRestitution` : le défaut devient bruyant

Avant : défaut 0,2 silencieux qui écrase la détente normale du contact
général — contamine toute mesure de restitution/éjection (revue WP6).
Après : bannière systématique `gcRestitution = X (deck | DEFAUT)`.
Aucun changement de comportement ; figé par ailleurs à 0,2 dans les trois
decks de bench_pulverisation.

## 5. `jointDeath = separation` : le piège devient visible

Le défaut reste `separation` — le changer sans données serait un pari, la
variante C de bench_pulverisation est faite pour trancher. En attendant,
toute percussion avec le défaut imprime la notice : sous l'indenteur un
joint écroui en compression ne meurt jamais, le relais contact roche/roche
ne s'engage pas (et `contactResidualMu` n'y a aucun accès). Décision sur le
défaut : APRÈS le banc C.

## Protocole de preuve (exécuté le 28/08, conteneur)

1. Binaire référence = WP6 (74a31bf) ; binaire réparé = ce commit.
2. Bit-identique : p2a/p2b/p3a/p3b (2D/3D × sans/avec bulkDamage, aucune
   des clés réparées) → toutes sorties `cmp`-identiques.
3. Démonstrations d'effet : §1-§3 ci-dessus, chacune avec run avant/après.
4. Suite fast complète sur le binaire réparé (résultat en annexe du commit).
5. Revue adverse (2 relecteurs indépendants) sur le diff avant commit.

## 6. Reprises post-revue (même jour — 2 relecteurs adverses indépendants)

La revue a validé les cinq réparations sur le fond et exigé sept reprises,
toutes appliquées le jour même :
1. **Métrique « tool KE loss » sous toolStop** : le zeroing détruisait par
   décret la KE restante — cliché `toolKEStop_` pris avant le setZero, la
   métrique lit le cliché (2D + 3D).
2. **Artefact anti-rebond de `toolImpulseCap`** : le plafond ∝ |v_outil|
   donne dv/dt ≤ C·v — avec un outil LIBRE le rebond est structurellement
   interdit (défaut 3D préexistant que le 2D rejoint). AVERTISSEMENT
   bruyant aux deux bannières ; le changement de formule (v de référence
   constante) reste une décision ouverte.
3. **Garde `bulkDamage` durcie** : `law = elastic` posé explicitement
   construisait une MatLaw et court-circuitait la pulvérisation EN SILENCE
   (le piège E3 exact) — désormais TOUTE clé `law` avec bulkDamage jette.
4. **Clé `meanTensionCapFactor` posée sous law** : avertissement bruyant
   (règle E3/E6) — et les 11 decks 3D à `law` du dépôt re-basés par
   `meanTensionCapFactor = 0` explicite (l'intention ne dépend plus du
   défaut) : imp3d_tri, v3d_fixed, v3d_adapt, indent3d_ye, indent3d_fin,
   indent3d_yan, indent3d_grad(+d000,d001), impact3d_dpdfh(+_gros).
5. **Bannière gcRestitution honnête** : suffixe « inerte sous contact =
   potential » (gcRest_ n'agit que dans la branche pénalité) — et déplacée
   APRÈS l'affectation de contactPot_ (bug d'ordre attrapé avant build).
6. **Notice jointDeath étendue à `toolShape = none`** (insert maillé — le
   piège mord autant) ; références de lignes remplacées par des ancres
   textuelles.
7. **Collatéral R1 déclaré** : `cut2d_panoplie.cfg` et `cut2d_repos.cfg`
   (shear + toolStop) gèlent désormais aussi la POSITION de l'outil après
   l'arrêt — les colonnes outil de leur history changent, la roche est
   inchangée (contact déjà coupé avant la réparation).
