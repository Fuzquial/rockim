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
