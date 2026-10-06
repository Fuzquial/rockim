# Résultats — suite de vérification rapide

`python3 tools/verify_suite.py --tier fast`, lancée sans limite de temps le 2026-10-06 de 22 h 19 à
22 h 36 (UTC), session cloud Linux, binaire `f638b0f`. Journal complet : `sim_out/verify_fast.log`
(hors git).

**50/51**, conforme au repère Linux de `A_LANCER.md` (le rejeu du 06/10 avait donné 49/51, maillage T1
introuvable : il est trouvé ici).

Seul échec, déjà connu :

```
FAIL  t1_toolcontact_penalty  28.7s  toolinj = 2.98519, attendu 4.42633 ± 0.0001;
      toolvb = 10.3848, attendu 8.2484 ± 0.001; broken = 3, attendu 2 ± 0
```
